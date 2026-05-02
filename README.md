# Phlex time-slicing pipeline (SHiP / UBT)

A Phlex-based pipeline that takes muon hits produced by the UBT Geant4
simulation (`run_ubt_alltile`) and produces per-hit reconstructed
timestamps in the SHiP-experiment "spill→coarse→fine" format
(12+32+12 bit packing with per-region Gaussian smearing). 
The pipeline runs with a conversion of CUDA muons (ask Matei for files
or simulation) or look at [https://github.com/matclim/UBTR-D](https://github.com/matclim/UBTR-D)

The pipeline is **chunked** — a single physical 1.2 s spill is split
into ~155 chunks of one million muons each, processed in parallel by
Phlex's TBB workers, and written to ROOT as one tree row per chunk.

```
┌─────────────────────────────────────────────────────────────────────────┐
│                                                                         │
│   muon_hits.root  ──►  Phlex pipeline  ──►  output.root                 │
│   (Geant4)             (this code)          (chunk-level reco hits)     │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

## Table of contents

1. [What the pipeline does](#what-the-pipeline-does)
2. [Files in this directory](#files-in-this-directory)
3. [Prerequisites](#prerequisites)
4. [Building](#building)
5. [Running](#running)
6. [Output format](#output-format)
7. [Diagnostic plots](#diagnostic-plots)
8. [Tuning the pipeline](#tuning-the-pipeline)
9. [Architecture and data flow](#architecture-and-data-flow)
10. [Known issues and gotchas](#known-issues-and-gotchas)
11. [Possible future work](#possible-future-work)

---

## What the pipeline does

Given an input ROOT file with one row per Geant4 event (a muon
propagated through the UBT tile detector, with vector branches
`tile_region`, `tile_edep`, `tile_tileX`, `tile_tileY` per event), the
pipeline:

1. **Assigns a continuous truth time to each muon** (drawn uniformly
   in the chunk's time slice, sorted to give a Poisson-process
   appearance). Per-chunk PRNG state is seeded by `(base_seed,
   spill_id, chunk_index)` so the simulation is deterministic and
   chunk-independent.

2. **Smears each tile hit's truth time with a region-dependent
   Gaussian.** The default is σ = 100 ps for the small ("Fine") tiles
   (regions 0, 1) and σ = 200 ps for the large ("Coarse") tiles
   (regions 2, 3, 4). Per-hit PRNG state is seeded by `(base_seed,
   spill_id, event_id, hit_index)`, independent of chunking.

3. **Packs each reconstructed time** into a 12-bit spill ID + 32-bit
   coarse counter (units 25 ns) + 12-bit fine counter (units
   25/4096 ns ≈ 6.10 ps). This matches the SHiP front-end electronics
   readout format.

The processing is split into chunks of ~10⁶ muons each so that peak
memory is bounded regardless of input file size. With 4 worker threads,
peak RSS for a 154M-muon file is ≈ 800 MB.

## Files in this directory

```
phlex_time_slicing/
├── README.md                             this file
│
├── CMakeLists.txt                        build everything
├── LinkDef.h                             tells rootcling which classes
│                                         to generate dictionaries for
├── time_slicing.hpp                      framework-independent data
│                                         products and algorithms
├── time_slicing.cpp                      ... their implementation
│                                         (does ROOT I/O for chunk loads)
│
├── muon_hits_source.cpp                  Phlex source plugin: emits
│                                         one SpillMetadata per spill
├── time_module.cpp                       Phlex algorithm plugin: unfold
│                                         spill->chunks, load+timestamp
│                                         each chunk, smear+pack hits
│
├── time-test-from-muon-hits.jsonnet      workflow configuration
│
├── flatten_phlex_output.py               post-processor: per-hit ntuple
│                                         (joins primary metadata,
│                                         strips colons from branch
│                                         names) -- requires UnROOT-
│                                         friendly Phlex output
├── diagnostic_plots.py                   five validation panels
│                                         (works on flattened output)
├── diagnostic_plots.C                    ROOT C++ macro version
│                                         (reads raw Phlex output)
└── diagnostic_plots.jl                   Julia version
                                          (currently has a known
                                          UnROOT.jl limitation, see
                                          "Known issues" below)
```

## Prerequisites

| Component                  | Version tested  | Notes                                         |
|----------------------------|-----------------|-----------------------------------------------|
| Phlex                      | post-2025       | built with `PHLEX_USE_FORM=ON`                |
| GCC                        | 15.x            | C++23 features required by Phlex headers      |
| ROOT                       | 6.x             | C++17 build is fine                           |
| Boost                      | 1.91 (or alias) | json component                                |
| TBB (oneAPI)               | any recent      | provided by Phlex's deps                      |
| spdlog, fmt                | any recent      | provided by Phlex's deps                      |
| Python                     | 3.10+           | uproot, numpy, scipy, matplotlib              |
| Julia                      | 1.10+ optional  | UnROOT, Plots, LsqFit                         |

If your distro packages Boost 1.91 but Phlex was built against 1.90,
the symlink workaround we used during development:

```bash
for sofile in /usr/lib/libboost_*.so.1.91.0; do
    base=$(basename "$sofile" .so.1.91.0)
    target="/usr/lib/${base}.so.1.90.0"
    [ ! -e "$target" ] && sudo ln -s "$sofile" "$target"
done
sudo ldconfig
```

## Building

This directory should be placed somewhere convenient and built with
CMake. Recommended workflow:

```bash
# Copy/clone these files into a project directory, then:
cd /path/to/your/project
mkdir -p build && cd build

# Point CMake at the Phlex install
cmake -DCMAKE_PREFIX_PATH=/path/to/phlex/install ..
make -j4
```

A successful build produces three shared libraries in `build/`:

```
libtime_slicing_dict.so       ROOT dictionary for our data products
libmuon_hits_source.so        the Phlex source plugin
libtime_module.so             the Phlex algorithm plugin
```

Plus the rootcling-generated `G__time_slicing_dict.cxx`,
`G__time_slicing_dict_rdict.pcm`, and `libtime_slicing_dict.rootmap`
files needed for ROOT to auto-load the dictionary when reading the
output file.

### Why the lib prefix and the dictionary library matter

- Phlex's plugin loader (`phlex/app/load_module.cpp`) explicitly looks
  for files named `lib<spec>.so` in `PHLEX_PLUGIN_PATH`. If the plugin
  is built without the `lib` prefix (e.g. via
  `set_target_properties(... PREFIX "")`), the loader won't find it
  and you get `Could not locate library with specification 'X'`. The
  CMakeLists.txt uses `add_library(... MODULE ...)` with no prefix
  override, which gives `lib<name>.so` automatically.

- `form_module` (Phlex's ROOT writer) refuses to serialize any class
  for which `TDictionary::GetDictionary(typeid)` returns null. Our
  data product types (`ChunkReconstructedHits` etc.) are not
  fundamental, so they need a dictionary. `libtime_slicing_dict.so`
  is built specifically to register dictionaries for those types.

- The dictionary library must exist as a real `.so` (not as object
  files baked into another library) so that ROOT's `TInterpreter`
  auto-load mechanism — which reads the `.rootmap` file alongside the
  `.so` — can find it later when consumer code reads the output file.

## Running

### Configure the runtime environment

```bash
cd build

# Plugins live here; Phlex's loader needs them on the path.
export PHLEX_PLUGIN_PATH=$PWD:/path/to/phlex/install/lib

# So ROOT can auto-load our dictionary at file-read time.
export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH
```

If `PHLEX_PLUGIN_PATH` is not set you'll get:

```
PHLEX_PLUGIN_PATH has not been set.
```

If `LD_LIBRARY_PATH` doesn't include the build dir, ROOT will print
warnings like `Error in <TInterpreter::TCling::AutoLoad>: failure
loading library libtime_slicing_dict.so` when reading the output.

### Run the pipeline

The workflow file expects `muon_hits.root` to be in the working
directory. Run from `build/` with that file in the same dir:

```bash
ln -s /path/to/your/muon_hits.root .
phlex -c ../time-test-from-muon-hits.jsonnet
```

Expected output (abbreviated):

```
Using configuration file: ../time-test-from-muon-hits.jsonnet
[info] Number of worker threads: 4
Registering FORM output module...
FORM output module registered successfully

=== FormOutputModule::save_data_products ===
... [one block per chunk, ~155 blocks] ...

--- spill 1  chunk 154  reco hits (108751 total, showing first 4) ---
evID      hit  reg    t_truth (ns)        t_reco  (ns)        coarse  fine ...

[info] Seen layers:
  job
   └ spill: 1
      └ chunk: 203
[info] CPU time: 27.65s  Real time: 9.33s  CPU efficiency: 296.23%
[info] Max. RSS: 775.620 MB
```

Wall time on a 4-core machine for 154M muons: ≈ 10 seconds.

The output is `output.root` in the working directory.

## Output format

`output.root` contains a single TTree named:

```
timing:reconstruct_hit_timestamps
```

(The colon comes from Phlex's plugin-label:algorithm-name naming
convention. ROOT C++ handles this fine; it's awkward in some other
tools — see "Known issues".)

The tree has 155 rows — one per non-empty chunk of the spill. Each row
holds two branches:

| Branch                | Type                                       | What it is                                  |
|-----------------------|--------------------------------------------|---------------------------------------------|
| `index`               | `phlex::data_cell_index`                   | Phlex internal: spill_id and chunk_index   |
| `chunk_reco_hits`     | `time_slicing::ChunkReconstructedHits`     | the actual hits, see below                  |

The `chunk_reco_hits` struct (see `time_slicing.hpp`) is a
column-of-vectors layout:

```cpp
struct ChunkReconstructedHits {
    std::uint32_t              spill_id, chunk_index, n_hits;
    std::vector<std::uint64_t> event_id;
    std::vector<std::int32_t>  hit_index;
    std::vector<std::int32_t>  region;
    std::vector<double>        t_truth_seconds;
    std::vector<double>        t_reco_seconds;
    std::vector<std::uint16_t> packed_spill_id;   // 12 bits used
    std::vector<std::uint32_t> coarse;            // 25 ns ticks
    std::vector<std::uint16_t> fine;              // 25/4096 ns ticks
    std::vector<double>        edep_MeV;
    std::vector<std::int32_t>  tileX, tileY;
};
```

So to access the third hit in the fifth chunk in interactive ROOT:

```cpp
TFile f("output.root");
auto* t = f.Get<TTree>("timing:reconstruct_hit_timestamps");
time_slicing::ChunkReconstructedHits* hits = nullptr;
t->SetBranchAddress("chunk_reco_hits", &hits);
t->GetEntry(5);
std::cout << hits->t_reco_seconds[3] << '\n';
```

ROOT will auto-load `libtime_slicing_dict.so` via the rootmap file
when you do `SetBranchAddress` (provided `LD_LIBRARY_PATH` includes
the build dir).

## Diagnostic plots

There are three diagnostic scripts. They produce the same five panels:

1. Residual `t_reco − t_truth` for **small** tiles — should be a
   Gaussian at 0 with σ ≈ 100 ps.
2. Same for **large** tiles — σ ≈ 200 ps.
3. **Overlay** of both, peak-normalized.
4. **All tiles combined** — looks like a mixture of two Gaussians.
5. **Encoder roundtrip residual**: |decoded_ns − t_reco·1e9| ≤ 3.05 ps
   everywhere, validating the (12+32+12) packing is lossless.

### ROOT C++ macro (recommended)

Reads `output.root` directly. Uses ROOT's dictionary so it handles the
composite struct natively.

```bash
cd build
root -l -b -q '../diagnostic_plots.C("output.root","plots")'
```

Output is in `build/plots/` as PNG and PDF.

### Python (requires the flatten step first)

Two-step: first flatten the per-chunk struct to a per-hit ntuple,
joining primary-particle metadata from `muon_hits.root`:

```bash
cd build
python3 ../flatten_phlex_output.py output.root muon_hits.root joined.root
python3 ../diagnostic_plots.py joined.root plots/
```

The intermediate `joined.root` has a `Hits` tree with simple POD
branches (no colons in names, one row per hit) which is convenient
for further analysis.

### Julia version

Currently has a UnROOT.jl limitation — see "Known issues" below.
If you want to use it, run flatten first, then point at `joined.root`
(updating the script's tree-finding logic accordingly).

## Tuning the pipeline

All knobs are in `time-test-from-muon-hits.jsonnet`:

```jsonnet
{
  driver: { cpp: 'generate_layers',
    layers: {
      spill: { parent: 'job',   total: 1,   starting_number: 1 },
      chunk: { parent: 'spill', total: 200, starting_number: 0 },
    },
  },

  sources: {
    muon_hits_source: {
      cpp: 'muon_hits_source',
      layer: 'spill',
      muon_hits_file: 'muon_hits.root',         // input file
      treename: 'Events',                        // input TTree name
      events_per_chunk: 1000000,                 // 1M muons/chunk
      spread_mode: 'by_count',                   // or 'fill_window'
      rate_hz: 3.33333333e11,                    // 4e11 / 1.2 s
      window_seconds: 1.2,                       // physical spill duration
      window_jitter_seconds: 0.0,
      base_seed: 12345,
    },
  },

  modules: {
    timing: {
      cpp: 'time_module',
      spill_layer: 'spill',
      chunk_layer: 'chunk',

      sigma_small_ps: 100,                       // regions 0/1
      sigma_large_ps: 200,                       // regions 2/3/4
      reso_base_seed: 99,

      print_first_n: 4,                          // hits printed per chunk
                                                 // by the observer
    },

    output: {
      cpp: 'form_module',
      products: ['chunk_reco_hits'],
    },
  },
}
```

### Spread modes

- **`by_count`** (default): the global muon window is `total_muons /
  rate_hz`, capped at `window_seconds`. Each chunk gets a contiguous
  time slice. Models a real Poisson process at the configured rate.
- **`fill_window`**: spread all muons uniformly across
  `window_seconds`. Useful for stress-testing chunk count vs window
  width.

### Memory tuning

If you run out of memory, lower `events_per_chunk`. Peak RSS scales
roughly as `n_threads × events_per_chunk × ~150 bytes`, so 1M
events/chunk × 4 threads ≈ 600 MB.

If you have plenty of memory and want fewer chunks, raise
`events_per_chunk`. The total chunk count is
`ceil(total_muons / events_per_chunk)`; the workflow's `chunk: total`
must be at least that. The default of 200 covers up to 200M muons.

### Reproducibility

All randomness flows from `base_seed`. Changing it gives an entirely
new realization; keeping it the same gives bit-identical output across
runs (the algorithms use counter-based seeding so parallelism does not
introduce nondeterminism).

## Architecture and data flow

```
                           ┌────────────────────┐
                           │   SpillMetadata    │  spill layer (1 cell)
                           │ (file path, total  │
                           │  muons, seeds)     │
                           └─────────┬──────────┘
                                     │
                          unfold_chunks (in time_module.cpp)
                                     │
                                     ▼
                           ┌────────────────────┐
                           │     ChunkSpec      │  chunk layer (~155 cells)
                           │ (row range, time   │
                           │  window, seeds)    │
                           └─────────┬──────────┘
                                     │
                       load_and_assign_truth (file I/O, seeds time)
                                     │
                                     ▼
                           ┌────────────────────┐
                           │    ChunkBundle     │  chunk layer
                           │  (event_ids,       │
                           │   timed_events,    │
                           │   hits_per_event)  │
                           └─────────┬──────────┘
                                     │
                       reconstruct_hit_timestamps (smear + pack)
                                     │
                                     ▼
                           ┌────────────────────┐
                           │ChunkReconstructed- │  chunk layer
                           │      Hits          │
                           │(column-vectors)    │
                           └─────────┬──────────┘
                                     │
                                form_module
                                     │
                                     ▼
                              output.root
```

### Layer hierarchy

```
job
 └── spill (1 cell — represents the physical 1.2 s spill)
      └── chunk (1 to N cells — implementation detail for memory bounding)
```

The chunking is **not physical**. There is one spill. Chunks exist
solely so that no single algorithm has to hold all 154M muons in
memory at once. Phlex's TBB graph parallelizes across chunks and the
peak in-flight memory is bounded by `n_threads × per_chunk_size`.

### Truth-time math

For mode `by_count`:

- The global muon window is `T_global = total_muons / rate_hz`,
  clamped at `window_seconds`. With 154M muons at 3.33×10¹¹ Hz this
  gives ~462 µs.
- The window is split into `N` contiguous time slices of width
  `T_global / N`. Chunk `i` claims `[i × T_global/N, (i+1) × T_global/N)`.
- Within its slice, each chunk draws `events_per_chunk` uniform doubles,
  sorts them, and assigns one to each muon in row order.

Because chunk slices tile the global window with no gaps, concatenating
all per-chunk truth times gives exactly the same distribution that one
single non-chunked pass would have produced (modulo PRNG-state
differences).

### Per-hit reco smearing

Each tile hit's seed is `mix(base_seed, spill_id, event_id, hit_index)`.
This is independent of `chunk_index`, so the same hit gets the same
smear regardless of which chunk it ends up in. A hit's reco time can
be pushed across a chunk boundary by smearing — that's expected and
physically correct (chunks are not physical).

## Known issues and gotchas

### 1. Tree name contains a colon

The output tree is named `timing:reconstruct_hit_timestamps`. The colon
is Phlex's separator between the workflow's plugin label (`timing`)
and the algorithm name (`reconstruct_hit_timestamps`). Most tools
handle it, but:

- **ROOT C++**: works natively. `f.Get<TTree>("timing:reconstruct_hit_timestamps")` is fine.
- **uproot (Python)**: works fine; the tree name is just a string.
- **UnROOT (Julia)**: works for the tree itself
  (`f["timing:reconstruct_hit_timestamps"]` returns a `TTree`), but
  see issue #2.

### 2. UnROOT.jl cannot read the composite ChunkReconstructedHits struct

The `chunk_reco_hits` branch is a `TBranchElement` holding the whole
struct serialized via the ROOT dictionary. UnROOT.jl's automatic type
detection (`streamerfor`) returns -1 for our custom type and throws
`BoundsError` because no streamer info is registered.

Two paths to a working Julia analysis:

- **Workaround**: run `flatten_phlex_output.py` first to produce a
  `joined.root` with simple POD branches; UnROOT handles those fine.
- **Proper fix**: split `ChunkReconstructedHits` into separate output
  products inside `time_module.cpp` — each a `std::vector<POD>`. Then
  Phlex writes one branch per column, all of which UnROOT can read.

### 3. form_module also writes intermediate products

When `phlex` runs you'll see `=== FormOutputModule::save_data_products ===`
blocks for `chunk_spec`, `chunk_bundle`, etc., even though only
`chunk_reco_hits` is in the workflow's `products: [...]` list. Those
saves print the message and immediately fail silently for unsupported
types. The actual output file only contains the listed products. This
is verbose but harmless.

### 4. Boost version mismatch

If Phlex was built against Boost 1.90 and your distro upgraded to
1.91, plugin loading fails with `cannot open shared object file:
libboost_json.so.1.90.0`. The `prerequisites` section above shows
the symlink workaround.

### 5. ROOT_USE_FILE downgrades C++ standard

If you add `include(${ROOT_USE_FILE})` to CMakeLists.txt, ROOT's
flags overwrite the C++ standard with whatever ROOT itself was built
against (typically C++17), breaking Phlex's C++23 headers. Don't use
`ROOT_USE_FILE`. Instead, use only the imported targets (`ROOT::Core`,
`ROOT::RIO`, `ROOT::Tree`) and call `ROOT_GENERATE_DICTIONARY` directly
— that's what the current CMakeLists.txt does.

### 6. PHLEX_USE_FORM build flag

If you see `Could not locate library with specification 'form_module'`
at runtime, your Phlex installation was built without
`PHLEX_USE_FORM=ON`. Reconfigure Phlex with:

```bash
cd /path/to/phlex/build
cmake -DPHLEX_USE_FORM=ON \
      -DCMAKE_INSTALL_PREFIX=/path/to/phlex/install \
      -DCMAKE_CXX_FLAGS="-Wno-error=cpp" \
      ..
make -j install
```

The `-Wno-error=cpp` is needed because ROOT's C++17 build emits
`#warning` for some Phlex C++23 features and Phlex builds with
`-Werror=cpp`.

## Possible future work

Things that would make this pipeline more robust or more useful but
aren't currently blocking:

- **Split outputs to fix UnROOT.jl** (see issue #2). The change is
  purely additive — a new transform that takes `ChunkReconstructedHits`
  and emits N parallel `std::vector<POD>` products.

- **Multi-spill workflows.** Currently the layer hierarchy has 1 spill;
  for a long beam-time simulation you'd set
  `spill: { parent: 'job', total: N }` and provide either `N` distinct
  input files or a single input with logical spill subdivision. The
  source plugin needs a small tweak to support this (currently it
  hardcodes one file).

- **Configurable per-region σ.** `sigma_small_ps` and `sigma_large_ps`
  cover the two-region case. If the SHiP front-end ends up with
  per-region resolution (5 distinct sigmas), generalize
  `resolution_config` and `sigma_for_region`.

- **Schema evolution.** The current dictionary is regenerated per
  build; old `output.root` files cannot be read by newer builds if
  `ChunkReconstructedHits`'s layout changes. For long-term archival,
  add a `ClassDef` macro to `ChunkReconstructedHits` with a version
  number and care about backward-compatible layouts.

- **Window jitter.** `window_jitter_seconds` is exposed as a knob but
  trivially equals 0 by default and isn't physics-motivated. Either
  use it (and document the model) or remove it.
