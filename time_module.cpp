// time_module.cpp
//
// Phlex algorithm module. Three algorithms in two layers:
//
//   unfold_chunks            (spill -> chunk)  SpillMetadata -> N x ChunkSpec
//   load_and_assign_truth    (chunk transform) ChunkSpec     -> ChunkBundle
//   reconstruct_hits         (chunk transform) ChunkBundle   -> ChunkRecoHits
//   print_reco_hits          (chunk observer)  ChunkRecoHits
//
// Configuration (from the workflow .jsonnet):
//   spill_layer       parent layer name      [spill]
//   chunk_layer       child layer name       [chunk]
//   sigma_small_ps    regions 0/1            [100]
//   sigma_large_ps    regions 2/3/4          [200]
//   reso_base_seed                           [99]
//   print_first_n     hits to print          [4]
//
// All inputs use the explicit product_query{...} form (no "_in"_layer
// UDLs in scope).

#include "phlex/configuration.hpp"
#include "phlex/module.hpp"
#include "phlex/model/data_cell_index.hpp"

#include "time_slicing.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <utility>

using namespace phlex;
using namespace time_slicing;

namespace {

// =========================================================================
//  Unfold functor
// =========================================================================
//
// State carried between calls is just "next chunk index". The functor
// holds the SpillMetadata to know when to stop and how to compute each
// chunk's parameters.
//
// Total number of chunks N = ceil(total_muons / events_per_chunk).
// Trailing empty chunks (row_start >= total_muons) are not generated.
//
class chunk_emitter {
public:
    explicit chunk_emitter(SpillMetadata const& meta)
        : meta_{meta}
    {
        if (meta_.events_per_chunk == 0) {
            n_chunks_ = 0;
        } else {
            n_chunks_ = static_cast<std::uint32_t>(
                (meta_.total_muons + meta_.events_per_chunk - 1)
                / meta_.events_per_chunk);
        }

        global_window_ = compute_global_muon_window(
            parse_spread_mode(meta_.spread_mode),
            meta_.total_muons,
            meta_.rate_hz,
            meta_.window_seconds);
    }

    // Required interface (mirrors test/unfold.cpp::iterate_through):
    std::uint32_t initial_value() const { return 0u; }

    bool predicate(std::uint32_t i) const { return i < n_chunks_; }

    std::pair<std::uint32_t, ChunkSpec>
    unfold(std::uint32_t i, data_cell_index const& /*lid*/) const
    {
        ChunkSpec s;
        s.spill_id       = meta_.spill_id;
        s.chunk_index    = i;
        s.n_chunks_total = n_chunks_;
        s.muon_hits_file = meta_.muon_hits_file;
        s.treename       = meta_.treename;
        s.row_start      = static_cast<std::uint64_t>(i)
                          * meta_.events_per_chunk;
        s.row_stop       = std::min(meta_.total_muons,
                                    s.row_start + meta_.events_per_chunk);
        compute_chunk_window(i, n_chunks_,
                             global_window_, s.t_start, s.t_stop);
        s.base_seed      = meta_.base_seed;
        return std::make_pair(i + 1, s);
    }

private:
    SpillMetadata meta_;
    std::uint32_t n_chunks_{0};
    double        global_window_{0.0};
};

} // namespace

PHLEX_REGISTER_ALGORITHMS(m, config)
{
    auto const spill_layer = config.get<std::string>("spill_layer", "spill");
    auto const chunk_layer = config.get<std::string>("chunk_layer", "chunk");

    resolution_config rc;
    rc.sigma_small_seconds = config.get<double>(       "sigma_small_ps", 100.0) * 1.0e-12;
    rc.sigma_large_seconds = config.get<double>(       "sigma_large_ps", 200.0) * 1.0e-12;
    rc.base_seed           = config.get<std::uint64_t>("reso_base_seed", 99);

    auto const print_first_n = config.get<std::size_t>("print_first_n", 4);

    // ---- Unfold: SpillMetadata -> N x ChunkSpec --------------------------
    m.unfold<chunk_emitter>(
         "unfold_chunks",
         &chunk_emitter::predicate,
         &chunk_emitter::unfold,
         chunk_layer,
         concurrency::unlimited)
     .input_family(product_query{.creator = "muon_hits",
                                 .layer   = spill_layer,
                                 .suffix  = "spill_metadata"})
     .output_product_suffixes("chunk_spec");

    // ---- Transform: ChunkSpec -> ChunkBundle (file I/O + truth times) ----
    m.transform(
         "load_and_assign_truth",
         load_and_assign_truth,
         concurrency::unlimited)
     .input_family(product_query{.creator = "unfold_chunks",
                                 .layer   = chunk_layer,
                                 .suffix  = "chunk_spec"})
     .output_product_suffixes("chunk_bundle");

    // ---- Transform: ChunkBundle -> ChunkReconstructedHits ----------------
    m.transform(
         "reconstruct_hit_timestamps",
         [rc](ChunkBundle const& chunk) -> ChunkReconstructedHits {
             return reconstruct_hits(chunk, rc);
         },
         concurrency::unlimited)
     .input_family(product_query{.creator = "load_and_assign_truth",
                                 .layer   = chunk_layer,
                                 .suffix  = "chunk_bundle"})
     .output_product_suffixes("chunk_reco_hits");

    // ---- Observer: print a few reconstructed hits per chunk --------------
    m.observe(
         "print_reco_hits",
         [print_first_n](ChunkReconstructedHits const& reco) {
             auto const n     = static_cast<std::size_t>(reco.n_hits);
             auto const shown = std::min(n, print_first_n);
             std::cout << "\n--- spill " << reco.spill_id
                       << "  chunk "  << reco.chunk_index
                       << "  reco hits (" << n
                       << " total, showing first " << shown << ") ---\n";
             std::cout << std::left
                       << std::setw(10) << "evID"
                       << std::setw(5)  << "hit"
                       << std::setw(7)  << "reg"
                       << std::setw(20) << "t_truth (ns)"
                       << std::setw(20) << "t_reco  (ns)"
                       << std::setw(13) << "coarse"
                       << std::setw(7)  << "fine"
                       << std::setw(12) << "edep (MeV)"
                       << '\n';
             for (std::size_t i = 0; i < shown; ++i) {
                 std::cout << std::left
                           << std::setw(10) << reco.event_id[i]
                           << std::setw(5)  << reco.hit_index[i]
                           << std::setw(7)  << reco.region[i]
                           << std::setw(20) << std::setprecision(12)
                                            << reco.t_truth_seconds[i] * 1.0e9
                           << std::setw(20) << std::setprecision(12)
                                            << reco.t_reco_seconds[i]  * 1.0e9
                           << std::setw(13) << reco.coarse[i]
                           << std::setw(7)  << reco.fine[i]
                           << std::setw(12) << std::setprecision(4)
                                            << reco.edep_MeV[i]
                           << '\n';
             }
         },
         concurrency::unlimited)
     .input_family(product_query{.creator = "reconstruct_hit_timestamps",
                                 .layer   = chunk_layer,
                                 .suffix  = "chunk_reco_hits"});
}
