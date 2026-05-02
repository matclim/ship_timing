#ifndef TIME_SLICING_HPP
#define TIME_SLICING_HPP

// ============================================================================
//  Framework-INDEPENDENT data product types and algorithms.
//
//  No #include of any phlex header here -- this file (and time_slicing.cpp)
//  could be reused outside of Phlex without modification, and the algorithms
//  are unit-testable with a plain main().
//
//  PIPELINE SHAPE (chunked, single-spill):
//
//       SpillMetadata               (spill layer, 1 cell)
//             |
//             | unfold_chunks  (one ChunkSpec per chunk)
//             v
//       ChunkSpec                   (chunk layer, N cells)
//             |
//             | load_and_assign_truth  (file I/O + Poisson process)
//             v
//       ChunkBundle                 (chunk layer)
//             |
//             | reconstruct_hits   (Gaussian smear + (12,32,12)-bit pack)
//             v
//       ChunkReconstructedHits      (chunk layer)
//             |
//             v
//       form_module                 (writes to ROOT)
// ============================================================================

#include <cstdint>
#include <string>
#include <vector>

namespace time_slicing {

// ---- Encoding constants ---------------------------------------------------
namespace tts {
inline constexpr double COARSE_NS = 25.0;
inline constexpr int    FINE_BITS = 12;
inline constexpr int    FINE_DIV  = 1 << FINE_BITS;          // 4096
inline constexpr double FINE_NS   = COARSE_NS / FINE_DIV;    // ~6.10 ps
}

// ============================================================================
//  Spill-level (1 cell per spill)
// ============================================================================

// All the information needed to deterministically split one spill into
// chunks. Produced once by the source provider; consumed by the unfold
// algorithm that emits one ChunkSpec per chunk.
struct SpillMetadata {
    std::string   muon_hits_file;        // path to muon_hits.root
    std::string   treename;              // typically "Events"
    std::uint32_t spill_id{};

    // Truth-time spreading parameters (the unfold needs these to compute
    // each chunk's time window):
    std::string   spread_mode;           // "by_count" or "fill_window"
    double        rate_hz{};             // for by_count
    double        window_seconds{};      // physical spill duration (1.2s ± jitter)

    // File geometry:
    std::uint64_t total_muons{};         // ttree row count
    std::uint64_t events_per_chunk{};    // rows per chunk

    // Reproducibility:
    std::uint64_t base_seed{};
};

// ============================================================================
//  Chunk-level data products
// ============================================================================

// One row in the unfold sequence. Self-contained: a chunk transform reading
// this product knows everything it needs to do its work without having to
// reach back to the SpillMetadata.
struct ChunkSpec {
    std::uint32_t spill_id{};
    std::uint32_t chunk_index{};
    std::uint32_t n_chunks_total{};

    // File access:
    std::string   muon_hits_file;
    std::string   treename;
    std::uint64_t row_start{};           // first row this chunk owns
    std::uint64_t row_stop{};            // one past the last row

    // Time window for this chunk's truth times. Continuous-time:
    // truth times will fall in [t_start, t_stop). Adjacent chunks are
    // exactly contiguous; per-hit Gaussian smearing later may push reco
    // times slightly across boundaries, which is the desired behaviour.
    double        t_start{};
    double        t_stop{};

    std::uint64_t base_seed{};
};

// One tile hit, as read out of muon_hits.root.
struct HitInfo {
    int    region{};            // 0/1 = small, 2/3/4 = large
    double edep_MeV{};
    int    tileX{};
    int    tileY{};
};

// One Geant4 event after time assignment (truth-level).
struct TruthEvent {
    std::uint64_t event_id{};
    std::uint32_t spill_id{};
    double        t_seconds{};
};

// Reconstructed (smeared, packed) hit. Written to ROOT.
struct ReconstructedHit {
    std::uint64_t event_id{};
    int           hit_index{};
    int           region{};
    double        t_truth_seconds{};
    double        t_reco_seconds{};
    std::uint16_t spill_id{};   // 12 bits used
    std::uint32_t coarse{};     // 32 bits, units 25 ns
    std::uint16_t fine{};       // 12 bits, units 25/4096 ns
    double        edep_MeV{};
    int           tileX{};
    int           tileY{};
};

// Bundle returned by load_and_assign_truth. Holds per-chunk hits and
// timed events; one product avoids a multi-output transform (Phlex's
// transform node is single-output by current API).
struct ChunkBundle {
    std::uint32_t                     spill_id{};
    std::uint32_t                     chunk_index{};
    std::vector<std::uint64_t>        event_ids;          // raw row indices
    std::vector<TruthEvent>           timed_events;       // event_ids[i] -> timed[i]
    std::vector<std::vector<HitInfo>> hits_per_event;     // lockstep with event_ids
};

// Column-store layout for the reconstructed hits. form_module's
// ROOT_TTREE backend can natively serialize std::vector<POD>, but
// not std::vector<UserStruct> (which would require a ROOT dictionary).
// So we flatten one ReconstructedHit per row into parallel vectors.
//
// This is also the schema the downstream join script wants -- one row
// per hit, columns for each field. Saves a Python flattening step.
struct ChunkReconstructedHits {
    std::uint32_t              spill_id{};
    std::uint32_t              chunk_index{};
    std::uint32_t              n_hits{};

    std::vector<std::uint64_t> event_id;
    std::vector<std::int32_t>  hit_index;
    std::vector<std::int32_t>  region;
    std::vector<double>        t_truth_seconds;
    std::vector<double>        t_reco_seconds;
    std::vector<std::uint16_t> packed_spill_id;   // 12 bits used
    std::vector<std::uint32_t> coarse;            // 25 ns ticks
    std::vector<std::uint16_t> fine;              // 25/4096 ns ticks
    std::vector<double>        edep_MeV;
    std::vector<std::int32_t>  tileX;
    std::vector<std::int32_t>  tileY;
};

// ============================================================================
//  Spread mode
// ============================================================================

enum class spread_mode { by_count, fill_window };

inline spread_mode parse_spread_mode(std::string const& s)
{
    if (s == "fill_window") return spread_mode::fill_window;
    return spread_mode::by_count;   // default
}

// ============================================================================
//  Resolution config (per-region 1-sigma timing resolution)
// ============================================================================

struct resolution_config {
    double        sigma_small_seconds{100.0e-12};   // regions 0/1
    double        sigma_large_seconds{200.0e-12};   // regions 2/3/4
    std::uint64_t base_seed          {99};
};

double sigma_for_region(int region, resolution_config const& cfg);

// ============================================================================
//  Chunk-time-window math (used by the unfold)
// ============================================================================

// Compute the global muon window (= span of the truth-time distribution).
// In by_count: total_muons / rate_hz, capped at window_seconds.
// In fill_window: window_seconds.
double compute_global_muon_window(spread_mode      mode,
                                  std::uint64_t   total_muons,
                                  double          rate_hz,
                                  double          window_seconds);

// Given chunk_index and total chunks, compute the [t_start, t_stop) for
// this chunk inside the global muon window.
void compute_chunk_window(std::uint32_t chunk_index,
                          std::uint32_t n_chunks_total,
                          double        global_muon_window,
                          double&       t_start,
                          double&       t_stop);

// ============================================================================
//  Algorithms
// ============================================================================

// Reads `[spec.row_start, spec.row_stop)` from spec.muon_hits_file, draws
// `events_per_chunk` uniform truth times in [spec.t_start, spec.t_stop),
// sorts, and returns everything bundled. The function opens its own TFile;
// concurrent calls on different ChunkSpecs are safe because each call uses
// an independent TFile instance.
ChunkBundle load_and_assign_truth(ChunkSpec const& spec);

// Take a ChunkBundle, smear each tile hit's time by a region-dependent
// Gaussian, and pack into the (12,32,12)-bit layout.
ChunkReconstructedHits reconstruct_hits(ChunkBundle const& chunk,
                                        resolution_config const& cfg);

} // namespace time_slicing

#endif
