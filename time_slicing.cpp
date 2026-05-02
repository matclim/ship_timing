// Implementation of the framework-INDEPENDENT algorithms declared in
// time_slicing.hpp. Pure C++ -- no Phlex includes. Unit-testable in
// isolation. ROOT is a dependency only because chunk loading reads
// muon_hits.root directly.

#include "time_slicing.hpp"

#include "TFile.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <random>
#include <stdexcept>

namespace time_slicing {

namespace {

// ---- Seed mixing ---------------------------------------------------------
constexpr std::uint64_t kGolden = 0x9E3779B97F4A7C15ULL;

std::uint64_t mix(std::uint64_t base, std::uint64_t a)
{
    return base ^ (a * kGolden);
}

std::uint64_t mix(std::uint64_t base, std::uint64_t a, std::uint64_t b)
{
    return mix(mix(base, a), b);
}

std::uint64_t mix(std::uint64_t base, std::uint64_t a,
                  std::uint64_t b, std::uint64_t c)
{
    return mix(mix(mix(base, a), b), c);
}

// ---- (12,32,12) packing --------------------------------------------------
void pack_into(double t_seconds, std::uint16_t spill_id_12bit,
               std::uint16_t& spill_out,
               std::uint32_t& coarse_out,
               std::uint16_t& fine_out)
{
    using namespace tts;
    double const t_ns = t_seconds * 1.0e9;
    double const ticks_d = std::llround(t_ns / FINE_NS);
    auto const ticks =
        (ticks_d < 0.0) ? 0ULL : static_cast<std::uint64_t>(ticks_d);
    spill_out  = spill_id_12bit & 0x0FFF;
    coarse_out = static_cast<std::uint32_t>(ticks / FINE_DIV);
    fine_out   = static_cast<std::uint16_t>(ticks % FINE_DIV);
}

} // namespace

// =========================================================================
//  Chunk-time-window math
// =========================================================================

double compute_global_muon_window(spread_mode    mode,
                                  std::uint64_t  total_muons,
                                  double         rate_hz,
                                  double         window_seconds)
{
    if (window_seconds <= 0.0)  return 0.0;

    if (mode == spread_mode::fill_window) {
        return window_seconds;
    }
    // by_count
    if (rate_hz <= 0.0 || total_muons == 0) return window_seconds;
    double const t = static_cast<double>(total_muons) / rate_hz;
    return std::min(t, window_seconds);
}

void compute_chunk_window(std::uint32_t chunk_index,
                          std::uint32_t n_chunks_total,
                          double        global_muon_window,
                          double&       t_start,
                          double&       t_stop)
{
    if (n_chunks_total == 0) {
        t_start = 0.0;
        t_stop  = global_muon_window;
        return;
    }
    double const dt = global_muon_window / static_cast<double>(n_chunks_total);
    t_start = static_cast<double>(chunk_index)     * dt;
    t_stop  = static_cast<double>(chunk_index + 1) * dt;
}

double sigma_for_region(int region, resolution_config const& cfg)
{
    if (region == 0 || region == 1) return cfg.sigma_small_seconds;
    return cfg.sigma_large_seconds;
}

// =========================================================================
//  load_and_assign_truth -- ROOT file I/O + truth-time draws.
//
//  Opens an INDEPENDENT TFile per call so this is safe to call concurrently
//  from multiple Phlex worker threads.
// =========================================================================
ChunkBundle load_and_assign_truth(ChunkSpec const& spec)
{
    ChunkBundle out;
    out.spill_id    = spec.spill_id;
    out.chunk_index = spec.chunk_index;

    // Defensive: if row_start is past EOF, return an empty bundle.
    if (spec.row_stop <= spec.row_start) return out;

    std::unique_ptr<TFile> file{
        TFile::Open(spec.muon_hits_file.c_str(), "READ")};
    if (!file || file->IsZombie())
        throw std::runtime_error(
            "load_and_assign_truth: cannot open '" + spec.muon_hits_file + "'");

    auto* tree = dynamic_cast<TTree*>(file->Get(spec.treename.c_str()));
    if (!tree)
        throw std::runtime_error(
            "load_and_assign_truth: tree '" + spec.treename
            + "' not found in " + spec.muon_hits_file);

    long long const n_total = tree->GetEntries();
    long long const r0 =
        std::min(static_cast<long long>(spec.row_start), n_total);
    long long const r1 =
        std::min(static_cast<long long>(spec.row_stop),  n_total);
    if (r1 <= r0) return out;

    std::vector<int>*    b_region = nullptr;
    std::vector<double>* b_edep   = nullptr;
    std::vector<int>*    b_tileX  = nullptr;
    std::vector<int>*    b_tileY  = nullptr;
    tree->SetBranchAddress("tile_region", &b_region);
    tree->SetBranchAddress("tile_edep",   &b_edep);
    tree->SetBranchAddress("tile_tileX",  &b_tileX);
    tree->SetBranchAddress("tile_tileY",  &b_tileY);

    std::size_t const n_rows = static_cast<std::size_t>(r1 - r0);
    out.event_ids.reserve(n_rows);
    out.hits_per_event.reserve(n_rows);

    for (long long e = r0; e < r1; ++e) {
        tree->GetEntry(e);
        out.event_ids.push_back(static_cast<std::uint64_t>(e));

        std::vector<HitInfo> hits;
        std::size_t const m = b_region->size();
        hits.reserve(m);
        for (std::size_t k = 0; k < m; ++k) {
            hits.push_back(HitInfo{
                (*b_region)[k],
                (*b_edep)  [k],
                (*b_tileX) [k],
                (*b_tileY) [k]
            });
        }
        out.hits_per_event.push_back(std::move(hits));
    }

    // -------- draw truth times in [t_start, t_stop), sort -------------
    std::uint64_t const seed = mix(spec.base_seed,
                                   static_cast<std::uint64_t>(spec.spill_id),
                                   static_cast<std::uint64_t>(spec.chunk_index));
    std::mt19937_64 rng{seed};
    std::uniform_real_distribution<double> u{spec.t_start, spec.t_stop};

    std::vector<double> ts;
    ts.reserve(n_rows);
    for (std::size_t i = 0; i < n_rows; ++i) ts.push_back(u(rng));
    std::sort(ts.begin(), ts.end());

    out.timed_events.reserve(n_rows);
    for (std::size_t i = 0; i < n_rows; ++i) {
        out.timed_events.push_back(
            TruthEvent{out.event_ids[i], spec.spill_id, ts[i]});
    }
    return out;
}

// =========================================================================
//  reconstruct_hits -- per-hit Gaussian smear + packing.
//  Outputs a column-of-vectors layout for ROOT-friendly serialization.
// =========================================================================
ChunkReconstructedHits reconstruct_hits(ChunkBundle const& chunk,
                                        resolution_config const& cfg)
{
    ChunkReconstructedHits out;
    out.spill_id    = chunk.spill_id;
    out.chunk_index = chunk.chunk_index;

    auto const n_events    = chunk.timed_events.size();
    auto const n_hit_lists = std::min(n_events, chunk.hits_per_event.size());

    std::size_t total_hits = 0;
    for (std::size_t i = 0; i < n_hit_lists; ++i)
        total_hits += chunk.hits_per_event[i].size();

    out.n_hits = static_cast<std::uint32_t>(total_hits);
    out.event_id        .reserve(total_hits);
    out.hit_index       .reserve(total_hits);
    out.region          .reserve(total_hits);
    out.t_truth_seconds .reserve(total_hits);
    out.t_reco_seconds  .reserve(total_hits);
    out.packed_spill_id .reserve(total_hits);
    out.coarse          .reserve(total_hits);
    out.fine            .reserve(total_hits);
    out.edep_MeV        .reserve(total_hits);
    out.tileX           .reserve(total_hits);
    out.tileY           .reserve(total_hits);

    std::uint16_t const sp12 =
        static_cast<std::uint16_t>(chunk.spill_id & 0x0FFF);

    for (std::size_t i = 0; i < n_hit_lists; ++i) {
        TruthEvent const& ev   = chunk.timed_events[i];
        auto       const& hits = chunk.hits_per_event[i];

        for (std::size_t j = 0; j < hits.size(); ++j) {
            HitInfo const& h = hits[j];

            std::uint64_t const seed = mix(
                cfg.base_seed,
                static_cast<std::uint64_t>(chunk.spill_id),
                static_cast<std::uint64_t>(ev.event_id),
                static_cast<std::uint64_t>(j));
            std::mt19937_64 rng{seed};

            double const sigma = sigma_for_region(h.region, cfg);
            std::normal_distribution<double> nd{0.0, sigma};
            double const t_reco = ev.t_seconds + nd(rng);

            std::uint16_t spill_packed{};
            std::uint32_t coarse{};
            std::uint16_t fine{};
            pack_into(t_reco, sp12, spill_packed, coarse, fine);

            out.event_id        .push_back(ev.event_id);
            out.hit_index       .push_back(static_cast<std::int32_t>(j));
            out.region          .push_back(h.region);
            out.t_truth_seconds .push_back(ev.t_seconds);
            out.t_reco_seconds  .push_back(t_reco);
            out.packed_spill_id .push_back(spill_packed);
            out.coarse          .push_back(coarse);
            out.fine            .push_back(fine);
            out.edep_MeV        .push_back(h.edep_MeV);
            out.tileX           .push_back(h.tileX);
            out.tileY           .push_back(h.tileY);
        }
    }
    return out;
}

} // namespace time_slicing
