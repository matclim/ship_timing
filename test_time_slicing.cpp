// Standalone unit test for the time_slicing algorithms.
// Compile as:
//   g++ -std=c++23 -O2 -I . time_slicing.cpp test_time_slicing.cpp -o test_time_slicing
// Verifies the two new spread modes (by_count, fill_window) and the per-hit
// reconstruction with realistic SHiP-style numbers.

#include "time_slicing.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cassert>

using namespace time_slicing;

static int n_pass = 0, n_fail = 0;
#define CHECK(cond) do { \
    if (cond) { ++n_pass; } \
    else { ++n_fail; std::fprintf(stderr, "FAIL: %s @ %s:%d\n", #cond, __FILE__, __LINE__); } \
} while (0)

// ----------------------------------------------------------------------------
// Test 1: by_count mode -- N muons should occupy a window of N/rate.
// With N = 10^6 and rate = 3.33e11, window should be ~3 us.
void test_by_count_window_size()
{
    std::printf("\n--- test_by_count_window_size ---\n");
    SpillTruth s;
    s.spill_id = 0;
    s.window_seconds = 1.2;
    s.event_ids.resize(1'000'000);
    for (std::size_t i = 0; i < s.event_ids.size(); ++i) s.event_ids[i] = i;

    spread_config cfg;
    cfg.mode           = spread_mode::by_count;
    cfg.rate_hz        = 4.0e11 / 1.2;            // ~3.33e11
    cfg.window_seconds = 1.2;
    cfg.base_seed      = 42;

    auto out = spread_events_in_time(s, cfg);
    CHECK(out.events.size() == s.event_ids.size());

    double t_min = out.events.front().t_seconds;
    double t_max = out.events.back ().t_seconds;
    double width = t_max - t_min;
    double expected = double(s.event_ids.size()) / cfg.rate_hz;  // ~3 us
    std::printf("  N=%zu  expected window=%.3e s  observed=%.3e s\n",
                s.event_ids.size(), expected, width);
    // For N=10^6 uniform draws in a window of width W, the empirical span
    // covers ~N/(N+1) of W, so it'll be just barely under. Check we're in
    // [0.95*expected, 1.05*expected] -- generous tolerance.
    CHECK(width < expected);
    CHECK(width > 0.99 * expected);
}

// ----------------------------------------------------------------------------
// Test 2: fill_window mode -- N muons should span (close to) the full window.
void test_fill_window_window_size()
{
    std::printf("\n--- test_fill_window_window_size ---\n");
    SpillTruth s;
    s.spill_id = 1;
    s.window_seconds = 1.2;
    s.event_ids.resize(10'000);   // smaller, so test runs fast
    for (std::size_t i = 0; i < s.event_ids.size(); ++i) s.event_ids[i] = i;

    spread_config cfg;
    cfg.mode           = spread_mode::fill_window;
    cfg.window_seconds = 1.2;
    cfg.base_seed      = 42;

    auto out = spread_events_in_time(s, cfg);
    double t_min = out.events.front().t_seconds;
    double t_max = out.events.back ().t_seconds;
    double width = t_max - t_min;
    std::printf("  N=%zu  fill_window=%.3e s  observed=%.3e s\n",
                s.event_ids.size(), cfg.window_seconds, width);
    CHECK(width <= cfg.window_seconds);
    // 10^4 uniform draws should span > 99% of the window
    CHECK(width > 0.99 * cfg.window_seconds);
}

// ----------------------------------------------------------------------------
// Test 3: truth times are sorted ascending.
void test_truth_times_sorted()
{
    std::printf("\n--- test_truth_times_sorted ---\n");
    SpillTruth s;
    s.spill_id = 7;
    s.window_seconds = 1.2;
    s.event_ids.resize(5000);
    for (std::size_t i = 0; i < s.event_ids.size(); ++i) s.event_ids[i] = i;

    for (auto m : {spread_mode::by_count, spread_mode::fill_window}) {
        spread_config cfg;
        cfg.mode           = m;
        cfg.rate_hz        = 4.0e11 / 1.2;
        cfg.window_seconds = 1.2;
        cfg.base_seed      = 42;
        auto out = spread_events_in_time(s, cfg);
        bool sorted_ok = true;
        for (std::size_t i = 1; i < out.events.size(); ++i)
            if (out.events[i].t_seconds < out.events[i-1].t_seconds) {
                sorted_ok = false; break;
            }
        std::printf("  mode=%d  sorted=%s\n",
                    int(m), sorted_ok ? "yes" : "no");
        CHECK(sorted_ok);
    }
}

// ----------------------------------------------------------------------------
// Test 4: the same seed gives the same times; different seeds differ.
void test_seed_reproducibility()
{
    std::printf("\n--- test_seed_reproducibility ---\n");
    SpillTruth s;
    s.spill_id = 3;
    s.window_seconds = 1.2;
    s.event_ids = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

    spread_config cfg;
    cfg.mode = spread_mode::by_count;
    cfg.rate_hz = 1.0e9;
    cfg.window_seconds = 1.2;
    cfg.base_seed = 42;

    auto a = spread_events_in_time(s, cfg);
    auto b = spread_events_in_time(s, cfg);
    bool same = true;
    for (std::size_t i = 0; i < a.events.size(); ++i)
        if (a.events[i].t_seconds != b.events[i].t_seconds) { same = false; break; }
    CHECK(same);

    cfg.base_seed = 43;
    auto c = spread_events_in_time(s, cfg);
    bool different = false;
    for (std::size_t i = 0; i < a.events.size(); ++i)
        if (a.events[i].t_seconds != c.events[i].t_seconds) { different = true; break; }
    CHECK(different);

    std::printf("  same-seed-same-output=%s  diff-seed-diff-output=%s\n",
                same ? "yes" : "no", different ? "yes" : "no");
}

// ----------------------------------------------------------------------------
// Test 5: per-hit Gaussian smear has the right RMS for each region.
void test_hit_smearing_rms()
{
    std::printf("\n--- test_hit_smearing_rms ---\n");
    // Build a degenerate "spill" where every event has truth t = 0,
    // then attach 1 hit per event. Look at the RMS of t_reco per region.
    int const N = 50000;
    SpillTimedEvents timed;
    timed.spill_id = 0;
    timed.events.reserve(N);
    for (int i = 0; i < N; ++i)
        timed.events.push_back(TruthEvent{std::uint64_t(i), 0, 0.0});

    auto hits_for = [&](int region) {
        SpillEventHits eh;
        eh.spill_id = 0;
        eh.hits_per_event.resize(N);
        for (int i = 0; i < N; ++i)
            eh.hits_per_event[i].push_back(HitInfo{region, 1.0, 0, 0});
        return eh;
    };

    resolution_config cfg;
    cfg.sigma_small_seconds = 100e-12;
    cfg.sigma_large_seconds = 200e-12;
    cfg.base_seed = 99;

    auto rms_of = [](SpillReconstructedHits const& r) {
        double s2 = 0.0;
        for (auto const& h : r.hits) s2 += h.t_reco_seconds * h.t_reco_seconds;
        return std::sqrt(s2 / double(r.hits.size()));
    };

    // small tile (region 0): expect rms ~ 100 ps
    auto small_rms = rms_of(reconstruct_hits(timed, hits_for(0), cfg));
    // large tile (region 3): expect rms ~ 200 ps
    auto large_rms = rms_of(reconstruct_hits(timed, hits_for(3), cfg));

    std::printf("  small region (0)  rms = %.3f ps  (expected 100)\n",
                small_rms * 1.0e12);
    std::printf("  large region (3)  rms = %.3f ps  (expected 200)\n",
                large_rms * 1.0e12);
    CHECK(small_rms > 95.0e-12 && small_rms < 105.0e-12);
    CHECK(large_rms > 190.0e-12 && large_rms < 210.0e-12);
}

// ----------------------------------------------------------------------------
// Test 6: the (coarse, fine) packing decodes to a value within half a fine
// tick of t_reco_seconds * 1e9.
void test_encoder_roundtrip()
{
    std::printf("\n--- test_encoder_roundtrip ---\n");
    int const N = 1000;
    SpillTimedEvents timed;
    timed.spill_id = 0;
    timed.events.reserve(N);
    SpillEventHits eh;
    eh.spill_id = 0;
    eh.hits_per_event.resize(N);
    for (int i = 0; i < N; ++i) {
        // Truth times: 1 ns, 2 ns, ..., N ns. Span is well below 1.2 s.
        timed.events.push_back(TruthEvent{std::uint64_t(i), 0, double(i+1) * 1.0e-9});
        eh.hits_per_event[i].push_back(HitInfo{0, 1.0, 0, 0});  // small tile
    }
    resolution_config cfg;
    cfg.sigma_small_seconds = 0.0;   // no smear, exact roundtrip
    cfg.sigma_large_seconds = 0.0;
    cfg.base_seed = 1;

    auto out = reconstruct_hits(timed, eh, cfg);
    double max_residual_ps = 0.0;
    for (auto const& h : out.hits) {
        double decoded_ns = double(h.coarse) * tts::COARSE_NS
                          + double(h.fine)   * tts::FINE_NS;
        double truth_ns = h.t_reco_seconds * 1.0e9;
        double residual_ps = std::fabs(decoded_ns - truth_ns) * 1000.0;
        if (residual_ps > max_residual_ps) max_residual_ps = residual_ps;
    }
    // Half a fine tick is FINE_NS/2 * 1000 ~ 3.05 ps
    std::printf("  max residual = %.3f ps (expected <= ~3.05 ps)\n", max_residual_ps);
    CHECK(max_residual_ps < 3.10);
}

// ----------------------------------------------------------------------------
int main()
{
    test_by_count_window_size();
    test_fill_window_window_size();
    test_truth_times_sorted();
    test_seed_reproducibility();
    test_hit_smearing_rms();
    test_encoder_roundtrip();

    std::printf("\n=== %d passed, %d failed ===\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
