// muon_hits_source.cpp
//
// Phlex provider plugin -- spill layer.
//
// Emits one SpillMetadata per spill cell. SpillMetadata is metadata only:
// the file path, treename, total row count, spill timing parameters, and
// seeds. It does NOT carry the actual hits -- those are read later, in
// parallel, by per-chunk transforms (load_and_assign_truth in
// time_module.cpp).
//
// The file is opened once at provider construction to read the row count.
// Subsequent chunk reads happen in their own TFile instances inside the
// chunk transforms (see time_slicing.cpp).
//
// Configuration (from the workflow .jsonnet):
//   layer                 layer to register on             [spill]
//   muon_hits_file        ROOT file with the muons         (required)
//   treename              tree name                        [Events]
//   events_per_chunk      rows per chunk                   [1000000]
//   window_seconds        physical spill duration [s]      [1.2]
//   window_jitter_seconds 1-sigma jitter on window         [0.0]
//   spread_mode           "by_count" | "fill_window"       [by_count]
//   rate_hz               for by_count                     [3.33333e11]
//   base_seed             root seed for everything         [12345]

#include "phlex/configuration.hpp"
#include "phlex/model/data_cell_index.hpp"
#include "phlex/source.hpp"

#include "time_slicing.hpp"

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <random>
#include <stdexcept>
#include <string>

using namespace phlex;
using time_slicing::SpillMetadata;

namespace {

class MetaProvider {
public:
    MetaProvider(std::string file, std::string tree,
                 std::uint64_t events_per_chunk,
                 std::string spread_mode,
                 double rate_hz,
                 double window_seconds,
                 double window_jitter_seconds,
                 std::uint64_t base_seed)
        : file_(std::move(file))
        , tree_(std::move(tree))
        , events_per_chunk_(events_per_chunk)
        , spread_mode_(std::move(spread_mode))
        , rate_hz_(rate_hz)
        , window_seconds_(window_seconds)
        , window_jitter_seconds_(window_jitter_seconds)
        , base_seed_(base_seed)
    {
        // One-time row count.
        std::unique_ptr<TFile> f{TFile::Open(file_.c_str(), "READ")};
        if (!f || f->IsZombie())
            throw std::runtime_error(
                "muon_hits_source: cannot open '" + file_ + "'");
        auto* t = dynamic_cast<TTree*>(f->Get(tree_.c_str()));
        if (!t)
            throw std::runtime_error(
                "muon_hits_source: tree '" + tree_
                + "' not found in " + file_);
        total_muons_ = static_cast<std::uint64_t>(t->GetEntries());
    }

    SpillMetadata build(std::uint32_t spill_id) const
    {
        SpillMetadata m;
        m.muon_hits_file   = file_;
        m.treename         = tree_;
        m.spill_id         = spill_id;
        m.spread_mode      = spread_mode_;
        m.rate_hz          = rate_hz_;
        m.window_seconds   = jittered_window_(spill_id);
        m.total_muons      = total_muons_;
        m.events_per_chunk = events_per_chunk_;
        m.base_seed        = base_seed_;
        return m;
    }

private:
    double jittered_window_(std::uint32_t spill_id) const
    {
        if (window_jitter_seconds_ <= 0.0) return window_seconds_;
        std::uint64_t seed =
            base_seed_ ^ (static_cast<std::uint64_t>(spill_id)
                          * 0x9E3779B97F4A7C15ULL);
        std::mt19937_64 rng{seed};
        std::normal_distribution<double> nd{
            window_seconds_, window_jitter_seconds_};
        double w = nd(rng);
        if (w < 1.0e-9) w = 1.0e-9;
        return w;
    }

    std::string   file_;
    std::string   tree_;
    std::uint64_t events_per_chunk_;
    std::string   spread_mode_;
    double        rate_hz_;
    double        window_seconds_;
    double        window_jitter_seconds_;
    std::uint64_t base_seed_;
    std::uint64_t total_muons_{};
};

} // namespace

PHLEX_REGISTER_PROVIDERS(m, config)
{
    auto const layer = config.get<std::string>("layer", "spill");

    auto provider = std::make_shared<MetaProvider>(
        config.get<std::string>("muon_hits_file"),
        config.get<std::string>("treename", "Events"),
        config.get<std::uint64_t>("events_per_chunk", 1000000),
        config.get<std::string>("spread_mode", "by_count"),
        config.get<double>("rate_hz", 4.0e11 / 1.2),
        config.get<double>("window_seconds", 1.2),
        config.get<double>("window_jitter_seconds", 0.0),
        config.get<std::uint64_t>("base_seed", 12345)
    );

    m.provide("provide_spill_metadata",
              [provider](data_cell_index const& id) -> SpillMetadata {
                  return provider->build(
                      static_cast<std::uint32_t>(id.number()));
              })
     .output_product({.creator = "muon_hits",
                      .layer   = layer,
                      .suffix  = "spill_metadata"});
}
