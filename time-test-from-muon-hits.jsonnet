// time-test-from-muon-hits.jsonnet
//
// End-to-end timing pipeline driven by run_ubt_alltile output.
//
// LAYER HIERARCHY:
//
//   job
//    └── spill   (1 cell)
//         └── chunk   (200 cells, 1M muons each)
//
// One physical 1.2 s spill, internally chunked for memory bounding.
// Total chunk capacity = 200 * 1_000_000 = 200 M, comfortably above
// our 154 M file. Trailing chunks past the file end are no-ops.
//
// The truth-time math is GLOBALLY CONSISTENT across chunks: the
// global muon window is divided into N contiguous time slices, one
// per chunk. Per-hit Gaussian smearing may push reco times across
// chunk boundaries, which is the desired physical behaviour.
//
// MODE switching: change spread_mode to "fill_window" to spread the
// 154 M muons across the whole 1.2 s spill instead of the natural
// N/rate ~462 us window.

{
  driver: {
    cpp: 'generate_layers',
    layers: {
      spill: { parent: 'job',   total: 1,   starting_number: 1 },
      chunk: { parent: 'spill', total: 200, starting_number: 0 },
    },
  },

  sources: {
    muon_hits_source: {
      cpp: 'muon_hits_source',
      layer: 'spill',
      muon_hits_file: 'muon_hits.root',
      treename: 'Events',
      events_per_chunk: 1000000,
      spread_mode: 'by_count',
      rate_hz: 3.33333333e11,           // 4e11 muons / 1.2 s
      window_seconds: 1.2,
      window_jitter_seconds: 0.0,
      base_seed: 12345,
    },
  },

  modules: {
    timing: {
      cpp: 'time_module',
      spill_layer: 'spill',
      chunk_layer: 'chunk',

      // Per-region 1-sigma timing resolution
      sigma_small_ps: 100,              // regions 0/1 (FineLeft/FineRight)
      sigma_large_ps: 200,              // regions 2/3/4 (Coarse*)
      reso_base_seed: 99,

      print_first_n: 4,
    },

    output: {
      cpp: 'form_module',
      products: ['chunk_reco_hits'],
    },
  },
}
