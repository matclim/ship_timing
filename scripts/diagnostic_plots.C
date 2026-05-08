// diagnostic_plots.C
//
// All diagnostic plots for the Phlex timing pipeline, reading the raw
// output.root directly. Each plot is saved as a standalone PNG and PDF.
//
// Styling rules:
//   - single 1D histogram per panel  -> solid fill
//   - several histograms per panel   -> no fill, linewidth = 3
//   - no transparency, no dashed/dotted lines anywhere
//
// Usage:
//   root -l -b -q 'diagnostic_plots.C("output.root","plots")'
// or interactively:
//   root -l
//   root [0] .L diagnostic_plots.C
//   root [1] diagnostic_plots("output.root", "plots")

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TKey.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TMath.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

// Encoder constants (matches tts:: in C++)
constexpr double COARSE_NS = 25.0;
constexpr int    FINE_BITS = 12;
constexpr int    FINE_DIV  = 1 << FINE_BITS;        // 4096
constexpr double FINE_NS   = COARSE_NS / FINE_DIV;  // ~6.10 ps

// ---- Margins -------------------------------------------------------------
// Generous left margin so y-axis title isn't clipped.
// Generous bottom margin for x-axis title.
// Right margin large enough for COLZ palette when present.
void setup_pad(TCanvas* c, bool has_palette = false)
{
    c->SetLeftMargin  (0.13);
    c->SetBottomMargin(0.13);
    c->SetTopMargin   (0.10);
    c->SetRightMargin (has_palette ? 0.15 : 0.05);
}

// Style title and axis labels uniformly.
void style_axes(TH1* h)
{
    h->GetXaxis()->SetTitleSize  (0.045);
    h->GetXaxis()->SetLabelSize  (0.040);
    h->GetXaxis()->SetTitleOffset(1.20);
    h->GetYaxis()->SetTitleSize  (0.045);
    h->GetYaxis()->SetLabelSize  (0.040);
    h->GetYaxis()->SetTitleOffset(1.40);
    h->GetZaxis()->SetTitleSize  (0.045);
    h->GetZaxis()->SetLabelSize  (0.040);
}

TTree* find_tree(TFile* file)
{
    TIter next(file->GetListOfKeys());
    while (auto* key = (TKey*)next()) {
        TObject* obj = key->ReadObj();
        if (auto* t = dynamic_cast<TTree*>(obj)) return t;
    }
    return nullptr;
}

std::string find_composite_branch(TTree* tree)
{
    std::string fallback;
    TIter biter(tree->GetListOfBranches());
    while (auto* b = (TBranch*)biter()) {
        std::string const cls  = b->GetClassName() ? b->GetClassName() : "";
        std::string const name = b->GetName();
        if (name == "index") continue;
        if (cls.find("ChunkReconstructedHits") != std::string::npos) return name;
        fallback = name;
    }
    return fallback;
}

void save_canvas(TCanvas* c, char const* outdir, char const* base)
{
    c->SaveAs(Form("%s/%s.png", outdir, base));
    c->SaveAs(Form("%s/%s.pdf", outdir, base));
}

// =========================================================================
//  Range discovery via a SINGLE FULL-PASS scan with coarse pre-binning.
//
//  This is robust: it never misses data because of an under-sampled pass 1.
//  We accept the cost of a full Draw because the panels we care about
//  (chunk index, coarse counter, tile coordinates) need accurate ranges.
// =========================================================================
struct Range1D { double lo, hi; bool found; };

Range1D scan_range_1d(TTree* tree, std::string const& expr,
                      double pre_lo, double pre_hi, int pre_nbins = 2000)
{
    auto const hname = std::string("htmp_range_") + std::to_string(rand());
    char binning[256];
    std::snprintf(binning, sizeof(binning), "(%d,%g,%g)",
                  pre_nbins, pre_lo, pre_hi);
    std::string const draw_expr = expr + ">>" + hname + binning;
    tree->Draw(draw_expr.c_str(), "", "goff");
    auto* h = (TH1D*)gDirectory->Get(hname.c_str());
    Range1D r{0, 0, false};
    if (h && h->GetEntries() > 0) {
        int const first = h->FindFirstBinAbove(0);
        int const last  = h->FindLastBinAbove(0);
        if (first > 0 && last > 0) {
            r.lo    = h->GetXaxis()->GetBinLowEdge(first);
            r.hi    = h->GetXaxis()->GetBinUpEdge(last);
            r.found = true;
        }
    }
    if (h) delete h;
    return r;
}

struct Range2D { double xlo, xhi, ylo, yhi; bool found; };

Range2D scan_range_2d(TTree* tree, std::string const& expr_y_then_x,
                      double pxlo, double pxhi,
                      double pylo, double pyhi,
                      int nbins = 100)
{
    auto const hname = std::string("htmp_range2d_") + std::to_string(rand());
    char binning[256];
    std::snprintf(binning, sizeof(binning), "(%d,%g,%g,%d,%g,%g)",
                  nbins, pxlo, pxhi, nbins, pylo, pyhi);
    std::string const draw_expr = expr_y_then_x + ">>" + hname + binning;
    tree->Draw(draw_expr.c_str(), "", "goff");
    auto* h = (TH2D*)gDirectory->Get(hname.c_str());
    Range2D r{0, 0, 0, 0, false};
    if (h && h->GetEntries() > 0) {
        // Find first/last filled bins along each axis
        int xfirst = -1, xlast = -1, yfirst = -1, ylast = -1;
        for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
            for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
                if (h->GetBinContent(ix, iy) > 0) {
                    if (xfirst == -1 || ix < xfirst) xfirst = ix;
                    if (ix > xlast) xlast = ix;
                    if (yfirst == -1 || iy < yfirst) yfirst = iy;
                    if (iy > ylast) ylast = iy;
                }
            }
        }
        if (xfirst > 0) {
            r.xlo = h->GetXaxis()->GetBinLowEdge(xfirst);
            r.xhi = h->GetXaxis()->GetBinUpEdge (xlast);
            r.ylo = h->GetYaxis()->GetBinLowEdge(yfirst);
            r.yhi = h->GetYaxis()->GetBinUpEdge (ylast);
            r.found = true;
        }
    }
    if (h) delete h;
    return r;
}

// =========================================================================
//  Two-pass histogram of (t_reco - t_truth) [ps]
// =========================================================================
TH1D* draw_residual_hist(TTree* tree,
                         std::string const& branch_path,
                         std::string const& region_cut,
                         char const* hist_name,
                         char const* hist_title,
                         long long sample_n = 200000)
{
    std::string const expr =
        "(" + branch_path + ".t_reco_seconds - "
            + branch_path + ".t_truth_seconds)*1e12";
    std::string cut;
    if (!region_cut.empty()) {
        std::string const r = branch_path + ".region";
        size_t pos = 0;
        while (pos < region_cut.size()) {
            size_t sep = region_cut.find("||", pos);
            std::string val = (sep == std::string::npos)
                              ? region_cut.substr(pos)
                              : region_cut.substr(pos, sep - pos);
            if (!cut.empty()) cut += "||";
            cut += "(" + r + "==" + val + ")";
            if (sep == std::string::npos) break;
            pos = sep + 2;
        }
    }

    // Range scan: small smearing values, sample-based is OK here.
    tree->Draw((expr + ">>htmp(200,-1e4,1e4)").c_str(),
               cut.c_str(), "goff", sample_n);
    auto* htmp = (TH1D*)gDirectory->Get("htmp");
    double mean = (htmp && htmp->GetEntries() > 0) ? htmp->GetMean() : 0.0;
    double rms  = (htmp && htmp->GetEntries() > 0) ? htmp->GetRMS()  : 100.0;
    double half_range = std::max(5.0 * rms, 50.0);
    if (htmp) delete htmp;

    char binning[128];
    std::snprintf(binning, sizeof(binning), "(121,%g,%g)",
                  mean - half_range, mean + half_range);
    std::string const hexpr = expr + ">>" + hist_name + binning;
    tree->Draw(hexpr.c_str(), cut.c_str(), "goff");

    auto* h = (TH1D*)gDirectory->Get(hist_name);
    if (!h) return nullptr;
    h->SetTitle(hist_title);
    h->GetXaxis()->SetTitle("t_{reco} - t_{truth}  [ps]");
    h->GetYaxis()->SetTitle("hits / bin");
    return h;
}

void fit_and_overlay_gaussian(TH1D* h, double& mu, double& sigma)
{
    auto* g = new TF1(Form("%s_gauss", h->GetName()),
                      "gaus",
                      h->GetXaxis()->GetXmin(),
                      h->GetXaxis()->GetXmax());
    g->SetParameter(0, h->GetMaximum());
    g->SetParameter(1, h->GetMean());
    g->SetParameter(2, h->GetRMS());
    h->Fit(g, "QN");
    mu    = g->GetParameter(1);
    sigma = g->GetParameter(2);
    g->SetLineColor(kBlack);
    g->SetLineWidth(3);
    g->Draw("same");
}

void make_residual_panel(TTree* tree, std::string const& branch,
                         std::string const& region_cut,
                         char const* hist_name,
                         char const* title,
                         int color,
                         char const* outdir,
                         char const* outbase,
                         double target_sigma_ps = -1.0)
{
    auto* c = new TCanvas(Form("c_%s", outbase), "", 800, 600);
    setup_pad(c);
    auto* h = draw_residual_hist(tree, branch, region_cut, hist_name, title);
    if (!h) return;
    style_axes(h);
    h->SetFillColor(color);
    h->SetLineColor(color);
    h->Draw("HIST");
    double mu, sigma;
    fit_and_overlay_gaussian(h, mu, sigma);
    save_canvas(c, outdir, outbase);
}

} // namespace

// ===========================================================================
void diagnostic_plots(char const* in_path = "output.root",
                      char const* outdir  = "plots")
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gSystem->mkdir(outdir, kTRUE);

    auto* file = TFile::Open(in_path, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: cannot open " << in_path << '\n';
        return;
    }

    TTree* tree = find_tree(file);
    if (!tree) {
        std::cerr << "ERROR: no TTree in " << in_path << '\n';
        file->Close(); return;
    }
    std::cout << "Reading " << in_path
              << "  using tree '" << tree->GetName() << "'\n";

    std::string const b = find_composite_branch(tree);
    if (b.empty()) {
        std::cerr << "ERROR: no composite branch found in tree.\n";
        file->Close(); return;
    }
    std::cout << "  using branch: " << b << '\n';
    Long64_t const n_chunks = tree->GetEntries();
    std::cout << "  number of chunk rows: " << n_chunks << '\n';

    // =========================================================================
    //  SMEARING VALIDATION
    // =========================================================================
    make_residual_panel(tree, b, "0||1", "h_small",
                        "Small tiles (regions 0, 1)",
                        kAzure+1, outdir, "residual_small", 100.0);
    make_residual_panel(tree, b, "2||3||4", "h_large",
                        "Large tiles (regions 2, 3, 4)",
                        kRed+1, outdir, "residual_large", 200.0);

    // Overlay
    {
        auto* c = new TCanvas("c_overlay", "", 800, 600);
        setup_pad(c);
        auto* hs = (TH1D*)gDirectory->Get("h_small");
        auto* hl = (TH1D*)gDirectory->Get("h_large");
        if (hs && hl) {
            auto* hsn = (TH1D*)hs->Clone("h_small_norm");
            auto* hln = (TH1D*)hl->Clone("h_large_norm");
            if (hsn->GetMaximum() > 0) hsn->Scale(1.0 / hsn->GetMaximum());
            if (hln->GetMaximum() > 0) hln->Scale(1.0 / hln->GetMaximum());

            hsn->SetFillStyle(0);
            hsn->SetLineColor(kAzure+1);
            hsn->SetLineWidth(3);
            hln->SetFillStyle(0);
            hln->SetLineColor(kRed+1);
            hln->SetLineWidth(3);

            hsn->SetMaximum(1.1);
            hsn->SetTitle("Small vs large tiles, residual distributions;"
                          "t_{reco} - t_{truth}  [ps];"
                          "hits / bin (peak = 1)");
            style_axes(hsn);
            hsn->Draw("HIST");
            hln->Draw("HIST SAME");

            auto* leg = new TLegend(0.65, 0.74, 0.93, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(hsn, "small tiles", "l");
            leg->AddEntry(hln, "large tiles", "l");
            leg->Draw();
            save_canvas(c, outdir, "residual_overlay");
        }
    }

    make_residual_panel(tree, b, "", "h_all",
                        "All tiles combined",
                        kGray+2, outdir, "residual_all", -1.0);

    // =========================================================================
    //  ENCODER VALIDATION
    // =========================================================================

    // Roundtrip
    {
        auto* c = new TCanvas("c_rt", "", 800, 600);
        setup_pad(c);
        char expr[1024];
        std::snprintf(expr, sizeof(expr),
            "(%s.coarse*%g + %s.fine*%g - %s.t_reco_seconds*1e9)*1000.0",
            b.c_str(), COARSE_NS, b.c_str(), FINE_NS, b.c_str());
        std::string const hexpr = std::string(expr) + ">>h_rt(81,-4,4)";
        tree->Draw(hexpr.c_str(), "", "goff");
        auto* h = (TH1D*)gDirectory->Get("h_rt");
        if (h) {
            h->SetTitle("Fine timestamp encoder residual;"
                        "t_{decoded} - t_{true} [ps];"
                        "hits / bin");
            style_axes(h);
            h->SetFillColor(kGreen+2);
            h->SetLineColor(kGreen+2);
            h->Draw("HIST");
            save_canvas(c, outdir, "encoder_roundtrip");
        }
    }

    // Fine counter (range is exactly [0, 4096))
    {
        auto* c = new TCanvas("c_fine", "", 800, 600);
        setup_pad(c);
        std::string const hexpr =
            b + ".fine>>h_fine(128,0,4096)";
        tree->Draw(hexpr.c_str(), "", "goff");
        auto* h = (TH1D*)gDirectory->Get("h_fine");
        if (h) {
            h->SetTitle("Fine counter;"
                        "fine value [0, 4096);"
                        "hits / bin");
            style_axes(h);
            h->SetFillColor(kOrange+1);
            h->SetLineColor(kOrange+1);
            h->SetMinimum(0.0);
            h->Draw("HIST");
            double const mean = h->Integral() / h->GetNbinsX();
            save_canvas(c, outdir, "fine_distribution");
        }
    }

    // Coarse counter -- full-pass range scan, then fill
    {
        auto* c = new TCanvas("c_coarse", "", 800, 600);
        setup_pad(c);
        // Coarse counter is uint32; pre-binning over [0, 2^32) is too coarse.
        // Use [0, 1e8] (covers up to 2.5 s; way more than our 1.2 s window).
        Range1D const r = scan_range_1d(tree, b + ".coarse",
                                         0.0, 1e8, 4000);
        if (!r.found) {
            std::cerr << "  coarse: range scan found nothing\n";
        } else {
            char binning[128];
            std::snprintf(binning, sizeof(binning),
                          "(120,%g,%g)", r.lo, r.hi * 1.02);
            std::string const hexpr = b + ".coarse>>h_coarse" + binning;
            tree->Draw(hexpr.c_str(), "", "goff");
            auto* h = (TH1D*)gDirectory->Get("h_coarse");
            if (h) {
                h->SetTitle("Coarse counter usage;"
                            "coarse value (units of 25 ns);"
                            "hits / bin");
                style_axes(h);
                h->SetFillColor(kViolet+2);
                h->SetLineColor(kViolet+2);
                h->Draw("HIST");
                save_canvas(c, outdir, "coarse_distribution");
            }
        }
    }

    // Truth timestamps in units of 1e-13 s (i.e., t_truth_seconds * 1e13).
    // Should be a flat distribution across the global muon window.
    {
        auto* c = new TCanvas("c_truth", "", 800, 600);
        setup_pad(c);
        // Range scan: very wide initial window to be safe.
        // 1.2 s window -> 1.2e13 in 1e-13 s units. Pre-bin generously.
        Range1D const r = scan_range_1d(tree,
            b + ".t_truth_seconds*1e13",
            0.0, 1.2e13, 4000);
        if (!r.found) {
            std::cerr << "  truth_timestamps: range scan found nothing\n";
        } else {
            char binning[128];
            std::snprintf(binning, sizeof(binning),
                          "(150,%g,%g)", r.lo, r.hi * 1.02);
            std::string const hexpr =
                b + ".t_truth_seconds*1e13>>h_truth" + binning;
            tree->Draw(hexpr.c_str(), "", "goff");
            auto* h = (TH1D*)gDirectory->Get("h_truth");
            if (h) {
                h->SetTitle("Truth timestamps;"
                            "t_{truth}  [10^{-13} s];"
                            "hits / bin");
                style_axes(h);
                h->SetFillColor(kAzure+2);
                h->SetLineColor(kAzure+2);
                h->SetMinimum(0.0);
                h->Draw("HIST");
                save_canvas(c, outdir, "truth_timestamps");
            }
        }
    }

    // =========================================================================
    //  ARCHITECTURE / GEOMETRY
    // =========================================================================

    // Truth time vs chunk index. The chunk index = tree row number
    // (Entry$). Time range we discover via full-pass scan; chunk range
    // is exactly [0, n_chunks).
    {
        auto* c = new TCanvas("c_chunk_truth", "", 900, 700);
        setup_pad(c, true);
        // t_truth in microseconds. Pre-bin generously: [0, 2 s] covers
        // anything reasonable, and we narrow down from there.
        Range1D const tr = scan_range_1d(tree,
            b + ".t_truth_seconds*1e6", 0.0, 2e6, 2000);
        if (!tr.found) {
            std::cerr << "  truth_vs_chunk: range scan found nothing\n";
        } else {
            char binning[256];
            std::snprintf(binning, sizeof(binning),
                          "(200,%g,%g, %lld,0,%lld)",
                          tr.lo, tr.hi * 1.02,
                          (long long)n_chunks, (long long)n_chunks);
            std::string const hexpr =
                "Entry$:" + b + ".t_truth_seconds*1e6"
                ">>h_chunk_truth" + binning;
            tree->Draw(hexpr.c_str(), "", "goff");
            auto* h = (TH2D*)gDirectory->Get("h_chunk_truth");
            if (h) {
                h->SetTitle("Chunk index vs truth time "
                            "(each chunk fills a contiguous time slice);"
                            "t_{truth}  [#mus];"
                            "chunk index");
                style_axes(h);
                h->Draw("COLZ");
                save_canvas(c, outdir, "truth_vs_chunk");
            }
        }
    }

    // Region counts. Range is exactly [-1, 4].
    {
        auto* c = new TCanvas("c_regions", "", 800, 600);
        setup_pad(c);
        std::string const hexpr =
            b + ".region>>h_regions(6,-1.5,4.5)";
        tree->Draw(hexpr.c_str(), "", "goff");
        auto* h = (TH1D*)gDirectory->Get("h_regions");
        if (h) {
            h->SetTitle("Hit count by region;"
                        "region;"
                        "n hits");
            style_axes(h);
            h->SetFillColor(kTeal+3);
            h->SetLineColor(kTeal+3);
            h->Draw("HIST");
            // Bin labels
            h->GetXaxis()->SetBinLabel(h->FindBin(-1.0), "-1 (none)");
            h->GetXaxis()->SetBinLabel(h->FindBin( 0.0), "0 (small)");
            h->GetXaxis()->SetBinLabel(h->FindBin( 1.0), "1 (small)");
            h->GetXaxis()->SetBinLabel(h->FindBin( 2.0), "2 (large)");
            h->GetXaxis()->SetBinLabel(h->FindBin( 3.0), "3 (large)");
            h->GetXaxis()->SetBinLabel(h->FindBin( 4.0), "4 (large)");
            h->GetXaxis()->SetLabelSize(0.045);
            save_canvas(c, outdir, "region_counts");

            std::printf("\n  region counts:\n");
            for (int r = -1; r <= 4; ++r) {
                int const bin = h->FindBin((double)r);
                std::printf("    region %2d: %lld\n",
                            r, (long long)h->GetBinContent(bin));
            }
        }
    }

    // tileX vs tileY heatmap -- full-pass range scan to handle any geometry
    {
        auto* c = new TCanvas("c_tilemap", "", 900, 750);
        setup_pad(c, true);
        // Scan: pre-bin a wide range. Tile coords are integers in some
        // unit; -1000..1000 is a safe initial window.
        Range2D const r = scan_range_2d(tree,
            b + ".tileY:" + b + ".tileX",
            -1000.0, 1000.0, -1000.0, 1000.0, 200);
        if (!r.found) {
            std::cerr << "  tilemap: range scan found nothing\n";
        } else {
            // Pad the range slightly
            double const xpad = (r.xhi - r.xlo) * 0.05 + 1.0;
            double const ypad = (r.yhi - r.ylo) * 0.05 + 1.0;
            double const xlo = r.xlo - xpad, xhi = r.xhi + xpad;
            double const ylo = r.ylo - ypad, yhi = r.yhi + ypad;

            char binning[256];
            std::snprintf(binning, sizeof(binning),
                          "(120,%g,%g, 120,%g,%g)", xlo, xhi, ylo, yhi);
            std::string const hexpr =
                b + ".tileY:" + b + ".tileX>>h_tilemap" + binning;
            tree->Draw(hexpr.c_str(), "", "goff");
            auto* h = (TH2D*)gDirectory->Get("h_tilemap");
            if (h) {
                h->SetTitle("Hit position by tile "
                            "(values from tile_tileX/tile_tileY in muon_hits.root);"
                            "tileX;"
                            "tileY");
                style_axes(h);
                h->Draw("COLZ");
                save_canvas(c, outdir, "tilemap");
            }
        }
    }

    // =========================================================================
    //  BIAS CHECK
    // =========================================================================

    // Per-region mean residual (TProfile)
    {
        auto* c = new TCanvas("c_bias", "", 800, 600);
        setup_pad(c);
        // First range-scan the residual to set Y bounds for the underlying TH2.
        Range1D const yr = scan_range_1d(tree,
            "(" + b + ".t_reco_seconds - " + b + ".t_truth_seconds)*1e12",
            -2000.0, 2000.0, 200);
        double ylo = -1000.0, yhi = 1000.0;
        if (yr.found) {
            ylo = std::min(yr.lo - 50, -50.0);
            yhi = std::max(yr.hi + 50,  50.0);
        }
        char binning[256];
        std::snprintf(binning, sizeof(binning),
                      "(6,-1.5,4.5, 200,%g,%g)", ylo, yhi);
        std::string const expr =
            "(" + b + ".t_reco_seconds - "
                + b + ".t_truth_seconds)*1e12 : "
                + b + ".region";
        std::string const hexpr = expr + ">>h_bias" + binning;
        tree->Draw(hexpr.c_str(), "", "goff");
        auto* h2 = (TH2D*)gDirectory->Get("h_bias");
        if (h2) {
            auto* prof = h2->ProfileX("h_bias_prof");
            prof->SetTitle("Per-region mean residual "
                           "(unbiased smearing #rightarrow 0);"
                           "region;"
                           "mean(t_{reco} - t_{truth})  [ps]");
            style_axes(prof);
            prof->SetMarkerStyle(20);
            prof->SetMarkerSize(1.5);
            prof->SetMarkerColor(kBlack);
            prof->SetLineColor(kBlack);
            prof->SetLineWidth(3);
            prof->Draw("PE");
            // Reference line at 0
            auto* zero = new TLine(prof->GetXaxis()->GetXmin(), 0.0,
                                   prof->GetXaxis()->GetXmax(), 0.0);
            zero->SetLineColor(kRed+1);
            zero->SetLineWidth(3);
            zero->Draw();
            save_canvas(c, outdir, "mean_residual_per_region");
        }
    }

    std::cout << "\nWrote panels to " << outdir << "/\n";
    file->Close();
}
