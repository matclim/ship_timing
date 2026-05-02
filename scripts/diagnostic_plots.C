// diagnostic_plots.C
//
// Validates the timing pipeline by reading the RAW Phlex output
// directly. The Phlex form_module writes one composite branch per data
// product, where the branch holds the entire ChunkReconstructedHits
// struct. ROOT walks the dictionary (libtime_slicing_dict.so) to
// access the struct's fields.
//
// In TTree::Draw expressions we refer to fields as
//
//     <branch_name>.<field>
//
// where <branch_name> may itself be qualified by the algorithm's
// product path (e.g. "chunk_reco_hits.region"). For Phlex this is
// just the suffix from the algorithm registration -- "chunk_reco_hits".
//
// Five panels:
//   1. residual t_reco - t_truth for SMALL tiles (regions 0,1)  -> sigma ~ 100 ps
//   2. same for LARGE tiles (regions 2,3,4)                     -> sigma ~ 200 ps
//   3. both overlaid
//   4. all tiles combined
//   5. encoder roundtrip residual                                -> |r| <= 3.05 ps
//
// Usage:
//   root -l -b -q 'diagnostic_plots.C("output.root","plots")'

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TBranchElement.h"
#include "TKey.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"

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

// ROOT supports tree names with colons in TTree::Get; we just call it directly.
TTree* find_tree(TFile* file)
{
    TIter next(file->GetListOfKeys());
    while (auto* key = (TKey*)next()) {
        TObject* obj = key->ReadObj();
        if (auto* t = dynamic_cast<TTree*>(obj)) {
            return t;
        }
    }
    return nullptr;
}

// Find the composite branch (the one whose class is our struct).
// Falls back to the last branch found, since typically the tree has
// just `index` (data_cell_index) + `chunk_reco_hits`.
std::string find_composite_branch(TTree* tree)
{
    std::string fallback;
    TIter biter(tree->GetListOfBranches());
    while (auto* b = (TBranch*)biter()) {
        std::string const class_name = b->GetClassName()
                                       ? b->GetClassName() : "";
        std::string const name = b->GetName();
        // Skip the bookkeeping `index` branch (data_cell_index).
        if (name == "index") continue;
        if (class_name.find("ChunkReconstructedHits") != std::string::npos)
            return name;
        fallback = name;
    }
    return fallback;
}

// Build a histogram of (t_reco - t_truth) [ps] using TTree::Draw with
// a lambda-style cut. `region_cut` is something like "==0||==1" applied
// to <branch>.region; pass empty for "all entries".
TH1D* draw_residual_hist(TTree* tree,
                         std::string const& branch_path,    // e.g. "chunk_reco_hits"
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
        // Build (region == 0 || region == 1) etc.
        std::string r = branch_path + ".region";
        // region_cut format example: "0||1" or "2||3||4"
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

    // Pass 1: scan first sample_n entries to set the histogram range.
    long long const n_scanned = tree->Draw(
        (expr + ">>htmp(200,-1e4,1e4)").c_str(),
        cut.c_str(), "goff", sample_n);
    auto* htmp = (TH1D*)gDirectory->Get("htmp");
    double mean = (htmp && htmp->GetEntries() > 0) ? htmp->GetMean() : 0.0;
    double rms  = (htmp && htmp->GetEntries() > 0) ? htmp->GetRMS()  : 100.0;
    double half_range = std::max(5.0 * rms, 50.0);
    if (htmp) delete htmp;

    // Pass 2: real fill.
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
    g->SetLineWidth(2);
    g->Draw("same");
}


void save_canvas(TCanvas* c, char const* outdir, char const* base)
{
    c->SaveAs(Form("%s/%s.png", outdir, base));
    c->SaveAs(Form("%s/%s.pdf", outdir, base));
}

} // namespace

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
        file->Close();
        return;
    }
    std::cout << "Reading " << in_path
              << "  using tree '" << tree->GetName() << "'\n";

    std::string const branch = find_composite_branch(tree);
    if (branch.empty()) {
        std::cerr << "ERROR: no composite branch found in tree.\n";
        file->Close();
        return;
    }
    std::cout << "  using branch: " << branch << '\n';

    // ---- Small tiles -------------------------------------------------------
    {
        auto* c = new TCanvas("c_small", "", 800, 600);
        auto* h = draw_residual_hist(tree, branch, "0||1",
                                     "h_small",
                                     "Small tiles (regions 0, 1)");
        if (h) {
            h->SetLineColor(kBlue+1);
            h->Draw("HIST");
            double mu, sigma;
            fit_and_overlay_gaussian(h, mu, sigma);
            auto* leg = new TLegend(0.58, 0.70, 0.92, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(h, Form("data, N=%lld", (long long)h->GetEntries()), "f");
            leg->AddEntry((TObject*)nullptr,
                          Form("fit: #mu=%.2f, #sigma=%.2f ps", mu, sigma), "");
            leg->AddEntry((TObject*)nullptr, "target #pm1#sigma = #pm100 ps", "");
            leg->Draw();
            save_canvas(c, outdir, "residual_small");
            std::printf("  small tiles: N=%lld, sigma = %.2f ps\n",
                        (long long)h->GetEntries(), sigma);
        }
    }

    // ---- Large tiles -------------------------------------------------------
    {
        auto* c = new TCanvas("c_large", "", 800, 600);
        auto* h = draw_residual_hist(tree, branch, "2||3||4",
                                     "h_large",
                                     "Large tiles (regions 2, 3, 4)");
        if (h) {
            h->SetLineColor(kRed+1);
            h->Draw("HIST");
            double mu, sigma;
            fit_and_overlay_gaussian(h, mu, sigma);
            auto* leg = new TLegend(0.58, 0.70, 0.92, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(h, Form("data, N=%lld", (long long)h->GetEntries()), "f");
            leg->AddEntry((TObject*)nullptr,
                          Form("fit: #mu=%.2f, #sigma=%.2f ps", mu, sigma), "");
            leg->AddEntry((TObject*)nullptr, "target #pm1#sigma = #pm200 ps", "");
            leg->Draw();
            save_canvas(c, outdir, "residual_large");
            std::printf("  large tiles: N=%lld, sigma = %.2f ps\n",
                        (long long)h->GetEntries(), sigma);
        }
    }

    // ---- Overlay -----------------------------------------------------------
    {
        auto* c = new TCanvas("c_overlay", "", 800, 600);
        auto* hs = (TH1D*)gDirectory->Get("h_small");
        auto* hl = (TH1D*)gDirectory->Get("h_large");
        if (hs && hl) {
            auto* hsn = (TH1D*)hs->Clone("h_small_norm");
            auto* hln = (TH1D*)hl->Clone("h_large_norm");
            if (hsn->GetMaximum() > 0) hsn->Scale(1.0 / hsn->GetMaximum());
            if (hln->GetMaximum() > 0) hln->Scale(1.0 / hln->GetMaximum());
            hsn->SetFillStyle(0);  hsn->SetLineColor(kBlue+1);  hsn->SetLineWidth(2);
            hln->SetFillStyle(0);  hln->SetLineColor(kRed+1);   hln->SetLineWidth(2);
            hsn->SetMaximum(1.1);
            hsn->SetTitle("Small vs large tiles, residual distributions;"
                          "t_{reco} - t_{truth}  [ps];"
                          "hits / bin (peak = 1)");
            hsn->Draw("HIST");
            hln->Draw("HIST SAME");
            auto* leg = new TLegend(0.62, 0.74, 0.92, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(hsn, "small tiles", "l");
            leg->AddEntry(hln, "large tiles", "l");
            leg->Draw();
            save_canvas(c, outdir, "residual_overlay");
        }
    }

    // ---- All tiles ---------------------------------------------------------
    {
        auto* c = new TCanvas("c_all", "", 800, 600);
        auto* h = draw_residual_hist(tree, branch, "",
                                     "h_all",
                                     "All tiles combined");
        if (h) {
            h->SetLineColor(kBlack);
            h->Draw("HIST");
            double mu, sigma;
            fit_and_overlay_gaussian(h, mu, sigma);
            auto* leg = new TLegend(0.62, 0.74, 0.92, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(h, Form("data, N=%lld", (long long)h->GetEntries()), "f");
            leg->AddEntry((TObject*)nullptr,
                          Form("fit: #mu=%.2f, #sigma=%.2f ps", mu, sigma), "");
            leg->Draw();
            save_canvas(c, outdir, "residual_all");
        }
    }

    // ---- Encoder roundtrip -------------------------------------------------
    {
        auto* c = new TCanvas("c_rt", "", 800, 600);
        char expr[1024];
        std::snprintf(expr, sizeof(expr),
            "(%s.coarse*%g + %s.fine*%g - %s.t_reco_seconds*1e9)*1000.0",
            branch.c_str(), COARSE_NS,
            branch.c_str(), FINE_NS,
            branch.c_str());

        std::string const hexpr = std::string(expr) + ">>h_rt(81,-4,4)";
        tree->Draw(hexpr.c_str(), "", "goff");
        auto* h = (TH1D*)gDirectory->Get("h_rt");
        if (h) {
            h->SetTitle("Encoder roundtrip residual (|r| #leq 3.05 ps);"
                        "decoded_{ns} - t_{reco}#cdot10^{9}  [ps];"
                        "hits / bin");
            h->SetLineColor(kGreen+2);
            h->Draw("HIST");
            auto* leg = new TLegend(0.55, 0.78, 0.92, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(h, Form("N=%lld", (long long)h->GetEntries()), "f");
            leg->AddEntry(lminus, "#pm half fine tick (3.05 ps)", "l");
            leg->Draw();
            save_canvas(c, outdir, "encoder_roundtrip");
        }
    }

    std::cout << "\nWrote 5 panels to " << outdir << "/\n"
              << "  residual_small.png/pdf      -- expected sigma ~100 ps\n"
              << "  residual_large.png/pdf      -- expected sigma ~200 ps\n"
              << "  residual_overlay.png/pdf    -- shape comparison\n"
              << "  residual_all.png/pdf        -- mixture of the two\n"
              << "  encoder_roundtrip.png/pdf   -- |residual| <= 3.05 ps\n";
    file->Close();
}
