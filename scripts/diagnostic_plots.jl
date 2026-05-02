#!/usr/bin/env julia
# diagnostic_plots.jl
#
# Validates the timing pipeline output by reading the RAW Phlex output
# directly. No flatten step needed.
#
# The raw output has one TTree row per chunk cell. Each row holds
# parallel column-vectors (event_id, region, t_truth_seconds, ...). The
# branch names contain Phlex's plugin/algorithm prefix, e.g.
#
#     timing:reconstruct_hit_timestamps/chunk_reco_hits/event_id
#
# We discover branches by their final ".../<field>" suffix, read them
# as jagged arrays, flatten them in memory, and produce the diagnostic
# panels:
#
#   1. residual t_reco - t_truth for SMALL tiles (regions 0,1)  -> sigma ~ 100 ps
#   2. same for LARGE tiles (regions 2,3,4)                     -> sigma ~ 200 ps
#   3. both overlaid
#   4. all tiles combined
#   5. encoder roundtrip residual                                -> |r| <= 3.05 ps
#
# Usage:
#     julia diagnostic_plots.jl output.root [outdir]
#
# Dependencies: UnROOT, Plots, LsqFit. Install once with:
#   julia> using Pkg; Pkg.add(["UnROOT", "Plots", "LsqFit"])

using UnROOT
using Plots
using LsqFit
using Statistics
using Printf

# ---- Region classification (matches sigma_for_region in time_slicing.cpp) ----
const SMALL_REGIONS = Set([0, 1])
const LARGE_REGIONS = Set([2, 3, 4])

# ---- Encoder constants (matches tts:: in C++) --------------------------------
const COARSE_NS = 25.0
const FINE_BITS = 12
const FINE_DIV  = 1 << FINE_BITS         # 4096
const FINE_NS   = COARSE_NS / FINE_DIV   # ~6.10 ps

# Fields we need from the per-hit column-vectors. Each must exist as a
# branch whose name ends in "/<field>" or "<sep><field>".
const HIT_FIELDS = ["event_id", "region",
                    "t_truth_seconds", "t_reco_seconds",
                    "coarse", "fine"]

# ============================================================================
#  Tree and branch discovery
# ============================================================================

"""
    find_hits_tree(file::ROOTFile)

Pick the TTree whose branches contain all of HIT_FIELDS as suffixes.
There should be exactly one such tree.
"""
function find_hits_tree(file::ROOTFile)
    candidates = String[]
    for k in keys(file)
        # Strip ROOT cycle suffix ";N"
        name = split(k, ";")[1]
        try
            t = file[name]
            if t isa UnROOT.TTree
                bnames = keys(t)
                if all(any(endswith(b, "/$f") || endswith(b, ".$f")
                           for b in bnames)
                       for f in HIT_FIELDS)
                    push!(candidates, name)
                end
            end
        catch
            # Not a tree, or a tree we can't introspect; skip.
        end
    end
    if isempty(candidates)
        error("No tree found whose branches cover all of $(HIT_FIELDS).\n" *
              "Trees in file: " *
              join(unique([split(k, ";")[1] for k in keys(file)]), ", "))
    end
    if length(candidates) > 1
        error("Multiple trees match (ambiguous): $(candidates).\n" *
              "Adjust this script to point at a specific one.")
    end
    return candidates[1]
end

"""
    locate_branches(tree)

Return a Dict{field -> branch name} for each field in HIT_FIELDS.
"""
function locate_branches(tree)
    bnames = keys(tree)
    result = Dict{String, String}()
    for field in HIT_FIELDS
        matches = [b for b in bnames
                   if endswith(b, "/$field") || endswith(b, ".$field")]
        if length(matches) == 1
            result[field] = matches[1]
        elseif isempty(matches)
            error("No branch ends in /$field. Available branches:\n  " *
                  join(sort(collect(bnames)), "\n  "))
        else
            error("Multiple branches end in /$field: $(matches).")
        end
    end
    return result
end

"""
    flatten_jagged(jagged_arr)

Take an array-of-vectors and concatenate to one flat vector. Returns
the type of the inner element.
"""
function flatten_jagged(jagged)
    total = sum(length, jagged; init=0)
    T = eltype(eltype(jagged))
    out = Vector{T}(undef, total)
    pos = 1
    @inbounds for v in jagged
        n = length(v)
        out[pos:pos+n-1] .= v
        pos += n
    end
    return out
end

# ============================================================================
#  Plot panels
# ============================================================================

gauss_model(x, p) = p[1] .* exp.(-0.5 .* ((x .- p[2]) ./ p[3]).^2)

function fit_gauss(values_ps; nbins::Int=121, half_range_min::Float64=50.0)
    s = std(values_ps)
    m = median(values_ps)
    half_range = max(5*s, half_range_min)
    edges  = range(m - half_range, m + half_range; length=nbins+1)
    counts = zeros(Int, nbins)
    @inbounds for v in values_ps
        i = searchsortedfirst(edges, v) - 1
        1 <= i <= nbins && (counts[i] += 1)
    end
    centers = collect(0.5 .* (edges[1:end-1] .+ edges[2:end]))
    A0  = float(maximum(counts))
    mu0 = float(mean(values_ps))
    sg0 = float(std(values_ps))
    p0  = [A0, mu0, sg0]
    fit_result = try
        curve_fit(gauss_model, centers, float.(counts), p0)
    catch err
        @warn "Gaussian fit failed; using sample stats" exception=err
        return (A0, mu0, sg0, edges, counts)
    end
    return (fit_result.param[1], fit_result.param[2], fit_result.param[3],
            edges, counts)
end

function panel_residual(values_ps, label; target_sigma=nothing)
    n = length(values_ps)
    if n == 0
        return plot(title="$label: NO DATA")
    end
    A, mu, sigma, edges, counts = fit_gauss(values_ps)
    centers = 0.5 .* (edges[1:end-1] .+ edges[2:end])

    p = bar(centers, counts;
            bar_width = step(edges),
            linewidth = 0,
            fillalpha = 0.4,
            color     = :steelblue,
            label     = "data, N=$n",
            xlabel    = "t_reco - t_truth  [ps]",
            ylabel    = "hits / bin",
            title     = label,
            grid      = true,
            legend    = :topright)

    fine = range(first(edges), last(edges); length=600)
    plot!(p, fine, gauss_model(fine, [A, mu, sigma]);
          color=:black, linewidth=1.8,
          label = @sprintf("fit: μ=%.2f ps, σ=%.2f ps", mu, sigma))

    if target_sigma !== nothing
        vline!(p, [-target_sigma, +target_sigma];
               color=:red, linestyle=:dash, alpha=0.5,
               label = @sprintf("target ±1σ = ±%.0f ps", target_sigma))
    end
    return p
end

function panel_overlay(small_ps, large_ps)
    s_l = isempty(large_ps) ? 200.0 : std(large_ps)
    half_range = max(5*s_l, 100.0)
    edges = range(-half_range, half_range; length=121)
    centers = 0.5 .* (edges[1:end-1] .+ edges[2:end])

    p = plot(xlabel = "t_reco - t_truth  [ps]",
             ylabel = "hits / bin (peak = 1)",
             title  = "Small vs large tiles: residual distributions",
             grid   = true,
             legend = :topright)

    if !isempty(small_ps)
        cs = zeros(Int, length(centers))
        @inbounds for v in small_ps
            i = searchsortedfirst(edges, v) - 1
            1 <= i <= length(cs) && (cs[i] += 1)
        end
        cs_n = cs ./ max(maximum(cs), 1)
        plot!(p, centers, cs_n;
              seriestype=:steppost, linewidth=1.7, color=:steelblue,
              label="small tiles  (N=$(length(small_ps)))")
    end
    if !isempty(large_ps)
        cl = zeros(Int, length(centers))
        @inbounds for v in large_ps
            i = searchsortedfirst(edges, v) - 1
            1 <= i <= length(cl) && (cl[i] += 1)
        end
        cl_n = cl ./ max(maximum(cl), 1)
        plot!(p, centers, cl_n;
              seriestype=:steppost, linewidth=1.7, color=:firebrick,
              label="large tiles  (N=$(length(large_ps)))")
    end
    return p
end

function panel_encoder_roundtrip(residual_ps)
    if isempty(residual_ps)
        return plot(title="Encoder roundtrip: NO DATA")
    end
    edges = range(-4, 4; length=81)
    centers = 0.5 .* (edges[1:end-1] .+ edges[2:end])
    counts = zeros(Int, length(centers))
    @inbounds for v in residual_ps
        i = searchsortedfirst(edges, v) - 1
        1 <= i <= length(counts) && (counts[i] += 1)
    end
    max_abs = maximum(abs, residual_ps)
    p = bar(centers, counts;
            bar_width = step(edges),
            linewidth = 0,
            fillalpha = 0.6,
            color     = :seagreen,
            label     = @sprintf("N=%d, max|r|=%.2f ps",
                                 length(residual_ps), max_abs),
            xlabel    = "decoded_ns - t_reco·1e9  [ps]",
            ylabel    = "hits / bin",
            title     = "Encoder roundtrip residual (|r| ≤ 3.05 ps)",
            grid      = true,
            legend    = :topright)
    vline!(p, [-3.05, +3.05]; color=:red, linestyle=:dash, alpha=0.5,
           label="± half fine tick")
    return p
end

# ============================================================================
function main(args)
    if length(args) < 1 || length(args) > 2
        println(stderr,
            "usage: julia diagnostic_plots.jl <phlex_output.root> [outdir]")
        exit(1)
    end
    in_path = args[1]
    outdir  = length(args) == 2 ? args[2] : "plots"
    isdir(outdir) || mkpath(outdir)

    println("Reading $in_path ...")
    file = ROOTFile(in_path)

    tree_name = find_hits_tree(file)
    println("  using tree: $tree_name")
    tree = file[tree_name]

    branch_map = locate_branches(tree)
    println("  found per-hit columns:")
    for (k, v) in sort(collect(branch_map); by=first)
        println("    $k  ->  $v")
    end

    # Read the per-hit columns as jagged arrays (one row per chunk).
    # LazyBranch loads only what we ask for.
    println("\nLoading per-hit columns (jagged: one entry per chunk)...")
    region_jag        = LazyBranch(file, branch_map["region"])
    t_truth_jag       = LazyBranch(file, branch_map["t_truth_seconds"])
    t_reco_jag        = LazyBranch(file, branch_map["t_reco_seconds"])
    coarse_jag        = LazyBranch(file, branch_map["coarse"])
    fine_jag          = LazyBranch(file, branch_map["fine"])

    println("  $(length(region_jag)) chunk rows, flattening...")
    region          = flatten_jagged(region_jag)
    t_truth_seconds = flatten_jagged(t_truth_jag)
    t_reco_seconds  = flatten_jagged(t_reco_jag)
    coarse          = flatten_jagged(coarse_jag)
    fine            = flatten_jagged(fine_jag)
    n = length(region)
    println("  $n hits total")

    # Derived columns
    residual_ps = (Float64.(t_reco_seconds) .- Float64.(t_truth_seconds)) .* 1.0e12
    decoded_t_ns = Float64.(coarse) .* COARSE_NS .+ Float64.(fine) .* FINE_NS
    decoded_resid_ps = (decoded_t_ns .- Float64.(t_reco_seconds) .* 1.0e9) .* 1000.0

    # Region classification
    is_small = [Int(r) in SMALL_REGIONS for r in region]
    is_large = [Int(r) in LARGE_REGIONS for r in region]
    n_small  = count(is_small)
    n_large  = count(is_large)
    n_other  = n - n_small - n_large
    println("  small tiles: $n_small")
    println("  large tiles: $n_large")
    n_other > 0 && println("  unclassified region: $n_other")

    small_ps = residual_ps[is_small]
    large_ps = residual_ps[is_large]

    # ---- Panels --------------------------------------------------------
    p1 = panel_residual(small_ps, "Small tiles (regions 0, 1)";
                        target_sigma = 100.0)
    p2 = panel_residual(large_ps, "Large tiles (regions 2, 3, 4)";
                        target_sigma = 200.0)
    p3 = panel_overlay(small_ps, large_ps)
    p4 = panel_residual(residual_ps, "All tiles combined")
    p5 = panel_encoder_roundtrip(decoded_resid_ps)

    for (name, fig) in (("residual_small",     p1),
                        ("residual_large",     p2),
                        ("residual_overlay",   p3),
                        ("residual_all",       p4),
                        ("encoder_roundtrip",  p5))
        for ext in ("png", "pdf")
            savefig(fig, joinpath(outdir, "$name.$ext"))
        end
    end

    println("\nWrote 5 panels to $outdir/")
    println("  residual_small.png/pdf      -- expected sigma ~100 ps")
    println("  residual_large.png/pdf      -- expected sigma ~200 ps")
    println("  residual_overlay.png/pdf    -- shape comparison")
    println("  residual_all.png/pdf        -- mixture of the two")
    println("  encoder_roundtrip.png/pdf   -- |residual| <= 3.05 ps")
end

main(ARGS)
