#!/usr/bin/env julia

using CSV
using DataFrames
using JSON
using LaTeXStrings
using Plots
using Statistics

include(joinpath(@__DIR__, "..", "src", "PhoenixFeedback.jl"))
using .PhoenixFeedback

# -- Okabe-Ito colorblind-safe palette (Wong 2011) ----------------------------
const CB_BLACK   = "#000000"
const CB_ORANGE  = "#E69F00"
const CB_SKYBLUE = "#56B4E9"
const CB_GREEN   = "#009E73"
const CB_YELLOW  = "#F0E442"
const CB_BLUE    = "#0072B2"
const CB_VERMIL  = "#D55E00"
const CB_PINK    = "#CC79A7"

# -- ApJ-compliant publication defaults ----------------------------------------
default(
    titlefontsize = 16,
    guidefontsize = 14,
    tickfontsize = 12,
    legendfontsize = 11,
    linewidth = 2,
    grid = false,
    dpi = 300,
    framestyle = :box,
    foreground_color_legend = nothing,
    background_color_legend = nothing,
)

root = normpath(joinpath(@__DIR__, ".."))
mkpath(joinpath(root, "results"))
mkpath(joinpath(root, "figures"))

# -- Cluster paths -------------------------------------------------------------
cluster_files = [
    ("Phoenix", joinpath(root, "data", "assumptions.json")),
    ("Perseus", joinpath(root, "data", "clusters", "perseus.json")),
    ("MS 0735", joinpath(root, "data", "clusters", "ms0735.json")),
    ("Hydra A", joinpath(root, "data", "clusters", "hydra_a.json")),
    ("Abell 2052", joinpath(root, "data", "clusters", "a2052.json")),
]

# -- Cluster colors ------------------------------------------------------------
cluster_colors = Dict(
    "Phoenix"    => CB_BLACK,
    "Perseus"    => CB_ORANGE,
    "MS 0735"    => CB_SKYBLUE,
    "Hydra A"    => CB_GREEN,
    "Abell 2052" => CB_VERMIL,
)

# -- Run all clusters ----------------------------------------------------------
all_summaries = Dict{String, Any}()
comparison_rows = NamedTuple[]
scatter_frames = Dict{String, DataFrame}()

for (tag, path) in cluster_files
    println("Processing $tag ($(relpath(path, root)))...")
    asm = load_assumptions(path)
    samples = run_multi_age_monte_carlo(asm; seed = 42)
    mdf = DataFrame(samples)
    scatter_frames[tag] = mdf
    summary = summarize_multi_age_samples(samples)
    summary["cluster_name"] = asm.cluster_name
    summary["redshift"] = asm.redshift
    summary["n_cavities"] = asm.n_cavities
    all_summaries[tag] = summary

    for key in ["sph_soundcross", "sph_buoyancy", "sph_refill",
                "ell_soundcross", "ell_buoyancy", "ell_refill"]
        s = summary[key]
        push!(comparison_rows, (
            cluster_name = asm.cluster_name,
            redshift = asm.redshift,
            model = s["label"],
            median_ratio = s["median"],
            p16 = s["p16"],
            p84 = s["p84"],
            p_exceeds_1 = s["p_exceeds_1"],
        ))
    end
end

# -- Save summaries ------------------------------------------------------------
open(joinpath(root, "results", "multi_cluster_summary.json"), "w") do io
    JSON.print(io, all_summaries, 4)
end

CSV.write(joinpath(root, "results", "multi_cluster_comparison.csv"), DataFrame(comparison_rows))

# -- P_cav vs L_cool scatter plot ----------------------------------------------
p_scatter = plot(;
    xlabel = L"L_{\rm cool}\ \mathrm{(erg\ s^{-1})}",
    ylabel = L"P_{\rm cav}\ \mathrm{(erg\ s^{-1})}",
    title = "Cavity Power vs. Cooling Luminosity",
    xscale = :log10,
    yscale = :log10,
    legend = :topleft,
    size = (2250, 1500),
)

for (tag, path) in cluster_files
    mdf = scatter_frames[tag]
    clr = cluster_colors[tag]
    scatter!(p_scatter,
        mdf.cooling_luminosity_erg_s,
        mdf.P_sph_buoyancy;
        label = tag,
        color = clr,
        markersize = 2,
        alpha = 0.25,
        markerstrokewidth = 0,
    )
end

savefig(p_scatter, joinpath(root, "figures", "pcav_lcool_correlation.png"))
savefig(p_scatter, joinpath(root, "figures", "pcav_lcool_correlation.pdf"))

# -- Summary output ------------------------------------------------------------
println()
println("===============================================================")
println("  Multi-Cluster Pipeline - Complete")
println("===============================================================")
println()
println("Per-cluster summaries: results/multi_cluster_summary.json")
println("Comparison table:      results/multi_cluster_comparison.csv")
println("Scatter plot:          figures/pcav_lcool_correlation.png")
println()

for (tag, _) in cluster_files
    s = all_summaries[tag]
    println("$tag ($(s["redshift"])):")
    for key in ["sph_soundcross", "sph_buoyancy", "sph_refill",
                "ell_soundcross", "ell_buoyancy", "ell_refill"]
        m = s[key]
        println("  $(rpad(m["label"], 30)) median=$(round(m["median"], digits=3))  " *
                "[$(round(m["p16"], digits=3)), $(round(m["p84"], digits=3))]  " *
                "P>=1=$(round(m["p_exceeds_1"], digits=3))")
    end
    println()
end
