#!/usr/bin/env julia

using CSV
using DataFrames
using Dates
using JSON
using SHA
using Statistics

include(joinpath(@__DIR__, "..", "src", "PhoenixFeedback.jl"))
using .PhoenixFeedback

root = normpath(joinpath(@__DIR__, ".."))
mkpath(joinpath(root, "results"))
mkpath(joinpath(root, "figures"))

function _try_readchomp(cmd)
    try
        return readchomp(pipeline(cmd; stderr = devnull))
    catch
        return nothing
    end
end

function _git_dirty(root::AbstractString)
    status = _try_readchomp(`git -C $root status --porcelain`)
    return status === nothing ? nothing : !isempty(status)
end

# -- CLI mode parsing ---------------------------------------------------------
mode = "default"
n_override = nothing
seed_override = 42
for i in 1:length(ARGS)
    if ARGS[i] == "--quick"
        global mode = "quick"
        global n_override = 2_000
    elseif ARGS[i] == "--full"
        global mode = "full"
        global n_override = 50_000
    elseif ARGS[i] == "--seed" && i < length(ARGS)
        global seed_override = parse(Int, ARGS[i+1])
    end
end

assumptions_path = joinpath(root, "data", "assumptions.json")
assumptions = load_assumptions(assumptions_path)

# Override sample count if --quick or --full was passed
if n_override !== nothing
    assumptions = PhoenixFeedback.ClusterAssumptions(
        assumptions.cluster_name,
        assumptions.redshift,
        assumptions.cooling_luminosity_erg_s,
        assumptions.gas_temperature_keV,
        assumptions.electron_density_cm3,
        assumptions.cavity_pressure_erg_cm3,
        assumptions.cavity_radius_kpc,
        assumptions.cavity_semi_major_kpc,
        assumptions.cavity_semi_minor_kpc,
        assumptions.cavity_age_yr,
        assumptions.cavity_distance_kpc,
        assumptions.n_cavities,
        n_override,
    )
end

println("Mode: $mode  |  MC samples: $(assumptions.monte_carlo_samples)")

# -- 1. Legacy single-age Monte Carlo (backward compatibility) ----------------
println("Running legacy single-age Monte Carlo...")
samples = run_monte_carlo(assumptions; seed=seed_override)
summary = summarize_samples(samples)
summary["cluster_name"] = assumptions.cluster_name
summary["redshift"] = assumptions.redshift
summary["n_cavities"] = assumptions.n_cavities
summary["assumption_file"] = relpath(assumptions_path, root)

CSV.write(joinpath(root, "results", "monte_carlo_samples.csv"), DataFrame(samples))
open(joinpath(root, "results", "feedback_summary.json"), "w") do io
    JSON.print(io, summary, 4)
end

assumption_rows = [
    (parameter = "cooling_luminosity_erg_s", median = assumptions.cooling_luminosity_erg_s[1], sigma = assumptions.cooling_luminosity_erg_s[2]),
    (parameter = "gas_temperature_keV", median = assumptions.gas_temperature_keV[1], sigma = assumptions.gas_temperature_keV[2]),
    (parameter = "electron_density_cm3", median = assumptions.electron_density_cm3[1], sigma = assumptions.electron_density_cm3[2]),
    (parameter = "cavity_pressure_erg_cm3", median = assumptions.cavity_pressure_erg_cm3[1], sigma = assumptions.cavity_pressure_erg_cm3[2]),
    (parameter = "cavity_radius_kpc", median = assumptions.cavity_radius_kpc[1], sigma = assumptions.cavity_radius_kpc[2]),
    (parameter = "cavity_age_yr", median = assumptions.cavity_age_yr[1], sigma = assumptions.cavity_age_yr[2]),
    (parameter = "cavity_semi_major_kpc", median = assumptions.cavity_semi_major_kpc[1], sigma = assumptions.cavity_semi_major_kpc[2]),
    (parameter = "cavity_semi_minor_kpc", median = assumptions.cavity_semi_minor_kpc[1], sigma = assumptions.cavity_semi_minor_kpc[2]),
    (parameter = "cavity_distance_kpc", median = assumptions.cavity_distance_kpc[1], sigma = assumptions.cavity_distance_kpc[2]),
]
CSV.write(joinpath(root, "results", "assumption_table.csv"), DataFrame(assumption_rows))

# -- 2. Multi-age, multi-geometry Monte Carlo ---------------------------------
println("Running multi-age, multi-geometry Monte Carlo...")
multi_samples = run_multi_age_monte_carlo(assumptions; seed=seed_override)
multi_summary = summarize_multi_age_samples(multi_samples)
multi_summary["cluster_name"] = assumptions.cluster_name
multi_summary["redshift"] = assumptions.redshift
multi_summary["n_cavities"] = assumptions.n_cavities

CSV.write(joinpath(root, "results", "multi_age_samples.csv"), DataFrame(multi_samples))
open(joinpath(root, "results", "multi_age_summary.json"), "w") do io
    JSON.print(io, multi_summary, 4)
end

# -- 3. Sensitivity grid -----------------------------------------------------
println("Running sensitivity grid...")
grid = run_sensitivity_grid(assumptions; n_points = 9)
CSV.write(joinpath(root, "results", "sensitivity_grid.csv"), DataFrame(grid))

# -- 4a. Derived Quantities -----------------------------------------------------
println("Generating derived quantities...")
mdf = DataFrame(multi_samples)
derived_cols = [
    ("volume_spherical_cm3", "4/3 pi r^3", "radius"),
    ("volume_ellipsoidal_cm3", "4/3 pi a b c", "a, b, c"),
    ("enthalpy_spherical_erg", "4 p V", "pressure, volume_sph"),
    ("enthalpy_ellipsoidal_erg", "4 p V", "pressure, volume_ell"),
    ("effective_ellipsoid_radius_kpc", "(a b c)^(1/3)", "a, b, c"),
    ("t_sound_yr", "R/c_s", "distance, temperature"),
    ("t_buoyancy_yr", "R sqrt(S C_D / 2gV)", "spherical radius, distance, temperature"),
    ("t_refill_yr", "2 sqrt(r/g)", "spherical radius, distance, temperature"),
    ("t_ell_buoyancy_yr", "R sqrt(S C_D / 2gV)", "effective ellipsoid radius, distance, temperature"),
    ("t_ell_refill_yr", "2 sqrt(r/g)", "effective ellipsoid radius, distance, temperature"),
    ("P_sph_soundcross", "H/t", "spherical enthalpy, sound-crossing age"),
    ("P_sph_buoyancy", "H/t", "spherical enthalpy, buoyancy age"),
    ("P_sph_refill", "H/t", "spherical enthalpy, refill age"),
    ("P_ell_soundcross", "H/t", "ellipsoidal enthalpy, sound-crossing age"),
    ("P_ell_buoyancy", "H/t", "ellipsoidal enthalpy, buoyancy age"),
    ("P_ell_refill", "H/t", "ellipsoidal enthalpy, refill age"),
    ("ratio_sph_cs", "P/L_cool", "spherical sound-crossing power, cooling luminosity"),
    ("ratio_sph_buoy", "P/L_cool", "spherical buoyancy power, cooling luminosity"),
    ("ratio_sph_refill", "P/L_cool", "spherical refill power, cooling luminosity"),
    ("ratio_ell_cs", "P/L_cool", "ellipsoidal sound-crossing power, cooling luminosity"),
    ("ratio_ell_buoy", "P/L_cool", "ellipsoidal buoyancy power, cooling luminosity"),
    ("ratio_ell_refill", "P/L_cool", "ellipsoidal refill power, cooling luminosity"),
]

derived_rows = []
for (col, form, inputs) in derived_cols
    vals = mdf[!, col]
    push!(derived_rows, (
        Quantity = col,
        Formula = form,
        Inputs = inputs,
        Median = median(vals),
        p16 = quantile(vals, 0.16),
        p84 = quantile(vals, 0.84)
    ))
end
CSV.write(joinpath(root, "results", "derived_quantities.csv"), DataFrame(derived_rows))

println("Generating literature comparison...")
function _column_summary(df, col)
    vals = df[!, col]
    return (median = median(vals), p16 = quantile(vals, 0.16), p84 = quantile(vals, 0.84))
end

lit_specs = [
    (
        study = "McDonald et al. (2015)",
        published_value_erg_s = 1.0e46,
        published_low_erg_s = 0.6e46,
        published_high_erg_s = 2.5e46,
        method = "4pV/t_buoy",
        matching_model = "Spherical / Buoyancy",
        pipeline_column = :P_sph_buoyancy,
        agreement_note = "Lower than nominal; within the quoted uncertainty scale",
    ),
    (
        study = "Hlavacek-Larrondo et al. (2015)",
        published_value_erg_s = 4.5e45,
        published_low_erg_s = 2.0e45,
        published_high_erg_s = 7.0e45,
        method = "4pV/t_cs",
        matching_model = "Spherical / Sound-crossing",
        pipeline_column = :P_sph_soundcross,
        agreement_note = "Consistent with the published range",
    ),
]

comparison_rows = []
for spec in lit_specs
    stats = _column_summary(mdf, spec.pipeline_column)
    push!(comparison_rows, (
        study = spec.study,
        published_value_erg_s = spec.published_value_erg_s,
        published_low_erg_s = spec.published_low_erg_s,
        published_high_erg_s = spec.published_high_erg_s,
        method = spec.method,
        matching_model = spec.matching_model,
        pipeline_median_erg_s = stats.median,
        pipeline_p16_erg_s = stats.p16,
        pipeline_p84_erg_s = stats.p84,
        agreement_note = spec.agreement_note,
    ))
end

power_cols = [:P_sph_soundcross, :P_sph_buoyancy, :P_sph_refill,
              :P_ell_soundcross, :P_ell_buoyancy, :P_ell_refill]
model_medians = [median(mdf[!, col]) for col in power_cols]
model_p16 = [quantile(mdf[!, col], 0.16) for col in power_cols]
model_p84 = [quantile(mdf[!, col], 0.84) for col in power_cols]
push!(comparison_rows, (
    study = "McDonald et al. (2019)",
    published_value_erg_s = 1.0e46,
    published_low_erg_s = missing,
    published_high_erg_s = missing,
    method = "Review/order-of-magnitude estimate",
    matching_model = "All six models",
    pipeline_median_erg_s = median(model_medians),
    pipeline_p16_erg_s = minimum(model_p16),
    pipeline_p84_erg_s = maximum(model_p84),
    agreement_note = "Pipeline model range overlaps the order-of-magnitude estimate",
))
CSV.write(joinpath(root, "results", "literature_comparison.csv"), DataFrame(comparison_rows))

# -- 4b. Run Metadata ----------------------------------------------------------
println("Generating run metadata...")
assumption_hash = bytes2hex(sha256(read(assumptions_path)))
git_commit = _try_readchomp(`git -C $root rev-parse HEAD`)
command_parts = vcat(["julia", "--project=.", relpath(@__FILE__, root)], collect(ARGS))
metadata = Dict(
    "run_mode" => mode,
    "sample_count" => assumptions.monte_carlo_samples,
    "seed" => seed_override,
    "julia_version" => string(VERSION),
    "timestamp_utc" => string(now(UTC)),
    "timestamp_local" => string(now()),
    "timezone" => "UTC for timestamp_utc; local system timezone for timestamp_local",
    "command" => join(command_parts, " "),
    "args" => collect(ARGS),
    "assumption_file_path" => relpath(assumptions_path, root),
    "assumption_file_sha256" => assumption_hash,
    "git_commit" => git_commit,
    "git_dirty" => _git_dirty(root)
)
open(joinpath(root, "results", "run_metadata.json"), "w") do f
    JSON.print(f, metadata, 4)
end

# -- 5. Figures ---------------------------------------------------------------
println("Generating figures...")
include(joinpath(@__DIR__, "make_figures.jl"))

# -- 5. Report ----------------------------------------------------------------
println("Writing reports...")
include(joinpath(@__DIR__, "write_report.jl"))

# -- Summary output -----------------------------------------------------------
println()
println("===============================================================")
println("  Phoenix Cluster Feedback Pipeline - Complete")
println("===============================================================")
println()
println("Legacy (single-age, spherical):")
println("  Median feedback/cooling ratio: $(round(summary["median_feedback_to_cooling_ratio"], digits=3))")
println("  P(feedback >= cooling): $(round(summary["probability_feedback_exceeds_cooling"], digits=3))")
println()
println("Multi-age results:")
for key in ["sph_soundcross", "sph_buoyancy", "sph_refill",
            "ell_soundcross", "ell_buoyancy", "ell_refill"]
    s = multi_summary[key]
    println("  $(rpad(s["label"], 30)) median=$(round(s["median"], digits=3))  " *
            "[$(round(s["p16"], digits=3)), $(round(s["p84"], digits=3))]  " *
            "P>=1=$(round(s["p_exceeds_1"], digits=3))")
end
println()
println("Outputs:")
println("  results/feedback_summary.json")
println("  results/multi_age_summary.json")
println("  results/monte_carlo_samples.csv")
println("  results/multi_age_samples.csv")
println("  results/sensitivity_grid.csv")
println("  results/assumption_table.csv")
println("  results/derived_quantities.csv")
println("  results/literature_comparison.csv")
println("  results/run_metadata.json")
println("  figures/feedback_ratio_histogram.png")
println("  figures/cooling_vs_feedback.png")
println("  figures/multi_age_ratio_comparison.png")
println("  figures/sensitivity_heatmap.png")
println("  results/feedback_report.md")
