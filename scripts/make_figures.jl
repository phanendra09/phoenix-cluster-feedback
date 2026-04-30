using CSV
using DataFrames
using Plots
using Statistics

# Apply publication-quality defaults
default(
    titlefontsize = 14,
    guidefontsize = 12,
    tickfontsize = 10,
    legendfontsize = 10,
    linewidth = 2,
    grid = false,
    dpi = 300,
    framestyle = :box
)

root = normpath(joinpath(@__DIR__, ".."))

# -- 1. Feedback ratio histogram (legacy) -------------------------------------
samples_path = joinpath(root, "results", "monte_carlo_samples.csv")
if isfile(samples_path)
    df = CSV.read(samples_path, DataFrame)

    histogram(
        df.feedback_to_cooling_ratio;
        bins = 60,
        xlabel = "AGN feedback power / cooling luminosity",
        ylabel = "Monte Carlo samples",
        title = "Phoenix Cluster feedback-cooling ratio",
        legend = false,
    )
    vline!([1.0], color = :red, linestyle = :dash, linewidth = 2.5)
    savefig(joinpath(root, "figures", "feedback_ratio_histogram.png"))
    savefig(joinpath(root, "figures", "feedback_ratio_histogram.pdf"))

    scatter(
        df.cooling_luminosity_erg_s,
        df.cavity_power_erg_s;
        xscale = :log10,
        yscale = :log10,
        xlabel = "Cooling luminosity (erg s^-1)",
        ylabel = "Cavity power (erg s^-1)",
        title = "Cooling losses versus AGN cavity power",
        markersize = 2,
        alpha = 0.35,
        markerstrokewidth = 0,
        legend = false,
    )
    savefig(joinpath(root, "figures", "cooling_vs_feedback.png"))
    savefig(joinpath(root, "figures", "cooling_vs_feedback.pdf"))
end

# -- 2. Multi-age comparison plot ---------------------------------------------
multi_path = joinpath(root, "results", "multi_age_samples.csv")
if isfile(multi_path)
    mdf = CSV.read(multi_path, DataFrame)

    ratio_cols = [
        (:ratio_sph_cs,     "Sph/Sound"),
        (:ratio_sph_buoy,   "Sph/Buoyancy"),
        (:ratio_sph_refill, "Sph/Refill"),
        (:ratio_ell_cs,     "Ell/Sound"),
        (:ratio_ell_buoy,   "Ell/Buoyancy"),
        (:ratio_ell_refill, "Ell/Refill"),
    ]

    p = plot(; xlabel = "Feedback / Cooling ratio",
               ylabel = "Density",
               title = "Feedback ratio by age definition & geometry",
               legend = :topright,
               size = (900, 500))

    for (col, lbl) in ratio_cols
        vals = mdf[!, col]
        # clip extreme outliers for clean plot
        clipped = filter(v -> v < quantile(vals, 0.99), vals)
        histogram!(p, clipped; bins = 60, alpha = 0.35, label = lbl, normalize = :pdf)
    end
    vline!(p, [1.0], color = :black, linestyle = :dash, linewidth = 2.5, label = "Balance")
    savefig(p, joinpath(root, "figures", "multi_age_ratio_comparison.png"))
    savefig(p, joinpath(root, "figures", "multi_age_ratio_comparison.pdf"))
end

# -- 3. Sensitivity heatmap ---------------------------------------------------
sens_path = joinpath(root, "results", "sensitivity_grid.csv")
if isfile(sens_path)
    sdf = CSV.read(sens_path, DataFrame)

    params = unique(sdf.varied_parameter)
    p_sens = plot(; xlabel = "Multiplicative factor",
                    ylabel = "Feedback / Cooling ratio",
                    title = "Parameter sensitivity analysis",
                    legend = :topleft,
                    size = (800, 500))

    for param in params
        sub = filter(r -> r.varied_parameter == param, sdf)
        plot!(p_sens, sub.factor, sub.feedback_ratio;
              label = replace(param, "_" => " "),
              linewidth = 2, marker = :circle, markersize = 4)
    end
    hline!(p_sens, [1.0], color = :black, linestyle = :dash, linewidth = 2.5, label = "Balance")
    savefig(p_sens, joinpath(root, "figures", "sensitivity_heatmap.png"))
    savefig(p_sens, joinpath(root, "figures", "sensitivity_heatmap.pdf"))
end
