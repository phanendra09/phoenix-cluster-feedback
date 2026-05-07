using CSV
using DataFrames
using JSON
using LaTeXStrings
using Plots
using Statistics

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
# ApJ full-width figure: 7.5 inches at 300 dpi = 2250 px
# ApJ single-column: 3.5 inches at 300 dpi = 1050 px
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

# -- 1. Feedback ratio histogram (legacy, single-age spherical) ----------------
samples_path = joinpath(root, "results", "monte_carlo_samples.csv")
if isfile(samples_path)
    df = CSV.read(samples_path, DataFrame)

    p1 = histogram(
        df.feedback_to_cooling_ratio;
        bins = 60,
        xlabel = L"P_{\rm cav} / L_{\rm cool}",
        ylabel = "Monte Carlo samples",
        title = "Phoenix Cluster Feedback-Cooling Ratio",
        legend = false,
        size = (2250, 1500),
        color = CB_SKYBLUE,
        alpha = 0.85,
    )
    vline!(p1, [1.0], color = CB_VERMIL, linestyle = :dash, linewidth = 3,
           label = "Energetic balance")
    savefig(p1, joinpath(root, "figures", "feedback_ratio_histogram.png"))
    savefig(p1, joinpath(root, "figures", "feedback_ratio_histogram.pdf"))

    p4 = scatter(
        df.cooling_luminosity_erg_s,
        df.cavity_power_erg_s;
        xscale = :log10,
        yscale = :log10,
        xlabel = L"Cooling luminosity $L_{\rm cool}$ (erg s$^{-1}$)",
        ylabel = L"Cavity power $P_{\rm cav}$ (erg s$^{-1}$)",
        title = "Cooling Losses vs. AGN Cavity Power",
        markersize = 2,
        alpha = 0.30,
        markerstrokewidth = 0,
        legend = false,
        size = (2250, 1500),
        color = CB_BLUE,
    )
    savefig(p4, joinpath(root, "figures", "cooling_vs_feedback.png"))
    savefig(p4, joinpath(root, "figures", "cooling_vs_feedback.pdf"))
end

# -- 2. Multi-age comparison plot ---------------------------------------------
multi_path = joinpath(root, "results", "multi_age_samples.csv")
if isfile(multi_path)
    mdf = CSV.read(multi_path, DataFrame)

    ratio_cols = [
        (:ratio_sph_cs,     "Sph / Sound",      CB_ORANGE),
        (:ratio_sph_buoy,   "Sph / Buoyancy",   CB_SKYBLUE),
        (:ratio_sph_refill, "Sph / Refill",      CB_GREEN),
        (:ratio_ell_cs,     "Ell / Sound",       CB_VERMIL),
        (:ratio_ell_buoy,   "Ell / Buoyancy",    CB_BLUE),
        (:ratio_ell_refill, "Ell / Refill",       CB_PINK),
    ]

    p2 = plot(; xlabel = L"P_{\rm cav} / L_{\rm cool}",
                  ylabel = "Probability density",
                  title = "Feedback Ratio by Age Definition & Geometry",
                  legend = :topright,
                  size = (2250, 1500))

    for (col, lbl, clr) in ratio_cols
        vals = mdf[!, col]
        clipped = filter(v -> v < quantile(vals, 0.99), vals)
        histogram!(p2, clipped; bins = 60, alpha = 0.40, label = lbl,
                   normalize = :pdf, color = clr)
    end
    vline!(p2, [1.0], color = CB_BLACK, linestyle = :dash, linewidth = 3,
           label = "Energetic balance")
    savefig(p2, joinpath(root, "figures", "multi_age_ratio_comparison.png"))
    savefig(p2, joinpath(root, "figures", "multi_age_ratio_comparison.pdf"))
end

# -- 3. Sensitivity analysis ---------------------------------------------------
sens_path = joinpath(root, "results", "sensitivity_grid.csv")
if isfile(sens_path)
    sdf = CSV.read(sens_path, DataFrame)

    param_colors = Dict(
        "cavity_radius_kpc"       => CB_ORANGE,
        "cavity_pressure_erg_cm3" => CB_SKYBLUE,
        "cavity_age_yr"           => CB_GREEN,
        "cooling_luminosity_erg_s" => CB_VERMIL,
    )

    param_labels = Dict(
        "cavity_radius_kpc"        => L"Cavity radius $r_{\rm cav}$",
        "cavity_pressure_erg_cm3"  => L"Cavity pressure $p_{\rm cav}$",
        "cavity_age_yr"            => L"Cavity age $t_{\rm cav}$",
        "cooling_luminosity_erg_s" => L"Cooling luminosity $L_{\rm cool}$",
    )

    params = unique(sdf.varied_parameter)
    p3 = plot(; xlabel = "Multiplicative factor",
                  ylabel = L"P_{\rm cav} / L_{\rm cool}",
                  title = "Parameter Sensitivity Analysis",
                  legend = :topleft,
                  size = (2250, 1500))

    for param in params
        sub = filter(r -> r.varied_parameter == param, sdf)
        clr = get(param_colors, param, CB_BLACK)
        lbl = get(param_labels, param, replace(param, "_" => " "))
        plot!(p3, sub.factor, sub.feedback_ratio;
              label = lbl, color = clr,
              linewidth = 2.5, marker = :circle, markersize = 5)
    end
    hline!(p3, [1.0], color = CB_BLACK, linestyle = :dash, linewidth = 3,
           label = "Energetic balance")
    savefig(p3, joinpath(root, "figures", "sensitivity_heatmap.png"))
    savefig(p3, joinpath(root, "figures", "sensitivity_heatmap.pdf"))
end

# -- 5. Cluster comparison bar chart -------------------------------------------
comp_path = joinpath(root, "data", "cluster_comparison.json")
if isfile(comp_path)
    comp_data = JSON.parsefile(comp_path)
    clusters = comp_data["clusters"]

    cluster_names = [c["name"] for c in clusters]
    ratios_buoyancy = [c["Pcav_Lcool_buoyancy"] for c in clusters]
    ratios_sound    = [c["Pcav_Lcool_sound"] for c in clusters]

    n_clusters = length(cluster_names)
    x_buoy = collect(1:n_clusters) .- 0.15
    x_snd  = collect(1:n_clusters) .+ 0.15

    p5 = bar(x_buoy, ratios_buoyancy;
             bar_width = 0.28,
             label = "Buoyancy age",
             xlabel = "Cluster",
             ylabel = L"P_{\rm cav} / L_{\rm cool}",
             title = "Feedback-Cooling Ratio: Phoenix in Context",
             legend = :topright,
             size = (2250, 1500),
             color = CB_SKYBLUE,
             alpha = 0.85,
             linewidth = 0,
             xticks = (collect(1:n_clusters), cluster_names),
    )
    bar!(p5, x_snd, ratios_sound;
         bar_width = 0.28,
         label = "Sound-crossing age",
         color = CB_ORANGE,
         alpha = 0.85,
         linewidth = 0,
    )
    hline!(p5, [1.0], color = CB_BLACK, linestyle = :dash, linewidth = 3,
           label = "Energetic balance")
    savefig(p5, joinpath(root, "figures", "cluster_comparison.png"))
    savefig(p5, joinpath(root, "figures", "cluster_comparison.pdf"))
end

# -- 6. Projection sensitivity figure ------------------------------------------
proj_path = joinpath(root, "results", "projection_sensitivity.csv")
if isfile(proj_path)
    pdf = CSV.read(proj_path, DataFrame)
    proj_models = ["Spherical / Sound-crossing", "Spherical / Buoyancy",
                   "Ellipsoidal / Sound-crossing", "Ellipsoidal / Buoyancy"]
    model_colors = Dict(
        "Spherical / Sound-crossing"   => CB_ORANGE,
        "Spherical / Buoyancy"         => CB_SKYBLUE,
        "Ellipsoidal / Sound-crossing" => CB_VERMIL,
        "Ellipsoidal / Buoyancy"       => CB_BLUE,
    )

    p6 = plot(; xlabel = "Projection factor  R_true / R_proj",
                  ylabel = L"P_{\rm cav} / L_{\rm cool}",
                  title = "Projection Sensitivity",
                  legend = :topright,
                  size = (2250, 1500))

    for mdl in proj_models
        sub = filter(r -> r.model == mdl, pdf)
        clr = get(model_colors, mdl, CB_BLACK)
        plot!(p6, sub.proj_factor, sub.median_ratio;
              label = mdl, color = clr, linewidth = 2.5,
              marker = :circle, markersize = 5,
              ribbon = (sub.median_ratio .- sub.p16,
                        sub.p84 .- sub.median_ratio),
              fillalpha = 0.15)
    end
    hline!(p6, [1.0], color = CB_BLACK, linestyle = :dash, linewidth = 3,
           label = "Energetic balance")
    savefig(p6, joinpath(root, "figures", "projection_sensitivity.png"))
    savefig(p6, joinpath(root, "figures", "projection_sensitivity.pdf"))
end

println("Figures saved: 6 figures (PNG + PDF each)")
