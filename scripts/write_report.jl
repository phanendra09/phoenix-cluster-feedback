using JSON

root = normpath(joinpath(@__DIR__, ".."))
report_path = joinpath(root, "results", "feedback_report.md")

# -- Load summaries -----------------------------------------------------------
legacy_path = joinpath(root, "results", "feedback_summary.json")
multi_path  = joinpath(root, "results", "multi_age_summary.json")

if !isfile(legacy_path)
    @warn "feedback_summary.json not found; skipping report."
    return
end

legacy = JSON.parsefile(legacy_path)

open(report_path, "w") do io
    println(io, "# Phoenix Cluster AGN Feedback Report")
    println(io)
    println(io, "**Cluster:** $(legacy["cluster_name"])  ")
    println(io, "**Redshift:** $(legacy["redshift"])  ")
    println(io, "**Monte Carlo samples:** $(legacy["n_samples"])  ")
    println(io)

    # -- Legacy result ----------------------------------------------------
    println(io, "## 1. Single-Age Spherical Analysis (Legacy)")
    println(io)
    ratio = legacy["median_feedback_to_cooling_ratio"]
    p16   = legacy["ratio_p16"]
    p84   = legacy["ratio_p84"]
    prob  = legacy["probability_feedback_exceeds_cooling"]
    println(io, "| Metric | Value |")
    println(io, "|--------|-------|")
    println(io, "| Median feedback/cooling ratio | $(round(ratio, sigdigits=4)) |")
    println(io, "| 16th-84th percentile interval | [$(round(p16, sigdigits=4)), $(round(p84, sigdigits=4))] |")
    println(io, "| P(feedback >= cooling) | $(round(prob, sigdigits=4)) |")
    println(io)

    # -- Multi-age result -------------------------------------------------
    if isfile(multi_path)
        multi = JSON.parsefile(multi_path)
        println(io, "## 2. Multi-Age, Multi-Geometry Analysis")
        println(io)
        println(io, "| Model | Median ratio | 16th pctl | 84th pctl | P(>=1) |")
        println(io, "|-------|-------------|-----------|-----------|-------|")

        for key in ["sph_soundcross", "sph_buoyancy", "sph_refill",
                     "ell_soundcross", "ell_buoyancy", "ell_refill"]
            s = multi[key]
            println(io, "| $(s["label"]) | $(round(s["median"], sigdigits=4)) " *
                        "| $(round(s["p16"], sigdigits=4)) " *
                        "| $(round(s["p84"], sigdigits=4)) " *
                        "| $(round(s["p_exceeds_1"], sigdigits=3)) |")
        end
        println(io)
        println(io, "Median cooling time: $(round(multi["median_cooling_time_yr"] / 1e6, sigdigits=4)) Myr")
        println(io)
    end

    # -- Interpretation ---------------------------------------------------
    println(io, "## 3. Interpretation")
    println(io)
    println(io, "Under the adopted literature-anchored assumptions, Phoenix Cluster AGN cavity ")
    println(io, "power is energetically comparable to the cooling luminosity, but most ")
    println(io, "age/geometry combinations remain below unity at the median. The conclusion ")
    println(io, "depends strongly on the adopted cavity age estimator and geometry model. ")
    println(io, "The buoyancy-time estimator yields the highest feedback ratios in this ")
    println(io, "implementation, while refill ages produce the lowest median ratios. ")
    println(io, "Sound-crossing results fall between these regimes. Ellipsoidal geometry ")
    println(io, "reduces cavity volume relative to the ")
    println(io, "spherical assumption when the projected semi-minor axis is smaller than ")
    println(io, "the semi-major axis.")
    println(io)
    println(io, "These results should be treated as energetic consistency tests under ")
    println(io, "transparent assumptions, not as definitive cavity analyses.")
end

println("Report written: results/feedback_report.md")
