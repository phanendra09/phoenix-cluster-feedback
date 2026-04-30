#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "PhoenixFeedback.jl"))
using .PhoenixFeedback

assumptions = load_assumptions(joinpath(@__DIR__, "..", "data", "assumptions.json"))
println("Configured cooling luminosity: $(assumptions.cooling_luminosity_erg_s[1]) erg/s")
println("Replace this script with aperture-specific X-ray luminosity extraction in v2.")
