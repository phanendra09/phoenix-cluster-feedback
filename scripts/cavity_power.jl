#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "PhoenixFeedback.jl"))
using .PhoenixFeedback

assumptions = load_assumptions(joinpath(@__DIR__, "..", "data", "assumptions.json"))
p = cavity_power_erg_s(
    assumptions.cavity_pressure_erg_cm3[1],
    assumptions.cavity_radius_kpc[1],
    assumptions.cavity_age_yr[1];
    n_cavities = assumptions.n_cavities,
)
println("Cavity power: $(p) erg/s")

