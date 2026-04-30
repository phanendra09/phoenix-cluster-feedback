#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "PhoenixFeedback.jl"))
using .PhoenixFeedback

assumptions = load_assumptions(joinpath(@__DIR__, "..", "data", "assumptions.json"))
t = cooling_time_seconds(assumptions.gas_temperature_keV[1], assumptions.electron_density_cm3[1])
println("Cooling time: $(t / YEAR_TO_SECONDS) yr")

