#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "PhoenixFeedback.jl"))
using .PhoenixFeedback

assumptions = load_assumptions(joinpath(@__DIR__, "..", "data", "assumptions.json"))
summary = summarize_samples(run_monte_carlo(assumptions))
for (key, value) in summary
    println("$key = $value")
end

