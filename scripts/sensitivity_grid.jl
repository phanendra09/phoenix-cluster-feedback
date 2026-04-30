#!/usr/bin/env julia
# Sensitivity grid: vary each key parameter independently and record
# how the feedback-to-cooling ratio changes.

using CSV
using DataFrames

include(joinpath(@__DIR__, "..", "src", "PhoenixFeedback.jl"))
using .PhoenixFeedback

root = normpath(joinpath(@__DIR__, ".."))
mkpath(joinpath(root, "results"))

assumptions_path = joinpath(root, "data", "assumptions.json")
assumptions = load_assumptions(assumptions_path)

grid = run_sensitivity_grid(assumptions; n_points = 9)
CSV.write(joinpath(root, "results", "sensitivity_grid.csv"), DataFrame(grid))

println("Sensitivity grid written: results/sensitivity_grid.csv  ($(length(grid)) rows)")
