#!/usr/bin/env julia
# Generate LaTeX input parameter table from data/assumptions.json
# Ensures the manuscript table is always consistent with the pipeline inputs.
#
# Usage: julia scripts/generate_input_table.jl
# Output: writes LaTeX fragment to stdout

using JSON

root = normpath(joinpath(@__DIR__, ".."))
raw = JSON.parsefile(joinpath(root, "data", "assumptions.json"))

function pct(block)
    sigma = Float64(block["sigma_fraction"])
    return string("\\pm ", round(Int, sigma * 100), "\\%")
end

function valfmt(val)
    if val >= 1e3
        return raw"$" * string(val) * raw"$"
    else
        return raw"$" * string(val) * raw"$"
    end
end

entries = [
    (raw"$\Lcool$ ($\erg\;\s^{-1}$)",
     raw"$8.2 \times 10^{45}$", pct(raw["cooling_luminosity_erg_s"]),
     raw"$r < 100\;\kpc$, 0.7--7.0~keV",
     raw"\citet{McDonald2012, McDonald2019}"),

    (raw"$kT$ ($\keV$)",
     raw"5.9", pct(raw["gas_temperature_keV"]),
     raw"10--30~kpc",
     raw"\citet{McDonald2019}"),

    (raw"$n_e$ ($\cm^{-3}$)",
     raw"0.12", pct(raw["electron_density_cm3"]),
     raw"$r < 10\;\kpc$",
     raw"\citet{McDonald2019}"),

    (raw"$\pcav$ ($\erg\;\cm^{-3}$)",
     raw"$1.5 \times 10^{-9}$", pct(raw["cavity_pressure_erg_cm3"]),
     raw"cavity location",
     raw"\citet{McDonald2015}"),

    (raw"$\rcav$ ($\kpc$)",
     raw"12", pct(raw["cavity_radius_kpc"]),
     raw"projected",
     raw"\citet{McDonald2015}"),

    (raw"$a$ ($\kpc$)",
     raw"14", pct(raw["cavity_geometry"]["semi_major_kpc"]),
     raw"major axis",
     raw"\citet{Hlavacek-Larrondo2015}"),

    (raw"$b$ ($\kpc$)",
     raw"9", pct(raw["cavity_geometry"]["semi_minor_kpc"]),
     raw"minor axis",
     raw"\citet{Hlavacek-Larrondo2015}"),

    (raw"$R$ ($\kpc$)",
     raw"25", pct(raw["cavity_distance_kpc"]),
     raw"projected distance",
     raw"\citet{McDonald2015}"),
]

for (param, value, unc, aperture, source) in entries
    println("$(param) & $(value) & $(unc) & $(aperture) & $(source) \\\\")
end