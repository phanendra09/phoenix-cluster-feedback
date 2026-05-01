module PhoenixFeedback

using JSON
using ProgressMeter
using Random
using Statistics

export ClusterAssumptions,
       load_assumptions,
       cooling_time_seconds,
       cavity_enthalpy_erg,
       cavity_power_erg_s,
       feedback_ratio,
       run_monte_carlo,
       summarize_samples,
       spherical_volume_cm3,
       ellipsoidal_volume_cm3,
       cavity_enthalpy_from_volume,
       cavity_power_from_enthalpy,
       sound_crossing_time_s,
       buoyancy_time_s,
       refill_time_s,
       run_multi_age_monte_carlo,
       summarize_multi_age_samples,
       run_sensitivity_grid

# -- Physical constants -------------------------------------------------------
const KEV_TO_ERG = 1.602176634e-9
const KPC_TO_CM = 3.0856775814913673e21
const YEAR_TO_SECONDS = 31_557_600.0
const PROTON_MASS_G = 1.67262192e-24

# -- Include sub-modules ------------------------------------------------------
include("geometry.jl")
include("cavity_age.jl")

# -- Data structures ----------------------------------------------------------

"""
Represents the full set of observational assumptions for the Phoenix Cluster
feedback analysis. Each physical quantity is a (median, sigma) tuple.
"""
struct ClusterAssumptions
    cluster_name::String
    redshift::Float64
    # Cooling
    cooling_luminosity_erg_s::Tuple{Float64, Float64}
    gas_temperature_keV::Tuple{Float64, Float64}
    electron_density_cm3::Tuple{Float64, Float64}
    # Cavity geometry
    cavity_pressure_erg_cm3::Tuple{Float64, Float64}
    # Spherical model
    cavity_radius_kpc::Tuple{Float64, Float64}
    # Ellipsoidal model (semi-axes a, b; c = sqrt(a*b) by default)
    cavity_semi_major_kpc::Tuple{Float64, Float64}
    cavity_semi_minor_kpc::Tuple{Float64, Float64}
    # Age inputs
    cavity_age_yr::Tuple{Float64, Float64}           # legacy single age
    cavity_distance_kpc::Tuple{Float64, Float64}     # projected distance from center
    n_cavities::Int
    monte_carlo_samples::Int
end

# -- JSON loading -------------------------------------------------------------

function _value_sigma(block)
    median_value = Float64(block["median"])
    sigma_fraction = Float64(block["sigma_fraction"])
    return (median_value, abs(median_value * sigma_fraction))
end

function load_assumptions(path::AbstractString)::ClusterAssumptions
    raw = JSON.parsefile(path)

    # Ellipsoidal geometry (fallback to spherical radius if not provided)
    geom = get(raw, "cavity_geometry", nothing)
    if geom !== nothing
        semi_major = _value_sigma(geom["semi_major_kpc"])
        semi_minor = _value_sigma(geom["semi_minor_kpc"])
    else
        # Fall back: treat spherical radius as both axes
        r = _value_sigma(raw["cavity_radius_kpc"])
        semi_major = r
        semi_minor = r
    end

    # Cavity distance from center for dynamical age calculations
    dist_block = get(raw, "cavity_distance_kpc", nothing)
    if dist_block !== nothing
        cavity_dist = _value_sigma(dist_block)
    else
        # Default: use spherical radius as distance proxy
        cavity_dist = _value_sigma(raw["cavity_radius_kpc"])
    end

    return ClusterAssumptions(
        raw["cluster_name"],
        Float64(raw["redshift"]),
        _value_sigma(raw["cooling_luminosity_erg_s"]),
        _value_sigma(raw["gas_temperature_keV"]),
        _value_sigma(raw["electron_density_cm3"]),
        _value_sigma(raw["cavity_pressure_erg_cm3"]),
        _value_sigma(raw["cavity_radius_kpc"]),
        semi_major,
        semi_minor,
        _value_sigma(raw["cavity_age_yr"]),
        cavity_dist,
        Int(raw["n_cavities"]),
        Int(raw["monte_carlo_samples"]),
    )
end

# -- Core physics (backward-compatible) ---------------------------------------

function cooling_time_seconds(temperature_keV::Real, electron_density_cm3::Real,
                              cooling_function_erg_cm3_s::Real = 2.0e-23)
    thermal_energy_density = 3.0 * electron_density_cm3 * temperature_keV * KEV_TO_ERG
    emissivity = electron_density_cm3^2 * cooling_function_erg_cm3_s
    return thermal_energy_density / emissivity
end

function cavity_enthalpy_erg(pressure_erg_cm3::Real, radius_kpc::Real; gamma = 4.0 / 3.0)
    radius_cm = radius_kpc * KPC_TO_CM
    volume_cm3 = 4.0 / 3.0 * pi * radius_cm^3
    return gamma / (gamma - 1.0) * pressure_erg_cm3 * volume_cm3
end

function cavity_power_erg_s(pressure_erg_cm3::Real, radius_kpc::Real,
                            age_yr::Real; n_cavities::Integer = 2)
    enthalpy = cavity_enthalpy_erg(pressure_erg_cm3, radius_kpc)
    age_seconds = age_yr * YEAR_TO_SECONDS
    return n_cavities * enthalpy / age_seconds
end

feedback_ratio(power_erg_s::Real, cooling_luminosity_erg_s::Real) =
    power_erg_s / cooling_luminosity_erg_s

# -- Sampling utility ---------------------------------------------------------

function _positive_normal(rng::AbstractRNG, mu::Float64, sigma::Float64)
    value = mu + sigma * randn(rng)
    while value <= 0.0
        value = mu + sigma * randn(rng)
    end
    return value
end

# -- Legacy Monte Carlo (backward-compatible) ---------------------------------

function run_monte_carlo(assumptions::ClusterAssumptions; seed::Integer = 42)
    rng = MersenneTwister(seed)
    n = assumptions.monte_carlo_samples
    rows = Vector{NamedTuple}(undef, n)

    for i in 1:n
        lcool = _positive_normal(rng, assumptions.cooling_luminosity_erg_s...)
        temp = _positive_normal(rng, assumptions.gas_temperature_keV...)
        ne = _positive_normal(rng, assumptions.electron_density_cm3...)
        pressure = _positive_normal(rng, assumptions.cavity_pressure_erg_cm3...)
        radius = _positive_normal(rng, assumptions.cavity_radius_kpc...)
        age = _positive_normal(rng, assumptions.cavity_age_yr...)

        tcool_s = cooling_time_seconds(temp, ne)
        pcav = cavity_power_erg_s(pressure, radius, age; n_cavities = assumptions.n_cavities)
        rows[i] = (
            cooling_luminosity_erg_s = lcool,
            temperature_keV = temp,
            electron_density_cm3 = ne,
            cooling_time_yr = tcool_s / YEAR_TO_SECONDS,
            cavity_pressure_erg_cm3 = pressure,
            cavity_radius_kpc = radius,
            cavity_age_yr = age,
            cavity_power_erg_s = pcav,
            feedback_to_cooling_ratio = feedback_ratio(pcav, lcool),
        )
    end

    return rows
end

# -- Multi-age, multi-geometry Monte Carlo ------------------------------------

"""
    run_multi_age_monte_carlo(assumptions; seed=42)

Run Monte Carlo sampling that computes feedback ratios for:
  - 3 age definitions: sound-crossing, buoyancy, refill
  - 2 geometry models: spherical, ellipsoidal
Returns a vector of NamedTuples with all combinations.
"""
function run_multi_age_monte_carlo(assumptions::ClusterAssumptions; seed::Integer = 42)
    rng = MersenneTwister(seed)
    n = assumptions.monte_carlo_samples
    rows = Vector{NamedTuple}(undef, n)

    prog = Progress(n; desc = "Multi-age MC: ", showspeed = true)
    for i in 1:n
        next!(prog)
        lcool    = _positive_normal(rng, assumptions.cooling_luminosity_erg_s...)
        temp     = _positive_normal(rng, assumptions.gas_temperature_keV...)
        ne       = _positive_normal(rng, assumptions.electron_density_cm3...)
        pressure = _positive_normal(rng, assumptions.cavity_pressure_erg_cm3...)
        radius   = _positive_normal(rng, assumptions.cavity_radius_kpc...)
        a_kpc    = _positive_normal(rng, assumptions.cavity_semi_major_kpc...)
        b_kpc    = _positive_normal(rng, assumptions.cavity_semi_minor_kpc...)
        dist     = _positive_normal(rng, assumptions.cavity_distance_kpc...)

        # Line-of-sight axis: geometric mean convention
        c_kpc = sqrt(a_kpc * b_kpc)

        # Volumes
        vol_sphere = spherical_volume_cm3(radius)
        vol_ellip  = ellipsoidal_volume_cm3(a_kpc, b_kpc, c_kpc)

        # Enthalpies (4pV for gamma = 4/3)
        H_sphere = cavity_enthalpy_from_volume(pressure, vol_sphere)
        H_ellip  = cavity_enthalpy_from_volume(pressure, vol_ellip)

        # Effective radius for ellipsoidal dynamics
        r_eff_ellip = cbrt(a_kpc * b_kpc * c_kpc)

        # Three dynamical ages (in seconds)
        # Sound-crossing age depends only on distance and temperature, not cavity
        # radius — so it is the same for both spherical and ellipsoidal geometries.
        t_cs     = sound_crossing_time_s(dist, temp)
        
        # Spherical ages
        t_sph_buoy   = buoyancy_time_s(dist, radius, temp)
        t_sph_refill = refill_time_s(radius, dist, temp)

        # Ellipsoidal ages
        t_ell_buoy   = buoyancy_time_s(dist, r_eff_ellip, temp)
        t_ell_refill = refill_time_s(r_eff_ellip, dist, temp)

        ncav = assumptions.n_cavities

        # Powers for spherical geometry
        P_sph_cs     = cavity_power_from_enthalpy(H_sphere, t_cs;         n_cavities = ncav)
        P_sph_buoy   = cavity_power_from_enthalpy(H_sphere, t_sph_buoy;   n_cavities = ncav)
        P_sph_refill = cavity_power_from_enthalpy(H_sphere, t_sph_refill; n_cavities = ncav)

        # Powers for ellipsoidal geometry
        P_ell_cs     = cavity_power_from_enthalpy(H_ellip, t_cs;         n_cavities = ncav)
        P_ell_buoy   = cavity_power_from_enthalpy(H_ellip, t_ell_buoy;   n_cavities = ncav)
        P_ell_refill = cavity_power_from_enthalpy(H_ellip, t_ell_refill; n_cavities = ncav)

        tcool_s = cooling_time_seconds(temp, ne)

        rows[i] = (
            cooling_luminosity_erg_s  = lcool,
            temperature_keV           = temp,
            electron_density_cm3      = ne,
            cooling_time_yr           = tcool_s / YEAR_TO_SECONDS,
            cavity_pressure_erg_cm3   = pressure,
            cavity_radius_kpc         = radius,
            volume_spherical_cm3      = vol_sphere,
            volume_ellipsoidal_cm3    = vol_ellip,
            enthalpy_spherical_erg    = H_sphere,
            enthalpy_ellipsoidal_erg  = H_ellip,
            semi_major_kpc            = a_kpc,
            semi_minor_kpc            = b_kpc,
            los_axis_kpc              = c_kpc,
            effective_ellipsoid_radius_kpc = r_eff_ellip,
            cavity_distance_kpc       = dist,
            t_sound_yr                = t_cs / YEAR_TO_SECONDS,
            t_buoyancy_yr             = t_sph_buoy / YEAR_TO_SECONDS,
            t_refill_yr               = t_sph_refill / YEAR_TO_SECONDS,
            t_ell_buoyancy_yr         = t_ell_buoy / YEAR_TO_SECONDS,
            t_ell_refill_yr           = t_ell_refill / YEAR_TO_SECONDS,
            # Spherical powers
            P_sph_soundcross          = P_sph_cs,
            P_sph_buoyancy            = P_sph_buoy,
            P_sph_refill              = P_sph_refill,
            # Ellipsoidal powers
            P_ell_soundcross          = P_ell_cs,
            P_ell_buoyancy            = P_ell_buoy,
            P_ell_refill              = P_ell_refill,
            # Feedback ratios - spherical
            ratio_sph_cs              = feedback_ratio(P_sph_cs, lcool),
            ratio_sph_buoy            = feedback_ratio(P_sph_buoy, lcool),
            ratio_sph_refill          = feedback_ratio(P_sph_refill, lcool),
            # Feedback ratios - ellipsoidal
            ratio_ell_cs              = feedback_ratio(P_ell_cs, lcool),
            ratio_ell_buoy            = feedback_ratio(P_ell_buoy, lcool),
            ratio_ell_refill          = feedback_ratio(P_ell_refill, lcool),
        )
    end

    return rows
end

# -- Summarize multi-age samples ----------------------------------------------

function _ratio_summary(ratios, label)
    return Dict{String, Any}(
        "label"       => label,
        "median"      => median(ratios),
        "p16"         => quantile(ratios, 0.16),
        "p84"         => quantile(ratios, 0.84),
        "p_exceeds_1" => count(r >= 1.0 for r in ratios) / length(ratios),
    )
end

function summarize_multi_age_samples(samples)
    return Dict(
        "n_samples"         => length(samples),
        "sph_soundcross"    => _ratio_summary([r.ratio_sph_cs     for r in samples], "Spherical / Sound-crossing"),
        "sph_buoyancy"      => _ratio_summary([r.ratio_sph_buoy   for r in samples], "Spherical / Buoyancy"),
        "sph_refill"        => _ratio_summary([r.ratio_sph_refill for r in samples], "Spherical / Refill"),
        "ell_soundcross"    => _ratio_summary([r.ratio_ell_cs     for r in samples], "Ellipsoidal / Sound-crossing"),
        "ell_buoyancy"      => _ratio_summary([r.ratio_ell_buoy   for r in samples], "Ellipsoidal / Buoyancy"),
        "ell_refill"        => _ratio_summary([r.ratio_ell_refill for r in samples], "Ellipsoidal / Refill"),
        "median_cooling_time_yr" => median([r.cooling_time_yr for r in samples]),
    )
end

# -- Legacy summarize (backward-compatible) -----------------------------------

function summarize_samples(samples::Vector{NamedTuple})
    ratios = [row.feedback_to_cooling_ratio for row in samples]
    powers = [row.cavity_power_erg_s for row in samples]
    tcool = [row.cooling_time_yr for row in samples]
    return Dict{String, Any}(
        "n_samples" => length(samples),
        "median_feedback_to_cooling_ratio" => median(ratios),
        "ratio_p16" => quantile(ratios, 0.16),
        "ratio_p84" => quantile(ratios, 0.84),
        "probability_feedback_exceeds_cooling" => count(r >= 1.0 for r in ratios) / length(ratios),
        "median_cavity_power_erg_s" => median(powers),
        "median_cooling_time_yr" => median(tcool),
    )
end

# -- Sensitivity grid ---------------------------------------------------------

"""
    run_sensitivity_grid(assumptions; n_points=5, seed=42)

Vary each key parameter by ±50% around its median while holding others fixed.
Returns a vector of NamedTuples for creating sensitivity tables.
"""
function run_sensitivity_grid(assumptions::ClusterAssumptions;
                               n_points::Integer = 7)
    results = NamedTuple[]
    base_p = assumptions.cavity_pressure_erg_cm3[1]
    base_r = assumptions.cavity_radius_kpc[1]
    base_age = assumptions.cavity_age_yr[1]
    base_lcool = assumptions.cooling_luminosity_erg_s[1]
    ncav = assumptions.n_cavities

    factors = range(0.5, 2.0, length = n_points)

    for f in factors
        # Vary cavity radius
        p_r = cavity_power_erg_s(base_p, base_r * f, base_age; n_cavities = ncav)
        push!(results, (varied_parameter = "cavity_radius_kpc",
                        factor = f,
                        value = base_r * f,
                        feedback_ratio = feedback_ratio(p_r, base_lcool)))

        # Vary pressure
        p_p = cavity_power_erg_s(base_p * f, base_r, base_age; n_cavities = ncav)
        push!(results, (varied_parameter = "cavity_pressure_erg_cm3",
                        factor = f,
                        value = base_p * f,
                        feedback_ratio = feedback_ratio(p_p, base_lcool)))

        # Vary cavity age
        p_a = cavity_power_erg_s(base_p, base_r, base_age * f; n_cavities = ncav)
        push!(results, (varied_parameter = "cavity_age_yr",
                        factor = f,
                        value = base_age * f,
                        feedback_ratio = feedback_ratio(p_a, base_lcool)))

        # Vary cooling luminosity
        p_l = cavity_power_erg_s(base_p, base_r, base_age; n_cavities = ncav)
        push!(results, (varied_parameter = "cooling_luminosity_erg_s",
                        factor = f,
                        value = base_lcool * f,
                        feedback_ratio = feedback_ratio(p_l, base_lcool * f)))
    end

    return results
end

end
