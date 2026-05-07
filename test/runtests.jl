include(joinpath(@__DIR__, "..", "src", "PhoenixFeedback.jl"))
using .PhoenixFeedback
using JSON
using Test

# -- Reference values computed deterministically from nominal parameters ------
function _deterministic_reference()
    p = 1.5e-9
    r = 12.0
    a, b = 14.0, 9.0
    c = sqrt(a * b)
    R = 25.0
    T = 5.9
    ncav = 2
    Lcool = 8.2e45

    vol_sph = spherical_volume_cm3(r)
    vol_ell = ellipsoidal_volume_cm3(a, b, c)
    r_eff = cbrt(a * b * c)
    H_sph = cavity_enthalpy_from_volume(p, vol_sph)
    H_ell = cavity_enthalpy_from_volume(p, vol_ell)
    t_cs = sound_crossing_time_s(R, T)
    t_buoy = buoyancy_time_s(R, r, T)
    t_refill = refill_time_s(r, R, T)

    P_cs = cavity_power_from_enthalpy(H_sph, t_cs; n_cavities = ncav)
    P_buoy = cavity_power_from_enthalpy(H_sph, t_buoy; n_cavities = ncav)

    return (H_sph = H_sph, H_ell = H_ell, r_eff = r_eff,
            t_cs_yr = t_cs / PhoenixFeedback.YEAR_TO_SECONDS,
            t_buoy_yr = t_buoy / PhoenixFeedback.YEAR_TO_SECONDS,
            t_refill_yr = t_refill / PhoenixFeedback.YEAR_TO_SECONDS,
            P_sph_cs = P_cs, P_sph_buoy = P_buoy,
            ratio_cs = feedback_ratio(P_cs, Lcool),
            ratio_buoy = feedback_ratio(P_buoy, Lcool))
end

@testset "PhoenixFeedback" begin

    # -- Legacy physics functions -----------------------------------------
    @testset "Core physics" begin
        @test cooling_time_seconds(8.0, 0.08) > 0.0
        @test cavity_enthalpy_erg(1.0e-9, 20.0) > 0.0
        @test cavity_power_erg_s(1.0e-9, 20.0, 3.0e7; n_cavities = 2) > 0.0
        @test feedback_ratio(2.0, 4.0) == 0.5
    end

    @testset "Numerical regression" begin
        ref = _deterministic_reference()
        # Sound-crossing time ~ 20 Myr for Phoenix
        @test ref.t_cs_yr ≈ 2.13e7 rtol = 0.1
        # Buoyancy time ~ 14 Myr for Phoenix
        @test ref.t_buoy_yr ≈ 1.37e7 rtol = 0.1
        # Spherical buoyancy power ~ 6e45 erg/s (McDonald 2015 scale)
        @test ref.P_sph_buoy ≈ 5.88e45 rtol = 0.2
        # Sound-crossing ratio ~ 0.5
        @test ref.ratio_cs ≈ 0.50 rtol = 0.2
        # Buoyancy ratio ~ 0.72
        @test ref.ratio_buoy ≈ 0.72 rtol = 0.2
        @test ref.H_sph > ref.H_ell  # spherical volume > ellipsoidal volume
        @test ref.r_eff ≈ 11.22 rtol = 0.02
    end

    @testset "Input validation" begin
        path = joinpath(@__DIR__, "..", "data", "assumptions.json")
        raw = JSON.parsefile(path)
        @test raw["redshift"] > 0.0
        @test raw["n_cavities"] > 0
        @test raw["monte_carlo_samples"] > 0
        for key in ["cooling_luminosity_erg_s", "gas_temperature_keV",
                     "electron_density_cm3", "cavity_pressure_erg_cm3",
                     "cavity_radius_kpc", "cavity_age_yr", "cavity_distance_kpc"]
            b = raw[key]
            @test b["median"] > 0.0
            @test 0.0 < b["sigma_fraction"] < 1.0
        end
        @test raw["cavity_geometry"]["semi_major_kpc"]["median"] > 0
        @test raw["cavity_geometry"]["semi_minor_kpc"]["median"] > 0
    end

    # -- Geometry ---------------------------------------------------------
    @testset "Geometry" begin
        v_sph = spherical_volume_cm3(10.0)
        @test v_sph > 0.0

        v_ell = ellipsoidal_volume_cm3(14.0, 9.0, 11.0)
        @test v_ell > 0.0

        # Sphere with r=10 should equal ellipsoid with a=b=c=10
        v_ell_sphere = ellipsoidal_volume_cm3(10.0, 10.0, 10.0)
        @test isapprox(v_sph, v_ell_sphere, rtol = 1e-10)

        H = cavity_enthalpy_from_volume(1e-9, v_sph)
        @test H > 0.0

        P = cavity_power_from_enthalpy(H, 1e15; n_cavities = 2)
        @test P > 0.0
    end

    # -- Cavity ages ------------------------------------------------------
    @testset "Cavity ages" begin
        t_cs = sound_crossing_time_s(25.0, 6.0)
        @test t_cs > 0.0

        t_buoy = buoyancy_time_s(25.0, 12.0, 6.0)
        @test t_buoy > 0.0

        t_ref = refill_time_s(12.0, 25.0, 6.0)
        @test t_ref > 0.0

        # All three ages should be physically reasonable (positive, finite)
        @test isfinite(t_cs)
        @test isfinite(t_buoy)
        @test isfinite(t_ref)
    end

    # -- Data loading -----------------------------------------------------
    @testset "Load assumptions" begin
        assumptions = load_assumptions(joinpath(@__DIR__, "..", "data", "assumptions.json"))
        @test assumptions.cluster_name == "SPT-CL J2344-4243 / Phoenix Cluster"
        @test assumptions.redshift ≈ 0.596
        @test assumptions.n_cavities == 2
        @test assumptions.cavity_semi_major_kpc[1] > 0.0
        @test assumptions.cavity_semi_minor_kpc[1] > 0.0
        @test assumptions.cavity_distance_kpc[1] > 0.0
    end

    # -- Legacy Monte Carlo -----------------------------------------------
    @testset "Legacy Monte Carlo" begin
        assumptions = load_assumptions(joinpath(@__DIR__, "..", "data", "assumptions.json"))
        samples = run_monte_carlo(assumptions; seed = 1)
        summary = summarize_samples(samples)
        @test summary["n_samples"] == assumptions.monte_carlo_samples
        @test summary["median_feedback_to_cooling_ratio"] > 0.0
    end

    # -- Multi-age Monte Carlo --------------------------------------------
    @testset "Multi-age Monte Carlo" begin
        assumptions = load_assumptions(joinpath(@__DIR__, "..", "data", "assumptions.json"))
        samples = run_multi_age_monte_carlo(assumptions; seed = 1)
        @test length(samples) == assumptions.monte_carlo_samples

        s = samples[1]
        @test s.ratio_sph_cs > 0.0
        @test s.ratio_sph_buoy > 0.0
        @test s.ratio_sph_refill > 0.0
        @test s.ratio_ell_cs > 0.0
        @test s.ratio_ell_buoy > 0.0
        @test s.ratio_ell_refill > 0.0
        @test s.volume_spherical_cm3 > 0.0
        @test s.volume_ellipsoidal_cm3 > 0.0
        @test s.enthalpy_spherical_erg > 0.0
        @test s.enthalpy_ellipsoidal_erg > 0.0
        @test s.effective_ellipsoid_radius_kpc > 0.0
        @test s.t_ell_buoyancy_yr > 0.0
        @test s.t_ell_refill_yr > 0.0

        msummary = summarize_multi_age_samples(samples)
        @test msummary["n_samples"] == assumptions.monte_carlo_samples
        for key in ["sph_soundcross", "sph_buoyancy", "sph_refill",
                     "ell_soundcross", "ell_buoyancy", "ell_refill"]
            @test haskey(msummary, key)
            @test msummary[key]["median"] > 0.0
        end
    end

    @testset "Monte Carlo reproducibility" begin
        assumptions = load_assumptions(joinpath(@__DIR__, "..", "data", "assumptions.json"))
        small = ClusterAssumptions(
            assumptions.cluster_name,
            assumptions.redshift,
            assumptions.cooling_luminosity_erg_s,
            assumptions.gas_temperature_keV,
            assumptions.electron_density_cm3,
            assumptions.cavity_pressure_erg_cm3,
            assumptions.cavity_radius_kpc,
            assumptions.cavity_semi_major_kpc,
            assumptions.cavity_semi_minor_kpc,
            assumptions.cavity_age_yr,
            assumptions.cavity_distance_kpc,
            assumptions.n_cavities,
            200,
        )

        same_a = run_multi_age_monte_carlo(small; seed = 7)
        same_b = run_multi_age_monte_carlo(small; seed = 7)
        different = run_multi_age_monte_carlo(small; seed = 8)

        @test same_a[1] == same_b[1]
        @test same_a[1] != different[1]
    end

    @testset "Legacy MC reproducibility" begin
        assumptions = load_assumptions(joinpath(@__DIR__, "..", "data", "assumptions.json"))
        small = ClusterAssumptions(
            assumptions.cluster_name, assumptions.redshift,
            assumptions.cooling_luminosity_erg_s,
            assumptions.gas_temperature_keV,
            assumptions.electron_density_cm3,
            assumptions.cavity_pressure_erg_cm3,
            assumptions.cavity_radius_kpc,
            assumptions.cavity_semi_major_kpc,
            assumptions.cavity_semi_minor_kpc,
            assumptions.cavity_age_yr,
            assumptions.cavity_distance_kpc,
            assumptions.n_cavities, 200,
        )
        a = run_monte_carlo(small; seed = 42)
        b = run_monte_carlo(small; seed = 42)
        c = run_monte_carlo(small; seed = 99)
        @test a[1] == b[1]
        @test a[1] != c[1]
    end

    @testset "Published-value recovery" begin
        # Deterministic (no uncertainty) check: with nominal Phoenix inputs,
        # pipeline should recover McDonald et al. (2015) and
        # Hlavacek-Larrondo et al. (2015) order-of-magnitude cavity powers.
        p = 1.5e-9; r = 12.0; R = 25.0; T = 5.9; ncav = 2
        H = cavity_enthalpy_from_volume(p, spherical_volume_cm3(r))
        t_buoy = buoyancy_time_s(R, r, T)
        t_cs   = sound_crossing_time_s(R, T)

        P_buoy = cavity_power_from_enthalpy(H, t_buoy; n_cavities = ncav)
        P_cs   = cavity_power_from_enthalpy(H, t_cs;   n_cavities = ncav)

        # McDonald 2015: P_cav ~ 1.0e46 (+1.5/-0.4) erg/s, buoyancy
        @test P_buoy ≈ 1.0e46 rtol = 0.5

        # Hlavacek-Larrondo 2015: P_cav ~ 2-7e45 erg/s, sound-crossing
        @test 2e45 <= P_cs <= 7e45
    end

    # -- Sensitivity grid -------------------------------------------------
    @testset "Sensitivity grid" begin
        assumptions = load_assumptions(joinpath(@__DIR__, "..", "data", "assumptions.json"))
        grid = run_sensitivity_grid(assumptions; n_points = 5)
        @test length(grid) == 5 * 4  # 5 points x 4 parameters
        @test all(r -> r.feedback_ratio > 0.0, grid)
    end

    @testset "Multi-age sensitivity grid" begin
        assumptions = load_assumptions(joinpath(@__DIR__, "..", "data", "assumptions.json"))
        grid = run_multi_age_sensitivity_grid(assumptions; n_points = 5)
        expected = 5 * 5 * 2
        @test length(grid) == expected
        @test all(r -> r.feedback_ratio > 0.0, grid)
        @test all(r -> r.model in ("Sph/Sound", "Sph/Buoyancy"), grid)
    end

    # -- Two-cavity Monte Carlo --------------------------------------------
    @testset "Two-cavity Monte Carlo" begin
        assumptions = load_assumptions(joinpath(@__DIR__, "..", "data", "assumptions.json"))
        cavities = load_cavities(joinpath(@__DIR__, "..", "data", "cavities.json"))
        @test length(cavities) == 2

        small = ClusterAssumptions(
            assumptions.cluster_name, assumptions.redshift,
            assumptions.cooling_luminosity_erg_s,
            assumptions.gas_temperature_keV,
            assumptions.electron_density_cm3,
            assumptions.cavity_pressure_erg_cm3,
            assumptions.cavity_radius_kpc,
            assumptions.cavity_semi_major_kpc,
            assumptions.cavity_semi_minor_kpc,
            assumptions.cavity_age_yr,
            assumptions.cavity_distance_kpc,
            assumptions.n_cavities, 100,
        )
        samples = run_two_cavity_monte_carlo(small, cavities; seed = 1)
        @test length(samples) == 100
        s = samples[1]
        @test s.ratio_sph_cs > 0.0
        @test s.ratio_sph_buoy > 0.0
        @test s.ratio_ell_cs > 0.0
        @test s.ratio_ell_buoy > 0.0
        @test s.proj_factor == 1.0

        # Each ratio should be similar to the representative-cavity case
        msummary = summarize_two_cavity_samples(samples)
        @test msummary["n_samples"] == 100
        @test haskey(msummary["sph_soundcross"], "median")

        # Reproducibility
        same = run_two_cavity_monte_carlo(small, cavities; seed = 7)
        same2 = run_two_cavity_monte_carlo(small, cavities; seed = 7)
        diff = run_two_cavity_monte_carlo(small, cavities; seed = 8)
        @test same[1] == same2[1]
        @test same[1] != diff[1]
    end

    @testset "Projection sensitivity" begin
        assumptions = load_assumptions(joinpath(@__DIR__, "..", "data", "assumptions.json"))
        cavities = load_cavities(joinpath(@__DIR__, "..", "data", "cavities.json"))
        small = ClusterAssumptions(
            assumptions.cluster_name, assumptions.redshift,
            assumptions.cooling_luminosity_erg_s,
            assumptions.gas_temperature_keV,
            assumptions.electron_density_cm3,
            assumptions.cavity_pressure_erg_cm3,
            assumptions.cavity_radius_kpc,
            assumptions.cavity_semi_major_kpc,
            assumptions.cavity_semi_minor_kpc,
            assumptions.cavity_age_yr,
            assumptions.cavity_distance_kpc,
            assumptions.n_cavities, 100,
        )
        results = run_projection_sensitivity(small, cavities; seed = 42,
                                              factors = [1.0, 1.5, 2.0])
        @test length(results) == 3
        @test results[1]["proj_factor"] == 1.0
        @test results[2]["proj_factor"] == 1.5
        @test results[3]["proj_factor"] == 2.0
        # Higher projection factor → lower ratio (larger distance → larger age → lower power)
        r1 = results[1]["sph_buoyancy"]["median"]
        r2 = results[2]["sph_buoyancy"]["median"]
        @test r2 < r1
    end

    # -- Integration: end-to-end data flow ---------------------------------
    @testset "Integration: end-to-end data flow" begin
        assumptions = load_assumptions(joinpath(@__DIR__, "..", "data", "assumptions.json"))
        small = ClusterAssumptions(
            assumptions.cluster_name,
            assumptions.redshift,
            assumptions.cooling_luminosity_erg_s,
            assumptions.gas_temperature_keV,
            assumptions.electron_density_cm3,
            assumptions.cavity_pressure_erg_cm3,
            assumptions.cavity_radius_kpc,
            assumptions.cavity_semi_major_kpc,
            assumptions.cavity_semi_minor_kpc,
            assumptions.cavity_age_yr,
            assumptions.cavity_distance_kpc,
            assumptions.n_cavities,
            100,
        )

        # Multi-age MC → summary → expected keys
        samples = run_multi_age_monte_carlo(small; seed = 42)
        @test length(samples) == 100

        msummary = summarize_multi_age_samples(samples)
        for key in ["sph_soundcross", "sph_buoyancy", "sph_refill",
                     "ell_soundcross", "ell_buoyancy", "ell_refill"]
            @test haskey(msummary, key)
            s = msummary[key]
            @test haskey(s, "label")
            @test haskey(s, "median")
            @test haskey(s, "p16")
            @test haskey(s, "p84")
            @test haskey(s, "p_exceeds_1")
            @test 0.0 <= s["p_exceeds_1"] <= 1.0
            @test s["p16"] <= s["median"] <= s["p84"]
        end
        @test haskey(msummary, "median_cooling_time_yr")
        @test msummary["median_cooling_time_yr"] > 0.0

        # Legacy MC flow
        legacy_samples = run_monte_carlo(small; seed = 42)
        @test length(legacy_samples) == 100
        legacy_summary = summarize_samples(legacy_samples)
        @test haskey(legacy_summary, "median_feedback_to_cooling_ratio")
        @test haskey(legacy_summary, "probability_feedback_exceeds_cooling")
        @test 0.0 <= legacy_summary["probability_feedback_exceeds_cooling"] <= 1.0

        # Sensitivity grid flow
        grid = run_sensitivity_grid(small; n_points = 5)
        @test length(grid) == 20

        # Multi-age sensitivity grid flow
        multi_grid = run_multi_age_sensitivity_grid(small; n_points = 5)
        @test length(multi_grid) == 5 * 5 * 2

        # Two-cavity + projection flow
        cavities = load_cavities(joinpath(@__DIR__, "..", "data", "cavities.json"))
        two_cav = run_two_cavity_monte_carlo(small, cavities; seed = 42)
        @test length(two_cav) == 100
        tsum = summarize_two_cavity_samples(two_cav)
        @test haskey(tsum, "sph_soundcross")

        proj = run_projection_sensitivity(small, cavities; seed = 42, factors = [1.0, 2.0])
        @test length(proj) == 2
    end
end
