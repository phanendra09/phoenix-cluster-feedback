include(joinpath(@__DIR__, "..", "src", "PhoenixFeedback.jl"))
using .PhoenixFeedback
using Test

@testset "PhoenixFeedback" begin

    # -- Legacy physics functions -----------------------------------------
    @testset "Core physics" begin
        @test cooling_time_seconds(8.0, 0.08) > 0.0
        @test cavity_enthalpy_erg(1.0e-9, 20.0) > 0.0
        @test cavity_power_erg_s(1.0e-9, 20.0, 3.0e7; n_cavities = 2) > 0.0
        @test feedback_ratio(2.0, 4.0) == 0.5
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

    # -- Sensitivity grid -------------------------------------------------
    @testset "Sensitivity grid" begin
        assumptions = load_assumptions(joinpath(@__DIR__, "..", "data", "assumptions.json"))
        grid = run_sensitivity_grid(assumptions; n_points = 5)
        @test length(grid) == 5 * 4  # 5 points x 4 parameters
        @test all(r -> r.feedback_ratio > 0.0, grid)
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
    end
end
