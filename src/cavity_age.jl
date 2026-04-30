# Cavity age estimators for AGN cavity power calculations.
# Three standard methods: sound-crossing, buoyancy rise, and refill time.
#
# References:
#   Birzan et al. (2004), ApJ 607, 800 -- Eqs. 1-3
#   Churazov et al. (2001), ApJ 554, 261 -- buoyancy formulation
#   McNamara & Nulsen (2007), ARA&A 45, 117 -- review
#   Hlavacek-Larrondo et al. (2015), ApJ 805, 35

"""
    sound_crossing_time_s(distance_kpc, temperature_keV)

Sound-crossing time (Birzan et al. 2004):

    t_cs = R / c_s

where c_s = sqrt(gamma * kT / (mu * m_p)) is the adiabatic sound speed in
the ICM with gamma = 5/3 and mu = 0.61 (fully ionized plasma).
"""
function sound_crossing_time_s(distance_kpc::Real, temperature_keV::Real)
    gamma_gas = 5.0 / 3.0
    mu = 0.61  # mean molecular weight for fully ionized ICM
    kT_erg = temperature_keV * KEV_TO_ERG
    c_s = sqrt(gamma_gas * kT_erg / (mu * PROTON_MASS_G))
    d_cm = distance_kpc * KPC_TO_CM
    return d_cm / c_s
end

"""
    buoyancy_time_s(distance_kpc, cavity_radius_kpc, temperature_keV)

Buoyancy rise time (Churazov et al. 2001; Birzan et al. 2004):

    t_buoy = R * sqrt(S * C_D / (2 * g * V))

where:
  S = pi * r^2     cross-sectional area of the cavity
  C_D = 0.75       drag coefficient (Churazov et al. 2001)
  V = 4/3 pi r^3   cavity volume (spherical approx for age estimate)
  g = 2 kT / (mu m_p R)  gravitational acceleration (isothermal sphere)

This simplifies to:
    t_buoy = R * sqrt(3 * C_D / (8 * g * r))
"""
function buoyancy_time_s(distance_kpc::Real, cavity_radius_kpc::Real,
                          temperature_keV::Real)
    C_D = 0.75  # drag coefficient (Churazov et al. 2001)
    mu = 0.61
    kT_erg = temperature_keV * KEV_TO_ERG
    R_cm = distance_kpc * KPC_TO_CM
    r_cm = cavity_radius_kpc * KPC_TO_CM

    # Gravitational acceleration assuming isothermal sphere
    g = 2.0 * kT_erg / (mu * PROTON_MASS_G * R_cm)

    # Cross-section and volume of spherical cavity
    S = pi * r_cm^2
    V = 4.0 / 3.0 * pi * r_cm^3

    return R_cm * sqrt(S * C_D / (2.0 * g * V))
end

"""
    refill_time_s(cavity_radius_kpc, distance_kpc, temperature_keV)

Refill time (McNamara & Nulsen 2007):

    t_refill = 2 * sqrt(r / g)

where r is the cavity radius and g is the gravitational acceleration
at the cavity position, estimated as g = 2 kT / (mu m_p R).

This is the dynamical time for the surrounding atmosphere to collapse
into the volume vacated by the rising cavity.
"""
function refill_time_s(cavity_radius_kpc::Real, distance_kpc::Real,
                        temperature_keV::Real)
    mu = 0.61
    kT_erg = temperature_keV * KEV_TO_ERG
    R_cm = distance_kpc * KPC_TO_CM
    r_cm = cavity_radius_kpc * KPC_TO_CM

    g = 2.0 * kT_erg / (mu * PROTON_MASS_G * R_cm)

    return 2.0 * sqrt(r_cm / g)
end
