# Cavity geometry calculations for AGN feedback analysis.
# Supports both spherical and ellipsoidal cavity models.

"""
    spherical_volume_cm3(radius_kpc)

Volume of a sphere in cm^3 given radius in kpc.
"""
function spherical_volume_cm3(radius_kpc::Real)
    r_cm = radius_kpc * KPC_TO_CM
    return 4.0 / 3.0 * pi * r_cm^3
end

"""
    ellipsoidal_volume_cm3(a_kpc, b_kpc, c_kpc)

Volume of an ellipsoid in cm^3 given semi-axes in kpc.
For projected cavities where line-of-sight axis is unknown,
the convention is c = sqrt(a * b) (geometric mean).
"""
function ellipsoidal_volume_cm3(a_kpc::Real, b_kpc::Real, c_kpc::Real)
    a_cm = a_kpc * KPC_TO_CM
    b_cm = b_kpc * KPC_TO_CM
    c_cm = c_kpc * KPC_TO_CM
    return 4.0 / 3.0 * pi * a_cm * b_cm * c_cm
end

"""
    cavity_enthalpy_from_volume(pressure_erg_cm3, volume_cm3; gamma=4/3)

Cavity enthalpy H = gamma/(gamma-1) * p * V.
For relativistic plasma gamma = 4/3, giving H = 4pV (the standard result).
"""
function cavity_enthalpy_from_volume(pressure_erg_cm3::Real, volume_cm3::Real;
                                     gamma::Real = 4.0 / 3.0)
    return gamma / (gamma - 1.0) * pressure_erg_cm3 * volume_cm3
end

"""
    cavity_power_from_enthalpy(enthalpy_erg, age_s; n_cavities=2)

Mechanical power = n_cavities * H / t_age.
"""
function cavity_power_from_enthalpy(enthalpy_erg::Real, age_s::Real;
                                    n_cavities::Integer = 2)
    return n_cavities * enthalpy_erg / age_s
end
