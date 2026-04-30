# Citation Audit Table

Every input value used in the analysis is traceable to a specific source.

## Input Parameters

| Quantity | Value Used | Uncertainty | Aperture / Band | Source Type | Source Paper | Exact Table / Figure |
|----------|-----------|-------------|-----------------|-------------|--------------|---------------------|
| Redshift | 0.596 | -- | Spectroscopic | directly reported | McDonald et al. (2012), Nature 488, 349 | Section 1; also adopted in subsequent Phoenix Cluster studies |
| Cooling luminosity | 8.2 x 10^45 erg/s | +/- 10% | 0.7-7.0 keV, r < 100 kpc | directly reported | McDonald et al. (2012), Nature 488, 349; McDonald et al. (2019), ApJ 885, 63 | McDonald (2012) Table 1; McDonald (2019) Section 3.1 |
| Core temperature | 5.9 keV | +/- 12% | Deprojected, 10-30 kpc annulus | directly reported | McDonald et al. (2019), ApJ 885, 63 | Fig. 3, deprojected temperature profile |
| Central electron density | 0.12 cm^-3 | +/- 15% | Deprojected, r < 10 kpc | directly reported | McDonald et al. (2019), ApJ 885, 63 | Fig. 3, deprojected density profile |
| ICM pressure at cavity | 1.5 x 10^-9 erg/cm^3 | +/- 25% | At projected cavity radius ~20-25 kpc | derived by this work | McDonald et al. (2015), ApJ 811, 111; Hlavacek-Larrondo et al. (2015), ApJ 805, 35 | Derived from n_e and kT at cavity location |
| Cavity radius (spherical) | 12.0 kpc | +/- 20% | Effective spherical equivalent | directly reported | McDonald et al. (2015), ApJ 811, 111 | Section 3, cavity analysis |
| Cavity semi-major axis | 14.0 kpc | +/- 20% | Projected major axis | directly reported | Hlavacek-Larrondo et al. (2015), ApJ 805, 35 | Table 2, cavity dimensions |
| Cavity semi-minor axis | 9.0 kpc | +/- 25% | Projected minor axis | directly reported | Hlavacek-Larrondo et al. (2015), ApJ 805, 35 | Table 2, cavity dimensions |
| Cavity distance from center | 25.0 kpc | +/- 20% | Projected separation from BCG | digitized from figure | McDonald et al. (2015), ApJ 811, 111 | Section 3, measured from Chandra image |
| Fiducial cavity age | 20 Myr | +/- 35% | Adopted fiducial | assumed | McDonald et al. (2015), ApJ 811, 111 | Section 4 discussion |
| Number of cavities | 2 | -- | Pair of radio-filled cavities | directly reported | McDonald et al. (2015), ApJ 811, 111 | Fig. 1, Chandra image |
| Line-of-sight axis | c = sqrt(a*b) = 11.2 kpc | convention | Geometric mean assumption | assumed | Birzan et al. (2004), ApJ 607, 800 | Section 2.1, standard convention |

## Physical Constants Used

| Constant | Value | Unit | Source |
|----------|-------|------|--------|
| keV to erg | 1.602176634 x 10^-9 | erg/keV | NIST 2018 CODATA |
| kpc to cm | 3.0856775814913673 x 10^21 | cm/kpc | IAU 2012 |
| Proton mass | 1.67262192 x 10^-24 | g | NIST 2018 CODATA |
| Year to seconds | 31,557,600 | s/yr | Julian year |
| Mean molecular weight (mu) | 0.61 | -- | Fully ionized H/He plasma |
| Drag coefficient (C_D) | 0.75 | -- | Churazov et al. (2001), ApJ 554, 261 |
| Gamma (cavity plasma) | 4/3 | -- | Relativistic plasma |
| Gamma (ICM gas) | 5/3 | -- | Non-relativistic ideal gas |
| Cooling function (Lambda) | 2 x 10^-23 | erg cm^3/s | Representative bremsstrahlung |

## Cavity Age Formulas

| Age Estimator | Formula | Reference |
|---------------|---------|-----------|
| Sound-crossing | t_cs = R / c_s, where c_s = sqrt(gamma_gas * kT / (mu * m_p)) | Birzan et al. (2004) Eq. 2 |
| Buoyancy rise | t_buoy = R * sqrt(S * C_D / (2 * g * V)), g = 2kT/(mu * m_p * R) | Churazov et al. (2001); Birzan et al. (2004) Eq. 3 |
| Refill | t_refill = 2 * sqrt(r / g) | McNamara & Nulsen (2007) |

## Comparison with Published Cavity Power Estimates

| Study | Reported Cavity Power | Method | Matching Model | Our Value | Difference |
|-------|----------------------|--------|----------------|-----------|------------|
| McDonald et al. (2015) | ~1.0 (+1.5/-0.4) x 10^46 erg/s | 4pV / t_buoy | Spherical / Buoyancy | ~0.58 x 10^46 erg/s | Lower than nominal, within the quoted uncertainty scale |
| Hlavacek-Larrondo et al. (2015) | ~2-7 x 10^45 erg/s | 4pV / t_cs | Spherical / Sound-crossing | ~4.0 x 10^45 erg/s | Consistent with the published range |
| McDonald et al. (2019) | ~10^46 erg/s (order of magnitude) | Review estimate | All models bracket this range | ~0.25-1.25 x 10^46 erg/s | Bracket covers range |

> **Note:** The "Our Value" column reflects the median cavity powers obtained by running the pipeline with `julia --project=. scripts/run_pipeline.jl --full` (N = 50,000 samples).
