# Paper Notes

## Working Title

AGN Feedback Versus Radiative Cooling in the Phoenix Cluster:
A Multi-Age, Multi-Geometry Cavity Analysis

## Scope

We compare Phoenix Cluster cooling losses and AGN mechanical feedback using
literature-based measurements with full uncertainty propagation. The analysis
reports feedback-to-cooling power ratios under three cavity age definitions
(sound-crossing, buoyancy, refill) and two geometries (spherical, ellipsoidal),
providing a comprehensive energetic consistency test.

## Key Citations

- McDonald et al. (2012): Discovery, cooling luminosity, SFR
- McDonald et al. (2015): Cavity detection, AGN properties
- McDonald et al. (2018): Phoenix thermodynamic structure, cool-core classification
- McDonald et al. (2019): Deprojected profiles, cooling anatomy
- Hlavacek-Larrondo et al. (2015): Cavity survey, SPT sample
- Russell et al. (2017): ALMA molecular gas
- Birzan et al. (2004): Cavity methodology, ellipsoidal convention
- McNamara & Nulsen (2007): AGN feedback review
- Churazov et al. (2001): Buoyancy age formulation
- Rafferty et al. (2006): Cavity power vs cooling in large sample
- Dunn & Fabian (2006): AGN heating efficiency
- Nulsen et al. (2005): Shock heating, MS 0735
- McNamara et al. (2005): MS 0735 extreme AGN outburst
- Fabian et al. (2006): Perseus sound waves / heating
- Fabian et al. (2017): Sound waves, turbulence heating
- Wise et al. (2007): Hydra A cavity analysis
- Voit et al. (2015): Precipitation-driven feedback
- Gaspari et al. (2013): Chaotic cold accretion

## Status

- [x] Input parameters anchored to cited values with aperture metadata
- [x] Multi-age analysis implemented (sound-crossing, buoyancy, refill)
- [x] Ellipsoidal geometry support added
- [x] Sensitivity analysis implemented
- [x] Literature comparison table created
- [x] BibTeX reference file created
- [x] Physical context (SFR, molecular gas, cooling rate) included
- [x] Conservative, conditional conclusions written
- [x] 50K-sample full Monte Carlo results regenerated
- [x] ApJ-compliant figures (colorblind-safe, LaTeX labels, 300 dpi)
- [x] Cluster comparison table (Phoenix, MS 0735, Perseus, Hydra A)
- [x] Expanded introduction with cool-core landscape context
- [x] Buoyancy time derivation added to §4.2
- [x] §9 reframed as Discussion with future work subsection
- [x] Systematic uncertainties expanded (§8.7 spherical approx, §8.8 cooling function)
- [x] ApJ-required metadata (keywords, ORCID, acknowledgments, data availability, facility, software)
- [x] References expanded to 19 entries
- [x] All 163 tests passing
- [ ] Replace approximate values with exact table entries if deeper archival analysis is done (requires archival data access)
- [x] Acknowledgments drafted (funding line item needs author input before submission)
- [x] Converted to AASTeX LaTeX (paper/manuscript.tex) — requires `\author`, `\affil`, `\email` fields to be filled by author
- [x] Two-cavity Monte Carlo (independent per-cavity sampling) implemented
- [x] Projection sensitivity analysis (R_true = [1.0, 1.2, 1.5, 2.0] × R_proj) implemented
- [x] Projection sensitivity figure (Fig 6) added
- [x] GitHub Actions CI configured (Julia tests + LaTeX build)

## Limitations (Acknowledged in Manuscript)

- Cavity geometry is approximate; real morphology is irregular
- Line-of-sight axis is assumed, not measured
- Age estimates are analytical, not from hydro simulations
- Cooling function is simplified (constant Lambda)
- No conduction, turbulence, or cosmic-ray heating
- ICM profiles are taken at representative radii, not full radial models
- Spherical approximation used in ellipsoidal buoyancy time (§8.7)
- Cooling function constant vs full APEC spectral model (§8.8)
