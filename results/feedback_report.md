# Phoenix Cluster AGN Feedback Report

**Cluster:** SPT-CL J2344-4243 / Phoenix Cluster  
**Redshift:** 0.596  
**Monte Carlo samples:** 50000  

## 1. Single-Age Spherical Analysis (Legacy)

| Metric | Value |
|--------|-------|
| Median feedback/cooling ratio | 0.4957 |
| 16th-84th percentile interval | [0.2231, 1.039] |
| P(feedback >= cooling) | 0.1728 |

## 2. Multi-Age, Multi-Geometry Analysis

| Model | Median ratio | 16th pctl | 84th pctl | P(>=1) |
|-------|-------------|-----------|-----------|-------|
| Spherical / Sound-crossing | 0.4899 | 0.231 | 0.944 | 0.138 |
| Spherical / Buoyancy | 0.702 | 0.2935 | 1.523 | 0.33 |
| Spherical / Refill | 0.386 | 0.2032 | 0.6661 | 0.034 |
| Ellipsoidal / Sound-crossing | 0.3859 | 0.2022 | 0.6921 | 0.0472 |
| Ellipsoidal / Buoyancy | 0.532 | 0.2506 | 1.062 | 0.184 |
| Ellipsoidal / Refill | 0.3169 | 0.1826 | 0.516 | 0.00526 |

Median cooling time: 373.7 Myr

## 3. Interpretation

Under the adopted literature-anchored assumptions, Phoenix Cluster AGN cavity 
power is energetically comparable to the cooling luminosity, but most 
age/geometry combinations remain below unity at the median. The conclusion 
depends strongly on the adopted cavity age estimator and geometry model. 
The buoyancy-time estimator yields the highest feedback ratios in this 
implementation, while refill ages produce the lowest median ratios. 
Sound-crossing results fall between these regimes. Ellipsoidal geometry 
reduces cavity volume relative to the 
spherical assumption when the projected semi-minor axis is smaller than 
the semi-major axis.

These results should be treated as energetic consistency tests under 
transparent assumptions, not as definitive cavity analyses.
