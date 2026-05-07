# Phoenix Cluster Feedback

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Language: Julia](https://img.shields.io/badge/Language-Julia%201.12-9558B2.svg)](https://julialang.org)
[![Monte Carlo: 50,000 samples](https://img.shields.io/badge/MC%20Samples-50%2C000-brightgreen.svg)](#key-results)
[![Reproducible: seed 42](https://img.shields.io/badge/Seed-42-orange.svg)](#reproducibility)

Reproducible Julia pipeline for evaluating whether AGN mechanical feedback in
the Phoenix Cluster (SPT-CL J2344-4243, z = 0.596) is energetically consistent with
offsetting its extreme radiative cooling flow. Built with literature-anchored
measurements, full uncertainty propagation, and transparent provenance tracking.

---

## Scientific Question

> **Is the central AGN feedback power in the Phoenix Cluster energetically
> comparable to the extreme cooling flow in its core?**

The Phoenix Cluster hosts the most extreme known cooling flow in the universe
(~3800 M☉/yr classical rate; McDonald et al. 2012), a massive starburst
(500-800 M☉/yr), and a powerful central AGN with X-ray cavities. This pipeline
quantifies whether the mechanical energy stored in those cavities can offset the
radiative losses.

---

## Pipeline Architecture

```mermaid
graph TD
    A["📄 Literature Inputs<br/>data/assumptions.json"] --> B("⚙️ PhoenixFeedback.jl<br/>Physics & MC Engine")
    
    B --> C{"🔷 Geometry Model"}
    C -->|Spherical| D["V = 4/3 π r³"]
    C -->|Ellipsoidal| E["V = 4/3 π a b c"]
    
    D --> F{"⏱️ Cavity Age Estimator"}
    E --> F
    
    F -->|Sound-Crossing| G["t_cs = R / c_s"]
    F -->|Buoyancy Rise| H["t_buoy with C_D = 0.75"]
    F -->|Refill| I["t_refill = 2 √(r/g)"]
    
    G --> J["🔥 Cavity Power<br/>P = 4pV / t"]
    H --> J
    I --> J
    
    J --> K["📊 Feedback Ratio<br/>R = P_cav / L_cool"]
    K --> L(("📁 Outputs"))

    style A fill:#1a1a2e,stroke:#e94560,color:#eee
    style B fill:#16213e,stroke:#0f3460,color:#eee
    style J fill:#0f3460,stroke:#e94560,color:#eee
    style K fill:#1a1a2e,stroke:#53d8fb,color:#eee
    style L fill:#533483,stroke:#e94560,color:#eee
```

## Monte Carlo Methodology

```mermaid
graph LR
    subgraph Inputs
        P1["L_cool ± 10%"]
        P2["kT ± 12%"]
        P3["n_e ± 15%"]
        P4["P_cav ± 25%"]
        P5["r_cav ± 20%"]
        P6["R_dist ± 20%"]
    end

    subgraph Sampling
        S["50,000 Draws<br/>Truncated Gaussians<br/>Seed = 42"]
    end

    subgraph Models
        M1["2 Geometries<br/>Spherical / Ellipsoidal"]
        M2["3 Age Methods<br/>Sound / Buoy / Refill"]
        M3["= 6 Model<br/>Combinations"]
    end

    subgraph Outputs
        O1["Median Ratio"]
        O2["16th-84th CI"]
        O3["P(R ≥ 1)"]
    end

    P1 --> S
    P2 --> S
    P3 --> S
    P4 --> S
    P5 --> S
    P6 --> S
    S --> M1
    M1 --> M2
    M2 --> M3
    M3 --> O1
    M3 --> O2
    M3 --> O3

    style S fill:#0f3460,stroke:#53d8fb,color:#eee
    style M3 fill:#533483,stroke:#e94560,color:#eee
```

---

## Key Results

### Summary Table

| Model | Geometry | Age Method | Median Ratio | 16th pctl | 84th pctl | P(≥1) |
|-------|----------|------------|:------------:|:---------:|:---------:|:-----:|
| **Sph / Sound** | Spherical | Sound-crossing | 0.490 | 0.231 | 0.944 | 13.8% |
| **Sph / Buoyancy** | Spherical | Buoyancy rise | **0.702** | 0.293 | 1.523 | **33.0%** |
| **Sph / Refill** | Spherical | Refill | 0.386 | 0.203 | 0.666 | 3.4% |
| **Ell / Sound** | Ellipsoidal | Sound-crossing | 0.386 | 0.202 | 0.692 | 4.7% |
| **Ell / Buoyancy** | Ellipsoidal | Buoyancy rise | 0.532 | 0.251 | 1.062 | 18.4% |
| **Ell / Refill** | Ellipsoidal | Refill | 0.317 | 0.183 | 0.516 | 0.5% |

> **Conclusion:** Phoenix Cluster cavity power is comparable to the cooling
> luminosity but remains below unity for most standard age/geometry assumptions,
> suggesting that cavities alone may not fully offset cooling unless additional
> heating channels contribute.

### Feedback Ratio Distribution (Single-Age, Spherical)

Distribution of the feedback-to-cooling power ratio from the legacy single-age
spherical Monte Carlo (N = 50,000). The dashed line marks the ratio = 1
balance point.

![Feedback Ratio Histogram](figures/feedback_ratio_histogram.png)

### Multi-Age Ratio Comparison

Comparison of feedback-to-cooling ratios across all six model configurations.
The buoyancy-rise estimator consistently produces the highest ratios, while
ellipsoidal geometry reduces all estimates by ~20-25%.

![Multi-Age Comparison](figures/multi_age_ratio_comparison.png)

### Parameter Sensitivity Analysis

One-at-a-time sensitivity scan showing how the feedback ratio responds to
perturbations in each input parameter. Cavity radius (cubic volume dependence)
and cavity age (inverse dependence) dominate the uncertainty budget.

![Sensitivity Heatmap](figures/sensitivity_heatmap.png)

### Cooling Luminosity vs Cavity Power

Scatter plot of cooling luminosity versus computed cavity power for individual
Monte Carlo samples (log-log scale), illustrating the correlation structure.

![Cooling vs Feedback](figures/cooling_vs_feedback.png)

### Cluster Comparison

Feedback-cooling ratio comparison across four well-studied cool-core clusters
(Phoenix, MS 0735, Perseus, Hydra A) for buoyancy-age and sound-crossing-age
estimates. Phoenix's ratio is comparable to MS 0735 under buoyancy-age
assumptions but its extreme absolute cooling rate makes the residual cooling
far more significant.

![Cluster Comparison](figures/cluster_comparison.png)

---

## Physical Context

```mermaid
graph LR
    subgraph PhoenixCluster["Phoenix Cluster at z = 0.596"]
        A["Central AGN"] -->|"Jets inflate<br/>X-ray cavities"| B["Cavity Pair"]
        B -->|"Mechanical<br/>work: 4pV"| C["ICM Heating"]
        D["Cooling Flow<br/>3800 Msun/yr"] -->|"Radiative<br/>losses"| E["Starburst<br/>500-800 Msun/yr"]
        C -.->|"Partially<br/>offsets?"| D
    end

    style A fill:#e94560,stroke:#1a1a2e,color:#fff
    style B fill:#0f3460,stroke:#53d8fb,color:#eee
    style D fill:#16213e,stroke:#53d8fb,color:#eee
    style E fill:#533483,stroke:#e94560,color:#eee
```

The Phoenix Cluster occupies a unique position in the cooling-flow landscape:

- **Cooling rate:** Classical mass deposition rate of ~3800 M☉/yr (McDonald et al. 2012)
- **Star formation:** 500-800 M☉/yr in the BCG, the highest known for any cluster
- **Molecular gas:** ~2 x 10¹⁰ M☉ detected in CO(3-2) by ALMA (Russell et al. 2017)
- **AGN cavities:** A pair of radio-filled X-ray cavities detected by Chandra (McDonald et al. 2015)

---

## Installation

### Prerequisites
- [Julia 1.10+](https://julialang.org/downloads/) (tested on Julia 1.12.6)

### Setup

```bash
# Clone the repository
git clone https://github.com/phanendra09/phoenix-cluster-feedback.git
cd phoenix-cluster-feedback

# Install dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Run tests
julia --project=. test/runtests.jl
```

## Usage

### Quick Mode (2,000 samples, ~30 seconds)
```bash
julia --project=. scripts/run_pipeline.jl --quick
```

### Full Mode (50,000 samples, ~3 minutes)
```bash
julia --project=. scripts/run_pipeline.jl --full
```

### Reproducible Run with Explicit Seed
```bash
julia --project=. scripts/run_pipeline.jl --full --seed 42
```

---

## Outputs

| File | Description |
|------|-------------|
| `results/feedback_summary.json` | Legacy single-age summary statistics |
| `results/multi_age_summary.json` | Multi-age, multi-geometry summary |
| `results/monte_carlo_samples.csv` | Legacy MC samples |
| `results/multi_age_samples.csv` | Full MC samples (6 ratio columns) |
| `results/sensitivity_grid.csv` | Sensitivity analysis grid |
| `results/derived_quantities.csv` | Derived physics audit table |
| `results/run_metadata.json` | Reproducibility metadata (seed, Julia version, timestamp) |
| `results/feedback_report.md` | Human-readable report |
| `figures/*.png` | Publication-quality plots |
| `figures/*.pdf` | Vector plots for LaTeX inclusion |

## Reproducibility

Every pipeline run generates `results/run_metadata.json` containing:

```json
{
    "run_mode": "full",
    "sample_count": 50000,
    "seed": 42,
    "julia_version": "1.12.6",
    "timestamp_utc": "2026-04-30T16:21:57.102",
    "assumption_file_path": "data/assumptions.json",
    "assumption_file_sha256": "2c2838..."
}
```

---

## Project Structure

```mermaid
graph TD
    subgraph "Source Code"
        S1["src/PhoenixFeedback.jl<br/>Core module"]
        S2["src/geometry.jl<br/>Volume calculations"]
        S3["src/cavity_age.jl<br/>Age estimators"]
    end

    subgraph "Pipeline Scripts"
        P1["scripts/run_pipeline.jl<br/>Main entry point"]
        P2["scripts/make_figures.jl"]
        P3["scripts/write_report.jl"]
        P4["scripts/sensitivity_grid.jl"]
    end

    subgraph "Data & Paper"
        D1["data/assumptions.json"]
        M1["paper/manuscript.md"]
        M2["paper/references.bib"]
        M3["paper/tables/literature_values.md"]
    end

    subgraph "Tests"
        T1["test/runtests.jl<br/>45 tests"]
    end

    P1 --> S1
    S1 --> S2
    S1 --> S3
    P1 --> P2
    P1 --> P3
    P1 --> P4
    D1 --> P1
    T1 --> S1

    style S1 fill:#0f3460,stroke:#53d8fb,color:#eee
    style P1 fill:#533483,stroke:#e94560,color:#eee
    style D1 fill:#1a1a2e,stroke:#e94560,color:#eee
    style T1 fill:#16213e,stroke:#53d8fb,color:#eee
```

```text
src/
  PhoenixFeedback.jl       Core module: physics, MC, sensitivity
  geometry.jl              Spherical and ellipsoidal cavity volumes
  cavity_age.jl            Sound-crossing, buoyancy, refill ages
scripts/
  run_pipeline.jl          End-to-end reproducible run (--quick / --full / --seed)
  make_figures.jl          Publication-quality figure generation
  sensitivity_grid.jl      Standalone sensitivity analysis
  write_report.jl          Markdown report generation
data/
  assumptions.json         Cited input assumptions with aperture metadata
  literature_targets.csv   Citation tracking checklist
paper/
  manuscript.md            Journal manuscript (ApJ/MNRAS ready)
  references.bib           BibTeX database (15 entries)
  tables/                  Literature comparison & citation audit tables
  notes.md                 Working notes & status checklist
test/
  runtests.jl              45-test suite covering physics, geometry, and MC
```

---

## Cavity Age Formulas

The three independent cavity age estimators implemented in this pipeline:

| Estimator | Formula | Reference |
|-----------|---------|-----------|
| **Sound-crossing** | t_cs = R / c_s, where c_s = √(γ kT / μ m_p) | Birzan et al. (2004), Eq. 2 |
| **Buoyancy rise** | t_buoy = R √(S C_D / 2gV), with C_D = 0.75 | Churazov et al. (2001); Birzan et al. (2004), Eq. 3 |
| **Refill** | t_refill = 2 √(r / g) | McNamara & Nulsen (2007) |

---

## Key References

| Paper | Topic | Used For |
|-------|-------|----------|
| McDonald et al. (2012), *Nature* 488, 349 | Discovery paper | Cooling luminosity, SFR, redshift |
| McDonald et al. (2015), *ApJ* 811, 111 | Cavity detection | Cavity radius, distance, age |
| McDonald et al. (2019), *ApJ* 885, 63 | Deprojected profiles | Temperature, density profiles |
| Hlavacek-Larrondo et al. (2015), *ApJ* 805, 35 | SPT cavity survey | Ellipsoidal dimensions |
| Russell et al. (2017), *ApJ* 836, 130 | ALMA molecular gas | Physical context |
| Birzan et al. (2004), *ApJ* 607, 800 | Cavity methodology | Age formulas, ellipsoidal convention |
| McNamara & Nulsen (2007), *ARA&A* 45, 117 | AGN feedback review | Refill time formula |
| Churazov et al. (2001), *ApJ* 554, 261 | Buoyancy formulation | Drag coefficient |

---

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{phoenix_feedback_2026,
  title  = {Phoenix Cluster AGN Feedback Pipeline},
  author = {Phanendra},
  url    = {https://github.com/phanendra09/phoenix-cluster-feedback},
  year   = {2026}
}
```

## License

This project is licensed under the [MIT License](LICENSE).
