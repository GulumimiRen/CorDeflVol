# CorDeflVol (Corneal Deflection Volume) - R Code

This repository provides an R script to compute **DeflVol**, the estimated **corneal deflection volume** at the **highest concavity (HC)** moment in the Corvis ST framework.

The implementation uses a **two-sphere spherical geometry model** and provides two mathematically equivalent algorithms:

- **Disk integration**: numerical trapezoidal integration of the rotated sagittal overlap region.
- **Closed-form intersection volume**: the analytical volume of intersection of two spheres.

Both values are computed for every row. The default reported result is written to `df$DeflVol`, which equals the closed-form result `df$DeflVol_closed`.

---

## Repository contents

- `CorDeflVol.R` - main R script for DeflVol calculation

---

## Requirements

Install the required R packages before running the script:

```r
install.packages("readr")
```

---

## Input data

Prepare a CSV file named `data.csv` in the working directory with the following columns:

- `CR` - Corneal Curvature Radius of the anterior surface
- `IR` - Inverse Curvature Radius
- `DeflAmp` - Deflection Amplitude at the HC moment

### Important notes

- Units must be consistent across `CR`, `IR`, and `DeflAmp` (for example, all in `mm`).
- The script expects the exact column names above.
- If any of the three input values is missing, the output for that row will be `NA`.
- If `CR <= 0` or `IR <= 0`, the output for that row will be `NA`.

If the input units are in `mm`, the resulting `DeflVol` will be in `mm^3`.

---

## Geometry convention

The script follows these deterministic geometry settings:

- Sphere-center distance:

  `d = CR + IR - DeflAmp`

- Coordinate convention:

  - `z_inverse = 0`
  - `z_cornea = d`

- Integration region for disk integration:

  The valid interval is determined from the overlap of the two spheres along the `z` axis.

- Cross-sectional radius:

  At each `z`, the effective radius is the minimum of the two sphere cross-sectional radii over the overlap region.

- Numerical disk integration:

  `DeflVol = pi * integral(rho(z)^2 dz)`

  computed numerically using the trapezoidal rule with default resolution `N = 4000`.

- Closed-form two-sphere intersection:

  `DeflVol = pi * (CR + IR - d)^2 * [d^2 + 2d(CR + IR) - 3(CR - IR)^2] / (12d)`

---

## How to run

1. Put `data.csv` in the same directory as `CorDeflVol.R`.
2. Open R or RStudio in that directory.
3. Run:

```r
source("CorDeflVol.R")
```

After execution, the script will:

- read `data.csv` into `df`
- compute disk-integration and closed-form DeflVol, plus their absolute/relative difference

If `data.csv` is not present, the script still loads the calculation functions and skips creating `df`.

The script also includes an optional output line:

```r
# write_csv(df, "data_with_DeflVol.csv")
```

You can uncomment it if you want to save the results as a new CSV file.

---

## Output

The computed values are stored in:

```r
df$DeflVol_integral   # numerical disk integration, N = 4000
df$DeflVol_closed     # closed-form two-sphere intersection volume
df$DeflVol            # default reported result, same as DeflVol_closed
df$DeflVol_abs_error       # absolute difference between the two algorithms
df$DeflVol_rel_error       # relative difference between the two algorithms
```

---

## Function usage

You can also source the file and call the calculation functions directly:

```r
source("CorDeflVol.R")

# One eye, closed-form algorithm
compute_deflvol_one(CR = 7.8, IR = 6.4, DeflAmp = 0.93, method = "closed")

# One eye, disk-integration algorithm
compute_deflvol_one(CR = 7.8, IR = 6.4, DeflAmp = 0.93, method = "integral", N = 4000)

# Data frame input, both algorithms
df_out <- compute_deflvol_dual(df, N = 4000)
```

---

## Notes

- The implementation uses numerical safety handling by clamping negative squared-radius terms to zero before square root evaluation.
- If the computed overlap interval is invalid or empty, the function returns `0` for that row.
- In physiological HC geometry, disk integration and the closed-form solution should agree up to numerical quadrature error. Increasing `N` reduces the disk-integration error.
