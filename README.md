# CorDeflVol

Compute **corneal deflection volume** (DeflVol) at the **highest concavity (HC)**
moment using a **two-sphere intersection model**. The default algorithm is a
closed-form analytical volume; an equivalent disk-integration method is also
available.

## Install

```r
install.packages("pak")
pak::pak("GulumimiRen/CorDeflVol")
```

From the package directory:

```r
install.packages(".", repos = NULL, type = "source")
```

### Install from release tarball

Download [CorDeflVol_0.1.0.tar.gz](https://github.com/GulumimiRen/CorDeflVol/releases/download/v0.1.0/CorDeflVol_0.1.0.tar.gz)
from the [v0.1.0 release](https://github.com/GulumimiRen/CorDeflVol/releases/tag/v0.1.0)
(Assets on that page).

With **pak** (local path or release URL):

```r
pak::pak("path/to/CorDeflVol_0.1.0.tar.gz")
```

```r
pak::pak("https://github.com/GulumimiRen/CorDeflVol/releases/download/v0.1.0/CorDeflVol_0.1.0.tar.gz")
```

With base R:

```r
install.packages("path/to/CorDeflVol_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Usage

```r
library(CorDeflVol)

CR <- 7.8
IR <- 6.4
DeflAmp <- 0.93

# scalar
deflvol(CR, IR, DeflAmp)

# vectorized
deflvol(c(7.8, 7.5), c(6.4, 6.2), c(0.93, 0.90))

# column-wise with dplyr
library(dplyr)
df <- df |> mutate(DeflVol = deflvol(CR, IR, DeflAmp))

# disk integration (optional)
deflvol(CR, IR, DeflAmp, method = "integral", N = 4000L)
```

## Inputs

| Parameter  | Description                                    | Unit |
|------------|------------------------------------------------|------|
| `CR`       | Corneal curvature radius                       | mm   |
| `IR`       | Inverse curvature radius                       | mm   |
| `DeflAmp`  | Deflection amplitude at HC                     | mm   |

Units must be consistent across all three inputs (mm).

Missing or non-positive radii return `NA`; non-overlapping geometry returns 0.

## Outputs

| Parameter | Description              | Unit |
|-----------|--------------------------|------|
| `DeflVol` | HC deflection volume     | μL   |

When inputs are in mm, `DeflVol` is in μL (1 mm³ = 1 μL).

## Geometry

Sphere-center distance: `d = CR + IR - DeflAmp`.

Closed-form two-sphere intersection volume:

`DeflVol = π (CR + IR − d)² [d² + 2d(CR + IR) − 3(CR − IR)²] / (12d)`

The integral method computes the same volume by trapezoidal disk integration
over the overlap region (default `N = 4000` segments).

## Citation

If you use CorDeflVol in your research, please cite the article describing the
DeflVol model:

> Ren Q, Deng Z, Chu Z, Pan Z and Wan P (2026). DeflVol: a model-derived volumetric
> descriptor of corneal deflection at highest concavity under Corvis ST air-puff
> loading. *Frontiers in Bioengineering and Biotechnology*, 14:1933733.
> https://doi.org/10.3389/fbioe.2026.1933733

After installing the package, retrieve the reference and its BibTeX entry in R:

```r
citation("CorDeflVol")
toBibtex(citation("CorDeflVol"))
```

```bibtex
@Article{Ren2026DeflVol,
  author  = {Qi Ren and Zhaozhe Deng and Zhe Chu and Zeyu Pan and Pengxia Wan},
  title   = {{DeflVol}: a model-derived volumetric descriptor of corneal deflection at highest concavity under {Corvis ST} air-puff loading},
  journal = {Frontiers in Bioengineering and Biotechnology},
  year    = {2026},
  volume  = {14},
  pages   = {1933733},
  doi     = {10.3389/fbioe.2026.1933733}
}
```

## License

MIT — see [LICENSE](LICENSE).
