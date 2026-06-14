# CorDeflVol

Compute **corneal deflection volume** (DeflVol) at the **highest concavity (HC)**
moment using a **two-sphere intersection model**. The default algorithm is a
closed-form analytical volume; an equivalent disk-integration method is also
available.

## Install

```r
# from GitHub
remotes::install_github("GulumimiRen/CorDeflVol")

# from the package directory
install.packages(".", repos = NULL, type = "source")
```

### Install from release tarball

Download [CorDeflVol_0.1.0.tar.gz](https://github.com/GulumimiRen/CorDeflVol/releases/download/v0.1.0/CorDeflVol_0.1.0.tar.gz)
from the [v0.1.0 release](https://github.com/GulumimiRen/CorDeflVol/releases/tag/v0.1.0)
(Assets on that page). The file is not stored in the git repository.

```r
install.packages("path/to/CorDeflVol_0.1.0.tar.gz", repos = NULL, type = "source")
```

Or download in R before installing:

```r
url <- "https://github.com/GulumimiRen/CorDeflVol/releases/download/v0.1.0/CorDeflVol_0.1.0.tar.gz"
dest <- tempfile(fileext = ".tar.gz")
download.file(url, dest, mode = "wb")
install.packages(dest, repos = NULL, type = "source")
```

## Usage

```r
library(CorDeflVol)

# scalar
deflvol(7.8, 6.4, 0.93)

# vectorized
deflvol(c(7.8, 7.5), c(6.4, 6.2), c(0.93, 0.90))

# column-wise with dplyr
library(dplyr)
df <- df |> mutate(DeflVol = deflvol(CR, IR, DeflAmp))

# disk integration (optional)
deflvol(7.8, 6.4, 0.93, method = "integral", N = 4000L)
```

## Inputs

| Parameter   | Description                                      |
|-------------|--------------------------------------------------|
| `CR`        | Corneal curvature radius                         |
| `IR`        | Inverse curvature radius                         |
| `defl_amp`  | Deflection amplitude at HC (same as `DeflAmp`)   |

Units must be consistent across all three inputs (e.g. mm). When inputs are in
mm, DeflVol is returned in μL (1 mm³ = 1 μL). Missing or non-positive radii
return `NA`; non-overlapping geometry returns 0.

## Geometry

Sphere-center distance: `d = CR + IR - defl_amp`.

Closed-form two-sphere intersection volume:

`DeflVol = π (CR + IR − d)² [d² + 2d(CR + IR) − 3(CR − IR)²] / (12d)`

The integral method computes the same volume by trapezoidal disk integration
over the overlap region (default `N = 4000` segments).

## License

MIT — see [LICENSE](LICENSE).
