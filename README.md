# alphahull: Generalization of the Convex Hull of a Sample of Points in the Plane

[![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![R build
status](https://github.com/beatrizpateiro/alphahull/workflows/R-CMD-check/badge.svg)](https://github.com/beatrizpateiro/alphahull/actions)
[![](https://www.r-pkg.org/badges/version/alphahull?color=green)](https://cran.r-project.org/package=alphahull)
[![](http://cranlogs.r-pkg.org/badges/grand-total/alphahull?color=green)](https://cran.r-project.org/package=alphahull)
[![](http://cranlogs.r-pkg.org/badges/last-month/alphahull?color=green)](https://cran.r-project.org/package=alphahull)


## Overview

Computation of the alpha-shape and alpha-convex hull of a given sample of points in the plane. The concepts of alpha-shape and alpha-convex hull generalize the definition of the convex hull of a finite set of points. The programming is based on the duality between the Voronoi diagram and Delaunay triangulation. The package also includes a function that returns the Delaunay mesh of a given sample of points and its dual Voronoi diagram in one single object.

## Installation

Get the released version from CRAN:

``` r
# Install the package
install.packages("alphahull")

# Load package
library(alphahull)
```

Alternatively, get the latest version from GitHub:

``` r
# Install the package
library(devtools)
install_github("beatriz.pateiro/alphahull")

# Load package
library(alphahull)
```
