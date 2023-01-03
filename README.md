![license](https://img.shields.io/badge/Licence-GPL--3-blue.svg) 

This package was originally forked from [`lidRplugins`](https://github.com/Jean-Romain/lidRplugins), and extends the `li2012()` tree segmentation function to adjust the constraints of the `R` parameter contained within the main [`lidR`](https://github.com/r-lidar/lidR) package. Thank you Jean-Romain for both of these packages and for your valuable input!

Here, I introduce a new idea for being able to adjust the `R` value of the `li2012()` tree segmentation algorithm in a couple of new ways:

1. Specifying `NULL` --> This will enable automatic window sizes to be created for each point using the `lmfxauto()` algorithm. These values get passed to the into the `filter_local_maxima` C++ function to find the local maximas;
2. Specifying a function for defining window sizes. This functionality is introduced in the original `lmf()` function from the `lidR` package to adjust the `ws` parameter;
3. Specifying a numeric value. This reverts to the default `li2012()` function usage.

Other algorithms, functions, and package exports from the `lidRplugins` package have been stripped away here so that the new `li2012_auto()` function may be focussed on. Testers of this implementation are welcome and feedback would be appreciated!  

### Parameter free tree detection :microscope:

`lmfxauto()` is a fast algorithm for individual tree detection with 0 parameters designed to process thousands of square kilometres without supervision. It is based on [`lmfauto()`](https://github.com/Jean-Romain/lidRplugins/blob/master/R/algo-lmfauto.R) from the `lidRplugins` package where the first step is changed to use `lmfx(5)` instead of `lmf(5)`.

## Installation

```r
remotes::install_github("mcoghill/lidR.li2012enhancement")
```

To install the package from github make sure you have a working development environment.

* **Windows**: Install [Rtools.exe](https://cran.r-project.org/bin/windows/Rtools/).  
* **Mac**: Install `Xcode` from the Mac App Store.
* **Linux**: Install the R development package, usually called `r-devel` or `r-base-dev`
