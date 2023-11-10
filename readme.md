## Repository for time delay of arrival type III burst implementation

The goal of this notebook is to demonstrate implementation of the Time delay of arrival (TDOA) method of triangulating the position of radio transients using public data only. 

This has recently been made much more tractable by the efforts of Krupar et al. ([2022](https://doi.org/10.25935/4tak-5225))  in archiving the STEREO native cadence radio spectrograms as CDF files at ([STEREO A](https://cdaweb.gsfc.nasa.gov/pub/data/stereo/ahead/l3/waves/) [STEREO B](https://cdaweb.gsfc.nasa.gov/pub/data/stereo/behind/l3/waves/) ).

The methodology is described in Badman+[2022](https://ui.adsabs.harvard.edu/abs/2022ApJ...938...95B/abstract). At its core, the procedure is quite simple. Assuming we have three spacecraft (e.g. A, B and C) located approximately co-planar (e.g. in the ecliptic plane) at 2D positions $(x,y)_{i={A,B,C}}$ which all receive the same radio signal at a given frequency $f$ at respective times $t_{i={A,B,C}}$, we can map these positions (assuming the radio source also is approximately co-planar) to a source origin position $(x_s(f),y_s(f))$ via an analytic construction. This construction essentially corresponds to intersecting two hyperbolae where the locations of the geometric foci are set by the positions of pairs of the spacecraft, and the hyperbola curvature are set by the relative time __delay__ of the measurement between those same two spacecraft  (see figure 11 in Badman+[2022](https://ui.adsabs.harvard.edu/abs/2022ApJ...938...95B/abstract) ).

To go from a real set of radio observations to a radio burst trajectory, there are four steps :

* Download radio spectrogram data and slice to isolate the feature of interest in all spacecraft
* Produce the position vectors of the individual spacecraft at the time of interest.
* Process the spectrograms to perform feature extraction of the feature of interest, resulting in a 1d curve f=f(t).
* Use the extracted feature from all three spacecraft and iterate the analytical procedure as function fo frequency to derive the source position as a function of frequency.

__Citations__

* Krupar, V., Q.N. Nguyen, X. Bonnin, B. Cecconi & M. Maksimovic (2022). STEREO/Waves/LFR-HFR L3 DF Data Collection (Version 1.0) [Data set]. PADC. https://doi.org/10.25935/4TAK-5225

* Samuel T. Badman et al 2022 ApJ 938 95 https://doi.org/10.3847/1538-4357/ac90c2

__Installation Instructions__

This repository requires a submodule (a fork of pyspedas : https://github.com/STBadman/pyspedas/tree/583fe444242c4e64923295579a476372dc3b0a9c which adds access to the cited STEREO dataset). To setup this submodule, in bash do : 

``` 
git submodule init
git submodule update
```

The repository has been tested using the bundled conda environment file. To create and use this :

```
cd <path/to/>Radio-Public
conda env create --file conda_env.yml
conda activate tdoa
```


