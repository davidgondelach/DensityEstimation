# DensityEstimation
Matlab code for global thermospheric density estimation using two-line element data.

This is a complete toolbox for Matlab that enables you to estimate the global thermospheric density using two-line element data. Three different reduced-order density models can be employed for the estimation. Details of the technique and models can be found in the journal paper, see https://doi.org/10.1029/2019SW002356 or https://arxiv.org/abs/1910.00695.


Copyright © 2020 by David Gondelach and Richard Linares

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3634245.svg)](https://doi.org/10.5281/zenodo.3634245)


### License
This code is licensed under the GNU General Public License version 3 - see the [LICENSE](LICENSE) file for details.


### Acknowledgments
The contributions by Dr. Piyush M. Mehta in the design and implementation of the code are acknowledged.

The MATLAB code for Jacchia-Bowman 2008 model was developed by Meysam Mahooti (copyright 2018) and was downloaded from https://www.mathworks.com/matlabcentral/fileexchange/56163-jacchia-bowman-atmospheric-density-model (version 2.0.0.0).

The MATLAB code for the SGP4 model and several time and reference frame routines was developed by David Vallado (and others) and was downloaded from https://celestrak.com/software/vallado-sw.php.


### References
The density modeling and estimation techniques are described in:
```
@article{gondelach2019realtime,
  author = {Gondelach, David J. and Linares, Richard},
  title = {Real-Time Thermospheric Density Estimation Via Two-Line-Element Data Assimilation},
  journal = {Space Weather},
  doi = {10.1029/2019SW002356},
  url = {https://doi.org/10.1029/2019SW002356}
}
```
see https://doi.org/10.1029/2019SW002356 or https://arxiv.org/abs/1910.00695.

### Installation instructions
1. Download the DensityEstimation Matlab code
2. Download and install SPICE Toolkit for Matlab: https://naif.jpl.nasa.gov/naif/toolkit_MATLAB.html
3. Set the path to the SPICE Toolkit directory in mainDensityEstimation.m
4. Download SPICE kernels (i.e. ephemeris files) from https://naif.jpl.nasa.gov/pub/naif/generic_kernels/ and put them in the folder Data. See links below.
5. Download space weather file from Celestrak and put in folder Data: https://www.celestrak.com/SpaceData/SW-All.txt
6. Download Earth orientation data file from Celestrak and put in folder Data: https://www.celestrak.com/SpaceData/EOP-All.txt
7. Download 2 space weather files needed for the JB2008 model and put in folder Data: http://sol.spacenvironment.net/jb2008/indices/SOLFSMY.TXT  and  http://sol.spacenvironment.net/jb2008/indices/DTCFILE.TXT 
8. If you want to use automatic download of TLE data, then specify your “www.space-track.org” username and password in runDensityEstimationTLE.m, line 95-96. (Alternatively, you can download the TLE data manually and put them in the folder TLEdata with naming convention: [NORADID].tle, e.g. 12388.tle )
9. For each object used for estimation, specify the ballistic coefficient (BC) in the text file: Data/BCdata.txt


### Run instructions
1. Open mainDensityEstimation.m
2. Specify the time window for density estimation by setting the start year, month and day and number of days
3. (Optionally) select the reduced-order density model and the reduction order (default order: r=10)
4. (Optionally) select the objects to use for density estimation
5. Run mainDensityEstimation


### Ephemeris files
Download the following ephemeris files and put them in the Data folder:
* latest_leapseconds.tls:  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/
* de430.bsp:  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/
* earthstns_itrf93_050714.bsp:  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/stations/
* pck00010.tpc, earth_fixed.tf, earth_070425_370426_predict.bpc, earth_720101_070426.bpc, earth_latest_high_prec.bpc:  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/


### Technical notes
This version of the toolbox does not include third-body and solar radiation pressure perturbations. The speed of the code has not been optimized; 10-day density estimation may required several hours of computation. The estimation can be speed up by reducing the degree of the gravity model (see mainDensityEstimation.m, line 79) at the cost of reduced density estimate accuracy.

MATLAB R2018b (Version 9.5) was used to develop the code.


### Version
The latest version is this toolbox may be found on: https://github.com/davidgondelach/DensityEstimation


David Gondelach, Jan 2020
email: davidgondelach@gmail.com
