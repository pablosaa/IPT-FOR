# IPT-FOR
## Integrated Profiling Technique for Atmospheric Parameters

IPT is a method to derive physically consistent profiles of temperature, humidity, and cloud liquid water content.
IPT has been developed at the _Meteorological Institute, University of Bonn, Germany_ by **Dr. Ulrich Löhnert** (see reference below) and in 2005 I have adapted the code from IDL to FORTRAN 90 mainly for improving speed performance and to have alternative code able to be run within non-proprietary software tools.

The IPT can be used by combining several remote sensing and in-situ instruments, for instance: a ground-based multichannel microwave radiometer, a cloud radar, a lidar-ceilometer, radiosonde measurement, and ground-level measurements of standard meteorological properties with statistics derived from results of a microphysical cloud model
.
All measurements are integrated within the framework of optimal estimation to guarantee a retrieved profile with maximum information content. It has been found that the liquid water content profiles obtained with the IPT are significantly more accurate than common methods that use the microwave-derived liquid water path to scale the radar reflectivity profile. 

### Structure ###
/IPT-FOR
        |-- /source  : Include all Fortran 90 source codes.
        |-- /bin     : Is where executable file will be created after running the Makefile inside this directory.
        |-- /databank: This directory should include all database needed as, for instance:
            |--/O2.dat
            |--/coeff.dat
            |--/oxygen_l93.dat
            |--/water_l93.dat
            |--/H2O.dat
            |--/coeff_fap_r98_r.dat
            etc.
        |-- /output   : In this directory will be placed all outputs from IPT


### Refence ###
Löhnert, U., S. Crewell, and C. Simmer, 2004: _"An Integrated Approach toward Retrieving Physically Consistent Profiles of Temperature, Humidity, and Cloud Liquid Water"_. J. Appl. Meteor., 43, 1295–1307, [doi:10.1175/1520-0450](https://doi.org/10.1175/1520-0450(2004)043<1295:AIATRP>2.0.CO;2)
