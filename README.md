# suncasa

Required packages

- CASA 5.4.0, which comes with
  * NumPy 1.14.3
  * Matplotlib 1.1.0
  
Additionally, we need
- SunPy 0.9.3
- AstroPy 2.0.12

Currently SunCASA requires an installation of CASA (Common Astronomy Software Applications; https://casa.nrao.edu/) to function. The latter is the arguably the most advanced general-purpose software for the new generation of radio interferometers including VLA, ALMA, and EOVSA. However CASA is platform dependent, and is not available on Windows. We provide packages of SunCASA for MacOS and certain Linux distros that already includes CASA 5.4.0. Please refer to this page for detailed instructions on installing SunCASA: http://www.ovsa.njit.edu/wiki/index.php/SunCASA_Installation 

Some latest examples for using SunCASA to reduce and visualize dynamic spectroscopic imaging data obtained from the Expanded Owens Valley Solar Array (EOVSA) are available at this link: http://www.ovsa.njit.edu/wiki/index.php/EOVSA_Data_Analysis_Tutorial. The same procedure has been tested to work on data from the Karl G. Jansky Very Large Array (VLA). More examples will be added in the near future.
