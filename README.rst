=======
SunCASA
=======
|Latest Version|

.. |Latest Version| image:: https://img.shields.io/pypi/v/suncasa.svg
   :target: https://pypi.python.org/pypi/suncasa/

SunCASA is an open-source CASA-based Python package for reducing, analyzing, and visualizing solar dynamic spectroscopic
imaging data at radio wavelengths. Our homepage `SunCASA`_ has more information about the project.

.. _SunCASA: https://github.com/suncasa/suncasa


Installation
============
Currently SunCASA requires an installation of `CASA`_ (Common Astronomy Software Applications) to
function. The latter is the arguably most advanced general-purpose software for the new generation of radio
interferometers including `VLA`_, `ALMA`_, and `EOVSA`_. However CASA is platform dependent, and is not available on Windows. The
stable version of SunCASA can be installed via pip, along with CASA 6 for MacOS and certain Linux distros.

The recommended way to install SunCASA is with `pip`_.
Once pip is installed, run the following command:

.. code:: bash

    $ pip install suncasa

For detailed installation instructions, see the `installation guide`_ in the EOVSA wiki.

.. _VLA: http://www.vla.nrao.edu/
.. _ALMA: https://almascience.nrao.edu/
.. _EOVSA: http://www.ovsa.njit.edu/
.. _CASA: https://casa.nrao.edu/
.. _pip: https://packaging.python.org/tutorials/installing-packages/
.. _installation guide: http://www.ovsa.njit.edu/wiki/index.php/SunCASA_Installation

Usage
=====
Here is a quick `example`_ for using SunCASA to reduce and visualize dynamic spectroscopic imaging data obtained from the
Expanded Owens Valley Solar Array (EOVSA). The same procedure has been tested to work
on data from the Karl G. Jansky Very Large Array (VLA). More example will be added in the near future.

.. _example: https://github.com/suncasa/suncasa-src/blob/master/examples/EOVSA_tutorial_RHESSI2021.ipynb



