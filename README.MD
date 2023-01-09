# suncasa

[![pypi](https://img.shields.io/pypi/v/suncasa.svg)](https://pypi.python.org/pypi/suncasa/)

[suncasa](https://github.com/suncasa/suncasa-src) is a Python wrapper around [CASA](https://casa.nrao.edu/) for
importing, calibrating, synthesis imaging and visualizing solar spectral imaging data. CASA is one of the leading
software tools for "supporting the data post-processing needs of the next generation of radio astronomical telescopes
such as ALMA and VLA", an international effort led by
the [National Radio Astronomy Observatory (NRAO)](https://public.nrao.edu/). The current version of CASA uses Python 3.6
to 3.8 on most of the platform. More information about the compatibility of CASA can be found
on [CASA Compatibility](https://casadocs.readthedocs.io/en/latest/notebooks/introduction.html#Compatibility/).
> **_NOTE:_** CASA is available ONLY on unix-based platforms, and therefore, so as for suncasa.

In principle, suncasa an casa can be installed easily using [pip](https://pypi.org/project/suncasa/), assuming the right python version is in use. However, YMMV if your do not have prerequisite environment configured.The followinig tutorial demonstrates the installation of **suncasa** in **python 3.8** environment on **MacOS 12** and **Ubuntu 18.04.6 LTS**.        

[comment]: <> (## Table of Contents)

[comment]: <> (* **[Installation]&#40;#installation&#41;**)

[comment]: <> (    * [Prerequisite OS Libraries]&#40;#Prerequisite-OS-Libraries&#41;)




[comment]: <> (## Table of Contents)

[comment]: <> (* **[Installation of the monolithic CASA]&#40;#installation-of-the-monolithic-casa&#41;**)

[comment]: <> (  * [a]&#40;a&#41;)

[comment]: <> (* **[Requirements]&#40;#requirements&#41;**)

[comment]: <> (* **[Installation]&#40;#installation&#41;**)

Full installation instructions of the standalone CASA package (as a Linux (.tar) or Mac (.dmg) files) including the build-in Python environment are available from
NRAO [Downloads](http://casa.nrao.edu/casa_obtaining.shtml)
page ([http://casa.nrao.edu/casa_obtaining.shtml](http://casa.nrao.edu/casa_obtaining.shtml))

The CASA 6.x series is also available as modular packages, giving users the flexibility to build CASA tools and tasks in
their own Python environment. This includes the casatools, casatasks, and casampi modules, allowing for core data
processing capabilities in parallel. 

<hr />

**If you want to use full funtionalities provided by the standalone CASA, please go through all steps.** 

**If you just want to install SUNCASA with the modular CASA option in the usual "Pythonic" manner, after Step 1 (Python installation), you can skip Step 2 and go directly to Step 3.**

<hr />

## Step 1: Python 3.8 Installation[](#installation) 
(Adapted from [this tutorial](https://casadocs.readthedocs.io/en/latest/notebooks/introduction.html#id1))

### Prerequisite OS Libraries[](#Prerequisite-OS-Libraries)

CASA requires certain libraries be installed in the users operating system. Some may already be present by default. In
case they are not, the following list should be checked before using CASA or if errors are encountered at runtime.
Commands and package names are for Red Hat Linux, but equivalents can be found for other Linux distributions.

**On Linux:**

A system administrator may be required to install OS libraries. One can
use [yum](https://access.redhat.com/documentation/en-us/red_hat_enterprise_linux/7/html/system_administrators_guide/ch-yum) (
RHEL) or [apt-get](https://help.ubuntu.com/community/AptGet/Howto) (Ubuntu) to install the missing libraries. For
example, for RHEL8 the following OS libraries should be installed:

```bash
sudo yum install libgfortran3
sudo yum install ImageMagick*
sudo yum install xorg-x11-server-Xvfb
sudo yum install compat-libgfortran-48
sudo yum install libnsl
sudo yum install libcanberra-gtk2
```

For modular CASA, one must supply their own Python environment. There are many, including ipython and Jupyter, here is a
basic example for
RHEL8 ([yum](https://access.redhat.com/documentation/en-us/red_hat_enterprise_linux/7/html/system_administrators_guide/ch-yum))
and Ubuntu ([apt-get](https://help.ubuntu.com/community/AptGet/Howto)), respectively.

```bash
sudo yum install python38-devel
```

```bash
sudo apt install python3.8 python3.8-venv python3-venv
```

**On Macintosh:**

Install [XQuartz](https://www.xquartz.org/)
and [Xcode](https://apps.apple.com/us/app/xcode/id497799835?mt=12). 

[comment]: <> (gcc &#40;[brew]&#40;https://brew.sh/&#41; install&#41;,)

Same as for Linux, one must supply their own Python environment for modular CASA. Here is an example for MacOS to install python 3.8 using PyEnv.

Clone pyenv to a local directory (e.g., `$HOME/.pyenv`):
```bash
git clone https://github.com/pyenv/pyenv.git $HOME/.pyenv
```
Set up your shell environment for Pyenv. Assuming you are using **bash**, add the configuration commands to `~/.bashrc` by running the following in your terminal:
~~~bash
bash
echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bashrc
echo 'command -v pyenv >/dev/null || export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bashrc
echo 'eval "$(pyenv init -)"' >> ~/.bashrc
~~~
Install the Python 3.8.

```bash
pyenv install 3.8.16
```
> **_Troubleshooting:_** If you encounter an error that "C compiler cannot create executables" then the simplest way to solve this is to reinstall Apple's [Xcode](https://apps.apple.com/us/app/xcode/id497799835?mt=12). Xcode is a tool created by Apple that includes all the C libraries and other tools that Python uses when it runs on MacOS. Xcode is a whopping 11 gigabytes, but you'll want to be up-to-date. You may want to run this while you're sleeping. Once it is done, just re-run the above `pyenv install 3.8.16` and it should now work.

Now you can set this version of Python to be the default version of Python MacOS uses.
```bash
pyenv global 3.8.16
```
[comment]: <> (When using the casampi package from modular CASA, additional MPI libraries are needed:)

[comment]: <> (*   `$: sudo yum install openmpi-devel`)

[comment]: <> (*   `$: sudo yum install mpich-devel`)

[comment]: <> (*   User installs mpi4py with:)

[comment]: <> (    *   `$: env MPICC=/usr/lib64/openmpi/bin/mpicc python -m pip install mpi4py --no-cache-dir`)

[comment]: <> (    *   Ensure mpirun is found: `which mpirun`)

[comment]: <> (    *   If not, try full path: `export PATH=/usr/lib64/openmpi/bin/mpirun\$PATH`)


[comment]: <> (Alternative method for _NRAO systems only_)

[comment]: <> (*   contact the helpdesk to install casa-toolset-3 &#40;which contains the previous libraries&#41;)

[comment]: <> (*   the run `export PATH=/opt/casa/03/bin:\$PATH`)

### Step 2: CASA (Standalone Monolithic Distribution) Installation [](#Monolithic-CASA-Distribution)

Download the standalone CASA package (as a Linux (.tar) or Mac (.dmg) files) including the build-in Python environment following
[http://casa.nrao.edu/casa_obtaining.shtml](http://casa.nrao.edu/casa_obtaining.shtml).

**On Linux:**

1. Download the .tar file and place it in a work directory (e.g. `$HOME/casa`).

2. From a Linux terminal window:
   ~~~ bash
   tar -xvf casa-xyz.tar.xz
   $HOME/casa/casa-xyz/bin/casa
   ~~~ 
3. Add casa to environment `PATH`. Assuming you are using *bash*, add the configuration commands to `~/.bashrc` by running the following in your terminal:

   ~~~ bash
   echo 'export CASAROOT="$HOME/casa/casa-xyz/"' >> ~/.bashrc
   echo 'export PATH="$PATH:$CASAROOT/bin"' >> ~/.bashrc
   ~~~

[comment]: <> (   Then, if you have `~/.profile`, `~/.bash_profile` or `~/.bash_login`, add the commands there as well. If you have)

[comment]: <> (   none of these, add them to `~/.profile`.)

[comment]: <> (    * to add to `~/.profile`:)

[comment]: <> (   ~~~ bash)

[comment]: <> (   echo 'export CASAROOT="$HOME/casa/casa-xyz/"' >> ~/.profile)

[comment]: <> (   echo 'export PATH="$PATH:$CASAROOT/bin"' >> ~/.profile)

[comment]: <> (   ~~~)

[comment]: <> (    * to add to `~/.bash_profile`:)

[comment]: <> (   ~~~ bash)

[comment]: <> (   echo 'export CASAROOT="$HOME/casa/casa-xyz/"' >> ~/.bash_profile)

[comment]: <> (   echo 'export PATH="$PATH:$CASAROOT/bin"' >> ~/.bash_profile)

[comment]: <> (   ~~~)

4. Restart your shell for the `PATH` changes to take effect.

  ```sh
  exec "$SHELL"
  ```    

You can now start casa.

  ```sh
  casa
  ```  

The one caveat is that CASA on Linux currently will not run if the Security-Enhanced Linux option of the linux operating
system is set to enforcing. For the non-root install to work, SElinux must be set to disabled or permissive (in _
/etc/selinux/config_) or you must run (as root): `$: setsebool -P allow_execheap=1`. Otherwise, you will encounter
errors like:

error while loading shared libraries: /opt/casa/casa-20.0.5653-001/lib/liblapack.so.3.1.1: cannot restore segment prot
after reloc: Permission denied_

> **_WARNING:_** By default, python 3.6 (and earlier versions of python 3) include the current working directory in the python path at startup. Any script in that directory with the same name as a standard python module or a CASA module will be detected and used by python instead of the code that is delivered with CASA. Protections have been included for files called “new.py” and “pickle.py”, but other scripts may cause problems with the CASA startup. For example, do not include a file named runpy.py in the working directory.

**On Macintosh:**

1. Download the .dmg disk image file

2. Double click on the disk image file (if your browser does not automatically open it).

3. Drag the CASA application to the _Applications_ folder of your hard disk.

4. Eject the CASA disk image.

5. Double click the CASA application to run it for the first time. If the OS does not allow you to install apps from
   non-Apple sources, please Change the settings in “System Preferences-> Security & Privacy -> General” and “Allow
   applications downloaded from: Mac App store and identified developers”.

6. Optional: Create symbolic links to the CASA version and its executables (Administrator privileges are required),
   which will allow you to run `casa`, `casaviewer`, `casaplotms`, etc. from any terminal command line. To do so, run

```console
 !create-symlinks
```   

### Step 3: Install SUNCASA (along with CASA modular version) [](#Modular-Packages)

Pip wheels for casatools, casatasks, as well as suncasa are available as Python 3 modules. This allows simple installation and import into standard Python environments. Make sure you have set up your machine with the _necessary prerequisite libraries_ first. Then a separate installation of desired modules (from a unix terminal window) as follows. First create a Python virtual environment
named `suncasaenv` under the `$HOME` directory:

```bash
cd $HOME
python3.8 -m venv suncasaenv
```

> **_NOTE:_**  We strongly recommend using a Python virtual environment to prevent breaking any packages within a pre-existing Python environment.

Then use [pip](https://pypi.org/project/suncasa/) to install suncasa within the newly-created virtual environment:

```bash
source suncasaenv/bin/activate
pip install --upgrade pip
pip install suncasa
```

> **_NOTE:_**  If this does not work, it could be due to unsuccessful installation of some dependencies. Running these commands should address this.

```bash
pip install casatasks
pip install casatools
pip install casadata
pip install PyQt5
pip install sunpy[all]
pip install suncasa
```

To exit the python venv, type `deactivate` from the terminal. However, the rest of this tutorial **assumes the venv is
active** (to reactivate, type `source $HOME/suncasaenv/bin/activate`)

You can update suncasa to its latest version by running:

```bash
pip install --upgrade suncasa
```

**Sanity check**

With the pip installation, suncasa as well as CASA may be used in a standard Pythonic manner. For example, suncasa
modules and CASA tasks can be invoked using “import”, while CASA tools first need to be instantiated:

```console
(suncasaenv) $ ipython
In [1]: import suncasa
In [2]: help(suncasa)
In [3]: import casatasks
In [4]: help(casatasks)
```

The use of python3 venv is a simple built-in method of containerizing the pip install such that multiple versions of suncasa can be kept on a single machine in different environments. In addition, suncasa is built and tested using
standard python 3.8 libraries (backward compatible) which can be replicated with a fresh venv, keeping the libraries needed for suncasa
isolated from other libraries which may already be installed on your machine.

**Add EOVSA to CASA observatories list**

The following command creates a configuration file to define the EOVSA instrument so that CASA is aware of
it. We plan to add EOVSA officially to the CASA distribution in the future, after which this step is unnecessary.

```bash
printf "import sys \nimport os \nimport sysconfig \nimport casatools \nimport casadata \nimport time \nlogfile='casalog-{}.log'.format(time.strftime('%%Y%%m%%d-%%H',time.localtime())) \ntelemetry_enabled = False \ncrashreporter_enabled = True \ntb=casatools.table() \nospathsep = os.path.sep \nlibpath = sysconfig.get_paths()['purelib'] \nobsdict = {'MJD': 57447.0, 'Name': 'EOVSA', 'Type': 'WGS84', 'Long': -118.287, \n            'Lat': 37.2332, 'Height': 1207.13, 'X': 0.0, 'Y': 0.0, 'Z': 0.0,  \n            'Source': 'Dale Gary'} \nobstable = os.path.join(casadata.datapath,'geodetic','Observatories') \ntb.open(obstable, nomodify=True) \nif 'EOVSA' not in tb.getcol('Name'): \n    print('Adding EOVSA to the Observatories') \n    tb.close() \n    tb.open(obstable, nomodify=False) \n    nrows = tb.nrows() \n    tb.addrows(1) \n    for k in obsdict.keys(): \n        tb.putcell(k, nrows, obsdict[k])     \ntb.close() \n" > $HOME/.casa/config.py
```

***Last but not least, add the `site-packages` path where the suncasa and its dependent packages installedd to the Monolithic CASA***
```bash
echo "sitepackagepath = '$HOME/suncasaenv/lib/python3.8/site-packages'" >> $HOME/.casa/config.py
echo "if sitepackagepath not in sys.path:" >> $HOME/.casa/config.py
echo "    sys.path.append(sitepackagepath)"     >> $HOME/.casa/config.py
```
>**_note:_** The `sitepackagepath` may varies for python versions on different platforms.

[comment]: <> (Now pick whichever subset of the available CASA packages you are interested in. Package dependencies are handled)

[comment]: <> (automatically by pip, with the exception of **casadata** which must be explicitly installed and updated by the user &#40;)

[comment]: <> (see [External Data]&#40;external-data.html#Updating-the-Data-Tables&#41;&#41;. The following packages are available:)

[comment]: <> (```bash)

[comment]: <> (pip install casatools==6.5.2.26)

[comment]: <> (pip install casatasks==6.5.2.26)

[comment]: <> (pip install casaplotms==1.8.7)

[comment]: <> (pip install casaviewer==1.6.6)

[comment]: <> (pip install casampi==0.5.01)

[comment]: <> (pip install casashell==6.5.2.26)

[comment]: <> (pip install casadata==2022.9.5)

[comment]: <> (pip install casaplotserver==1.4.6)

[comment]: <> (pip install almatasks==1.5.2)

[comment]: <> (pip install casatestutils==6.5.2.26)

[comment]: <> (```)

**Note for Mac M1 users**: For macOS 12 on an ARM-based M1 chip, users will need to install the wheels of CASA version
11 for x86 architecture. For that, we recommend to use the following command to pip install the CASA wheels:

```bash
arch -x86\_64 python3.8 -m pip install ...
```

[comment]: <> (The use of python3 venv is a simple built-in method of containerizing the pip install such that multiple versions of)

[comment]: <> (CASA 6.x can be kept on a single machine in different environments. In addition, CASA is built and tested using)

[comment]: <> (standard &#40;python 3.6&#41; libraries which can be replicated with a fresh venv, keeping the libraries needed for CASA)

[comment]: <> (isolated from other libraries which may already be installed on your machine.)

[comment]: <> (With the pip installation, CASA may be used in a standard Pythonic manner. Examples can be found)

[comment]: <> (in [this Jupyter Notebook]&#40;https://colab.research.google.com/github/casangi/examples/blob/master/casa6/CASA6_demo.ipynb&#41;)

[comment]: <> (.)
