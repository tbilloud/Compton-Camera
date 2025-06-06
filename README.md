# Timepix3 Compton camera simulation and data processing

![Screenshot of a reconstruction](doc/img.png)

A single python script to:
- Simulate particle transport via Geant4 and detector response via Allpix2
- Reconstruct cones from:
  - Geant4/Gate 'hits'
  - Gate 'singles'
  - Allpix² 'pixel hits' (WIP)
  - measured data (WIP)
- Validate simulation/measurement of gamma point sources
- Reconstruct 3D image with basic back-projection
- Visualize 3D images with napari (Linux only)

Requires:
- Linux, MacOS, or Windows + WSL
- 20 GB of disk space
- Python >= 3.10 and pip
- Allpix²
- Optional: CUDA, napari

Tested with:
- OS: 
  - Ubuntu 22.04 / 24.04
  - MacOS Sequoia (15.4.1)
  - Windows + WSL Ubuntu 22.04 (some QT issues, solved with ChatGPT)
- Python:
  - 3.10, 3.11, 3.12 (issue with OpenSSL with 3.9 and opengate-core with 3.13)
- OpenGate:
  - 10.0.1
  - 10.0.2 (problem with opengate/Allpix interface, see TODOs)
- Allpix²:
  - 3.1.0
- GPU: 
  - cupy-cuda115/128 + GeForce RTX 2080 Ti

Python3 installs by default:
opengate==10.0.0 if 3.9, 3.10
opengate==10.0.2 if 3.11

## [Installation](#install)

If not already installed, install a suitable python version and pip.
MacOS: use pyenv for example

### 1) Download or clone: 
```
git clone https://github.com/tbilloud/ComptonCamera
```

### 2) Create a virtual environment
```
cd ComptonCamera
python3 -m venv venv
source venv/bin/activate
```

### 3) Install required python packages
Ubuntu: `pip install -r requirements-OS.txt`  
MacOS: `pip install -r requirements-macos.txt`

### 4) Install Allpix²
 
#### Prerequisites Ubuntu 
```
sudo apt-get install libboost-all-dev # installs BOOST
sudo apt-get install libeigen3-dev # installs Eigen3
wget https://root.cern/download/root_v6.32.10.Linux-ubuntu22.04-x86_64-gcc11.4.tar.gz
tar -xzvf root_v6.32.10.Linux-ubuntu22.04-x86_64-gcc11.4.tar.gz
source root/bin/thisroot.sh # installs ROOT6
```
Note: adapt lines 3-4 to your system (https://root.cern/install/all_releases/)

#### Prerequisites MacOS 
```
xcode-select --install # install Command Line Developer Tools, if not already there
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)" # installs brew
brew install cmake # installs CMake
brew install boost # installs BOOST
brew install eigen # installs Eigen3
brew install root # installs ROOT6
```
Notes:
- Only tested with brew. Other ways might work too (macports, or from source).
- As of April 2025, brew installs ROOT 6.34.08 built with C++17 (needed by Allpix²)


#### Then install Allpix² without Geant4:  
https://allpix-squared.docs.cern.ch/docs/02_installation/  
```
cd allpix
git clone https://gitlab.cern.ch/allpix-squared/allpix-squared
rm -rf .git # remove the git folder to avoid conflicts
cd allpix-squared
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=../install-noG4 -DBUILD_GeometryBuilderGeant4=OFF -DBUILD_DepositionCosmics=OFF -DBUILD_DepositionGeant4=OFF -DBUILD_DepositionGenerator=OFF -DBUILD_GDMLOutputWriter=OFF -DBUILD_VisualizationGeant4=OFF ..
make -j4
make install
```

### 5) Set environment (each time you open a new terminal)
```
export GLIBC_TUNABLES=glibc.rtld.optional_static_tls=2000000
```

If using PyCharm:
- Expand the box with the name of the main script on the top (near ▶)
- Edit Configurations > Path to ".env files" > /path/Compton-Camera/.env
- Replace the absolute project paths with your own

If environment is not set properly, you might get errors like:  
`ImportError: .../libG4geometry-cf4c216c.so: cannot allocate memory in static TLS block`
`ERROR: ld.so: object ... from LD_PRELOAD cannot be preloaded (...): ignored.'`

`QObject::moveToThread: Current thread (0x5dfc4786c040) is not the object's thread...`  
=> for this error see section [QT/opengate conflict](#qtopengate-conflict) below.

### 6) Optional: Install GPU tools
To use the GPU-based functions (point source validation, reconstruction):
a) Install CUDA (https://developer.nvidia.com/cuda-downloads)
b) Install the Cupy package suited to your CUDA version, e.g.  
`pip install cupy-cuda115`

## [Installation without simulation packages](#install-offline)

Use cases:
- For analyzing measured data, avoiding opengate and Allpix² saves time and disk space.
- For visualizing with napari on MacOS, a 2nd virtual environment is needed (see QT/opengate conflict below).

In opengate, running code without Monte Carlo (geant4) simulation is called 'offline'.

If not already done, clone the repository:
```
git clone https://github.com/tbilloud/ComptonCamera
cd ComptonCamera
```

Create a virtual environment without opengate and Allpix²:
```
python3 -m venv venv-offline
source venv-offline/bin/activate
pip install -r requirements-offline.txt
```

Try the offline example:
```
python3 main_offline.py
```


## [Getting started](#getting-started)
Run the tests:  
```
python3 main_basic.py
python3 main_allpix.py`
```
The 1st time you run a simulation, Gate10 will install Geant4 datasets, which can take a while. This is done only once.

In main.py, the 1st part (code until block 'ANALYSIS AND RECONSTRUCTION') is the Gate 10 simulation. See user manual:  
https://opengate-python.readthedocs.io/en/master/user_guide/user_guide_intro.html

Then, several functions are available. Step-by-step:
1) Print basic info about Gate output files 
- Gate hits with analyze_hits(). Gate hits are like Geant4 hits and are different from pixel hits.
- singles, which are group of hits per pixel, with analyze_singles()
2) Simulate pixel hits:
- from Gate singles with gSingles2pixelHits()
- from Allpix² output with gHits2allpix2pixelHits()
3) Reconstruct cones:
- from Gate4 hits with gHits2cones_byEvtID()
- from pixel hits (WIP)
4) Check cones from a point sources:
- validate_psource() plots cone projections. It's slow, ~1 sec per cone.
- with GPU acceleration with validate_psource_gpu()
5) Reconstruct 3D image with:
- simple backprojection with backprojection(). It's slow, ~1 sec per cone.
- PU-accelerated backpropagation with backprojection_gpu()

For Linux users, potting functions using napari are available:
- scroll between cones with plot_stack_napari()
- show reconstructed source and detector geometry in 3D with plot_reconstruction_napari()
Does not work on MacOS yet, see QT/opengate conflict below.


## [Allpix²](#allpix2)

Allpix² is a C++ software for precise simulation of semiconductor pixel detectors.
It simulates the transport of charge carriers in semiconductor sensors and their signal induction.
It is used primarily for detector R&D in particle physics.  
https://cern.ch/allpix-squared


Allpix² can read the hits root file from Gate.
Combined with Gate10, the entire simulation can be done with a single python file, using the function
gHits2allpix2pixelHits() after the sim.run() in the main.py script. It does the following:
1) run Gate10 and creates the hits root file
2) generate the three .conf files needed by Allpix²
3) run Allpix² and creates the output files data.txt and modules.root in the sub-folder 'allpix'
4) read data.txt and return a pandas dataframe with the pixel hits

An Allpix² simulation needs 3 configuration (.conf) files:
- detector geometry
- detector model
- simulation parameters

The main configuration file (simulation parameters) is a 'simulation chain' made of several components:
- global parameters
- electric field
- charge deposition
- charge propagation
- charge transfer
- digitization
https://allpix-squared.docs.cern.ch/docs/03_getting_started/06_simulation_chain/

For each component, several modules are available.
- Charge propagation:
  - ProjectionPropagation: fast but only silicon sensors, linear electric field, one carrier type at a time
  - GenericPropagation
  - TransientPropagation
- Charge transfer:
  - SimpleTransfer: no ToA (as of v3.1.0, 2025-01-08)
  - CapacitiveTransfer
  - PulseTransfer
  - InducedTransfer
- Digitization:
  - DefaultDigitizer
  - CSADigitizer


## [QT/opengate conflict](#qt)
Napari and OpenGate use the same QT backends (PyQt5), which causes conflicts.
Using Qt-based code (e.g. napari) and gate in the same script leads to warning/crashes.  
### Ubuntu
```
WARNING: QObject::moveToThread: Current thread (0x57ad941535d0) is not the object's thread (0x57ad94c1ef50).
Cannot move to target thread (0x57ad941535d0)
```
Solution:
```
mv /path-to-virtual-environment/lib/python3.XX/site-packages/opengate_core/plugins /path-to-virtual-environment/lib/python3.XX/site-packages/opengate_core/plugins.bak
```
-> replace XX with your python version, e.g. 10

### Macos
```
objc[16117]: Class QT_... is implemented in both .../opengate_core/.dylibs/QtCore and .../QtCore (0x16c7d1278) ... One of the duplicates must be removed or renamed.
objc[16117]: Class KeyV... is implemented in both .../opengate_core/.dylibs/QtCore and .../QtCore (0x16c7d12a0) ... One of the duplicates must be removed or renamed.
objc[16117]: Class RunL... is implemented in both .../opengate_core/.dylibs/QtCore and .../QtCore (0x16c7d12f0) ... One of the duplicates must be removed or renamed.
```
Solution: TODO !
