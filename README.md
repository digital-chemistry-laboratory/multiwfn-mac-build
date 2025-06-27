# Introduction

[Multiwfn](http://sobereva.com/multiwfn/) is one of the most versatile software packages for electronic wave function analysis. Sadly, support for MacOS was dropped after version 3.7. This is an independent CMake-based build recipe which I have found to work reasonably well on an M1-based Mac running MacOS Monterey 12.4. Others have reproduced it on Intel Macs. No guarantees given for anything.

# Homebrew tap

The easiest way to install is with the [Homebrew tap](https://github.com/digital-chemistry-laboratory/homebrew-multiwfn).

# Requirements

The build has been tested with Homebrew-installed:
- [GCC](https://formulae.brew.sh/formula/gcc) 
- [OpenBLAS](https://formulae.brew.sh/formula/openblas)
- [CMake](https://formulae.brew.sh/formula/cmake)
- [arb](https://formulae.brew.sh/formula/arb)
- [flint](https://formulae.brew.sh/formula/flint)

Apple's Accelerate for BLAS/LAPACK also seems to work. No guarantees are given for other combinations of compilers and linear algebra backends. The build recipe might break with future versions of Multiwfn.

# Building from source distribution

The [source_dist](https://github.com/digital-chemistry-laboratory/multiwfn-mac-build/tree/source_dist) branch is updated daily based on the published Multiwfn [source files](http://sobereva.com/multiwfn/download.html). We clone this branch and build with cmake:

```zsh
$ git clone --branch source_dist https://github.com/digital-chemistry-laboratory/multiwfn-mac-build.git
$ cd multiwfn-mac-build
$ cmake -B build
$ cmake --build build
$ cp build/multiwfn .
```

Alternatively, use `cmake -B build -DCMAKE_BUILD_TYPE=Release` to fully optimize the build, but it will take additonal time.

We then need to add the Multiwfn path to `~/.zshrc`. Adapt the specific path to where you have installed Multiwfn on your system. See the Multiwfn manual for more specific instructions.

```zsh
export PATH=$PATH:$HOME/bin/multiwfn-mac-build
export Multiwfnpath=$HOME/bin/multiwfn-mac-build
```

## Building with OpenMP

To build with OpenMP, set the flag `cmake -B build -DWITH_OpenMP=ON`. Then you need to specify the number of processors in `settings.ini` (see below).

## Building with OpenBLAS

To link against OpenBLAS rather than Apple's Accelerate, on Apple Silicon the following environment variables have to be set before calling CMake.

```
export LDFLAGS="-L/opt/homebrew/opt/openblas/lib"
export PKG_CONFIG_PATH="/opt/homebrew/opt/openblas/lib/pkgconfig"
```

## Install

Use `cmake --install build` to install the `multiwfn` executable and `settings.ini`. Change install prefix with the flag `--prefix <dir>`

## Building with GUI

The GUI version can be built on Intel machines and seem to work with bugs on Apple Silicon machnines. The required libraries can be installed with homebrew:

```shell
brew install brew install libx11 libxt mesa openmotif
```

### Install DISLIN

The DISLIN package needs to be installed on the system. Follow these steps:
1. Download the right distribution from the DISLIN download [page](https://www.dislin.de/darwin.html)
2. Unpack the .tar.gz file and go into the folder
3. Set the preferred install directory for DISLIN
4. Run the installation script
5. Patch the library install name

An example of how to do this is:
```shell
wget https://www.dislin.de/downloads/macOS/dislin-11.5.macOS.intel.64.tar.gz
# for Apple Silicon: wget https://www.dislin.de/downloads/macOS/dislin-11.5.macOS.arm.64.tar.gz
tar -xvf dislin-11.5.macOS.intel.64.tar.gz
export DISLIN=$HOME/bin/dislin
cd dislin-11.5
./INSTALL
cd ..
install_name_tool -id $DISLIN/libdislin_d.dylib $DISLIN/libdislin_d.dylib
```

### Building

To build with the GUI, specify the following options to CMake,

```shell
$ cmake -B build -DWITH_GUI=ON -DDISLIN_DIR=$DISLIN
```

where `$DISLIN` points to the same directory as during the installation. These options can be combined with OpenMP and the release build type for optimal performance (see above).

### Running

To avoid memory issues, make sure to follow the instructions in the Multiwfn manual and set the `OMP_STACKSIZE` before executing the program.

```shell
export OMP_STACKSIZE=64000000
multiwfn
```