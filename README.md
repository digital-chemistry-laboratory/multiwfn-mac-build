# Introduction

[Multiwfn](http://sobereva.com/multiwfn/) is one of the most versatile software packages for electronic wave function analysis. Sadly, support for MacOS was dropped after version 3.7. This is an independent CMake-based build recipe which I have found to work reasonably well on an M1-based Mac running MacOS Monterey 12.4. No guarantees given for anything.

# Homebrew tap

The easiest way to install is with the [Homebrew tap](https://github.com/kjelljorner/homebrew-multiwfn).

# Requirements

The build has been tested with Homebrew-installed:
- GFortran 
- OpenBLAS
- CMake

No guarantees are given for other combinations of compilers and linear algebra backends. The build recipe might break with future versions of Multiwfn.

# Building from source distribution

The [source_dist](https://github.com/kjelljorner/multiwfn-mac-build/tree/source_dist) branch is updated daily based on the published Multiwfn [source files](http://sobereva.com/multiwfn/download.html). We clone this branch and build with cmake:

```zsh
$ git clone --branch source_dist https://github.com/kjelljorner/multiwfn-mac-build.git
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

## Install

Use `cmake --install build` to install the `multiwfn` executable and `settings.ini`. Change install prefix with the flag `--prefix <dir>`

## Building with GUI

To build with GUI version, the following additional packages and libraries can be installed via Homebrew and configured correctly, respectively. 
- openmotif
- libxt
- XQuartz
- x11lib
Commandline tools can be installed via `xcode-select --install`.
[DISLIN](https://www.dislin.de) can be downloaded in official website and configured with _Official User Manual_. 

The packages LIB files should linked to your building envirment. Add the library PATH to your bash profile such as `.zprofile` or `.zshrc` or `$export` it.

```zsh
$ vi ~/.zprofile #or .zshrc
```

#Add these to `.zprofile` or `.zshrc` file in home folder. Do not just copy the following fields. All libraries should be confiugred correctly.

```
PATH=$PATH:/other/executable/path/folder

DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/usr/lib:/usr/local/lib:/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib:/your/path/to/libxt/1.2.1/lib:/your/path/to/dislin/lib:/your/path/to/openmotif/lib

#In order to let cmake find libdislin, LIBRARY_PATH should be configured, LD_LIBRARY_PATH can be configured additonaly.

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$DYLD_LIBRARY_PATH
LIBRARY_PATH=$DYLD_LIBRARY_PATH

export DYLD_LIBRARY_PATH
export LD_LIBRARY_PATH
export LIBRARY_PATH
export PATH

```
Restart the terminal or use `$source ~/.zprofile or .zshrc` , then build the GUI version multiwfn.

```zsh
$ git clone --branch source_dist https://github.com/kjelljorner/multiwfn-mac-build.git
$ cd multiwfn-mac-build
$ cmake -B build -DWITH_OpenMP=ON -DWITH_GUI=ON
$ cd build
$ make -j6    # or just use make
$ cp multiwfn ../.
```
Then we can try execute the multiwfn to find if we build it correctly with GUI.

After that, we can install the multiwfn just as previous steps and configured it into PATH.
