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

Currently, it is not possible to build Multiwfn with a GUI on machines with an Apple Silicon processor, due to lack of a [DISLIN](https://www.dislin.de) distribution. If someone has an Intel Mac and would like to work out a recipe, please see the open issues on the question.
