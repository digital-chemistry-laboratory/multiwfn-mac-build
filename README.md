# Introduction

[Multiwfn](http://sobereva.com/multiwfn/) is one of the most versatile software packages for electronic wave function analysis. Sadly, support for MacOS was dropped after version 3.7. This is an independent CMake-based build recipe which I have found to work reasonably well on an M1-based Mac running MacOS Monterey 12.4. No guarantees given for anything.

# Requirements

The build has been tested with Homebrew-installed:
- GFortran 
- OpenBLAS
- CMake

No guarantees are given for other combinations of compilers and linear algebra backends. The build recipe might break with future versions of Multiwfn.

# Building from source distribution

Here we get the source distribution of Multiwfn 3.8dev with `wget`. The latest versions can be found on the Multiwfn [download page](http://sobereva.com/multiwfn/). Clone the repository and run the following commands in the repository directory.

```zsh
$ wget http://sobereva.com/multiwfn/misc/Multiwfn_3.8_dev_src_Linux.zip
$ unzip Multiwfn_3.8_dev_src_Linux.zip && mv Multiwfn_3.8_dev_src_Linux/* . && rmdir Multiwfn_3.8_dev_src_Linux
$ patch < util.patch
$ cmake -B build
$ cmake --build build
$ cp build/multiwfn .
```

Alternatively, use `cmake -B build -DCMAKE_BUILD_TYPE=Release` to fully optimize the build, but it will take additonal time.

We then need to add the Multiwfn path to the `~/.zshrc`. Adapt the specific path to where you have installed Multiwfn on your system. See the Multiwfn manual for more specific instructions.

```zsh
export PATH=$PATH:$HOME/bin/Multiwfn
export Multiwfnpath=$HOME/bin/Multiwfn
```

# Settings

The source distribution of Multiwfn doesn't come with the file `settings.ini`. You can [download](http://sobereva.com/multiwfn/) the binary distribution for Linux and copy just the file `settings.ini` into your Multiwfn folder.
