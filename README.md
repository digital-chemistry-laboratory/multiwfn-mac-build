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

```console
$ wget http://sobereva.com/multiwfn/misc/Multiwfn_3.8_dev_src_Linux.zip
$ unzip Multiwfn_3.8_dev_src_Linux.zip && mv Multiwfn_3.8_dev_src_Linux/* . && rmdir Multiwfn_3.8_dev_src_Linux
$ sed -i '' 's/if ((.not.present(irewind)).or.(present(irewind).and.irewind==1)) rewind(fileid)/if (.not. present(irewind)) then\n  rewind(fileid)\nelse\n  if (irewind==1) rewind(fileid)\nend if/g' util.f90
$ cmake -B build
$ cmake --build build
$ cp build/multiwfn .
```

`sed` is used to patch one of the source files so that it can be compiled with GFortran.

We then need to add the Multiwfn path to the `~/.zshrc`. Adapt the specific path to where you have installed Multiwfn on your system. See the Multiwfn manual for more specific instructions.

```zsh
export PATH=$PATH:$HOME/bin/Multiwfn
export Multiwfnpath=$HOME/bin/Multiwfn
```

# Settings

The source distribution of Multiwfn doesn't come with the file `settings.ini`. You can [download](http://sobereva.com/multiwfn/) the binary distribution for Linux and copy just the file `settings.ini` into your Multiwfn folder.