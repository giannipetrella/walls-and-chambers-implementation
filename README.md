# Walls and chambers for quiver moduli

This repository is part of the paper "Finding the walls for quiver moduli"
by Hans Franzen, Gianni Petrella and Rachel Webb,
[arXiv:2506.20568](https://arxiv.org/abs/2506.20568).

This code can be used, among other things, to verify every example in the paper.

## Instructions

You will need a Sagemath and/or a Julia installation. Get each from their official websites.

### How to set up the Julia environment

Open a shell in a directory where you want these files, and launch the Julia REPL
by running

```bash
julia --project=.
```

Install QuiverTools by running the following commands, in the given order.

```julia
julia> using Pkg

julia> Pkg.add(url="https://github.com/QuiverTools/QuiverTools.jl", rev="f1647cc");
```

Now install Oscar and Memoize by running
```julia
julia> Pkg.add("Oscar"); Pkg.add("Memoize");
```

To save raw computation results, install Serialization.jl via

```julia
julia> Pkg.add("Serialization")
```

To run the script, close the Julia REPL with `Ctrl+D` and run the following command in the shell:

```bash
julia --project=. examples.jl
```

This is configured by default to reproduce the computations for the Segre cubic.
You can use this code to experiment on any other example by importing it in an interactive
session; do so by launching a REPL (with `julia --project=.`) and running

```julia
julia> include("walls-and-chambers.jl")
```

### How to set up the Sagemath environment

In the directory where you want these files,
install QuiverTools by running the following command in a shell:

```bash
sage --pip install git+https://github.com/QuiverTools/QuiverTools.git
```

You can then experiment with the code after loading it with

```sagemath
sage: load("git_equivalence.py")
```
