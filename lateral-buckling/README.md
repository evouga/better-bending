# Lateral Buckling Experiment

This code will run our implementation of [Romero et al.]'s lateral buckling experiment.

## Building and Running

Uses the standard CMake toolchain. Requires:
- LibShell (https://github.com/evouga/libshell)
- Polyscope (https://github.com/nmwsharp/polyscope)

These dependencies are provided as Git submodules at the repository root. From the repository root, clone with submodules:

```bash
git clone --recurse-submodules <repository-url>
```

Or, if you already cloned the repository:

```bash
git submodule update --init --recursive
```

The build expects `libshell/` and `polyscope/` as siblings of this directory (see the [repository README](../README.md)).

## Details

The code outputs the results in log.txt. Each row tests one w/L, Gamma\* pair and is of the format:

```
w/L: BAC lateral displacement, QB(PL) lateral displacement, QB(CR) lateral displacement.
```

The full simulation will take several hours to finish.