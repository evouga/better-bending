# Better Bending Source Code

This repository contains source code for our paper **Better Bending: Analysis, Construction and Verification of Discrete Bending Models for Kirchhoff-Love Shells**.

- [`energy-benchmark/`](energy-benchmark/) — the energy convergence and creased vee benchmarks (Section 8 and 9)
- [`lateral-buckling/`](lateral-buckling/) — our reproduction of [Romero et al.]'s lateral buckling experiment (Section 10.3)

Both projects depend on [LibShell](https://github.com/evouga/libshell) and [Polyscope](https://github.com/nmwsharp/polyscope). These dependencies are included as Git submodules at the repository root.

## Getting started

Clone the repository and fetch submodules in one step:

```bash
git clone --recurse-submodules <repository-url>
```

If you have already cloned without submodules, initialize them with:

```bash
git submodule update --init --recursive
```

## Building

From the repository root, configure and build all projects (libshell, then both experiments):

```bash
cmake -S . -B build
cmake --build build
```

This produces the `energy-benchmark` and `lateral-buckling` executables under `build/`.

Each subfolder can also be built on its own with the standard CMake workflow after libshell has been built and installed to `libshell/lib/`.
