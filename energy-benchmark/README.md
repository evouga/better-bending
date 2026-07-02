# Bending Model Energy Verification Benchmark

This code will run the full test suite of energy scaling and creased vee experiments in Sections 8 and 9 of our paper.

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

The code outputs the results in log.txt. Experiments can be turned on/off in main.cpp.

The full simulation will take many hours to finish. A machine with at least 32 GB of RAM is strongly recommended.