# Bending Model Verification Benchmark

This code will run the full test suite of experiments probing the energy convergence behavior of discrete shell models.

## Building and Running

Uses the standard CMake toolchain. Requires:
- LibShell (https://github.com/evouga/libshell)
- Polyscope (https://github.com/nmwsharp/polyscope)

The simplest setup is to place these in sibling directories under the same parent as the folder containing this README.

## Details

The code outputs the results in log.txt. Experiments can be turned on/off in main.cpp.

The full simulation will take many hours to finish. A machine with at least 32 GB of RAM is strongly recommended.