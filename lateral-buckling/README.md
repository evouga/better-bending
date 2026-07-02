# Lateral Buckling Experiment

This code will run the lateral buckling experiment for the BAC, QB(PL), and QB(CR) methods.

## Building and Running

Uses the standard CMake toolchain. Requires:
- LibShell (https://github.com/evouga/libshell)
- Polyscope (https://github.com/nmwsharp/polyscope)

The simplest setup is to place these in sibling directories under the same parent as the folder containing this README.

## Details

The code outputs the results in log.txt. Each row tests one w/L, Gamma\* pair and is of the format:

```
w/L: BAC lateral displacement, QB(PL) lateral displacement, QB(CR) lateral displacement.
```

The full simulation will take several hours to finish.