# ROOT-Sim Random Number Generators Library

[![Continuous Integration](https://github.com/ROOT-Sim/random-number-generators/actions/workflows/build_and_test.yml/badge.svg)](https://github.com/ROOT-Sim/random-number-generators/actions/workflows/build_and_test.yml)
[![GitHub issues](https://img.shields.io/github/issues/ROOT-Sim/random-number-generators)](https://github.com/ROOT-Sim/random-number-generators/issues)
[![REUSE status](https://api.reuse.software/badge/github.com/ROOT-Sim/random-number-generators)](https://api.reuse.software/info/github.com/ROOT-Sim/random-number-generators)

*Brought to you by the [High Performance Computing & Simulation (HPCS)](https://hpcs.ing.uniroma2.it/) Research Group*

----------------------------------------------------------------------------------------

## Overview

The ROOT-Sim Random Number Generators (RNG) Library provides a collection of Piece-Wise Deterministic Random Number Generators designed for use in simulation models. It serves as a core component for the [ROOT-Sim](https://github.com/ROOT-Sim/cROOT-Sim) framework, delivering robust, high-performance, and reproducible random streams for Logical Processes.

## Features

This library supports multiple common statistical distributions and relies on a fast underlying Pseudo-Random Number Generator (PRNG), such as Xoroshiro, to ensure performance and statistical quality.

Available distributions and utilities include:

*   **Uniform Distributions:** Standard continuous random numbers and discrete 64-bit unsigned integers.
*   **Normal Distributions:** Standard normally distributed random values.
*   **Exponential Distributions:** Exponentially distributed random numbers for inter-arrival times.
*   **Gamma Distributions:** Support for the Gamma distribution shape parameters.
*   **Poisson Distributions:** Discrete probability distributions for event occurrences.
*   **Zipf Distributions:** Skewed distributions often used to model popularity or access patterns.
*   **Range Selections:** Uniform and non-uniform selections across defined minimum and maximum bounds.

## Integration

To utilize the random number generation facilities in your model, include the main header:

```c
#include <ROOT-Sim/random.h>
```

This project relies on CMake and requires a C11-compliant compiler. To build and link the library, follow the standard CMake build process:

```bash
mkdir build
cd build
cmake ..
make
```

## License

This project is licensed under the GPL-3.0-only license, with some components utilizing the CC0-1.0 license. See the `LICENSES` directory or the REUSE configuration for more precise details regarding file-level licensing.
