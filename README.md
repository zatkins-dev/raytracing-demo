# raytracing-demo

This program implements raytracing between two concentric cylinders.
One ray is cast from each inward face of the outer cylinder to the center of each other inward face of the outer cylinder.
We compute collisions with the inner, blocking cylinder.

## Requirements
- gnuplot
- [NVIDIA HPC SDK 22.11](https://developer.nvidia.com/nvidia-hpc-sdk-downloads)

## Running
The [run](scripts/run) script will run a full suite of tests.
If you would prefer to run a single test, use the command line options.
