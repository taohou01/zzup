# Updating zigzag barcodes and representatives for computing vineyards for dynamic point clouds

This project implements the algorithms for computing vines and vineyards for dynamic point clouds using representative-based zigzag update algorihtms. The approaches are described in the following paper `[1]`:

[Updating Barcodes and Representatives for Zigzag Persistence](https://arxiv.org/pdf/2112.02352.pdf)

by Tamal K. Dey and Tao Hou. The programming language used is C++.

## Group Information

This project is developed by [Tao Hou](https://taohou01.github.io) under the [CGTDA](https://www.cs.purdue.edu/homes/tamaldey/CGTDAwebsite/) research group at Purdue University lead by [Dr. Tamal K. Dey](https://www.cs.purdue.edu/homes/tamaldey/).

## About the Implementation

Given a dynamic point cloud as input, the software first builds a sequence of simplex-wise operations on the zigzag filtration, which consists of the following as in `[1]`:

- forward switches
- backward switches
- outward switches
- outward expansions
- inward contractions

It then performs updates on the zigzag barcodes for the operations and builds *vines and vineyards* `[2]` for the input data. More explanations on the implementation are provided in the Section [Implementation Details](https://github.com/taohou01/zzup/edit/main/README.md#implementation-details)

## Building

The building utilized [`cmake`](https://cmake.org/) software, and all building problems should be solvable by inspecting and changing the `CMakeLists.txt` file. The building has two dependencies: one is `boost`, which is quite standard (see `CMakeLists.txt`); the other is [`phat`](https://github.com/blazs/phat) used in the [`FastZigzag`](https://github.com/taohou01/fzz) implementation for computing zigzag barcodes from scratch for timing comparisons (see [Implementation Details](https://github.com/taohou01/zzup/edit/main/README.md#implementation-details)). For [`phat`](https://github.com/blazs/phat), users should first download the codes themselves and then include the header files of phat into the project by adding

```
include_directories( "path/to/phat/include" ) 
```

into CMakeLists.txt.

Commands for building are quite standard:

```
cd [dir-to-zzup]
mkdir build
cd build
cmake ..
make
```

The software is developed and tested under MacOS. Compiling and running under other platforms (e.g., Linux) should not have any problems as only standard features of C++ are utilized.

## Usage

The software runs with following command:

```
./dpc_vine [OPTIONS] input_dpc_file
```

A sample input file specifying a dynamic point cloud (DPC) is provided as `sample_in.txt` with the source codes. A `input_dpc_file` starts with a line specifying the number (`N`) of points and each `N+1` lines that follow specify the positions of the points at a time `i`, with `i` always starting from `0`. The first of the `N+1` lines is `# t i` and each remaining line starts with an id for the point followed by its position (currently only in *2D*).

## Implementation Details

## References

1. Tamal K. Dey and Tao Hou. Updating Barcodes and Representatives for Zigzag Persistence.
2. David Cohen-Steiner, Herbert Edelsbrunner, and Dmitriy Morozov. Vines and vineyards by
updating persistence in linear time. *The Twenty-Second Annual Symposium
on Computational Geometry*.
