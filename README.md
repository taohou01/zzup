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

A sample input file specifying a dynamic point cloud (DPC) is provided as `sample_in.txt` with the source codes. An `input_dpc_file` starts with a line specifying the number `N` of points and each `N+1` lines that follow specify the positions of the points at a time `i`, with `i` always starting from `0`. The first of the `N+1` lines is `# t i` and each remaining line starts with an id for the point followed by its position (currently only in *2D*).

The following options control the behavior of the program: 
- `-d` specifies the maximum dimension for the Rips simplices considered (default to 2); 
- `-s` specifies the starting time considered for the input DPC (default to 0); 
- `-e` specifies the ending time considered for the input DPC (default to maximum time in input file);
- `--fzz` makes the software also invoke [`FastZigzag`](https://github.com/taohou01/fzz) for computing barcodes from scratch for the purpose of timing comparison;
- `-E` prints some status messages concerning the edge operations performed during the run;
- `--sop` prints some status messages concerning the simplex-wise operations performed during the run (automatically turns on `-E`).

The name of the output file encoding the vineyard produces is of the following form

```
[INPUT_FILENAME]_d_[X]_t_[Y]_[Z]_vines.txt
```
where `[X]` is the maximum dimension and `[Y]`/`[Z]` are the starting/ending time. For example, the following command

```
./dpc_vine -d 3 -s 1 -e 5 ../sample_in.txt 
```

generates an output file `sample_in_d_3_t_1_5_vines.txt`.

As in `[1]`, the vineyard specified in an output file consists of 2D points (birth-death time for the persistence bars) sweeping through a third dimension (distance threshold for Rips) forming vines (piece-wise linear lines connected by points). The first line of the output file specifies the maximum distance of points in the DPC. The reason for specifying such max distance is to provide a reference for manullay choosing a finite capping value for the third dimension for those vines starting from infinite distance.

The remaining lines of an ouput file specifies the vines in the vineyard, each of which could start with `s`, `c`, or `e`. 

- A line starting with `s` indicates the start of a vine; the following five numbers specify the (unique) id, birth time (x-value), death time (y-value), distance threshold (z-value), and dimension for the vine and its starting point. 
- A line starting with `c` denotes the next point for the vine (identified by its id), connecting to the previous point (specified by a line starting with `s` or `c`); the following four numbers are the same as the first four numbers (excluding dimension) in a line starting with `s`. 
- The format and meaning of a line starting with `e` are the same as one starting with `c` with the only difference that we now know for sure that the vine ends. 

Notice that the id's for the vines are only guaranteed to be unique; the id's may not start with `0`, nor being consecutive or monotonically increasing. 

A python script `disp_vines.py` is provided for displaying the vines for an output file.

## Implementation Details

The implementation can be roughly broken into two parts:

- The *front end* (`dpc.h/cpp`) generates the *distance-time curves* (see `[1]`) from the DPC and detects all the *critical events* (local min/max and intersections) of the curves. It then generates *edge-level* operations from the critical events (inserting/removing consecutive addition/deletion of an edge in the middle or switching two addition/deletion of edges). From the edge-level operations, it generates the *simplex-wise* opeartions on the zigzag filration (aka. the switches, expansions, and contractions) and performs the barcode updates using the back end. During the barcode updating, points for vines are being output.

- The *back end* (`dynamic_zigzag.h/cpp`) implements the representative-based update algorithms described in `[1]` for the five operations which are used by the front end.

We also notice the following:

- Since we sweep the distance threhold in decreasing order, the starting filtration is an *up-down* zigzag filtration which corresponds to the infinite distance. The representatives for the up-down filtration are computed from an algorithm similar to the one described in Appendix A of `[3]` ([link](https://arxiv.org/pdf/2105.00518.pdf)).

- For efficient implementation, when manipulating chains during the update of representatives, contents of a chain (which is an increasing array of simplex id's) never change once the chain is initially set. Hence, with [`std::share_ptr`](https://en.cppreference.com/w/cpp/memory/shared_ptr), copying a chain is nothing but copying the chain's smart pointer (which is constant time). Whenever summing two chains, the resulting chain always points to a newly allocated memory location. Doing this can avoid a lot of unnecessary chain copies.

- [`FastZigzag`](https://github.com/taohou01/fzz) is used for computing zigzag barcodes from scratch (for the purpose of timing comparison). Codes of [`FastZigzag`](https://github.com/taohou01/fzz) are wrapped into a class (`fzz.h/cpp`). [phat](https://github.com/blazs/phat) is used for computing the standard (non-zigzag) persistence in [`FastZigzag`](https://github.com/taohou01/fzz).


## References

1. Tamal K. Dey and Tao Hou. Updating Barcodes and Representatives for Zigzag Persistence.
2. David Cohen-Steiner, Herbert Edelsbrunner, and Dmitriy Morozov. Vines and vineyards by
updating persistence in linear time. *The Twenty-Second Annual Symposium
on Computational Geometry*.
3. Tamal K. Dey and Tao Hou. Computing Optimal Persistent Cycles for Levelset Zigzag on Manifold-like Complexes.
