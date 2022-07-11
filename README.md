# Updating zigzag barcodes and representatives for computing vineyards for dynamic point clouds

This project implements the algorithms for computing vines and vineyards for dynamic point clouds using representative-based zigzag update algorihtms. The approaches are described in the following paper `[1]`:

[Updating Barcodes and Representatives for Zigzag Persistence](https://arxiv.org/pdf/2112.02352.pdf)

by Tamal K. Dey and Tao Hou.

## Group Information

This project is developed by [Tao Hou](https://taohou01.github.io) under the [CGTDA](https://www.cs.purdue.edu/homes/tamaldey/CGTDAwebsite/) research group at Purdue University lead by [Dr. Tamal K. Dey](https://www.cs.purdue.edu/homes/tamaldey/).

## About the Implementation

Given a dynamic point cloud as input, the software first builds a sequence of simplex-wise operations on the zigzag filtration, which consist of the following as in `[1]`:

- backward switches
- forward switches
- outward switches
- outward expansions
- inward contractions

It then perform updates on the zigzag barcodes for the operations and builds *vines and vineyards* `[2]` for the input data. More explanations on the implementation are provided in Section [Implementation Details](https://github.com/taohou01/zzup/edit/main/README.md#implementation-details)

## Implementation Details

## References

1. Tamal K. Dey and Tao Hou. Updating Barcodes and Representatives for Zigzag Persistence.
2. David Cohen-Steiner, Herbert Edelsbrunner, and Dmitriy Morozov. Vines and vineyards by
updating persistence in linear time. *The Twenty-Second Annual Symposium
on Computational Geometry*.
