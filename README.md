# Updating zigzag barcodes and representatives for computing vineyards for dynamic point clouds

This project implements the algorithms for computing vines and vineyards for dynamic point clouds using representative-based zigzag update algorihtms. The approaches are described in the following paper:

[Updating Barcodes and Representatives for Zigzag Persistence](https://arxiv.org/pdf/2112.02352.pdf)

by Tamal K. Dey and Tao Hou.

## Group Information

This project is developed by [Tao Hou](https://taohou01.github.io) under the [CGTDA](https://www.cs.purdue.edu/homes/tamaldey/CGTDAwebsite/) research group at Purdue University lead by [Dr. Tamal K. Dey](https://www.cs.purdue.edu/homes/tamaldey/).

## About the Implementation

Given a dynamic point cloud as input, the software first builds a sequence of simplex-wise operations on the zigzag filtration, which consist of the following as in the [paper](https://arxiv.org/pdf/2112.02352.pdf):

- backward switches
- forward switches
- outward switches
- outward expansions
- inward contractions
