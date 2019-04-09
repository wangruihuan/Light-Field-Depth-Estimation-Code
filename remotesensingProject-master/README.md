# remotesensingProject


This repository implements some depth map estimation algorithm using 3D light fields, originally proposed by Kim et al. in [Scene Reconstruction from High Spatio-Angular Resolution Light Fields](https://www.disneyresearch.com/publication/scene-reconstruction-from-high-spatio-angular-resolution-light-fields/) (2013). 


This program was made in the framework of the course [Remote sensing: from sensor to large-scale geospatial data exploitation](https://mvaisat.wp.imt.fr/) of the Master 2 MVA (ENS Paris-Saclay, Telecom ParisTech).


The Doxygen documentation for this project can be found in the branch `gh-pages` and is published under [this website](https://14chanwa.github.io/remotesensingProject/).


<p align="center">
<img src="https://raw.githubusercontent.com/14chanwa/remotesensingProject/master/report/images/SkysatLR18_240_img/1521805051081_dmap_050.png" width="800">
</p>
<p align="center"><em>Some sample satellite image and the computed disparity map</em></p>


## Dependancies


* OpenCV 3.x should be installed and findable.

* The program should be compiled with C++11 standards. In particular, this program makes use of `<experimental/filesystem>` and its corresponding library `stdc++fs`.

* OpenMP.


## Minimal working example (with the right paths)


Assume we are in the folder containing `README.md`. The following commands builds the library and the tests and run the library's "Hello World!".

```
mkdir build
cd build
cmake ../RSLightFields
make
./test_read_tiff 0
./test_read_tiff 1
```

