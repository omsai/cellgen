# Cell Images Generator #

Python script to generate test data for [CellProfiler-Analyst](https://github.com/CellProfiler/CellProfiler-Analyst).

Generates 3 channel monochrome cell images with nuclei and centromeres.
Each image varies the amount of overlap between the centromeres.

# Installation #

## Native ##

The Python script depends on scipy and Python Imaging Library (PIL).

It also calls a command-line program, `[fast_delaunay](https://github.com/thouis/fast-poisson-disk)`,
to generate uniform spatial distributions.
To compile `fast_delaunay` one needs the GNU Triangulated Surface Library (libgts).
With libgts installed, one can build `fast_delaunay` using:

``` sh
make -f fast-poisson-disk.makefile
```

Then run the Python script with to generate the `*.tif` images with:

``` sh
python cellgen.py path/to/output/directory
```

Or leave off the path to output the images in the current directory:

``` sh
python cellgen.py
```

## Docker ##

If you prefer to use docker to install the dependencies and run the program,
create the docker image with:

``` sh
docker build -t cellgen .
```

Then run the program with:

``` sh
docker run -v /path/to/output/directory:/data -it cellgen 
```

To create the `*.tif` images in the current directory:

``` sh
docker run -v $PWD:/data -it cellgen 
```
