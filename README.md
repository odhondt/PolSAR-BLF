## PolSAR-BLF: Iterative Bilateral Filtering of Polarimetric SAR data.

- This C++ code implements the method described in the paper:

[1] D'Hondt, O., Guillaso, S. and Hellwich, O. **Iterative Bilateral Filtering of Polarimetric SAR Data**,
_Selected Topics in Applied Earth Observations and Remote Sensing, IEEE Journal of, 2013, 6, 1628-1639_

- If you use this code in a publication or public presentation, please cite the above reference.

- This code is based on the C++ CImg library http://cimg.sourceforge.net/ and uses Eigen http://eigen.tuxfamily.org for matrix computations.

### Installation

- To get the package you can download the zip file from the project page or you can clone the directory

`git clone https://github.com/odhondt/PolSAR-BLF`

- CImg is included in this distribution as it consists in a single header file `CImg.h`.

- Eigen has to be downloaded from http://eigen.tuxfamily.org/

- To indicate the location of Eigen to the compiler, please edit the file `Makefile` and set the  

- If there is a problem with multi-core, it is possible to disable openmp by simply removing the option -fopenmp.

- The filter can then be compiled by typing `make olinux` (normally `make linux` is sufficient, but it does not enable compiler optimization resulting in a much slower execution.)

### Usage

- To display help and default parameter values, simply type:

`./PolSAR-BLF -h`

- The default parameters will give a result identical to the ones of the publication [1] with affine invariant distance.

- The input format is .cimg which is a binary file with a simple header.

- The simplest way to use our code is to first convert your data to the .rat format with the RADAR-TOOLS software (http://sourceforge.net/projects/radartools.berlios/). Then you can use the functions `import_rat` and `export_rat` that can be found in the directory `io_rat` to convert to and from cimg. Invoking these programs with `-h` option will display help. 

- **IMPORTANT:** The input data must be in the form of a 3x3 polarimetric covariance or coherency matrix with a minimum number of looks equal to the dimension of the matrix so that the matrices are full rank.

- Please note that the display functions are assuming the input is a coherency matrix and show the Pauli RGB base. If using a covariance (lexicographic) basis, please use your own display function. We will implement optional covariance visualization in the future versions.

- If you want to use the CImg data format, you can find information about how the data is stored in memory here: http://cimg.sourceforge.net/reference/group__cimg__storage.html

- For polarimetric covariances storage we use the following convention: the data is stored in a CImgList object which is a list of two images. The first image correspond to the real part and the second image to the imaginary part. To save space, only the upper and diagonal elements are stored: each pixel is a vector with components (Re(C11), Re(C12), Re(C13), Re(C22), Re(C23), Re(C33)) in the first image and (Im(C11), Im(C12), Im(C13), Im(C22), Im(C23), Im(C33)) in the second one. An object list is indexed as follows: `mylist(pos, x, y, z, c)` where `pos` is the image number (0 or 1), `x` and `y` are the pixel spatial coordinates, `z=0` (we deal with 2D images), and `c` the channel varying from 0 to 5.

- The Eigen format is used for matrix computations and the size of the matrices is hard-coded for a faster runtime. That is why for now the code only runs with 3x3 matrices. We will implement functionalities for arbitrary sized matrices in future versions. 

- This software has been tested for linux and mac. 

