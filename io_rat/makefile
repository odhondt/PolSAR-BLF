###############################################################################
###############################################################################
####                           DEFINE PARAMETERS                           ####
###############################################################################
###############################################################################

# PREFIX is the location of the CSAR application, where binaries will be stored
PREFIX = ./.


###############################################################################
###############################################################################
#####                        DO NOT CHANGE FROM HERE                       ####
###############################################################################
###############################################################################
CXX = g++
BIN = $(PREFIX)
LIB = $(PREFIX)/lib

CIMG_FILE     = CImg-1.5.1
EIGEN_FILE    = eigen-eigen-43d9075b23ef
FFTW3_FILE    = fftw-3.3.2
GDAL_FILE     = gdal-1.9.1
RAPIDXML_FILE = rapidxml-1.13
TRIANGLE_FILE = triangle

#--------------------------------------------------
# Set compilation flags allowing to customize CImg
#--------------------------------------------------
X11PATH = /usr/X11R6

# Flags to enable strict code standards
#CIMG_ANSI_CFLAGS = -ansi -ffast-math -pedantic
CIMG_ANSI_CFLAGS = -ffast-math

# Flags to enable code debugging.
#CIMG_DEBUG_CFLAGS = 
CIMG_DEBUG_CFLAGS = -O3 -Dcimg_verbosity=3 -Dcimg_strict_warnings #-g


# Flags to enable OpenMP support.
#CIMG_OPENMP_CFLAGS = 
CIMG_OPENMP_CFLAGS = -Dcimg_use_openmp -fopenmp

# Flags to enable OpenCV support.
CIMG_OPENCV_CFLAGS = 
CIMG_OPENCV_LDFLAGS = 
#CIMG_OPENCV_CFLAGS = -Dcimg_use_opencv -I/usr/include/opencv
#CIMG_OPENCV_LDFLAGS = -lcv -lhighgui

# Flags to enable the use of the X11 library.
# (X11 is used by CImg to handle display windows)
# !!! For 64bits systems : replace -L$(X11PATH)/lib by -L$(X11PATH)/lib64 !!!
CIMG_X11_CFLAGS = -I$(X11PATH)/include
CIMG_X11_LDFLAGS = -L$(X11PATH)/lib -lpthread -lX11

# Flags to enable fast image display, using the XSHM library (when using X11).
# !!! Seems to randomly crash when used on MacOSX, so do not use it on MacOSX !!!
# remove it for linux system
CIMG_XSHM_CFLAGS = 
CIMG_XSHM_LDFLAGS = 
#CIMG_XSHM_CFLAGS = -Dcimg_use_xshm
#CIMG_XSHM_LDFLAGS = -lXext

# Flags to enable GDI32 display (Windows native).
CIMG_GDI32_CFLAGS =
CIMG_GDI32_LDFLAGS =
#CIMG_GDI32_CFLAGS = -mwindows
#CIMG_GDI32_LDFLAGS = -lgdi32

# Flags to enable screen mode switching, using the XRandr library (when using X11).
# ( http://www.x.org/wiki/Projects/XRandR )
# !!! Not supported by the X11 server on MacOSX, so do not use it on MacOSX !!!
CIMG_XRANDR_CFLAGS = 
CIMG_XRANDR_LDFLAGS = 
#CIMG_XRANDR_CFLAGS = -Dcimg_use_xrandr
#CIMG_XRANDR_LDFLAGS = -lXrandr

# Flags to enable native support for PNG image files, using the PNG library.
# ( http://www.libpng.org/ )
#CIMG_PNG_CFLAGS = 
#CIMG_PNG_LDFLAGS = 
CIMG_PNG_CFLAGS = -Dcimg_use_png
CIMG_PNG_LDFLAGS = -lpng -lz

# Flags to enable native support for JPEG image files, using the JPEG library.
# ( http://www.ijg.org/ )
#CIMG_JPEG_CFLAGS = 
#CIMG_JPEG_LDFLAGS = 
CIMG_JPEG_CFLAGS = -Dcimg_use_jpeg
CIMG_JPEG_LDFLAGS = -ljpeg

# Flags to enable native support for TIFF image files, using the TIFF library.
# ( http://www.libtiff.org/ )
#CIMG_TIFF_CFLAGS = 
#CIMG_TIFF_LDFLAGS = 
CIMG_TIFF_CFLAGS = -Dcimg_use_tiff
CIMG_TIFF_LDFLAGS = -ltiff

# Flags to enable native support for MINC2 image files, using the MINC2 library.
# ( http://en.wikibooks.org/wiki/MINC/Reference/MINC2.0_Users_Guide )
CIMG_MINC2_CFLAGS =
CIMG_MINC2_LDFLAGS =
#CIMG_MINC2_CFLAGS = -Dcimg_use_minc2
#CIMG_MINC2_LDFLAGS = -lminc2

# Flags to enable native support for EXR image files, using the OpenEXR library.
# ( http://www.openexr.com/ )
CIMG_EXR_CFLAGS =
CIMG_EXR_LDFLAGS =
#CIMG_EXR_CFLAGS = -Dcimg_use_openexr -I/usr/include/OpenEXR
#CIMG_EXR_LDFLAGS = -lIlmImf -lHalf

# Flags to enable native support for various video files, using the FFMPEG library.
# ( http://www.ffmpeg.org/ )
CIMG_FFMPEG_CFLAGS =
CIMG_FFMPEG_LDFLAGS = 
#CIMG_FFMPEG_CFLAGS = -Dcimg_use_ffmpeg -D__STDC_CONSTANT_MACROS -I/usr/include/libavcodec -I/usr/include/libavformat -I/usr/include/libswscale -I/usr/include/ffmpeg
#CIMG_FFMPEG_LDFLAGS = -lavcodec -lavformat -lswscale

# Flags to enable native support for compressed .cimgz files, using the Zlib library.
# ( http://www.zlib.net/ )
# CIMG_ZLIB_CFLAGS = 
# CIMG_ZLIB_LDFLAGS = 
CIMG_ZLIB_CFLAGS = -Dcimg_use_zlib
CIMG_ZLIB_LDFLAGS = -lz

# Flags to enable native support of most classical image file formats, using the Magick++ library.
# ( http://www.imagemagick.org/Magick++/ )
CIMG_MAGICK_CFLAGS =
CIMG_MAGICK_LDFLAGS =
#CIMG_MAGICK_CFLAGS = -Dcimg_use_magick `Magick++-config --cppflags` `Magick++-config --cxxflags`
#CIMG_MAGICK_LDFLAGS = `Magick++-config --ldflags` `Magick++-config --libs`

# Flags to enable faster Discrete Fourier Transform computation, using the FFTW3 library
# ( http://www.fftw.org/ )
CIMG_FFTW3_CFLAGS = 
#CIMG_FFTW3_CFLAGS = -Dcimg_use_fftw3 -I$(LIB)/fftw3/include
#ifeq ($(MSYSTEM),MINGW32)
#CIMG_FFTW3_LDFLAGS = -lfftw3-3
#else
 CIMG_FFTW3_LDFLAGS =
#CIMG_FFTW3_LDFLAGS = -L$(LIB)/fftw3/lib -lfftw3
#endif

# Flags to enable the use of LAPACK routines for matrix computation
# ( http://www.netlib.org/lapack/ )
CIMG_LAPACK_CFLAGS =
CIMG_LAPACK_LDFLAGS =
#CIMG_LAPACK_CFLAGS = -Dcimg_use_lapack
#CIMG_LAPACK_LDFLAGS = -lblas -lg2c -llapack

# Flags to enable the use of GDAL routines for input/output file.
# ( http://www.gdal.org )
# be sure that GDAL is installed in the CSAR/lib/GDAL.
CSAR_GDAL_CFLAGS =
CSAR_GDAL_LDFLAGS =
#CSAR_GDAL_CFLAGS = -Dcsar_use_gdal -I$(LIB)/GDAL/include
#CSAR_GDAL_LDFLAGS = -L$(LIB)/GDAL/lib -lgdal

# Flags to enable the use of EIGEN routines
# be sure that eigen is installed in the CSAR/lib/EIGEN directory, use the absolute path.
CSAR_EIGEN_CFLAGS = 
#CSAR_EIGEN_CFLAGS = -Dcsar_use_eigen -I$(LIB)/Eigen

# Flags to enable the use of CImg routines (do not comment)
# be sure that CImg is installed in the CSAR/lib/CIMG directory (at least the CImg.h)
#CSAR_CIMG_CFLAGS = -I$(LIB)/CImg 
CSAR_CIMG_CFLAGS = 

# Flags to enable the use of RapidXML routines
# be sure that rapidxml is installed in the CSAR/lib/RAPIDXML directory
CSAR_RAPIDXML_CFLAGS = 
#CSAR_RAPIDXML_CFLAGS = -Dcsar_use_rapidxml -I$(LIB)/rapidxml

CSAR_TRIANGLE_CFLAGS = 
CSAR_TRIANGLE_LDFLAGS = 
#CSAR_TRIANGLE_CFLAGS = -Dcsar_use_triangle  -DTRILIBRARY -DANSI_DECLARATORS -I$(LIB)/triangle
#CSAR_TRIANGLE_LDFLAGS = $(LIB)/triangle/triangle.o

# Flags to enable the use of csar.h, csar_plugin_cimg.h ... (do not comment)
CSAR_CFLAGS = -I. 
CSAR_LDFLAGS = -lm


CONF_CFLAGS = \
$(CSAR_CFLAGS) \
$(CSAR_CIMG_CFLAGS) \
$(CSAR_EIGEN_CFLAGS) \
$(CSAR_TRIANGLE_CFLAGS) \
$(CSAR_RAPIDXML_CFLAGS) \
$(CSAR_GDAL_CFLAGS) \
$(CIMG_LAPACK_CFLAGS) \
$(CIMG_FFTW3_CFLAGS) \
$(CIMG_MAGICK_CFLAGS) \
$(CIMG_ZLIB_CFLAGS) \
$(CIMG_FFMPEG_CFLAGS) \
$(CIMG_EXR_CFLAGS) \
$(CIMG_MINC2_CFLAGS) \
$(CIMG_TIFF_CFLAGS) \
$(CIMG_JPEG_CFLAGS) \
$(CIMG_PNG_CFLAGS) \
$(CIMG_XRANDR_CFLAGS) \
$(CIMG_GDI32_CFLAGS) \
$(CIMG_XSHM_CFLAGS) \
$(CIMG_X11_CFLAGS) \
$(CIMG_OPENCV_CFLAGS) \
$(CIMG_OPENMP_CFLAGS) \
$(CIMG_DEBUG_CFLAGS) \
$(CIMG_ANSI_CFLAGS)

CONF_LDFLAGS = \
$(CSAR_LDFLAGS) \
$(CSAR_TRIANGLE_LDFLAGS) \
$(CSAR_GDAL_LDFLAGS) \
$(CIMG_LAPACK_LDFLAGS) \
$(CIMG_FFTW3_LDFLAGS) \
$(CIMG_MAGICK_LDFLAGS) \
$(CIMG_ZLIB_LDFLAGS) \
$(CIMG_FFMPEG_LDFLAGS) \
$(CIMG_EXR_LDFLAGS) \
$(CIMG_MINC2_LDFLAGS) \
$(CIMG_TIFF_LDFLAGS) \
$(CIMG_JPEG_LDFLAGS) \
$(CIMG_PNG_LDFLAGS) \
$(CIMG_XRANDR_LDFLAGS) \
$(CIMG_GDI32_LDFLAGS) \
$(CIMG_XSHM_LDFLAGS) \
$(CIMG_X11_LDFLAGS) \
$(CIMG_OPENCV_LDFLAGS)

all: 	_io_

####################################################################################################
#                                           I/O DIRECTORY                                          #
####################################################################################################

_io_:  cimg2rat \
	import_rat

cimg2rat: cimg2rat.cpp
	$(CXX) $(CONF_CFLAGS) cimg2rat.cpp -o cimg2rat $(CONF_LDFLAGS)

import_rat: import_rat.cpp
	$(CXX) $(CONF_CFLAGS) import_rat.cpp -o import_rat $(CONF_LDFLAGS)

