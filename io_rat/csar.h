 /*
 #
 #  File        : csar.h
 #
 #  Copyright   : Stéphane Guillaso
 #
 #  License     : CeCILL v2.0
 #                ( http://www.cecill.info/licences/Licence_CeCILL_V2-en.html )
 #
 #  This software is governed by the CeCILL  license under French law and
 #  abiding by the rules of distribution of free software.  You can  use,
 #  modify and/ or redistribute the software under the terms of the CeCILL
 #  license as circulated by CEA, CNRS and INRIA at the following URL
 #  "http://www.cecill.info".
 #
 #  As a counterpart to the access to the source code and  rights to copy,
 #  modify and redistribute granted by the license, users are provided only
 #  with a limited warranty  and the software's author,  the holder of the
 #  economic rights,  and the successive licensors  have only  limited
 #  liability.
 #
 #  In this respect, the user's attention is drawn to the risks associated
 #  with loading,  using,  modifying and/or developing or reproducing the
 #  software by the user in light of its specific status of free software,
 #  that may mean  that it is complicated to manipulate,  and  that  also
 #  therefore means  that it is reserved for developers  and  experienced
 #  professionals having in-depth computer knowledge. Users are therefore
 #  encouraged to load and test the software's suitability as regards their
 #  requirements in conditions enabling the security of their systems and/or
 #  data to be ensured and,  more generally, to use and operate it in the
 #  same conditions as regards security.
 #
 #  The fact that you are presently reading this means that you have had
 #  knowledge of the CeCILL license and that you accept its terms.
 #
*/
/**
 * 
 * \file csar.h
 * \brief C++ template SAR Image Processing Toolkit
 * \author Stéphane Guillaso
 * \version 1.0.0
 * \date 8/09/2012
 *	csar.h is a plugin for the CImg Library that allows to process SAR 
 *	(Synthetic Aperture Radar Images
 *	provided in different format. csar.h ... 

*/

/**
 * \mainpage The CSAR Software
 * \section intro_sec Introduction
 * This is the introduction
 *
 * \note
 * This plugin is only dedicated to CImgList having a specific form: the number
 * of list correspond to store
 * the real/imaginary part.
 * Add example.
 */


/*
 #
 #  File            : csar.h
 #                    ( C++ header file )
 #
 #  Description     : The C++ Template SAR Image Processing Toolkit.
 #                    This file is the main component of the csar software project.
 #                    ( http://www.sourceforge.net/project/c-sar )
 #
 #  Project manager : Stephane Guillaso.
 #                    ( http://www.cv.tu-berlin.de/menue/mitarbeiter_doktoranden/stephane_guillaso/parameter/en/ )
 #
 #                    A complete list of contributors is available in file 'README.txt'
 #                    distributed within the "csar" package.
 #
 #  Licenses        : TO BE DEFINE
 # 
*/
//TODO: Define the license of for the soft.
#ifndef csar_h
#define csar_h


//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
// TODO: reorganize the cimg_for
#define cimg_forCZXY(img,c,z,x,y) cimg_forY(img,y) cimg_forX(img,x) cimg_forZ(img, z) cimg_forC(img,c)
#define cimg_forCXY(img,c,x,y) cimg_forY(img,y) cimg_forX(img,x) cimg_forC(img,c)
#define cimg_forCZ(img,c,z)  cimg_forZ(img,z) cimg_forC(img,c)
#define cimg_forCZXY(img,c,z,x,y) cimg_forY(img,y) cimg_forX(img,x) cimg_forZ(img, z) cimg_forC(img,c)
#define cimg_forCXY(img,c,x,y) cimg_forY(img,y) cimg_forX(img,x) cimg_forC(img,c)
#define cimg_forZXY(img,z,x,y) cimg_forY(img,y) cimg_forX(img,x) cimg_forZ(img,z)
//#define cimg_forCZYX(img,c,z,y,x) cimg_forX(img,x) cimg_forY(img,y) cimg_forZ(img, z) cimg_forC(img,c)
//#define cimg_forCZXY(img,c,z,x,y) cimg_forY(img,y) cimg_forX(img,x) cimg_forZ(img, z) cimg_forC(img,c)
//#define cimg_forCZX(img,c,z,x)  cimg_forX(img,x) cimg_forZ(img, z) cimg_forC(img,c)


// define new fields for handling SAR data and block processing
#define BLOCK_SIZE 256
#define MAX_BLOCK_SIZE 512
#define MAX_FILE_SIZE 52428800
#define DO_TILING 1
#define NO_TILING 0

#define TILE_READ 1
#define TILE_WRITE 2

#define RADEG 57.2958




//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
// field to read/write cimg file
unsigned int _cimg_xdim, _cimg_ydim, _cimg_zdim, _cimg_cdim, _cimg_cplx, _cimg_type;
std::FILE *_cimg_fid;
unsigned long _cimg_seek1, _cimg_seek2;
//std::string _cimg_type;
const char* _cimg_fname;
std::ifstream cimg_ifile;

//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
// field to read/write rat file
unsigned int _rat_xdim, _rat_ydim, _rat_zdim, _rat_vdim, _rat_xdr, _rat_type, _rat_ndim, _rat_cplx, _rat_mode;
std::FILE *_rat_fid;
unsigned long long _rat_header;
bool _rat_lex, _rat_3D;
int _imode;
std::string _mode;

//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
// field to read/write file using gdal library
#ifdef csar_use_gdal
unsigned int _gdal_xdim, _gdal_ydim, _gdal_cdim, _gdal_type, _gdal_cplx;
GDALDataset *_gdal_p;
GDALRasterBand *_gdal_rb;
#endif



//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
// field to use tiling
unsigned long *_bp;
int _ov, *_bs, _bn;



//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
// Overloaded Operators
//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
template<typename t>
CImgList<T>& operator+=(const t val){
  (*this)(0) += val;
  if ((*this).width() == 2) (*this)(1) += val;
  return *this;
}

template<typename t>
CImgList<T> operator+(const t val){
  return CImgList<T>(*this, false)+=val;
}

template<typename t>
CImgList<T>& operator+=(const CImgList<t>& img) {
  cimglist_for(*this,l){
    (*this)(l) += img(l);
  }
  return *this;
}


template<typename t>
CImgList<T> operator+(const CImgList<t>& img) const {
  return CImgList<T>(*this,false)+=img;
}


template<typename t>
CImgList<T>& operator-=(const CImgList<t>& img) {
  cimglist_for(*this,l){
    (*this)(l) -= img(l);
  }
  return *this;
}

template<typename t>
CImgList<T>& operator-=(const t val){
  (*this)(0) -= val;
  if ((*this).width() == 2) (*this)(1) -= val;
  return *this;
}

template<typename t>
CImgList<T> operator-(const t val){
  return CImgList<T>(*this, false)-=val;
}

template<typename t>
CImgList<T> operator-(const CImgList<t>& img) const {
  return CImgList<T>(*this,false)-=img;
}


template<typename t>
CImgList<T>& operator*=(const t val){
  cimglist_for(*this, l) (*this)(0) *= val;
  return *this;
}

template<typename t>
CImgList<T> operator*(const t val) const {
  return CImgList<T>(*this, false)*=val;
}

template<typename t>
CImgList<T>& operator/=(const t val){
  cimglist_for(*this, l) (*this)(0) /= val;
  return *this;
}

template<typename t>
CImgList<T> operator/(const t val) const {
  return CImgList<T>(*this, false) /= val;
}

//! In-place multiplicator operator (for list of real and complex images)
template<typename t>
CImgList<T>& operator*=(const CImgList<t>& img){
  if ((*this).size() == 1){
    CImgList<T> res(img);
    res(0) = (*this)(0).get_mul(img(0));
    *this = res;
    return *this;
  }
  if ((*this).size() > 2 || img.size() > 2){
    std::cerr << "Only for real and complex data" << std::endl;
    exit(1);
  }
  if (img(0).width() == 1 && img(0).height() == 1 && img(0).depth() == 1 && img(0).spectrum() == 1) {
    CImgList<T> res(*this);
    res(0) = (*this)(0) * img(0, 0, 0, 0, 0) - (*this)(1) * img(1, 0, 0, 0, 0);
    res(1) = (*this)(1) * img(0, 0, 0, 0, 0) + (*this)(0) * img(1, 0, 0, 0, 0);
    *this = res;
    return *this;
  }
  CImgList<T> res(img);
  res(0) = (*this)(0).get_mul(img(0)) - (*this)(1).get_mul(img(1));
  res(1) = (*this)(1).get_mul(img(0)) + (*this)(0).get_mul(img(1));
  *this = res;
  return *this;
}

template<typename t>
CImgList<T> operator*(const CImgList<t>& img) const{
  return CImgList<T>(*this,false)*=img;
}


//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
/** @name ARRAY CREATION */
#pragma mark -
#pragma mark -
#pragma mark ARRAY CREATION

//@{
#pragma mark - INDGEN
//! The \c indgen functions returns real array  which is set to the value of its one-dimensional subscript
//! If complex is set to true, then create complex array with imaginary part sets to zero
CImgList<T>& indgen(long dx, bool complex = false){
  return (*this).indgen(dx, 1, 1, 1, complex);
}
CImgList<T> get_indgen(long dx, bool complex = false){
  return CImgList<T>(*this, false).indgen(dx, complex);
}

CImgList<T>& indgen(long dx, long dy, bool complex = false){
  return (*this).indgen(dx, dy, 1, 1, complex);
}
CImgList<T> get_indgen(long dx, long dy, bool complex = false){
  return CImgList<T>(*this, false).indgen(dx, dy, complex);
}

CImgList<T>& indgen(long dx, long dy, long dz, bool complex = false){
  return (*this).indgen(dx, dy, dz, 1, complex);
}
CImgList<T> get_indgen(long dx, long dy, long dz, bool complex = false){
  return CImgList<T>(*this, false).indgen(dx, dy, dz, complex);
}

CImgList<T>& indgen(long dx, long dy, long dz, long dc, bool complex = false){
  CImgList<T> res(1, dx, dy, dz, dc);
  if (!complex) res.assign(1, dx, dy, dz, dc); else res.assign(2, dx, dy, dz, dc);
#pragma omp parallel for
	cimg_forXYZC(res(0), x, y, z, c) res(0, x, y, z, c) = x + dx * y + dx * dy * z + dx * dy * dz * c;
  *this = res;
  return *this;
}

CImgList<T> get_indgen(long dx, long dy, long dz, long dc, bool complex = false){
  return CImgList<T>(*this, false).indgen(dx, dy, dz, dc, complex);
}
//@}


//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
/** @name ARRAY MANIPULATION */
#pragma mark -
#pragma mark -
#pragma mark ARRAY MANIPULATION
//@{

//! Get type of CSAR data
int get_type(const char *fname){
	
}
	
//! Extract channel
CImgList<T>& channel(int ch){
	CImgList<T> res;
	res.assign((*this).width(), (*this)(0).width(), (*this)(0).height());
	res(0) = (*this)(0).get_channel(ch);
	if ((*this).width() == 2) res(1) = (*this)(1).get_channel(ch);
	*this = res;
	return *this;
}

CImgList<T> get_channel(int ch){
	return CImgList<T>(*this, false).channel(ch);
}

#ifdef CSAR_INT_ONE_PARAM
void get_channel(CImgList<T> **i, CImgList<T> **o, void *pP){
	int_one_param_struct *p = static_cast<int_one_param_struct *>(pP);
	*o[0] = i[0]->get_channel((*p).p);
}
#endif

//@}

//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
//! @name IMAGE PROCESSING
#pragma mark -
#pragma mark -
#pragma mark IMAGE PROCESSING
//@{
/** @name Image Geometry Transformations */
//@{
#pragma mark * Image Geometry Transformations
#pragma mark - REBIN

CImgList<T>& rebin(long mx, bool sample = false){
  return (*this).rebin(mx, 1, 1, 1, sample);
}

CImgList<T> get_rebin(long mx, bool sample = false){
  return CImgList<T>(*this,false).rebin(mx, sample);
}

CImgList<T>& rebin(long mx, long my, bool sample = false){
  return (*this).rebin(mx, my, 1, 1, sample);
}

CImgList<T> get_rebin(long mx, long my, bool sample = false){
  return CImgList<T>(*this,false).rebin(mx, my, sample);
}

CImgList<T>& rebin(long mx, long my, long mz, bool sample = false){
  return (*this).rebin(mx, my, mz, 1, sample);
}

CImgList<T> get_rebin(long mx, long my, long mz, bool sample = false){
  return CImgList<T>(*this,false).rebin(mx, my, mz, sample);
}

//! rebin (like in idl, use bilinear/ neighborhood averaging or nearest neighbor sampling (sample = true)
/**
 Return the resized array or vector of the specified dimensions (only 1D and 2D)
 \param mx: new X-axis dimension (should be integer multiples or factor of original dimensions)
 \param my: new Y-axis dimension (should be integer multiples or factor of original dimensions)
 \param mz: new Z-axis dimension (should be integer multiples or factor of original dimensions)
 \param mc: new C-axis dimension (should be integer multiples or factor of original dimensions)
 \param sample: use nearest neighbor instead of bilinear/neighborhood averaging.
 */
CImgList<T>& rebin(long mx, long my, long mz, long mc, bool sample = false){
  
  long nx = (*this)(0).width();
  long ny = (*this)(0).height();
  long nz = (*this)(0).depth();
  long nc = (*this)(0).spectrum();
  
  // test factor
  int factor = nx < mx ? mx % nx : nx % mx;
  if (factor != 0){
    std::cerr << "Result dimensions must be integer factor of original dimensions" << std::endl;
    exit(1);
  }
  factor = ny < my ? my % ny : ny % my;
  if (factor != 0){
    std::cerr << "Result dimensions must be integer factor of original dimensions" << std::endl;
    exit(1);
  }
  factor = nz < mz ? mz % nz : nz % mz;
  if (factor != 0){
    std::cerr << "Result dimensions must be integer factor of original dimensions" << std::endl;
    exit(1);
  }
  factor = nc < mc ? mc % nc : nc % mc;
  if (factor != 0){
    std::cerr << "Result dimensions must be integer factor of original dimensions" << std::endl;
    exit(1);
  }
  
  CImgList<T> res;
  int CPLX = (*this).width();
  
  // case 'x'
  if (nx < mx) { // expensions
    double f = (double) nx / (double) mx;
    res.assign(CPLX, mx, ny, nz, nc);
    cimg_forXYZC(res(0), kx, ky, kz, kc){
      double p = f * (double) kx;
      T x_p  = (*this)(0, (long) std::floor(p), ky, kz, kc);
      T x_p1 = p < (nx - 1)? (*this)(0, (long) std::floor(p) +1, ky, kz, kc) : x_p;
      res(0, kx, ky, kz, kc) = sample ? x_p : (T) x_p + (p - std::floor(p)) * (x_p1 - x_p);
      if (CPLX == 2) {
        x_p  = (*this)(1, (long) std::floor(p), ky, kz, kc);
        x_p1 = p < (nx - 1)? (*this)(1, (long) std::floor(p) +1, ky, kz, kc) : x_p;
        res(1, kx, ky, kz, kc) = sample ? x_p : (T) x_p + (p - std::floor(p)) * (x_p1 - x_p);        
      }
    }
  }
  if (nx > mx) { // compression
    double f = (double) mx / (double) nx;
    long f1 = (long) 1./f;
    res.assign(CPLX, mx, ny, nz, nc);
    cimg_forXYZC(res(0), kx, ky, kz, kc){
      T x_p = (T) 0;
      for (int k = 0; k < f1; k++) {
        x_p += (*this)(0, f1 * kx + k , ky, kz, kc) / (T) f1;
      }
      res(0, kx, ky, kz, kc) = sample ? (*this)(0, f1 * kx, ky, kz, kc) : x_p;
      if (CPLX == 2){
        x_p = (T) 0;
        for (int k = 0; k < f1; k++) {
          x_p += (*this)(1, f1 * kx + k , ky, kz, kc) / (T) f1;
        }
        res(1, kx, ky, kz, kc) = sample ? (*this)(1, f1 * kx, ky, kz, kc) : x_p;
      }
    }
  }
  if (mx == nx) {
    res = *this;
  }
  *this = res; nx = (*this)(0).width();
  
  // case 'y'
  if (ny < my) { // expensions
    double f = (double) ny / (double) my;
    res.assign(CPLX, nx, my, nz, nc);
    cimg_forXYZC(res(0), kx, ky, kz, kc){
      double p = f * (double) ky;
      T x_p  = (*this)(0, kx, (long) std::floor(p), kz, kc);
      T x_p1 = p < (ny - 1)? (*this)(0, kx, (long) std::floor(p) +1, kz, kc) : x_p;
      res(0, kx, ky, kz, kc) = sample ? x_p : (T) x_p + (p - std::floor(p)) * (x_p1 - x_p);
      if (CPLX == 2) {
        x_p  = (*this)(1, kx, (long) std::floor(p), kz, kc);
        x_p1 = p < (ny - 1)? (*this)(1, kx, (long) std::floor(p) +1, kz, kc) : x_p;
        res(1, kx, ky, kz, kc) = sample ? x_p : (T) x_p + (p - std::floor(p)) * (x_p1 - x_p);        
      }
    }
  }
  if (ny > my) { // compression
    double f = (double) my / (double) ny;
    long f1 = (long) 1./f;
    res.assign(CPLX, nx, my, nz, nc);
    cimg_forXYZC(res(0), kx, ky, kz, kc){
      T x_p = (T) 0;
      for (int k = 0; k < f1; k++) {
        x_p += (*this)(0, kx , f1 * ky + k, kz, kc) / (T) f1;
      }
      res(0, kx, ky, kz, kc) = sample ? (*this)(0, kx, f1 * ky, kz, kc) : x_p;
      if (CPLX == 2){
        x_p = (T) 0;
        for (int k = 0; k < f1; k++) {
          x_p += (*this)(1, kx , f1 * ky + k, kz, kc) / (T) f1;
        }
        res(1, kx, ky, kz, kc) = sample ? (*this)(1, kx, f1 * ky, kz, kc) : x_p;
      }
    }
  }
  if (my == ny) {
    res = *this;
  }
  *this = res; ny = (*this)(0).height();
  
  // case 'z'
  if (nz < mz) { // expensions
    double f = (double) nz / (double) mz;
    res.assign(CPLX, nx, ny, mz, nc);
    cimg_forXYZC(res(0), kx, ky, kz, kc){
      double p = f * (double) kz;
      T x_p  = (*this)(0, kx, ky, (long) std::floor(p), kc);
      T x_p1 = p < (nz - 1)? (*this)(0, kx, ky, (long) std::floor(p) +1, kc) : x_p;
      res(0, kx, ky, kz, kc) = sample ? x_p : (T) x_p + (p - std::floor(p)) * (x_p1 - x_p);
      if (CPLX == 2) {
        x_p  = (*this)(1, kx, ky, (long) std::floor(p), kc);
        x_p1 = p < (nz - 1)? (*this)(1, kx, ky, (long) std::floor(p) +1, kc) : x_p;
        res(1, kx, ky, kz, kc) = sample ? x_p : (T) x_p + (p - std::floor(p)) * (x_p1 - x_p);        
      }
    }
  }
  if (nz > mz) { // compression
    double f = (double) mz / (double) nz;
    long f1 = (long) 1./f;
    res.assign(CPLX, nx, ny, mz, nc);
    cimg_forXYZC(res(0), kx, ky, kz, kc){
      T x_p = (T) 0;
      for (int k = 0; k < f1; k++) {
        x_p += (*this)(0, kx , ky, f1 * kz + k, kc) / (T) f1;
      }
      res(0, kx, ky, kz, kc) = sample ? (*this)(0,  kx, ky, f1 * kz, kc) : x_p;
      if (CPLX == 2){
        x_p = (T) 0;
        for (int k = 0; k < f1; k++) {
          x_p += (*this)(1,  kx , ky, f1 * kz + k, kc) / (T) f1;
        }
        res(1, kx, ky, kz, kc) = sample ? (*this)(1,  kx, ky, f1 *kz, kc) : x_p;
      }
    }
  }
  if (mz == nz) {
    res = *this;
  }
  *this = res; nz = (*this)(0).depth();
  
  // case 'c'
  if (nc < mc) { // expensions
    double f = (double) nc / (double) mc;
    res.assign(CPLX, nx, ny, nz, mc);
    cimg_forXYZC(res(0), kx, ky, kz, kc){
      double p = f * (double) kc;
      T x_p  = (*this)(0, kx, ky, kz, (long) std::floor(p));
      T x_p1 = p < (nc - 1)? (*this)(0, kx,  ky, kz, (long) std::floor(p) +1) : x_p;
      res(0, kx, ky, kz, kc) = sample ? x_p : (T) x_p + (p - std::floor(p)) * (x_p1 - x_p);
      if (CPLX == 2) {
        x_p  = (*this)(1, kx, ky, kz, (long) std::floor(p));
        x_p1 = p < (nc - 1)? (*this)(1, kx, ky, kz, (long) std::floor(p) +1) : x_p;
        res(1, kx, ky, kz, kc) = sample ? x_p : (T) x_p + (p - std::floor(p)) * (x_p1 - x_p);        
      }
    }
  }
  if (nc > mc) { // compression
    double f = (double) mc / (double) nc;
    long f1 = (long) 1./f;
    res.assign(CPLX, nx, ny, nz, mc);
    cimg_forXYZC(res(0), kx, ky, kz, kc){
      T x_p = (T) 0;
      for (int k = 0; k < f1; k++) {
        x_p += (*this)(0, kx, ky, kz, f1 * kc + k ) / (T) f1;
      }
      res(0, kx, ky, kz, kc) = sample ? (*this)(0, kx, ky, kz, f1 * kc) : x_p;
      if (CPLX == 2){
        x_p = (T) 0;
        for (int k = 0; k < f1; k++) {
          x_p += (*this)(1, kx , ky, kz, f1 * kc + k) / (T) f1;
        }
        res(1, kx, ky, kz, kc) = sample ? (*this)(1, kx, ky, kz, f1 * kc) : x_p;
      }
    }
  }
  if (mc == nc) {
    res = *this;
  }
  *this = res;
  return *this;  
}

/*CImgList<T>& rebin(long m, bool sample = false){
  
  // test if 1D 
  if ((*this)(0).height() > 1 || (*this)(0).depth() > 1 || (*this)(0).spectrum() > 1){
    std::cerr << "only 1 dimension, see help for other rebin function" << std::endl;
    exit(1);
  }
  long n = (*this)(0).width();
  
  if (n == m) return *this;
  
  // test if expension or compressing
  bool expension = false;
  if (n < m) expension = true;
  CImgList<T> res;
  
  // case expension
  if (expension){
    if ((m % n) != 0) {
      std::cerr << "Result dimensions must be integer factor of original dimensions" << std::endl;
      exit(1);
    }
    double f = (double) n / (double) m;
    res.assign((*this).width(), m);
    cimg_forX(res(0), i){
      double p = f * (double) i;
      T x_p = (*this)(0, (long) std::floor(p));
      T x_p1 = p < (n-1) ? (*this)(0, (long) std::floor(p) + 1) : x_p;
      res(0, i) = sample ? x_p : (T) x_p + (p - std::floor(p)) * (x_p1 - x_p);
      if ((*this).width() == 2) {
        x_p = (*this)(1, (long) std::floor(p));
        x_p1 = p < (n-1) ? (*this)(0, (long) std::floor(p) + 1) : x_p;
        res(1, i) = sample ? x_p : (T) x_p + (p - std::floor(p)) * (x_p1 - x_p);
      }
    }
  }
  *this = res;
  return *this;
  
}
CImgList<T> get_rebin(long m, bool sample = false){
  return CImgList<T>(*this, false).rebin(m,sample);
}*/
//@}
//@}


//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
//! \name Mathematical Functions
//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
#pragma mark -
#pragma mark -
#pragma mark MATHEMATICAL FUNCTIONS

#pragma mark - SQR Function
//! Compute the square value of each pixel (takes into account the complex case)
CImgList<T>& sqr(){
  CImgList<T> res(*this);
  res = (*this) * (*this);
  *this = res;
  return *this;
}
CImgList<T> get_sqr(){
  return CImgList<T>(*this, false).sqr();
}

#pragma mark - SQRT Function
//! Compute the square root value of each pixel (only for real image)
CImgList<T>& sqrt(){
  if ((*this).size() > 1) {
    std::cerr << "Only for real" << std::endl;
    exit(1);
  }
  CImgList<T> res(*this);
  res(0) = (*this)(0).get_sqrt();
  *this = res;
  return *this;
}
CImgList<T> get_sqrt(){
  return CImgList<T>(*this,false).sqrt();
}

#pragma mark - ATAN Function
//! Compute the arctangent of each pixel value, if complex, compute the argument of a complex number.
CImgList<T>& atan(){
  CImgList<T> res(1);
  if ((*this).size() == 1) res(0) = (*this)(0).get_atan();
  else                     res(0) = (*this)(1).get_atan2((*this)(0));
  *this = res;
  return *this;
}
CImgList<T> get_atan(){
  return CImgList<T>(*this, false).atan();
}

/*#pragma mark - HADAMAR Function
//! Compute the Hadamar product.
CImgList<T>& hadamar(CImgList<T> img){
  
}
*/

#pragma mark - Real 2 Complex Phase
//! Transform a real list into a complex phase of the form: \c res = exp(complex(0, ph));
CImgList<T>& phase2complex(){
  if ((*this).size() >= 2) {
    std::cerr << "Work only for real CImgList image" << std::endl;
    exit(1);
  }
  CImgList<T> res(2, (*this)(0).width(), (*this)(0).height(), (*this)(0).depth(), (*this)(0).spectrum());
  res(0) = (*this)(0).get_cos();
  res(1) = (*this)(0).get_sin();
  *this = res;
  return *this;
}
CImgList<T> get_phase2complex(){
  return CImgList<T>(*this, false).phase2complex();
}

//! Transform a real list into a complex phase of the form: \c res = exp(complex(0, ph));
CImgList<T>& complexphase(){
  if ((*this).size() >= 2) {
    std::cerr << "Work only for real CImgList image" << std::endl;
    exit(1);
  }
  CImgList<T> res(2, (*this)(0).width(), (*this)(0).height(), (*this)(0).depth(), (*this)(0).spectrum());
  res(0) = (*this)(0).get_cos();
  res(1) = (*this)(0).get_sin();
  *this = res;
  return *this;
}
CImgList<T> get_complexphase(){
  return CImgList<T>(*this, false).complexphase();
}

#pragma mark - Conjugate of complex number
//! Compute the conjugate of a complex number
CImgList<T>& conj(){
  if ((*this).size() == 2)
    (*this)(1) = - (*this)(1);
  return *this;
}
CImgList<T> get_conj(){
  return CImgList<T>(*this, false).conj();
}

#pragma mark - Average [MEAN]
//! Return the average value as CImgList complex
CImgList<double> mean(){
  CImgList<double> res((*this).width(), 1, 1, 1, 1);
  res(0, 0, 0, 0, 0) = (*this)(0).mean();
  if((*this).width() == 2) res(1, 0, 0, 0, 0) = (*this)(1).mean();
  return res;
}

double kurtosis(){
  long X = (*this)(0).width();
  long Y = (*this)(0).height();
  long Z = (*this)(0).depth();
  long C = (*this)(0).spectrum();
  CImg<T> tmp(X, Y, Z, C), tmp2;
  cimg_forXYZC((*this)(0), x, y, z, c) tmp(x,y,z,c) = (*this)(0, x, y, z, c);
  double mean = tmp.mean();
  double std  = std::sqrt(tmp.variance());
  tmp = (tmp - mean) / std; 
  tmp.sqr().sqr();
  return tmp.mean() - 3;
}

#pragma mark - Median [MEDIAN]
//! Return the median value as T variable. For complex values, do only to the real part
T median() {
  return (*this).median((*this)(0).max()+ (T) 1);
}
T median(T missing){
  T res;
  unsigned long nr = 0;
  cimg_forXYZC( (*this)(0), x, y, z, c) if ((*this)(0, x, y, z, c) < missing) nr++;
  CImg<T> med(nr);
  nr = 0;
  cimg_forXYZC((*this)(0), x, y, z, c){
    if ((*this)(0, x, y, z, c) < missing) {
      med(nr++) = (*this)(0, x, y, z, c);
    }
  }
  return med.median();
}

//! Return statistic of the image: min, max, mean, stddev
/**
 * im(0, 0, 0, 0, 0) = min for the band 0
 * im(1, 0, 0, 0, 1) = max for the second band
 * im(2, 0, 0, 0, 2) = mean ...
 * im(3, 0, 0, 0, 0) = stddev ...
 *
 * Only in tile processing, image should be opened using im.load_cimg_header(fname);
 */
void get_image_stat(CImgList<T>& o){
  o.assign(1, 4, 1, 1, _cimg_cdim, 0.0);
  get_tiles();
  double av = 0., mmin, mmax;
  for (int k=0; k<get_bn(); k++) {
    std::cout << k << " step00" << std::endl;
    read_cimg(k);
    std::cout << k << " step01" << std::endl;
    CImg<T> tmp((*this)(0).width(), (*this)(0).height());
    std::cout << k << " step1" << std::endl;
    tmp = ((*this)(0).get_sqr() + (*this)(1).get_sqr()).get_sqrt();
    std::cout << k << " step2" << std::endl;
    o(0, 2, 0, 0, 0) += tmp.sum() / (double) _cimg_xdim / (double) _cimg_ydim;
    std::cout << k << " step3" << std::endl;
    Tdouble mmin, mmax;
    std::cout << k << " step4" << std::endl;
    mmin = (Tdouble) tmp.min_max(mmax);
    std::cout << k << " step5" << std::endl;
    if (k == 0){
      std::cout << k << " step61" << std::endl;
      o(0, 0, 0, 0, 0) = mmin;
      std::cout << k << " step62" << std::endl;
      o(0, 1, 0, 0, 0) = mmax;
    }
    else {
      std::cout << k << " step63" << std::endl;
      if (o(0, 0, 0, 0, 0) > mmin) o(0, 0, 0, 0, 0) = mmin;
      std::cout << k << " step64" << std::endl;
      if (o(0, 1, 0, 0, 0) < mmax) o(0, 1, 0, 0, 0) = mmax;
    }

  }
  std::cout << o(0, 0) << " " << o(0, 1) << " " << o(0, 2) << std::endl;
}

Tdouble sum() const{
  Tdouble res;
  res = (*this)(0).sum();
  return res;
}
#pragma mark - Eigen decomposition
//! eigen decomposition
#ifdef csar_use_eigen
//void eig_ord(Eigen::MatrixXcf cov, Eigen::VectorXf v, Eigen::MatrixXcf w){
//  using namespace Eigen;
//  SelfAdjointEigenSolver<MatrixXcf> es;
//  es.compute(cov);
//  v = es.eigenvectors();
//  std::cout << v << std::endl;
//}
#endif

#pragma mark - Poly Fit
//! Poly fit
CImgList<T>& polyfit(CImg<T> x, CImg<T> y, int deg){
  int n = x.width();
  int m = deg + 1;
  // construct work array
  CImgList<T> res(1, m); // define the output coefficient
  CImg<T> A(m,m); // least square matrix
  CImg<T> b(m);   // will contain sum weights*y*x^i
  CImg<T> z(x);
  z.fill(1.0);
  A(0,0) = n;
  b(0) = y.sum();
  for (int p=1; p<= (2*deg); p++){
    cimg_forX(x, xx) z(xx) *= x(xx);
    if (p<m) b(p) = (y.get_mul(z)).sum();
    double sum = z.sum();
    int beg, end;
    beg = (p-deg) > 0 ? p-deg : 0;
    end = deg < p ? deg : p;
    for (int j=beg; j<=end; j++) {
      A(j,p-j) = sum;
    }
  }
  A.invert();
  for (int k=0; k<m; k++) res(0)(k) = (b.get_mul(A.get_column(k))).sum();
  *this = res;
  return *this;
}

//! Poly
CImgList<T> get_poly(CImgList<T> x){
  CImgList<T> res(x);
  int n = (*this)(0).width() - 1;
  if (n==0) {
    cimg_forX(res,xx) res(xx) = (*this)(0)(0);
    return res;
  }
  cimg_forX(res,xx) res(xx) = (*this)(0)(n);
  for (int i=n-1; i == 0; i--) {
    cimg_forX(res, xx) res(xx) = res(xx) * x(xx) + (*this)(0)(i);
  }
  
  return res;
}

#pragma mark - value_locate
CImgList<T> get_value_locate(CImgList<T> value, int inc = 1){
  CImgList<T> res;
  res.assign(1, value(0).width());
  // check if increasing or decreasing. The vector is assumed monotonic.
  int N = (*this)(0).width();
  
  // case increasing
  if (inc == 1){
    int j = 0;
    cimg_forX(res(0), x){
      if (value(0)(x) < (*this)(0)(0)) res(0)(x) = -1;
      else 
        if (value(0)(x) >= (*this)(0)(N-1)) res(0)(x) = N-1;
        else {
          for (int j=0; j<N-1; j++) {
            if ((*this)(0)(j) <= value(0)(x) && value(0)(x) < (*this)(0)(j+1)) {
              res(0)(x) = j;
            }
          }
        }
    } // TODO: decreasing case
  }
  return res;
}

//! function interpol (named after idl) works only for long t1 and t2 and monotonic increasing
CImgList<T> get_interpol(CImgList<long> t1, CImgList<long> t2){
  CImgList<long> s = t1.get_value_locate(t2);
  int N = t1(0).width();
  cimg_forX(s(0), x){
    if (s(0)(x) < 0) s(0)(x) = 0;
    if (s(0)(x) > (N -2)) s(0)(x) = N-2;    
  }
  CImgList<T> p(1, t2(0).width());
  cimg_forX(p(0), x){
    T diff = (*this)(0)(s(0)(x)+1) - (*this)(0)(s(0)(x));
    p(0)(x) = (t2(0)(x) - t1(0)(s(0)(x))) * diff / (t1(0)(s(0)(x)+1) - t1(0)(s(0)(x))) + (*this)(0)(s(0)(x));
  }
  return p;
}

#pragma mark - ZERO PADDING
//! zero padding (use FFT)
/**
 Perform a zero padding (increase or decrease the image size) by inserting (and removing) zero in the 
 Fourier domain, work only for 2D images (real or complex)
 \param f : f > 0 : increase image, f < 0 decrease image
 \param axis : direction -> 'x' for X-axis, 'y' for Y-axis
 */

CImgList<T>& zero_padding(int f, const char axis = 'x'){
  // check dimensions
  if ((*this)(0).depth() != 1 || (*this)(0).spectrum() != 1) {
    std::cerr << "Only 2D image (real or complex)" << std::endl;
    exit(1);
  }
  CImgList<T> res, tmp;
  res.assign(2);
  CImg<T> zero;
  int X = (*this)(0).width();
  int Y = (*this)(0).height();
  int N;
  if (f<0) {
    f = -f;
    switch (cimg::uncase(axis)) {
      case 'x':
        N = X / f;
        tmp = (*this).get_fft_2D(-1,1);
        res(0) = tmp(0).get_columns(0, N/2 - 1);
        res(0).append(tmp(0).get_columns(f * N - N/2, X-1));
        res(1) = tmp(1).get_columns(0, N/2 - 1);
        res(1).append(tmp(1).get_columns(f * N - N/2, X-1));
        res.fft_2D(1, 1);
        break;
      case 'y':
        N = Y /f;
        tmp = (*this).get_fft_2D(-1, 2);
        res(0) = tmp(0).get_rows(0, N/2-1);
        res(0).append(tmp(0).get_rows( f* N - N/2, Y-1), 'y');
        res(1) = tmp(1).get_rows(0, N/2-1);
        res(1).append(tmp(1).get_rows( f* N - N/2, Y-1), 'y');
        res.fft_2D(1, 2);
        break;
    }
  } else {
    switch (cimg::uncase(axis)) {
      case 'x':{ // Along the X-axis.
        zero.assign((f-1) * X, Y).fill((T) 0);
        tmp = (*this).get_fft_2D(-1, 1);
        res(0) = tmp(0).get_columns(0, X/2-1);
        res(0).append(zero);
        res(0).append(tmp(0).get_columns(X/2,X-1));
        res(1) = tmp(1).get_columns(0, X/2-1);
        res(1).append(zero);
        res(1).append(tmp(1).get_columns(X/2,X-1));
        res.fft_2D(1, 1);
      } break;
      case 'y':{ // Along the Y-axis
        zero.assign(X, (f-1) * Y).fill((T) 0);
        tmp = (*this).get_fft_2D(-1, 2);
        res(0) = tmp(0).get_rows(0, Y/2 - 1);
        res(0).append(zero, 'y');
        res(0).append(tmp(0).get_rows(Y/2, Y-1), 'y');
        res(1) = tmp(1).get_rows(0, Y/2 - 1);
        res(1).append(zero, 'y');
        res(1).append(tmp(1).get_rows(Y/2, Y-1), 'y');
        res.fft_2D(1, 2);
      } break;
    }
  }
  *this = res;
  return *this;
}

CImgList<T> get_zero_padding(int f, const char axis = 'x'){
  CImgList<T>(*this, false).zero_padding(f, axis);
}

/*CImgList<T>& zero_padding(int fx, int fy = 1, bool reverse = false){
  // check dimensions
  if ((*this)(0).depth() != 1 || (*this)(0).spectrum() != 1) {
    std::cerr << "Only 2D image (real or complex)" << std::endl;
    exit(1);
  }
  CImgList<T> res, tmp;
  CImg<T> zero;
  res.assign((*this).width());
  // dimensions of original image
  int X = (*this)(0).width();
  int Y = (*this)(0).height();
  int NX = X, NY = Y;
  // case fx
  if (fx > 1){
    zero.assign((fx-1) * X, Y).fill((T)0);
    tmp = (*this).get_fft_2D(-1, 1);
    NX = X * fx; // new image size
    res(0) = tmp(0).get_columns(0, X/2-1);
    res(0).append(zero);
    res(0).append(tmp(0).get_columns(X/2,X-1));
    res(1) = tmp(1).get_columns(0, X/2-1);
    res(1).append(zero);
    res(1).append(tmp(1).get_columns(X/2,X-1));
    res.fft_2D(1, 1);
    *this = res;
    return *this;
  }
  // case fy
  if (fy > 1) {
    exit(1);
  }
}
*/
// zero padding, work only in 1D
/*CImgList<T>& zero_padding(long fx, bool backward = false){
  CImgList<T> res, tmp;
  tmp.assign((*this).width(), (*this)(0).width());
  cimg_forX(tmp(0), x) tmp(0,x) = (*this)(0,x);
  tmp.fft(-1);
  if (!backward){
    res.assign((*this).width(), (*this)(0).width() * fx);
    res(0).fill(0.0); if ((*this).width() == 2) res(1).fill(0.0);
    for (int k=0; k<((*this)(0).width() / 2); k++) {
      res(0, k) = tmp(0, k);
      if (tmp.width() == 2) res(1, k) = tmp(1, k);
    }
    for (int k=(*this)(0).width() / 2; k < (*this)(0).width(); k++) {
      res(0, (over - 1) * (*this)(0).width() + k) = tmp(0, k);
      if (tmp.width() == 2) res(1, (over - 1) * (*this)(0).width() + k) = tmp(1, k);
    }
  } else {
    res.assign((*this).width(), (*this)(0).width() / fx);
    for (int k=0; k < ((*this)(0).width() / fx) / 2; k++) {
      
    }
  }

}*/

#pragma mark - CONGRID
//! congrid (Like in idl, use cubic form)

CImgList<T>& congrid(long nx, long ny, std::string mode = "cubic"){
	return get_congrid(nx, ny, mode).move_to(*this);
}

CImgList<T> get_congrid(long nx, long ny, std::string mode = "cubic"){
  CImg<float> srx, sry;
  float dimx, dimy;
  dimx = (float) (*this)(0).width();
  dimy = (float) (*this)(0).height();
  srx.assign(nx), sry.assign(ny);
  cimg_forX(srx, x) srx(x) = (float) x * dimx / (float) nx;
  cimg_forX(sry, y) sry(y) = (float) y * dimy / (float) ny;
  CImgList<T> res;
  res.assign((*this).width(), nx, ny);
  cimg_forXY(res(0), xx, yy){
    res(0)(xx, yy) = mode == "cubic" ? (*this)(0).cubic_atXY(srx(xx), sry(yy)) : (*this)(0).linear_atXY(srx(xx), sry(yy));
  }
  return res;
}

//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------------------------------


#pragma mark - Interpolate
CImgList<T>& interpolate(CImg<T> x, std::string mode = "linear"){
  return get_interpolate(x,mode).move_to(*this);
}
CImgList<T> get_interpolate(CImg<T> x, std::string mode = "linear"){
  CImgList<T> res;
  res.assign((*this).width(), x.width(), x.height());
#pragma omp parallel for  
	cimg_forX(res(0), xx) {
    res(0)(xx) = mode == "linear" ? (*this)(0).linear_atX(x(xx)) : (*this)(0).cubic_atX(x(xx), 0, 0, 0, 0.0);
    if (res.width() == 2) res(1,xx) = mode == "linear" ? (*this)(1).linear_atX(x(xx)) : (*this)(1).cubic_atX(x(xx),0,0,0,0.0);
  }
  return res;
}

CImgList<T>& interpolate(CImg<T> x, CImg<T> y, std::string mode = "linear"){
  return get_interpolate(x, y, mode).move_to(*this);
}
CImgList<T> get_interpolate(CImg<T> x, CImg<T> y, std::string mode = "linear"){
  CImgList<T> res;
  res.assign((*this).width(), x.width(), x.height());
#pragma omp parallel for
	cimg_forXY(res(0), xx, yy){
    res(0)(xx, yy) = mode == "linear" ? (*this)(0).linear_atXY(x(xx,yy), y(xx,yy), 0, 0, 0.0) : (*this)(0).cubic_atXY(x(xx,yy), y(xx,yy), 0, 0, 0.0);
    if (res.width() == 2) res(1,xx,yy) =  mode == "linear" ? (*this)(1).linear_atXY(x(xx,yy), y(xx,yy), 0, 0, 0.0) : (*this)(1).cubic_atXY(x(xx,yy), y(xx,yy), 0, 0, 0.0);
  }
  return res;
}

#ifdef csar_use_triangle
#pragma mark - TRIANGULATE
//! \c delaunay_triangulate: construct a Delaunay triangulation. Need triangle.o 
// TODO: check input/output
CImgList<T>& delaunay_triangulate(CImgList<T> x, CImgList<T> y){
  struct triangulateio In, Out, vorOut;
  // set everything to 0 and NULL
  std::memset((void*) &In,     0, sizeof(struct triangulateio));
  std::memset((void*) &Out,    0, sizeof(struct triangulateio));
  std::memset((void*) &vorOut, 0, sizeof(struct triangulateio));
  
  In.numberofpoints = x(0).width();
  In.pointlist      = new REAL [2 * x(0).width()];
  int i, j;
  for (i=j=0; i<x(0).width(); i++) { In.pointlist[j++] = x(0)(i); In.pointlist[j++] = y(0)(i);}
  Out.pointlist = (REAL *) NULL;
  vorOut.pointlist = (REAL *) NULL;
  // perform triangulation:
  triangulate("znIQB", &In, &Out, &vorOut);
  // assign triangles:
  CImgList<T> res;
  res.assign(1, Out.numberoftriangles, 3);
  for (i=j=0; i<Out.numberoftriangles; i++) {
    res(0)(i, 0) = Out.trianglelist[j++];
    res(0)(i, 1) = Out.trianglelist[j++];
    res(0)(i, 2) = Out.trianglelist[j++];
  }
  if(In.pointlist) delete[] In.pointlist;
  if(Out.pointlist)    std::free(Out.pointlist);
  if(Out.trianglelist) std::free(Out.trianglelist);
  *this = res;
  return *this;
}

//! \c trigrid (names after IDL function)
// TODO change to be not so close from GDL function
CImgList<T>& trigrid(CImgList<T> x, CImgList<T> y, CImgList<T> z, CImgList<T> tr, double missing = 0.0){
  std::vector <double> gs(2), limits(4);
  limits[0] = x(0).min();
  limits[1] = y(0).min();
  limits[2] = x(0).max();
  limits[3] = y(0).max();
  gs[0]     = (limits[2] - limits[0])/50.;
  gs[1]     = (limits[3] - limits[1])/50.;
  return (*this).trigrid(x, y, z, tr, gs, limits, missing);
}


CImgList<T>& trigrid(CImgList<T> x, CImgList<T> y, CImgList<T> z, CImgList<T> tr, std::vector<double> gs, std::vector<double> limits, double missing = 0.0){
  
  double x_min = limits[0];
  double y_min = limits[1];
  double x_size = limits[2] - limits[0];
  double y_size = limits[3] - limits[1];
  
  
  // determine grid spacing
  double x_step = gs[0];
  double y_step = gs[1];
  double x_dim = x_size / x_step;
  double y_dim = y_size / y_step;
  
  // setup the return array
  CImgList<T> res;
  res.assign(1, x_dim+1, y_dim+1);
  res(0).fill(missing);
  
  for (int k=0; k<tr(0).width(); k++){
    
    int i0 = (int) tr(0)(k, 0);
    int i1 = (int) tr(0)(k, 1);
    int i2 = (int) tr(0)(k, 2);
    
    
    // get the 3 corners of the triangle
    double x0 = x(0)(i0);
    double x1 = x(0)(i1);
    double x2 = x(0)(i2);
    double y0 = y(0)(i0);
    double y1 = y(0)(i1);
    double y2 = y(0)(i2);
    
    // *** plane interpolation *** 
    double diff_x10 = x(0)(i1) - x(0)(i0);
    double diff_x21 = x(0)(i2) - x(0)(i1);
    double diff_y10 = y(0)(i1) - y(0)(i0);
    double diff_y21 = y(0)(i2) - y(0)(i1);
    double diff_z10 = z(0)(i1) - z(0)(i0);
    double diff_z21 = z(0)(i2) - z(0)(i1);
        
    double C = (diff_x21 * diff_z10 - diff_x10 * diff_z21) / (diff_x21 * diff_y10 - diff_x10 * diff_y21);
    double B = (diff_z10 - C * diff_y10) / diff_x10;
    double A = z(0)(i0) - B * x(0)(i0) - C * y(0)(i0);
    
    // sort according increasing 'x'
    if (x0 > x1) {
      cimg::swap(x0, x1);
      cimg::swap(y0, y1);
    }
    if (x1 > x2) {
      cimg::swap(x1, x2);
      cimg::swap(y1, y2);
    }
    if (x0 > x1) {
      cimg::swap(x0, x1);
      cimg::swap(y0, y1);
    }
    
    // determine the 3 line equations
    double a0 = (y1 - y0) / (x1 - x0); double b0 = y1 - a0 * x1;
    double a1 = (y2 - y1) / (x2 - x1); double b1 = y2 - a1 * x2;
    double a2 = (y2 - y0) / (x2 - x0); double b2 = y0 - a2 * x0;
    
    // test the first segment (x0 and x1)
    double x0k = (x0 - x_min)/x_step;
    double x1k = (x1 - x_min)/x_step;
    x0k = (x0k / std::floor(x0k)) == 1 ? std::floor(x0k) : std::floor(x0k) + 1;
    x1k = std::floor(x1k);
    
    for (long kx = (long) x0k; kx <= (long) x1k; kx++) {
      double y0k = a0 * (x_min + kx * x_step) + b0;
      double y1k = a2 * (x_min + kx * x_step) + b2;
      if (y0k > y1k) cimg::swap(y0k, y1k);
      y0k = (y0k - y_min)/y_step;
      y1k = (y1k - y_min)/y_step;
      y0k = (y0k / std::floor(y0k)) == 1 ? std::floor(y0k) : std::floor(y0k) + 1;
      y1k = std::floor(y1k);
      for (long ky=(long)y0k; ky<=(long) y1k; ky++) {
        if(kx>=0 && ky >= 0 && kx<=x_dim && ky<= y_dim) res(0)(kx,ky) = A + B * (x_min + kx * x_step) + C * (y_min + ky * y_step);
      }
    }
    
    // test the second segment x1 -> x2
    x1k = (x1 - x_min)/x_step;
    double  x2k = (x2 - x_min)/x_step;
    x1k = (x1k / std::floor(x1k)) == 1 ? std::floor(x1k) : std::floor(x1k) + 1;
    x2k = std::floor(x2k);
    for (long kx = (long) x1k; kx <= (long) x2k; kx++) {
      double y1k = a1 * (x_min + kx * x_step) + b1;
      double y2k = a2 * (x_min + kx * x_step) + b2;
      if (y1k > y2k) cimg::swap(y1k, y2k);
      y1k = (y1k - y_min)/y_step;
      y2k = (y2k - y_min)/y_step;
      y1k = (y1k / std::floor(y1k)) == 1 ? std::floor(y1k) : std::floor(y1k) + 1;
      y2k = std::floor(y2k);
      for (long ky=(long)y1k; ky<=(long)y2k; ky++) {
        if(kx>=0 && ky >= 0 && kx<=x_dim && ky<= y_dim) res(0)(kx,ky) = A + B * (x_min + kx * x_step) + C * (y_min + ky * y_step);
      }
    }
  }
  
  *this = res;
  return *this;

}
#endif
#pragma mark - MAX
//! \c get_max_pos: get the position of the maximum value of an image
//T max_pos(long& posx, long& posy);

#pragma mark - MIN
//! \c get_min_pos: get the position of the minimum value of an image
T min_pos(long& posx){
  T min_val = 9.999999e99;
  posx = 0;
  cimg_forX((*this)(0), x) if ((*this)(0)(x) < min_val) {
    min_val = (*this)(0)(x);
    posx = x;
  }
  return min_val;
}

//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
//! \name Geometric / Spatial Manipulation
//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
#pragma mark -
#pragma mark -
#pragma mark GEOMETRIC / SPATIAL MANIPULATION
CImgList<T> get_channel(const int c0) const{
  CImgList<T> res((*this).size(), (*this)(0).width(), (*this)(0).height(), (*this)(0).depth(), (*this)(0).spectrum());
  res(0) = (*this)(0).get_channel(c0);
  if((*this).size() == 2 )
    res(1) = (*this)(1).get_channel(c0);
  return res;
}

//---------------------------------------------------------------------------------------------------------------------
// ---> ROTATE
//---------------------------------------------------------------------------------------------------------------------
//! Transpose the a file
/**
 \param oname Filename, as a C-string corresponding to the transposed file.
 
 <b>Usage:</b>
 \code
 CImgList<float> im;
 im.load_cimg_header("input.cimg");
 im.transpose("output.cimg");
 im.close_cimg_file();
 \endcode
 The \c output.cimg contains now the transposed file.
 \note
 Equivalent to
 \code
 rotate(oname, 4);
 \endcode
 */
void transpose(const char* oname){
  rotate(oname, 4);
}
//! \c rotate (named after IDL routine). Rotate only in multiples of 90 degrees. Work only for file. Use cimg command otherwise
void rotate (const char* oname, int dir){
  dir = dir % 8;
  if (dir == 0) {
    std::cerr << "You selected a direction that do nothing, 1 <= dir <= 7" << std::endl;
    exit(1);
  }
  std::vector<long> dim(2);
  dim[0] = _cimg_xdim; dim[1] = _cimg_ydim;
  // create output object
  CImgList<T> res;
  if (dir == 1 || dir == 3 || dir == 4 || dir == 6) cimg::swap(dim[0], dim[1]);
  res.write_cimg_header(oname, _cimg_cplx, dim[0], dim[1], _cimg_zdim, _cimg_cdim);
  get_tiles(0, TILE_READ);
  res.get_tiles(0, TILE_WRITE);
  res.assign(_cimg_cplx);
  
  // get number of blocks
  //long n_i_blocks = get_bn();
  //long n_o_blocks = o.get_bn();
  
  int mode;
  switch (dir) {
    case 2:  mode = 2; break;
    case 5:  mode = 1; break;
    case 7:  mode = 2; break;
    default: mode = 0; break;
  }
  
  if (mode == 0) {  // rotate
    std::vector<long> idim(4), odim(4);
    idim[0] =    _cimg_xdim;   idim[1] =    _cimg_ydim;   idim[2] =    _cimg_zdim;   idim[3] =    _cimg_cdim;
    odim[0] = res.cimg_xdim(); odim[1] = res.cimg_ydim(); odim[2] = res.cimg_zdim(); odim[3] = res.cimg_cdim();
    int sk = _bn + res.get_bn();
    std::vector< CImgList<T> > tmpIm(_bn);
    std::vector<std::string>   tmpName(_bn);
    for (int k=0; k<_bn; k++) {
      progress(k, sk);
      odim[0] = res.get_bs(k);
      tmpName[k] = temporary_filename();
      //tmpIm.assign(_cimg_cplx, odim[0], odim[1], odim[2], odim[3]);
      (*this).read_cimg(k);
      
      switch (dir) {
        case 1:
          (*this)(0).rotate(90); if (_cimg_cplx == 2) (*this)(1).rotate(90);
          break;
        case 3:
          (*this)(0).rotate(270); if (_cimg_cplx == 2) (*this)(1).rotate(270);
          break;
        case 4:
          (*this)(0).permute_axes("yxzc"); if(_cimg_cplx == 2) (*this)(1).permute_axes("yxzc");
          break;
        case 6:
          (*this)(0).rotate(180).permute_axes("yxzc"); if (_cimg_cplx == 2) (*this)(1).rotate(180).permute_axes("yxzc");
          break;
        default:
          break;
      }
      (*this).save_cimg(tmpName[k].c_str());
      (tmpIm[k]).load_cimg_header((tmpName[k]).c_str());
      (tmpIm[k]).get_tiles();
    }
    for (int k=0; k<res.get_bn(); k++) {
      progress(k+_bn, sk);
      //create output data
      //res.assign(res.cimg_cplx(), res.cimg_xdim(), res.get_bs[k], res.cimg_zdim(), res.cimg_cdim());
      CImgList<T> dummy(_cimg_cplx);
      for (int l=0; l<_bn; l++) {
        int ll;
        ll = dir == 1 || dir == 6 ? _bn - l - 1 : l;
        (tmpIm[ll]).read_cimg(k);
        dummy(0).append((tmpIm[ll])(0), 'x');
        if (_cimg_cplx == 2) dummy(1).append((tmpIm[ll])(1), 'x');
      }
      res(0) = dummy(0); if (_cimg_cplx == 2) res(1) = dummy(1);
      res.write_cimg(k);
    }
    progress();
    for (int k=0; k<_bn; k++) std::remove((tmpName[k]).c_str());
  } else { // flip
    if (mode == 2) res.reverse_tile();
    for (int k=0; k<_bn; k++) {
      progress(k);
      (*this).read_cimg(k);
      switch (dir) {
        case 2:
          (*this)(0).rotate(180); if (_cimg_cplx == 2) (*this)(1).rotate(180);
          break;
        case 5:
          (*this)(0).rotate(270).permute_axes("yxzc"); if (_cimg_cplx == 2) (*this)(1).rotate(270).permute_axes("yxzc");
          break;
        case 7:
          (*this)(0).rotate(90).permute_axes("yxzc"); if (_cimg_cplx == 2) (*this)(1).rotate(90).permute_axes("yxzc");
          break;
        default:
          break;
      }
      res(0) = (*this)(0); if (_cimg_cplx == 2) res(1) = (*this)(1);
      if (mode == 2) res.write_cimg(_bn - 1 - k); else res.write_cimg(k);
    }
    progress();
  }
  res.close_cimg_file();
}


//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
//! \name I/O General Functions
//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
#pragma mark -
#pragma mark -
#pragma mark GENERAL FUNCTIONS

int get_covmat_dim(){
	T delta = 1 + 8 * _cimg_cdim;
	delta = std::sqrt(delta);
	return (int) (delta - 1) / 2;
}

//! Cut out/Crop a rectangular part of an image [work only in a tile processing mode]
#pragma mark - Cut out/Crop
void cut_region(const char* ifname, const char* ofname, int xmin, int xmax, int ymin, int ymax){
  // test region - X
  
  CImgList<T> i;
  i.load_cimg_header(ifname);
  
  
  if (xmin < 0)                    xmin = 0;
  if (xmax < 1)                    xmax = 1;
  if (xmin >= (i.cimg_xdim() - 1)) xmin = i.cimg_xdim() - 2;
  if (xmax >= i.cimg_xdim())       xmax = i.cimg_xdim() - 1;
  if (xmin >= xmax){
    std::cerr << "Check out your \"xmin\"/\"xmax\" input variables" << std::endl;
    exit(1);
  }
  // test region - Y
  if (ymin < 0)                    ymin = 0;
  if (ymax < 1)                    ymax = 1;
  if (ymin >= (i.cimg_ydim() - 1)) ymin = i.cimg_ydim() - 2;
  if (ymax >= i.cimg_ydim())       ymax = i.cimg_ydim() - 1;
  if (ymin >= ymax){
    std::cerr << "Check out your \"ymin\"/\"ymax\" input variables" << std::endl;
    exit(1);
  }
  
  // define output file
  CImgList<T> o;
  o.write_cimg_header(ofname, i.cimg_cplx(), xmax - xmin + 1, ymax - ymin + 1, i.cimg_zdim(), i.cimg_cdim());
  o.get_tiles(0, TILE_WRITE);
  for (int k=0; k<o.get_bn(); k++) {
    o.progress(k);
    int yymin = ymin + o.get_bp(k);
    o.load_cimg(ifname, 0, i.cimg_cplx() - 1, xmin, yymin, 0, 0, xmax, yymin + o.get_bs(k) - 1, i.cimg_zdim(), i.cimg_cdim());
    o.write_cimg(k);
  }
  //o.close_cimg_file();
  o.progress();  
}

//! Transform a complex image to absolute values
#pragma mark - Complex -> Absolute
CImgList<T>& complex2abs(){
  CImgList<T> res(1, (*this)(0).width(), (*this)(0).height(), (*this)(0).depth(), (*this)(0).spectrum(), 0.0);
  res(0) = ((*this)(0).get_sqr() + (*this)(1).get_sqr()).get_sqrt();
  (*this) = res;
  return *this;
}
CImgList<T> get_complex2abs(){return CImgList<T>(*this, false).complex2abs();}
void        get_complex2abs(CImgList<T> **i, CImgList<T> **o, void *p){*o[0] = i[0]->get_complex2abs();}

//! Construct a scattering vector
//#pragma mark - 
#pragma mark - slcs -> vector (construct scattering vector)
//CImgList<T>& construct_vec(int n, ...)

//---------------------------------------------------------------------------------------------------------------------
// ! Pressuming
//#pragma mark - 
#pragma mark - Presumming

#define CSAR_PRESUMMING_DEFAULT 0 
  //!< presumming default mode
#define CSAR_PRESUMMING_PHASE   1 
  //!< presumming in case of a phase data: phase -> complex phase -> presumming -> extract phase from complex

void set_output_dimensions(int ix, int iy, int fx, int fy, int *ox, int *oy, int *rebx){
	*ox   = (int) std::floor((float) ix / (float) fx);
	*oy   = (int) std::floor((float) iy / (float) fy);
	*rebx = *ox * fx - 1;
}

CImgList<T>& presumming(int fx, int fy, int mode = CSAR_PRESUMMING_DEFAULT){
	return get_presumming(fx, fy, mode).move_to(*this);
}

CImgList<T> get_presumming(int fx, int fy, int mode = CSAR_PRESUMMING_DEFAULT){
  //FIXME: Test the function with resize from CIMG
	int S = size();
	int X = (*this)(0).width();
	int Y = (*this)(0).height();
	int Z = (*this)(0).depth();
	int C = (*this)(0).spectrum();
	CImgList<T>res (S, X/fx, Y/fy, Z, C, 0.0);
	cimg_forXYZC(res(0), x, y, z, c){
		for (int kx=0; kx<fx; kx++)
			for (int ky=0; ky<fy; ky++){
				            res(0, x, y, z, c) += (*this)(0)(x*fx+kx, y*fy+ky, z, c);
				if (S == 2) res(1, x, y, z, c) += (*this)(1)(x*fx+kx, y*fy+ky, z, c);
		}
		            res(0, x, y, z, c) /= ((T) fx * (T) fy);
		if (S == 2) res(1, x, y, z, c) /= ((T) fx * (T) fy);
	}
	return res;
}

#ifdef CSAR_PRESUMMING_TILING
void get_presumming(CImgList<T> **i, CImgList<T> **o, void *pP){
	presumming_param_struct *p = static_cast<presumming_param_struct *>(pP);
	*o[0] = i[0]->get_presumming((*p).fx, (*p).fy);
}
#endif



//---------------------------------------------------------------------------------------------------------------------
//! Transform scattering vector into its covariance matrix
//#pragma mark -
#pragma mark - Vector -> Matrix

CImgList<T>& k2m(){
	int S = size();
	int X = (*this)(0).width();
	int Y = (*this)(0).height();
	int Z = (*this)(0).depth();
	int C = (*this)(0).spectrum();
	
	CImgList<T> res(S, X, Y, Z, (C*(C+1))/2, 0.0);
	
	switch (C) {
    case 2:
      res(0).get_shared_channel(0)  = (*this)(0).get_channel(0).get_mul((*this)(0).get_channel(0));
      res(0).get_shared_channel(0) += (*this)(1).get_channel(0).get_mul((*this)(1).get_channel(0));
      res(1).get_shared_channel(0)  = (*this)(1).get_channel(0).get_mul((*this)(0).get_channel(0));
      res(1).get_shared_channel(0) -= (*this)(0).get_channel(0).get_mul((*this)(1).get_channel(0));

      res(0).get_shared_channel(1)  = (*this)(0).get_channel(0).get_mul((*this)(0).get_channel(1));
      res(0).get_shared_channel(1) += (*this)(1).get_channel(0).get_mul((*this)(1).get_channel(1));
      res(1).get_shared_channel(1)  = (*this)(1).get_channel(0).get_mul((*this)(0).get_channel(1));
      res(1).get_shared_channel(1) -= (*this)(0).get_channel(0).get_mul((*this)(1).get_channel(1));
      
      res(0).get_shared_channel(2)  = (*this)(0).get_channel(1).get_mul((*this)(0).get_channel(1));
      res(0).get_shared_channel(2) += (*this)(1).get_channel(1).get_mul((*this)(1).get_channel(1));
      res(1).get_shared_channel(2)  = (*this)(1).get_channel(1).get_mul((*this)(0).get_channel(1));
      res(1).get_shared_channel(2) -= (*this)(0).get_channel(1).get_mul((*this)(1).get_channel(1));
      break;
    case 3:
      res(0).get_shared_channel(0)  = (*this)(0).get_channel(0).get_mul((*this)(0).get_channel(0));
      res(0).get_shared_channel(0) += (*this)(1).get_channel(0).get_mul((*this)(1).get_channel(0));
      res(1).get_shared_channel(0)  = (*this)(1).get_channel(0).get_mul((*this)(0).get_channel(0));
      res(1).get_shared_channel(0) -= (*this)(0).get_channel(0).get_mul((*this)(1).get_channel(0));

      res(0).get_shared_channel(1)  = (*this)(0).get_channel(0).get_mul((*this)(0).get_channel(1));
      res(0).get_shared_channel(1) += (*this)(1).get_channel(0).get_mul((*this)(1).get_channel(1));
      res(1).get_shared_channel(1)  = (*this)(1).get_channel(0).get_mul((*this)(0).get_channel(1));
      res(1).get_shared_channel(1) -= (*this)(0).get_channel(0).get_mul((*this)(1).get_channel(1));

      res(0).get_shared_channel(2)  = (*this)(0).get_channel(0).get_mul((*this)(0).get_channel(2));
      res(0).get_shared_channel(2) += (*this)(1).get_channel(0).get_mul((*this)(1).get_channel(2));
      res(1).get_shared_channel(2)  = (*this)(1).get_channel(0).get_mul((*this)(0).get_channel(2));
      res(1).get_shared_channel(2) -= (*this)(0).get_channel(0).get_mul((*this)(1).get_channel(2));

      res(0).get_shared_channel(3)  = (*this)(0).get_channel(1).get_mul((*this)(0).get_channel(1));
      res(0).get_shared_channel(3) += (*this)(1).get_channel(1).get_mul((*this)(1).get_channel(1));
      res(1).get_shared_channel(3)  = (*this)(1).get_channel(1).get_mul((*this)(0).get_channel(1));
      res(1).get_shared_channel(3) -= (*this)(0).get_channel(1).get_mul((*this)(1).get_channel(1));

      res(0).get_shared_channel(4)  = (*this)(0).get_channel(1).get_mul((*this)(0).get_channel(2));
      res(0).get_shared_channel(4) += (*this)(1).get_channel(1).get_mul((*this)(1).get_channel(2));
      res(1).get_shared_channel(4)  = (*this)(1).get_channel(1).get_mul((*this)(0).get_channel(2));
      res(1).get_shared_channel(4) -= (*this)(0).get_channel(1).get_mul((*this)(1).get_channel(2));

      res(0).get_shared_channel(5)  = (*this)(0).get_channel(2).get_mul((*this)(0).get_channel(2));
      res(0).get_shared_channel(5) += (*this)(1).get_channel(2).get_mul((*this)(1).get_channel(2));
      res(1).get_shared_channel(5)  = (*this)(1).get_channel(2).get_mul((*this)(0).get_channel(2));
      res(1).get_shared_channel(5) -= (*this)(0).get_channel(2).get_mul((*this)(1).get_channel(2));
      break;
    case 4:
      res(0).get_shared_channel(0)  = (*this)(0).get_channel(0).get_mul((*this)(0).get_channel(0));
      res(0).get_shared_channel(0) += (*this)(1).get_channel(0).get_mul((*this)(1).get_channel(0));
      res(1).get_shared_channel(0)  = (*this)(1).get_channel(0).get_mul((*this)(0).get_channel(0));
      res(1).get_shared_channel(0) -= (*this)(0).get_channel(0).get_mul((*this)(1).get_channel(0));
      
      res(0).get_shared_channel(1)  = (*this)(0).get_channel(0).get_mul((*this)(0).get_channel(1));
      res(0).get_shared_channel(1) += (*this)(1).get_channel(0).get_mul((*this)(1).get_channel(1));
      res(1).get_shared_channel(1)  = (*this)(1).get_channel(0).get_mul((*this)(0).get_channel(1));
      res(1).get_shared_channel(1) -= (*this)(0).get_channel(0).get_mul((*this)(1).get_channel(1));
      
      res(0).get_shared_channel(2)  = (*this)(0).get_channel(0).get_mul((*this)(0).get_channel(2));
      res(0).get_shared_channel(2) += (*this)(1).get_channel(0).get_mul((*this)(1).get_channel(2));
      res(1).get_shared_channel(2)  = (*this)(1).get_channel(0).get_mul((*this)(0).get_channel(2));
      res(1).get_shared_channel(2) -= (*this)(0).get_channel(0).get_mul((*this)(1).get_channel(2));
      
      res(0).get_shared_channel(3)  = (*this)(0).get_channel(0).get_mul((*this)(0).get_channel(3));
      res(0).get_shared_channel(3) += (*this)(1).get_channel(0).get_mul((*this)(1).get_channel(3));
      res(1).get_shared_channel(3)  = (*this)(1).get_channel(0).get_mul((*this)(0).get_channel(3));
      res(1).get_shared_channel(3) -= (*this)(0).get_channel(0).get_mul((*this)(1).get_channel(3));
      
      res(0).get_shared_channel(4)  = (*this)(0).get_channel(1).get_mul((*this)(0).get_channel(1));
      res(0).get_shared_channel(4) += (*this)(1).get_channel(1).get_mul((*this)(1).get_channel(1));
      res(1).get_shared_channel(4)  = (*this)(1).get_channel(1).get_mul((*this)(0).get_channel(1));
      res(1).get_shared_channel(4) -= (*this)(0).get_channel(1).get_mul((*this)(1).get_channel(1));
      
      res(0).get_shared_channel(5)  = (*this)(0).get_channel(1).get_mul((*this)(0).get_channel(2));
      res(0).get_shared_channel(5) += (*this)(1).get_channel(1).get_mul((*this)(1).get_channel(2));
      res(1).get_shared_channel(5)  = (*this)(1).get_channel(1).get_mul((*this)(0).get_channel(2));
      res(1).get_shared_channel(5) -= (*this)(0).get_channel(1).get_mul((*this)(1).get_channel(2));
      
      res(0).get_shared_channel(6)  = (*this)(0).get_channel(1).get_mul((*this)(0).get_channel(3));
      res(0).get_shared_channel(6) += (*this)(1).get_channel(1).get_mul((*this)(1).get_channel(3));
      res(1).get_shared_channel(6)  = (*this)(1).get_channel(1).get_mul((*this)(0).get_channel(3));
      res(1).get_shared_channel(6) -= (*this)(0).get_channel(1).get_mul((*this)(1).get_channel(3));
      
      res(0).get_shared_channel(7)  = (*this)(0).get_channel(2).get_mul((*this)(0).get_channel(2));
      res(0).get_shared_channel(7) += (*this)(1).get_channel(2).get_mul((*this)(1).get_channel(2));
      res(1).get_shared_channel(7)  = (*this)(1).get_channel(2).get_mul((*this)(0).get_channel(2));
      res(1).get_shared_channel(7) -= (*this)(0).get_channel(2).get_mul((*this)(1).get_channel(2));
      
      res(0).get_shared_channel(8)  = (*this)(0).get_channel(2).get_mul((*this)(0).get_channel(3));
      res(0).get_shared_channel(8) += (*this)(1).get_channel(2).get_mul((*this)(1).get_channel(3));
      res(1).get_shared_channel(8)  = (*this)(1).get_channel(2).get_mul((*this)(0).get_channel(3));
      res(1).get_shared_channel(8) -= (*this)(0).get_channel(2).get_mul((*this)(1).get_channel(3));
      
      res(0).get_shared_channel(9)  = (*this)(0).get_channel(3).get_mul((*this)(0).get_channel(3));
      res(0).get_shared_channel(9) += (*this)(1).get_channel(3).get_mul((*this)(1).get_channel(3));
      res(1).get_shared_channel(9)  = (*this)(1).get_channel(3).get_mul((*this)(0).get_channel(3));
      res(1).get_shared_channel(9) -= (*this)(0).get_channel(3).get_mul((*this)(1).get_channel(3));
      break;
    default:
      int j = 0;
      for (int c1 = 0; c1 < C; c1++)
        for (int c2 = c1; c2 < C; c2++){
          res(0).get_shared_channel(j)    = (*this)(0).get_channel(c1).get_mul((*this)(0).get_channel(c2));
          res(0).get_shared_channel(j)   += (*this)(1).get_channel(c1).get_mul((*this)(1).get_channel(c2));
          res(1).get_shared_channel(j)    = (*this)(1).get_channel(c1).get_mul((*this)(0).get_channel(c2));
          res(1).get_shared_channel(j++) -= (*this)(0).get_channel(c1).get_mul((*this)(1).get_channel(c2));
        }
      break;
  }
  (*this) = res;
  return *this;
}

CImgList<T> get_k2m(){
  return CImgList<T> (*this, false).k2m();
}

#ifdef CSAR_K2M_TILING
void get_k2m(CImgList<T> **i, CImgList<T> **o, void *pP){
  k2m_param_struct *p = static_cast<k2m_param_struct *> (pP);
  *o[0] = i[0]->get_k2m();
  if ((*p).fx != 1 || (*p).fy != 1) o[0]->presumming((*p).fx, (*p).fy);
}
#endif



/*
//---------------------------------------------------------------------------------------------------------------------
// ! Calculate the span of a vector
//#pragma mark -
#pragma mark - Vector -> Span

CImgList<T>& vec2span(){
	CImgList<T> res;
	res.assign(1, (*this)(0).width(), (*this)(0).height());
	res(0).fill(0.0);
	for(int k = 0; k<(*this)(0).spectrum(); k++){
		res(0) += ((*this)(0).get_channel(k).get_sqr() + (*this)(1).get_channel(k).get_sqr()).get_sqrt();
	}
	res(0) /= (T)(*this)(0).spectrum();
	(*this) = res;
	return *this;
}

CImgList<T> get_vec2span(){
	return CImgList<T>(*this, false).vec2span();
}
// Tiling processing
void get_vec2span(CImgList<T> **i, CImgList<T> **o, void *p){
	*o[0] = i[0]->get_vec2span();
}
*/
//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
//! \name Filters / TRANSFORM
//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
#pragma mark -
#pragma mark -
#pragma mark FILTERS / TRANSFORM 

// Get 2D image integral
CImgList<long double> get_intlImage(){
  return (*this).get_image_2D_integral();
}

CImgList<long double> get_image_2D_integral(){
  int S = (*this).size();
  int X = (*this)(0).width();
  int Y = (*this)(0).height();
  int Z = (*this)(0).depth();
  int C = (*this)(0).spectrum();
  CImgList<long double> s(S, X, Y, Z, C, 0.0);
  CImgList<long double> ii(S, X, Y, Z, C, 0.0);
  for (int n=0; n<S; n++){
    for(int c = 0; c < C; c++)
    for(int z = 0; z < Z; z++)
    for (int y = 0; y < Y; y++)
 		#pragma omp parallel for
    for (int x = 0; x < X; x++){
    	if(y>0) s(n,x,y,z,c) = s(n, x, y-1, z, c) + (*this)(n,x,y,z,c);
      else    s(n,x,y,z,c) = (*this)(n,x,y,z,c);
    }
    for(int c = 0; c < C; c++)
    for(int z = 0; z < Z; z++)
 		#pragma omp parallel for
    for (int y = 0; y < Y; y++)
    for (int x = 0; x < X; x++){
      if(x>0) ii(n,x,y,z,c) = ii(n, x-1, y, z, c) + s(n,x,y,z,c);
      else    ii(n,x,y,z,c) = s(n,x,y,z,c);
    }
  }
  return ii;
}

//-----------------------------------------------------------------------------

#pragma mark - Averaging box car filter

CImgList<T>& speck_mean(int win){
  return (*this).speck_mean(win, win);
}

CImgList<T>& speck_mean(int xwin, int ywin){
  CImg<T> box(xwin, ywin, 1, 1, 1.0);
  for (int i=0; i<size(); i++) {
    (*this)(i).convolve(box);
    (*this)(i) /= (T) (xwin*ywin);
  }
  return *this;
}

CImgList<T> get_speck_mean(int win){
  return (*this).get_speck_mean(win, win);
}

CImgList<T> get_speck_mean(int xwin, int ywin){
  return CImgList<T>(*this, false).speck_mean(xwin, ywin);
}

//-----------------------------------------------------------------------------
CImgList<T>& fast_speck_mean(int xwin, int ywin){
	if (xwin <= 1 && ywin <= 1) return *this; // check if both window sizes have been set
	unsigned long S = (*this).size();
	unsigned long X = (*this)(0).width();
	unsigned long Y = (*this)(0).height();
	unsigned long Z = (*this)(0).depth();
	unsigned long C = (*this)(0).spectrum();
	if (xwin % 2 == 0) xwin += 1;
	if (ywin % 2 == 0) ywin += 1;
	CImgList<T> res(S, X, Y, Z, C, 0.0); // declare output
	CImgList<long double> II(S, X, Y, Z, C, 0.0);  // declade 2D integrale image
	II = (*this).get_image_2D_integral();
	int xH = (xwin - 1) / 2;
	int yH = (ywin - 1) / 2;
	for (int i = 0; i < S; i++){
    for(int c = 0; c < C; c++)
    	for(int z = 0; z < Z; z++)
  			#pragma omp parallel for
    		for(int y = 0; y < Y; y++)
    			for(int x = 0; x < X; x++){
    				int x1 = x - xH - 1;
    				int y1 = y - yH - 1;
    				int x2 = cimg::min(X - 1, x + xH);
    				int y2 = cimg::min(Y - 1, y + yH);
    				if (x1 >= 0 && y1 >= 0) res(i, x, y, z, c) = II(i, x1, y1, z, c) + II(i, x2, y2, z, c) - II(i, x1, y2, z, c) - II(i, x2, y1, z, c);
    				else {
        			if (x1 < 0 && y1 < 0)  res(i, x, y, z, c) = II(i, x2, y2, z, c);
        			if (x1 < 0 && y1 >= 0) res(i, x, y, z, c) = II(i, x2, y2, z, c) - II(i, x2, y1, z, c);
        			if (x1 >= 0 && y1 < 0) res(i, x, y, z, c) = II(i, x2, y2, z, c) - II(i, x1, y2, z, c);
    				}
    }
    res(i) /= (long double)(xwin * ywin);
	}
	(*this) = res;
	return *this;
}

CImgList<T>& fast_speck_mean(int win){return (*this).fast_speck_mean(win, win);}
CImgList<T> get_fast_speck_mean(int xwin, int ywin){return CImgList<T>(*this, false).fast_speck_mean(xwin, ywin);}
CImgList<T> get_fast_speck_mean(int win){           return CImgList<T>(*this, false).fast_speck_mean(win);}

#ifdef CSAR_SPECK_MEAN_TILING
void get_fast_speck_mean(CImgList<T> **i, CImgList<T> **o, void *pP){
	speck_mean_param_struct *p = static_cast<speck_mean_param_struct *>(pP);
  if (i[0]->cimg_zdim() == 1 && i[0]->cimg_cdim() == 1)i[0]->complex2abs();
	*o[0] = i[0]->get_fast_speck_mean((*p).xwin, (*p).ywin);
}
#endif

//! \c fft: compute fast fourier transform: need fftw (use IDL notation)
#ifdef cimg_use_fftw3
CImgList<T>& fft(int dir2, bool debug = false){
  if ((*this)(0).height() != 1 || (*this)(0).depth() != 1 || (*this)(0).spectrum() != 1) {
    std::cerr << "Input array has to be one dimensions: only 'x' is allowed" << std::endl;
    exit(1);
  }
  if ((*this).width() != 2){
    std::cerr << "Input array has to be complex: only 'x' is allowed" << std::endl;
    exit(1);
  }
  fftw_complex *i, *o;
  long k;
  fftw_plan p;
  i = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (*this)(0).width());
  o = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (*this)(0).width());
#pragma omp parallel for
  cimg_forX((*this)(0), x){
    i[x][0] = (*this)(0)(x);
    i[x][1] = (*this)(1)(x);
    if (debug) std::cout << "FFT: this->i " << x << " -> " << i[x][0] << " " << i[x][1] << std::endl;
  }
  p = fftw_plan_dft_1d((*this)(0).width(), i, o, dir2 == -1? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
#pragma omp parallel for
  cimg_forX((*this)(0), x){
    (*this)(0)(x) = dir2 == -1 ? o[x][0] / (T) (*this)(0).width() : o[x][0];
    (*this)(1)(x) = dir2 == -1 ? o[x][1] / (T) (*this)(0).width() : o[x][1];
    if (debug) std::cout << "FFT: i->this " << x << " -> " << i[x][0] << " " << i[x][1] << std::endl;
  }
  fftw_free(i);
  fftw_free(o);
  return *this;
}

CImgList<T>& fft_2D(int dir2, int mode = 0){
  if ((*this)(0).depth() != 1 || (*this)(0).spectrum() != 1) {
    std::cerr << "Input array has to be 1D or 2D elements" << std::endl;
    exit(1);
  }
  if ((*this)(0).height() == 1) return (*this).fft(dir2);
  fftw_complex *i, *o;
  long k, l, dk, dl;
  fftw_plan p;
  int *n, dir3;
  i = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (*this)(0).width() * (*this)(0).height());
  o = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (*this)(0).width() * (*this)(0).height());
  
  if (mode == 2){
    (*this)(0).transpose();
    if ((*this).width() == 2) (*this)(1).transpose();
  }
  long xxx = (*this)(0).width();
  long yyy = (*this)(0).height();
  
#pragma omp parallel for
  for (long kkk = 0; kkk < (xxx * yyy); kkk++){
    i[kkk][0] = (*this)(0, kkk % xxx, kkk / xxx);
    i[kkk][1] = (*this).width() == 2 ? (*this)(1, kkk % xxx, kkk / xxx) : 0.0;
  }
  
  n = (int*) malloc((unsigned) (1) * sizeof(int));
  dir3 = dir2 == -1 ? FFTW_FORWARD : FFTW_BACKWARD;
  dl = (*this)(0).height();
  dk = (*this)(0).width();
  switch (mode) {
    case 0:
      p = fftw_plan_dft_2d(dl, dk, i, o, dir3, FFTW_ESTIMATE);
      break;
    case 1:
      n[0] = dk;
      p = fftw_plan_many_dft(1, n, dl, i, NULL, 1, dk, o, NULL, 1, dk, dir3, FFTW_ESTIMATE);
      dl = 1;
      break;
    case 2:
      n[0] = dk;
      p = fftw_plan_many_dft(1, n, dl, i, NULL, 1, dk, o, NULL, 1, dk, dir3, FFTW_ESTIMATE);
      dl = 1;
      break;
    default:
      std::cerr << "Mode not defined!" << std::endl;
      exit(1);
      break;
  }
  fftw_execute(p);
  fftw_destroy_plan(p);
  free((int *) n);
  if ((*this).width() == 1)(*this).assign(2, xxx, yyy);
  
#pragma omp parallel for
  for (long kkk=0; kkk < (xxx * yyy); kkk++) {
    (*this)(0, kkk % xxx, kkk / xxx) = o[kkk][0]; 
    (*this)(1, kkk % xxx, kkk / xxx) = o[kkk][1];
  }
  if (dir2 == -1){
    (*this)(0) /= (float)(dl * dk);
    (*this)(1) /= (float)(dl * dk);
  }
  if (mode == 2) {
    (*this)(0).transpose();
    (*this)(1).transpose();
  }
  fftw_free(i);
  fftw_free(o);
  return *this;
}

CImgList<T> get_fft_2D(int dir2, int mode = 0){
  return CImgList<T>(*this, false).fft_2D(dir2, mode);
}

#endif
//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
//! \name GEODESY FUNCTIONS
//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
#pragma mark -
#pragma mark -
#pragma mark GEODESY FUNCTIONS

CImgList<T>& ellip2cart_WGS84(){
  double a = 6378137.0;
  double b = 6356752.31425;
  double f = (a - b) / a;
  b = a * (1.0 - f);
  double e_sq = (a*a - b*b) / (a*a);
  double v;
  cimg_forXY((*this)(0), x, y){
    double lon = (double) (*this)(0)(x, y, 0, 0);
    double lat = (double) (*this)(0)(x, y, 0, 1);
    double hgt = (double) (*this)(0)(x, y, 0, 2);
    v = a / std::sqrt(1.0 - e_sq * std::sin(lat) * std::sin(lat));
    (*this)(0)(x, y, 0, 0) = (T) (v + hgt) * std::cos(lat) * std::cos(lon);
    (*this)(0)(x, y, 0, 1) = (T) (v + hgt) * std::cos(lat) * std::sin(lon);
    (*this)(0)(x, y, 0, 2) = (T) (v * (1.0 - e_sq) + hgt) * std::sin(lat);
  }
  return *this;
}

CImgList<T>& cart2ellip_WGS84(){
  double a = 6378137.0;
  double b = 6356752.31425;
  double f = (a - b) / a;
  b = a * (1.0 - f);
  double e_sq = (a*a - b*b) / (a*a);
  double e_prima_sq = (a*a - b*b) / (b*b);
  double v, s, theta, dummy1, dummy2, e0, e1, e2;
  cimg_forXY((*this)(0), x, y){
    double X = (double) (*this)(0)(x, y, 0, 0);
    double Y = (double) (*this)(0)(x, y, 0, 1);
    double Z = (double) (*this)(0)(x, y, 0, 2);
    s = std::sqrt(X*X + Y*Y);
    theta = std::atan2(Z*a, s*b);
    e0     = std::atan2(Y,X);
    dummy1 = std::sin(theta) * std::sin(theta) * std::sin(theta);
    dummy2 = std::cos(theta) * std::cos(theta) * std::cos(theta);
    e1     = std::atan2(Z + e_prima_sq * b * dummy1, s - e_sq * a * dummy2);
    dummy1 = std::sin(e1) * std::sin(e1);
    v      = a / std::sqrt(1 - e_sq * dummy1);
    e2     = (s / std::cos(e1)) - v;
    (*this)(0)(x, y, 0, 0) = (T) e0;
    (*this)(0)(x, y, 0, 1) = (T) e1;
    (*this)(0)(x, y, 0, 2) = (T) e2;
  }
  return *this;
}


//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
//! \name I/O CImg File Format
//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
#pragma mark -
#pragma mark -
#pragma mark IO CIMG FILE FORMAT

//#pragma mark - Return Attributes

//! Return the number of image columns stored in the header file (mainly for tile processing).
/**
Return the image width, i.e. the image dimensions along the X-axis (Horizontal), as 
indicating in the file header. This function is mainly use for tiling processing.
\note
To use cimg_xdim() method read the cimg header file with the \c read_cimg_header method.
Otherwise, the routine will fail.
*/

unsigned int cimg_xdim() const {return _cimg_xdim;}

//! Return the number of image rows stored in  header file (mainly for tile processing).
/**
Return the image height, i.e. the image dimensions along the Y-axis (Vertical), as 
indicating in the file header. This function is mainly use for tiling processing.
\note
To use cimg_ydim() method read the cimg header file with the \c read_cimg_header method.
\code
CImgList<float> im;
im.read_cimg_header("SAR_product.cimg");
std::cout << im.cimg_ydim() << std::endl;
\endcode
*/
unsigned int cimg_ydim() const {return _cimg_ydim;}

//! Return the number of image slices stored in  header file (mainly for tile processing).
/**
Return the image width, i.e. the image dimensions along the Z-axis (Depth), as 
indicating in the file header. This function is mainly use for tiling processing.
\note
To use cimg_zdim() method read the cimg header file with the \c read_cimg_header method.
\code
CImgList<float> im;
im.read_cimg_header("SAR_product.cimg");
std::cout << im.cimg_zdim() << std::endl;
\endcode
*/
unsigned int cimg_zdim() const {return _cimg_zdim;}

//! Return the number of image channels stored in  header file (mainly for tile processing).
/**
Return the image width, i.e. the image dimensions along the C-axis (Spectrum), as 
indicating in the file header. This function is mainly use for tiling processing.
\note
To use cimg_cdim() method read the cimg header file with the \c read_cimg_header method.
\code
CImgList<float> im;
im.read_cimg_header("SAR_product.cimg");
std::cout << im.cimg_cdim() << std::endl;
\endcode
*/
unsigned int cimg_cdim() const {return _cimg_cdim;}


//! Return 1 if the data is real, 2 if complex (corresponds to the number of list).
/**
Return 1 if real data.

Return 2 if complex data.
*/
unsigned int cimg_cplx() const {return _cimg_cplx;}

//! Set a new height for the whole image (mainly for tile processing).
/**
\param y height size.
To be defined
*/
void ydim_cimg(unsigned int y){_cimg_ydim = y;}

unsigned int cimg_type() const {return _cimg_type;}

//---------------------------------------------------------------------------------------------------------------------
// READ A CIMG IMAGE FILE
//---------------------------------------------------------------------------------------------------------------------
//#pragma mark -
#pragma mark - Read CIMG image file

#define CIMG_REAL 1
#define CIMG_COMPLEX 2

//! read cimg header information using filename (for tiling processing)
/**
Read cimg header information.
\param fname Filename, as a C-string.
*/
void load_cimg_header(const char *fname){
  std::FILE *fid;
  fid = cimg::fopen(fname, "rb");
  _cimg_fname = fname;
  _cimg_seek1 = 0;
  _cimg_seek2 = 0;
	char tmp[256] = {0}, str_pixeltype[256] = {0}, strendian[256] = {0};
	unsigned int j, err, N = 0, W, H, D, C, csiz;
	int i;
	// read first line
	j=0;
	while ((i=std::fgetc(fid))!= '\n' && i != EOF && j<256){tmp[j++] = (char)i;tmp[j]   = 0;}
  err = std::sscanf(tmp,"%u%*c%255[A-Za-z_]%*c%255[sA-Za-z_ ]",&N,str_pixeltype,strendian);
  if (err < 2){cimg::fclose(fid);std::cerr << "load_cimg_header() : CImg header not found in file" << std::endl;exit(1);}
  _cimg_cplx = (int) N;
	// GET TYPE
	_cimg_type = 2;
	if (!cimg::strcasecmp(str_pixeltype, "float"))  _cimg_type = 4;
	if (! cimg::strcasecmp(str_pixeltype, "double")) _cimg_type = 5; 
  // READ THE SECOND LINE
  j=0;
  while ((i=std::fgetc(fid))!= '\n' && i != EOF && j<256 ){tmp[j++] = (char)i;tmp[j]   = 0;}
  err = std::sscanf(tmp,"%u %u %u %u #%u",&W,&H,&D,&C,&csiz);
  _cimg_seek1 = std::ftell(fid);
  cimg::fclose(fid);
  //cimg::fclose(_cfid);
  if (err < 4) exit(1);
  
  // assign some dimension information
  _cimg_xdim = (int) W;
  _cimg_ydim = (int) H;
  _cimg_zdim = (int) D;
  _cimg_cdim = (int) C;
  
  // get the position for the first image (the real part)
  if (_cimg_cplx == 1) return;
  
  fid = cimg::fopen(fname, "rb");
  j=0;
  std::fseek(fid, _cimg_seek1 + W * H * D * C * sizeof(T), SEEK_SET);
  while ((i=std::fgetc(fid))!= '\n' && i != EOF && j<256){tmp[j++] = (char)i;tmp[j]   = 0;}
  err = std::sscanf(tmp, "%u %u %u %u #%u", &W, &H, &D, &C, &csiz);
  _cimg_seek2 = std::ftell(fid);
  cimg::fclose(fid);
  if (err < 4) exit(1);
}

//! read cimg header information using file pointer (for tiling processing and for CSAR mode)
/*void  load_cimg_header(std::FILE *fid){
	//_cimg_fid =  fid;
  _cimg_fid = cimg::fopen(_cimg_fname, "rb");
	char tmp[256] = {0}, str_pixeltype[256] = {0}, strendian[256] = {0};
	unsigned int j, err, N = 0, W, H, D, C, csiz;
	int i;
	// read first line
	j=0;
	while ((i=std::fgetc(_cimg_fid))!= '\n' && i != EOF && j<256){tmp[j++] = (char)i;tmp[j]   = 0;}
  err = std::sscanf(tmp,"%u%*c%255[A-Za-z_]%*c%255[sA-Za-z_ ]",&N,str_pixeltype,strendian);
  if (err < 2){close_cimg_file();std::cerr << "load_cimg_header() : CImg header not found in file" << std::endl;exit(1);}
  _cimg_cplx = (int) N;
  // READ THE SECOND LINE
  j=0;
  while ((i=std::fgetc(_cimg_fid))!= '\n' && i != EOF && j<256 ){tmp[j++] = (char)i;tmp[j]   = 0;}
  err = std::sscanf(tmp,"%u %u %u %u #%u",&W,&H,&D,&C,&csiz);
  //cimg::fclose(_cfid);
  if (err < 4) exit(1);

  // assign some dimension information
  _cimg_xdim = (int) W;
  _cimg_ydim = (int) H;
  _cimg_zdim = (int) D;
  _cimg_cdim = (int) C;
  
  // get the position for the first image (the real part)
  _cimg_seek1 = std::ftell(_cimg_fid);
  _cimg_seek2 = 0;
  cimg::fclose(_cimg_fid);
  if (_cimg_cplx == 1) return;
  }
  
  j=0;
  std::fseek(_cimg_fid, _cimg_seek1 + W * H * D * C * sizeof(T), SEEK_SET);
  while ((i=std::fgetc(_cimg_fid))!= '\n' && i != EOF && j<256){tmp[j++] = (char)i;tmp[j]   = 0;}
  err = std::sscanf(tmp, "%u %u %u %u #%u", &W, &H, &D, &C, &csiz);
  _cimg_seek2 = std::ftell(_cimg_fid);
  cimg::fclose(_cimg_fid);
}
*/


//! read a cimg image file in a tiling procedure.
CImgList<T>& read_cimg(const unsigned int k){return (*this).read_cimg(_bp[k], _bp[k] + _bs[k] - 1);}

CImgList<T>& read_cimg(const unsigned int y0, const unsigned int y1){
//  std::ifstream ifile(_cimg_fname, std::ios::out | std::ios::binary);
  unsigned long height = y0 >= 0 && y1 >= 0 && y1 >= y0 ? y1 - y0 + 1 : -1;
	if (height == -1) {std::cerr << "Wrong y0 and y1: y1 > y0 > 0" << std::endl; exit(1);}
	(*this).assign(_cimg_cplx, _cimg_xdim, height, _cimg_zdim, _cimg_cdim);
  (*this).get_cimg(height, y0, 0);
  if (_cimg_cplx == 2) (*this).get_cimg(height, y0, 1);
  return *this;
/*  
  //(*this)(0).fill(0.0);
  //(*this)(1).fill(0.0);
  //return *this;
  unsigned long pos = _cimg_seek1 + _cimg_xdim * y0 * sizeof(T);
  ifile.seekg(pos);
  T *data; 
  data = (T*) std::malloc(sizeof(T) * height * _cimg_xdim);
  ifile.read(reinterpret_cast<char *> (data), sizeof(T) * height * _cimg_xdim);
  for (long k=0; k<(height * _cimg_xdim); k++) (*this)(0, k % _cimg_xdim, k / _cimg_xdim) = data[k];
//  T *data; data = (T*) std::malloc(sizeof(T) * height * _cimg_xdim);
//  ifile.read(reinterpret_cast<char*> (data), height * _cimg_xdim * sizeof(T));
  //ifile.read((char*) &data,  sizeof(T));
  //std::cout << "pos 3: " << pos << " - > " << data[0] << " " << data[1] << " " << sizeof(T) <<  std::endl;
//#pragma omp parallel for
	//cimg_forXY((*this)(0), x, y) (*this)(0, x, y) = data[x + y * _cimg_xdim];
  //std::free(data);
  ifile.close();
  if (_cimg_cplx == 1) return *this;
  ifile.open(_cimg_fname, std::ios::out | std::ios::binary);
  //std::cout << "pos 4: " << pos << std::endl;
  pos = _cimg_seek2 + _cimg_xdim * y0 * sizeof(T);
  //std::cout << "pos 5: " << pos << std::endl;
  ifile.seekg(pos);
  ifile.read(reinterpret_cast<char *> (data), sizeof(T) * height * _cimg_xdim);
  for (long k=0; k<(height * _cimg_xdim); k++) (*this)(1, k % _cimg_xdim, k / _cimg_xdim) = data[k];
  //std::cout << "imag: " << data[11] << " " << data[12] << " " << data[13] << " " << std::endl;
  //std::cout << "pos 6: " << pos << std::endl;
  ////data = (T*) std::malloc(sizeof(T) * height * _cimg_xdim);
  //std::cout << "pos 7: " << pos << std::endl;
//  ifile.read(reinterpret_cast<char*> (data), height * _cimg_xdim * sizeof(T));
  //ifile.read((char*) &data,  sizeof(T));
  //std::cout << "pos 8: " << pos << " - > " << data[10] << " " << data[11] << " " << sizeof(T) <<  std::endl;
//#pragma omp parallel for
  //cimg_forXY((*this)(1), x, y) (*this)(1, x, y) = data[x + y * _cimg_xdim];
  std::free(data);
  ifile.close();
  return *this;*/
}

CImgList<T>& get_cimg(long height, unsigned long y0, int roc){
	for (int z=0; z<_cimg_zdim; z++) {
    for (int c=0; c<_cimg_cdim; c++) {
      
      std::ifstream ifile(_cimg_fname, std::ios::in | std::ios::binary);
      ifile.seekg((roc == 0 ? _cimg_seek1 : _cimg_seek2) + _cimg_xdim * (y0 + _cimg_ydim * (z + c * _cimg_zdim)) * sizeof(T));
      T *data; 
      data = (T*) std::malloc(sizeof(T) * height * _cimg_xdim);
      ifile.read(reinterpret_cast<char *> (data), sizeof(T) * height * _cimg_xdim);
      #pragma omp parallel for
      for (long k=0; k<(height * _cimg_xdim); k++) (*this)(roc, k % _cimg_xdim, k / _cimg_xdim, z, c) = data[k];
      std::free(data);
      ifile.close();
    }
  }
}

//! read a cimg image file from line y0 to line y1
/*CImgList<T>& read_cimg(const unsigned int y0, const unsigned int y1){
	//std::fseek(_cimg_fid, 0, SEEK_SET);
	return (*this).load_cimg(_cimg_fname, 0, _cimg_cplx - 1, 0, y0, 0, 0, _cimg_xdim, y1, _cimg_zdim, _cimg_cdim);
}

/*CImgList<T>& read_cimg(const unsigned int y0, const unsigned int y1){
  unsigned long height = y0 >= 0 && y1 >= 0 && y1 >= y0 ? y1 - y0 + 1 : -1;
	if (height == -1) {std::cerr << "Wrong y0 and y1: y1 > y0 > 0" << std::endl; close_cimg_file(); exit(1);}
  _cimg_fid = cimg::fopen(_cimg_fname, "rb");
  if (_cimg_zdim == 1 && _cimg_cdim == 1) return (*this).cimg2Dimg(height, y0);
  (*this).assign(_cimg_cplx, _cimg_xdim, height, _cimg_zdim, _cimg_cdim, 0.0);
	for (int z=0; z<_cimg_zdim; z++) {
    for (int c=0; c<_cimg_cdim; c++) {
      unsigned long pos = _cimg_seek1 + _cimg_xdim * (y0 + _cimg_ydim * (z + c * _cimg_zdim)) * sizeof(T);
      T *data;
      data = (T *) std::malloc(sizeof(T) * height * _cimg_xdim);
      _cimg_fid = cimg::fopen(_cimg_fname, "rb");
      std::fseek(_cimg_fid, pos, SEEK_SET);
      cimg::fread(data, height * _cimg_xdim, _cimg_fid);
      cimg::fclose(_cimg_fid);
#pragma omp parallel for
      for ( long y = 0; y < height; y++)
      	for( long x = 0; x < _cimg_xdim; x++)
      		(*this)(0, x, y, z, c) = data[x + y * _cimg_xdim];
      std::free(data);
    }
  }
  if (_cimg_cplx == 1) return *this;
  for (int z=0; z<_cimg_zdim; z++) {
    for (int c=0; c<_cimg_cdim; c++) {
      unsigned long pos = _cimg_seek2 + _cimg_xdim * (y0 + _cimg_ydim * (z + c * _cimg_zdim)) * sizeof(T);
      T *data;
      data = (T *) std::malloc(sizeof(T) * height * _cimg_xdim);
      _cimg_fid = cimg::fopen(_cimg_fname, "rb");
      std::fseek(_cimg_fid, pos, SEEK_SET);
      cimg::fread(data, height * _cimg_xdim, _cimg_fid);
      cimg::fclose(_cimg_fid);
#pragma omp parallel for
      for ( long y = 0; y < height; y++)
        for( long x = 0; x < _cimg_xdim; x++)
          (*this)(1, x, y, z, c) = data[x + y * _cimg_xdim];
      std::cout << "data imaginary: " << data[11] << " " << data[12] << " " << data[13] << std::endl;
      std::free(data);
    }
  }
  return *this;
}*/

//---------------------------------------------------------------------------------------------------------------------
// Read a 2D cimg image
/*CImgList<T>& cimg2Dimg(long height, unsigned long y0){
	(*this).assign(_cimg_cplx, _cimg_xdim, height);
  unsigned long pos1 = _cimg_seek1 + _cimg_xdim * y0 * sizeof(T);
  std::fseek(_cimg_fid, pos1, SEEK_SET);
  T *data; data = (T*) std::malloc(sizeof(T) * height * _cimg_xdim);
  cimg::fread(data, height * _cimg_xdim, _cimg_fid);
#pragma omp parallel for
	cimg_forXY((*this)(0), x, y) (*this)(0, x, y) = data[x + y * _cimg_xdim];
  std::free(data);
  if (_cimg_cplx == 1) {
    cimg::fclose(_cimg_fid);
    return *this;
  }
  unsigned long long pos2 = _cimg_seek2 + _cimg_xdim * y0 * sizeof(T);
  std::fseek(_cimg_fid, pos2, SEEK_SET);
  data = (T*) std::malloc(sizeof(T) * height * _cimg_xdim);
  cimg::fread(data, height * _cimg_xdim, _cimg_fid);
#pragma omp parallel for
  cimg_forXY((*this)(1), x, y) (*this)(1, x, y) = data[x + y * _cimg_xdim];
  std::free(data);
  cimg::fclose(_cimg_fid);
  return *this;
}*/


//---------------------------------------------------------------------------------------------------------------------
// WRITE A CIMG FILE HEADER
//---------------------------------------------------------------------------------------------------------------------
//#pragma mark - 
#pragma mark - Write CIMG file header

//! write cimg header
void write_cimg_header(const char *fname, const unsigned int nl, const unsigned int nx, const unsigned int ny = 1, 
                                          const unsigned int nz = 1, const unsigned int nc = 1){
	_cimg_fid = cimg::fopen(fname, "wb");
	write_cimg_header(_cimg_fid, nl, nx, ny, nz, nc);
}

//! write cimg header using file pointer (work only for the CSAR framework)
void write_cimg_header(std::FILE *fid, const unsigned int nl, const unsigned int nx, const unsigned int ny = 1, 
                                       const unsigned int nz = 1, const unsigned int nc = 1){
	_cimg_fid = fid;
	std::fprintf(_cimg_fid, "%u %s\n", nl, pixel_type());
	std::fprintf(_cimg_fid, "%u %u %u %u\n", nx, ny, nz, nc);
	_cimg_seek1 = std::ftell(_cimg_fid);
	_cimg_seek2 = 0;
	if (nl == 2){
		const unsigned long siz = (unsigned long) nx * ny * nz * nc * sizeof(T);
		std::fseek(_cimg_fid, _cimg_seek1 + siz, SEEK_SET);
		std::fprintf(_cimg_fid, "%u %u %u %u\n", nx, ny, nz, nc);
		_cimg_seek2 = std::ftell(_cimg_fid);
	}
	_cimg_xdim = nx; _cimg_ydim = ny; _cimg_zdim = nz; _cimg_cdim = nc;
	_ov = 0;
}

//---------------------------------------------------------------------------------------------------------------------
// WRITE INTO CIMG IMAGE 
//---------------------------------------------------------------------------------------------------------------------
//#pragma mark - 
#pragma mark - Write CIMG file header

//! write cimg file (only in a tiling mode
void write_cimg(const unsigned int k){(*this).write_cimg(_bp[k], _bp[k] + _bs[k] - 1);}

//! write cimg file between ymin and ymax
void write_cimg(const unsigned int y0, const unsigned int y1){
	unsigned long height = y0 >= 0 && y1 >= 0 && y1 >= y0 ? y1 - y0 + 1 : _cimg_ydim;
	unsigned long ymin = 0;
	ymin = y0==0 ? 0 : _ov;
	unsigned long ymax = ymin + height - 1;
	cimg_forCZ((*this)(0), c, z){
		unsigned long pos = _cimg_seek1 + _cimg_xdim * (y0 + _cimg_ydim * (z + c * _cimg_zdim)) * sizeof(T);
		std::fseek(_cimg_fid, pos, SEEK_SET);
		T *data;
		data = (T *) std::malloc(sizeof(T) * height * (*this)(0).width());
		for (unsigned long y = ymin; y <= ymax; y++)
			for(unsigned long x = 0; x < (*this)(0).width(); x++)
				data[x + (y - ymin) * (*this)(0).width()] = (*this)(0, x, y, z, c);
		cimg::fwrite(data, height * (*this)(0).width(), _cimg_fid);
		std::free(data);
	}
	if (size() == 2) {
		cimg_forCZ((*this)(0), c, z){
			unsigned long pos = _cimg_seek2 + _cimg_xdim * (y0 + _cimg_ydim * (z + c * _cimg_zdim)) * sizeof(T);
			std::fseek(_cimg_fid, pos, SEEK_SET);
			T *data;
			data = (T *) std::malloc(sizeof(T) * height * (*this)(0).width());
			for (unsigned long y = ymin; y <= ymax; y++)
				for(unsigned long x = 0; x < (*this)(1).width(); x++)
					data[x + (y - ymin) * (*this)(1).width()] = (*this)(1, x, y, z, c);
			cimg::fwrite(data, height * (*this)(1).width(), _cimg_fid);
			std::free(data);
		}
	}
}

//---------------------------------------------------------------------------------------------------------------------
// CLOSE CIMG FILE
//---------------------------------------------------------------------------------------------------------------------
//#pragma mark - 
#pragma mark - Close CIMG file

//! close cimg file
void close_cimg_file(){cimg::fclose(_cimg_fid);}
void close_cimg_file2(){cimg_ifile.close();}


//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
//!  \name I/O Rat File Format
//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
#pragma mark -
#pragma mark -
#pragma mark IO RAT FILE FORMAT

// Define internal variables
#define RAT_XDR 1
#define RAT_NOXDR 0
#define RAT_MAGICNUMBER 844382546
#define RAT_REAL 1
#define RAT_COMPLEX 2

int rat_xdim() const {return (int) _rat_xdim;}
int rat_ydim() const {return (int) _rat_ydim;}
int rat_zdim() const {return (int) _rat_zdim;}
int rat_vdim() const {return (int) _rat_vdim;}
int rat_type() const {return (int) _rat_type;}
int rat_cplx() const {return (int) _rat_cplx;}
int rat_mode() const {return (int) _rat_mode;}
int rat_ndim() const {return (int) _rat_ndim;}



//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
//#pragma mark - 
#pragma mark - Read RAT header

// load the rat header
void load_rat_header(const char *fname, int imode = 101, bool mode_lex = false, bool mode_3D = false){
	_rat_fid = cimg::fopen(fname, "rb");
	load_rat_header(imode, mode_lex, mode_3D);
}

void load_rat_header(int imode = 101, bool mode_lex = false, bool mode_3D = false){
	_rat_xdim = 1;
	_rat_ydim = 1;
	_rat_zdim = 1;
	_rat_vdim = 1;
	_rat_xdr  = RAT_NOXDR;
  _rat_mode = 0;
	_imode = imode;
  _rat_lex = mode_lex;
  _rat_3D  = mode_3D;
	cimg::fread(&_rat_ndim, 1, _rat_fid);
	// check version
	
  if (_rat_ndim == RAT_MAGICNUMBER){
		std::fseek(_rat_fid, 8, SEEK_SET);
		// read number of dimensions
		cimg::fread(&_rat_ndim, 1, _rat_fid);
		if (_rat_ndim < 1 || _rat_ndim > 4){std::cerr << "Only 1 -> 4 dimensions are supported!" << std::endl; close_rat_file(); exit(1);}
		std::fseek(_rat_fid, 16, SEEK_SET);
		if (_rat_ndim == 4) cimg::fread(&_rat_vdim, 1, _rat_fid);
		if (_rat_ndim >= 3) cimg::fread(&_rat_zdim, 1, _rat_fid);
		cimg::fread(&_rat_xdim, 1, _rat_fid);
		if (_rat_ndim >= 2) cimg::fread(&_rat_ydim, 1, _rat_fid);
		
		// read the idl type of data number
		std::fseek(_rat_fid, 48, SEEK_SET);
		cimg::fread(&_rat_type, 1, _rat_fid);
    _rat_mode = _rat_type; // FIXME: check if there is some rat type code (definitions.pro?)
		// set the header size (in bytes)
		_rat_header = 1000;		
	} else { // ------------------------------ if not version 2, maybe version 1?
		// check if xdr?
		if (_rat_ndim < 0 || _rat_ndim > 9) _rat_xdr = RAT_XDR; if(_rat_xdr) cimg::invert_endianness(_rat_ndim);
		if (_rat_ndim < 0 || _rat_ndim > 9) {std::cerr << "Wrong RAT file format, exit now!" << std::endl; close_rat_file(); exit(1);}
		if (_rat_ndim < 1 || _rat_ndim > 4) {std::cerr << "Only 1 -> 4 dimensions are supported!" << std::endl; close_rat_file(); exit(1);}
		if (_rat_ndim == 4) {cimg::fread(&_rat_vdim, 1, _rat_fid); if (_rat_xdr) cimg::invert_endianness(_rat_vdim);}
		if (_rat_ndim >= 3) {cimg::fread(&_rat_zdim, 1, _rat_fid); if (_rat_xdr) cimg::invert_endianness(_rat_zdim);}
    cimg::fread(&_rat_xdim, 1, _rat_fid); if (_rat_xdr) cimg::invert_endianness(_rat_xdim); 
		if (_rat_ndim >= 2) {cimg::fread(&_rat_ydim, 1, _rat_fid); if (_rat_xdr) cimg::invert_endianness(_rat_ydim);}
    cimg::fread(&_rat_type, 1, _rat_fid); if (_rat_xdr) cimg::invert_endianness(_rat_type);
    cimg::fread(&_rat_mode, 1, _rat_fid); if (_rat_xdr) cimg::invert_endianness(_rat_mode);
    if (_rat_mode == 0) _rat_mode = _rat_type;
		_rat_header = 104 + _rat_ndim * 4; if (_rat_xdr) _rat_header += 4;
	}
  if (_rat_mode == 200 || _rat_mode == 220) _rat_lex = true;
	_rat_cplx = _rat_type == 6 || _rat_type == 9 ? RAT_COMPLEX : RAT_REAL;
}

//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
//#pragma mark - 
#pragma mark - Read RAT image file

// LOAD FULL RAT IMAGE
CImgList<T>& read_rat(){
	return (*this).read_rat(0, _rat_ydim - 1);
}

// LOAD RAT TILE
CImgList<T>& read_rat(int k){
	return (*this).read_rat(_bp[k], _bp[k] + _bs[k] - 1);
}

// LOAD A PORTION OF RAT IMAGE
// only between y0 and y1 (for the moment)
CImgList<T>& read_rat(const unsigned int y0, const unsigned int y1){
	long height = y0 >= 0 && y1 >= 0 && y1 >= y0 ? y1 - y0 + 1 : -1;
	if (height == -1) {std::cerr << "Wrong y0 and y1: y1 > y0 > 0" << std::endl; close_rat_file(); exit(1);}
	unsigned long long idl_ta[16] = {0ULL, 1ULL, 2ULL, 4ULL, 4ULL, 8ULL, 8ULL, 0ULL, 0ULL, 16ULL, 0ULL, 0ULL, 2ULL, 4ULL, 8ULL, 8ULL};
	unsigned long long offset = (unsigned long long) (y0 * _rat_xdim * _rat_zdim * _rat_vdim * idl_ta[_rat_type]);
	std::fseek(_rat_fid, _rat_header + offset, SEEK_SET);
	int var = _rat_type == 6 || _rat_type == 9 ? RAT_COMPLEX : RAT_REAL;
	unsigned long size = (unsigned long) _rat_cplx * height * _rat_xdim * _rat_zdim * _rat_vdim;
	// read data
	T *data;
	data = (T *) std::malloc(sizeof(T) * size);
	cimg::fread(data, size, _rat_fid); if (_rat_xdr) cimg::invert_endianness(data, size);
  if (_rat_ndim == 2) return (*this).rat2img(data, height);
  if (_rat_ndim == 3) {
    if (_rat_3D)  return (*this).rat2tomo(data, height);
    if (_rat_lex) return (*this).rat2lex(data, height);
    return (*this).rat2vec(data, height);
  }
  if (_rat_ndim == 4) {
    if (_rat_lex) return (*this).rat2cov(data, height);
    return (*this).rat2covmat(data, height);
  }

	if (_imode == 101) return (*this).rat2img(data, height);
	if (_imode == 200) return (*this).rat2lex(data, height);
  if (_imode == 210) return (*this).rat2vec(data, height);
  if (_imode == 220) return (*this).rat2cov(data, height);
  if (_imode == 221) return (*this).rat2covmat(data, height);
  //	if (_imode == 209) return (*this).load_rat_vec(data, height);
  //	if (_imode == 210) return (*this).load_rat_vec(data, height);
  //	if (_imode == 220) return (*this).load_rat_cov(data, height);
  //	if (_imode == 221) return (*this).load_rat_covmat(data, height);
  //	if (_imode == 222) return (*this).load_rat_covmat(data, height);
	exit(1);
}

CImgList<T>& rat2tomo(T *data, long height){
  std::cerr << "Not yet implemented! Exit Now!" << std::endl;
  exit(1);
}

//---------------------------------------------------------------------------------------------------------------------
// convert a slc/sla/phase/intensity/... image
CImgList<T>& rat2img(T *data, long height){
	(*this).assign(_rat_cplx, _rat_xdim, height);
#pragma omp parallel for
	cimg_forXY((*this)(0), x, y){
		(*this)(0, x, y) = data[0 + size() * (x + y * _rat_xdim)];
		if (size() == 2) (*this)(1, x, y) = data[1 + size() * (x + y * _rat_xdim)];
	}
	std::free(data);
	return *this;
}

//---------------------------------------------------------------------------------------------------------------------
// load scattering vector, lexicographic basis in rat internal structure
CImgList<T>& rat2lex(T *data, long height){
	(*this).assign(_rat_cplx, _rat_xdim, height, 1, _rat_zdim);
	CImgList<T> im(_rat_cplx, _rat_xdim, height, 1, _rat_zdim, 0.0);
#pragma omg parallel for
	cimg_forCXY((*this)(0), c,x,y){
    im(0, x, y, 0, c) = data[0 + c * size() + x * size() * _rat_zdim + y * size() * _rat_zdim * _rat_xdim];
		if(size() == 2) im(1, x, y, 0, c) = data[1 + c * size() + x * size() * _rat_zdim + y * size() * _rat_zdim * _rat_xdim];
	}
	std::free(data);
	if (_rat_zdim == 3){
		(*this)(0).get_shared_channel(0) = im(0).get_channel(0); if (size() == 2) (*this)(1).get_shared_channel(0) = im(1).get_channel(0); // HH
		(*this)(0).get_shared_channel(1) = im(0).get_channel(2); if (size() == 2) (*this)(1).get_shared_channel(1) = im(1).get_channel(2); // HV
		(*this)(0).get_shared_channel(2) = im(0).get_channel(1); if (size() == 2) (*this)(1).get_shared_channel(2) = im(1).get_channel(1); // VV
		return *this;
	}
	(*this)(0).get_shared_channel(0) = im(0).get_channel(0); if (size() == 2) (*this)(1).get_shared_channel(0) = im(1).get_channel(0); // HH
	(*this)(0).get_shared_channel(1) = im(0).get_channel(2); if (size() == 2) (*this)(1).get_shared_channel(1) = im(1).get_channel(2); // HV
	(*this)(0).get_shared_channel(2) = im(0).get_channel(3); if (size() == 2) (*this)(1).get_shared_channel(2) = im(1).get_channel(3); // VH
	(*this)(0).get_shared_channel(3) = im(0).get_channel(1); if (size() == 2) (*this)(1).get_shared_channel(3) = im(1).get_channel(1); // VV
	return *this;
}

//---------------------------------------------------------------------------------------------------------------------
// load scattering vector, general form
CImgList<T>& rat2vec(T *data, long height){
  (*this).assign(_rat_cplx, _rat_xdim, height, 1, _rat_zdim);
#pragma omp parallel for
  cimg_forCXY((*this)(0), c, x, y){
    (*this)(0, x, y, 0, c) = data[0 + c * size() + x * size() * _rat_zdim + y * size() * _rat_zdim * _rat_xdim];
    if (size() == 2) (*this)(1, x, y, 0, c) = data[1 + c * size() + x * size() * _rat_zdim + y * size() * _rat_zdim * _rat_xdim];
  }
}

//---------------------------------------------------------------------------------------------------------------------
// Load a polarimetric covariance matric [C] (3x3 or 4x4) in the rat internal data structure.
CImgList<T>& rat2cov(T *data, long height){
	int var = _rat_type == 6 || _rat_type == 9 ? RAT_COMPLEX : RAT_REAL;
	(*this).assign(_rat_cplx, _rat_xdim, height, 1, (_rat_zdim * (_rat_zdim+1))/2);
	CImgList<T> im(_rat_cplx, _rat_xdim, height, _rat_zdim, _rat_vdim);
#pragma omp parallel for
	for (int y=0; y<height; y++)
		for (int x=0; x<_rat_xdim; x++)
			for (int z = 0; z<_rat_zdim; z++)
				for (int c = 0; c < _rat_vdim; c++){
					im(0, x, y, z, c) = data[0 + size()*(c + _rat_vdim * (z + _rat_zdim * (x + y * _rat_xdim)))];
					im(1, x, y, z, c) = data[1 + size()*(c + _rat_vdim * (z + _rat_zdim * (x + y * _rat_xdim)))];
        }
	std::free(data);
	if (_rat_zdim == 3){
		(*this)(0).get_shared_channel(0) = im(0).get_channel(0).get_slice(0); if (size() == 2) (*this)(1).get_shared_channel(0) = im(1).get_channel(0).get_slice(0);
		(*this)(0).get_shared_channel(1) = im(0).get_channel(2).get_slice(0); if (size() == 2) (*this)(1).get_shared_channel(1) = im(1).get_channel(2).get_slice(0);
		(*this)(0).get_shared_channel(2) = im(0).get_channel(1).get_slice(0); if (size() == 2) (*this)(1).get_shared_channel(2) = im(1).get_channel(1).get_slice(0);
		(*this)(0).get_shared_channel(3) = im(0).get_channel(2).get_slice(2); if (size() == 2) (*this)(1).get_shared_channel(3) = im(1).get_channel(2).get_slice(2);
		(*this)(0).get_shared_channel(4) = im(0).get_channel(1).get_slice(2); if (size() == 2) (*this)(1).get_shared_channel(4) = im(1).get_channel(1).get_slice(2);
		(*this)(0).get_shared_channel(5) = im(0).get_channel(1).get_slice(1); if (size() == 2) (*this)(1).get_shared_channel(5) = im(1).get_channel(1).get_slice(1);
		return *this;
	}
	std::cerr << "get 4x4 covariance matrix is not yet implemented! Exit now!" << std::endl;
	exit(1);
	return *this;
}


//---------------------------------------------------------------------------------------------------------------------
// load general covariance matrix (not the lexicographic from RAT, usable for coherency matrix)
CImgList<T>& rat2covmat(T *data, long height){
	(*this).assign(_rat_cplx, _rat_xdim, height, 1, (_rat_zdim * (_rat_zdim+1))/2);
#pragma omp parallel for
  for (int y=0; y<height;y++)
		for (int x=0; x<_rat_xdim; x++){
			int cc = -1;
			for (int z=0; z<_rat_zdim; z++)
				for(int c=z; c<_rat_vdim; c++){
          (*this)(0, x, y, 0, ++cc) = data[0 + size()*(c + _rat_vdim * (z + _rat_zdim * (x + y * _rat_xdim)))];
					if(size() == 2) (*this)(1, x, y, 0,   cc) = data[1 + size()*(c + _rat_vdim * (z + _rat_zdim * (x + y * _rat_xdim)))];
				}
		}
	return *this;
}


// load scattering vector
CImgList<T>& load_rat_vec(T *data, long height){
	(*this).assign(_rat_type == 6 || _rat_type == 9 ? RAT_COMPLEX : RAT_REAL, _rat_xdim, height, 1, _rat_zdim);
	cimg_forCXY((*this)(0), c,x,y){
    (*this)(0, x, y, 0, c) = data[0 + c * size() + x * size() * _rat_zdim + y * size() * _rat_zdim * _rat_xdim];
		if(size() == 2) (*this)(1, x, y, 0, c) = data[1 + c * size() + x * size() * _rat_zdim + y * size() * _rat_zdim * _rat_xdim];
	}
	std::free(data);
	return *this;
}


CImgList<T>& load_rat_covmat(T *data, long height){
	(*this).assign(_rat_type == 6 || _rat_type == 9 ? RAT_COMPLEX : RAT_REAL, _rat_xdim, height, 1, (_rat_zdim * (_rat_zdim+1))/2);
	for (int y=0; y<height;y++)
		for (int x=0; x<_rat_xdim; x++){
			int cc = -1;
			for (int z=0; z<_rat_zdim; z++)
				for(int c=z; c<_rat_vdim; c++){
          (*this)(0, x, y, 0, ++cc) = data[0 + size()*(c + _rat_vdim * (z + _rat_zdim * (x + y * _rat_xdim)))];
					if(size() == 2) (*this)(1, x, y, 0,   cc) = data[1 + size()*(c + _rat_vdim * (z + _rat_zdim * (x + y * _rat_xdim)))];
				}
		}
	return *this;
}
/*
 HH HH* HH VV* HH HV*
 VV HH* VV VV* VV HV*
 HV HH* HV VV* HV HV*
 */

// load covariance matrix
CImgList<T>& load_rat_2D_cov(T *data, long height){
	(*this).assign(_rat_type == 6 || _rat_type == 9 ? RAT_COMPLEX : RAT_REAL, _rat_xdim, height, 1, (_rat_zdim * (_rat_zdim+1))/2);
	//CImgList<T> im(_rat_type == 6 || _rat_type == 9 ? RAT_COMPLEX : RAT_REAL, _rat_xdim, height, _rat_zdim, rat_vdim, 0.0);
	for (int y=0; y<height;y++)
		for (int x=0; x<_rat_xdim; x++){
			int cc = -1;
			for (int z=0; z<_rat_zdim; z++)
				for(int c=z; c<_rat_vdim; c++){
          (*this)(0, x, y, 0, ++cc) = data[0 + size()*(c + _rat_vdim * (z + _rat_zdim * (x + y * _rat_xdim)))];
					if(size() == 2) (*this)(1, x, y, 0,   cc) = data[1 + size()*(c + _rat_vdim * (z + _rat_zdim * (x + y * _rat_xdim)))];
				}
		}
	return *this;
	//return *this;
	//cimg_forCZXY(im(0), c, z, x, y) {
	//	                im(0, x, y, z, c) = data[0 + size()*(c + _rat_vdim * (z + _rat_zdim * (x + y * _rat_xdim)))];
	//	if(size() == 2) im(1, x, y, z, c) = data[1 + size()*(c + _rat_vdim * (z + _rat_zdim * (x + y * _rat_xdim)))];
	//}
	//(*this)(0).get_shared_channel(0) = im(0).get_channel(0).get_plane(0);
	//(*this)(1).get_shared_channel(0) = im(1).get_channel(0).get_plane(0);
	//(*this)(0).get_shared_channel(1) = im(0).get_channel(1).get_plane(0);
	//(*this)(1).
	
}


//#pragma mark - 
#pragma mark - Write RAT header


void write_rat_header(const char* fname, std::string mode, 
                      const unsigned int nx, const unsigned int ny = 1,
                      const unsigned int nz = 1, const unsigned int nv = 1){
  _rat_fid = cimg::fopen(fname, "wb");
  write_rat_header(mode, nx, ny, nz, nv);
}

void write_rat_header(std::FILE *rat_fid, std::string mode, 
                      const unsigned int nx, const unsigned int ny = 1,
                      const unsigned int nz = 1, const unsigned int nv = 1){
  _rat_fid = rat_fid;
  write_rat_header(mode, nx, ny, nz, nv);
}


void write_rat_header(std::string mode, 
                      const unsigned int nx, const unsigned int ny = 1,
                      const unsigned int nz = 1, const unsigned int nv = 1){
  _rat_ndim = ny == 1 ? 1 : (nz == 1 && nv == 1 ? 2 : nv == 1 ? 3 : 4);
  if (mode.compare(1, 1, "f") == 0) _rat_type = 4;
  if (mode.compare(1, 1, "c") == 0) _rat_type = 6;
  if (mode.compare(1, 1, "d") == 0) _rat_type = 5;
  if (mode.compare(1, 1, "z") == 0) _rat_type = 9;
  _mode = mode;
	cimg::fwrite(&_rat_ndim, 1, _rat_fid); // write number of dimension
	if(_rat_ndim == 4) cimg::fwrite(&nv, 1, _rat_fid);
	if(_rat_ndim >= 3) cimg::fwrite(&nz, 1, _rat_fid);
	cimg::fwrite(&nx, 1, _rat_fid);
	if(_rat_ndim >= 2) cimg::fwrite(&ny, 1, _rat_fid);
	cimg::fwrite(&_rat_type, 1, _rat_fid);
	_rat_header = 104 + _rat_ndim * 4;
	_rat_cplx = _rat_type == 6 || _rat_type == 9 ? RAT_COMPLEX : RAT_REAL;
	_rat_xdim = nx;
	_rat_ydim = ny;
	_rat_zdim = nz;
	_rat_vdim = nv;	
}


#pragma mark - Write RAT image file

//! write the whole image
void write_rat(){(*this).write_rat(0, _rat_ydim - 1);}

//! write a tile
void write_rat(int k){(*this).write_rat(_bp[k], _bp[k] + _bs[k] - 1);}

// write a portion of RAT image only between y0 and y1
void write_rat(const unsigned int y0, const unsigned int y1){
	long height = y0 >= 0 && y1 >= 0 && y1 >= y0 ? y1 - y0 + 1: -1;
	if (height == -1) {std::cerr << "Wrong y0 and y1: y1 > y0 > 0" << std::endl; close_rat_file(); exit(1);}
	unsigned long long idl_ta[16] = {0ULL, 1ULL, 2ULL, 4ULL, 4ULL, 8ULL, 8ULL, 0ULL, 0ULL, 16ULL, 0ULL, 0ULL, 2ULL, 4ULL, 8ULL, 8ULL};
	unsigned long long offset = (unsigned long long) (y0 * _rat_xdim * _rat_zdim * _rat_vdim * idl_ta[_rat_type]);
	std::fseek(_rat_fid, _rat_header + offset, SEEK_SET);
	unsigned long size = (unsigned long) _rat_cplx * height * _rat_xdim * _rat_zdim * _rat_vdim;
	T *data;
	data = (T *) std::malloc(sizeof(T) * size);
//	if (_imode == 101) (*this).img2rat(data, height);
//  if (_imode == 210) (*this).vec2rat(data, height);
//	if (_imode == 221) (*this).covmat2rat(data, height);
  if (_mode == "-f2Dimg") (*this).cimg_f2Dimg(data, height);
  if (_mode == "-d2Dimg") (*this).cimg_f2Dimg(data, height);
  if (_mode == "-c2Dimg") (*this).cimg_c2Dimg(data, height);
  if (_mode == "-f2Dvec") (*this).cimg_f2Dvec(data, height);
  if (_mode == "-d2Dvec") (*this).cimg_d2Dvec(data, height);
  if (_mode == "-c2Dvec") (*this).cimg_c2Dvec(data, height);
  if (_mode == "-c2Dmat") (*this).cimg_c2Dmat(data, height);
  if (_mode == "-f2Dmat") (*this).cimg_f2Dmat(data, height);
  if (_mode == "-f3Dimg") (*this).cimg_f3Dimg(data, height);
  
	cimg::fwrite(data, size, _rat_fid);
	std::free(data);
	return;
}

void cimg_f2Dimg(T* data, long height){
  cimg_forXY((*this)(0), x, y)
  data[x + y * _rat_xdim] = (*this)(0, x, y);
}
void cimg_c2Dimg(T* data, long height){
#pragma omp parallel for
  cimg_forCXY((*this)(0), c, x, y){
    data[0 + 2 * (c + _rat_zdim * (x + y * _rat_xdim))] = (*this)(0, x, y, 0, c);  
    data[1 + 2 * (c + _rat_zdim * (x + y * _rat_xdim))] = (*this)(1, x, y, 0, c);  
  }
}

void cimg_f2Dvec(T* data, long height){
#pragma omp parallel for
  cimg_forCXY((*this)(0), c, x, y)
  data[c + _rat_zdim * (x + y * _rat_xdim)] = (*this)(0, x, y, 0, c);  
}

void cimg_d2Dvec(T* data, long height){
//#pragma omp parallel for
  cimg_forCXY((*this)(0), c, x, y){
    data[c + _rat_zdim * (x + y * _rat_xdim)] = (*this)(0, x, y, 0, c);  
  }
}

void cimg_c2Dvec(T* data, long height){
#pragma omp parallel for
  cimg_forCXY((*this)(0), c, x, y){
    data[0 + 2 * (c + _rat_zdim * (x + y * _rat_xdim))] = (*this)(0, x, y, 0, c);  
    data[1 + 2 * (c + _rat_zdim * (x + y * _rat_xdim))] = (*this)(1, x, y, 0, c);  
  }
}

void cimg_c2Dmat(T* data, long height){
#pragma omp parallel for
  cimg_forXY((*this)(0), x, y){
    int cc = -1;
    for (int z=0; z<_rat_zdim; z++) {
      for (int v=z; v<_rat_vdim; v++) {
        data[0 + 2 * (v + _rat_vdim * (z + _rat_zdim * (x + y * _rat_xdim)))] = (*this)(0, x, y, 0, ++cc);
        //std::cout << z << " " << v << " " << cc << std::endl;
        if (v != z) data[0 + 2 * (z + _rat_vdim * (v + _rat_zdim * (x + y * _rat_xdim)))] = (*this)(0, x, y, 0,   cc);
        data[1 + 2 * (v + _rat_vdim * (z + _rat_zdim * (x + y * _rat_xdim)))] = (*this)(1, x, y, 0, cc);
        if (v != z) data[1 + 2 * (z + _rat_vdim * (v + _rat_zdim * (x + y * _rat_xdim)))] = - (*this)(1, x, y, 0,   cc);
      }
    }
  }
}

void cimg_f2Dmat(T* data, long height){
#pragma omp parallel for
  cimg_forXY((*this)(0), x, y){
    int cc = -1;
    for (int z=0; z<_rat_zdim; z++) {
      for (int v=0; v<_rat_vdim; v++) {
        data[v + _rat_vdim * (z + _rat_zdim * (x + y * _rat_xdim))] = (*this)(0, x, y, 0, ++cc);
        if (v != z) data[z + _rat_vdim * (v + _rat_zdim * (x + y * _rat_xdim))] = (*this)(0, x, y, 0,   cc);
      }
    }
  }
}


void cimg_f3Dimg(T* data, long height){
#pragma omp parallel for
  cimg_forZXY((*this)(0), z, x, y)
  data[z + _rat_zdim * (x + y * _rat_xdim)] = (*this)(0, x, y, z, 0);
}
/*
// convert cimg -> slc/sla/phase/intensity/...
void img2rat(T *data, long height){
#pragma omp parallel for
	cimg_forXY((*this)(0), x, y){
		data[0 + size() * (x + y * _rat_xdim)] = (*this)(0, x, y);
		if(size() == 2) data[1 + size() * (x + y * _rat_xdim)] = (*this)(1, x, y);
	}
}

// convert cimg -> pauli / general scattering vector
CImgList<T>& vec2rat(T *data, long height){
#pragma omp parallel for
  cimg_forCXY((*this)(0), c, x, y){
    data[0 + c * size() + x * size() * _rat_zdim + y * size() * _rat_zdim * _rat_xdim] = (*this)(0, x, y, 0, c);
    if (size() == 2) data[1 + c * size() + x * size() * _rat_zdim + y * size() * _rat_zdim * _rat_xdim] = (*this)(1, x, y, 0, c);
  }
}

// convert cimg -> coherency / general covariance matrix
void covmat2rat(T *data, long height){
	int S = _rat_cplx;
	int X = _rat_xdim;
	int Y = height;
	int Z = _rat_zdim;
	int V = _rat_vdim;
	CImgList<T> res(S, X, Y, Z, V, 0.0);
#pragma omp parallel for
	for (int y = 0; y < Y; y ++)
		for (int x = 0; x < X; x++){
			int cc = -1;
			for (int z = 0; z < Z; z++)
				for (int v = z; v < V; v++){
          data[0 + S * (v + V * (z + Z * (x + y * X)))] = (*this)(0, x, y, 0, ++cc);
					if (S == 2) data[1 + S * (v + V * (z + Z * (x + y * X)))] = (*this)(1, x, y, 0,   cc);
					if (v != z) {
            data[0 + S * (z + V * (v + Z * (x + y * X)))] = (*this)(0, x, y, 0,   cc);
            if (S == 2) data[1 + S * (z + V * (v + Z * (x + y * X)))] = -(*this)(1, x, y, 0,   cc);
					} 
				}
		}
}

*/
//#pragma mark - 
#pragma mark - Close RAT file
void close_rat_file(){cimg::fclose(_rat_fid);}


//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
//!  \name I/O using GDAL Library
//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
#pragma mark -
#pragma mark -
#pragma mark IO using GDAL Library
#ifdef csar_use_gdal
#define GDAL_REAL    1
#define GDAL_COMPLEX 2
int gdal_xdim() const {return (int) _gdal_xdim;}
int gdal_ydim() const {return (int) _gdal_ydim;}
int gdal_cdim() const {return (int) _gdal_cdim;}
int gdal_type() const {return (int) _gdal_type;}
int gdal_cplx() const {return (int) _gdal_cplx;}

void load_gdal_header(const char *fname){
  GDALAllRegister();
  _gdal_p    = (GDALDataset *) GDALOpen (fname, GA_ReadOnly);
  if (_gdal_p == NULL) exit(1);
  _gdal_xdim = _gdal_p->GetRasterXSize();
  _gdal_ydim = _gdal_p->GetRasterYSize();
  _gdal_cdim = _gdal_p->GetRasterCount();
  _gdal_rb   = _gdal_p->GetRasterBand(1);
  _gdal_type = _gdal_rb->GetRasterDataType();
  _gdal_cplx = _gdal_rb->GetRasterDataType() > 6 && _gdal_rb->GetRasterDataType() <= 11 ? GDAL_COMPLEX : GDAL_REAL;
}

//! Load full image using GDAL library
CImgList<T>& read_gdal(){
  return (*this).read_gdal(0, _gdal_ydim - 1);
}

//! Load tiled image using GDAL library
CImgList<T>& read_gdal(int k){
  return (*this).read_gdal(_bp[k], _bp[k] + _bs[k] - 1);
}

//! Load portion of an image using GDAL library
CImgList<T>& read_gdal(const unsigned y0, const unsigned y1){
	long height = y0 >= 0 && y1 >= 0 && y1 >= y0 ? y1 - y0 + 1 : -1;
	if (height == -1) {std::cerr << "Wrong y0 and y1: y1 > y0 > 0" << std::endl; close_rat_file(); exit(1);}
  (*this).assign(_gdal_cplx, _gdal_xdim, height, 1, _gdal_cdim);
  for (int c=0; c<_gdal_cdim; c++) {
    _gdal_rb = _gdal_p->GetRasterBand(c+1);
    switch (_gdal_rb->GetRasterDataType()) {
      case 2:
        uint16_t *data_2;
        data_2 = (uint16_t *) CPLMalloc(sizeof(uint16_t) * _gdal_xdim * height);
        _gdal_rb->RasterIO(GF_Read, 0, y0, _gdal_xdim, height, data_2, _gdal_xdim, height, _gdal_rb->GetRasterDataType(), 0, 0);
#pragma omp parallel for
        cimg_forXY((*this)(0), x, y) (*this)(0, x, y, 0, c) = (T) data_2[x + y * _gdal_xdim];
        CPLFree(data_2);
        break;    
      case 8:
        short *data_s;
        data_s = (short *) CPLMalloc(sizeof(short) * 2 * _gdal_xdim * height);
        _gdal_rb->RasterIO(GF_Read, 0, y0, _gdal_xdim, height, data_s, _gdal_xdim, height, _gdal_rb->GetRasterDataType(), 0, 0);
//#pragma omp parallel for
        cimg_forXY((*this)(0), x, y){
          (*this)(0, x, y, 0, c) = (T) data_s[2 * (x + y * _gdal_xdim)];
          (*this)(1, x, y, 0, c) = (T) data_s[2 * (x + y * _gdal_xdim) + 1];
        }
        CPLFree(data_s);
        break;
      default:
        T * data;
        data = (T*) CPLMalloc(sizeof(T) * _gdal_cplx * _gdal_xdim * height);
        _gdal_rb->RasterIO(GF_Read, 0, y0, _gdal_xdim, height, data, _gdal_xdim, height, _gdal_rb->GetRasterDataType(), 0, 0);
#pragma omp parallel for
        cimg_forXY((*this)(0), x, y){
          (*this)(0, x, y, 0, c) = (T) data[_gdal_cplx * (x + y * _gdal_xdim)];
          if (_gdal_cplx == 2) (*this)(1, x, y, 0, c) = (T) data[2 * (x + y * _gdal_xdim) + 1];
        }
        CPLFree(data);
        break;
    } 
  }
  return *this;
}

//! Read a band of the image (tiling processing)
CImgList<T>& read_gdal_band(int k, int band){
  return (*this).read_gdal_band(_bp[k], _bp[k] + _bs[k] - 1, band);
}

//! Read a band of a portion of the image using GDAL library
CImgList<T>& read_gdal_band(const unsigned y0, const unsigned y1, int b){
	long height = y0 >= 0 && y1 >= 0 && y1 >= y0 ? y1 - y0 + 1 : -1;
	if (height == -1) {std::cerr << "Wrong y0 and y1: y1 > y0 > 0" << std::endl; close_rat_file(); exit(1);}
  (*this).assign(_gdal_cplx, _gdal_xdim, height, 1, 1);
  _gdal_rb = _gdal_p->GetRasterBand(b);
  switch (_gdal_rb->GetRasterDataType()) {
    case 2:
      uint16_t *data_2;
      data_2 = (uint16_t *) CPLMalloc(sizeof(uint16_t) * _gdal_xdim * height);
      _gdal_rb->RasterIO(GF_Read, 0, y0, _gdal_xdim, height, data_2, _gdal_xdim, height, _gdal_rb->GetRasterDataType(), 0, 0);
#pragma omp parallel for
      cimg_forXY((*this)(0), x, y) (*this)(0, x, y) = (T) data_2[x + y * _gdal_xdim];
      CPLFree(data_2);
      break;    
    case 8:
      short *data_s;
      data_s = (short *) CPLMalloc(sizeof(short) * 2 * _gdal_xdim * height);
      _gdal_rb->RasterIO(GF_Read, 0, y0, _gdal_xdim, height, data_s, _gdal_xdim, height, _gdal_rb->GetRasterDataType(), 0, 0);
#pragma omp parallel for
      cimg_forXY((*this)(0), x, y){
        (*this)(0, x, y) = (T) data_s[2 * (x + y * _gdal_xdim)];
        (*this)(1, x, y) = (T) data_s[2 * (x + y * _gdal_xdim) + 1];
      }
      CPLFree(data_s);
      break;
    default:
      T * data;
      data = (T*) CPLMalloc(sizeof(T) * _gdal_cplx * _gdal_xdim * height);
      _gdal_rb->RasterIO(GF_Read, 0, y0, _gdal_xdim, height, data, _gdal_xdim, height, _gdal_rb->GetRasterDataType(), 0, 0);
#pragma omp parallel for
      cimg_forXY((*this)(0), x, y){
        (*this)(0, x, y) = (T) data[_gdal_cplx * (x + y * _gdal_xdim)];
        if (_gdal_cplx == 2) (*this)(1, x, y) = (T) data[2 * (x + y * _gdal_xdim) + 1];
      }
      CPLFree(data);
      break;
  }
  return *this;
}


//! Close File using GDAL library
void close_gdal_file(){GDALClose( (GDALDatasetH) _gdal_p);}
#endif

//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
//!  \name Tile Processing
//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
#pragma mark -
#pragma mark -
#pragma mark TILE PROCESSING


//! Return the block size for the block \c k
/**
\param k iterator over block loop
*/
int get_bs(int k){return _bs[k];}

//! Return the Y-axis position for the block \c k
unsigned long get_bp(int k){return _bp[k];}

//! Return the number of blocks
int get_bn(){     return _bn;}


#pragma mark - Tile Processing
//! Perform the tile processing.
/**
This method calculate the number of block for each input/output images. It takes into account the fact that data could have different size
but they should have the same ratio.
\param f pointer over the function ...
\param n_idata number of input data
\param n_odata number of output data
\param iList array of CImgList input data
\param oList array of CImgList output data
\param pParam pointer over a structure that could be used for the function \c f
\param ov overlap. Useful for filter, should be half of the size of the filter in Y-Axis
\param bs Block size. if not set use the library defined value: \c BLOCK_SIZE
*/
void tile_processing(void (CImgList<T>::*f)(CImgList<T> **, CImgList<T> **, void*), int n_idata, int n_odata, CImgList<T> **iList, CImgList<T> **oList, void *pParam, int ov = 0, int bs = BLOCK_SIZE){
  int n_data = n_idata + n_odata;
  // initialize data 
  int ydim[n_data], n_ydim[n_data], block_size[n_data];
  for (int k=0; k<n_idata; k++) ydim[k] = iList[k]->cimg_ydim();
  for (int k=0; k<n_odata; k++) ydim[k+n_idata] = oList[k]->cimg_ydim();
  // calculate the ymin
  int ymin = ydim[0];
  for (int k=0; k<n_data; k++) ymin = ymin > ydim[k] ? ydim[k] : ymin;
  // calculate the yratio
  int y_ratio[n_data];
  for (int k=0; k<n_data; k++) y_ratio[k] = std::floor((float)ydim[k] / (float) ymin);
  // calculate new ydim min
  int n_ymin = ydim[0] / y_ratio[0];
  for (int k=0; k<n_data; k++) n_ymin = n_ymin > ydim[k] / y_ratio[k] ? ydim[k] : n_ymin;
  // new dim array
  for (int k=0; k<n_data; k++) n_ydim[k] = n_ymin * y_ratio[k];
  // calculate the new blocksize
  //int bs   = BLOCK_SIZE;
  bool is_done = false;
  while (!is_done) {
    int ymax = -1;
    for (int k=0; k<n_data; k++) {
      block_size[k] = bs * y_ratio[k];
      if (ymax < block_size[k]) ymax = block_size[k];
    }
    if (ymax > MAX_BLOCK_SIZE) bs /= 2; else is_done = true;
  }
  // generate tiles for input/output object
  for (int k=0; k<n_idata; k++){
    iList[k]->ydim_cimg(n_ydim[k]);
    iList[k]->get_tiles(ov, TILE_READ, block_size[k]);
  }
  for (int k=0; k<n_odata; k++){
    oList[k]->ydim_cimg(n_ydim[k+n_idata]);
    oList[k]->get_tiles(ov, TILE_WRITE, block_size[k+n_idata]);
  }

  progress("");
  for (int kl=0; kl<iList[0]->get_bn(); kl++){
  	progress(kl, iList[0]->get_bn());
  	//- read input data
  	for (int k = 0; k<n_idata; k++) iList[k]->read_cimg(kl);
  	//- call the function
		(this->*f)(iList, oList, pParam);
		//- write output data
		for (int k=0; k<n_odata; k++) oList[k]->write_cimg(kl);
  }
  progress();
  std::cout << std::endl;
  return;
}

//! Generate tile in reading and writing mode.
/**
By default generate tiles with no overlap in reading mode.
Normaly, this method is called by \c tile_processing but could be manually called
if needed.
*/
void get_tiles(int ov = 0, int mode = TILE_READ, int bs = 128){
  
  int bls;
  
  _ov = ov;
  _bn = 1;
  
  // calcualte tile
  if (bs >= _cimg_ydim){
    bs  = _cimg_ydim;
    bls = _cimg_ydim;
  } else {
    if (_ov == 0) {
      int pos = bs;
      bls     = bs;
      while (pos < _cimg_ydim) {
        pos += bs;
        _bn++;
        bls = _cimg_ydim - (_bn - 1)* bs;
      }
    } else {
      bool ok = false;
      while (!ok) {
        ok = true;
        int pos1 = 0;
        int pos2 = bs - _ov;
        while (pos2 < _cimg_ydim) {
          pos1 = pos2 - _ov;
          pos2 = pos1 + bs;
          _bn++;
          if (pos2 >= _cimg_ydim) break;
          pos2 -= _ov;
        }
        bls = _cimg_ydim - pos1;
        // if last block too small:
        if (bls < _ov) {_ov++; ok = false;}
      }
    }
  }
  _bs = new int [_bn];
  _bp = new unsigned long  [_bn];
  
  switch (mode) {
    case TILE_READ:
      for (int k=0; k<_bn; k++) { 
      	_bs[k] = bs; 
      	_bp[k] = (unsigned long) k * (bs - 2 * _ov);
      }
      _bs[_bn-1] = bls;
      return;
      break;
    case TILE_WRITE:
      _bs[0] = bs - _ov;
      if (_bn == 1) _bs[0] = bls;
      if (_bn >  2) for (int k=1; k<_bn-1; k++) _bs[k] = bs - 2 * _ov;
      _bp[0] = 0;
      if (_bn > 1) {
        _bs[_bn-1] = bls - _ov;
        for (int k=1; k<_bn; k++) _bp[k] = _bp[k-1] + _bs[k-1];
      }
      return;
      break;
  }
}

void reverse_tile(){
  // reverse block size (just swap first and last one)
  cimg::swap(_bs[0], _bs[_bn-1]);
  // reverse block position (in that case, it is assuming no overlap)
  _bp[0] = 0;
  for (int k=1; k<_bn; k++) _bp[k] = _bp[k-1] + _bs[k-1];
}

//---------------------------------------------------------------------------------------------------------------------
// ---> PROGRESS BAR 
//---------------------------------------------------------------------------------------------------------------------

//#pragma mark - 
#pragma mark - Progress Bar
const char *tile_message;
bool display_progress;
unsigned long tile_prev_percent;

// progress bar in tiling processing
void progress(const char message[256]){
  std::cout << message << "[                                                  ]";
  display_progress = true;
  fflush (stdout);
  tile_message = message;
  tile_prev_percent = 0;
}

void progress (const unsigned long ind){
  double dpercent = ((double) ind + 1.0) * 100.0 / ((double) ((*this).get_bn()));
  int percent = std::floor(dpercent/2.0 + 0.5);
  if (percent == tile_prev_percent) return;
  tile_prev_percent = percent;
  char output[128] = "[                                                  ]";
  for(int k = 1; k <= percent; k++) output[k] = '#';
  for(int k = 0; k <= 51; k++) std::cout <<  "\b";
  std::cout << output;
  fflush(stdout);
}

void progress(const unsigned long ind, const unsigned long nloop){
  double dpercent = ((double) ind + 1.0) * 100.0 / ((double) nloop);
  int percent = std::floor(dpercent/2.0 + 0.5);
  if (percent == tile_prev_percent) return;
  tile_prev_percent = percent;
  char output[128] = "[                                                  ]";
  for(int k = 1; k <= percent; k++) output[k] = '#';
  for(int k = 0; k <= 51; k++) std::cout <<  "\b";
  std::cout << output;
  fflush(stdout);
}

void progress(void){
  std::cout << " - Done! \n" << std::endl;
}

//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
//!  \name File operations
//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------------------
// ---> SET TEMPORARY FILE NAME
//---------------------------------------------------------------------------------------------------------------------
//! Set temporary file name (useful to rotate a file)
const char* temporary_filename(){
  std::string res;
  std::ostringstream oss;
  oss << (long) std::floor(cimg::rand() * 1e10);
  res = "workfile_" + oss.str() + ".cimg";
  return  res.c_str();
}


//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
//!  \name SAR Interferometry
//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
//! Calculate theoretical offset from dem and orbit information for squint 0 geometry. Could be time consumming.
/**
 Calculate the theoretical offset in range and azimuth given an external DEM and orbit data.
 \param demx  Cartesian coordinates of the DEM - X axis
 \param demy  Cartesian coordinates of the DEM - Y axis
 \param demz  Cartesian coordinates of the DEM - Z axis
 \param orb_m Master Orbit
 \param orb_s Slave  Orbit
 \param offset_rg Offset in the <b>range</b> direction
 \param offset_az offset in the azimuth direction
 \param spha Synthetic interferometric phase between master and slave.
 \param ratio Ratio between range and azimuth
 \param sps_rg Range pixel spacing (from the slave image)
 \param k_lp Loop index if inside tile processing (let to zero if not)
 \note
 \code
 Here is the code
 \endcode
*/
void get_offset_theo_squint0(CImgList<double>& demx, CImgList<double>& demy, CImgList<double>& demz,
                             CImgList<double> orb_m, CImgList<double> orb_s,
                             CImgList<double>&offset_rg, CImgList<double>&offset_az, CImgList<double>&spha, long ratio, double sps_rg, long k_lp = 0){
  
  // define local variable
  
  long mlen = (long) orb_m(0).height() / ratio * ratio;
  long slen = (long) orb_s(0).height() / ratio * ratio;
  
  CImgList<double> moX(1, mlen), moY(1, mlen), moZ(1, mlen);
  for (long yy=0; yy<mlen; yy++){
    moX(0, yy) = orb_m(0, 1, yy);
    moY(0, yy) = orb_m(0, 2, yy);
    moZ(0, yy) = orb_m(0, 3, yy);
  }
  moX *= 1e-6; moY *= 1e-6; moZ *= 1e-6;
  CImgList<double>soX(1, slen), soXe(1, 100), soY(1, slen), soYe(1, 100), soZ(1, slen), soZe(1, 100), dumX(1, 100), dumY(1, 100), dumZ(1, 100);
  for (long yy=0; yy<slen; yy++) {
    soX(0, yy) = orb_s(0, 1, yy);
    soY(0, yy) = orb_s(0, 2, yy);
    soZ(0, yy) = orb_s(0, 3, yy);
  }
  soX *= 1e-6; soY *= 1e-6; soZ *= 1e-6;
  cimg_forX(soXe(0), xx){
    soXe(0, xx) = soX(0, 0); dumX(0, xx) = soX(0, slen-1);
    soYe(0, xx) = soY(0, 0); dumY(0, xx) = soY(0, slen-1);
    soZe(0, xx) = soZ(0, 0); dumZ(0, xx) = soZ(0, slen-1);
  }
  soXe(0).append(soX(0)).append(dumX(0));
  soYe(0).append(soY(0)).append(dumY(0));
  soZe(0).append(soZ(0)).append(dumZ(0));
  
  
  std::vector<double> dpos(3);
  CImgList<double> rs0, rsk, dem1d(1, demx(0).width(), 1, 1, 3);
  rs0.assign(1, soXe(0).width());
  rsk.assign(1, 181);
  long az, az_est, az2, jjj;
  double azp;
  
  
  // slave orbit extension
  //long slen = (long) orb_s(0).height() / ratio * ratio;
  //std::cout << "step 08" << std::endl;
  
  //CImgList<double> soXe(1, 100), soYe(1, 100), soZe(1, 100), dumX(1, 100), dumY(1, 100), dumZ(1, 100);
  //std::cout << "step 09" << std::endl;
  
 // cimg_forX(soXe(0), xx) {
 //   soXe(0, xx) = orb_s(0, 1, 0);
 //   soYe(0, xx) = orb_s(0, 2, 0);
 //   soZe(0, xx) = orb_s(0, 3, 0);
 //   dumX(0, xx) = orb_s(0, 1, slen - 1);
 //   dumY(0, xx) = orb_s(0, 2, slen - 1);
 //   dumZ(0, xx) = orb_s(0, 3, slen - 1);
 // }
 // std::cout << "step 10" << std::endl;
  
 // soXe(0).append(orb_s(0).get_row(1)).append(dumX(0));
 // std::cout << "step 11" << std::endl;
 // 
  //soYe(0).append(orb_s(0).get_row(2)).append(dumY(0));
  //std::cout << "step 12" << std::endl;
  
  //soZe(0).append(orb_s(0).get_row(3)).append(dumZ(0));
  //std::cout << "step 13" << std::endl;
  
  
  
  // double loop to estimate offset and generate synthetic phase
  offset_rg.assign(1, demx(0).width(), demx(0).height());
  offset_az.assign(1, demx(0).width(), demx(0).height());
  spha.assign(     1, demx(0).width(), demx(0).height());
  
  for (int ky=0; ky<demx(0).height(); ky++) {
    
   //if (k_lp > 0) jjj = k_lp * demx.get_bs(0) + ky; else jjj = 0;
    jjj = k_lp *demx.get_bs(0) + ky;
    dem1d(0).get_shared_channel(0) = demx(0).get_row(ky);
    dem1d(0).get_shared_channel(1) = demy(0).get_row(ky);
    dem1d(0).get_shared_channel(2) = demz(0).get_row(ky);

    dpos[0] = dem1d(0, 0, 0, 0, 0);
    dpos[1] = dem1d(0, 0, 0, 0, 1);
    dpos[2] = dem1d(0, 0, 0, 0, 2);
    rs0(0) = ((soXe(0) - dpos[0]).get_sqr() + (soYe(0) - dpos[1]).get_sqr() + (soZe(0) - dpos[2]).get_sqr()).get_sqrt();
    rs0.min_pos(az);
   az_est = az > 100 ? (az < (slen+100) ? az : slen + 100 ) : 100;
    for (int kx=0; kx<demx(0).width(); kx++) {
      dpos[0] = dem1d(0, kx, 0, 0, 0);
      dpos[1] = dem1d(0, kx, 0, 0, 1);
      dpos[2] = dem1d(0, kx, 0, 0, 2);
      rsk(0) = ((soXe(0).get_columns(az_est-90, az_est+90) - dpos[0]).get_sqr() + (soYe(0).get_columns(az_est-90, az_est+90) - dpos[1]).get_sqr() + (soZe(0).get_columns(az_est-90, az_est+90) - dpos[2]).get_sqr()).get_sqrt();
      double rslave = rsk.min_pos(az);
      if (az >= 10 && az < 170){
        if (rsk(0,az) == rsk(0, az+2) || rsk(0, az) == rsk(0, az-2))
          offset_az(0, kx, ky) = +9999999.0;
        else {
          azp = 0.5 * (rsk(0)(az-1) - rsk(0)(az+1)) / (rsk(0)(az-1) + rsk(0)(az+1) - 2*rsk(0)(az));
          offset_az(0, kx, ky) = ((double) az + azp - 190.) + (double) az_est - (double) ratio * (double)jjj;
        }
      } else offset_az(0, kx, ky) = +9999999.0;
      az_est += (az - 90);
      az_est = az_est > 100 ? (az_est < (slen+100) ? az_est : slen + 100 ) : 100;
      // FIXME: blablabl
      double m0 = moX(0, ratio * jjj) - dpos[0]; m0 *= m0;
      double m1 = moY(0, ratio * jjj) - dpos[1]; m1 *= m1;
      double m2 = moZ(0, ratio * jjj) - dpos[2]; m2 *= m2;
      double rmaster = std::sqrt(m0 + m1 + m2);
      std::cout << rmaster << " - " << rslave << std::endl;
      offset_rg(0, kx, ky) = (rslave - rmaster) * 1e6 / sps_rg;
      spha(0, kx, ky) = (rslave - rmaster) * 1e6;
    }
  }
}

#endif // csar

