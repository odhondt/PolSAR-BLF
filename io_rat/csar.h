/*
 #
 #  File        : csar.h
 #
 #  Description : A CImg plugin
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
//
//
//
//
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

unsigned long *_bp;
int _ov, *_bs, _bn;

int get_covmat_dim(){
	T delta = 1 + 8 * _cimg_cdim;
	delta = std::sqrt(delta);
	return (int) (delta - 1) / 2;
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

////---------------------------------------------------------------------------------------------------------------------
#endif // csar
//
