/*
 #
 #  File        : SARToolPlugin.h
 #                ( C++ header file - CImg plug-in )
 #
 #  Description : Header file providing basic functionalities for the 
 #                processing of multi-channel SAR images. Please note that
 #                this code is to be used for research purpose and is by no
 #                means a finished product.
 #                This plugin is using the CImg library.
 #                ( http://cimg.sourceforge.net )
 #
 #  Copyright   : Olivier D'Hondt
 #                (https://sites.google.com/site/dhondtolivier/)
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

// IMPORTANT: SARToolPlugin is using eigen library for matrix computations. Please install eigen and put:
// #include <iostream>
// #include <Eigen/Core>
// at the beginning of your code.

/**************************** FILE IO **********************************/

// import for RAT files
CImgList<T>& load_rat(const char* filename)
{
    
  if (filename == NULL){
    std::cerr<<"load_rat: empty filename\n";
    exit(1);
  }

  char *dot = strrchr((char*) filename, '.');
  if (strcmp(dot, ".rat") != 0){
    std::cerr<<"load_rat (%s): extension .rat expected\n";
    exit(1);
  }

  FILE *stream = fopen(filename, "rb");
  std::cout<<"load_rat: opening file "<<filename<<'\n';
  if (stream == 0){
    std::cerr<<"load_rat: could not open %s\n";
    exit(1);
  }

  // Reading header data
  unsigned int Dim;
  unsigned int* SizeVec;
  unsigned int Var;
  unsigned int Type;
  unsigned int Dummy;
  char Info[80];

  if((int)fread(&Dim, sizeof(Dim), 1, stream) != 1){ 
    std::cerr<<"load_rat: problem reading file header.\n";
    std::exit(1);
  }
  cimg::invert_endianness(Dim);
  if(Dim>4){
    std::cerr<<"Cannot (yet) handle files with more than 4 dimensions\n";
    exit(1);
  }
  SizeVec = new unsigned int[Dim];

  std::cout<<"Dim: "<<Dim<<'\n';
  std::cout<<"Size Vec: ";
  for(int i = 0; i< (int) Dim; i++){
    if((int)fread(SizeVec + i, sizeof(SizeVec[i]), 1, stream) != 1){ 
      std::cerr<<"load_rat: problem reading file header.\n";
      exit(1);
    }
    cimg::invert_endianness(SizeVec[i]);
    std::cout<<SizeVec[i]<<" ";
  }
  std::cout<<'\n';

    if(((int)fread(&Var, sizeof(Var), 1, stream) != 1)||
    ((int)fread(&Type, sizeof(Type), 1, stream) != 1)||
    ((int)fread(&Dummy, sizeof(Dummy), 1, stream) != 1)||
    ((int)fread(&Dummy, sizeof(Dummy), 1, stream) != 1)||
    ((int)fread(&Dummy, sizeof(Dummy), 1, stream) != 1)||
    ((int)fread(&Dummy, sizeof(Dummy), 1, stream) != 1)||
    ((int)fread(Info, sizeof(Info), 1, stream) != 1)){
        std::cerr<<"load_rat: problem reading file.\n";
        exit(1);
      }
    cimg::invert_endianness(Var);
    cimg::invert_endianness(Type);
    
    std::cout<<"Var: "<<Var<<'\n';
    std::cout<<"Type: "<<Type<<'\n';
    std::size_t SizeVar; // Size of variables
    int RealOrClx; // 2 if data is complex, 1 if real;
    RealOrClx = 1;
    
    switch(Var){
      case 1: RealOrClx = 1; SizeVar = sizeof(char); break; 
      case 2: RealOrClx = 1; SizeVar = sizeof(int); break;
      case 3: RealOrClx = 1; SizeVar = sizeof(long int); break;
      case 4: RealOrClx = 1; SizeVar = sizeof(float); break;
      case 5: RealOrClx = 1; SizeVar = sizeof(double); break;
      case 6: RealOrClx = 2; SizeVar = sizeof(float); break;
      case 9: RealOrClx = 2; SizeVar = sizeof(double); break;
      default: std::cerr << "load_rat: variable format not recognized.\n";
    }

    switch(Dim){
      case 1: assign(RealOrClx, SizeVec[0], 1, 1, 1, 0.0); break;
      case 2: assign(RealOrClx, SizeVec[0], SizeVec[1], 1, 1, 0.0); break;
      case 3: assign(RealOrClx, SizeVec[1], SizeVec[2], 1, SizeVec[0], 0.0); break;
      case 4: assign(RealOrClx, SizeVec[2], SizeVec[3], 1, SizeVec[0] * SizeVec[1], 0.0); break;
      default: std::cerr << "load_rat: data is limited to 4 dimensions.\n";
    }


// TODO: test the code on other types (int, long...)
// Here, cimg::invert_endianness fails with -O3 compilation
// I replaced it by a simple for loop that swaps bytes

int DimX = (*this)(0).width();
int DimY = (*this)(0).height();
int DimC = (*this)(0).spectrum();

#define _cimg_load_rat(Ts) \
  for(int y = 0; y < DimY; y++) \
    for(int x = 0; x < DimX; x++) \
      for(int c = 0; c < DimC; c++) \
        for(int i = 0; i < RealOrClx; i++){ \
          Ts val; \
          fread(&val, sizeof(val), 1, stream); \
          Ts retVal; \
          char *ToConvert = ( char* ) & val; \
          char *retVec = ( char* ) & retVal; \
          int SSiz = (int) sizeof(Ts); \
          for(int n = 0; n < SSiz; n++){ \
            retVec[n] = ToConvert[SSiz - n - 1]; \
          } \
          (*this)(i, x, y, 0, c) = (T) retVal; \
        } \

    switch(Var){
      case 1: _cimg_load_rat(char); break; 
      case 2: _cimg_load_rat(int); break;
      case 3: _cimg_load_rat(long int); break;
      case 4: _cimg_load_rat(float); break;
      case 5: _cimg_load_rat(double); break;
      case 6: _cimg_load_rat(float); break;
      case 9: _cimg_load_rat(double); break;
      default: std::cerr << "load_rat: variable format not recognized.\n";
    }
  return *this;

}


// Temporary function to write .rat files
// TODO: write generic function
void save_rat_float_hack(const char* filename)
{
    
  FILE *myfile = fopen(filename, "wb");
  std::cout<<"saving file "<<filename<<'\n';
  if (myfile == 0){
    std::cerr<<"load_rat: could not open %s\n";
    exit(1);
  }

  // Writing header data
  unsigned int Dim;
  unsigned int* SizeVec;
  unsigned int Var;
  unsigned int Type;
  unsigned int Dummy=0;
  char Info[80];
  
  Dim = 4;
  cimg::invert_endianness(Dim);
  fwrite(&Dim, sizeof(Dim), 1, myfile); 
  
  SizeVec = new unsigned int[Dim];
    
  SizeVec[0] = 3;
  SizeVec[1] = 3;
  SizeVec[2] = (*this)(0).width();
  SizeVec[3] = (*this)(0).height();
  cimg::invert_endianness(SizeVec[0]);
  cimg::invert_endianness(SizeVec[1]);
  cimg::invert_endianness(SizeVec[2]);
  cimg::invert_endianness(SizeVec[3]);
//  for(int i = 0; i< (int) Dim; i++)
//    fwrite(SizeVec + i, sizeof(SizeVec[i]), 1, myfile); 
    fwrite(SizeVec, sizeof(unsigned int), 1, myfile); 
    fwrite(SizeVec+1, sizeof(unsigned int), 1, myfile); 
    fwrite(SizeVec+2, sizeof(unsigned int), 1, myfile); 
    fwrite(SizeVec+3, sizeof(unsigned int), 1, myfile); 
  
  Var = 6; // Variable of type complex float
  Type = 221; // Coherency matrix (see definitions.pro in RAT)
//  Type = 100;
  cimg::invert_endianness(Var);
  fwrite(&Var, sizeof(Var), 1, myfile);
  cimg::invert_endianness(Type);
  fwrite(&Type, sizeof(Type), 1, myfile);

    fwrite(&Dummy, sizeof(Dummy), 1, myfile);
    fwrite(&Dummy, sizeof(Dummy), 1, myfile);
    fwrite(&Dummy, sizeof(Dummy), 1, myfile);
    fwrite(&Dummy, sizeof(Dummy), 1, myfile);
    fwrite(Info, sizeof(Info), 1, myfile);
    

    int DimX = (*this)(0).width();
    int DimY = (*this)(0).height();
  for(int y = 0; y < DimY; y++)
    for(int x = 0; x < DimX; x++){
      CImgList<float> Mat = (*this).get_tensor_at(x, y);
//      Mat.display();
      for(int s = 0; s < 3; s++){
        for(int t = 0; t < 3; t++){
          float val_re = Mat(0, s, t);
          cimg::invert_endianness(val_re);
          fwrite(&val_re, sizeof(val_re), 1, myfile);
          float val_im = Mat(1, s, t);
          cimg::invert_endianness(val_im);
          fwrite(&val_im, sizeof(val_im), 1, myfile);
        }
      }
    }
  fclose(myfile);

}


/**************************** DISPLAY **********************************/

// Displays the (rescaled) amplitude of complex data
CImgList<T>& dispclxsar(float lambda = 2.5)
{ 
  CImg<T> ImgDisp((*this)(0).width(), (*this)(0).height(), (*this)(0).depth(), 1, 0.0);
  if(size()>2)
    std::cout<<"dispclxsar: only valid for single channel SAR complex images. Ignoring images with list index > 1.";
  for(int i = 0; i < (int) size(); i++)
    ImgDisp += (*this)(i).get_sqr();
  ImgDisp.sqrt();
  T s = lambda*ImgDisp.mean();
  ImgDisp.cut(0, s).display();
  return *this;
}

// Displays the (rescaled) log intensity of complex data
CImgList<T>& dispclxlsar(float lambda = 2, char* MODE = "int")
{ 
  CImg<T> ImgDisp((*this)(0).width(), (*this)(0).height(), (*this)(0).depth(), 1, 0.0);
  if(size()>2)
    std::cout<<"dispclxlsar: only valid for single channel SAR complex images. Ignoring images with list index > 1.";
  for(int i = 0; i < size(); i++){
    ImgDisp += (*this)[i].get_sqr();
  }
  T s;
  if(strcmp(MODE,"amp")==0){
    ImgDisp.sqrt();
    s = lambda*ImgDisp.variance();
  }else{
    s = lambda*std::sqrt(ImgDisp.variance());
  }
  T m = ImgDisp.min();
  if(m<0)
    ImgDisp -= m;

  (ImgDisp+1e-6).log().cut(cimg::min(0, std::log(s)),cimg::max(0, std::log(s))).display();
  return *this;
}

// Displays the (rescaled) amplitude data
CImgList<T>& dispsar(float lambda = 2.5)
{ 
  CImg<T> ImgDisp((*this)(0).width(), (*this)(0).height(), (*this)(0).depth(), 1, 0.0);
  for(int i = 0; i < (int) size(); i++)
    ImgDisp += (*this)(i);
  T m = ImgDisp.min();
  if(m<0)
    ImgDisp -= m;
  ImgDisp.sqrt();
  T s = lambda*ImgDisp.mean();

  ImgDisp.cut(0, s).display();
  return *this;
}


// Displays (rescaled) intensity in logarithmic range
CImgList<T>& displsar(float lambda = 2)
{ 
  CImg<T> ImgDisp((*this)(0).width(), (*this)(0).height(), (*this)(0).depth(), 1, 0.0);
  
  ImgDisp = (*this)(0);
  T s = lambda*std::sqrt(ImgDisp.variance());
  T m = ImgDisp.min();
  if(m<0)
    ImgDisp -= m;

  (ImgDisp+1e-6).log().cut(cimg::min(0, std::log(s)),cimg::max(0, std::log(s))).display();
  return *this;
}

// Displays phase with HSV color scale
// assuming the instance object is an angle image
CImgList<T>& dispclxarg(float MaxDeg = 240, float Saturation = 1.0)
{ 
  CImg<T> ImgArg((*this)(0).width(), (*this)(0).height(), (*this)(0).depth(), 1, 0.0);
  CImg<T> ImgDisp((*this)(0).width(), (*this)(0).height(), (*this)(0).depth(), 3, 1.0);
  
  ImgArg = (*this)(0);
  ImgArg = ((0.5 * ImgArg / M_PI) + 0.5);
  ImgDisp.get_shared_channel(0) = (1-ImgArg) * MaxDeg;
  ImgDisp.get_shared_channel(1).fill(Saturation);
  ImgDisp.HSVtoRGB();
  // Trick: using CImgDisplay to have an unnormalized display
  CImgDisplay disp(cimg_fitscreen(ImgDisp.width(),ImgDisp.height(),ImgDisp.depth()), 0, 0);
  ImgDisp.display(disp,true);
  return *this;
}

// Displays complex coherence between interferometric pairs with HSV color scale
// Color is coding the argument and saturation is coding the magnitude
CImgList<T>& dispclxcoh(float MaxDeg = 240, const char* MODE = "sat")
{ 
  CImg<T> ImgArg((*this)(0).width(), (*this)(0).height(), (*this)(0).depth(), 1, 0.0);
  CImg<T> ImgDisp((*this)(0).width(), (*this)(0).height(), (*this)(0).depth(), 3, 1.0);
  
  ImgArg = ((*this).get_clxarg())(0);
  ImgArg = ((0.5 * ImgArg / M_PI) + 0.5);
  ImgDisp.get_shared_channel(0) = (1-ImgArg) * MaxDeg;
  int dispmode = 1;
  if(strcmp(MODE, "val")==0)
    dispmode = 2;
  ImgDisp.get_shared_channel(dispmode) = ((*this).get_clxmag())(0).get_cut(0.0, 1.0);
  ImgDisp.HSVtoRGB();
  // Trick: using CImgDisplay to have an unnormalized display
  CImgDisplay disp(cimg_fitscreen(ImgDisp.width(),ImgDisp.height(),ImgDisp.depth()), 0, 0);
  ImgDisp.display(disp,true);
  return *this;
}

// Display without normalization
CImgList<T>& dispnn()
{
  CImgDisplay disp(cimg_fitscreen((*this)(0).width(),(*this)(0).height(),(*this)(0).depth()), 0, 0);
  (*this).display(disp,true);
  return *this;
}


// Displays a color image of complex covariance image
// IND_NORM imposes independent normalization of the channels like in RAT!!!!
void dispcolcov3(float lambda = 2.0, bool lex = false, bool IND_NORM = true)
{ 
  CImgList<T> ImgDisp(1, (*this)(0).width(), (*this)(0).height(), (*this)(0).depth(), 3, 0.0);
  if(size()!=2)
    std::cout<<"dispcov3: only valid for complex images of 3x3 covariance matrices.";
  else{
    if(!lex){
      // |S_hh + S_vv| -> B
      ImgDisp(0).get_shared_channel(2) = (*this)(0).get_channel(0);
      // |S_hh - S_vv| -> R
      ImgDisp(0).get_shared_channel(0) = (*this)(0).get_channel(3);
      // |S_hv| -> G
      ImgDisp(0).get_shared_channel(1) = (*this)(0).get_channel(5);
    }else{
      // |S_vv| -> R
      ImgDisp(0).get_shared_channel(0) = (*this)(0).get_channel(5);
      // |S_hv| -> G
      ImgDisp(0).get_shared_channel(1) = (*this)(0).get_channel(3);
      // |S_hh| -> B
      ImgDisp(0).get_shared_channel(2) = (*this)(0).get_channel(0);
    }
    ImgDisp(0)-=ImgDisp(0).min();
    ImgDisp(0).sqrt();
    T s = lambda*ImgDisp(0).mean();
    if(IND_NORM){
      s = lambda*ImgDisp(0).get_channel(0).mean();
      ImgDisp(0).get_shared_channel(0).cut(0, s);
      ImgDisp(0).get_shared_channel(0) = ImgDisp(0).get_channel(0) * 255.0 / s;
      //      ImgDisp(0).get_shared_channel(0).normalize(0, 255);
      s = lambda*ImgDisp(0).get_channel(1).mean();
      ImgDisp(0).get_shared_channel(1).cut(0, s);
      ImgDisp(0).get_shared_channel(1) = ImgDisp(0).get_channel(1) * 255.0 / s;
      //      ImgDisp(0).get_shared_channel(1).normalize(0, 255);
      s = lambda*ImgDisp(0).get_channel(2).mean();
      ImgDisp(0).get_shared_channel(2).cut(0, s);
      ImgDisp(0).get_shared_channel(2) = ImgDisp(0).get_channel(2) * 255.0 / s;
      //      ImgDisp(0).get_shared_channel(2).normalize(0, 255);
    }
    ImgDisp(0).display();
  }
}

// Returns a color image of complex covariance image
// IND_NORM imposes independent normalization of the channels like in RAT!!!!
CImgList<T> get_colcov3(float lambda = 2.5, bool lex = false, bool IND_NORM = true)
{ 
  CImgList<T> ImgDisp(1, (*this)(0).width(), (*this)(0).height(), (*this)(0).depth(), 3, 0.0);
  if(size()!=2)
    std::cout<<"dispcov3: only valid for complex images of 3x3 covariance matrices.";
  else{
    if(!lex){
      // |S_hh + S_vv| -> B
      ImgDisp(0).get_shared_channel(2) = (*this)(0).get_channel(0);
      // |S_hh - S_vv| -> R
      ImgDisp(0).get_shared_channel(0) = (*this)(0).get_channel(3);
      // |S_hv| -> G
      ImgDisp(0).get_shared_channel(1) = (*this)(0).get_channel(5);
    }else{
      // |S_vv| -> R
      ImgDisp(0).get_shared_channel(0) = (*this)(0).get_channel(5);
      // |S_hv| -> G
      ImgDisp(0).get_shared_channel(1) = (*this)(0).get_channel(3);
      // |S_hh| -> B
      ImgDisp(0).get_shared_channel(2) = (*this)(0).get_channel(0);
    }
    ImgDisp(0)-=ImgDisp(0).min();
    ImgDisp(0).sqrt();
//    ImgDisp(0) = ImgDisp(0).get_mul(ImgDisp(0).get_sqrt());
    T s = lambda*ImgDisp(0).mean();
    if(IND_NORM){
      s = lambda*ImgDisp(0).get_channel(0).mean();
      ImgDisp(0).get_shared_channel(0).cut(0, s);
      ImgDisp(0).get_shared_channel(0) = ImgDisp(0).get_channel(0) * 255.0 / s;
//      ImgDisp(0).get_shared_channel(0).normalize(0, 255);
      s = lambda*ImgDisp(0).get_channel(1).mean();
      ImgDisp(0).get_shared_channel(1).cut(0, s);
      ImgDisp(0).get_shared_channel(1) = ImgDisp(0).get_channel(1) * 255.0 / s;
//      ImgDisp(0).get_shared_channel(1).normalize(0, 255);
      s = lambda*ImgDisp(0).get_channel(2).mean();
      ImgDisp(0).get_shared_channel(2).cut(0, s);
      ImgDisp(0).get_shared_channel(2) = ImgDisp(0).get_channel(2) * 255.0 / s;
//      ImgDisp(0).get_shared_channel(2).normalize(0, 255);
    }
  }
  return ImgDisp;
}

CImgList<T> get_coljet(bool NORM = false) 
{
  if((*this).size()>1){
    std::cerr<<"get_coljet(): only valid for images of real numbers (list size should be = 1).\n";
  }
  CImg<> CMap;
  CMap = CImg<>::jet_LUT256();
  //CMap.display();
  CImgList<T> Disp(1, (*this)(0).width(),(*this)(0).height(),(*this)(0).depth(),3,0.0);
  CImg<float> ImNorm((*this)(0), "xyz", 0);
  if((*this)(0).spectrum() == 1){
    ImNorm = (*this)(0);
  }else{
    cimg_forC((*this)(0), c){
      ImNorm += (*this)(0).get_channel(c);
    }
    ImNorm /= (*this)(0).spectrum();
  }
  if(NORM)
    ImNorm.normalize(0.0,255.0); 
  cimg_forXYZ(ImNorm, x, y, z){
    int Col = ImNorm(x, y, z);
    int Idx;
    Idx = (Col>255?255:Col);
    Idx = (Col<0?0:Col);
    Disp(0, x, y, z, 0) = CMap(0, Idx, 0, 0);
    Disp(0, x, y, z, 1) = CMap(0, Idx, 0, 1);
    Disp(0, x, y, z, 2) = CMap(0, Idx, 0, 2);
  }
  return Disp;
}

CImgList<T> get_sarcut(float lambda = 2.0) const {
  if((*this).size()>1) {
    std::cerr<<"get_cutsar(): only valid for images of real numbers (list size should be = 1).\n";
  }
  CImgList<T> Disp(*this);

  float avg = Disp(0).mean();
  Disp(0).cut(0, lambda*avg);

  return Disp;
}


/**************************** COMPUTATIONS **********************************/


//CImgList<T>& cropXY(int x0, int y0, int x1, int y1)
//{
//  CImgList<T> Res(*this);
//  cimg_for
//}
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
CImgList<T> operator-(const CImgList<t>& img) const {
  return CImgList<T>(*this,false)-=img;
}

template<typename t>
CImgList<T>& operator*=(const t val) {
  cimglist_for(*this,l){
    (*this)(l) *= val;
  }
  return *this;
}


template<typename t>
CImgList<T> operator*(const t val) const {
  return CImgList<T>(*this,false)*=val;
}


template<typename t>
CImgList<T>& operator/=(const t val) {
  cimglist_for(*this,l){
    (*this)(l) /= val;
  }
  return *this;
}


template<typename t>
CImgList<T> operator/(const t val) const {
  return CImgList<T>(*this,false)/=val;
}



// Computes complex squared magnitude
CImgList<T> get_clxmag() const{
  if((*this).size() != 2)
    std::cerr<<"clxmag: Error, should be applied to CImgList with dimension 2 representing a complex type.\n";
  CImgList<T> res(1, (*this)(0).width(), (*this)(0).height(), (*this)(0).depth(), 1, 0.0);
  res(0) = (*this)(0).get_sqr() + (*this)(1).get_sqr(); 
  return res;
} 

// Computes complex phase (argument)
CImgList<T> get_clxarg() const{
  if((*this).size() != 2)
    std::cerr<<"clxarg: Error, should be applied to CImgList with dimension 2 representing a complex type.\n";
  CImgList<T> res(1, (*this)(0).width(), (*this)(0).height(), (*this)(0).depth(), 1, 0.0);
  res(0) = (*this)(1).get_atan2((*this)(0)); 
  return res;
} 

// Computes complex conjgate
CImgList<T>& clxcj(){
  if((*this).size() != 2)
    std::cerr<<"clxcj: Error, should be applied to CImgList with dimension 2 representing a complex type.\n";
  cimg_for((*this)(1),ptrd,T) { const T val = *ptrd; *ptrd = -val; };
  return *this;
}

CImgList<T> get_clxcj(){
  return CImgList<T>(*this, false).clxcj();
}


// pointwise product of two complex images
CImgList<T>& clxmul(const CImgList<T>& iml){
  if((*this).size() != 2 || iml.size() != 2)
    std::cerr<<"clxmul: Error, should be applied to CImgList with dimension 2 representing a complex type.\n";
   CImgList<T> res(iml);
   res(0) = (*this)(0).get_mul(iml(0)) - (*this)(1).get_mul(iml(1));
   res(1) = (*this)(1).get_mul(iml(0)) + (*this)(0).get_mul(iml(1));
   *this=res;
   return *this;
}

CImgList<T> get_clxmul(const CImgList<T>& iml){
  return CImgList<T>(*this, false).clxmul(iml);
}

// pointwise division of two complex images
CImgList<T>& clxdiv(const CImgList<T>& iml){
  if((*this).size() != 2 || iml.size() != 2)
    std::cerr<<"clxdiv: Error, should be applied to CImgList with dimension 2 representing a complex type.\n";
   CImgList<T> res(iml);
   CImgList<T> mag = iml.get_clxmag();
   res(0) = (*this)(0).get_mul(iml(0)) + (*this)(1).get_mul(iml(1)); 
   res(1) = (*this)(1).get_mul(iml(0)) - (*this)(0).get_mul(iml(1));
   res(0).div(mag(0)); 
   res(1).div(mag(1)); 
   *this=res;
   return *this;
}

CImgList<T> get_clxdiv(const CImgList<T>& iml){
  return CImgList<T>(*this, false).clxdiv(iml);
}

// pointwise product of a complex image by the conjugate of the other
CImgList<T>& clxmulcj(const CImgList<T>& iml){
  if((*this).size() != 2 || iml.size() != 2)
    std::cerr<<"clxmul: Error, should be applied to CImgList with dimension 2 representing a complex type.\n";
   CImgList<T> res(iml);
   res(0) = (*this)(0).get_mul(iml(0)) + (*this)(1).get_mul(iml(1));
   res(1) = (*this)(1).get_mul(iml(0)) - (*this)(0).get_mul(iml(1));
   *this = res;
   return *this;
}

CImgList<T> get_clxmulcj(const CImgList<T>& iml){
  return CImgList<T>(*this, false).clxmulcj(iml);
}

// pointwise product of a complex image by a real one
CImgList<T>& rlmul(const CImgList<T>& iml){
  if((*this).size() > 2 || iml.size() != 1)
    std::cerr<<"clxmul: Error, should be applied to multiply a complex by a real number.\n";
  (*this)(0).mul(iml(0));
  (*this)(1).mul(iml(0));
  return *this;
}

CImgList<T> get_rlmul(const CImgList<T>& iml){
  return CImgList<T>(*this, false).rlmul(iml);
}

// pointwise division of a complex image by a real one
CImgList<T>& rldiv(const CImgList<T>& iml){
  if((*this).size() > 2 || iml.size() != 1)
    std::cerr<<"clxmul: Error, should be applied to divide a complex by a real number.\n";
  (*this)(0).div(iml(0));
  (*this)(1).div(iml(0));
  return *this;
}

CImgList<T> get_rldiv(const CImgList<T>& iml){
  return CImgList<T>(*this, false).rldiv(iml);
}


// Computes the 2D integral image for fast boxcar filtering
CImgList<long double> get_intlImage()
{
  int S = (*this).size();
  int W = (*this)(0).width();
  int H = (*this)(0).height();
  int C = (*this)(0).spectrum();
  CImgList<long double> s(S, W, H, 1, C, 0.0);
  CImgList<long double> ii(S, W, H, 1, C, 0.0);
  for(int n = 0; n < (int) size(); n++){
    cimg_forXYC(ii(n), x, y, v){
      if(y > 0)
        s(n ,x, y, 0, v) = s(n, x, y - 1, 0, v) + (*this)(n, x, y, 0, v);
      else
        s(n, x, y, 0, v) = (*this)(n, x, y, 0, v);
    }
    cimg_forXYC(ii(n), x, y, v){
      if(x > 0)
        ii(n, x, y, 0, v) = ii(n, x - 1, y, 0, v) + s(n, x, y, 0, v);
      else
        ii(n, x, y, 0, v) = s(n, x, y, 0, v);
    }
  }
  return ii;
}

// Integral image of selected channel with index idx
CImgList<long double> get_intlImage(int idx)
{
  int S = (*this).size();
  int W = (*this)(0).width();
  int H = (*this)(0).height();
  CImgList<long double> s(S, W, H, 1, 1, 0.0);
  CImgList<long double> ii(S, W, H, 1, 1, 0.0);
  for(int n = 0; n < (int) size(); n++){
    cimg_forXY(ii(n), x, y){
      if(y > 0)
        s(n ,x, y) = s(n, x, y - 1) + (*this)(n, x, y, 0, idx);
      else
        s(n, x, y) = (*this)(n, x, y, 0, idx);
    }
    cimg_forXY(ii(n), x, y){
      if(x > 0)
        ii(n, x, y) = ii(n, x - 1, y) + s(n, x, y);
      else
        ii(n, x, y) = s(n, x, y);
    }
  }
  return ii;
}



//CImgList<T>& Test()
//{
//using Eigen::MatrixXd;
//  MatrixXd m(2,2);
//
//  m(0,0) = 3;
//  m(1,0) = 2.5;
//  m(0,1) = -1;
//  m(1,1) = m(1,0) + m(0,1);
//  std::cout << m << std::endl;
//
//  return *this;
//}

// Simple boxcar filter
CImgList<T>& boxcar(int BoxSize)
{
  CImg<T> Box(BoxSize, BoxSize, 1, 1, 1.0);
  for(int i = 0; i < size(); i++){
    (*this)(i).convolve(Box);
    (*this)(i) /= (T) (BoxSize*BoxSize);
  }

  return *this;
}

CImgList<T> get_boxcar(int BoxSize)
{
  return CImgList<T>(*this, false).boxcar(BoxSize);
}

// Fast boxcar filter using 2D integral image
CImgList<T>& fastboxcar(int BoxSize) 
{

  int Siz = (*this).size();
  int DimX = (*this)(0).width();
  int DimY = (*this)(0).height();
  int DimC = (*this)(0).spectrum();

  if(BoxSize % 2 == 0){ 
    std::cout<<"SARToolPlugin fastboxcar: Window size has to be odd, adding 1 to the original size\n";
    BoxSize += 1;
  } 

  CImgList<long double> res(Siz, DimX, DimY, 1, DimC, 0.0);
  if(BoxSize > 1){
    CImgList<long double> II(Siz, DimX, DimY, 1, DimC, 0.0);
    II = (*this).get_intlImage();

    int H = (BoxSize - 1) / 2;

    for(int i = 0; i < Siz; i++){
      cimg_forXYC((*this)(i), x, y, v){
        int x1 = x - H - 1;  
        int y1 = y - H - 1;  
        int x2 = cimg::min(DimX - 1, x + H);  
        int y2 = cimg::min(DimY - 1, y + H);

        if(x1 >= 0 && y1 >= 0){
          res(i, x, y, 0, v) = II(i, x1, y1, 0, v) + II(i, x2, y2, 0, v) - II(i, x1, y2, 0, v) - II(i, x2, y1, 0, v);
        }
        else{
          if(x1 < 0 && y1 < 0){
            res(i, x, y, 0, v) = II(i, x2, y2, 0, v);
          }
          if(x1 < 0 && y1 >= 0){
            res(i, x, y, 0, v) = II(i, x2, y2, 0, v) - II(i, x2, y1, 0, v);
          }
          if(x1 >= 0 && y1 < 0){
            res(i, x, y, 0, v) = II(i, x2, y2, 0, v) - II(i, x1, y2, 0, v);
          }
        }
      }
      res(i) /= (long double) (BoxSize*BoxSize);
    }

//  (*this) = res;
  }
//  return res;
  return res.move_to(*this);
}

CImgList<T> get_fastboxcar(int BoxSize) const
{
  return CImgList<T>(*this, false).fastboxcar(BoxSize);
}

// Fast boxcar filter using 2D integral image
// overloaded functions to filter only one channel and return only 
// this image
CImgList<T>& fastboxcar(int BoxSize, int idx)
{

  int Siz = (*this).size();
  int DimX = (*this)(0).width();
  int DimY = (*this)(0).height();

  if(BoxSize % 2 == 0){ 
    std::cout<<"SARToolPlugin fastboxcar: Window size has to be odd, adding 1 to the original size\n";
    BoxSize += 1;
  } 

  CImgList<long double> res(Siz, DimX, DimY, 1, 1, 0.0);
  if(BoxSize > 1){
    CImgList<long double> II(Siz, DimX, DimY, 1, 1, 0.0);
    II = (*this).get_intlImage(idx);

    int H = (BoxSize - 1) / 2;

    for(int i = 0; i < Siz; i++){
      cimg_forXY((*this)(i), x, y){
        int x1 = x - H - 1;  
        int y1 = y - H - 1;  
        int x2 = cimg::min(DimX - 1, x + H);  
        int y2 = cimg::min(DimY - 1, y + H);

        if(x1 >= 0 && y1 >= 0){
          res(i, x, y) = II(i, x1, y1) + II(i, x2, y2) - II(i, x1, y2) - II(i, x2, y1);
        }
        else{
          if(x1 < 0 && y1 < 0){
            res(i, x, y) = II(i, x2, y2);
          }
          if(x1 < 0 && y1 >= 0){
            res(i, x, y) = II(i, x2, y2) - II(i, x2, y1);
          }
          if(x1 >= 0 && y1 < 0){
            res(i, x, y) = II(i, x2, y2) - II(i, x1, y2);
          }
        }
      }
      res(i) /= (long double) (BoxSize*BoxSize);
    }

  (*this) = res;
  }

  return res.move_to(*this);
//  return (*this);
}

CImgList<T> get_fastboxcar(int BoxSize, int idx) const
{
  return CImgList<T>(*this, false).fastboxcar(BoxSize, idx);
}



// Computes normalized complex coherence between two complex images
CImgList<T>& clxcoh(const CImgList<T> iml, int BoxSize = 1)
{
  CImgList<T> P1(1, (*this)(0).width(), (*this)(0).height(), (*this)(0).depth(), 1, 0.0);
  CImgList<T> P2(1, (*this)(0).width(), (*this)(0).height(), (*this)(0).depth(), 1, 0.0);
  P1 = (*this).get_clxmag().get_fastboxcar(BoxSize); 
  P2 = iml.get_clxmag().get_fastboxcar(BoxSize); 

  (*this).clxmulcj(iml).fastboxcar(BoxSize);
  CImgList<T> Q(1, (*this)(0).width(), (*this)(0).height(), (*this)(0).depth(), 1, 0.0);
  Q[0] = (P1[0].get_mul(P2[0])).sqrt();
  
  (*this).rldiv(Q);
  return *this;
}

CImgList<T> get_clxcoh(CImgList<T>& ClxImg, int BoxSize = 1)
{
  return CImgList<T>(*this, false).clxcoh(ClxImg, BoxSize); 
}


    //! Return a new image list corresponding to the complex tensor located at (x, y, z) of the current vector-valued list.
    CImgList<T> get_tensor_at(const unsigned int x, const unsigned int y=0, const unsigned int z=0) const 
{
  int S = (*this)(0).spectrum();
  if((*this).size()>2){
    std::cerr<<"get_tensor_at() for lists cannot be applied for lists of size > 2!\n"; 
    exit(1);
  }
  if(S!=3 && S!=6 && S!=21){
    std::cerr<<"get_tensor_at() for lists is handling only 2x2, 3x3 and 6x6 tensors!\n"; 
    exit(1);
  }
  CImgList<T> Res;
  if (S==21){
    Res.assign(2, 6, 6, 1, 1, 0); 
    for(int l = 0; l < 2; l++){
      int v = 0;
      cimg_forXY(Res(l), s, t){
        if(s + t <  6)
          Res(l, s + t, t) = (*this)(l, x, y, z, v++);
        if(s!=t){
          Res(0, t, s) = Res(0, s, t);  
          Res(1, t, s) = -Res(1, s, t);  
        }
      }
    }
  }
  if (S==6){
    Res.assign((*this).size(), 3, 3, 1, 1, 0);
      Res(0, 0, 0) = (*this)(0, x, y, z ,0); 
      Res(0, 1, 1) = (*this)(0, x, y, z, 3); 
      Res(0, 2, 2) = (*this)(0, x, y, z, 5); 
      
      Res(0, 1, 0) = Res(0, 0, 1) = (*this)(0, x, y, z, 1); 
      Res(0, 2, 0) = Res(0, 0, 2) = (*this)(0, x, y, z, 2); 
      Res(0, 2, 1) = Res(0, 1, 2) = (*this)(0, x, y, z, 4); 
      
      Res(1, 1, 0) = - (*this)(1, x, y, z, 1); 
      Res(1, 2, 0) = - (*this)(1, x, y, z, 2); 
      Res(1, 2, 1) = - (*this)(1, x, y, z, 4); 

      Res(1, 0, 1) = (*this)(1, x, y, z, 1); 
      Res(1, 0, 2) = (*this)(1, x, y, z, 2); 
      Res(1, 1, 2) = (*this)(1, x, y, z, 4); 
  }
  if (S==3){
    Res.assign((*this).size(), 2, 2, 1, 1, 0);
      Res(0, 0, 0) = (*this)(0, x, y, z,0); 
      Res(0, 1, 1) = (*this)(0, x, y, z, 2); 

      Res(0, 1, 0) = (*this)(0, x, y, z, 1); 
      Res(1, 1, 0) = - (*this)(1, x, y, z, 1); 
  }

  return Res;
}

//TODO: write for 2x2 and 6x6 tensors!!!!
    CImgList<T>& set_tensor_at(const CImgList<T> &Ts, const unsigned int x, const unsigned int y=0, const unsigned int z=0)
{
  int S = (*this)(0).spectrum();
  if((*this).size()>2){
    std::cerr<<"set_tensor_at() for lists cannot be applied for lists of size > 2!\n"; 
    exit(1);
  }
  if(S!=3 && S!=6 && S!=21){
    std::cerr<<"set_tensor_at() for lists is handling only 2x2, 3x3 and 6x6 tensors!\n"; 
    exit(1);
  }

  if (S==6){
    for(int l = 0; l < (int) (*this).size(); l++){
      (*this)(l, x, y, z ,0) = Ts(l, 0, 0); 
      (*this)(l, x, y, z, 1) = Ts(l, 1, 0); 
      (*this)(l, x, y, z, 2) = Ts(l, 2, 0); 
      (*this)(l, x, y, z, 3) = Ts(l, 1, 1); 
      (*this)(l, x, y, z, 4) = Ts(l, 2, 1); 
      (*this)(l, x, y, z, 5) = Ts(l, 2, 2); 
    }
  }
  return (*this);
}
   //! Return a new eigen matrix corresponding to the complex tensor located at (x, y, z) of the current vector-valued list.
//TODO: write for 2x2 and 6x6 tensors!!!!
//inline Eigen::Matrix3cf get_eigenmat_at(const unsigned int x, const unsigned int y=0, const unsigned int z=0) 
Eigen::Matrix3cf get_eigenmat_at(const unsigned int x, const unsigned int y=0, const unsigned int z=0) 
{
  Eigen::Matrix3cf Res;
  Res(0, 0) = std::complex<float> ((*this)(0, x, y, z ,0), 0.0); 
  Res(1, 1) = std::complex<float> ((*this)(0, x, y, z, 3), 0.0); 
  Res(2, 2) = std::complex<float> ((*this)(0, x, y, z, 5), 0.0); 

  Res(0, 1) = std::complex<float> ((*this)(0, x, y, z, 1), (*this)(1, x, y, z, 1)); 
  Res(0, 2) = std::complex<float> ((*this)(0, x, y, z, 2), (*this)(1, x, y, z, 2)); 
  Res(1, 2) = std::complex<float> ((*this)(0, x, y, z, 4), (*this)(1, x, y, z, 4)); 

  Res(1, 0) = std::conj(Res(0, 1)); 
  Res(2, 0) = std::conj(Res(0, 2)); 
  Res(2, 1) = std::conj(Res(1, 2)); 
  return Res;
}

// Same for any matrix dimension
// The name is the same, the argument is different (overloading would be impossible)
// Matrix dimension has to match CImg storage size (only upper elements are stored in the CImg)
// Please pre-allocate matrix M before using (thus we avoid control on the dimension).
void get_eigenmat_at(Eigen::MatrixXcf &M, const unsigned int x, const unsigned int y=0, const unsigned int z=0) const
{
  int d = M.rows(); 
  if(M.cols() != d) std::cerr<<"Eigen matrix must be square.\n";
  
  int cnt = 0;
  for(int i = 0; i < d; i++ ) 
    for(int j = i; j < d; j++) {
      M(i, j) = std::complex<float> ((*this)(0, x, y, z, cnt), (*this)(1, x, y, z, cnt)); 
      if(i != j)
        M(j, i) = std::complex<float> ((*this)(0, x, y, z, cnt), -(*this)(1, x, y, z, cnt)); 
      cnt++;
    }
}



//TODO: write for 2x2 and 6x6 tensors!!!!
//inline CImgList<T>& set_eigenmat_at(const Eigen::Matrix3cf &Ts, const unsigned int x, const unsigned int y=0, const unsigned int z=0)
CImgList<T>& set_eigenmat_at(const Eigen::Matrix3cf &Ts, const unsigned int x, const unsigned int y=0, const unsigned int z=0)
{
  (*this)(0, x, y, z ,0) = Ts(0, 0).real(); 
  (*this)(0, x, y, z, 1) = Ts(0, 1).real(); 
  (*this)(0, x, y, z, 2) = Ts(0, 2).real(); 
  (*this)(0, x, y, z, 3) = Ts(1, 1).real(); 
  (*this)(0, x, y, z, 4) = Ts(1, 2).real(); 
  (*this)(0, x, y, z, 5) = Ts(2, 2).real(); 

  (*this)(1, x, y, z ,0) = Ts(0, 0).imag(); 
  (*this)(1, x, y, z, 1) = Ts(0, 1).imag(); 
  (*this)(1, x, y, z, 2) = Ts(0, 2).imag(); 
  (*this)(1, x, y, z, 3) = Ts(1, 1).imag(); 
  (*this)(1, x, y, z, 4) = Ts(1, 2).imag(); 
  (*this)(1, x, y, z, 5) = Ts(2, 2).imag(); 
  return (*this);
}

// Same for any dimension with dynamic eigen matrices
// Please pre-allocate matrix M before using (thus we avoid control on the dimension).
CImgList<T>& set_eigenmat_at(const Eigen::MatrixXcf &M, const unsigned int x, const unsigned int y=0, const unsigned int z=0)
{
  int d = M.rows(); 
  if(M.cols() != d) std::cerr<<"Eigen matrix must be square.\n";
  
  int cnt = 0;
  for(int i = 0; i < d; i++ ) 
    for(int j = i; j < d; j++) {
      (*this)(0, x, y, z, cnt) = std::real(M(i, j)); 
      (*this)(1, x, y, z, cnt) = std::imag(M(i, j)); 
      cnt++;
    }

  return *this;
}




CImgList<T>& cov3_elt(const CImgList<T> &VecK, int ch1, int ch2)
{
  int ch_cov;
  if(ch1 == 0 && ch2 == 0) ch_cov = 0;
  if(ch1 == 0 && ch2 == 1) ch_cov = 1;
  if(ch1 == 0 && ch2 == 2) ch_cov = 2;
  if(ch1 == 1 && ch2 == 1) ch_cov = 3;
  if(ch1 == 1 && ch2 == 2) ch_cov = 4;
  if(ch1 == 2 && ch2 == 2) ch_cov = 5;

  if(ch1 != ch2){
    (*this)(0).get_shared_channel(ch_cov) = VecK(0).get_channel(ch1).get_mul(VecK(0).get_channel(ch2));
    (*this)(0).get_shared_channel(ch_cov) += VecK(1).get_channel(ch1).get_mul(VecK(1).get_channel(ch2));
    (*this)(1).get_shared_channel(ch_cov) = VecK(1).get_channel(ch1).get_mul(VecK(0).get_channel(ch2));
    (*this)(1).get_shared_channel(ch_cov) -= VecK(0).get_channel(ch1).get_mul(VecK(1).get_channel(ch2));
  }else{
    (*this)(0).get_shared_channel(ch_cov) = VecK(0).get_channel(ch1).get_sqr() + VecK(1).get_channel(ch1).get_sqr();
    (*this)(1).get_shared_channel(ch_cov) = CImg<T> (VecK(0), "xyz", 0);
  }
  return (*this);
}

 //! Return a new image list corresponding to the polarimetric covariance matrix. To be applied to a 3 elements scattering vector.
CImgList<T> get_cov3(int WinSiz = 3) const 
{
  int S = (*this)(0).spectrum();
  int DimX = (*this)(0).width();
  int DimY = (*this)(0).height();
  int DimZ = (*this)(0).depth();
  if((*this).size()!=2){
    std::cerr<<"get_cov3() cannot be applied for lists of size = 2!\n"; 
    exit(1);
  }
  if(S!=3){
    std::cerr<<"get_cov3() is only for complex 3 elements vectors!\n"; 
    exit(1);
  }
  CImgList<T> Res(2, DimX, DimY, DimZ, 6, 0);
  Res.cov3_elt((*this), 0, 0);
  Res.cov3_elt((*this), 0, 1);
  Res.cov3_elt((*this), 0, 2);
  Res.cov3_elt((*this), 1, 1);
  Res.cov3_elt((*this), 1, 2);
  Res.cov3_elt((*this), 2, 2);

  if(WinSiz > 1){
    Res.fastboxcar(WinSiz);
  }
  return Res;
}


CImgList<T>& k3_lex2pauli()
{
  CImgList<T> Buff((*this));
  (*this)(0).get_shared_channel(0) = (1.0/std::sqrt(2.0)) * (Buff(0).get_channel(0) + Buff(0).get_channel(2));  
  (*this)(1).get_shared_channel(0) = (1.0/std::sqrt(2.0)) * (Buff(1).get_channel(0) + Buff(1).get_channel(2));  

  (*this)(0).get_shared_channel(1) = (1.0/std::sqrt(2.0)) * (Buff(0).get_channel(0) - Buff(0).get_channel(2));  
  (*this)(1).get_shared_channel(1) = (1.0/std::sqrt(2.0)) * (Buff(1).get_channel(0) - Buff(1).get_channel(2));  

  (*this)(0).get_shared_channel(2) = Buff(0).get_channel(1);  
  (*this)(1).get_shared_channel(2) = Buff(1).get_channel(1);  

  return (*this);
} 

CImgList<T>& T3toC3()
{
  using namespace Eigen;
  Matrix3cf U_pl;
  U_pl(0, 0) = 1.0; U_pl(0, 1) = 1.0; U_pl(0, 2) = 0.0;
  U_pl(1, 0) = 0.0; U_pl(1, 1) = 0.0; U_pl(1, 2) = std::sqrt(2.0);
  U_pl(2, 0) = 1.0; U_pl(2, 1) = -1.0; U_pl(2, 2) = 0.0;
  U_pl /= std::sqrt(2.0);
  cimg_forXY((*this)(0), x, y){
    Matrix3cf T3 = (*this).get_eigenmat_at(x,y);
    Matrix3cf C3 = U_pl * T3 * U_pl.transpose();
    (*this).set_eigenmat_at(C3, x, y);
  }
  return (*this); 
}


CImgList<T>& get_T3toC3()
{
  return CImgList<T>(*this, false).T3toC3(); 
}

CImgList<T>& C3toT3()
{
  using namespace Eigen;
  Matrix3cf U_lp;
  U_lp(0, 0) = 1.0; U_lp(0, 1) = 0.0; U_lp(0, 2) = 1.0;
  U_lp(1, 0) = 1.0; U_lp(1, 1) = 0.0; U_lp(1, 2) = -1.0;
  U_lp(2, 0) = 0.0; U_lp(2, 1) = std::sqrt(2.0); U_lp(2, 2) = 0.0;
  U_lp /= std::sqrt(2.0);
  cimg_forXY((*this)(0), x, y){
    Matrix3cf C3 = (*this).get_eigenmat_at(x,y);
    Matrix3cf T3 = U_lp * C3 * U_lp.transpose();
    (*this).set_eigenmat_at(T3, x, y);
  }
  return (*this); 
}


CImgList<T>& get_C3toT3()
{
  return CImgList<T>(*this, false).C3toT3(); 
}

CImgList<T>& sum4look()
{
  int L = (*this).size();
  int W = (*this)(0).width();
  int H = (*this)(0).height();
  int D = (*this)(0).depth();
  int S = (*this)(0).spectrum();
  CImgList<T> Res(L, cimg::max(1, W/2), cimg::max(1, H/2), cimg::max(1, D/2), S, 0.0);

  cimglist_for(Res, l){
    CImg<T> I(2,2);
    CImg<T> Im ((*this)(l));

    cimg_forZC(Im, z, c){
      cimg_for2x2(Im, x, y, z, c, I, T){
        if(x%2 == 0 && y%2 == 0 && z%2 == 0 && x/2 < Res(l).width() && y/2 < Res(l).height() && z/2 < Res(l).depth())
          Res(l, x/2, y/2, z/2, c) = I.sum() / 4.0;
      }
    }
  }
  (*this).assign(Res);
  return *this;
  
}

CImgList<T> get_sum4look()
{
  return CImgList<T>(*this, false).sum4look(); 
}

// Computes the coefficient of variation of an image.
// If this is a multichannel (covariance) image, we compute CV for
// diagonal elements and make the average.
CImgList<T> get_coeff_var(int WinSiz = 7) const {
  int w = (*this)(0).width();
  int h = (*this)(0).height();
  int nc = (*this)(0).spectrum();
  CImg<T> res(w, h, 1, 1, 0.0);
  
  // finding the dimension of the matrix from 
  // the number of channels of the image (upper diagonal)
  int d = 0; 
  int nelt = nc; 
  while(nelt > 0) nelt-=(++d); // removing # of elements from each line of the triangular matrix and counting the number of removals (i.e. lines) 
  if(nelt!=0) std::cerr<<"Dimensions do not correspond to a square matrix.\n";

  // Calculating CV on diagonal channels
  int i = 0;
  int cnt = 0;
  while(i<nc){
    CImgList<T> ImgC((*this)(0).get_channel(i), (*this)(1).get_channel(i));
    CImg<T> sqrMeanImg(ImgC.get_fastboxcar(WinSiz)(0).get_sqr());
    CImg<T> varImg(ImgC.get_clxmag().get_fastboxcar(WinSiz)(0)-sqrMeanImg);
    res += varImg.get_div(sqrMeanImg);  
    i += (d - cnt++); // index of nex diagonal elt
  }
  res /= d;
  return CImgList(res);
}

CImgList<T>& coeff_var(int WinSiz=7) {
  return get_coeff_var(WinSiz).move_to(*this);
} 

// Computes span image of a N-dim covariance matrix
CImgList<T> get_span() const {
  int w = (*this)(0).width();
  int h = (*this)(0).height();
  int nc = (*this)(0).spectrum();
  CImg<T> res(w, h, 1, 1, 0.0);
  
  // finding the dimension of the matrix from 
  // the number of channels of the image (upper diagonal)
  int d = 0; 
  int nelt = nc; 
  while(nelt > 0) nelt-=(++d); // removing # of elements from each line of the triangular matrix and counting the number of removals (i.e. lines) 
  if(nelt!=0) std::cerr<<"Dimensions do not correspond to a square matrix.\n";

  // Summing diagonal channels.
  int i = 0;
  int cnt = 0;
  while(i<nc){
    CImgList<T> ImgC((*this)(0).get_channel(i), (*this)(1).get_channel(i));
    res += ImgC(0);  
    i += (d - cnt++); // index of nex diagonal elt
  }
  return CImgList(res);
}

CImgList<T>& span() {
  return get_span().move_to(*this);
} 

// Computes eigenvalues of a 3x3 image of hermitian matrices.
// TODO: assign and fill eigenvector image
void get_eigT3x3(CImgList<T> &val, CImgList<T> &vec)
{
  using namespace Eigen;
  //  using Eigen::Matrix3cf;
  int DimX = (*this)(0).width();
  int DimY = (*this)(0).height();
  int DimZ = (*this)(0).depth();
  int DimV = (*this)(0).spectrum();
  val.assign(1, DimX, DimY, DimZ, 3);
//  vec.assign(2, DimX, DimY, DimZ, 9);
  cimg_forXYZ((*this)(0), x, y, z){
    Matrix<std::complex<T>, 3, 3> M;
    M(0, 0) = std::complex<float>( (*this)(0,x,y,z,0) , (*this)(1,x,y,z,0) );
    M(1, 1) = std::complex<float>( (*this)(0,x,y,z,3) , (*this)(1,x,y,z,3) );
    M(2, 2) = std::complex<float>( (*this)(0,x,y,z,5) , (*this)(1,x,y,z,5) );
    M(0, 1) = std::complex<float>( (*this)(0,x,y,z,1) , (*this)(1,x,y,z,1) ); M(1, 0) = std::conj(M(0, 1));
    M(0, 2) = std::complex<float>( (*this)(0,x,y,z,2) , (*this)(1,x,y,z,2) ); M(2, 0) = std::conj(M(0, 2));
    M(1, 2) = std::complex<float>( (*this)(0,x,y,z,4) , (*this)(1,x,y,z,4) ); M(2, 1) = std::conj(M(1, 2));

 //    std::cout<<M<<'\n';
    SelfAdjointEigenSolver<Matrix<std::complex<T>, 3, 3> > es(M);
//    std::cout<<es.eigenvalues().reverse()<<'\n';
    Matrix<T, 3, 1> Va(es.eigenvalues().reverse());
    val(0, x, y, z, 0) = Va(0);
    val(0, x, y, z, 1) = Va(1);
    val(0, x, y, z, 2) = Va(2);
  }
//  val(0).get_channel(0).get_cut(0, 2.5*val(0).get_channel(0).mean()).display();
//  val(0).get_channel(1).get_cut(0, 2.5*val(0).get_channel(1).mean()).display();
//  val(0).get_channel(2).get_cut(0, 2.5*val(0).get_channel(2).mean()).display();
  //  return Res;
}

// Computes the H-alpha decomposition for polarimetric T matrices.
void get_HAlpha(CImgList<T> &H, CImgList<T>& Alpha)
{
  using namespace Eigen;
  int DimX = (*this)(0).width();
  int DimY = (*this)(0).height();
  int DimZ = (*this)(0).depth();
  H.assign(1, DimX, DimY, DimZ, 1, 0.0);
  Alpha.assign(1, DimX, DimY, DimZ, 1, 0.0);

#pragma omp parallel for
  cimg_forXYZ((*this)(0), x, y, z){
    Matrix<std::complex<T>, 3, 3> M;
    M(0, 0) = std::complex<T>( (*this)(0,x,y,z,0) , (*this)(1,x,y,z,0) );
    M(1, 1) = std::complex<T>( (*this)(0,x,y,z,3) , (*this)(1,x,y,z,3) );
    M(2, 2) = std::complex<T>( (*this)(0,x,y,z,5) , (*this)(1,x,y,z,5) );
    M(0, 1) = std::complex<T>( (*this)(0,x,y,z,1) , (*this)(1,x,y,z,1) ); M(1, 0) = std::conj(M(0, 1));
    M(0, 2) = std::complex<T>( (*this)(0,x,y,z,2) , (*this)(1,x,y,z,2) ); M(2, 0) = std::conj(M(0, 2));
    M(1, 2) = std::complex<T>( (*this)(0,x,y,z,4) , (*this)(1,x,y,z,4) ); M(2, 1) = std::conj(M(1, 2));

    SelfAdjointEigenSolver<Matrix<std::complex<T>, 3, 3> > es(M);
    Matrix<T, 3, 1> Va(es.eigenvalues());
    Matrix<std::complex<T>, 3, 3> Ve(es.eigenvectors());
    Va /= Va.sum();

    float Hpix = -Va(0) * std::log(Va(0)) / std::log(3.0);
    Hpix += -Va(1) * std::log(Va(1)) / std::log(3.0);
    Hpix += -Va(2) * std::log(Va(2)) / std::log(3.0);
    H(0, x, y, z) = Hpix;
    float Apix = std::acos(std::abs(Ve(0,0))) * Va(0); 
    Apix += std::acos(std::abs(Ve(0,1))) * Va(1); 
    Apix += std::acos(std::abs(Ve(0,2))) * Va(2); 
    Alpha(0, x, y, z) = Apix;
//    if(x==635 && y==465){ 
//        std::cout<<"\nSingularity detected\n";
//        std::cout<< "x: "<<x<<'\n';
//        std::cout<< "y: "<<y<<'\n';
//        std::cout<< "T: "<<M<<'\n';
//        std::cout<< "Va: "<<Va<<'\n';
//        std::cout<< "Ve: "<<Ve<<'\n';
//        
//        std::cin.ignore();
//      }

  }

}

CImgList<T> get_scatterplot2D(int NBinX, int NBinY)
{
  CImg<T> ImX((*this)(0));
  CImg<T> ImY((*this)(1));
  
  ImX.normalize(0, NBinX-1);
  ImY.normalize(0, NBinY-1);
  CImg<int> SPlot(NBinX, NBinY, 1, 1, 0); 
  CImgList<unsigned char> Disp(1, NBinX, NBinY, 1, 3, 0); 
  cimg_forXY(ImX, x, y){
    SPlot((int)ImX(x, y), (int)ImY(x, y)) += 1; 
  }

  CImg<> CMap;
  CMap = CImg<>::jet_LUT256();
  CMap(0, 0, 0, 0) = 255;
  CMap(0, 0, 0, 1) = 255;
  CMap(0, 0, 0, 2) = 255;

  SPlot.normalize(0,255); 
  cimg_forXY(Disp(0), x, y){
    Disp(0, x, y, 0, 0) = CMap(0, SPlot(x, y), 0, 0);   
    Disp(0, x, y, 0, 1) = CMap(0, SPlot(x, y), 0, 1);   
    Disp(0, x, y, 0, 2) = CMap(0, SPlot(x, y), 0, 2);   
  }
//  Disp(0).mirror('y').display();
  Disp(0).mirror('y');

  return Disp;
}

// Makes a scatter plot of H and alpha parameters with curves for the limits of validity
CImgList<T> get_scatterHAlpha(int SizX, int SizY)
{
  CImg<T> ImX((*this)(0));
  CImg<T> ImY((*this)(1));
  
  const unsigned char black[] = {0,0,0};

  CImg<int> SPlot(SizX, SizY, 1, 1, 0); 
  CImgList<unsigned char> Disp(1, SizX, SizY, 1, 3, 0); 
  
  cimg_forXY(ImX, x, y){
    if(ImX(x, y) > 0.0 && ImX(x, y) < 1.0 
        && ImY(x, y) > 0.0 && ImY(x, y) < 0.5*M_PI)
    SPlot(floor(ImX(x, y) * (SizX-1)), floor(ImY(x, y) * (SizY-1) / (0.5*M_PI))) += 1; 
  }

    CImg<> CMap;
  CMap = CImg<>::jet_LUT256();
  CMap(0, 0, 0, 0) = 255;
  CMap(0, 0, 0, 1) = 255;
  CMap(0, 0, 0, 2) = 255;

  SPlot.normalize(0,255); 
  cimg_forXY(Disp(0), x, y){
    Disp(0, x, y, 0, 0) = CMap(0, SPlot(x, y), 0, 0);   
    Disp(0, x, y, 0, 1) = CMap(0, SPlot(x, y), 0, 1);   
    Disp(0, x, y, 0, 2) = CMap(0, SPlot(x, y), 0, 2);   
  }
   
  float h1p, a1p, h2p, a2p;
  a2p = 1.0;
  h2p = -logf(0.5) / logf(3.0); 
  for(int j = 0; j < SizX; j++){
    float m = (float) j / (float) (SizX-1); 
    float p1, p2, p3, h1, a1;
    float h2, a2;
    // Curve I
    p1 = 1.0 / (1.0 + 2.0*m);
    p2 = m / (1.0 + 2.0*m); p3 = p2;
    if(m > 0.0)
      h1 = (1.0 / logf(3.0)) * (-p1*logf(p1) - p2*logf(p2) - p3*logf(p3));
    else
      h1 = 0.0;
    a1 = p2 + p3;
    if(j > 0){
      int x0 = floor (h1p*(SizX-1)); int y0 = floor (a1p*(SizY-1));
      int x1 = floor (h1*(SizX-1)); int y1 = floor (a1*(SizY-1));
      Disp(0).draw_line(x0, y0, x1, y1, black);
    }
    h1p = h1; a1p = a1;
    
    // Curve II
    if(m > 0.5){
      p1 = (2*m - 1.0) / (1.0 + 2.0*m); 
      p2 = 1 / (1.0 + 2.0*m); 
      p3 = p2; 

      h2 = (1.0 / logf(3.0)) * (-p1*logf(p1) - p2*logf(p2) - p3*logf(p3));
      a2 = p2 + p3;
      int x0 = floor (h2p*(SizX-1)); int y0 = floor (a2p*(SizY-1));
      int x1 = floor (h2*(SizX-1)); int y1 = floor (a2*(SizY-1));
      Disp(0).draw_line(x0, y0, x1, y1, black);
      h2p = h2; a2p = a2;
    }

  }

  Disp(0).mirror('y');

  return Disp;
}

// Makes a colormap corresponding to an HSV representation of alpha entropy
CImgList<T> get_coldiagramHAlpha(int SizX, int SizY)
{
  CImg<T> ImX((*this)(0));
  CImg<T> ImY((*this)(1));
  
  const float black[] = {0,0,0};

  CImgList<float> Disp(1, SizX, SizY, 1, 3, 0); 
  

  cimg_forXY(Disp(0), x, y){
    float cx = (float) x / float (SizX-1);
    float cy = (float) y / float (SizY-1);
    Disp(0, x, y, 0, 0) =  (360 - 240) * cy + 240;  
    Disp(0, x, y, 0, 1) = 1.0 - cx; 
//    Disp(0, x, y, 0, 2) = 0.05 * cx + 0.95;
    Disp(0, x, y, 0, 2) = 1.0;
  }
  Disp(0).HSVtoRGB();
   
  float h1p, a1p, h2p, a2p;
  a2p = 1.0;
  h2p = -logf(0.5) / logf(3.0); 

  for(int j = 0; j < SizX; j++){
    float m = (float) j / (float) (SizX-1); 
    float p1, p2, p3, h1, a1;
    float h2, a2;
    // Curve I
    p1 = 1.0 / (1.0 + 2.0*m);
    p2 = m / (1.0 + 2.0*m); p3 = p2;
    if(m > 0.0)
      h1 = (1.0 / logf(3.0)) * (-p1*logf(p1) - p2*logf(p2) - p3*logf(p3));
    else
      h1 = 0.0;
    a1 = p2 + p3;
    if(j > 0){
      int x0 = floor (h1p*(SizX-1)); int y0 = floor (a1p*(SizY-1));
      int x1 = floor (h1*(SizX-1)); int y1 = floor (a1*(SizY-1));
      Disp(0).draw_line(x0, y0, x1, y1, black);
    }
    h1p = h1; a1p = a1;
    
    // Curve II
    if(m > 0.5){
      p1 = (2*m - 1.0) / (1.0 + 2.0*m); 
      p2 = 1 / (1.0 + 2.0*m); 
      p3 = p2; 

      h2 = (1.0 / logf(3.0)) * (-p1*logf(p1) - p2*logf(p2) - p3*logf(p3));
      a2 = p2 + p3;
      int x0 = floor (h2p*(SizX-1)); int y0 = floor (a2p*(SizY-1));
      int x1 = floor (h2*(SizX-1)); int y1 = floor (a2*(SizY-1));

      Disp(0).draw_line(x0, y0, x1, y1, black);
      h2p = h2; a2p = a2;
    }
  }

  Disp(0).mirror('y');

  return Disp;
}



// Experimental colormap for H-Alpha parameters.
CImgList<T> get_colmapHAlpha()
{
 if((*this).size()!=2){
    std::cerr<<"get_colmapHAlpha(): To be used with H Alpha parameters.\n";
  }
  CImg<> CMap;
  CImgList<T> Disp(1, (*this)(0).width(),(*this)(0).height(),(*this)(0).depth(),3,1.0);
  Disp(0).get_shared_channel(1) = 1.0 - (*this)(0);
  Disp(0).get_shared_channel(2) = 0.05*((*this)(0)) + 0.95;
  Disp(0).get_shared_channel(0) = (360 - 240) * ((*this)(1)/(0.5*M_PI)) + 240;
  Disp(0).HSVtoRGB();
  return Disp;

}

//std::complex<float> logfn(std::complex<float> x, int)
//{
//  return std::log(x);
//}


CImgList<T>& logT3x3()
{
  using namespace Eigen;

  cimg_forXYZ((*this)(0), x, y, z){
    Matrix<std::complex<T>, 3, 3> M;
//    MatrixXcf M(3,3);
//    Matrix3cf M;
//    M(0, 0) = std::complex<float>( (*this)(0,x,y,z,0) , (*this)(1,x,y,z,0) );
//    M(1, 1) = std::complex<float>( (*this)(0,x,y,z,3) , (*this)(1,x,y,z,3) );
//    M(2, 2) = std::complex<float>( (*this)(0,x,y,z,5) , (*this)(1,x,y,z,5) );
//    M(0, 1) = std::complex<float>( (*this)(0,x,y,z,1) , (*this)(1,x,y,z,1) ); M(1, 0) = std::conj(M(0, 1));
//    M(0, 2) = std::complex<float>( (*this)(0,x,y,z,2) , (*this)(1,x,y,z,2) ); M(2, 0) = std::conj(M(0, 2));
//    M(1, 2) = std::complex<float>( (*this)(0,x,y,z,4) , (*this)(1,x,y,z,4) ); M(2, 1) = std::conj(M(1, 2));
    M(0, 0) = std::complex<T>( (*this)(0,x,y,z,0) , (*this)(1,x,y,z,0) );
    M(1, 1) = std::complex<T>( (*this)(0,x,y,z,3) , (*this)(1,x,y,z,3) );
    M(2, 2) = std::complex<T>( (*this)(0,x,y,z,5) , (*this)(1,x,y,z,5) );
    M(0, 1) = std::complex<T>( (*this)(0,x,y,z,1) , (*this)(1,x,y,z,1) ); M(1, 0) = std::conj(M(0, 1));
    M(0, 2) = std::complex<T>( (*this)(0,x,y,z,2) , (*this)(1,x,y,z,2) ); M(2, 0) = std::conj(M(0, 2));
    M(1, 2) = std::complex<T>( (*this)(0,x,y,z,4) , (*this)(1,x,y,z,4) ); M(2, 1) = std::conj(M(1, 2));


    SelfAdjointEigenSolver<Matrix<std::complex<T>, 3, 3> > es(M);
//    ComplexEigenSolver<Matrix<std::complex<T>, 3, 3> > es(M);
    Matrix<T , 3, 3> Va;
//    Matrix<std::complex<T> , 3, 3> Va;
    Va.setZero();
    Matrix<std::complex<T>, 3, 3> Ve(es.eigenvectors());

//    std::cout<<"MatOrig\n"<<M<<"\n\n";
    
    Va(0,0) = std::log(es.eigenvalues()(0));
    Va(1,1) = std::log(es.eigenvalues()(1));
    Va(2,2) = std::log(es.eigenvalues()(2));
//    Matrix<std::complex<T>, 3, 3> L;
//    MatrixXcf L(3,3);
//    M.matrixFunction(logfn); 
//    Matrix3cf L = M.log(); 
    M = Ve*(Va*Ve.adjoint());
   
//    std::cout<<"Vec\n"<<Ve<<"\n\n";
//    std::cout<<"Vec sq norm\n"<<Ve.squaredNorm()<<"\n\n";
//    std::cout<<"Val\n"<<es.eigenvalues()<<"\n\n";
//    std::cout<<"LogVal\n"<<Va<<"\n\n";
//    std::cout<<"LogMat\n"<<M<<"\n\n";
//    std::cin.ignore();

    (*this)(0,x,y,z,0) = std::real(M(0, 0)); (*this)(1,x,y,z,0) = 0.0;
    (*this)(0,x,y,z,3) = std::real(M(1, 1)); (*this)(1,x,y,z,3) = 0.0;
    (*this)(0,x,y,z,5) = std::real(M(2, 2)); (*this)(1,x,y,z,5) = 0.0;
    (*this)(0,x,y,z,1) = std::real(M(0, 1)); (*this)(1,x,y,z,1) = std::imag(M(0, 1)); 
    (*this)(0,x,y,z,2) = std::real(M(0, 2)); (*this)(1,x,y,z,2) = std::imag(M(0, 2)); 
    (*this)(0,x,y,z,4) = std::real(M(1, 2)); (*this)(1,x,y,z,4) = std::imag(M(1, 2)); 

  }
  (*this).display();
    return (*this);
}

CImgList<T> get_logT3x3(){
  return CImgList<T>(*this, false).logT3x3();
}

CImgList<T>& expT3x3()
{
  using namespace Eigen;
  cimg_forXYZ((*this)(0), x, y, z){
    Matrix<std::complex<T>, 3, 3> M;
    M(0, 0) = std::complex<T>( (*this)(0,x,y,z,0) , (*this)(1,x,y,z,0) );
    M(1, 1) = std::complex<T>( (*this)(0,x,y,z,3) , (*this)(1,x,y,z,3) );
    M(2, 2) = std::complex<T>( (*this)(0,x,y,z,5) , (*this)(1,x,y,z,5) );
    M(0, 1) = std::complex<T>( (*this)(0,x,y,z,1) , (*this)(1,x,y,z,1) ); M(1, 0) = std::conj(M(0, 1));
    M(0, 2) = std::complex<T>( (*this)(0,x,y,z,2) , (*this)(1,x,y,z,2) ); M(2, 0) = std::conj(M(0, 2));
    M(1, 2) = std::complex<T>( (*this)(0,x,y,z,4) , (*this)(1,x,y,z,4) ); M(2, 1) = std::conj(M(1, 2));

//    SelfAdjointEigenSolver<Matrix<std::complex<T>, 3, 3> > es(M);
    ComplexEigenSolver<Matrix<std::complex<T>, 3, 3> > es(M);
    Matrix<std::complex<T> , 3, 3> Va;
//    Matrix<T , 3, 3> Va;
    Va.setZero();
    Matrix<std::complex<T>, 3, 3> Ve(es.eigenvectors());

    Va(0,0) = std::exp(es.eigenvalues()(0));
    Va(1,1) = std::exp(es.eigenvalues()(1));
    Va(2,2) = std::exp(es.eigenvalues()(2));
    M = Ve*(Va*Ve.adjoint());
   
    (*this)(0,x,y,z,0) = std::real(M(0, 0)); (*this)(1,x,y,z,0) = 0.0;
    (*this)(0,x,y,z,3) = std::real(M(1, 1)); (*this)(1,x,y,z,3) = 0.0;
    (*this)(0,x,y,z,5) = std::real(M(2, 2)); (*this)(1,x,y,z,5) = 0.0;
    (*this)(0,x,y,z,1) = std::real(M(0, 1)); (*this)(1,x,y,z,1) = std::imag(M(0, 1)); 
    (*this)(0,x,y,z,2) = std::real(M(0, 2)); (*this)(1,x,y,z,2) = std::imag(M(0, 2)); 
    (*this)(0,x,y,z,4) = std::real(M(1, 2)); (*this)(1,x,y,z,4) = std::imag(M(1, 2)); 

  }
    return (*this);
}

CImgList<T> get_expT3x3(){
  return CImgList<T>(*this, false).expT3x3();
}





// Computes the Frobenius norm of a complex matrix
//inline T FrobNorm()
T FrobNorm()
{
  T Norm2 = 0.0;
  cimglist_for(*this, l){
    cimg_for((*this)[l], ptrs, T){
      Norm2 += cimg::sqr(*ptrs); 
//      Norm2 += (*this)(l, x, y)*(*this)(l, x, y); 
    }
  }
  return std::sqrt(Norm2);
}

// Computes the Frobenius norm of an image of complex tensors 
CImgList<T> get_frobnorm() const
{
  CImg<T> Img((*this)(0), "xy", 0.0);
  cimg_forXY(Img, x, y){
    CImgList<T> Tens = (*this).get_tensor_at(x, y);
    Img(x, y) = Tens.FrobNorm();
  }
  return CImgList<T> (Img);
}

// Computes the squared Frobenius norm of a complex matrix
//inline T FrobNorm()
T SqrFrobNorm()
{
  T Norm2 = 0.0;
  cimglist_for(*this, l){
    cimg_for((*this)[l], ptrs, T){
      Norm2 += cimg::sqr(*ptrs); 
    }
  }
  return Norm2;
}

// Computes the squared Frobenius norm of an image of complex tensors 
CImgList<T> get_sqrfrobnorm() const
{
  CImg<T> Img((*this)(0), "xy", 0.0);
  cimg_forXY(Img, x, y){
    CImgList<T> Tens = (*this).get_tensor_at(x, y);
    Img(x, y) = Tens.SqrFrobNorm();
  }
  return CImgList<T> (Img);
}


//#endif
