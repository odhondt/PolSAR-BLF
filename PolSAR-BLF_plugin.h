/*
 #
 #  File        : PolSAR-BLF_plugin.h
 #                ( C++ header file - CImg plug-in )
 #
 #  Description : Header file for the method described in the paper 
 #                'Iterative Bilateral Filtering of Polarimetric SAR Data'. 
 #                This plugin is using the CImg library.
 #                ( http://cimg.sourceforge.net )
 #
 #  Copyright   : Olivier D'Hondt
 #                (https://sites.google.com/site/dhondtolivier/)
 #                (https://github.com/odhondt)
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

Eigen::Matrix3cf LogMat(const Eigen::Matrix3cf &A)
{
  using namespace::Eigen;
  SelfAdjointEigenSolver<Matrix3cf> es(A);
  Vector3f Va = es.eigenvalues();
  Matrix3cf Ve = es.eigenvectors(); 

  // Recomposing the matrix
  Va = Va.array().log();
  Matrix3f M = Va.asDiagonal(); 
//  return Ve*M*(Ve.adjoint());
  return Ve*M*Ve.adjoint();
}

float SymKullbackDist(const Eigen::Matrix3cf &A, const Eigen::Matrix3cf &B) 
{
  using namespace::Eigen;

  float dist = std::real( (A.inverse()*B + B.inverse()*A).trace() );
  
  dist = 0.5 * dist - 3.0;
  return (dist>0.0)?dist:0.0;
}

// --- BILATERAL FILTER WITH DIFFERENT MATRIX SIMILARITIES ---


// Bilateral filter with affine invariant similarity 
CImgList<T>& polsar_blf_ai(float gammaS = 2.0, float gammaR = 2.0, bool TRICK=false, bool FLATW=false, CImgDisplay *disp = 0) 
{
  CImgList<T> Res((*this)(0), (*this)(1));
  Res(0).fill(1.0);
  Res(1).fill(1.0);

  float gr2 = gammaR * gammaR;
  float gs2 = gammaS * gammaS;
 
   // Calculating the window size for the gaussian weights
  int H = ceil(std::sqrt(3.0)*gammaS); 

  // spatial Gaussian weights
  CImg<T> WIm(2*H+1, 2*H+1, 1, 1, 1.0);
  if(!FLATW)
    cimg_forXY(WIm, x, y){
      int s = x - H;
      int t = y - H;
      WIm(x, y) = std::exp(- (s*s + t*t) / gs2 );
    }else{
      cimg_forXY(WIm, x, y){
        int s = x - H;
        int t = y - H;
        if((s*s + t*t) <= (H+1)*(H+1))
          WIm(x, y) = 1.0;
        else
          WIm(x, y) = 0.0;
      }
    }


  // Making an image of matrices to speed-up computations
  using namespace Eigen;
  CImg<Matrix3cf> Data((*this)(0), "xy");
  GeneralizedSelfAdjointEigenSolver<Matrix3cf> es;
  cimg_forXY((*this)(0), x, y){
    Data(x, y) = (*this).get_eigenmat_at(x, y);
  }


  Matrix3cf TsFilt;

  int cnt = 0; 
  if(!TRICK){
#pragma omp parallel for firstprivate(es, TsFilt)
    cimg_forXY(Res[0], x, y){
#pragma omp atomic
      ++cnt;
      float SumWeight = 0;
      TsFilt.setZero();

      const int tmin = cimg::max(y-H, 0),
            tmax = cimg::min(y+H, Res[0].height()-1),
            smin = cimg::max(x-H, 0),
            smax = cimg::min(x+H, Res[0].width()-1);

      for(int t = tmin; t <= tmax; ++t)
        for(int s = smin; s <= smax; ++s) {
          es.compute(Data(s, t), Data(x, y), EigenvaluesOnly);
          const float D = es.eigenvalues().array().log().abs2().sum();
          float W = 0.0;
          if(!std::isnan(D)) W = std::exp(- D / gr2) * WIm(s - x + H, t - y + H);
          TsFilt += W*Data(s, t);
          SumWeight += W;
        }

      if(SumWeight > 1.0e-10) TsFilt /= SumWeight;
      else TsFilt = Data(x, y);
      Res.set_eigenmat_at(TsFilt, x, y);

      if(cnt%9000 == 0 && disp) {
#pragma omp critical
        disp->display(Res.get_colcov3()).set_title("Filtering.");
      }
    }
  }else{
//#pragma omp parallel for
#pragma omp parallel for firstprivate(es, TsFilt)
    cimg_forXY(Res[0], x, y){
#pragma omp atomic
      ++cnt;
      TsFilt.setZero();
      // Computing the weights
      float SumWeight = 0.0, WRMax = 0.0;

      const int tmin = cimg::max(y-H, 0),
            tmax = cimg::min(y+H, Res[0].height()-1),
            smin = cimg::max(x-H, 0),
            smax = cimg::min(x+H, Res[0].width()-1);

      for(int t = tmin; t <= tmax; ++t) {
        for(int s = smin; s <= smax; ++s) {
          if(s != x || t != y) {
            es.compute(Data(s, t), Data(x, y), EigenvaluesOnly);
            const float D = es.eigenvalues().array().log().abs2().sum();
            float WR = 0.0;
            if(!std::isnan(D)){
              WR = std::exp(- D / gr2);
              if(WR > WRMax && WR < 1.0)
                WRMax = WR;
            }
            const float W = WR * WIm(s - x + H, t - y + H);
            TsFilt += W * Data(s, t);
            SumWeight += W;
          }
        }
      }
      // Adding central coefficient
      TsFilt += WRMax*Data(x, y);
      SumWeight += WRMax;
      if(SumWeight > 1.0e-10) TsFilt /= SumWeight;
      else TsFilt = Data(x, y) ;
      Res.set_eigenmat_at(TsFilt, x, y);

      if(cnt%9000 == 0 && disp) {
#pragma omp critical
        disp->display(Res.get_colcov3()).set_title("Filtering.");
      }
    }

  }

  return Res.move_to(*this);
}

CImgList<T> get_polsar_blf_ai(float gammaS = 2.0, float gammaR = 2.0, bool TRICK=false, bool FLATW=false, CImgDisplay *disp = 0) 
{
  return CImgList<T>(*this, false).polsar_blf_ai(gammaS, gammaR, TRICK, FLATW, disp);
}

// Bilateral filter with log Euclidean similarity
CImgList<T>& polsar_blf_le(float gammaS = 2.0, float gammaR = 2.0, bool TRICK=false, bool FLATW=false, CImgDisplay *disp = 0) 
{

  CImgList<T> Res((*this)(0), (*this)(1));
  Res(0).fill(1.0);
  Res(1).fill(1.0);
  float gr2 = gammaR * gammaR;
  float gs2 = gammaS * gammaS;

  // Calculating the window size for the gaussian weights
  int H = ceil(std::sqrt(3.0)*gammaS); 

  CImg<T> WIm(2*H+1, 2*H+1, 1, 1, 1.0);
  if(!FLATW) {
    cimg_forXY(WIm, x, y) {
      int s = x - H;
      int t = y - H;
      WIm(x, y) = std::exp(- (s*s + t*t) / gs2 );
    }
  }else{
      cimg_forXY(WIm, x, y) {
        int s = x - H;
        int t = y - H;
        if((s*s + t*t) <= (H+1)*(H+1))
          WIm(x, y) = 1.0;
        else
          WIm(x, y) = 0.0;
      }
    }

  // Making an image of matrices and pre-computing log to speed-up computations
  using namespace Eigen;
  CImg<Matrix3cf> Data((*this)(0), "xy");
  CImg<Matrix3cf> logData((*this)(0), "xy");
#pragma omp parallel for
  cimg_forXY((*this)(0), x, y){
    Data(x, y) = (*this).get_eigenmat_at(x, y);
    logData(x, y) = LogMat(Data(x, y));
  }


  Matrix3cf TsFilt;
  int cnt = 0; 
  if(!TRICK){
#pragma omp parallel for firstprivate(TsFilt)
    cimg_forXY(Res[0], x, y){
#pragma omp atomic
      cnt++;
      float SumWeight = 0;

      Matrix3cf TsFilt;
      TsFilt.setZero();

      const int tmin = cimg::max(y-H, 0),
            tmax = cimg::min(y+H, Res[0].height()-1),
            smin = cimg::max(x-H, 0),
            smax = cimg::min(x+H, Res[0].width()-1);

      // Computing the weights
      for(int t = tmin; t <= tmax; ++t)
        for(int s = smin; s <= smax; ++s) {
          //const float D = (LogMat(Data(x, y)) - LogMat(Data(s, t))).squaredNorm();
          const float D = (logData(x, y) - logData(s, t)).squaredNorm();
          float W = 0.0;
          if(!std::isnan(D)) W = std::exp(- D / gr2) * WIm(s - x + H, t - y + H);
          TsFilt += W*Data(s, t);
          SumWeight += W;
        }
      //      TsFilt /= SumWeight;
      if(SumWeight > 1.0e-10) TsFilt /= SumWeight;
      else TsFilt = Data(x, y);
      Res.set_eigenmat_at(TsFilt, x, y);

      if(cnt%9000 == 0 && disp) {
#pragma omp critical
        disp->display(Res.get_colcov3()).set_title("Filtering.");
      }
    }
  }
  else{
    std::cout<<"prout\n";
#pragma omp parallel for firstprivate(TsFilt)
    cimg_forXY(Res[0], x, y){
      cnt++;
      float SumWeight = 0.0, WRMax = 0.0;
      TsFilt.setZero();

      const int tmin = cimg::max(y-H, 0),
            tmax = cimg::min(y+H, Res[0].height()-1),
            smin = cimg::max(x-H, 0),
            smax = cimg::min(x+H, Res[0].width()-1);

      // Computing the weights
      for(int t = tmin; t <= tmax; ++t) {
        for(int s = smin; s <= smax; ++s) {
          if(s != x || t != y) {
            //float D = (LogMat(Data(x, y)) - LogMat(Data(s, t))).squaredNorm();
            const float D = (logData(x, y) - logData(s, t)).squaredNorm();
            float WR = 0.0;
            if(!std::isnan(D)){
              WR = std::exp(- D / gr2);
              if(WR > WRMax && WR < 1.0)
                WRMax = WR;
            }
            const float W = WR * WIm(s - x + H, t - y + H);
            TsFilt += W * Data(s, t);
            SumWeight += W;
          }
        }
      }
      TsFilt += WRMax*Data(x, y);
      SumWeight += WRMax;
      if(SumWeight > 1.0e-10) TsFilt /= SumWeight;
      else TsFilt = Data(x, y);
      Res.set_eigenmat_at(TsFilt, x, y);
      
      if(cnt%9000 == 0 && disp) {
#pragma omp critical
        disp->display(Res.get_colcov3()).set_title("Filtering.");
      }
    }

  }

  return Res.move_to(*this);
}

CImgList<T> get_polsar_blf_le(float gammaS = 2.0, float gammaR = 2.0, bool TRICK=false, bool FLATW=false, CImgDisplay *disp = 0) 
{
  return CImgList<T>(*this, false).polsar_blf_le(gammaS, gammaR, TRICK, FLATW, disp);
}

// Bilateral with kullback leibler similarity ********
CImgList<float>& polsar_blf_kl(float gammaS = 2.0, float gammaR = 2.0,bool TRICK=false, bool FLATW=false, CImgDisplay *disp = 0) 
{

  CImgList<T> Res((*this)(0), (*this)(1));
  Res(0).fill(1.0);
  Res(1).fill(1.0);

  float gr2 = gammaR * gammaR;
  float gs2 = gammaS * gammaS;

  // Calculating the window size for the gaussian weights
  int H = ceil(std::sqrt(3.0)*gammaS); 

  // spatial Gaussian weights
  CImg<T> WIm(2*H+1, 2*H+1, 1, 1, 1.0);
  if(!FLATW) { 
    cimg_forXY(WIm, x, y) {
      int s = x - H;
      int t = y - H;
      WIm(x, y) = std::exp(- (s*s + t*t) / gs2 );
    }
  } else {
    cimg_forXY(WIm, x, y) {
      int s = x - H;
      int t = y - H;
      if((s*s + t*t) <= (H+1)*(H+1))
        WIm(x, y) = 1.0;
      else
        WIm(x, y) = 0.0;
    }
  }

//  WIm.display();

  // Making an image of matrices to speed-up computations
  using namespace Eigen;
  CImg<Matrix3cf> Data((*this)(0), "xy");
  cimg_forXY((*this)(0), x, y){
    Data(x, y) = (*this).get_eigenmat_at(x, y);
  }

  Matrix3cf TsFilt;
  
  CImg<float> DImg((*this)(0), "xy", 0.0);

  int cnt = 0; 
  if(!TRICK){
#pragma omp parallel for firstprivate(TsFilt)
    cimg_forXY(Res[0], x, y){
      ++cnt;
      float SumWeight = 0;
      TsFilt.setZero();
      const int tmin = cimg::max(y-H, 0),
            tmax = cimg::min(y+H, Res[0].height()-1),
            smin = cimg::max(x-H, 0),
            smax = cimg::min(x+H, Res[0].width()-1);

      for(int t = tmin; t <= tmax; ++t)
        for(int s = smin; s <= smax; ++s) {
          const float D = SymKullbackDist(Data(x, y), Data(s, t));
          float W = 0.0;
          if(!std::isnan(D)) W = std::exp(- D / gr2) * WIm(s - x + H, t - y + H);
          TsFilt += W*Data(s, t);
          SumWeight += W;
        }

      if(SumWeight > 1e-10) TsFilt /= SumWeight;
      else TsFilt = Data(x, y);
      Res.set_eigenmat_at(TsFilt, x, y);

      if(cnt%9000 == 0 && disp) {
#pragma omp critical
        disp->display(Res.get_colcov3()).set_title("Filtering.");
      }
    }
  }else{
#pragma omp parallel for firstprivate(TsFilt)
    cimg_forXY(Res[0], x, y){
#pragma omp atomic
      ++cnt;
      TsFilt.setZero();
      // Computing the weights
      float SumWeight = 0.0, WRMax = 0.0;

      const int tmin = cimg::max(y-H, 0),
            tmax = cimg::min(y+H, Res[0].height()-1),
            smin = cimg::max(x-H, 0),
            smax = cimg::min(x+H, Res[0].width()-1);
      CImg<float> Img(2*H+1,2*H+1,1,1,0);
      for(int t = tmin; t <= tmax; ++t) {
        for(int s = smin; s <= smax; ++s) {

          if(s != x || t != y) {
            const float D = SymKullbackDist(Data(x, y), Data(s, t));
            float WR = 0.0;

            if(!std::isnan(D) && D >= 0.0) {
              WR = std::exp(- D / gr2);
              if(WR > WRMax && WR < 1.0) WRMax = WR;
            }

            const float W = WR * WIm(s - x + H, t - y + H);
            Img(s - x + H, t - y + H) = 1;
            TsFilt += W * Data(s, t);
            SumWeight += W;
          }
        }
      }

      // Adding central coefficient
      TsFilt += WRMax*Data(x, y);
      SumWeight += WRMax;
      if(SumWeight > 1.0e-10) TsFilt /= SumWeight;
      else TsFilt = Data(x, y) ;

      Res.set_eigenmat_at(TsFilt, x, y);
      
      if(cnt%9000 == 0 && disp) {
#pragma omp critical
        disp->display(Res.get_colcov3()).set_title("Filtering.");
      }
    }

  }

  return Res.move_to(*this);
}
CImgList<T> get_polsar_blf_kl(int WinSiz = 11, float gammaR = 2.0, bool TRICK=false, bool FLATW=false, CImgDisplay *disp = 0) 
{
  return CImgList<T>(*this, false).polsar_blf_kl(WinSiz, gammaR, TRICK, FLATW, disp);
}

