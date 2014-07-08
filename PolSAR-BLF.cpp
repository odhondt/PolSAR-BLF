/*
 #
 #  File        : PolSAR-BLF.cpp
 #
 #  Description : Main file for the method described in the paper 
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

#include<iostream>
#include<queue>
#include<string>
#include<math.h>
#include <Eigen/Core>
#include <Eigen/Dense> 

#define cimglist_plugin1 "SARtoolPlugin.h"
#define cimglist_plugin2 "PolSAR-BLF_plugin.h"
#include "CImg.h"


//#define EIGEN_DONT_PARALLELIZE

using namespace cimg_library;
using namespace std;

// Get basename of a file
string GetBase (const string& Str)
{
  string Tmp;
  string Res;
  size_t Slash;
  size_t Dot;
  Slash=Str.find_last_of("/");
  Dot=Str.find_last_of(".");

  if(Dot - Slash > 0)
    Tmp = Str.substr(Slash + 1, Dot - Slash -1);
  else Tmp = Str;

  Res = Tmp;

  return Res;
    
}

string GetExt (const string& Str)
{
  string Tmp;
  string Res;
  size_t Slash;
  size_t Dot;
  Slash=Str.find_last_of("/");
  Dot=Str.find_last_of(".");

  if(Dot - Slash > 0)
    Tmp = Str.substr(Dot, Str.size() - 1);
  else Tmp = Str;

  Res = Tmp;

  return Res;
    
}


string Int2String(int n)
{
  stringstream ss;
  string str;
  ss << n;
  ss >> str;
  return str;
}

string Float2String(float n)
{
  stringstream ss;
  string str;
  ss << n;
  ss >> str;
  return str;
}

int String2Int(string s)
{
  stringstream ss;
  int n;
  ss << s;
  ss >> n;
  return n;
}

int main(int argc, char **argv)
{

  cimg_usage("Iterative bilateral filtering of PolSAR data.");
  cimg_help("Options:");
  const char *file_in = cimg_option("-i", (char*)0, "Input coherency/covariance matrix file");
  const float gammaR = cimg_option("-gr", 1.33, "Radiometric scale.");
  const float gammaS = cimg_option("-gs", 2.8, "Spatial scale.");
  const int METHOD = cimg_option("-meth", 1, "Filtering method: 1 -> Affine Invariant, 2 -> Log-Euclidean, 3 -> Kullback-Leibler");
  const int NIter = cimg_option("-nit", 4, "Number of iterations.");
  const bool TRICK = cimg_option("-trick", true, "Trick for correcting the central weight (cf ipol non-local means article).");
  const bool FLATW = cimg_option("-fl", false, "Using flat weights instead of gaussian ones.");
  const bool PREVIEW = cimg_option("-p", false, "Save a preview of the filtered image (png).");
  const bool DISP = cimg_option("-d", false, "Display images (using Pauli RGB, valid for choherency matrix).");

  
  if(!file_in){
    cerr<<"Must specify one input file with -i option...\n";
    exit(0);
  }


  string FILE_IN = file_in;
  string BASE = GetBase(file_in);
  string EXT = GetExt(file_in);
 
  CImgList<float> Cov3((FILE_IN).c_str());
  CImgList<float> Cov3Filt(Cov3);

  if(Cov3(0).spectrum()!=6){
    cerr<<"Not implemented! Matrices with size != 3 cannot be processed yet.\n";
    exit(0);
  }

  CImgDisplay *disp = 0;
  if(DISP)
    disp = new CImgDisplay;

  string TMP;
  string F_NAME; 
  cout<<"\n---Starting filtering \n";

  // PolSAR-BLF with affine invariant distance
  if(METHOD==1){
    TMP = "_polsar_blf_ai_gs_" + Float2String(gammaS);
    TMP += "_gr_" + Float2String(gammaR) + "_niter_" + Int2String(NIter);
    if(FLATW) TMP += "_flatw"; 
    if(!TRICK) TMP += "_notrick";
    for(int i = 0; i < NIter; i++) {
      cout<<"Iteration "<<i+1<<'\n';
      Cov3Filt.polsar_blf_ai(gammaS, gammaR, TRICK, FLATW, disp);
    }
  }
  
  // PolSAR-BLF with log euclidean distance
  if(METHOD==2){
    TMP = "_polsar_blf_le_gs_" + Float2String(gammaS);
    TMP += "_gr_" + Float2String(gammaR) + "_niter_" + Int2String(NIter);
    if(FLATW) TMP += "_flatw"; 
    if(!TRICK) TMP += "_notrick";
    for(int i = 0; i < NIter; i++) {
      cout<<"Iteration "<<i+1<<'\n';
      Cov3Filt.polsar_blf_le(gammaS, gammaR, TRICK, FLATW, disp);
    }
  }

  // PolSAR-BLF with symmetric kullback-leibler distance
  if(METHOD==3){
    TMP = "_polsar_blf_kl_gs_" + Float2String(gammaS);
    TMP += "_gr_" + Float2String(gammaR) + "_niter_" + Int2String(NIter);
    if(FLATW) TMP += "_flatw"; 
    if(!TRICK) TMP += "_notrick";
    for(int i = 0; i < NIter; i++) {
      cout<<"Iteration "<<i+1<<'\n';
      Cov3Filt.polsar_blf_kl(gammaS, gammaR, TRICK, FLATW, disp);
    }
  }

  std::cout<<"\n--- Saving Data: "<<'\n'; 
  F_NAME = BASE + TMP + ".cimg"; 
  std::cout<<F_NAME<<'\n'; 
  Cov3Filt.save(F_NAME.c_str());

  if(PREVIEW) {
    string P_NAME;
    P_NAME = BASE + TMP + ".png";
      Cov3Filt.get_colcov3().save(P_NAME.c_str());
  }

 }
