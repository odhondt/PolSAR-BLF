/*
 #
 #  File        : cimg2rat.cpp
 #
 #  Description : See below
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


/*
 * File: cimg2rat.cpp
 *       ( C++ source file)
 *
 * Description: A file converter from RAT (version 1.0 - Display) to CImg format.
 *
 */


#include<iostream>

#define CSAR_PROGNAME "cimg2rat"
#define CSAR_USAGE "Export a \"cimg\" file format into a \"rat\" file format" 
#define CSAR_VERSION  "1.0.0"

#include <csar_soft.h>
using namespace cimg_library;
using namespace std;


template <typename T>
int cimg2rat(const char *ifname, const char *ofname, string mode){
  CImgList<T> im;
  im.load_cimg_header(ifname);
  if (mode == "-f2Dimg" || mode == "-d2Dimg" || mode == "-c2Dimg")                      im.write_rat_header(ofname, mode, im.cimg_xdim(), im.cimg_ydim());
  if (mode == "-f2Dvec" || mode == "-d2Dvec" || mode == "-c2Dvec" || mode == "-c2Dlex") im.write_rat_header(ofname, mode, im.cimg_xdim(), im.cimg_ydim(), im.cimg_cdim());
  if (mode == "-f2Dmat" || mode == "-c2Dmat" || mode == "-c2Dcov")                      im.write_rat_header(ofname, mode, im.cimg_xdim(), im.cimg_ydim(), im.get_covmat_dim(), im.get_covmat_dim());
  if (mode == "-f3Dimg")                                                                im.write_rat_header(ofname, mode, im.cimg_xdim(), im.cimg_ydim(), im.cimg_zdim());
  im.get_tiles();
  for (int k=0; k<im.get_bn(); k++) {
    im.progress(k);
    im.read_cimg(k);
//    im.read_cimg(k);
    im.write_rat(k);
  }
//  im.close_cimg_file();
//  im.close_cimg_file2();
  im.close_rat_file();
  im.progress();
  return 0;
}

int main (int argc, char **argv) {
  
  // check the command line
  cimg_usage(CSAR_USAGE);
  const char *ifname       = cimg_option("-i",      (char*)0, "[Input]  Input data                           - CIMG File\n");
  const char *ofname       = cimg_option("-o",      (char*)0, "[Output] Output data                          - RAT File\n");
  
  const char *mode_f2Dimg  = cimg_option("-f2Dimg", (char*)0, "[Param]  Convert float   2D image             - e.g. SAR Amplitude");
  const char *mode_f2Dvec  = cimg_option("-f2Dvec", (char*)0, "[Param]  Convert float   2D image of vector   - e.g. To be defined");
  const char *mode_f2Dmat  = cimg_option("-f2Dmat", (char*)0, "[Param]  Convert float   2D image of matrix   - e.g. To be defined\n");
  const char *mode_d2Dimg  = cimg_option("-d2Dimg", (char*)0, "[Param]  Convert double  2D image             - e.g. TBD");
  const char *mode_d2Dvec  = cimg_option("-d2Dvec", (char*)0, "[Param]  Convert double  2D image of vector   - e.g. TBD");
  const char *mode_c2Dimg  = cimg_option("-c2Dimg", (char*)0, "[Param]  Convert complex 2D image             - e.g. SAR slcs");
  const char *mode_c2Dvec  = cimg_option("-c2Dvec", (char*)0, "[Param]  Convert complex 2D image of vector   - e.g. SAR Scattering Vector");
  const char *mode_c2Dlex  = cimg_option("-c2Dlex", (char*)0, "[Param]  Convert complex 2D image of vector   - RAT lexicographic basis");
  const char *mode_c2Dmat  = cimg_option("-c2Dmat", (char*)0, "[Param]  Convert complex 2D image of matrix   - e.g. SAR Coherency Matrix\n");
  const char *mode_c2Dcov  = cimg_option("-c2Dcov", (char*)0, "[Param]  Convert complex 2D covariance matrix - RAT covariance matrix (lexicographix basis)");
  const char *mode_f3Dimg  = cimg_option("-f3Dimg", (char*)0, "[Param]  Convert float   3D image             - e.g. Tomograms\n");

  const int   mode_version = cimg_option("-version",       1, "[Param]  RAT format version 1 or 2            - default: 1\n");
  
  const char
  *const is_help1 = cimg_option("-h", (char*)0, 0),
  *const is_help2 = cimg_option("--h", (char*)0, 0),
  *const is_help3 = cimg_option("-help", (char*)0, 0),
  *const is_help4 = cimg_option("--help", (char*)0, 0);
  if (is_help1 || is_help2 || is_help3 || is_help4) return 0;
	
	// ----- TEST INPUT FILE -----
	if (!ifname) io_error_disp(CSAR_I_FNAME);
  if (!ofname) io_error_disp(CSAR_O_FNAME);
  disp_begin_prog(string(ifname), string(ofname));
  
  if (mode_f2Dimg) return cimg2rat<float>(ifname, ofname, string(mode_f2Dimg));
  if (mode_f2Dvec) return cimg2rat<float>(ifname, ofname, string(mode_f2Dvec));
  if (mode_f2Dmat) return cimg2rat<float>(ifname, ofname, string(mode_f2Dmat));
  if (mode_d2Dimg) return cimg2rat<double>(ifname, ofname, string(mode_d2Dimg));
  if (mode_d2Dvec) return cimg2rat<double>(ifname, ofname, string(mode_d2Dvec));
  if (mode_c2Dimg) return cimg2rat<float>(ifname, ofname, string(mode_c2Dimg));
  if (mode_c2Dvec) return cimg2rat<float>(ifname, ofname, string(mode_c2Dvec));
  if (mode_c2Dlex) return cimg2rat<float>(ifname, ofname, string(mode_c2Dlex));
  if (mode_c2Dmat) return cimg2rat<float>(ifname, ofname, string(mode_c2Dmat));
  if (mode_c2Dcov) return cimg2rat<float>(ifname, ofname, string(mode_c2Dcov));
  if (mode_f3Dimg) return cimg2rat<float>(ifname, ofname, string(mode_f3Dimg));
  
	
    
  
  
  //const char *mode_slc    = cimg_option("-slc",    (char*)0, "Convert a SAR complex image");
  //const char *mode_amp    = cimg_option("-amp",    (char*)0, "Convert a SAR amplitude image");
  //const char *mode_pha    = cimg_option("-pha",    (char*)0, "Convert a SAR phase image");
  //const char *mode_int    = cimg_option("-int",    (char*)0, "Convert a SAR intensity image");
  //const char *mode_lex    = cimg_option("-lex",    (char*)0, "Convert a polarimetric scattering vector, lexicographic basis");
  //const char *mode_pauli  = cimg_option("-pauli",  (char*)0, "Convert a polarimetric scattering vector, Pauli basis");
  //const char *mode_vec    = cimg_option("-vec",    (char*)0, "Convert a general scattering vector");
  //const char *mode_cov    = cimg_option("-cov",    (char*)0, "Convert a polarimetric covariance matrix [C]");
  //const char *mode_coh    = cimg_option("-coh",    (char*)0, "Convert a polarimetric coherency matrix [T]");
  //const char *mode_covmat = cimg_option("-covmat", (char*)0, "Convert an arbitrary covariance matrix NxN");
  //const char *mode_floatvec  = cimg_option("-floatvec",  (char*)0, "[Option]  Convert 2D float array as vector");
  //const char *mode_doubleimg = cimg_option("-doubleimg", (char*)0, "Convert double 2D array");
  //const char *mode_doublevec = cimg_option("-doublevec", (char*)0, "Convert double 2D array as vector");
  
  
  
  /*int mode = -1;
  int idl_type = 0;
	if (mode_slc)    { mode = 101; idl_type = 6;}
  if (mode_amp)    { mode = 101; idl_type = 4;}
  if (mode_pauli)  { mode = 210; idl_type = 6;}
  if (mode_vec)    { mode = 210; idl_type = 6;}
  if (mode_coh)    { mode = 221; idl_type = 6;}
  if (mode_covmat) { mode = 221; idl_type = 6;}
  if (mode_floatvec)  {cimg2rat<float>(ifname, ofname, 210, 4);  return 0;}
  if (mode_doubleimg) {cimg2rat<double>(ifname, ofname, 101, 5); return 0;}
  if (mode_doublevec) {cimg2rat<double>(ifname, ofname, 210, 5); return 0;}
	if (mode != -1) cimg2rat<float>(ifname, ofname, mode, idl_type);*/
  std::cerr <<"Not yet implemented: see cimg2rat -h for more information!\n" << std::endl;
  return 0;
}
