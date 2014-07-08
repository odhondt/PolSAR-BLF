/*
 #
 #  File        : csar_soft.h
 #
 #  Description : See below
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


/*
 * File: import_rat.cpp
 *       ( C++ source file)
 *
 * Description: Import RAT file.
 *

*/


#define CSAR_PROGNAME "import_rat"
#define CSAR_USAGE    "Import RAT file"
#define CSAR_VERSION  "1.0.0"

#include <csar_soft.h>
using namespace cimg_library;
using namespace std;

template <typename T>
int import_rat(const char* ifname, const char* ofname, bool mode_lex=false, bool mode_3D=false){
  CImgList<T> im;
  im.load_rat_header(ifname, mode_lex, mode_3D);
  // simple case of sar image
  if (im.rat_ndim() == 2) im.write_cimg_header(ofname, im.rat_cplx(), im.rat_xdim(), im.rat_ydim());
  if (im.rat_ndim() == 3){
    if (!mode_3D) im.write_cimg_header(ofname, im.rat_cplx(), im.rat_xdim(), im.rat_ydim(), 1, im.rat_zdim());
    else          im.write_cimg_header(ofname, im.rat_cplx(), im.rat_xdim(), im.rat_ydim(), im.rat_zdim());
  }
  if (im.rat_ndim() == 4) {
    if (im.rat_zdim() != im.rat_vdim()) {
      cerr << "Not yet implemented! Exit now!" << endl;
      im.close_rat_file();
      return 0;
    }
    im.write_cimg_header(ofname, im.rat_cplx(), im.rat_xdim(), im.rat_ydim(), 1, (im.rat_zdim() * (im.rat_zdim()+1))/2);
  }
  im.get_tiles();
  for (int k=0; k<im.get_bn(); k++) {
    im.progress(k);
    im.read_rat(k);
    im.write_cimg(k);
  }
  im.close_rat_file();
  im.close_cimg_file();
  im.progress();
  return 0;
}

int main(int argc, char **argv){

  // check the command line
  cimg_help("\nUsage: import_rat -i input_data.rat -o output_data.cimg [-lex] [-3D] \n");
  const char* ifname = cimg_option("-i",     (char*)0, "[Input]  RAT format data file name    - *.rat  file\n");
  const char* ofname = cimg_option("-o",     (char*)0, "[Output] CSAR format data file name   - *.cimg file\n");
  cimg_help("Options:");
  const char* mlex   = cimg_option("-lex",   (char*)0, "[Option] import a vector/matrix using lexicographic basis (special coding in RAT)");
  const char* m3D    = cimg_option("-3D",    (char*)0, "[Option] import a 3D image (mainly tomogram)\n");
  
  cimg_help("\"import_rat\" will try to guess the kind of data. If generating using RAT, it will contain a code that will be used.\n");
  
  const char
  *const is_help1 = cimg_option("-h",     (char*)0, 0),
  *const is_help2 = cimg_option("--h",    (char*)0, 0),
  *const is_help3 = cimg_option("-help",  (char*)0, 0),
  *const is_help4 = cimg_option("--help", (char*)0, 0);
  if (is_help1 || is_help2 || is_help3 || is_help4) return 0;
  
  disp_begin_prog();
  
  //----------------------------------------------------------------------------------------------------
  // ---> Test input information (command line)
  //----------------------------------------------------------------------------------------------------
  if (!ifname) io_error_disp(CSAR_I_FNAME);
  if (!ofname) io_error_disp(CSAR_O_FNAME);
  CImgList<float> in;
  in.load_rat_header(ifname);
  in.close_rat_file();
  bool mode_lex = false, mode_3D = false;
  if (mlex) mode_lex = true;
  if (m3D)  mode_3D  = true;
  if (in.rat_type() == 2) return import_rat<short> (ifname, ofname, mode_lex, mode_3D);
  if (in.rat_type() == 3) return import_rat<long>  (ifname, ofname, mode_lex, mode_3D);
  if (in.rat_type() == 4) return import_rat<float> (ifname, ofname, mode_lex, mode_3D);
  if (in.rat_type() == 5) return import_rat<double>(ifname, ofname, mode_lex, mode_3D);
  if (in.rat_type() == 6) return import_rat<float> (ifname, ofname, mode_lex, mode_3D);
  if (in.rat_type() == 9) return import_rat<double>(ifname, ofname, mode_lex, mode_3D);
  
        
  
  
  std::cerr <<"Not yet implemented: see import_rat -h for more information!\n" << std::endl;
  return 0;
  
}
