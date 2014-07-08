/*
 #
 #  File        : csar_soft.h
 #  
 #  Description : A stupid file totally obsolete but necessary.
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

#ifndef CSAR_PROGNAME
#define CSAR_PROGNAME "routines_name"
#endif
#ifndef CSAR_USAGE
#define CSAR_USAGE "define an usage"
#endif
#ifndef CSAR_VERSION
#define CSAR_VERSION  "x.x.x"
#endif

#define CSAR_I_FNAME 0
#define CSAR_O_FNAME 1

#ifndef _csar_soft_h
#define _csar_soft_h

typedef unsigned short uint16_t;


#include<fstream>
#include<iostream>
#include<complex>
#include<cstdarg>
#include<cstring>
#include<ctime>
#include<new>
#include<string>
#include<vector>
#include<cmath>

#define cimg_display 0
#define cimglist_plugin "csar.h"
#include <CImg.h>

typedef struct {int p;} int_one_param_struct;

//----------------------------------------------------------------------------------------------------
// DISPLAY I/O ERROR MESSAGE AND EXIT
//----------------------------------------------------------------------------------------------------
#pragma mark -
#pragma mark Display i/o error message.
//! Display standard i/o error message and exit
void io_error_disp(int mode){
  if (mode == 0)
    std::cerr << "["<< CSAR_PROGNAME << "] runtime error! No input filename set (type \'" << CSAR_PROGNAME << " -h\' to get help)." << std::endl << std::endl;
  if (mode == 1)    
    std::cerr << "["<< CSAR_PROGNAME << "] runtime error! No output filename set (type \'" << CSAR_PROGNAME << " -h\' to get help)." << std::endl << std::endl;
  exit(1);
}

//! Display a specific message and exit.
void io_error_disp(const char* mess){
  std::cerr << "[" << CSAR_PROGNAME << "] runtime error! " << mess << " (type \'" << CSAR_PROGNAME << " -h\' to get help)." << std::endl << std::endl;
  exit(1);
}

void io_miss_disp(const char* mess){
  std::cerr << "[" << CSAR_PROGNAME << "] runtime error! " << mess << " file name is missing! (type \'" << CSAR_PROGNAME << " -h\' to get help)." << std::endl << std::endl;
  exit(1);
}

//----------------------------------------------------------------------------------------------------
// DISPLAY BEGIN OF THE PROGRAM
//----------------------------------------------------------------------------------------------------
#pragma mark -
#pragma mark Display begin of a program.
void disp_begin_prog(){
  char buffer[512];
  int n = std::sprintf(buffer, "*****       %s - %s - v%s       *****", CSAR_PROGNAME, CSAR_USAGE, CSAR_VERSION);
  std::cout << std::endl;
  for (int c=0; c<(n); c++) std::cout << "-"; std::cout << std::endl;
  std::cout << std::endl;
  std::cout << buffer << std::endl;
  std::cout << std::endl;
  for (int c=0; c<(n); c++) std::cout << "-"; std::cout << "\n" << std::endl;
}

void disp_begin_prog(const char *ifname, const char *ofname){
  char buffer[512];
  int n = std::sprintf(buffer, "*****       %s [%s] - v%s       *****", CSAR_USAGE, CSAR_PROGNAME, CSAR_VERSION);
  std::cout << std::endl;
  for (int c=0; c<(n); c++) std::cout << "-"; std::cout << std::endl;
  std::cout << std::endl;
  std::cout << buffer << std::endl;
  std::cout << std::endl;
  for (int c=0; c<(n); c++) std::cout << "-"; std::cout << "\n" << std::endl;
	std::cout << "  INPUT  File:    " << ifname  << std::endl;
	std::cout << "  OUTPUT File:    " << ofname << "\n" << std::endl;  
}

void disp_begin_prog(std::string ifname, std::string ofname){
  char buffer[512];
  int n = std::sprintf(buffer, "*****       %s [%s] - v%s       *****", CSAR_USAGE, CSAR_PROGNAME, CSAR_VERSION);
  std::cout << std::endl;
  for (int c=0; c<(n); c++) std::cout << "-"; std::cout << std::endl;
  std::cout << std::endl;
  std::cout << buffer << std::endl;
  std::cout << std::endl;
  for (int c=0; c<(n); c++) std::cout << "-"; std::cout << "\n" << std::endl;
	std::cout << "  INPUT  File:    " << ifname  << std::endl;
	std::cout << "  OUTPUT File:    " << ofname << "\n" << std::endl;  
}

//----------------------------------------------------------------------------------------------------
// READ XML FILE AND TRANSFERT IT INTO A BUFFER (USEFUL FOR RAPIDXML)
//----------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------
// CONVERT A VARIABLE INTO A STRING
//----------------------------------------------------------------------------------------------------
#pragma mark -
#pragma mark convert to type into string 

// see:http://cpp.developpez.com/faq/cpp/?page=strings#STRINGS_convertform for more details (in french)
template<typename T>
std::string to_string( const T &input){
  std::ostringstream oss;
  oss << input;
  return oss.str();
}




#endif
