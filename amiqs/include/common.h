#ifndef COMMON_H
#define COMMON_H

#include <vector>
#include <cmath>
#include <iostream>
#include <ctime>   // localtime
#include <sstream> // stringstream
#include <iomanip> // put_time
#include <string>  // string



#define _VERSION_ "v1.0"

#define _AMIQS_SUCCESS_ 0 /**< integer returned after successful call of a function */

const double _PI_ = 3.141592653589793238462643383279;
                
const double _KAV_      = pow(_PI_,2)/18./1.202056903159594;
const double _FD_       = 9.*1.202056903159594/pow(_PI_,2);

const double _ASYMFAC_  = 1.26*pow(10,-3);

const double _CABDIAG_  = 257./711.;
const double _CABOFF_   = 20./711.;


const double _G0_       = 0.013296432287963416; //0.0088; //0.013296432287963416; //0.0091;
const double _G1_       = 0.00673722695812018; //0.00673722695812018; //0.0051;
const double _G2_       = -0.00343289; //-0.0026; //-0.00343289;
const double _S0_       = 0.0428199999972177;
const double _S1_       = 0.008850333885126667;
const double _S2_       = -0.01636907336566425;

const double _VEV_      = 246.22/sqrt(2);
const double _MPL_      = 1.2209*pow(10,19)*sqrt(45./4./pow(_PI_,3)/106.75);
const double _TSP_      = 131.7;
const double _TC_       = 160.0;
const double _XC_       = _TSP_/_TC_;
const double _MPLTSP_   = _MPL_/_TSP_;

const double _ALPHA_W_  = 0.425/4./_PI_;

const double _YBobs_    = 8.65*pow(10,-11);
const double _YBerr_    = _YBobs_/100.;


//the following are UBUNTU/LINUX, and MacOS ONLY terminal color codes.
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */





#endif 