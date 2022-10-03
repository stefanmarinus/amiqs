#include <fstream>
#include <complex>
#include <chrono>
#include <ctime>  
#include <thread>
#include <gsl/gsl_spline2d.h>
#include <stdio.h>
#include <string.h>

#include "include/LSODA.h"
#include "include/common.h"

#include "Source/Yukawas.cpp"
#include "Source/rates.cpp"
#include "Source/MatMul.cpp"
#include "Source/kinetic_eq.cpp"
#include "Source/Ana_Approx.cpp"
#include "Source/read_ini.cpp"


using namespace std;
typedef complex<double> complx;




/*
define the Number of rows
*/
#define N 3
#define M 2


/*
Function to find maximal matrix element
mat[][] : 2D array to find max element
*/
double findMax(complx mat[N][M]){
     /*
     Initializing max element
     */
     double maxElement = -99999999999.9;

     /*
     checking each element of matrix
     if it is greater than maxElement,
     update maxElement
     */
     for (int i = 0; i < N; i++) {
          for (int j = 0; j < M; j++) {
               if (abs(mat[i][j]) > maxElement) {
                    maxElement = abs(mat[i][j]);
               }
          }
     }
  return maxElement;
}


/*
Function to find minimal matrix element
mat[][] : 2D array to find max element
*/
double findMin(complx mat[N][M]){
     /*
     Initializing min element
     */
     double minElement = 99999999999999.9;

     /*
     checking each element of matrix
     if it is smalles than minElement,
     update maxElement
     */
     for (int i = 0; i < N; i++) {
          for (int j = 0; j < M; j++) {
               if (abs(mat[i][j]) < minElement) {
                    minElement = abs(mat[i][j]);
               }
          }
     }
     return minElement;
}


/*
check if input argument is int or not
*/
bool is_number(const string& s){
    return !s.empty() && find_if(s.begin(), 
        s.end(), [](unsigned char c) { return !isdigit(c); }) == s.end();
}


/*
Function translating the inputs of Lambda & mu to
physical parameters M1 & DM/M
*/
double Lambda_to_M1(double Lambda, double mu){
     return log10(pow(10,Lambda) - pow(10,mu));
}
double mu_to_DM_M(double Lambda, double mu){
     return log10(2.*pow(10,mu)/(pow(10,Lambda) - pow(10,mu)));
}


/*log scale input but return in linear scale*/
double D12(double M1, double DM_M){
     // vector<double> param_var      = *(vector<double> *)param_ini;
     // double M1                     = pow(10,param_var[0]);
     // double DM_M                   = pow(10,param_var[1]);
     double M2                     = pow(10,M1) * (pow(10,DM_M) + 1.);

     return pow(M2,2) - pow(pow(10,M1),2);
}


/*
Calculate some transition times which define different washout regimes
*/
vector<double> transition_times(double parametrization, void* param_ini){
     /*
     declare variables used
     */
     vector<double> param_var;
     double M1, DM_M;
     double d12;
     double y_sum;



     param_var = *(vector<double> *)param_ini;
     if (parametrization == 0){
          M1= param_var[1]; DM_M= param_var[2];
          Yukawa.YukawaLN(param_ini, YLN);
     }
     else{
          M1= Lambda_to_M1(param_var[0], param_var[1]); DM_M= mu_to_DM_M(param_var[0], param_var[1]);
          Yukawa.Yukawa_anaLN(param_ini, YLN);
     }


     d12       = D12(M1, DM_M);


     y_sum     = sqrt( pow(abs(YLN[0][0]),2) + pow(abs(YLN[1][0]),2) + pow(abs(YLN[2][0]),2) );


     if (y_sum > 0){
          // double x_transition = pow((16.*pow(rates.lnc_0(pow(10,M1),1.),2) + pow(_KAV_,2))*(5*pow(_TSP_,5.)*(9.*rates.lnc_0(pow(10,M1),1.)+rates.lnc_1(pow(10,M1),1.)*_FD_)*pow(y_sum,2))
          //                     /(4.*_MPL_*pow(d12,2)*rates.lnc_0(pow(10,M1),1.)*(9.*rates.lnc_0(pow(10,M1),1.)+2.*rates.lnc_1(pow(10,M1),1.)*_FD_)*pow(_KAV_,2)),1./5.);
          // double x_0          = x_transition;
          double x_0          = pow(5.*pow(y_sum,2.)* (pow(_G0_*_MPLTSP_,2.) + 4.*pow(_KAV_/8.*_MPLTSP_,2.) )/(_G0_*_MPLTSP_*pow(_KAV_/2.*d12/pow(_TSP_,2.)*_MPLTSP_,2.)),1./5.);
          double x_osc        =  _PI_/3.*_G0_*pow(y_sum,2)*pow(2/_KAV_,1./3.)*pow(_MPL_*_MPL_/d12,1./3.);
          double s0           = _S0_;
          double s1           = _S1_;
          double x_lnv        = pow(1./3. * pow(10,M1)*pow(10,M1) * s0 * _MPL_/pow(_TSP_,3.) * pow(y_sum,2.), -1./3.); //pow(1./3. * pow(10,M1)*pow(10,M1) * rates.lnv_0(pow(10,M1), 0.001) * _MPL_/pow(_TSP_,3.) * pow(y_sum,2.) ,-1./3.);
          double x_lnv_int    = pow(2./3.* ( _G1_ * s0 + _G0_*s1*_FD_*pow(10,M1)*pow(10,M1)*_MPL_/pow(_TSP_,3.) *pow(y_sum,2.)/(9.*_G0_ + 2.*_FD_*_G1_)) ,-1./3.);
          double ww           = _MPL_/_TSP_ * (2. * (3. *_G0_ * _G1_ * _FD_ + pow(_G1_*_FD_,2)))/(3. * (3. * _G0_ + 2. * _G1_ * _FD_))*1/2.*pow(y_sum,2);

          double sw_e         = _MPLTSP_*1./2.*_FD_*_G1_*pow(abs(YLN[0][0]),2);
          double sw_mu        = _MPLTSP_*1./2.*_FD_*_G1_*pow(abs(YLN[1][0]),2);
          double sw_tau       = _MPLTSP_*1./2.*_FD_*_G1_*pow(abs(YLN[2][0]),2);

          return {x_0, x_osc, x_lnv, x_lnv_int, ww, sw_e, sw_mu, sw_tau};
     }
     else{
          return {0, 0, 0, 0, 0, 0, 0, 0};
     }
}


/*
function which prints the information of the transition times and the corresponding washout regime
*/
void print_transition_times(double parametrization, double lnv_rates, void* param_ini){
     /*
     declare variables
     */
     vector<double> x_vec;
     double x_0, x_osc, x_lnv, x_lnv_int, ww;


     x_vec = transition_times(parametrization, param_ini);
     x_0= x_vec[0]; x_osc= x_vec[1]; x_lnv= x_vec[2]; x_lnv_int= x_vec[3]; ww= x_vec[4];


     if (x_0>0 && x_osc>0 && x_lnv>0 && x_lnv_int>0 && ww>0){

          cout << " [AMIQS] x_0:                       " << x_0 << endl;
          cout << " [AMIQS] x_osc:                     " << x_osc << endl;
          cout << " [AMIQS] x_M^ov:                    " << x_lnv << endl;
          cout << " [AMIQS] x_M^int:                   " << x_lnv_int << endl;



          if(lnv_rates == 1){
               if(x_0>1){
                    if(x_lnv >1){
                         if (ww < 1){
                         cout << "[AMIQS] >> WEAK LNV + OVERDAMPED WEAK WASHOUT <<   " << endl;
                         }
                         else{
                         cout << "[AMIQS] >> WEAK LNV + OVERDAMPED <<                 " << endl;
                         }
                    }
                    else{
                         if(ww <1){
                         cout << "[AMIQS] >> STRONG LNV + OVERDAMPED WEAK WASHOUT <<  " << endl;
                         }
                         else{
                         cout << "[AMIQS] >> STRONG LNV + OVERDAMPED <<               " << endl;
                         }
                    }
               }
               else{
                    if(x_lnv_int >1){
                         if (ww <1){
                         cout << "[AMIQS] >> WEAK LNV + OSCILLATION WEAK WASHOUT <<  " << endl;
                         }
                         else if( x_osc <1){
                         cout << "[AMIQS] >> WEAK LNV + FAST OSCILLATIONS <<          " << endl;
                         }
                         else{
                         cout << "[AMIQS] >> WEAK LNV + INTERMEDIATE <<               " << endl;
                         }
                    }
                    else{
                         if (ww < 1){
                         cout << "[AMIQS] >> STRONG LNV + OSCILLATION WEAK WASHOUT <<               " << endl;
                         }
                         else if(x_osc < 1){
                         cout << "[AMIQS] >> STRONG LNV + FAST OSCILLATIONS <<               " << endl;
                         }
                         else{
                         cout << "[AMIQS] >> STRONG LNV + INTERMEDIATE <<               " << endl;
                         }
                    }
               }
          }
          else{
               if(x_0 > 1){
                    if (ww<1){
                         cout << "[AMIQS] >> PURE LNC OVERDAMPED WEAK WASHOUT <<               " << endl;
                    }
                    else{
                         cout << "[AMIQS] >> PURE LNC OVERDAMPED <<                " << endl;
                    }
               }
               else{
                    if (ww <1){
                         cout << "[AMIQS] >> PURE LNC OSCILLATION WEAK WASHOUT <<               " << endl;
                    }
                    else if( x_osc <1){
                         cout << "[AMIQS] >> PURE LNC FAST OSCILLATIONS <<              " << endl;
                    }
                    else{
                         cout << "[AMIQS] >> PURE LNC INTERMEDIATE <<            " << endl;
                    }
               }
          }
     }
}


/*
these two functions define the x_value at which the oscilations reach
the frequency of 10^5 or pass 1000 complete oscillations
This will be the point at which the kinetic equations will be averaged -> x_ave
*/
double x_osc_num(double num, double d12){
     return  pow(num*2.*_PI_/_KAV_*6.*pow(_TSP_,3.)/_MPL_/d12,1./3.);
}
double x_osc_freq(double freq, double d12){
     double factor = 2.*_PI_/_KAV_*6.*pow(_TSP_,3)/_MPL_/d12;
     if (factor > 1/pow(freq,3.)){
          double num_to_reach_freq = 1./18.*(sqrt(3)*sqrt((pow(pow(freq,3.)*factor+2.,2)*(4.*pow(freq,3)*factor - 1.))/(pow(freq,6.)*pow(factor,2.))) + 9.);
          // cout << "[AMIQS] # osc. to reach max freq.   " << num_to_reach_freq << endl;
          return x_osc_num(num_to_reach_freq, d12);
     }
     else{
          return 1;
     }
}


/*
returns integer number corresponding to a concrete washout regime
0 -> Overdamped + wLNV
1 -> Overdamped + sLNV
2 -> Intermediate 
3 -> Fast oscillations
*/
bool regime_check(int regime, double parametrization, void* param_ini){
     vector<double> var;
     var = transition_times(parametrization, param_ini);

     /*
     overdamped wLNV
     */
     if (regime == 0){
          if (var[0] > 5 && var[2] > 10 && var[5] > 5 && var[6] > 5 && var[7] > 5){
               return true;
          }
          else{
               return false;
          }
     }
     /*
     overdamped sLNV
     */
     else if (regime == 1){
          if (var[0] > 5 && var[2] < 0.5 && var[5] > 5 && var[6] > 5 && var[7] > 5){
               return true;
          }
          else{
               return false;
          }
     }
     /*
     Intermediate
     */
     else if (regime == 2){
          if (var[0] < 0.01 && var[1] > 10 && var[5] < 0.2 && var[6] > 5 && var[7] > 5){
               return true;
          }
          else{
               return false;
          }
     }
     /*
     Fast oscillations
     */
     else if (regime == 3){
          if (var[0] < 0.01 && var[1] < 0.1 && var[5] < 0.2 && var[6] > 5 && var[7] > 5){
               return true;
          }
          else{
               return false;
          }
     }
     else{
          return false;
     }
}


/*
function joining all the stuff to call and solve the kinetic equaiton 
*/
double amiqs(double parametrization, void* setting_ini, void* param_ini, double info, double safe){
     vector<double> param_var           = *(vector<double> *)param_ini;
     vector<double> setting_var         = *(vector<double> *)setting_ini;
     double hierarchy, M1, DM_M, yukawa, theta, delta, phi;
     double zre, zim;
     double Lambda, mu, ye1, ymu1, ytau1, ye2, ymu2, ytau2, phie2, phimu2, phitau2;
     double lnv_rates, nl_approx, rates_approx, c_matrix_approx, g1_approx, sph_approx;
     double asym_inst, asym_smooth;


     lnv_rates= setting_var[0]; nl_approx= setting_var[1]; rates_approx= setting_var[2]; c_matrix_approx= setting_var[3]; g1_approx= setting_var[4], sph_approx= setting_var[5];

     if (parametrization == 0){
          hierarchy= param_var[0]; M1 = param_var[1]; DM_M = param_var[2]; yukawa = param_var[3]; theta = param_var[4]; delta = param_var[5]; phi = param_var[6];
     }
     else if (parametrization == 1){
          Lambda= param_var[0]; mu= param_var[1]; 
          ye1= param_var[2]; ymu1= param_var[3]; ytau1= param_var[4]; 
          ye2= param_var[5]; ymu2= param_var[6]; ytau2= param_var[7]; 
          phie2= param_var[8]; phimu2= param_var[9]; phitau2= param_var[10];

          M1        = Lambda_to_M1(Lambda, mu);
          DM_M      = mu_to_DM_M(Lambda, mu);
     }
     else if (parametrization == 2){
          hierarchy= param_var[0]; M1 = param_var[1]; DM_M = param_var[2]; zre = param_var[3]; zim = param_var[4]; delta = param_var[5]; phi = param_var[6]; 
     }

     /*
     Number of equations and time span to solve
     */
     int neq = 12;  double tmin = 0.00002; double t = tmin; double tout;

     if (safe == 0){
          tout = 1.;
     }
     else{
          tout = 1.001*tmin;
     }

     /*
     initial condition of the system
     */
     vector<double> y    = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
     /*
     set solver state to true
     */
     int istate          = 1;
     /*
     initialize LSODA solver
     */
     LSODA lsoda;
     /*
     initialize aux result vector
     */
     vector<double> myres;
     /*
     initialize the output vector
     */
     vector<double> yout;
     /*
     set averaging scales
     */
     double osc_num_max  = pow(10,3);
     double osc_freq_max = pow(10,4);

     /*
     print some info
     */
     if (info == 1){
          if (parametrization == 0 || parametrization == 2){
               if (hierarchy == 0){
                    cout << "[AMIQS] Normal Hierarchy                 " << endl;
               }
               else if (hierarchy == 1){
                    cout << "[AMIQS] Inverted hierarchy               " << endl;
               }
               else{
                    cout << RED << "[AMIQS]" << RESET << " Hierarchy not correctly read. Choose between 0 (NH) and 1 (IH)" << endl;
                    exit(0);
               }
          }
          if (lnv_rates == 0){
               cout << "[AMIQS] LN Conserving                    " << endl;
          }
          else if (lnv_rates == 1){
               cout << "[AMIQS] LN Conserving+Violating          " << endl;
          }
          else{
               cout << RED << "[AMIQS]" << RESET << " Interactions not correctly read. Choose between 0 (HC) and 1 (HC+HV)" << endl;
               exit(0);
          }
     }


     /*
     get the squared mass difference for given input
     */
     double d12          = D12(M1, DM_M);

     /*
     get the Yukawa matrices for given input and check if LN constraint is fulfilled
     */
     if(parametrization == 0){
          Yukawa.YukawaLN(param_ini, YLN);
          Yukawa.Yukawa(YLN, Y);

          // YLN[0][0] = 1.99746340* pow(10,-7) * polar(1.,1.07082355); YLN[1][0] = 2.08930761 * pow(10,-6.) * polar(1.,-0.21441471); YLN[2][0] = 1.86841946* pow(10,-6.) * polar(1.,-0.8711756);
          // YLN[0][1] = 1.99751319*pow(10,-7.) * polar(1.,2.64156543); YLN[1][1] = 2.08931149 * pow(10,-6.) * polar(1.,1.356393); YLN[2][1] = 1.86839934 * pow(10,-6.) * polar(1.,0.69962497);

          // Y[0][0] = YLN[0][0]; Y[1][0] = YLN[1][0]; Y[2][0] = YLN[2][0];
          // Y[0][1] = YLN[0][1]; Y[1][1] = YLN[1][1]; Y[2][1] = YLN[2][1];

     }
     else if (parametrization == 1){
          Yukawa.Yukawa_anaLN(param_ini, YLN);
          Yukawa.Yukawa_ana(YLN, Y);
     }
     else{
          Yukawa.Yukawa_CILN(param_ini, YLN);
          Yukawa.Yukawa_CI(YLN, Y);
     }
     // Yukawa.Yukawa(YLN, Y); 
     Yukawa.YukawaHC(Y,YHC); Yukawa.YukawaT(Y,YT); Yukawa.YukawaC(Y,YC);


     /*
     print yukawa matrices
     no fa falta pro bo...
     */
     if (info == 1){
          if (parametrization < 2){
               cout << "" << endl;
               for(int i = 0; i<3; i++){
                    for(int j = 0; j<2; j++){
                         cout << "YLN[" << i << "][" << j << "] = " << YLN[i][j] << endl;
                    }
               }
          }
          cout << "" << endl;
          for(int i = 0; i<3; i++){
               for(int j = 0; j<2; j++){
                    cout << "Y[" << i << "][" << j << "] = " << Y[i][j] << endl;
               }
          }
     cout << "" << endl;
     }
     


     double y_sum = sqrt( pow(abs(YLN[0][0]),2) + pow(abs(YLN[1][0]),2) + pow(abs(YLN[2][0]),2) );
     if (y_sum == 0){
          Yukawa.Yukawa_reset(YLN);
          Yukawa.Yukawa_reset(Y);
          y.clear();
          yout.clear();
          cout << "[AMIQS] Asymmetry set to            " << 0 << endl;
          return _AMIQS_SUCCESS_;
     }

    /*
    Check for complete thermalization of the system -> Asymmetry will be 0
    */
    double ylnmin       = findMin(YLN);
    double ratio        = _KAV_*d12/pow(_TSP_,2) / pow(y_sum,2) / rates.lnc_0(pow(10,M1), 0.001);
    if (ylnmin > 5*pow(10,-6)){
          cout << "[AMIQS] All Yukawas thermalize      " << endl;
          cout << "[AMIQS] Asymmetry set to            " << 0 << endl;
          Yukawa.Yukawa_reset(YLN);
          Yukawa.Yukawa_reset(Y);
          y.clear();
          yout.clear();
          return _AMIQS_SUCCESS_;
    }

    /*
    check for unflavoured strong washout outside the overdamped regime 
    Only valid for helicity conserving interactions
    */
    if (lnv_rates == 0){
          if (1./2.*_MPLTSP_*rates.lnc_0(pow(10,M1),0.001)*pow(min(abs(YLN[0][0]), min(abs(YLN[1][0]), abs(YLN[2][0]))),2) > 20  && ratio > 10){
               cout << "[AMIQS] SW outside ov. regime       " << endl;
               cout << "[AMIQS] Asymmetry set to            " << 0 << endl;
               /*
               clear vectors and return final asymmetry
               */
               Yukawa.Yukawa_reset(YLN);
               Yukawa.Yukawa_reset(Y);
               y.clear();
               yout.clear();
               return _AMIQS_SUCCESS_;
          }
    }
    /*
    in the case of including helcity violating interactions check if both rates thermalize
    */
     if (lnv_rates == 1){
          if (1./2.*_MPLTSP_*rates.lnc_0(pow(10,M1),0.001)*pow(min(abs(YLN[0][0]), min(abs(YLN[1][0]), abs(YLN[2][0]))),2) > 20 &&
              1./2.*_MPLTSP_*rates.lnv_0(pow(10,M1),0.001)*pow(min(abs(YLN[0][0]), min(abs(YLN[1][0]), abs(YLN[2][0]))),2) > 20){
               cout << "[AMIQS] HC & HV rates thermalize    " << endl;
               cout << "[AMIQS] Asymmetry set to            " << 0 << endl;
               /*
               clear vectors and return final asymmetry
               */
               Yukawa.Yukawa_reset(YLN);
               Yukawa.Yukawa_reset(Y);
               y.clear();
               yout.clear();
               return _AMIQS_SUCCESS_;
        }
   }

     /* 
     check up to which x we need to evolve the full system before averaging the oscillations 
     x_transition is time at which the system changes from the overdamped to fast oscilltion,
     it is derived by looking at the small (d12) eigenvalue within the overdamped regime. When
     this eigenvalue becomes comparable to the other ones we leave the overdamp regime.
     Then x_ave has to be added to this, which defines the timespan for which the oscillations
     reach a frequency of 10^5 or 1000 complete oscillations are done.
     */
     double x_transition = pow((16.*pow(rates.lnc_0(pow(10,M1),1.),2) + pow(_KAV_,2))*(5*pow(_TSP_,5.)*(9.*rates.lnc_0(pow(10,M1),1.)+rates.lnc_1(pow(10,M1),1.)*_FD_)*pow(y_sum,2))
                         /(4.*_MPL_*pow(d12,2)*rates.lnc_0(pow(10,M1),1.)*(9.*rates.lnc_0(pow(10,M1),1.)+2.*rates.lnc_1(pow(10,M1),1.)*_FD_)*pow(_KAV_,2)),1./5.);
     double x_ave        =  x_transition + min(x_osc_freq(osc_freq_max,d12), x_osc_num(osc_num_max,d12));
     if (info == 1){
          cout << CYAN <<  "[AMIQS]" << RESET << " Averaging x                 " << x_ave << endl;
     }

     /*
     call this function to print washout regime information
     if you dont need this info just comment the function call
     */
     if (info == 1){
          if (parametrization != 2){
               print_transition_times(parametrization, lnv_rates, param_ini);
          }
     }

     /*
     store doubles in a pointer which can be passed to the solver
     */
     vector<double> ini  = {pow(10,M1), d12, x_ave, lnv_rates, nl_approx, rates_approx, c_matrix_approx, g1_approx};
     void *ptr           = &ini;


     /****************************************************************************************************
                                             LSODA Integrator
     ****************************************************************************************************/
     double abs_tol      = 1e-13;
     double rel_tol      = 1e-12;
     while(tout < 1){
          lsoda.lsoda_update(fun, neq, y, yout, &t, tout, &istate, ptr, rel_tol, abs_tol);
          myres.push_back(t);
          myres.push_back(_ASYMFAC_ * (yout[9] + yout[10] + yout[11]));
          for(int i = 0; i< 12;i++ ){
               y[i] = yout[i + 1];
               /*
               ese if(){} pues no hace falta pero lo he puesto por tener los plots mas bonitos.. 
               */
               if(t > 1.1*x_ave){
                    y[2] = 0; y[3] = 0; y[6] = 0; y[7] = 0;
               }
               myres.push_back(y[i]);
          }
          if (info == 1){
               cout << " x = " << " " << t << "\r" << flush;
          }
          tout = tout * 1.01;
     }
     lsoda.lsoda_update(fun, neq, y, yout, &t, 1, &istate, ptr, rel_tol, abs_tol);
     cout << "[AMIQS] Integration finished" << "\r" << flush;



     
     if (safe == 1){
          myres.push_back(t);
          myres.push_back(_ASYMFAC_ * (yout[9] + yout[10] + yout[11]));
          for(int i = 0; i< 12;i++ ){
               y[i] = yout[i + 1];
               myres.push_back(y[i]);
          }
          int dim_vector = 14;

          int row_size   = int(myres.size()/dim_vector);

          double myres_matrix [row_size][dim_vector];


          for (int i = 0; i < row_size; ++i) {
               for (int j=0; j < dim_vector; j++){
                    myres_matrix[i][j] =  myres[dim_vector*i+j];
               }
          }
          /*
          write to file the calculated asymmetry for given input parameters
          */
          ofstream file;
          if(lnv_rates == 0){
               file.open("../output/asym_results_lnc.dat");
          }
          else{
               file.open("../output/asym_results_lnv.dat");  
          }
          
          char hierarchy_string [10];
          char lnv_rates_string [10];
          if (parametrization == 0){
               if(hierarchy == 0){
                    strcpy (hierarchy_string,"NH");
               }
               else{
                    strcpy (hierarchy_string,"IH");
               }
          }
          if(lnv_rates == 0){
               strcpy (lnv_rates_string,"LNC");
          }
          else{
               strcpy (lnv_rates_string,"LNC+LNV");
          }

          if (parametrization == 0){
               file << "# " << "Hierarchy:" << " " << hierarchy_string << " " << "Interactions:" << " " <<  lnv_rates_string << " " <<
               "M1:" << " " << M1 << " " << "DM/M:" << " " << DM_M << " " << "y:" << " " << yukawa << " " << "theta:" << " " << theta << " " << 
               "delta:" << " " << delta << " " << "phi:" << " " << phi << endl;
          }
          else if (parametrization == 1){
          file << "# " << "Interactions:" << " " <<  lnv_rates_string << " " <<
          "Lambda:" << " " << Lambda << " " << "mu:" << " " << mu << " " << "ye1:" << " " << ye1 << " " << "ymu1:" << " " << ymu1 << " " << 
          "ytau1:" << " " << ytau1 << " " << "ye2:" << " " << ye2 << " " << "ymu2:" << " " << ymu2 << " " << "ytau2:" << " " << ytau2 << " " << 
          " " << "phie2:" << " " << phie2 << " " << "phimu2:" << " " << phimu2 << " " << "phitau2:" << " " << phitau2 << endl;
          }
          else{
               file << "# " << "Hierarchy:" << " " << hierarchy_string << " " << "Interactions:" << " " <<  lnv_rates_string << " " <<
               "M1:" << " " << M1 << " " << "DM/M:" << " " << DM_M << " " << "zre:" << " " << zre << " " << "zim:" << " " << zim << " " << 
               "delta:" << " " << delta << " " << "phi:" << " " << phi << endl;
          }
          

          file << "# " << "t" << " "  << "asym" << " " << "mu_e" << " " << "mu_mu" << " " << "mu_tau" << " " << "r11" << " " << "r22" << " " << "Re(r12)" << " " << "Im(r12)" << " " << "rb11" << " " << "rb22" << " " <<
          "Re(rb12)" << " " << "Im(r12)" << " " << "Smooth-Sphaleron-Asym" << endl;
          file << 0 << " "<< 0 << " "<< 0 << " " << 0 << " " << 0 << " "<< 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 <<  endl;
          for (int i = 0; i < row_size ; i++){
               double t = myres_matrix[i][0]; double asym = myres_matrix[i][1];
               double r11 = myres_matrix[i][2]; double r22 = myres_matrix[i][3]; double Rer12 = myres_matrix[i][4];
               double Imr12  = myres_matrix[i][5]; double rb11 = myres_matrix[i][6]; double rb22 = myres_matrix[i][7];
               double Rerb12  = myres_matrix[i][8]; double Imrb12 = myres_matrix[i][9]; double mue  = myres_matrix[i][10];
               double mumu = myres_matrix[i][11]; double mutau = myres_matrix[i][12]; double sp = myres_matrix[i][13];
               file << t - tmin << " "<< asym << " "<< mue << " " << mumu << " " << mutau << " "<< r11 << " " << r22 << " " << Rer12 << " " << Imr12 << " " << rb11 << " " << rb22 << " " << Rerb12 << " " <<
               Imrb12 << " " << sp <<  endl;
          }
          file.close();

          cout << "\n" << endl;
     }


     if (t < 0.95){
          /*
          If the solver failed return a 0 asymmetry as this is due to the exponential washout causing the numerical instability
          */
          cout << "[AMIQS] Solver failed at time x     " << t << endl;
          asym_inst      = 0;
          asym_smooth    = 0;
          cout << RED << "[AMIQS]" << RESET << " Asymmetry set to            " << asym_inst << endl;
     }
     else{
          /*
          assign final asymmetry generates 1) Instantenous Sphaleron freeze-out and 2) gradual, smooth Sphaleron freeze out
          */
          asym_inst      = _ASYMFAC_ * (yout[9] + yout[10] + yout[11]);
          asym_smooth    = yout[12];
          if (info == 1){
               cout << GREEN <<    "[AMIQS]" << RESET <<    " Instantaneous Asymmetry:    " << asym_inst << endl;
               cout << GREEN <<    "[AMIQS]" << RESET <<    " Smooth Asymmetry:           " << asym_smooth << endl;
               cout <<             "[AMIQS]" <<             " Deviation:                  " << asym_smooth/asym_inst << endl;
               if (parametrization == 0){
                    cout <<             "[AMIQS]" <<             " Input Parameters:           " << "M1=" << M1 << " " << "DM/M="<< DM_M << " " << "y=" << yukawa << " " << "th="<< theta << " " << "d=" << delta << " " << "phi=" << phi << endl;
               }
               else if (parametrization == 1){
                    cout <<             "[AMIQS]" <<             " Input Parameters:           " << "Lambda=" << Lambda << " " << "mu="<< mu << " " << "ye1=" << ye1 << " " << "ymu1="<< ymu1 << " " << "ytau1=" << ytau1 << " " << "ye2=" << ye2
                    << " " << "ymu2=" << ymu2 << " " << "ytau2=" << ytau2 << " " << "phie2=" << phie2 << " " << "phitau2=" << phitau2 << endl;
               }
               else{
                    cout <<             "[AMIQS]" <<             " Input Parameters:           " << "M1=" << M1 << " " << "DM/M="<< DM_M << " " << "zre=" << zre << " " << "zim="<< zim << " " << "d=" << delta << " " << "phi=" << phi << endl; 
               }
          }
     }
     double asym_final = asym_smooth;

     /*
     clear vectors and return final asymmetry
     */
     Yukawa.Yukawa_reset(YLN);
     Yukawa.Yukawa_reset(Y);
     y.clear();
     yout.clear();
     ini.clear();
     lsoda._freevectors();

     if (sph_approx == 0){
          return asym_final;
     }
     else{
          return asym_inst;
     }
}


/*
create matrix with random input parameters for testing purpose
*/
vector<vector<double>> rand_test(int regime, int test_num){
     srand (static_cast <unsigned> (time(0)));
     double low_angle, high_angle, low_y_alpha, high_y_alpha, low_y_beta, high_y_beta, high_yp, low_yp, low_Lambda, high_Lambda, low_mu, high_mu;
     double parametrization;

     double Lambda, mu, ye1, ymu1, ytau1, ye2, ymu2, ytau2, phie2, phimu2, phitau2;

     vector<double> param_ini;
     vector<double> vec_Lambda(test_num), vec_mu(test_num);
     vector<double> vec_ye1(test_num), vec_ymu1(test_num), vec_ytau1(test_num);
     vector<double> vec_ye2(test_num), vec_ymu2(test_num), vec_ytau2(test_num);
     vector<double> vec_phie2(test_num), vec_phimu2(test_num), vec_phitau2(test_num);

     /*
     testing routine using the analytical Lepton Number parametrization
     */
     parametrization= 1;

     /*
     three phases always between 0 and 2*PI
     */
     low_angle= 0.0; high_angle= 2.0*_PI_;
     /*
     depending on the regime generate a weak flavour alpha or strong flavour alpha
     0, 1 -> Overdamped with weak/strong LNV do not need flavour effects
     2, 3 -> Intermediate and Fast oscillations need flavour effects for sizable asymmetry
     */
     if (regime > 1){
          low_y_alpha= -8.5; high_y_alpha= -8.0; 
          low_mu= -8.0; high_mu= -1.0;
     } 
     else{
          low_y_alpha= -6.0; high_y_alpha= -5.0; 
          low_mu= -11.5; high_mu= -7.0;
     }
     /*
     perturbative parameter always in weak washout
     */
     low_yp= -9; high_yp= -8.5;
     /*
     flavour beta always in strong washout
     */
     low_y_beta= -6.0; high_y_beta= -5.0; 
     low_Lambda= -0.8; high_Lambda= 1.8;

     static int index_run;
     index_run= 0;
     while (index_run < test_num){
          Lambda       =  low_Lambda       + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(high_Lambda-low_Lambda))); 
          mu           =  low_mu           + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(high_mu-low_mu)));
          
          ye1          =  low_y_alpha      + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(high_y_alpha-low_y_alpha)));
          ymu1         =  low_y_beta       + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(high_y_beta-low_y_beta)));
          ytau1        =  low_y_beta       + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(high_y_beta-low_y_beta)));
          
          ye2          =  low_yp           + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(high_yp-low_yp)));
          ymu2         =  low_yp           + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(high_yp-low_yp)));
          ytau2        =  low_yp           + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(high_yp-low_yp)));

          phie2        =  low_angle        + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(high_angle-low_angle))); 
          phimu2       =  low_angle        + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(high_angle-low_angle))); 
          phitau2      =  low_angle        + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(high_angle-low_angle)));

          param_ini           = {Lambda, mu, ye1, ymu1, ytau1, ye2, ymu2, ytau2, phie2, phimu2, phitau2};
          void *ptr_param_ini = &param_ini;
          if (regime_check(regime, parametrization, ptr_param_ini)==true){
               vec_Lambda[index_run]= Lambda; vec_mu[index_run]= mu;
               vec_ye1[index_run]= ye1; vec_ymu1[index_run]= ymu1; vec_ytau1[index_run]= ytau1;
               vec_ye2[index_run]= ye2; vec_ymu2[index_run]= ymu2; vec_ytau2[index_run]= ytau2;
               vec_phie2[index_run]= phie2; vec_phimu2[index_run]= phimu2; vec_phitau2[index_run]= phitau2;
               index_run++;
               cout << " [AMIQS - TEST] Collecting testing points " << index_run << "/" << test_num << "\r" << flush;
          }
     }

     cout << " " << endl;
     return {vec_Lambda, vec_mu, vec_ye1, vec_ymu1, vec_ytau1, vec_ye2, vec_ymu2, vec_ytau2, vec_phie2, vec_phimu2, vec_phitau2};
}


/*
simplified testing routine performing 3 non-trivial checks
1) Comparison of numerical to analytical solution within the Overdamped weak LNV regime
2) Comparison of numerical to analytical solution within the Overdamped strong LNV regime
3) Comparison of numerical to analytical solution within the flavoured Intermediate weak/strong LNV regime
4) Comparison of numerical to analytical solution within the flavoured fast osc weak/strong LNV regime
*/
double amiqs_testing(void* setting_ini, int test_num){
     vector<vector<double>> rand_ini;
     vector<double> param_ini;
     vector<double> Lambda; vector<double> mu; 
     vector<double> ye1; vector<double> ymu1; vector<double> ytau1; 
     vector<double> ye2; vector<double> ymu2; vector<double> ytau2; 
     vector<double> phie2; vector<double> phimu2; vector<double> phitau2; 

     double progress, pos;
     int barWidth= 70;

     double parametrization;

     static double index_fail, index_succes;

     double dev;

     double ana, num;

     /*
     test using the analytical LN parametrization
     */
     parametrization= 1;

     /*
     acceptance ratio for ana vs num
     */
     dev= 5.;

     /*
     set running index to 0 which keep track of succes & failure of individual tests
     */
     index_fail= 0, index_succes= 0;



     /*
     Test the overdamped wLNV regime
     */
     progress = 0.0;
     rand_ini= rand_test(0, test_num);
     cout << " [AMIQS - TEST]" << " overdamped + wLNV " << endl;
     Lambda= rand_ini[0]; mu= rand_ini[1]; 
     ye1= rand_ini[2]; ymu1= rand_ini[3]; ytau1= rand_ini[4];
     ye2= rand_ini[5]; ymu2= rand_ini[6]; ytau2= rand_ini[7];
     phie2= rand_ini[8]; phimu2= rand_ini[9]; phitau2= rand_ini[10];
     for (int i=0; i<test_num; i++){
          param_ini = {Lambda[i], mu[i], ye1[i], ymu1[i], ytau1[i], ye2[i], ymu2[i], ytau2[i], phie2[i], phimu2[i], phitau2[i]};
          void *ptr_param_ini = &param_ini;
          ana= ana_approx.overdamped_wLNV(parametrization, ptr_param_ini);
          num= amiqs(parametrization, setting_ini, ptr_param_ini, 0, 0);
          if (num/ana < dev && num/ana > 1./dev){
               index_succes++;
          }
          else{
               index_fail++;
          }


          cout << " [";
          progress = ((double)i)/((double)test_num) + 1./(double)test_num; 
          pos = barWidth * progress - 1.;
          for (int k = 0; k < barWidth; ++k) {
               if (k < int(pos) ) cout << "=";
               else if (k == int(pos)) cout << ">";
               else cout << " ";
          }
          cout << "] " << int(progress * 100.0) << " %\r";
          cout.flush();
         
     }
     cout << endl;



     /*
     Test the overdamped sLNV regime
     */
     progress = 0.0;
     rand_ini= rand_test(1, test_num);
     cout << " [AMIQS - TEST]" << " overdamped + sLNV " << endl;
     Lambda= rand_ini[0]; mu= rand_ini[1]; 
     ye1= rand_ini[2]; ymu1= rand_ini[3]; ytau1= rand_ini[4];
     ye2= rand_ini[5]; ymu2= rand_ini[6]; ytau2= rand_ini[7];
     phie2= rand_ini[8]; phimu2= rand_ini[9]; phitau2= rand_ini[10];
     for (int i=0; i<test_num; i++){
          param_ini = {Lambda[i], mu[i], ye1[i], ymu1[i], ytau1[i], ye2[i], ymu2[i], ytau2[i], phie2[i], phimu2[i], phitau2[i]};
          void *ptr_param_ini = &param_ini;
          ana= ana_approx.overdamped_sLNV(parametrization, ptr_param_ini);
          num= amiqs(parametrization, setting_ini, ptr_param_ini, 0, 0);
          if (num/ana < dev && num/ana > 1./dev){
               index_succes++;
          }
          else{
               index_fail++;
          }

          cout << " [";
          progress = ((double)i)/((double)test_num) + 1./(double)test_num; 
          pos = barWidth * progress - 1.;
          for (int k = 0; k < barWidth; ++k) {
               if (k < int(pos) ) cout << "=";
               else if (k == int(pos)) cout << ">";
               else cout << " ";
          }
          cout << "] " << int(progress * 100.0) << " %\r";
          cout.flush();
     }
     cout << endl;



     /*
     Test the flavour weak intermediate regime
     */
     progress= 0;
     rand_ini= rand_test(2, test_num);
     cout << " [AMIQS - TEST]" << " flavoured intermediate " << endl;
     Lambda= rand_ini[0]; mu= rand_ini[1]; 
     ye1= rand_ini[2]; ymu1= rand_ini[3]; ytau1= rand_ini[4];
     ye2= rand_ini[5]; ymu2= rand_ini[6]; ytau2= rand_ini[7];
     phie2= rand_ini[8]; phimu2= rand_ini[9]; phitau2= rand_ini[10];
     for (int i=0; i<test_num; i++){
          param_ini = {Lambda[i], mu[i], ye1[i], ymu1[i], ytau1[i], ye2[i], ymu2[i], ytau2[i], phie2[i], phimu2[i], phitau2[i]};
          void *ptr_param_ini = &param_ini;
          ana= ana_approx.intermediate_fw(parametrization, ptr_param_ini);
          num= amiqs(parametrization, setting_ini, ptr_param_ini, 0, 0);
          if (num/ana < dev && num/ana > 1./dev){
               index_succes++;
          }
          else{
               index_fail++;
          }

          cout << " [";
          progress = ((double)i)/((double)test_num) + 1./(double)test_num; 
          pos = barWidth * progress - 1.;
          for (int k = 0; k < barWidth; ++k) {
               if (k < int(pos) ) cout << "=";
               else if (k == int(pos)) cout << ">";
               else cout << " ";
          }
          cout << "] " << int(progress * 100.0) << " %\r";
          cout.flush();
     }
     cout << endl;



     /*
     Test the flavour weak fast osc regime
     */
     progress= 0;
     rand_ini= rand_test(3, test_num);
     cout << " [AMIQS - TEST]" << " flavoured fast osc " << endl;
     Lambda= rand_ini[0]; mu= rand_ini[1]; 
     ye1= rand_ini[2]; ymu1= rand_ini[3]; ytau1= rand_ini[4];
     ye2= rand_ini[5]; ymu2= rand_ini[6]; ytau2= rand_ini[7];
     phie2= rand_ini[8]; phimu2= rand_ini[9]; phitau2= rand_ini[10];
     for (int i=0; i<test_num; i++){
          param_ini = {Lambda[i], mu[i], ye1[i], ymu1[i], ytau1[i], ye2[i], ymu2[i], ytau2[i], phie2[i], phimu2[i], phitau2[i]};
          void *ptr_param_ini = &param_ini;
          ana= ana_approx.fast_osc_fw(parametrization, ptr_param_ini);
          num= amiqs(parametrization, setting_ini, ptr_param_ini, 0, 0);
          if (num/ana < dev && num/ana > 1./dev){
               index_succes++;
          }
          else{
               index_fail++;
          }

          cout << " [";
          progress = ((double)i)/((double)test_num) + 1./(double)test_num; 
          pos = barWidth * progress - 1.;
          for (int k = 0; k < barWidth; ++k) {
               if (k < int(pos) ) cout << "=";
               else if (k == int(pos)) cout << ">";
               else cout << " ";
          }
          cout << "] " << int(progress * 100.0) << " %\r";
          cout.flush();
     }
     cout << endl;

     /***********************************************************
     Done with everything -- do not put anything below this line
     ***********************************************************/
     cout << " " << endl;
     if (index_fail/4/test_num < 0.15){
          cout << GREEN <<  "[AMIQS - TEST]" << RESET <<    " Succesfully passed                      " << endl;
          return _AMIQS_SUCCESS_;     
     }
     else{
          cout << RED <<   "[AMIQS - TEST]" << RESET <<    " Failed for " << index_fail/4/test_num*100 << " % of the points (expected <15%) "<< endl;
          exit(0);
     }


}



int main(int argc, const char* argv[]){
     (void)argc; (void)argv;
     auto start_main = chrono::system_clock::now();
     cout.precision(4);
     vector<double> setting_ini, param_ini;



     /*
     set the rates
     */
     set_rates();

     if(argc < 2){
          cout << RED <<  "[AMIQS]" << RESET <<    " Need to pass either:\n i) .ini file \n ii) test X where X > 0 " << endl;
          exit(0);
     }

     
     string arg1(argv[1]);
     bool found = arg1.find(".ini") != string::npos;
     if (found == true){
          vector<double> vars      = read_ini(argv[1]);
          if (vars[0] == 0){
               param_ini = {  vars[1], vars[2], vars[3], vars[4], vars[5], vars[6], vars[7]  };
          }
          else if(vars[0] == 1){
               param_ini = {  vars[8], vars[9], vars[10], vars[11], vars[12], vars[13], 
                              vars[14], vars[15], vars[16], vars[17], vars[18]  };
          }
          else if(vars[0] == 2){
               param_ini = {  vars[1], vars[2], vars[3], vars[27], vars[28], vars[6], vars[7]  };
          }
          else{
               cout << RED <<  "[AMIQS]" << RESET <<    " Need to pass correct parametrization in the .ini file; 0 or 1 or 2. " << endl;
               exit(0);
          }
          void *ptr_param_ini      = &param_ini;
          setting_ini              = {vars[19], vars[20], vars[21], vars[22], vars[23], vars[24]};
          void *ptr_setting_ini    = &setting_ini;
          amiqs(vars[0], ptr_setting_ini, ptr_param_ini, vars[25], vars[26]);
     }
     else{
          if (argv[1] == string("test")){
               /*
               fixed settings, dont change
               the analytical formulas are derived under these assumptions
               */
               setting_ini              = {1, 1, 1, 1, 0, 1};
               void *ptr_setting_ini    = &setting_ini;
               /*
               amount of testing points is fixed to some random number
               */
              double test_num= 10.;
               amiqs_testing(ptr_setting_ini, test_num);
          }
          else{
               cout << RED <<  "[AMIQS]" << RESET <<    " Need to pass either:\n i) .ini file \n ii) test " << endl;
               exit(0);       
          }
     }



     /***********************************************************
     Done with everything -- do not put anything below this line
     ***********************************************************/
     auto end_main = chrono::system_clock::now();
     chrono::duration<double> elapsed_seconds_main = end_main-start_main;
     cout <<"\n[AMIQS]" << " time: " << scientific << setprecision(4) << elapsed_seconds_main.count() << " seconds" << endl;
     /*
     Si arribes ací.... Xe que bo =)
     */
     // cout << CYAN << "Si arribes ací.... Xe que bo =)" << RESET << endl;
     return _AMIQS_SUCCESS_;
}
