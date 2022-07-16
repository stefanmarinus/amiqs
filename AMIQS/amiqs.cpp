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
bool is_number(const std::string& s){
    return !s.empty() && std::find_if(s.begin(), 
        s.end(), [](unsigned char c) { return !std::isdigit(c); }) == s.end();
}


/*log scale input but return in linear scale*/
double D12(void* param_ini){
     vector<double> param_var      = *(vector<double> *)param_ini;
     double M1                     = pow(10,param_var[0]);
     double DM_M                   = pow(10,param_var[1]);
     double M2                     = M1 * (DM_M + 1.);

     return pow(M2,2) - pow(M1,2);
}

/*
Calculate some transition times which define different washout regimes
*/
vector<double> transition_times(int hierarchy, double hv_interactions, void* param_ini){
     /*
     declare variables used
     */
     vector<double> param_var;
     double M1;
     double d12;
     double y_sum;



     param_var = *(vector<double> *)param_ini;
     M1        = param_var[0]; 


     d12       = D12(param_ini);

     Yukawa.YukawaLN(hierarchy, param_ini, YLN);

     y_sum = sqrt( pow(abs(YLN[0][0]),2) + pow(abs(YLN[1][0]),2) + pow(abs(YLN[2][0]),2) );


     if (y_sum > 0){
          double x_transition = pow((16.*pow(rates.hc0(pow(10,M1),1.),2) + pow(_KAV_,2))*(5*pow(_TSP_,5.)*(9.*rates.hc0(pow(10,M1),1.)+rates.hc1(pow(10,M1),1.)*_FD_)*pow(y_sum,2))
                              /(4.*_MPL_*pow(d12,2)*rates.hc0(pow(10,M1),1.)*(9.*rates.hc0(pow(10,M1),1.)+2.*rates.hc1(pow(10,M1),1.)*_FD_)*pow(_KAV_,2)),1./5.);
          double x_0          = x_transition;
          double x_osc        =  _PI_/3.*rates.hc0(pow(10,M1), 0.001)*pow(y_sum,2)*pow(2/_KAV_,1./3.)*pow(_MPL_*_MPL_/d12,1./3.);
          double s0           = rates.hv0(pow(10,M1), 0.001) * pow(_TSP_/0.001, 2.) * pow(1./pow(10,M1),2.);
          double s1           = rates.hv1(pow(10,M1), 0.001) * pow(_TSP_/0.001, 2.) * pow(1./pow(10,M1),2.);
          double x_lnv        = pow(1./3. * pow(10,M1)*pow(10,M1) * s0 * _MPL_/pow(_TSP_,3.) * pow(y_sum,2.), -1./3.); //pow(1./3. * pow(10,M1)*pow(10,M1) * rates.hv0(pow(10,M1), 0.001) * _MPL_/pow(_TSP_,3.) * pow(y_sum,2.) ,-1./3.);
          double x_lnv_int    = pow(2./3.* (rates.hc1(pow(10,M1), 0.001)* s0 + rates.hc0(pow(10,M1), 0.001)*s1*_FD_*pow(10,M1)*pow(10,M1)*_MPL_/pow(_TSP_,3.) *pow(y_sum,2.)/(9.*rates.hc0(pow(10,M1), 0.001) + 2.*_FD_*rates.hc1(pow(10,M1), 0.001))) ,-1./3.);
          double ww           = _MPL_/_TSP_ * (2. * (3. *rates.hc0(pow(10,M1), 0.001) * rates.hc1(pow(10,M1), 0.001) * _FD_ + pow(rates.hc1(pow(10,M1), 0.001)*_FD_,2)))/(3. * (3. * rates.hc0(pow(10,M1), 0.001) + 2. * rates.hc1(pow(10,M1), 0.001) * _FD_))*1/2.*pow(y_sum,2);

          double sw_e         = _MPLTSP_*1./2.*rates.hc1(pow(10,M1),1.)*pow(abs(YLN[0][0]),2);
          double sw_mu        = _MPLTSP_*1./2.*rates.hc1(pow(10,M1),1.)*pow(abs(YLN[1][0]),2);
          double sw_tau       = _MPLTSP_*1./2.*rates.hc1(pow(10,M1),1.)*pow(abs(YLN[2][0]),2);

          return {x_0, x_osc, x_lnv, x_lnv_int, ww, sw_e, sw_mu, sw_tau};
     }
     else{
          return {0, 0, 0, 0, 0, 0, 0, 0};
     }
}

/*
function which prints the information of the transition times and the corresponding washout regime
*/
void print_transition_times(int hierarchy, double hv_interactions, void* param_ini){
     /*
     declare variables
     */
     vector<double> x_vec;
     double x_0, x_osc, x_lnv, x_lnv_int, ww;


     x_vec = transition_times(hierarchy, hv_interactions, param_ini);
     x_0= x_vec[0]; x_osc= x_vec[1]; x_lnv= x_vec[2]; x_lnv_int= x_vec[3]; ww= x_vec[4];


     if (x_0>0 && x_osc>0 && x_lnv>0 && x_lnv_int>0 && ww>0){

          cout << " [AMIQS] x_0:                       " << x_0 << endl;
          cout << " [AMIQS] x_osc:                     " << x_osc << endl;
          cout << " [AMIQS] x_M^ov:                    " << x_lnv << endl;
          cout << " [AMIQS] x_M^int:                   " << x_lnv_int << endl;



          if(hv_interactions == 1){
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
function joining all the stuff to call and solve the kinetic equaiton 
*/
double amiqs(double hierarchy, double hv_interactions, void* setting_ini, void* param_ini, double info, double safe){
     vector<double> param_var      = *(vector<double> *)param_ini;
     vector<double> setting_var       = *(vector<double> *)setting_ini;
     double M1, DM_M, yukawa, theta, delta, phi;
     double nl_approx, rates_approx, c_matrix_approx, g1_approx, sph_approx;
     double asym_inst, asym_smooth;

     M1 = param_var[0]; DM_M = param_var[1]; yukawa = param_var[2]; theta = param_var[3]; delta = param_var[4]; phi = param_var[5];
     nl_approx= setting_var[0]; rates_approx= setting_var[1]; c_matrix_approx= setting_var[2]; g1_approx= setting_var[3], sph_approx= setting_var[4];
     /*
     Number of equations and time span to solve
     */
     int neq = 12;  double tmin = 0.00002; double t = tmin; double tout = tmin*1.001;

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
     double osc_freq_max = pow(10,5);

     /*
     print some info
     */
     if (info == 1){
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
          if (hv_interactions == 0){
               cout << "[AMIQS] LN Conserving                    " << endl;
          }
          else if (hv_interactions == 1){
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
     double d12          = D12(param_ini);

     /*
     get the Yukawa matrices for given input and check if LN constraint is fulfilled
     */
     Yukawa.YukawaLN(hierarchy, param_ini, YLN);
     Yukawa.Yukawa(YLN, Y); Yukawa.YukawaHC(Y,YHC); Yukawa.YukawaT(Y,YT); Yukawa.YukawaC(Y,YC);


     /*
     print yukawa matrices
     no fa falta pro bo...
     */
     if (info == 1){
          cout << "" << endl;
          for(int i = 0; i<3; i++){
               for(int j = 0; j<2; j++){
                    cout << "YLN[" << i << "][" << j << "] = " << YLN[i][j] << endl;
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
    double ratio        = _KAV_*d12/pow(_TSP_,2) / pow(y_sum,2) / rates.hc0(pow(10,M1), 0.001);
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
    if (hv_interactions == 0){
          if (1./2.*_MPLTSP_*rates.hc0(pow(10,M1),0.001)*pow(min(abs(YLN[0][0]), min(abs(YLN[1][0]), abs(YLN[2][0]))),2) > 20  && ratio > 10){
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
     if (hv_interactions == 1){
          if (1./2.*_MPLTSP_*rates.hc0(pow(10,M1),0.001)*pow(min(abs(YLN[0][0]), min(abs(YLN[1][0]), abs(YLN[2][0]))),2) > 20 &&
              1./2.*_MPLTSP_*rates.hv0(pow(10,M1),0.001)*pow(min(abs(YLN[0][0]), min(abs(YLN[1][0]), abs(YLN[2][0]))),2) > 20){
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
     double x_transition = pow((16.*pow(rates.hc0(pow(10,M1),1.),2) + pow(_KAV_,2))*(5*pow(_TSP_,5.)*(9.*rates.hc0(pow(10,M1),1.)+rates.hc1(pow(10,M1),1.)*_FD_)*pow(y_sum,2))
                         /(4.*_MPL_*pow(d12,2)*rates.hc0(pow(10,M1),1.)*(9.*rates.hc0(pow(10,M1),1.)+2.*rates.hc1(pow(10,M1),1.)*_FD_)*pow(_KAV_,2)),1./5.);
     double x_ave        =  x_transition + min(x_osc_freq(osc_freq_max,d12), x_osc_num(osc_num_max,d12));
     if (info == 1){
          cout << CYAN <<  "[AMIQS]" << RESET << " Averaging x                 " << x_ave << endl;
     }

     /*
     call this function to print washout regime information
     if you dont need this info just comment the function call
     */
     if (info == 1){
          print_transition_times(hierarchy, hv_interactions, param_ini);
     }

     /*
     store doubles in a pointer which can be passed to the solver
     */
     vector<double> ini  = {pow(10,M1), d12, x_ave, hv_interactions, nl_approx, rates_approx, c_matrix_approx, g1_approx};
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
     
     if (safe == 1){
          /*
          write to file the calculated asymmetry for given input parameters
          */
          ofstream file;
          if (rates_approx == 1 && c_matrix_approx == 1){
               if(hv_interactions == 0){
                    file.open("../output/asym_approx_results_lnc.dat");
               }
               else{
                    file.open("../output/asym_approx_results_lnv.dat");
               }
          }
          else{
               if(hv_interactions == 0){
                    file.open("../output/asym_results_lnc.dat");
               }
               else{
                    file.open("../output/asym_results_lnv.dat");  
               }
          }
          char hierarchy_string [10];
          char hv_interactions_string [10];
          if(hierarchy == 0){
               // hierarchy_string = "NH";
               strcpy (hierarchy_string,"NH");
          }
          else{
               strcpy (hierarchy_string,"IH");
          }
          if(hv_interactions == 0){
               strcpy (hv_interactions_string,"LNC");
          }
          else{
               strcpy (hv_interactions_string,"LNC+LNV");
          }

          file << "# " << "Hierarchy:" << " " << hierarchy_string << " " << "Interactions:" << " " <<  hv_interactions_string << " " <<
          "M1:" << " " << M1 << " " << "DM/M" << " " << DM_M << " " << "y" << " " << yukawa << " " << "theta" << " " << theta << " " << 
          "delta" << " " << delta << " " << "phi" << " " << phi << endl;
          
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
               cout <<             "[AMIQS]" <<             " Input Parameters:           " << "M1=" << M1 << " " << "DM/M="<< DM_M << " " << "y=" << yukawa << " " << "th="<< theta << " " << "d=" << delta << " " << "phi=" << phi << endl;
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
vector<vector<double>> rand_test(double len_vec){
     srand (static_cast <unsigned> (time(0)));
     double low_angle, high_angle, low_y, high_y, low_dmm, high_dmm, low_m1, high_m1;
     vector<double> vec_m1(len_vec), vec_dmm(len_vec), vec_y(len_vec), vec_theta(len_vec), vec_delta(len_vec), vec_phi(len_vec);

     low_angle= 0.0; high_angle= 2.0*_PI_; low_y= -8.0; high_y= -4.0; low_dmm= -14.0; high_dmm= -1.0; low_m1= -0.9; high_m1= 1.9;

     for (int i = 0; i<len_vec; i++){
          vec_m1[i]      =  low_m1      + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(high_m1-low_m1))); 
          vec_dmm[i]     =  low_dmm     + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(high_dmm-low_dmm)));
          vec_y[i]       =  low_y       + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(high_y-low_y)));
          vec_theta[i]   =  low_angle   + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(high_angle-low_angle))); 
          vec_delta[i]   =  low_angle   + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(high_angle-low_angle))); 
          vec_phi[i]   =  low_angle   + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(high_angle-low_angle)));
     }

    return {vec_m1, vec_dmm, vec_y, vec_theta, vec_delta, vec_phi};
}

/*
simplified testing routine performing 3 non-trivial checks
1) Comparison of numerical to analytical solution within the Overdamped weak LNV regime
2) Comparison of numerical to analytical solution within the Overdamped strong LNV regime
3) Check numerical for vanishing asymmetry if there is no weak mode
*/
double amiqs_testing(double rand_num, void* setting_ini){
     /*
     declare variables used here
     */
     vector<vector<double>> rand_ini;
     vector<double> param_ini;
     vector<double> regimes;
     vector<double> M1; vector<double> DM_M; vector<double> yukawa; vector<double> theta; vector<double> delta; vector<double> phi;
     double ov_num_wLNV, ov_ana_wLNV, ov_num_sLNV, ov_ana_sLNV, white_num;
     double x_0, x_lnv, x_lnv_int, sw_e, sw_mu, sw_tau;
     double hv_interactions, dev;
     int hierarchy;
     static double index_fail, index_succes, index_test;
     static int index_run;

     
     /*
     get the random input to test, some dummy number of variables
     and store the result in individual vectors
     */
     rand_ini= rand_test(pow(10,5));
     M1= rand_ini[0]; DM_M= rand_ini[1]; yukawa= rand_ini[2]; theta= rand_ini[3]; delta= rand_ini[4]; phi= rand_ini[5];

     /*
     set LNV interactions to true
     */
     hv_interactions= 1;

     /*
     acceptance ratio for ana vs num
     */
     dev= 4.;

     /*
     set running index to 0 which keep track of succes & failure of individual tests
     */
     index_fail= 0, index_succes= 0, index_test= 0, index_run= 0;
     


     /*
     start the test
     */
     while (index_test < rand_num){
          /*
          half of the test for NH rest for IH
          */
          if (index_test < round(rand_num/2)){
               hierarchy= 0;
          }
          else{
               hierarchy= 1;
          }

          /*
          define i. input parameters
          */
          param_ini = {M1[index_run], DM_M[index_run], yukawa[index_run], theta[index_run], delta[index_run], phi[index_run]};
          void *ptr_param_ini = &param_ini;

          /*
          get all relevant transition times to check in which washout regime we are
          */
          regimes= transition_times(hierarchy, hv_interactions, ptr_param_ini);
          x_0= regimes[0]; x_lnv= regimes[2]; x_lnv_int= regimes[3]; sw_e= regimes[5], sw_mu= regimes[6], sw_tau= regimes[7];

          /*
          if point falls inside the overdamped + wLNV regime -> test
          */
          if (x_0 > 1 && x_lnv > 1 && sw_e > 1 && sw_mu > 1 && sw_tau > 1) {
               ov_ana_wLNV    = ana_approx.overdamped_wLNV(hierarchy, hv_interactions, ptr_param_ini);
               if (abs(ov_ana_wLNV) > pow(10,-12)){
                    index_test++;
                    ov_num_wLNV    = amiqs(hierarchy, hv_interactions, setting_ini, ptr_param_ini, 0, 0);
                    cout << "[AMIQS]" << " Testing routine " << index_test << "/" << rand_num << "\r" << flush;
                    if (ov_num_wLNV/ov_ana_wLNV < dev && ov_num_wLNV/ov_ana_wLNV > 1./dev){
                         index_succes++;
                    }
                    else{
                         index_fail++;
                    }
               }
          }

          /*
          if point falls inside the overdamped + sLNV regime -> test
          */
          if (x_0 > 1 && x_lnv < 1 && sw_e > 1 && sw_mu > 1 && sw_tau > 1) {
               ov_ana_sLNV    = ana_approx.overdamped_sLNV(hierarchy, hv_interactions, ptr_param_ini);
               if (abs(ov_ana_sLNV) > pow(10,-12)){
                    index_test++;
                    ov_num_sLNV    = amiqs(hierarchy, hv_interactions, setting_ini, ptr_param_ini, 0, 0);
                    cout << "[AMIQS]" << " Testing routine " << index_test << "/" << rand_num << "\r" << flush;
                    if (ov_num_sLNV/ov_ana_sLNV < dev && ov_num_sLNV/ov_ana_sLNV > 1./dev){
                         index_succes++;
                    }
                    else{
                         index_fail++;
                    }
               }
          /*
          for points putside the overdamped, with sLNV and strong flavour washout the asymmetry should be zero
          */
          if (x_0 < 1 && x_lnv_int < 1 && sw_e > 1 && sw_mu > 1 && sw_tau > 1){
               index_test++;
               white_num = amiqs(hierarchy, hv_interactions, setting_ini, ptr_param_ini, false, false);
               cout << "[AMIQS]" << " Testing routine " << index_test << "/" << rand_num << "\r" << std::flush;
               if (abs(white_num) < pow(10,-13)){
                    index_succes++;
               }
               else{
                    index_fail++;
               }
          }
          }
          index_run++;
     }




     /***********************************************************
     Done with everything -- do not put anything below this line
     ***********************************************************/
     if (index_fail/index_test < 0.15){
          cout << GREEN <<  "[AMIQS]" << RESET <<    " Test succesfully passed                      " << endl;
          // cout << GREEN <<  "[AMIQS]" << RESET <<    " Test passed with " << index_succes/index_test*100 << " % accurary (expected > 85%) " << endl;
          return _AMIQS_SUCCESS_;     
     }
     else{
          cout << RED <<   "[AMIQS]" << RESET <<    " Test failed for " << index_fail/index_test*100 << " % of the points (expected <15%) "<< endl;
          cout << "Try to increase the number of test points to generate more statistical significance." << endl;
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
     bool found = arg1.find(".ini") != std::string::npos;
     if (found == true){
          vector<double> vars = read_ini(argv[1]);
          param_ini                = {vars[2], vars[3], vars[4], vars[5], vars[6], vars[7]};
          void *ptr_param_ini      = &param_ini;
          setting_ini              = {vars[8], vars[9], vars[10], vars[11], vars[12]};
          void *ptr_setting_ini    = &setting_ini;
          amiqs(vars[0], vars[1], ptr_setting_ini, ptr_param_ini, vars[13], vars[14]);
     }
     else{
          if (argv[1] == string("test")){
               if (argc > 2){
                    if (is_number(argv[2])){
                         /*
                         fixed settings, dont change
                         the analytical formulas are derived under these assumptions
                         */
                         setting_ini              = {1, 1, 1, 0, 1};
                         void *ptr_setting_ini    = &setting_ini;
                         amiqs_testing(atof(argv[2]), ptr_setting_ini);
                    }
                    else{
                         cout << RED <<  "[AMIQS]" << RESET <<    " Need to pass either:\n i) .ini file \n ii) test X where X > 0 " << endl;
                         exit(0);       
                    }
               }
               else{
                    cout << RED <<  "[AMIQS]" << RESET <<    " Need to pass either:\n i) .ini file \n ii) test X where X > 0 " << endl;
                    exit(0);
               }
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
