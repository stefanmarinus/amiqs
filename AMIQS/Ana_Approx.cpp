#include <complex>
#include "common.h"
// #include "Yukawas.cpp"

using namespace std;

class Ana_Approx{
private:
    double d12, ye1, ymu1, ytau1, ye2, ymu2, ytau2, phie2, phimu2, phitau2, y_sum;
    double g0           = 0.0091;
    double g1           = 0.0051;
    double s0           = 0.0434;
    double s1           = 0.0086;
    double kappa        = 1.096144454102158;
    double ch           = 0.4561442592067353;
    double omega        = _MPLTSP_*ch/8.;
    double vth          = pow(rate_scaled(g0),2)/(pow(rate_scaled(g0),2)+pow(2*omega,2));
    double xmov(double M1, double y_sum){
        return pow(1./3. * pow(10,M1)*pow(10,M1) * s0 * _MPL_/pow(_TSP_,3.) * pow(y_sum,2.), -1./3.);
    }      
    double d12_scaled(double d12){
        return ch/2.*_MPL_/pow(_TSP_,3.)*d12;
    }
    double rate_scaled(double x){
        return x*_MPLTSP_;
    }
    /*log scale input but return in linear scale*/
    double D12(void* param_ini){
        vector<double> param_var      = *(vector<double> *)param_ini;
        double M1                     = pow(10,param_var[0]);
        double DM_M                   = pow(10,param_var[1]);
        double M2                     = M1 * (DM_M + 1.);

        return pow(M2,2) - pow(M1,2);
    }
    vector<double>  ini(int hierarchy, void* param_ini){
        d12                 = D12(param_ini);
        Yukawa.YukawaLN(hierarchy, param_ini, YLN);
        ye1= abs(YLN[0][0]); ymu1= abs(YLN[1][0]); ytau1= abs(YLN[2][0]);
        ye2= abs(YLN[0][1]); ymu2= abs(YLN[1][1]); ytau2= abs(YLN[2][1]);
        phie2= arg(YLN[0][1]) - arg(YLN[0][0]); phimu2= arg(YLN[1][1]) - arg(YLN[1][0]); phitau2= arg(YLN[2][1]) - arg(YLN[2][0]);
        y_sum = sqrt( pow(abs(YLN[0][0]),2) + pow(abs(YLN[1][0]),2) + pow(abs(YLN[2][0]),2) );

        vector<double> ana_ini;
        ana_ini = {d12, ye1, ymu1, ytau1, ye2, ymu2, ytau2, phie2, phimu2, phitau2, y_sum};
        return ana_ini;
    }

public:
    double overdamped_wLNV(int hierarchy, double hv_interactions, void* param_ini){
        vector<double> param_var      = *(vector<double> *)param_ini;
        double M1                     = param_var[0];
        vector<double> ana_ini = ini(hierarchy, param_ini);
        d12= ana_ini[0]; 
        ye1= ana_ini[1]; ymu1= ana_ini[2]; ytau1= ana_ini[3]; 
        ye2= ana_ini[4]; ymu2= ana_ini[5]; ytau2= ana_ini[6]; 
        phie2= ana_ini[7]; phimu2= ana_ini[8]; phitau2= ana_ini[9];
        y_sum= ana_ini[10];


        double fac_lnc      = 4.*kappa*d12_scaled(d12)/(6*rate_scaled(g0) + kappa*rate_scaled(g1))*vth;
        double y_lnc        = 1/pow(y_sum,2.)*(ye1*ye2*sin(phie2)*(1/pow(ye1,2.)  - 3./pow(y_sum,2.)) + ymu1*ymu2*sin(phimu2)*(1/pow(ymu1,2.)  - 3./pow(y_sum,2.)) + ytau1*ytau2*sin(phitau2)*(1/pow(ytau1,2.)  - 3./pow(y_sum,2.)));
        
        double fac_lnv      = 48./5.*kappa*rate_scaled(s0)*d12_scaled(d12)/(6*rate_scaled(g0) + kappa*rate_scaled(g1))*vth*pow(pow(10,M1)/_TSP_,2.);
        double y_lnv        = 1/pow(y_sum,2.)*(ye1*ye2*sin(phie2) + ymu1*ymu2*sin(phimu2) + ytau1*ytau2*sin(phitau2));

        double res          = _ASYMFAC_* (-fac_lnc*y_lnc + fac_lnv*y_lnv);


        Yukawa.Yukawa_reset(YLN);
        return res;
    }

    double overdamped_sLNV(int hierarchy, double hv_interactions, void* param_ini){
        vector<double> param_var      = *(vector<double> *)param_ini;
        double M1                     = param_var[0];
        vector<double> ana_ini = ini(hierarchy, param_ini);
        d12= ana_ini[0]; 
        ye1= ana_ini[1]; ymu1= ana_ini[2]; ytau1= ana_ini[3]; 
        ye2= ana_ini[4]; ymu2= ana_ini[5]; ytau2= ana_ini[6]; 
        phie2= ana_ini[7]; phimu2= ana_ini[8]; phitau2= ana_ini[9];
        y_sum= ana_ini[10];

        double fac_lnv      = 48./5.*kappa*pow(rate_scaled(s0),2.)/(6*rate_scaled(g0)*rate_scaled(s0) + kappa*(rate_scaled(g0)*rate_scaled(s1) + rate_scaled(g1)*rate_scaled(s0)))*vth*pow(pow(10,M1)/_TSP_,2)*d12_scaled(d12);
        double y_lnv        = 1/pow(y_sum,2.)*(ye1*ye2*sin(phie2) + ymu1*ymu2*sin(phimu2) + ytau1*ytau2*sin(phitau2));

        double res          = _ASYMFAC_ * fac_lnv * y_lnv * pow(xmov(M1, y_sum),5.);
        
        
        Yukawa.Yukawa_reset(YLN);
        return res;
    }
};

/*
initialize globally Ana_Approx class
*/
Ana_Approx ana_approx;