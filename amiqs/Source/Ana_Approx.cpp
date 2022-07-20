class Ana_Approx{
private:
    double res;
    double hierarchy, Lambda, mu, M1, DM_M, d12, ye1, ymu1, ytau1, ye2, ymu2, ytau2, phie2, phimu2, phitau2, y_sum;
    double g0     = _G0_;
    double g1     = _G1_;
    double s0     = _S0_ ;
    double s1     = _S1_ ;
    double kappa        = _FD_;
    double ch           = _KAV_;
    double omega        = _MPLTSP_*ch/8.;
    double vth          = pow(rate_scaled(g0),2)/(pow(rate_scaled(g0),2)+pow(2*omega,2));
    double x_m_ov(double M1, double y_sum){
        return pow(1./3. * pow(10,M1)*pow(10,M1) * s0 * _MPL_/pow(_TSP_,3.) * pow(y_sum,2.), -1./3.);
    }      
    double d12_scaled(double d12){
        return ch/2.*_MPL_/pow(_TSP_,3.)*d12;
    }
    double rate_scaled(double x){
        return x*_MPLTSP_;
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
    initial conditions
    */
    vector<double> ini(double parametrization, void* param_ini){
        vector<double> param_var  = *(vector<double> *)param_ini;
        

        if (parametrization == 0){
            hierarchy= param_var[0]; M1 = param_var[1]; DM_M = param_var[2];
            d12                 = D12(M1, DM_M);
            Yukawa.YukawaLN(param_ini, YLN);
            phie2= arg(YLN[0][1]) - arg(YLN[0][0]); phimu2= arg(YLN[1][1]) - arg(YLN[1][0]); phitau2= arg(YLN[2][1]) - arg(YLN[2][0]);
        }
        else{
            Lambda= param_var[0]; mu= param_var[1]; 
            M1= Lambda_to_M1(Lambda, mu); DM_M= mu_to_DM_M(Lambda, mu);
            d12 = D12(M1, DM_M);
            Yukawa.Yukawa_anaLN(param_ini, YLN);
            phie2= arg(YLN[0][1]); phimu2= arg(YLN[1][1]); phitau2= arg(YLN[2][1]);
        }
        ye1= abs(YLN[0][0]); ymu1= abs(YLN[1][0]); ytau1= abs(YLN[2][0]);
        ye2= abs(YLN[0][1]); ymu2= abs(YLN[1][1]); ytau2= abs(YLN[2][1]);
        
        y_sum = sqrt( pow(abs(YLN[0][0]),2) + pow(abs(YLN[1][0]),2) + pow(abs(YLN[2][0]),2) );

        vector<double> ana_ini;
        
        ana_ini = {d12, ye1, ymu1, ytau1, ye2, ymu2, ytau2, phie2, phimu2, phitau2, y_sum};
        return ana_ini;
    }

public:
    double overdamped_wLNV(double parametrization, void* param_ini){
        vector<double> param_var      = *(vector<double> *)param_ini;
        if (parametrization == 0){
            M1                     = param_var[1];
        }
        else{
            M1 = Lambda_to_M1(param_var[0], param_var[1]);
        }
        vector<double> ana_ini = ini(parametrization, param_ini);
        d12= ana_ini[0]; 
        ye1= ana_ini[1]; ymu1= ana_ini[2]; ytau1= ana_ini[3]; 
        ye2= ana_ini[4]; ymu2= ana_ini[5]; ytau2= ana_ini[6]; 
        phie2= ana_ini[7]; phimu2= ana_ini[8]; phitau2= ana_ini[9];
        y_sum= ana_ini[10];


        double fac_lnc      = 4.*kappa*d12_scaled(d12)/(6*rate_scaled(g0) + kappa*rate_scaled(g1))*vth;
        double y_lnc        = 1/pow(y_sum,2.)*(ye1*ye2*sin(phie2)*(1/pow(ye1,2.)  - 3./pow(y_sum,2.)) + ymu1*ymu2*sin(phimu2)*(1/pow(ymu1,2.)  - 3./pow(y_sum,2.)) + ytau1*ytau2*sin(phitau2)*(1/pow(ytau1,2.)  - 3./pow(y_sum,2.)));
        
        double fac_lnv      = 48./5.*kappa*rate_scaled(s0)*d12_scaled(d12)/(6*rate_scaled(g0) + kappa*rate_scaled(g1))*vth*pow(pow(10,M1)/_TSP_,2.);
        double y_lnv        = 1/pow(y_sum,2.)*(ye1*ye2*sin(phie2) + ymu1*ymu2*sin(phimu2) + ytau1*ytau2*sin(phitau2));

        res          = _ASYMFAC_* (-fac_lnc*y_lnc + fac_lnv*y_lnv);


        Yukawa.Yukawa_reset(YLN);
        return res;
    }

    double overdamped_sLNV(double parametrization, void* param_ini){
        vector<double> param_var      = *(vector<double> *)param_ini;
        if (parametrization == 0){
            M1                     = param_var[1];
        }
        else{
            M1 = Lambda_to_M1(param_var[0], param_var[1]);
        }
        vector<double> ana_ini = ini(parametrization, param_ini);
        d12= ana_ini[0]; 
        ye1= ana_ini[1]; ymu1= ana_ini[2]; ytau1= ana_ini[3]; 
        ye2= ana_ini[4]; ymu2= ana_ini[5]; ytau2= ana_ini[6]; 
        phie2= ana_ini[7]; phimu2= ana_ini[8]; phitau2= ana_ini[9];
        y_sum= ana_ini[10];

        double fac_lnv      = 48./5.*kappa*pow(rate_scaled(s0),2.)/(6*rate_scaled(g0)*rate_scaled(s0) + kappa*(rate_scaled(g0)*rate_scaled(s1) + rate_scaled(g1)*rate_scaled(s0)))*vth*pow(pow(10,M1)/_TSP_,2)*d12_scaled(d12);
        double y_lnv        = 1/pow(y_sum,2.)*(ye1*ye2*sin(phie2) + ymu1*ymu2*sin(phimu2) + ytau1*ytau2*sin(phitau2));

        res          = _ASYMFAC_ * fac_lnv * y_lnv * pow(x_m_ov(M1, y_sum),5.);
        
        
        Yukawa.Yukawa_reset(YLN);
        return res;
    }

    double intermediate_fw(double parametrization, void* param_ini){
        vector<double> param_var      = *(vector<double> *)param_ini;
        if (parametrization == 0){
            M1                     = param_var[1];
        }
        else{
            M1 = Lambda_to_M1(param_var[0], param_var[1]);
        }
        vector<double> ana_ini = ini(parametrization, param_ini);
        d12= ana_ini[0]; 
        ye1= ana_ini[1]; ymu1= ana_ini[2]; ytau1= ana_ini[3]; 
        ye2= ana_ini[4]; ymu2= ana_ini[5]; ytau2= ana_ini[6]; 
        phie2= ana_ini[7]; phimu2= ana_ini[8]; phitau2= ana_ini[9];
        y_sum= ana_ini[10];



        double x_0= pow(5*(pow(rate_scaled(g0),2.) + 4*pow(omega,2.)) /(rate_scaled(g0)*pow(d12_scaled(d12),2.))*pow(y_sum,2.),1./5.);



        double gamma_m_slow     = 1./3. * (g1*s0 + g0*s1)*_FD_*pow(pow(10,M1),2)*_MPL_*pow(_TSP_,3) * pow(y_sum,2)/(3.*g0 + _FD_*g1);

        double fac_lnc_LN       = - 2./3.*d12_scaled(d12)*vth*2.*kappa*rate_scaled(g0)/(2*rate_scaled(g0) + kappa*rate_scaled(g1));
        double fac_lnc_alpha    = 2./3.*d12_scaled(d12)*vth*kappa;
        double y_lnc            = 1./pow(y_sum,4) *( pow(ye1,2) *(ymu1 * ymu2 * sin(phimu2) + ytau1 * ytau2 * sin(phitau2)) - (ye1 * ye2 * sin(phie2)*( pow(ymu1,2) + pow(ytau1,2))));


        res                     = _ASYMFAC_* pow(x_0,3.) *( fac_lnc_LN*y_lnc*exp(-gamma_m_slow*( 1 - pow(x_0,3.))) + fac_lnc_alpha*y_lnc);


        Yukawa.Yukawa_reset(YLN);
        return res;
    }


    double fast_osc_fw(double parametrization, void* param_ini){
        vector<double> param_var      = *(vector<double> *)param_ini;
        if (parametrization == 0){
            M1                     = param_var[1];
        }
        else{
            M1 = Lambda_to_M1(param_var[0], param_var[1]);
        }
        vector<double> ana_ini = ini(parametrization, param_ini);
        d12= ana_ini[0]; 
        ye1= ana_ini[1]; ymu1= ana_ini[2]; ytau1= ana_ini[3]; 
        ye2= ana_ini[4]; ymu2= ana_ini[5]; ytau2= ana_ini[6]; 
        phie2= ana_ini[7]; phimu2= ana_ini[8]; phitau2= ana_ini[9];
        y_sum= ana_ini[10];


        double x_0= pow(5*(pow(rate_scaled(g0),2.) + 4*pow(omega,2.)) /(rate_scaled(g0)*pow(d12_scaled(d12),2.))*pow(y_sum,2.),1./5.);

        double integral_J_200= 1.436463044073161/pow(d12_scaled(d12),2./3.);

        double gamma_m_slow     = 1./3. * (g1*s0 + g0*s1)*_FD_*pow(pow(10,M1),2)*_MPL_*pow(_TSP_,3) * pow(y_sum,2)/(3.*g0 + _FD_*g1);

        double fac_lnc_LN  = -2.*pow(rate_scaled(g0),3.)*kappa/(2*rate_scaled(g0)+rate_scaled(g1)*kappa)*integral_J_200;
        double fac_lnc_alpha = pow(rate_scaled(g0),2)*kappa*integral_J_200;
        double y_lnc            =  pow(ye1,2) *(ymu1 * ymu2 * sin(phimu2) + ytau1 * ytau2 * sin(phitau2)) - (ye1 * ye2 * sin(phie2)*( pow(ymu1,2) + pow(ytau1,2)));

        res                     = _ASYMFAC_* ( fac_lnc_LN*y_lnc*exp(-gamma_m_slow*( 1 - pow(x_0,3.))) + fac_lnc_alpha*y_lnc);


        Yukawa.Yukawa_reset(YLN);
        return res;    
    }
};

/*
initialize globally Ana_Approx class
*/
Ana_Approx ana_approx;