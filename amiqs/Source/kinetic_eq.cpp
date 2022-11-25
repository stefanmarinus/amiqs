double yaux[11];
double sph;

/*
Right Hand Side of the kinetic equation
*/
void RHS(double t, double* y, void* data, double (&yaux)[11]){
        vector<double> ini = *(vector<double> *)data;
        double M1, d12, x_ave, lnv_rates, nl_approx, rates_approx, c_matrix_approx, g1_approx, mthM_approx;
        M1= ini[0]; d12= ini[1]; x_ave= ini[2]; lnv_rates= ini[3]; 
        nl_approx= ini[4]; rates_approx= ini[5]; c_matrix_approx= ini[6], g1_approx= ini[7]; mthM_approx= ini[8];

        double g0, g1, g2, s0, s1 ,s2, c_aa, c_ab, g1_r_to_mu, thmlnv, VNM;
        if (rates_approx == 1){
                g0     = _G0_;
                g1     = _G1_;
                g2     = _G2_;
                s0     = _S0_   * pow(M1*t/_TSP_,2.);
                s1     = _S1_   * pow(M1*t/_TSP_,2.);
                s2     = _S2_   * pow(M1*t/_TSP_,2.);
        }
        else{
                g0      = rates.lnc_0(M1,t);
                g1      = rates.lnc_1(M1,t);
                g2      = rates.lnc_2(M1,t);
                s0      = rates.lnv_0(M1,t);
                s1      = rates.lnv_1(M1,t);
                s2      = rates.lnv_2(M1,t);
        }
        if (nl_approx == 1){
                g2      = 0.;
                s2      = 0.;
        }
        if (c_matrix_approx == 1){
                c_aa    = 1./2.;
                c_ab    = 0.;
        }
        else{
                c_aa    = _CABDIAG_;
                c_ab    = _CABOFF_;
        }
        if (g1_approx == 1){
                g1_r_to_mu = 0;
        }
        else{
                g1_r_to_mu = 1;
        }
        if (mthM_approx == 1){
                thmlnv  = 0;
                VNM     = 0;
        }
        else{
                thmlnv  = 1;
                VNM     = thermal_mass_lnv.mthM(M1,t);
        }

        
        double y00re, y00im, y01re, y01im, y10re, y10im, y11re, y11im, y20re, y20im, y21re, y21im;
        y00re = real(Y[0][0]); y00im = imag(Y[0][0]); y01re = real(Y[0][1]); y01im = imag(Y[0][1]);
        y10re = real(Y[1][0]); y10im = imag(Y[1][0]); y11re = real(Y[1][1]); y11im = imag(Y[1][1]);
        y20re = real(Y[2][0]); y20im = imag(Y[2][0]); y21re = real(Y[2][1]); y21im = imag(Y[2][1]);

        if (t > 1.1*x_ave){
                y[2] = 0; y[3] = 0; y[6] = 0; y[7] = 0;
        }

        yaux[0]         = -0.25*(_MPLTSP_*(4*c_aa*g1_r_to_mu*g1*y[8]*pow(y00im,2.) - 4*c_aa*s1*lnv_rates*y[8]*pow(y00im,2.) + 
                        4*c_ab*g1_r_to_mu*g1*y[9]*pow(y00im,2.) - 4*c_ab*s1*lnv_rates*y[9]*pow(y00im,2.) + 
                        4*c_ab*g1_r_to_mu*g1*y[10]*pow(y00im,2.) - 4*c_ab*s1*lnv_rates*y[10]*pow(y00im,2.) - 
                        4*c_aa*g2*y[8]*y[0]*pow(y00im,2.) + 4*c_aa*s2*lnv_rates*y[8]*y[0]*pow(y00im,2.) - 
                        4*c_ab*g2*y[9]*y[0]*pow(y00im,2.) + 4*c_ab*s2*lnv_rates*y[9]*y[0]*pow(y00im,2.) - 
                        4*c_ab*g2*y[10]*y[0]*pow(y00im,2.) + 4*c_ab*s2*lnv_rates*y[10]*y[0]*pow(y00im,2.) + 
                        4*c_aa*g1_r_to_mu*g1*y[8]*pow(y00re,2.) - 4*c_aa*s1*lnv_rates*y[8]*pow(y00re,2.) + 
                        4*c_ab*g1_r_to_mu*g1*y[9]*pow(y00re,2.) - 4*c_ab*s1*lnv_rates*y[9]*pow(y00re,2.) + 
                        4*c_ab*g1_r_to_mu*g1*y[10]*pow(y00re,2.) - 4*c_ab*s1*lnv_rates*y[10]*pow(y00re,2.) - 
                        4*c_aa*g2*y[8]*y[0]*pow(y00re,2.) + 4*c_aa*s2*lnv_rates*y[8]*y[0]*pow(y00re,2.) - 
                        4*c_ab*g2*y[9]*y[0]*pow(y00re,2.) + 4*c_ab*s2*lnv_rates*y[9]*y[0]*pow(y00re,2.) - 
                        4*c_ab*g2*y[10]*y[0]*pow(y00re,2.) + 4*c_ab*s2*lnv_rates*y[10]*y[0]*pow(y00re,2.) + 
                        _KAV_*y[3]*y00im*y01im - 4*c_aa*g2*y[8]*y[2]*y00im*y01im + 
                        4*c_aa*s2*lnv_rates*y[8]*y[2]*y00im*y01im - 4*c_ab*g2*y[9]*y[2]*y00im*y01im + 
                        4*c_ab*s2*lnv_rates*y[9]*y[2]*y00im*y01im - 4*c_ab*g2*y[10]*y[2]*y00im*y01im + 
                        4*c_ab*s2*lnv_rates*y[10]*y[2]*y00im*y01im - 4*c_aa*g2*y[8]*y[3]*y00re*y01im - 
                        4*c_aa*s2*lnv_rates*y[8]*y[3]*y00re*y01im - 4*c_ab*g2*y[9]*y[3]*y00re*y01im - 
                        4*c_ab*s2*lnv_rates*y[9]*y[3]*y00re*y01im - 4*c_ab*g2*y[10]*y[3]*y00re*y01im - 
                        4*c_ab*s2*lnv_rates*y[10]*y[3]*y00re*y01im - _KAV_*y[2]*y00re*y01im + 
                        4*c_aa*g2*y[8]*y[3]*y00im*y01re + 4*c_aa*s2*lnv_rates*y[8]*y[3]*y00im*y01re + 
                        4*c_ab*g2*y[9]*y[3]*y00im*y01re + 4*c_ab*s2*lnv_rates*y[9]*y[3]*y00im*y01re + 
                        4*c_ab*g2*y[10]*y[3]*y00im*y01re + 4*c_ab*s2*lnv_rates*y[10]*y[3]*y00im*y01re + 
                        _KAV_*y[2]*y00im*y01re + _KAV_*y[3]*y00re*y01re - 4*c_aa*g2*y[8]*y[2]*y00re*y01re + 
                        4*c_aa*s2*lnv_rates*y[8]*y[2]*y00re*y01re - 4*c_ab*g2*y[9]*y[2]*y00re*y01re + 
                        4*c_ab*s2*lnv_rates*y[9]*y[2]*y00re*y01re - 4*c_ab*g2*y[10]*y[2]*y00re*y01re + 
                        4*c_ab*s2*lnv_rates*y[10]*y[2]*y00re*y01re + 4*c_ab*g1_r_to_mu*g1*y[8]*pow(y10im,2.) - 
                        4*c_ab*s1*lnv_rates*y[8]*pow(y10im,2.) + 4*c_aa*g1_r_to_mu*g1*y[9]*pow(y10im,2.) - 
                        4*c_aa*s1*lnv_rates*y[9]*pow(y10im,2.) + 4*c_ab*g1_r_to_mu*g1*y[10]*pow(y10im,2.) - 
                        4*c_ab*s1*lnv_rates*y[10]*pow(y10im,2.) - 4*c_ab*g2*y[8]*y[0]*pow(y10im,2.) + 
                        4*c_ab*s2*lnv_rates*y[8]*y[0]*pow(y10im,2.) - 4*c_aa*g2*y[9]*y[0]*pow(y10im,2.) + 
                        4*c_aa*s2*lnv_rates*y[9]*y[0]*pow(y10im,2.) - 4*c_ab*g2*y[10]*y[0]*pow(y10im,2.) + 
                        4*c_ab*s2*lnv_rates*y[10]*y[0]*pow(y10im,2.) + 4*c_ab*g1_r_to_mu*g1*y[8]*pow(y10re,2.) - 
                        4*c_ab*s1*lnv_rates*y[8]*pow(y10re,2.) + 4*c_aa*g1_r_to_mu*g1*y[9]*pow(y10re,2.) - 
                        4*c_aa*s1*lnv_rates*y[9]*pow(y10re,2.) + 4*c_ab*g1_r_to_mu*g1*y[10]*pow(y10re,2.) - 
                        4*c_ab*s1*lnv_rates*y[10]*pow(y10re,2.) - 4*c_ab*g2*y[8]*y[0]*pow(y10re,2.) + 
                        4*c_ab*s2*lnv_rates*y[8]*y[0]*pow(y10re,2.) - 4*c_aa*g2*y[9]*y[0]*pow(y10re,2.) + 
                        4*c_aa*s2*lnv_rates*y[9]*y[0]*pow(y10re,2.) - 4*c_ab*g2*y[10]*y[0]*pow(y10re,2.) + 
                        4*c_ab*s2*lnv_rates*y[10]*y[0]*pow(y10re,2.) + _KAV_*y[3]*y10im*y11im - 
                        4*c_ab*g2*y[8]*y[2]*y10im*y11im + 4*c_ab*s2*lnv_rates*y[8]*y[2]*y10im*y11im - 
                        4*c_aa*g2*y[9]*y[2]*y10im*y11im + 4*c_aa*s2*lnv_rates*y[9]*y[2]*y10im*y11im - 
                        4*c_ab*g2*y[10]*y[2]*y10im*y11im + 4*c_ab*s2*lnv_rates*y[10]*y[2]*y10im*y11im - 
                        4*c_ab*g2*y[8]*y[3]*y10re*y11im - 4*c_ab*s2*lnv_rates*y[8]*y[3]*y10re*y11im - 
                        4*c_aa*g2*y[9]*y[3]*y10re*y11im - 4*c_aa*s2*lnv_rates*y[9]*y[3]*y10re*y11im - 
                        4*c_ab*g2*y[10]*y[3]*y10re*y11im - 4*c_ab*s2*lnv_rates*y[10]*y[3]*y10re*y11im - 
                        _KAV_*y[2]*y10re*y11im + 4*c_ab*g2*y[8]*y[3]*y10im*y11re + 
                        4*c_ab*s2*lnv_rates*y[8]*y[3]*y10im*y11re + 4*c_aa*g2*y[9]*y[3]*y10im*y11re + 
                        4*c_aa*s2*lnv_rates*y[9]*y[3]*y10im*y11re + 4*c_ab*g2*y[10]*y[3]*y10im*y11re + 
                        4*c_ab*s2*lnv_rates*y[10]*y[3]*y10im*y11re + _KAV_*y[2]*y10im*y11re + _KAV_*y[3]*y10re*y11re - 
                        4*c_ab*g2*y[8]*y[2]*y10re*y11re + 4*c_ab*s2*lnv_rates*y[8]*y[2]*y10re*y11re - 
                        4*c_aa*g2*y[9]*y[2]*y10re*y11re + 4*c_aa*s2*lnv_rates*y[9]*y[2]*y10re*y11re - 
                        4*c_ab*g2*y[10]*y[2]*y10re*y11re + 4*c_ab*s2*lnv_rates*y[10]*y[2]*y10re*y11re + 
                        4*c_ab*g1_r_to_mu*g1*y[8]*pow(y20im,2.) - 4*c_ab*s1*lnv_rates*y[8]*pow(y20im,2.) + 
                        4*c_ab*g1_r_to_mu*g1*y[9]*pow(y20im,2.) - 4*c_ab*s1*lnv_rates*y[9]*pow(y20im,2.) + 
                        4*c_aa*g1_r_to_mu*g1*y[10]*pow(y20im,2.) - 4*c_aa*s1*lnv_rates*y[10]*pow(y20im,2.) - 
                        4*c_ab*g2*y[8]*y[0]*pow(y20im,2.) + 4*c_ab*s2*lnv_rates*y[8]*y[0]*pow(y20im,2.) - 
                        4*c_ab*g2*y[9]*y[0]*pow(y20im,2.) + 4*c_ab*s2*lnv_rates*y[9]*y[0]*pow(y20im,2.) - 
                        4*c_aa*g2*y[10]*y[0]*pow(y20im,2.) + 4*c_aa*s2*lnv_rates*y[10]*y[0]*pow(y20im,2.) + 
                        4*c_ab*g1_r_to_mu*g1*y[8]*pow(y20re,2.) - 4*c_ab*s1*lnv_rates*y[8]*pow(y20re,2.) + 
                        4*c_ab*g1_r_to_mu*g1*y[9]*pow(y20re,2.) - 4*c_ab*s1*lnv_rates*y[9]*pow(y20re,2.) + 
                        4*c_aa*g1_r_to_mu*g1*y[10]*pow(y20re,2.) - 4*c_aa*s1*lnv_rates*y[10]*pow(y20re,2.) - 
                        4*c_ab*g2*y[8]*y[0]*pow(y20re,2.) + 4*c_ab*s2*lnv_rates*y[8]*y[0]*pow(y20re,2.) - 
                        4*c_ab*g2*y[9]*y[0]*pow(y20re,2.) + 4*c_ab*s2*lnv_rates*y[9]*y[0]*pow(y20re,2.) - 
                        4*c_aa*g2*y[10]*y[0]*pow(y20re,2.) + 4*c_aa*s2*lnv_rates*y[10]*y[0]*pow(y20re,2.) + 
                        _KAV_*y[3]*y20im*y21im - 4*c_ab*g2*y[8]*y[2]*y20im*y21im + 
                        4*c_ab*s2*lnv_rates*y[8]*y[2]*y20im*y21im - 4*c_ab*g2*y[9]*y[2]*y20im*y21im + 
                        4*c_ab*s2*lnv_rates*y[9]*y[2]*y20im*y21im - 4*c_aa*g2*y[10]*y[2]*y20im*y21im + 
                        4*c_aa*s2*lnv_rates*y[10]*y[2]*y20im*y21im - 4*c_ab*g2*y[8]*y[3]*y20re*y21im - 
                        4*c_ab*s2*lnv_rates*y[8]*y[3]*y20re*y21im - 4*c_ab*g2*y[9]*y[3]*y20re*y21im - 
                        4*c_ab*s2*lnv_rates*y[9]*y[3]*y20re*y21im - 4*c_aa*g2*y[10]*y[3]*y20re*y21im - 
                        4*c_aa*s2*lnv_rates*y[10]*y[3]*y20re*y21im - _KAV_*y[2]*y20re*y21im + 
                        (4*(g2 + s2*lnv_rates)*(c_ab*(y[8] + y[9]) + c_aa*y[10])*y[3]*y20im + _KAV_*y[2]*y20im + 
                        _KAV_*y[3]*y20re - 4*(g2 - s2*lnv_rates)*(c_ab*(y[8] + y[9]) + c_aa*y[10])*y[2]*y20re)*y21re + 
                        4*g0*((-1 + y[0])*pow(y00im,2.) - pow(y00re,2.) + y[3]*y00re*y01im + y[2]*y00re*y01re + 
                        y00im*(y[2]*y01im - y[3]*y01re) - pow(y10im,2.) - pow(y10re,2.) + y[2]*y10im*y11im + 
                        y[3]*y10re*y11im - y[3]*y10im*y11re + y[2]*y10re*y11re - pow(y20im,2.) - pow(y20re,2.) + 
                        y[0]*(pow(y00re,2.) + pow(y10im,2.) + pow(y10re,2.) + pow(y20im,2.) + pow(y20re,2.)) + 
                        y[2]*y20im*y21im + y[3]*y20re*y21im - y[3]*y20im*y21re + y[2]*y20re*y21re) + 
                        4*s0*lnv_rates*((-1 + y[0])*pow(y00im,2.) - pow(y00re,2.) - y[3]*y00re*y01im + y[2]*y00re*y01re + 
                        y00im*(y[2]*y01im + y[3]*y01re) - pow(y10im,2.) - pow(y10re,2.) + y[2]*y10im*y11im - 
                        y[3]*y10re*y11im + y[3]*y10im*y11re + y[2]*y10re*y11re - pow(y20im,2.) - pow(y20re,2.) + 
                        y[0]*(pow(y00re,2.) + pow(y10im,2.) + pow(y10re,2.) + pow(y20im,2.) + pow(y20re,2.)) + 
                        y[2]*y20im*y21im - y[3]*y20re*y21im + y[3]*y20im*y21re + y[2]*y20re*y21re)))    - //stef, new mthM term
                        2*_MPLTSP_*thmlnv*VNM*(y[2]*(y00re*y01im - y00im*y01re + y10re*y11im - y10im*y11re + y20re*y21im - 
                        y20im*y21re) + y[3]*(y00im*y01im + y00re*y01re + y10im*y11im + y10re*y11re + y20im*y21im + 
                        y20re*y21re));

        yaux[1]         = (_MPLTSP_*(_KAV_*(y[2]*(-(y00re*y01im) + y00im*y01re - y10re*y11im + y10im*y11re - y20re*y21im + 
                        y20im*y21re) + y[3]*(y00im*y01im + y00re*y01re + y10im*y11im + y10re*y11re + y20im*y21im + 
                        y20re*y21re)) - 4*(-(c_aa*g2*y[8]*y[2]*y00im*y01im) + 
                        c_aa*s2*lnv_rates*y[8]*y[2]*y00im*y01im - c_ab*g2*y[9]*y[2]*y00im*y01im + 
                        c_ab*s2*lnv_rates*y[9]*y[2]*y00im*y01im - c_ab*g2*y[10]*y[2]*y00im*y01im + 
                        c_ab*s2*lnv_rates*y[10]*y[2]*y00im*y01im - c_aa*g2*y[8]*y[3]*y00re*y01im - 
                        c_aa*s2*lnv_rates*y[8]*y[3]*y00re*y01im - c_ab*g2*y[9]*y[3]*y00re*y01im - 
                        c_ab*s2*lnv_rates*y[9]*y[3]*y00re*y01im - c_ab*g2*y[10]*y[3]*y00re*y01im - 
                        c_ab*s2*lnv_rates*y[10]*y[3]*y00re*y01im + c_aa*g1_r_to_mu*g1*y[8]*pow(y01im,2.) - 
                        c_aa*s1*lnv_rates*y[8]*pow(y01im,2.) + c_ab*g1_r_to_mu*g1*y[9]*pow(y01im,2.) - 
                        c_ab*s1*lnv_rates*y[9]*pow(y01im,2.) + c_ab*g1_r_to_mu*g1*y[10]*pow(y01im,2.) - 
                        c_ab*s1*lnv_rates*y[10]*pow(y01im,2.) - c_aa*g2*y[8]*y[1]*pow(y01im,2.) + 
                        c_aa*s2*lnv_rates*y[8]*y[1]*pow(y01im,2.) - c_ab*g2*y[9]*y[1]*pow(y01im,2.) + 
                        c_ab*s2*lnv_rates*y[9]*y[1]*pow(y01im,2.) - c_ab*g2*y[10]*y[1]*pow(y01im,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[1]*pow(y01im,2.) + c_aa*g2*y[8]*y[3]*y00im*y01re + 
                        c_aa*s2*lnv_rates*y[8]*y[3]*y00im*y01re + c_ab*g2*y[9]*y[3]*y00im*y01re + 
                        c_ab*s2*lnv_rates*y[9]*y[3]*y00im*y01re + c_ab*g2*y[10]*y[3]*y00im*y01re + 
                        c_ab*s2*lnv_rates*y[10]*y[3]*y00im*y01re - c_aa*g2*y[8]*y[2]*y00re*y01re + 
                        c_aa*s2*lnv_rates*y[8]*y[2]*y00re*y01re - c_ab*g2*y[9]*y[2]*y00re*y01re + 
                        c_ab*s2*lnv_rates*y[9]*y[2]*y00re*y01re - c_ab*g2*y[10]*y[2]*y00re*y01re + 
                        c_ab*s2*lnv_rates*y[10]*y[2]*y00re*y01re + c_aa*g1_r_to_mu*g1*y[8]*pow(y01re,2.) - 
                        c_aa*s1*lnv_rates*y[8]*pow(y01re,2.) + c_ab*g1_r_to_mu*g1*y[9]*pow(y01re,2.) - 
                        c_ab*s1*lnv_rates*y[9]*pow(y01re,2.) + c_ab*g1_r_to_mu*g1*y[10]*pow(y01re,2.) - 
                        c_ab*s1*lnv_rates*y[10]*pow(y01re,2.) - c_aa*g2*y[8]*y[1]*pow(y01re,2.) + 
                        c_aa*s2*lnv_rates*y[8]*y[1]*pow(y01re,2.) - c_ab*g2*y[9]*y[1]*pow(y01re,2.) + 
                        c_ab*s2*lnv_rates*y[9]*y[1]*pow(y01re,2.) - c_ab*g2*y[10]*y[1]*pow(y01re,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[1]*pow(y01re,2.) - c_ab*g2*y[8]*y[2]*y10im*y11im + 
                        c_ab*s2*lnv_rates*y[8]*y[2]*y10im*y11im - c_aa*g2*y[9]*y[2]*y10im*y11im + 
                        c_aa*s2*lnv_rates*y[9]*y[2]*y10im*y11im - c_ab*g2*y[10]*y[2]*y10im*y11im + 
                        c_ab*s2*lnv_rates*y[10]*y[2]*y10im*y11im - c_ab*g2*y[8]*y[3]*y10re*y11im - 
                        c_ab*s2*lnv_rates*y[8]*y[3]*y10re*y11im - c_aa*g2*y[9]*y[3]*y10re*y11im - 
                        c_aa*s2*lnv_rates*y[9]*y[3]*y10re*y11im - c_ab*g2*y[10]*y[3]*y10re*y11im - 
                        c_ab*s2*lnv_rates*y[10]*y[3]*y10re*y11im + c_ab*g1_r_to_mu*g1*y[8]*pow(y11im,2.) - 
                        c_ab*s1*lnv_rates*y[8]*pow(y11im,2.) + c_aa*g1_r_to_mu*g1*y[9]*pow(y11im,2.) - 
                        c_aa*s1*lnv_rates*y[9]*pow(y11im,2.) + c_ab*g1_r_to_mu*g1*y[10]*pow(y11im,2.) - 
                        c_ab*s1*lnv_rates*y[10]*pow(y11im,2.) - c_ab*g2*y[8]*y[1]*pow(y11im,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[1]*pow(y11im,2.) - c_aa*g2*y[9]*y[1]*pow(y11im,2.) + 
                        c_aa*s2*lnv_rates*y[9]*y[1]*pow(y11im,2.) - c_ab*g2*y[10]*y[1]*pow(y11im,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[1]*pow(y11im,2.) + c_ab*g2*y[8]*y[3]*y10im*y11re + 
                        c_ab*s2*lnv_rates*y[8]*y[3]*y10im*y11re + c_aa*g2*y[9]*y[3]*y10im*y11re + 
                        c_aa*s2*lnv_rates*y[9]*y[3]*y10im*y11re + c_ab*g2*y[10]*y[3]*y10im*y11re + 
                        c_ab*s2*lnv_rates*y[10]*y[3]*y10im*y11re - c_ab*g2*y[8]*y[2]*y10re*y11re + 
                        c_ab*s2*lnv_rates*y[8]*y[2]*y10re*y11re - c_aa*g2*y[9]*y[2]*y10re*y11re + 
                        c_aa*s2*lnv_rates*y[9]*y[2]*y10re*y11re - c_ab*g2*y[10]*y[2]*y10re*y11re + 
                        c_ab*s2*lnv_rates*y[10]*y[2]*y10re*y11re + c_ab*g1_r_to_mu*g1*y[8]*pow(y11re,2.) - 
                        c_ab*s1*lnv_rates*y[8]*pow(y11re,2.) + c_aa*g1_r_to_mu*g1*y[9]*pow(y11re,2.) - 
                        c_aa*s1*lnv_rates*y[9]*pow(y11re,2.) + c_ab*g1_r_to_mu*g1*y[10]*pow(y11re,2.) - 
                        c_ab*s1*lnv_rates*y[10]*pow(y11re,2.) - c_ab*g2*y[8]*y[1]*pow(y11re,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[1]*pow(y11re,2.) - c_aa*g2*y[9]*y[1]*pow(y11re,2.) + 
                        c_aa*s2*lnv_rates*y[9]*y[1]*pow(y11re,2.) - c_ab*g2*y[10]*y[1]*pow(y11re,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[1]*pow(y11re,2.) - c_ab*g2*y[8]*y[2]*y20im*y21im + 
                        c_ab*s2*lnv_rates*y[8]*y[2]*y20im*y21im - c_ab*g2*y[9]*y[2]*y20im*y21im + 
                        c_ab*s2*lnv_rates*y[9]*y[2]*y20im*y21im - c_aa*g2*y[10]*y[2]*y20im*y21im + 
                        c_aa*s2*lnv_rates*y[10]*y[2]*y20im*y21im - c_ab*g2*y[8]*y[3]*y20re*y21im - 
                        c_ab*s2*lnv_rates*y[8]*y[3]*y20re*y21im - c_ab*g2*y[9]*y[3]*y20re*y21im - 
                        c_ab*s2*lnv_rates*y[9]*y[3]*y20re*y21im - c_aa*g2*y[10]*y[3]*y20re*y21im - 
                        c_aa*s2*lnv_rates*y[10]*y[3]*y20re*y21im + c_ab*g1_r_to_mu*g1*y[8]*pow(y21im,2.) - 
                        c_ab*s1*lnv_rates*y[8]*pow(y21im,2.) + c_ab*g1_r_to_mu*g1*y[9]*pow(y21im,2.) - 
                        c_ab*s1*lnv_rates*y[9]*pow(y21im,2.) + c_aa*g1_r_to_mu*g1*y[10]*pow(y21im,2.) - 
                        c_aa*s1*lnv_rates*y[10]*pow(y21im,2.) - c_ab*g2*y[8]*y[1]*pow(y21im,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[1]*pow(y21im,2.) - c_ab*g2*y[9]*y[1]*pow(y21im,2.) + 
                        c_ab*s2*lnv_rates*y[9]*y[1]*pow(y21im,2.) - c_aa*g2*y[10]*y[1]*pow(y21im,2.) + 
                        c_aa*s2*lnv_rates*y[10]*y[1]*pow(y21im,2.) + 
                        (c_ab*(y[8] + y[9]) + c_aa*y[10])*
                        (g2*y[3]*y20im + s2*lnv_rates*y[3]*y20im - g2*y[2]*y20re + s2*lnv_rates*y[2]*y20re)*y21re + 
                        (c_ab*(y[8] + y[9]) + c_aa*y[10])*(g1_r_to_mu*g1 - s1*lnv_rates - g2*y[1] + s2*lnv_rates*y[1])*pow(y21re,2.) + 
                        g0*(y[3]*(y00re*y01im - y00im*y01re + y10re*y11im - y10im*y11re + y20re*y21im - 
                        y20im*y21re) + y[2]*(y00im*y01im + y00re*y01re + y10im*y11im + y10re*y11re + 
                        y20im*y21im + y20re*y21re) + 
                        (-1 + y[1])*(pow(y01im,2.) + pow(y01re,2.) + pow(y11im,2.) + pow(y11re,2.) + 
                        pow(y21im,2.) + pow(y21re,2.))) + 
                        s0*lnv_rates*(y[3]*(-(y00re*y01im) + y00im*y01re - y10re*y11im + y10im*y11re - y20re*y21im + 
                        y20im*y21re) + y[2]*(y00im*y01im + y00re*y01re + y10im*y11im + y10re*y11re + 
                        y20im*y21im + y20re*y21re) + 
                        (-1 + y[1])*(pow(y01im,2.) + pow(y01re,2.) + pow(y11im,2.) + pow(y11re,2.) + 
                        pow(y21im,2.) + pow(y21re,2.))))))/4. + //stef, new mthM term 
                        2*_MPLTSP_*thmlnv*VNM*(y[2]*(y00re*y01im - y00im*y01re + y10re*y11im - y10im*y11re + y20re*y21im - 
                        y20im*y21re) + y[3]*(y00im*y01im + y00re*y01re + y10im*y11im + y10re*y11re + y20im*y21im + 
                        y20re*y21re));

        yaux[2]         = -0.125*(_MPLTSP_*(_KAV_*(4*d12*y[3]*pow(t,2.) + 
                        pow(_TSP_,2.)*((y[0] - y[1])*(y00re*y01im - y00im*y01re + y10re*y11im - y10im*y11re + 
                        y20re*y21im - y20im*y21re) + 
                        y[3]*(-pow(y00im,2.) - pow(y00re,2.) + pow(y01im,2.) + pow(y01re,2.) - 
                        pow(y10im,2.) - pow(y10re,2.) + pow(y11im,2.) + pow(y11re,2.) - pow(y20im,2.) - 
                        pow(y20re,2.) + pow(y21im,2.) + pow(y21re,2.)))) + 
                        4*pow(_TSP_,2.)*(-(c_aa*g2*y[8]*y[2]*pow(y00im,2.)) + 
                        c_aa*s2*lnv_rates*y[8]*y[2]*pow(y00im,2.) - c_ab*g2*y[9]*y[2]*pow(y00im,2.) + 
                        c_ab*s2*lnv_rates*y[9]*y[2]*pow(y00im,2.) - c_ab*g2*y[10]*y[2]*pow(y00im,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[2]*pow(y00im,2.) - c_aa*g2*y[8]*y[2]*pow(y00re,2.) + 
                        c_aa*s2*lnv_rates*y[8]*y[2]*pow(y00re,2.) - c_ab*g2*y[9]*y[2]*pow(y00re,2.) + 
                        c_ab*s2*lnv_rates*y[9]*y[2]*pow(y00re,2.) - c_ab*g2*y[10]*y[2]*pow(y00re,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[2]*pow(y00re,2.) + 2*c_aa*g1_r_to_mu*g1*y[8]*y00im*y01im - 
                        2*c_aa*s1*lnv_rates*y[8]*y00im*y01im + 2*c_ab*g1_r_to_mu*g1*y[9]*y00im*y01im - 
                        2*c_ab*s1*lnv_rates*y[9]*y00im*y01im + 2*c_ab*g1_r_to_mu*g1*y[10]*y00im*y01im - 
                        2*c_ab*s1*lnv_rates*y[10]*y00im*y01im - c_aa*g2*y[8]*y[0]*y00im*y01im + 
                        c_aa*s2*lnv_rates*y[8]*y[0]*y00im*y01im - c_ab*g2*y[9]*y[0]*y00im*y01im + 
                        c_ab*s2*lnv_rates*y[9]*y[0]*y00im*y01im - c_ab*g2*y[10]*y[0]*y00im*y01im + 
                        c_ab*s2*lnv_rates*y[10]*y[0]*y00im*y01im - c_aa*g2*y[8]*y[1]*y00im*y01im + 
                        c_aa*s2*lnv_rates*y[8]*y[1]*y00im*y01im - c_ab*g2*y[9]*y[1]*y00im*y01im + 
                        c_ab*s2*lnv_rates*y[9]*y[1]*y00im*y01im - c_ab*g2*y[10]*y[1]*y00im*y01im + 
                        c_ab*s2*lnv_rates*y[10]*y[1]*y00im*y01im - c_aa*g2*y[8]*y[2]*pow(y01im,2.) + 
                        c_aa*s2*lnv_rates*y[8]*y[2]*pow(y01im,2.) - c_ab*g2*y[9]*y[2]*pow(y01im,2.) + 
                        c_ab*s2*lnv_rates*y[9]*y[2]*pow(y01im,2.) - c_ab*g2*y[10]*y[2]*pow(y01im,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[2]*pow(y01im,2.) + 2*c_aa*g1_r_to_mu*g1*y[8]*y00re*y01re - 
                        2*c_aa*s1*lnv_rates*y[8]*y00re*y01re + 2*c_ab*g1_r_to_mu*g1*y[9]*y00re*y01re - 
                        2*c_ab*s1*lnv_rates*y[9]*y00re*y01re + 2*c_ab*g1_r_to_mu*g1*y[10]*y00re*y01re - 
                        2*c_ab*s1*lnv_rates*y[10]*y00re*y01re - c_aa*g2*y[8]*y[0]*y00re*y01re + 
                        c_aa*s2*lnv_rates*y[8]*y[0]*y00re*y01re - c_ab*g2*y[9]*y[0]*y00re*y01re + 
                        c_ab*s2*lnv_rates*y[9]*y[0]*y00re*y01re - c_ab*g2*y[10]*y[0]*y00re*y01re + 
                        c_ab*s2*lnv_rates*y[10]*y[0]*y00re*y01re - c_aa*g2*y[8]*y[1]*y00re*y01re + 
                        c_aa*s2*lnv_rates*y[8]*y[1]*y00re*y01re - c_ab*g2*y[9]*y[1]*y00re*y01re + 
                        c_ab*s2*lnv_rates*y[9]*y[1]*y00re*y01re - c_ab*g2*y[10]*y[1]*y00re*y01re + 
                        c_ab*s2*lnv_rates*y[10]*y[1]*y00re*y01re - c_aa*g2*y[8]*y[2]*pow(y01re,2.) + 
                        c_aa*s2*lnv_rates*y[8]*y[2]*pow(y01re,2.) - c_ab*g2*y[9]*y[2]*pow(y01re,2.) + 
                        c_ab*s2*lnv_rates*y[9]*y[2]*pow(y01re,2.) - c_ab*g2*y[10]*y[2]*pow(y01re,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[2]*pow(y01re,2.) - c_ab*g2*y[8]*y[2]*pow(y10im,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[2]*pow(y10im,2.) - c_aa*g2*y[9]*y[2]*pow(y10im,2.) + 
                        c_aa*s2*lnv_rates*y[9]*y[2]*pow(y10im,2.) - c_ab*g2*y[10]*y[2]*pow(y10im,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[2]*pow(y10im,2.) - c_ab*g2*y[8]*y[2]*pow(y10re,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[2]*pow(y10re,2.) - c_aa*g2*y[9]*y[2]*pow(y10re,2.) + 
                        c_aa*s2*lnv_rates*y[9]*y[2]*pow(y10re,2.) - c_ab*g2*y[10]*y[2]*pow(y10re,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[2]*pow(y10re,2.) + 2*c_ab*g1_r_to_mu*g1*y[8]*y10im*y11im - 
                        2*c_ab*s1*lnv_rates*y[8]*y10im*y11im + 2*c_aa*g1_r_to_mu*g1*y[9]*y10im*y11im - 
                        2*c_aa*s1*lnv_rates*y[9]*y10im*y11im + 2*c_ab*g1_r_to_mu*g1*y[10]*y10im*y11im - 
                        2*c_ab*s1*lnv_rates*y[10]*y10im*y11im - c_ab*g2*y[8]*y[0]*y10im*y11im + 
                        c_ab*s2*lnv_rates*y[8]*y[0]*y10im*y11im - c_aa*g2*y[9]*y[0]*y10im*y11im + 
                        c_aa*s2*lnv_rates*y[9]*y[0]*y10im*y11im - c_ab*g2*y[10]*y[0]*y10im*y11im + 
                        c_ab*s2*lnv_rates*y[10]*y[0]*y10im*y11im - c_ab*g2*y[8]*y[1]*y10im*y11im + 
                        c_ab*s2*lnv_rates*y[8]*y[1]*y10im*y11im - c_aa*g2*y[9]*y[1]*y10im*y11im + 
                        c_aa*s2*lnv_rates*y[9]*y[1]*y10im*y11im - c_ab*g2*y[10]*y[1]*y10im*y11im + 
                        c_ab*s2*lnv_rates*y[10]*y[1]*y10im*y11im - c_ab*g2*y[8]*y[2]*pow(y11im,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[2]*pow(y11im,2.) - c_aa*g2*y[9]*y[2]*pow(y11im,2.) + 
                        c_aa*s2*lnv_rates*y[9]*y[2]*pow(y11im,2.) - c_ab*g2*y[10]*y[2]*pow(y11im,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[2]*pow(y11im,2.) + 2*c_ab*g1_r_to_mu*g1*y[8]*y10re*y11re - 
                        2*c_ab*s1*lnv_rates*y[8]*y10re*y11re + 2*c_aa*g1_r_to_mu*g1*y[9]*y10re*y11re - 
                        2*c_aa*s1*lnv_rates*y[9]*y10re*y11re + 2*c_ab*g1_r_to_mu*g1*y[10]*y10re*y11re - 
                        2*c_ab*s1*lnv_rates*y[10]*y10re*y11re - c_ab*g2*y[8]*y[0]*y10re*y11re + 
                        c_ab*s2*lnv_rates*y[8]*y[0]*y10re*y11re - c_aa*g2*y[9]*y[0]*y10re*y11re + 
                        c_aa*s2*lnv_rates*y[9]*y[0]*y10re*y11re - c_ab*g2*y[10]*y[0]*y10re*y11re + 
                        c_ab*s2*lnv_rates*y[10]*y[0]*y10re*y11re - c_ab*g2*y[8]*y[1]*y10re*y11re + 
                        c_ab*s2*lnv_rates*y[8]*y[1]*y10re*y11re - c_aa*g2*y[9]*y[1]*y10re*y11re + 
                        c_aa*s2*lnv_rates*y[9]*y[1]*y10re*y11re - c_ab*g2*y[10]*y[1]*y10re*y11re + 
                        c_ab*s2*lnv_rates*y[10]*y[1]*y10re*y11re - c_ab*g2*y[8]*y[2]*pow(y11re,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[2]*pow(y11re,2.) - c_aa*g2*y[9]*y[2]*pow(y11re,2.) + 
                        c_aa*s2*lnv_rates*y[9]*y[2]*pow(y11re,2.) - c_ab*g2*y[10]*y[2]*pow(y11re,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[2]*pow(y11re,2.) - c_ab*g2*y[8]*y[2]*pow(y20im,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[2]*pow(y20im,2.) - c_ab*g2*y[9]*y[2]*pow(y20im,2.) + 
                        c_ab*s2*lnv_rates*y[9]*y[2]*pow(y20im,2.) - c_aa*g2*y[10]*y[2]*pow(y20im,2.) + 
                        c_aa*s2*lnv_rates*y[10]*y[2]*pow(y20im,2.) - c_ab*g2*y[8]*y[2]*pow(y20re,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[2]*pow(y20re,2.) - c_ab*g2*y[9]*y[2]*pow(y20re,2.) + 
                        c_ab*s2*lnv_rates*y[9]*y[2]*pow(y20re,2.) - c_aa*g2*y[10]*y[2]*pow(y20re,2.) + 
                        c_aa*s2*lnv_rates*y[10]*y[2]*pow(y20re,2.) + 2*c_ab*g1_r_to_mu*g1*y[8]*y20im*y21im - 
                        2*c_ab*s1*lnv_rates*y[8]*y20im*y21im + 2*c_ab*g1_r_to_mu*g1*y[9]*y20im*y21im - 
                        2*c_ab*s1*lnv_rates*y[9]*y20im*y21im + 2*c_aa*g1_r_to_mu*g1*y[10]*y20im*y21im - 
                        2*c_aa*s1*lnv_rates*y[10]*y20im*y21im - c_ab*g2*y[8]*y[0]*y20im*y21im + 
                        c_ab*s2*lnv_rates*y[8]*y[0]*y20im*y21im - c_ab*g2*y[9]*y[0]*y20im*y21im + 
                        c_ab*s2*lnv_rates*y[9]*y[0]*y20im*y21im - c_aa*g2*y[10]*y[0]*y20im*y21im + 
                        c_aa*s2*lnv_rates*y[10]*y[0]*y20im*y21im - c_ab*g2*y[8]*y[1]*y20im*y21im + 
                        c_ab*s2*lnv_rates*y[8]*y[1]*y20im*y21im - c_ab*g2*y[9]*y[1]*y20im*y21im + 
                        c_ab*s2*lnv_rates*y[9]*y[1]*y20im*y21im - c_aa*g2*y[10]*y[1]*y20im*y21im + 
                        c_aa*s2*lnv_rates*y[10]*y[1]*y20im*y21im - c_ab*g2*y[8]*y[2]*pow(y21im,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[2]*pow(y21im,2.) - c_ab*g2*y[9]*y[2]*pow(y21im,2.) + 
                        c_ab*s2*lnv_rates*y[9]*y[2]*pow(y21im,2.) - c_aa*g2*y[10]*y[2]*pow(y21im,2.) + 
                        c_aa*s2*lnv_rates*y[10]*y[2]*pow(y21im,2.) + 
                        (c_ab*(y[8] + y[9]) + c_aa*y[10])*(2*g1_r_to_mu*g1 - 2*s1*lnv_rates - g2*(y[0] + y[1]) + s2*lnv_rates*(y[0] + y[1]))*
                        y20re*y21re - (g2 - s2*lnv_rates)*(c_ab*(y[8] + y[9]) + c_aa*y[10])*y[2]*pow(y21re,2.) + 
                        g0*((-2 + y[0] + y[1])*(y00im*y01im + y00re*y01re + y10im*y11im + y10re*y11re + y20im*y21im + 
                        y20re*y21re) + y[2]*(pow(y00im,2.) + pow(y00re,2.) + pow(y01im,2.) + 
                        pow(y01re,2.) + pow(y10im,2.) + pow(y10re,2.) + pow(y11im,2.) + pow(y11re,2.) + 
                        pow(y20im,2.) + pow(y20re,2.) + pow(y21im,2.) + pow(y21re,2.))) + 
                        s0*lnv_rates*((-2 + y[0] + y[1])*(y00im*y01im + y00re*y01re + y10im*y11im + y10re*y11re + 
                        y20im*y21im + y20re*y21re) + 
                        y[2]*(pow(y00im,2.) + pow(y00re,2.) + pow(y01im,2.) + pow(y01re,2.) + 
                        pow(y10im,2.) + pow(y10re,2.) + pow(y11im,2.) + pow(y11re,2.) + pow(y20im,2.) + 
                        pow(y20re,2.) + pow(y21im,2.) + pow(y21re,2.))))))/pow(_TSP_,2.) + //stef, new mthM term
                        _MPLTSP_*thmlnv*VNM*((y[0] - y[1])*(y00re*y01im - y00im*y01re + y10re*y11im - y10im*y11re + y20re*y21im - 
                        y20im*y21re) + y[3]*(pow(y00im,2.) + pow(y00re,2.) - pow(y01im,2.) - pow(y01re,2.) + 
                        pow(y10im,2.) + pow(y10re,2.) - pow(y11im,2.) - pow(y11re,2.) + pow(y20im,2.) + 
                        pow(y20re,2.) - pow(y21im,2.) - pow(y21re,2.)));

        yaux[3]         = (_MPLTSP_*(-4*(g0 + s0*lnv_rates - (g2 - s2*lnv_rates)*(c_aa*y[8] + c_ab*(y[9] + y[10])))*y[3]*
                        (pow(y00im,2.) + pow(y00re,2.)) + 
                        _KAV_*((4*d12*y[2]*pow(t,2.))/pow(_TSP_,2.) + 
                        (y[0] - y[1])*(y00im*y01im + y00re*y01re + y10im*y11im + y10re*y11re + y20im*y21im + 
                        y20re*y21re) + y[2]*(-pow(y00im,2.) - pow(y00re,2.) + pow(y01im,2.) + pow(y01re,2.) - 
                        pow(y10im,2.) - pow(y10re,2.) + pow(y11im,2.) + pow(y11re,2.) - pow(y20im,2.) - 
                        pow(y20re,2.) + pow(y21im,2.) + pow(y21re,2.))) - 
                        4*(2*c_aa*g1_r_to_mu*g1*y[8]*y00re*y01im + 2*c_aa*s1*lnv_rates*y[8]*y00re*y01im + 
                        2*c_ab*g1_r_to_mu*g1*y[9]*y00re*y01im + 2*c_ab*s1*lnv_rates*y[9]*y00re*y01im + 2*c_ab*g1_r_to_mu*g1*y[10]*y00re*y01im + 
                        2*c_ab*s1*lnv_rates*y[10]*y00re*y01im - c_aa*g2*y[8]*y[0]*y00re*y01im - 
                        c_aa*s2*lnv_rates*y[8]*y[0]*y00re*y01im - c_ab*g2*y[9]*y[0]*y00re*y01im - 
                        c_ab*s2*lnv_rates*y[9]*y[0]*y00re*y01im - c_ab*g2*y[10]*y[0]*y00re*y01im - 
                        c_ab*s2*lnv_rates*y[10]*y[0]*y00re*y01im - c_aa*g2*y[8]*y[1]*y00re*y01im - 
                        c_aa*s2*lnv_rates*y[8]*y[1]*y00re*y01im - c_ab*g2*y[9]*y[1]*y00re*y01im - 
                        c_ab*s2*lnv_rates*y[9]*y[1]*y00re*y01im - c_ab*g2*y[10]*y[1]*y00re*y01im - 
                        c_ab*s2*lnv_rates*y[10]*y[1]*y00re*y01im - c_aa*g2*y[8]*y[3]*pow(y01im,2.) + 
                        c_aa*s2*lnv_rates*y[8]*y[3]*pow(y01im,2.) - c_ab*g2*y[9]*y[3]*pow(y01im,2.) + 
                        c_ab*s2*lnv_rates*y[9]*y[3]*pow(y01im,2.) - c_ab*g2*y[10]*y[3]*pow(y01im,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[3]*pow(y01im,2.) - 2*c_aa*g1_r_to_mu*g1*y[8]*y00im*y01re - 
                        2*c_aa*s1*lnv_rates*y[8]*y00im*y01re - 2*c_ab*g1_r_to_mu*g1*y[9]*y00im*y01re - 
                        2*c_ab*s1*lnv_rates*y[9]*y00im*y01re - 2*c_ab*g1_r_to_mu*g1*y[10]*y00im*y01re - 
                        2*c_ab*s1*lnv_rates*y[10]*y00im*y01re + c_aa*g2*y[8]*y[0]*y00im*y01re + 
                        c_aa*s2*lnv_rates*y[8]*y[0]*y00im*y01re + c_ab*g2*y[9]*y[0]*y00im*y01re + 
                        c_ab*s2*lnv_rates*y[9]*y[0]*y00im*y01re + c_ab*g2*y[10]*y[0]*y00im*y01re + 
                        c_ab*s2*lnv_rates*y[10]*y[0]*y00im*y01re + c_aa*g2*y[8]*y[1]*y00im*y01re + 
                        c_aa*s2*lnv_rates*y[8]*y[1]*y00im*y01re + c_ab*g2*y[9]*y[1]*y00im*y01re + 
                        c_ab*s2*lnv_rates*y[9]*y[1]*y00im*y01re + c_ab*g2*y[10]*y[1]*y00im*y01re + 
                        c_ab*s2*lnv_rates*y[10]*y[1]*y00im*y01re - c_aa*g2*y[8]*y[3]*pow(y01re,2.) + 
                        c_aa*s2*lnv_rates*y[8]*y[3]*pow(y01re,2.) - c_ab*g2*y[9]*y[3]*pow(y01re,2.) + 
                        c_ab*s2*lnv_rates*y[9]*y[3]*pow(y01re,2.) - c_ab*g2*y[10]*y[3]*pow(y01re,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[3]*pow(y01re,2.) - c_ab*g2*y[8]*y[3]*pow(y10im,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[3]*pow(y10im,2.) - c_aa*g2*y[9]*y[3]*pow(y10im,2.) + 
                        c_aa*s2*lnv_rates*y[9]*y[3]*pow(y10im,2.) - c_ab*g2*y[10]*y[3]*pow(y10im,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[3]*pow(y10im,2.) - c_ab*g2*y[8]*y[3]*pow(y10re,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[3]*pow(y10re,2.) - c_aa*g2*y[9]*y[3]*pow(y10re,2.) + 
                        c_aa*s2*lnv_rates*y[9]*y[3]*pow(y10re,2.) - c_ab*g2*y[10]*y[3]*pow(y10re,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[3]*pow(y10re,2.) + 2*c_ab*g1_r_to_mu*g1*y[8]*y10re*y11im + 
                        2*c_ab*s1*lnv_rates*y[8]*y10re*y11im + 2*c_aa*g1_r_to_mu*g1*y[9]*y10re*y11im + 
                        2*c_aa*s1*lnv_rates*y[9]*y10re*y11im + 2*c_ab*g1_r_to_mu*g1*y[10]*y10re*y11im + 
                        2*c_ab*s1*lnv_rates*y[10]*y10re*y11im - c_ab*g2*y[8]*y[0]*y10re*y11im - 
                        c_ab*s2*lnv_rates*y[8]*y[0]*y10re*y11im - c_aa*g2*y[9]*y[0]*y10re*y11im - 
                        c_aa*s2*lnv_rates*y[9]*y[0]*y10re*y11im - c_ab*g2*y[10]*y[0]*y10re*y11im - 
                        c_ab*s2*lnv_rates*y[10]*y[0]*y10re*y11im - c_ab*g2*y[8]*y[1]*y10re*y11im - 
                        c_ab*s2*lnv_rates*y[8]*y[1]*y10re*y11im - c_aa*g2*y[9]*y[1]*y10re*y11im - 
                        c_aa*s2*lnv_rates*y[9]*y[1]*y10re*y11im - c_ab*g2*y[10]*y[1]*y10re*y11im - 
                        c_ab*s2*lnv_rates*y[10]*y[1]*y10re*y11im - c_ab*g2*y[8]*y[3]*pow(y11im,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[3]*pow(y11im,2.) - c_aa*g2*y[9]*y[3]*pow(y11im,2.) + 
                        c_aa*s2*lnv_rates*y[9]*y[3]*pow(y11im,2.) - c_ab*g2*y[10]*y[3]*pow(y11im,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[3]*pow(y11im,2.) - 2*c_ab*g1_r_to_mu*g1*y[8]*y10im*y11re - 
                        2*c_ab*s1*lnv_rates*y[8]*y10im*y11re - 2*c_aa*g1_r_to_mu*g1*y[9]*y10im*y11re - 
                        2*c_aa*s1*lnv_rates*y[9]*y10im*y11re - 2*c_ab*g1_r_to_mu*g1*y[10]*y10im*y11re - 
                        2*c_ab*s1*lnv_rates*y[10]*y10im*y11re + c_ab*g2*y[8]*y[0]*y10im*y11re + 
                        c_ab*s2*lnv_rates*y[8]*y[0]*y10im*y11re + c_aa*g2*y[9]*y[0]*y10im*y11re + 
                        c_aa*s2*lnv_rates*y[9]*y[0]*y10im*y11re + c_ab*g2*y[10]*y[0]*y10im*y11re + 
                        c_ab*s2*lnv_rates*y[10]*y[0]*y10im*y11re + c_ab*g2*y[8]*y[1]*y10im*y11re + 
                        c_ab*s2*lnv_rates*y[8]*y[1]*y10im*y11re + c_aa*g2*y[9]*y[1]*y10im*y11re + 
                        c_aa*s2*lnv_rates*y[9]*y[1]*y10im*y11re + c_ab*g2*y[10]*y[1]*y10im*y11re + 
                        c_ab*s2*lnv_rates*y[10]*y[1]*y10im*y11re - c_ab*g2*y[8]*y[3]*pow(y11re,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[3]*pow(y11re,2.) - c_aa*g2*y[9]*y[3]*pow(y11re,2.) + 
                        c_aa*s2*lnv_rates*y[9]*y[3]*pow(y11re,2.) - c_ab*g2*y[10]*y[3]*pow(y11re,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[3]*pow(y11re,2.) - c_ab*g2*y[8]*y[3]*pow(y20im,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[3]*pow(y20im,2.) - c_ab*g2*y[9]*y[3]*pow(y20im,2.) + 
                        c_ab*s2*lnv_rates*y[9]*y[3]*pow(y20im,2.) - c_aa*g2*y[10]*y[3]*pow(y20im,2.) + 
                        c_aa*s2*lnv_rates*y[10]*y[3]*pow(y20im,2.) - c_ab*g2*y[8]*y[3]*pow(y20re,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[3]*pow(y20re,2.) - c_ab*g2*y[9]*y[3]*pow(y20re,2.) + 
                        c_ab*s2*lnv_rates*y[9]*y[3]*pow(y20re,2.) - c_aa*g2*y[10]*y[3]*pow(y20re,2.) + 
                        c_aa*s2*lnv_rates*y[10]*y[3]*pow(y20re,2.) + 2*c_ab*g1_r_to_mu*g1*y[8]*y20re*y21im + 
                        2*c_ab*s1*lnv_rates*y[8]*y20re*y21im + 2*c_ab*g1_r_to_mu*g1*y[9]*y20re*y21im + 
                        2*c_ab*s1*lnv_rates*y[9]*y20re*y21im + 2*c_aa*g1_r_to_mu*g1*y[10]*y20re*y21im + 
                        2*c_aa*s1*lnv_rates*y[10]*y20re*y21im - c_ab*g2*y[8]*y[0]*y20re*y21im - 
                        c_ab*s2*lnv_rates*y[8]*y[0]*y20re*y21im - c_ab*g2*y[9]*y[0]*y20re*y21im - 
                        c_ab*s2*lnv_rates*y[9]*y[0]*y20re*y21im - c_aa*g2*y[10]*y[0]*y20re*y21im - 
                        c_aa*s2*lnv_rates*y[10]*y[0]*y20re*y21im - c_ab*g2*y[8]*y[1]*y20re*y21im - 
                        c_ab*s2*lnv_rates*y[8]*y[1]*y20re*y21im - c_ab*g2*y[9]*y[1]*y20re*y21im - 
                        c_ab*s2*lnv_rates*y[9]*y[1]*y20re*y21im - c_aa*g2*y[10]*y[1]*y20re*y21im - 
                        c_aa*s2*lnv_rates*y[10]*y[1]*y20re*y21im - c_ab*g2*y[8]*y[3]*pow(y21im,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[3]*pow(y21im,2.) - c_ab*g2*y[9]*y[3]*pow(y21im,2.) + 
                        c_ab*s2*lnv_rates*y[9]*y[3]*pow(y21im,2.) - c_aa*g2*y[10]*y[3]*pow(y21im,2.) + 
                        c_aa*s2*lnv_rates*y[10]*y[3]*pow(y21im,2.) - 
                        (c_ab*(y[8] + y[9]) + c_aa*y[10])*(2*g1_r_to_mu*g1 + 2*s1*lnv_rates - (g2 + s2*lnv_rates)*(y[0] + y[1]))*y20im*y21re - 
                        (g2 - s2*lnv_rates)*(c_ab*(y[8] + y[9]) + c_aa*y[10])*y[3]*pow(y21re,2.) + 
                        g0*((-2 + y[0] + y[1])*y00re*y01im - 
                        (-2 + y[0] + y[1])*(y00im*y01re - y10re*y11im + y10im*y11re - y20re*y21im + y20im*y21re) + 
                        y[3]*(pow(y01im,2.) + pow(y01re,2.) + pow(y10im,2.) + pow(y10re,2.) + pow(y11im,2.) + 
                        pow(y11re,2.) + pow(y20im,2.) + pow(y20re,2.) + pow(y21im,2.) + pow(y21re,2.))) + 
                        s0*lnv_rates*(-((-2 + y[0] + y[1])*y00re*y01im) + 
                        (-2 + y[0] + y[1])*(y00im*y01re - y10re*y11im + y10im*y11re - y20re*y21im + y20im*y21re) + 
                        y[3]*(pow(y01im,2.) + pow(y01re,2.) + pow(y10im,2.) + pow(y10re,2.) + pow(y11im,2.) + 
                        pow(y11re,2.) + pow(y20im,2.) + pow(y20re,2.) + pow(y21im,2.) + pow(y21re,2.))))))/8. + //stef, new mthM term
                        _MPLTSP_*thmlnv*VNM*((y[0] - y[1])*(y00im*y01im + y00re*y01re + y10im*y11im + y10re*y11re + y20im*y21im + 
                        y20re*y21re) + y[2]*(-pow(y00im,2.) - pow(y00re,2.) + pow(y01im,2.) + pow(y01re,2.) - 
                        pow(y10im,2.) - pow(y10re,2.) + pow(y11im,2.) + pow(y11re,2.) - pow(y20im,2.) - 
                        pow(y20re,2.) + pow(y21im,2.) + pow(y21re,2.)));


        yaux[4]         = -0.25*(_MPLTSP_*(-4*c_aa*g1_r_to_mu*g1*y[8]*pow(y00im,2.) + 4*c_aa*s1*lnv_rates*y[8]*pow(y00im,2.) - 
                        4*c_ab*g1_r_to_mu*g1*y[9]*pow(y00im,2.) + 4*c_ab*s1*lnv_rates*y[9]*pow(y00im,2.) - 
                        4*c_ab*g1_r_to_mu*g1*y[10]*pow(y00im,2.) + 4*c_ab*s1*lnv_rates*y[10]*pow(y00im,2.) + 
                        4*c_aa*g2*y[8]*y[4]*pow(y00im,2.) - 4*c_aa*s2*lnv_rates*y[8]*y[4]*pow(y00im,2.) + 
                        4*c_ab*g2*y[9]*y[4]*pow(y00im,2.) - 4*c_ab*s2*lnv_rates*y[9]*y[4]*pow(y00im,2.) + 
                        4*c_ab*g2*y[10]*y[4]*pow(y00im,2.) - 4*c_ab*s2*lnv_rates*y[10]*y[4]*pow(y00im,2.) - 
                        4*c_aa*g1_r_to_mu*g1*y[8]*pow(y00re,2.) + 4*c_aa*s1*lnv_rates*y[8]*pow(y00re,2.) - 
                        4*c_ab*g1_r_to_mu*g1*y[9]*pow(y00re,2.) + 4*c_ab*s1*lnv_rates*y[9]*pow(y00re,2.) - 
                        4*c_ab*g1_r_to_mu*g1*y[10]*pow(y00re,2.) + 4*c_ab*s1*lnv_rates*y[10]*pow(y00re,2.) + 
                        4*c_aa*g2*y[8]*y[4]*pow(y00re,2.) - 4*c_aa*s2*lnv_rates*y[8]*y[4]*pow(y00re,2.) + 
                        4*c_ab*g2*y[9]*y[4]*pow(y00re,2.) - 4*c_ab*s2*lnv_rates*y[9]*y[4]*pow(y00re,2.) + 
                        4*c_ab*g2*y[10]*y[4]*pow(y00re,2.) - 4*c_ab*s2*lnv_rates*y[10]*y[4]*pow(y00re,2.) + 
                        _KAV_*y[7]*y00im*y01im + 4*c_aa*g2*y[8]*y[6]*y00im*y01im - 
                        4*c_aa*s2*lnv_rates*y[8]*y[6]*y00im*y01im + 4*c_ab*g2*y[9]*y[6]*y00im*y01im - 
                        4*c_ab*s2*lnv_rates*y[9]*y[6]*y00im*y01im + 4*c_ab*g2*y[10]*y[6]*y00im*y01im - 
                        4*c_ab*s2*lnv_rates*y[10]*y[6]*y00im*y01im - 4*c_aa*g2*y[8]*y[7]*y00re*y01im - 
                        4*c_aa*s2*lnv_rates*y[8]*y[7]*y00re*y01im - 4*c_ab*g2*y[9]*y[7]*y00re*y01im - 
                        4*c_ab*s2*lnv_rates*y[9]*y[7]*y00re*y01im - 4*c_ab*g2*y[10]*y[7]*y00re*y01im - 
                        4*c_ab*s2*lnv_rates*y[10]*y[7]*y00re*y01im + _KAV_*y[6]*y00re*y01im + 
                        4*c_aa*g2*y[8]*y[7]*y00im*y01re + 4*c_aa*s2*lnv_rates*y[8]*y[7]*y00im*y01re + 
                        4*c_ab*g2*y[9]*y[7]*y00im*y01re + 4*c_ab*s2*lnv_rates*y[9]*y[7]*y00im*y01re + 
                        4*c_ab*g2*y[10]*y[7]*y00im*y01re + 4*c_ab*s2*lnv_rates*y[10]*y[7]*y00im*y01re - 
                        _KAV_*y[6]*y00im*y01re + _KAV_*y[7]*y00re*y01re + 
                        4*c_aa*g2*y[8]*y[6]*y00re*y01re - 4*c_aa*s2*lnv_rates*y[8]*y[6]*y00re*y01re + 
                        4*c_ab*g2*y[9]*y[6]*y00re*y01re - 4*c_ab*s2*lnv_rates*y[9]*y[6]*y00re*y01re + 
                        4*c_ab*g2*y[10]*y[6]*y00re*y01re - 4*c_ab*s2*lnv_rates*y[10]*y[6]*y00re*y01re - 
                        4*c_ab*g1_r_to_mu*g1*y[8]*pow(y10im,2.) + 4*c_ab*s1*lnv_rates*y[8]*pow(y10im,2.) - 
                        4*c_aa*g1_r_to_mu*g1*y[9]*pow(y10im,2.) + 4*c_aa*s1*lnv_rates*y[9]*pow(y10im,2.) - 
                        4*c_ab*g1_r_to_mu*g1*y[10]*pow(y10im,2.) + 4*c_ab*s1*lnv_rates*y[10]*pow(y10im,2.) + 
                        4*c_ab*g2*y[8]*y[4]*pow(y10im,2.) - 4*c_ab*s2*lnv_rates*y[8]*y[4]*pow(y10im,2.) + 
                        4*c_aa*g2*y[9]*y[4]*pow(y10im,2.) - 4*c_aa*s2*lnv_rates*y[9]*y[4]*pow(y10im,2.) + 
                        4*c_ab*g2*y[10]*y[4]*pow(y10im,2.) - 4*c_ab*s2*lnv_rates*y[10]*y[4]*pow(y10im,2.) - 
                        4*c_ab*g1_r_to_mu*g1*y[8]*pow(y10re,2.) + 4*c_ab*s1*lnv_rates*y[8]*pow(y10re,2.) - 
                        4*c_aa*g1_r_to_mu*g1*y[9]*pow(y10re,2.) + 4*c_aa*s1*lnv_rates*y[9]*pow(y10re,2.) - 
                        4*c_ab*g1_r_to_mu*g1*y[10]*pow(y10re,2.) + 4*c_ab*s1*lnv_rates*y[10]*pow(y10re,2.) + 
                        4*c_ab*g2*y[8]*y[4]*pow(y10re,2.) - 4*c_ab*s2*lnv_rates*y[8]*y[4]*pow(y10re,2.) + 
                        4*c_aa*g2*y[9]*y[4]*pow(y10re,2.) - 4*c_aa*s2*lnv_rates*y[9]*y[4]*pow(y10re,2.) + 
                        4*c_ab*g2*y[10]*y[4]*pow(y10re,2.) - 4*c_ab*s2*lnv_rates*y[10]*y[4]*pow(y10re,2.) + 
                        _KAV_*y[7]*y10im*y11im + 4*c_ab*g2*y[8]*y[6]*y10im*y11im - 
                        4*c_ab*s2*lnv_rates*y[8]*y[6]*y10im*y11im + 4*c_aa*g2*y[9]*y[6]*y10im*y11im - 
                        4*c_aa*s2*lnv_rates*y[9]*y[6]*y10im*y11im + 4*c_ab*g2*y[10]*y[6]*y10im*y11im - 
                        4*c_ab*s2*lnv_rates*y[10]*y[6]*y10im*y11im - 4*c_ab*g2*y[8]*y[7]*y10re*y11im - 
                        4*c_ab*s2*lnv_rates*y[8]*y[7]*y10re*y11im - 4*c_aa*g2*y[9]*y[7]*y10re*y11im - 
                        4*c_aa*s2*lnv_rates*y[9]*y[7]*y10re*y11im - 4*c_ab*g2*y[10]*y[7]*y10re*y11im - 
                        4*c_ab*s2*lnv_rates*y[10]*y[7]*y10re*y11im + _KAV_*y[6]*y10re*y11im + 
                        4*c_ab*g2*y[8]*y[7]*y10im*y11re + 4*c_ab*s2*lnv_rates*y[8]*y[7]*y10im*y11re + 
                        4*c_aa*g2*y[9]*y[7]*y10im*y11re + 4*c_aa*s2*lnv_rates*y[9]*y[7]*y10im*y11re + 
                        4*c_ab*g2*y[10]*y[7]*y10im*y11re + 4*c_ab*s2*lnv_rates*y[10]*y[7]*y10im*y11re - 
                        _KAV_*y[6]*y10im*y11re + _KAV_*y[7]*y10re*y11re + 
                        4*c_ab*g2*y[8]*y[6]*y10re*y11re - 4*c_ab*s2*lnv_rates*y[8]*y[6]*y10re*y11re + 
                        4*c_aa*g2*y[9]*y[6]*y10re*y11re - 4*c_aa*s2*lnv_rates*y[9]*y[6]*y10re*y11re + 
                        4*c_ab*g2*y[10]*y[6]*y10re*y11re - 4*c_ab*s2*lnv_rates*y[10]*y[6]*y10re*y11re - 
                        4*c_ab*g1_r_to_mu*g1*y[8]*pow(y20im,2.) + 4*c_ab*s1*lnv_rates*y[8]*pow(y20im,2.) - 
                        4*c_ab*g1_r_to_mu*g1*y[9]*pow(y20im,2.) + 4*c_ab*s1*lnv_rates*y[9]*pow(y20im,2.) - 
                        4*c_aa*g1_r_to_mu*g1*y[10]*pow(y20im,2.) + 4*c_aa*s1*lnv_rates*y[10]*pow(y20im,2.) + 
                        4*c_ab*g2*y[8]*y[4]*pow(y20im,2.) - 4*c_ab*s2*lnv_rates*y[8]*y[4]*pow(y20im,2.) + 
                        4*c_ab*g2*y[9]*y[4]*pow(y20im,2.) - 4*c_ab*s2*lnv_rates*y[9]*y[4]*pow(y20im,2.) + 
                        4*c_aa*g2*y[10]*y[4]*pow(y20im,2.) - 4*c_aa*s2*lnv_rates*y[10]*y[4]*pow(y20im,2.) - 
                        4*c_ab*g1_r_to_mu*g1*y[8]*pow(y20re,2.) + 4*c_ab*s1*lnv_rates*y[8]*pow(y20re,2.) - 
                        4*c_ab*g1_r_to_mu*g1*y[9]*pow(y20re,2.) + 4*c_ab*s1*lnv_rates*y[9]*pow(y20re,2.) - 
                        4*c_aa*g1_r_to_mu*g1*y[10]*pow(y20re,2.) + 4*c_aa*s1*lnv_rates*y[10]*pow(y20re,2.) + 
                        4*c_ab*g2*y[8]*y[4]*pow(y20re,2.) - 4*c_ab*s2*lnv_rates*y[8]*y[4]*pow(y20re,2.) + 
                        4*c_ab*g2*y[9]*y[4]*pow(y20re,2.) - 4*c_ab*s2*lnv_rates*y[9]*y[4]*pow(y20re,2.) + 
                        4*c_aa*g2*y[10]*y[4]*pow(y20re,2.) - 4*c_aa*s2*lnv_rates*y[10]*y[4]*pow(y20re,2.) + 
                        _KAV_*y[7]*y20im*y21im + 4*c_ab*g2*y[8]*y[6]*y20im*y21im - 
                        4*c_ab*s2*lnv_rates*y[8]*y[6]*y20im*y21im + 4*c_ab*g2*y[9]*y[6]*y20im*y21im - 
                        4*c_ab*s2*lnv_rates*y[9]*y[6]*y20im*y21im + 4*c_aa*g2*y[10]*y[6]*y20im*y21im - 
                        4*c_aa*s2*lnv_rates*y[10]*y[6]*y20im*y21im - 4*c_ab*g2*y[8]*y[7]*y20re*y21im - 
                        4*c_ab*s2*lnv_rates*y[8]*y[7]*y20re*y21im - 4*c_ab*g2*y[9]*y[7]*y20re*y21im - 
                        4*c_ab*s2*lnv_rates*y[9]*y[7]*y20re*y21im - 4*c_aa*g2*y[10]*y[7]*y20re*y21im - 
                        4*c_aa*s2*lnv_rates*y[10]*y[7]*y20re*y21im + _KAV_*y[6]*y20re*y21im + 
                        (4*(g2 + s2*lnv_rates)*(c_ab*(y[8] + y[9]) + c_aa*y[10])*y[7]*y20im - _KAV_*y[6]*y20im + 
                        _KAV_*y[7]*y20re + 4*(g2 - s2*lnv_rates)*(c_ab*(y[8] + y[9]) + c_aa*y[10])*y[6]*y20re)*
                        y21re + 4*s0*lnv_rates*((-1 + y[4])*pow(y00im,2.) - pow(y00re,2.) + y[7]*y00re*y01im + 
                        y[6]*y00re*y01re + y00im*(y[6]*y01im - y[7]*y01re) - pow(y10im,2.) - 
                        pow(y10re,2.) + y[6]*y10im*y11im + y[7]*y10re*y11im - y[7]*y10im*y11re + 
                        y[6]*y10re*y11re - pow(y20im,2.) - pow(y20re,2.) + 
                        y[4]*(pow(y00re,2.) + pow(y10im,2.) + pow(y10re,2.) + pow(y20im,2.) + pow(y20re,2.)) + 
                        y[6]*y20im*y21im + y[7]*y20re*y21im - y[7]*y20im*y21re + y[6]*y20re*y21re) + 
                        4*g0*((-1 + y[4])*pow(y00im,2.) - pow(y00re,2.) - y[7]*y00re*y01im + 
                        y[6]*y00re*y01re + y00im*(y[6]*y01im + y[7]*y01re) - pow(y10im,2.) - 
                        pow(y10re,2.) + y[6]*y10im*y11im - y[7]*y10re*y11im + y[7]*y10im*y11re + 
                        y[6]*y10re*y11re - pow(y20im,2.) - pow(y20re,2.) + 
                        y[4]*(pow(y00re,2.) + pow(y10im,2.) + pow(y10re,2.) + pow(y20im,2.) + pow(y20re,2.)) + 
                        y[6]*y20im*y21im - y[7]*y20re*y21im + y[7]*y20im*y21re + y[6]*y20re*y21re))) - //stef, new mthM term
                        2*_MPLTSP_*thmlnv*VNM*(y[6]*(-(y00re*y01im) + y00im*y01re - y10re*y11im + y10im*y11re - y20re*y21im + 
                        y20im*y21re) + y[7]*(y00im*y01im + y00re*y01re + y10im*y11im + y10re*y11re + y20im*y21im + 
                        y20re*y21re));


        yaux[5]         = (_KAV_*_MPLTSP_*(y[6]*(y00re*y01im - y00im*y01re + y10re*y11im - y10im*y11re + y20re*y21im - 
                        y20im*y21re) + y[7]*(y00im*y01im + y00re*y01re + y10im*y11im + y10re*y11re + 
                        y20im*y21im + y20re*y21re)))/4. + 
                        _MPLTSP_*(-(c_aa*g2*y[8]*y[6]*y00im*y01im) + c_aa*s2*lnv_rates*y[8]*y[6]*y00im*y01im - 
                        c_ab*g2*y[9]*y[6]*y00im*y01im + c_ab*s2*lnv_rates*y[9]*y[6]*y00im*y01im - 
                        c_ab*g2*y[10]*y[6]*y00im*y01im + c_ab*s2*lnv_rates*y[10]*y[6]*y00im*y01im + 
                        c_aa*g2*y[8]*y[7]*y00re*y01im + c_aa*s2*lnv_rates*y[8]*y[7]*y00re*y01im + 
                        c_ab*g2*y[9]*y[7]*y00re*y01im + c_ab*s2*lnv_rates*y[9]*y[7]*y00re*y01im + 
                        c_ab*g2*y[10]*y[7]*y00re*y01im + c_ab*s2*lnv_rates*y[10]*y[7]*y00re*y01im + 
                        c_aa*g1_r_to_mu*g1*y[8]*pow(y01im,2.) - c_aa*s1*lnv_rates*y[8]*pow(y01im,2.) + c_ab*g1_r_to_mu*g1*y[9]*pow(y01im,2.) - 
                        c_ab*s1*lnv_rates*y[9]*pow(y01im,2.) + c_ab*g1_r_to_mu*g1*y[10]*pow(y01im,2.) - 
                        c_ab*s1*lnv_rates*y[10]*pow(y01im,2.) - c_aa*g2*y[8]*y[5]*pow(y01im,2.) + 
                        c_aa*s2*lnv_rates*y[8]*y[5]*pow(y01im,2.) - c_ab*g2*y[9]*y[5]*pow(y01im,2.) + 
                        c_ab*s2*lnv_rates*y[9]*y[5]*pow(y01im,2.) - c_ab*g2*y[10]*y[5]*pow(y01im,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[5]*pow(y01im,2.) - c_aa*g2*y[8]*y[7]*y00im*y01re - 
                        c_aa*s2*lnv_rates*y[8]*y[7]*y00im*y01re - c_ab*g2*y[9]*y[7]*y00im*y01re - 
                        c_ab*s2*lnv_rates*y[9]*y[7]*y00im*y01re - c_ab*g2*y[10]*y[7]*y00im*y01re - 
                        c_ab*s2*lnv_rates*y[10]*y[7]*y00im*y01re - c_aa*g2*y[8]*y[6]*y00re*y01re + 
                        c_aa*s2*lnv_rates*y[8]*y[6]*y00re*y01re - c_ab*g2*y[9]*y[6]*y00re*y01re + 
                        c_ab*s2*lnv_rates*y[9]*y[6]*y00re*y01re - c_ab*g2*y[10]*y[6]*y00re*y01re + 
                        c_ab*s2*lnv_rates*y[10]*y[6]*y00re*y01re + c_aa*g1_r_to_mu*g1*y[8]*pow(y01re,2.) - 
                        c_aa*s1*lnv_rates*y[8]*pow(y01re,2.) + c_ab*g1_r_to_mu*g1*y[9]*pow(y01re,2.) - 
                        c_ab*s1*lnv_rates*y[9]*pow(y01re,2.) + c_ab*g1_r_to_mu*g1*y[10]*pow(y01re,2.) - 
                        c_ab*s1*lnv_rates*y[10]*pow(y01re,2.) - c_aa*g2*y[8]*y[5]*pow(y01re,2.) + 
                        c_aa*s2*lnv_rates*y[8]*y[5]*pow(y01re,2.) - c_ab*g2*y[9]*y[5]*pow(y01re,2.) + 
                        c_ab*s2*lnv_rates*y[9]*y[5]*pow(y01re,2.) - c_ab*g2*y[10]*y[5]*pow(y01re,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[5]*pow(y01re,2.) - c_ab*g2*y[8]*y[6]*y10im*y11im + 
                        c_ab*s2*lnv_rates*y[8]*y[6]*y10im*y11im - c_aa*g2*y[9]*y[6]*y10im*y11im + 
                        c_aa*s2*lnv_rates*y[9]*y[6]*y10im*y11im - c_ab*g2*y[10]*y[6]*y10im*y11im + 
                        c_ab*s2*lnv_rates*y[10]*y[6]*y10im*y11im + c_ab*g2*y[8]*y[7]*y10re*y11im + 
                        c_ab*s2*lnv_rates*y[8]*y[7]*y10re*y11im + c_aa*g2*y[9]*y[7]*y10re*y11im + 
                        c_aa*s2*lnv_rates*y[9]*y[7]*y10re*y11im + c_ab*g2*y[10]*y[7]*y10re*y11im + 
                        c_ab*s2*lnv_rates*y[10]*y[7]*y10re*y11im + c_ab*g1_r_to_mu*g1*y[8]*pow(y11im,2.) - 
                        c_ab*s1*lnv_rates*y[8]*pow(y11im,2.) + c_aa*g1_r_to_mu*g1*y[9]*pow(y11im,2.) - 
                        c_aa*s1*lnv_rates*y[9]*pow(y11im,2.) + c_ab*g1_r_to_mu*g1*y[10]*pow(y11im,2.) - 
                        c_ab*s1*lnv_rates*y[10]*pow(y11im,2.) - c_ab*g2*y[8]*y[5]*pow(y11im,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[5]*pow(y11im,2.) - c_aa*g2*y[9]*y[5]*pow(y11im,2.) + 
                        c_aa*s2*lnv_rates*y[9]*y[5]*pow(y11im,2.) - c_ab*g2*y[10]*y[5]*pow(y11im,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[5]*pow(y11im,2.) - c_ab*g2*y[8]*y[7]*y10im*y11re - 
                        c_ab*s2*lnv_rates*y[8]*y[7]*y10im*y11re - c_aa*g2*y[9]*y[7]*y10im*y11re - 
                        c_aa*s2*lnv_rates*y[9]*y[7]*y10im*y11re - c_ab*g2*y[10]*y[7]*y10im*y11re - 
                        c_ab*s2*lnv_rates*y[10]*y[7]*y10im*y11re - c_ab*g2*y[8]*y[6]*y10re*y11re + 
                        c_ab*s2*lnv_rates*y[8]*y[6]*y10re*y11re - c_aa*g2*y[9]*y[6]*y10re*y11re + 
                        c_aa*s2*lnv_rates*y[9]*y[6]*y10re*y11re - c_ab*g2*y[10]*y[6]*y10re*y11re + 
                        c_ab*s2*lnv_rates*y[10]*y[6]*y10re*y11re + c_ab*g1_r_to_mu*g1*y[8]*pow(y11re,2.) - 
                        c_ab*s1*lnv_rates*y[8]*pow(y11re,2.) + c_aa*g1_r_to_mu*g1*y[9]*pow(y11re,2.) - 
                        c_aa*s1*lnv_rates*y[9]*pow(y11re,2.) + c_ab*g1_r_to_mu*g1*y[10]*pow(y11re,2.) - 
                        c_ab*s1*lnv_rates*y[10]*pow(y11re,2.) - c_ab*g2*y[8]*y[5]*pow(y11re,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[5]*pow(y11re,2.) - c_aa*g2*y[9]*y[5]*pow(y11re,2.) + 
                        c_aa*s2*lnv_rates*y[9]*y[5]*pow(y11re,2.) - c_ab*g2*y[10]*y[5]*pow(y11re,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[5]*pow(y11re,2.) - c_ab*g2*y[8]*y[6]*y20im*y21im + 
                        c_ab*s2*lnv_rates*y[8]*y[6]*y20im*y21im - c_ab*g2*y[9]*y[6]*y20im*y21im + 
                        c_ab*s2*lnv_rates*y[9]*y[6]*y20im*y21im - c_aa*g2*y[10]*y[6]*y20im*y21im + 
                        c_aa*s2*lnv_rates*y[10]*y[6]*y20im*y21im + c_ab*g2*y[8]*y[7]*y20re*y21im + 
                        c_ab*s2*lnv_rates*y[8]*y[7]*y20re*y21im + c_ab*g2*y[9]*y[7]*y20re*y21im + 
                        c_ab*s2*lnv_rates*y[9]*y[7]*y20re*y21im + c_aa*g2*y[10]*y[7]*y20re*y21im + 
                        c_aa*s2*lnv_rates*y[10]*y[7]*y20re*y21im + c_ab*g1_r_to_mu*g1*y[8]*pow(y21im,2.) - 
                        c_ab*s1*lnv_rates*y[8]*pow(y21im,2.) + c_ab*g1_r_to_mu*g1*y[9]*pow(y21im,2.) - 
                        c_ab*s1*lnv_rates*y[9]*pow(y21im,2.) + c_aa*g1_r_to_mu*g1*y[10]*pow(y21im,2.) - 
                        c_aa*s1*lnv_rates*y[10]*pow(y21im,2.) - c_ab*g2*y[8]*y[5]*pow(y21im,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[5]*pow(y21im,2.) - c_ab*g2*y[9]*y[5]*pow(y21im,2.) + 
                        c_ab*s2*lnv_rates*y[9]*y[5]*pow(y21im,2.) - c_aa*g2*y[10]*y[5]*pow(y21im,2.) + 
                        c_aa*s2*lnv_rates*y[10]*y[5]*pow(y21im,2.) - 
                        (c_ab*(y[8] + y[9]) + c_aa*y[10])*((g2 + s2*lnv_rates)*y[7]*y20im + (g2 - s2*lnv_rates)*y[6]*y20re)*
                        y21re + (c_ab*(y[8] + y[9]) + c_aa*y[10])*(g1_r_to_mu*g1 - s1*lnv_rates - g2*y[5] + s2*lnv_rates*y[5])*
                        pow(y21re,2.) - s0*lnv_rates*(y[7]*(y00re*y01im - y00im*y01re + y10re*y11im - y10im*y11re + 
                        y20re*y21im - y20im*y21re) + y[6]*
                        (y00im*y01im + y00re*y01re + y10im*y11im + y10re*y11re + y20im*y21im + y20re*y21re) + 
                        (-1 + y[5])*(pow(y01im,2.) + pow(y01re,2.) + pow(y11im,2.) + pow(y11re,2.) + 
                        pow(y21im,2.) + pow(y21re,2.))) - 
                        g0*(y[7]*(-(y00re*y01im) + y00im*y01re - y10re*y11im + y10im*y11re - y20re*y21im + 
                        y20im*y21re) + y[6]*(y00im*y01im + y00re*y01re + y10im*y11im + y10re*y11re + 
                        y20im*y21im + y20re*y21re) + (-1 + y[5])*
                        (pow(y01im,2.) + pow(y01re,2.) + pow(y11im,2.) + pow(y11re,2.) + pow(y21im,2.) + 
                        pow(y21re,2.)))) + //stef, new mthM term
                        2*_MPLTSP_*thmlnv*VNM*(y[6]*(-(y00re*y01im) + y00im*y01re - y10re*y11im + y10im*y11re - y20re*y21im + 
                        y20im*y21re) + y[7]*(y00im*y01im + y00re*y01re + y10im*y11im + y10re*y11re + y20im*y21im + 
                        y20re*y21re));

        
        yaux[6]         = -0.125*(_MPLTSP_*(_KAV_*(4*d12*y[7]*pow(t,2.) + 
                        pow(_TSP_,2.)*(-((y[4] - y[5])*
                        (y00re*y01im - y00im*y01re + y10re*y11im - y10im*y11re + y20re*y21im - y20im*y21re)) + 
                        y[7]*(-pow(y00im,2.) - pow(y00re,2.) + pow(y01im,2.) + pow(y01re,2.) - 
                        pow(y10im,2.) - pow(y10re,2.) + pow(y11im,2.) + pow(y11re,2.) - pow(y20im,2.) - 
                        pow(y20re,2.) + pow(y21im,2.) + pow(y21re,2.)))) + 
                        4*pow(_TSP_,2.)*(c_aa*g2*y[8]*y[6]*pow(y00im,2.) - 
                        c_aa*s2*lnv_rates*y[8]*y[6]*pow(y00im,2.) + c_ab*g2*y[9]*y[6]*pow(y00im,2.) - 
                        c_ab*s2*lnv_rates*y[9]*y[6]*pow(y00im,2.) + c_ab*g2*y[10]*y[6]*pow(y00im,2.) - 
                        c_ab*s2*lnv_rates*y[10]*y[6]*pow(y00im,2.) + c_aa*g2*y[8]*y[6]*pow(y00re,2.) - 
                        c_aa*s2*lnv_rates*y[8]*y[6]*pow(y00re,2.) + c_ab*g2*y[9]*y[6]*pow(y00re,2.) - 
                        c_ab*s2*lnv_rates*y[9]*y[6]*pow(y00re,2.) + c_ab*g2*y[10]*y[6]*pow(y00re,2.) - 
                        c_ab*s2*lnv_rates*y[10]*y[6]*pow(y00re,2.) - 2*c_aa*g1_r_to_mu*g1*y[8]*y00im*y01im + 
                        2*c_aa*s1*lnv_rates*y[8]*y00im*y01im - 2*c_ab*g1_r_to_mu*g1*y[9]*y00im*y01im + 
                        2*c_ab*s1*lnv_rates*y[9]*y00im*y01im - 2*c_ab*g1_r_to_mu*g1*y[10]*y00im*y01im + 
                        2*c_ab*s1*lnv_rates*y[10]*y00im*y01im + c_aa*g2*y[8]*y[4]*y00im*y01im - 
                        c_aa*s2*lnv_rates*y[8]*y[4]*y00im*y01im + c_ab*g2*y[9]*y[4]*y00im*y01im - 
                        c_ab*s2*lnv_rates*y[9]*y[4]*y00im*y01im + c_ab*g2*y[10]*y[4]*y00im*y01im - 
                        c_ab*s2*lnv_rates*y[10]*y[4]*y00im*y01im + c_aa*g2*y[8]*y[5]*y00im*y01im - 
                        c_aa*s2*lnv_rates*y[8]*y[5]*y00im*y01im + c_ab*g2*y[9]*y[5]*y00im*y01im - 
                        c_ab*s2*lnv_rates*y[9]*y[5]*y00im*y01im + c_ab*g2*y[10]*y[5]*y00im*y01im - 
                        c_ab*s2*lnv_rates*y[10]*y[5]*y00im*y01im + c_aa*g2*y[8]*y[6]*pow(y01im,2.) - 
                        c_aa*s2*lnv_rates*y[8]*y[6]*pow(y01im,2.) + c_ab*g2*y[9]*y[6]*pow(y01im,2.) - 
                        c_ab*s2*lnv_rates*y[9]*y[6]*pow(y01im,2.) + c_ab*g2*y[10]*y[6]*pow(y01im,2.) - 
                        c_ab*s2*lnv_rates*y[10]*y[6]*pow(y01im,2.) - 2*c_aa*g1_r_to_mu*g1*y[8]*y00re*y01re + 
                        2*c_aa*s1*lnv_rates*y[8]*y00re*y01re - 2*c_ab*g1_r_to_mu*g1*y[9]*y00re*y01re + 
                        2*c_ab*s1*lnv_rates*y[9]*y00re*y01re - 2*c_ab*g1_r_to_mu*g1*y[10]*y00re*y01re + 
                        2*c_ab*s1*lnv_rates*y[10]*y00re*y01re + c_aa*g2*y[8]*y[4]*y00re*y01re - 
                        c_aa*s2*lnv_rates*y[8]*y[4]*y00re*y01re + c_ab*g2*y[9]*y[4]*y00re*y01re - 
                        c_ab*s2*lnv_rates*y[9]*y[4]*y00re*y01re + c_ab*g2*y[10]*y[4]*y00re*y01re - 
                        c_ab*s2*lnv_rates*y[10]*y[4]*y00re*y01re + c_aa*g2*y[8]*y[5]*y00re*y01re - 
                        c_aa*s2*lnv_rates*y[8]*y[5]*y00re*y01re + c_ab*g2*y[9]*y[5]*y00re*y01re - 
                        c_ab*s2*lnv_rates*y[9]*y[5]*y00re*y01re + c_ab*g2*y[10]*y[5]*y00re*y01re - 
                        c_ab*s2*lnv_rates*y[10]*y[5]*y00re*y01re + c_aa*g2*y[8]*y[6]*pow(y01re,2.) - 
                        c_aa*s2*lnv_rates*y[8]*y[6]*pow(y01re,2.) + c_ab*g2*y[9]*y[6]*pow(y01re,2.) - 
                        c_ab*s2*lnv_rates*y[9]*y[6]*pow(y01re,2.) + c_ab*g2*y[10]*y[6]*pow(y01re,2.) - 
                        c_ab*s2*lnv_rates*y[10]*y[6]*pow(y01re,2.) + c_ab*g2*y[8]*y[6]*pow(y10im,2.) - 
                        c_ab*s2*lnv_rates*y[8]*y[6]*pow(y10im,2.) + c_aa*g2*y[9]*y[6]*pow(y10im,2.) - 
                        c_aa*s2*lnv_rates*y[9]*y[6]*pow(y10im,2.) + c_ab*g2*y[10]*y[6]*pow(y10im,2.) - 
                        c_ab*s2*lnv_rates*y[10]*y[6]*pow(y10im,2.) + c_ab*g2*y[8]*y[6]*pow(y10re,2.) - 
                        c_ab*s2*lnv_rates*y[8]*y[6]*pow(y10re,2.) + c_aa*g2*y[9]*y[6]*pow(y10re,2.) - 
                        c_aa*s2*lnv_rates*y[9]*y[6]*pow(y10re,2.) + c_ab*g2*y[10]*y[6]*pow(y10re,2.) - 
                        c_ab*s2*lnv_rates*y[10]*y[6]*pow(y10re,2.) - 2*c_ab*g1_r_to_mu*g1*y[8]*y10im*y11im + 
                        2*c_ab*s1*lnv_rates*y[8]*y10im*y11im - 2*c_aa*g1_r_to_mu*g1*y[9]*y10im*y11im + 
                        2*c_aa*s1*lnv_rates*y[9]*y10im*y11im - 2*c_ab*g1_r_to_mu*g1*y[10]*y10im*y11im + 
                        2*c_ab*s1*lnv_rates*y[10]*y10im*y11im + c_ab*g2*y[8]*y[4]*y10im*y11im - 
                        c_ab*s2*lnv_rates*y[8]*y[4]*y10im*y11im + c_aa*g2*y[9]*y[4]*y10im*y11im - 
                        c_aa*s2*lnv_rates*y[9]*y[4]*y10im*y11im + c_ab*g2*y[10]*y[4]*y10im*y11im - 
                        c_ab*s2*lnv_rates*y[10]*y[4]*y10im*y11im + c_ab*g2*y[8]*y[5]*y10im*y11im - 
                        c_ab*s2*lnv_rates*y[8]*y[5]*y10im*y11im + c_aa*g2*y[9]*y[5]*y10im*y11im - 
                        c_aa*s2*lnv_rates*y[9]*y[5]*y10im*y11im + c_ab*g2*y[10]*y[5]*y10im*y11im - 
                        c_ab*s2*lnv_rates*y[10]*y[5]*y10im*y11im + c_ab*g2*y[8]*y[6]*pow(y11im,2.) - 
                        c_ab*s2*lnv_rates*y[8]*y[6]*pow(y11im,2.) + c_aa*g2*y[9]*y[6]*pow(y11im,2.) - 
                        c_aa*s2*lnv_rates*y[9]*y[6]*pow(y11im,2.) + c_ab*g2*y[10]*y[6]*pow(y11im,2.) - 
                        c_ab*s2*lnv_rates*y[10]*y[6]*pow(y11im,2.) - 2*c_ab*g1_r_to_mu*g1*y[8]*y10re*y11re + 
                        2*c_ab*s1*lnv_rates*y[8]*y10re*y11re - 2*c_aa*g1_r_to_mu*g1*y[9]*y10re*y11re + 
                        2*c_aa*s1*lnv_rates*y[9]*y10re*y11re - 2*c_ab*g1_r_to_mu*g1*y[10]*y10re*y11re + 
                        2*c_ab*s1*lnv_rates*y[10]*y10re*y11re + c_ab*g2*y[8]*y[4]*y10re*y11re - 
                        c_ab*s2*lnv_rates*y[8]*y[4]*y10re*y11re + c_aa*g2*y[9]*y[4]*y10re*y11re - 
                        c_aa*s2*lnv_rates*y[9]*y[4]*y10re*y11re + c_ab*g2*y[10]*y[4]*y10re*y11re - 
                        c_ab*s2*lnv_rates*y[10]*y[4]*y10re*y11re + c_ab*g2*y[8]*y[5]*y10re*y11re - 
                        c_ab*s2*lnv_rates*y[8]*y[5]*y10re*y11re + c_aa*g2*y[9]*y[5]*y10re*y11re - 
                        c_aa*s2*lnv_rates*y[9]*y[5]*y10re*y11re + c_ab*g2*y[10]*y[5]*y10re*y11re - 
                        c_ab*s2*lnv_rates*y[10]*y[5]*y10re*y11re + c_ab*g2*y[8]*y[6]*pow(y11re,2.) - 
                        c_ab*s2*lnv_rates*y[8]*y[6]*pow(y11re,2.) + c_aa*g2*y[9]*y[6]*pow(y11re,2.) - 
                        c_aa*s2*lnv_rates*y[9]*y[6]*pow(y11re,2.) + c_ab*g2*y[10]*y[6]*pow(y11re,2.) - 
                        c_ab*s2*lnv_rates*y[10]*y[6]*pow(y11re,2.) + c_ab*g2*y[8]*y[6]*pow(y20im,2.) - 
                        c_ab*s2*lnv_rates*y[8]*y[6]*pow(y20im,2.) + c_ab*g2*y[9]*y[6]*pow(y20im,2.) - 
                        c_ab*s2*lnv_rates*y[9]*y[6]*pow(y20im,2.) + c_aa*g2*y[10]*y[6]*pow(y20im,2.) - 
                        c_aa*s2*lnv_rates*y[10]*y[6]*pow(y20im,2.) + c_ab*g2*y[8]*y[6]*pow(y20re,2.) - 
                        c_ab*s2*lnv_rates*y[8]*y[6]*pow(y20re,2.) + c_ab*g2*y[9]*y[6]*pow(y20re,2.) - 
                        c_ab*s2*lnv_rates*y[9]*y[6]*pow(y20re,2.) + c_aa*g2*y[10]*y[6]*pow(y20re,2.) - 
                        c_aa*s2*lnv_rates*y[10]*y[6]*pow(y20re,2.) - 2*c_ab*g1_r_to_mu*g1*y[8]*y20im*y21im + 
                        2*c_ab*s1*lnv_rates*y[8]*y20im*y21im - 2*c_ab*g1_r_to_mu*g1*y[9]*y20im*y21im + 
                        2*c_ab*s1*lnv_rates*y[9]*y20im*y21im - 2*c_aa*g1_r_to_mu*g1*y[10]*y20im*y21im + 
                        2*c_aa*s1*lnv_rates*y[10]*y20im*y21im + c_ab*g2*y[8]*y[4]*y20im*y21im - 
                        c_ab*s2*lnv_rates*y[8]*y[4]*y20im*y21im + c_ab*g2*y[9]*y[4]*y20im*y21im - 
                        c_ab*s2*lnv_rates*y[9]*y[4]*y20im*y21im + c_aa*g2*y[10]*y[4]*y20im*y21im - 
                        c_aa*s2*lnv_rates*y[10]*y[4]*y20im*y21im + c_ab*g2*y[8]*y[5]*y20im*y21im - 
                        c_ab*s2*lnv_rates*y[8]*y[5]*y20im*y21im + c_ab*g2*y[9]*y[5]*y20im*y21im - 
                        c_ab*s2*lnv_rates*y[9]*y[5]*y20im*y21im + c_aa*g2*y[10]*y[5]*y20im*y21im - 
                        c_aa*s2*lnv_rates*y[10]*y[5]*y20im*y21im + c_ab*g2*y[8]*y[6]*pow(y21im,2.) - 
                        c_ab*s2*lnv_rates*y[8]*y[6]*pow(y21im,2.) + c_ab*g2*y[9]*y[6]*pow(y21im,2.) - 
                        c_ab*s2*lnv_rates*y[9]*y[6]*pow(y21im,2.) + c_aa*g2*y[10]*y[6]*pow(y21im,2.) - 
                        c_aa*s2*lnv_rates*y[10]*y[6]*pow(y21im,2.) - 
                        (c_ab*(y[8] + y[9]) + c_aa*y[10])*
                        (2*g1_r_to_mu*g1 - 2*s1*lnv_rates - g2*(y[4] + y[5]) + s2*lnv_rates*(y[4] + y[5]))*y20re*y21re + 
                        (g2 - s2*lnv_rates)*(c_ab*(y[8] + y[9]) + c_aa*y[10])*y[6]*pow(y21re,2.) + 
                        g0*((-2 + y[4] + y[5])*(y00im*y01im + y00re*y01re + y10im*y11im + y10re*y11re + 
                        y20im*y21im + y20re*y21re) + 
                        y[6]*(pow(y00im,2.) + pow(y00re,2.) + pow(y01im,2.) + pow(y01re,2.) + 
                        pow(y10im,2.) + pow(y10re,2.) + pow(y11im,2.) + pow(y11re,2.) + pow(y20im,2.) + 
                        pow(y20re,2.) + pow(y21im,2.) + pow(y21re,2.))) + 
                        s0*lnv_rates*((-2 + y[4] + y[5])*(y00im*y01im + y00re*y01re + y10im*y11im + y10re*y11re + 
                        y20im*y21im + y20re*y21re) + 
                        y[6]*(pow(y00im,2.) + pow(y00re,2.) + pow(y01im,2.) + pow(y01re,2.) + 
                        pow(y10im,2.) + pow(y10re,2.) + pow(y11im,2.) + pow(y11re,2.) + pow(y20im,2.) + 
                        pow(y20re,2.) + pow(y21im,2.) + pow(y21re,2.))))))/pow(_TSP_,2.) + //stef, new mthM term
                        _MPLTSP_*thmlnv*VNM*(-((y[4] - y[5])*(y00re*y01im - y00im*y01re + y10re*y11im - y10im*y11re + 
                        y20re*y21im - y20im*y21re)) + y[7]*
                        (pow(y00im,2.) + pow(y00re,2.) - pow(y01im,2.) - pow(y01re,2.) + pow(y10im,2.) + 
                        pow(y10re,2.) - pow(y11im,2.) - pow(y11re,2.) + pow(y20im,2.) + pow(y20re,2.) - 
                        pow(y21im,2.) - pow(y21re,2.)));


        yaux[7]         = (_MPLTSP_*(-4*(g0 + s0*lnv_rates + (g2 - s2*lnv_rates)*(c_aa*y[8] + c_ab*(y[9] + y[10])))*y[7]*
                        (pow(y00im,2.) + pow(y00re,2.)) + 
                        _KAV_*((4*d12*y[6]*pow(t,2.))/pow(_TSP_,2.) + 
                        (y[4] - y[5])*(y00im*y01im + y00re*y01re + y10im*y11im + y10re*y11re + y20im*y21im + 
                        y20re*y21re) + y[6]*(-pow(y00im,2.) - pow(y00re,2.) + pow(y01im,2.) + 
                        pow(y01re,2.) - pow(y10im,2.) - pow(y10re,2.) + pow(y11im,2.) + pow(y11re,2.) - 
                        pow(y20im,2.) - pow(y20re,2.) + pow(y21im,2.) + pow(y21re,2.))) + 
                        4*(-2*c_aa*g1_r_to_mu*g1*y[8]*y00re*y01im - 2*c_aa*s1*lnv_rates*y[8]*y00re*y01im - 
                        2*c_ab*g1_r_to_mu*g1*y[9]*y00re*y01im - 2*c_ab*s1*lnv_rates*y[9]*y00re*y01im - 2*c_ab*g1_r_to_mu*g1*y[10]*y00re*y01im - 
                        2*c_ab*s1*lnv_rates*y[10]*y00re*y01im + c_aa*g2*y[8]*y[4]*y00re*y01im + 
                        c_aa*s2*lnv_rates*y[8]*y[4]*y00re*y01im + c_ab*g2*y[9]*y[4]*y00re*y01im + 
                        c_ab*s2*lnv_rates*y[9]*y[4]*y00re*y01im + c_ab*g2*y[10]*y[4]*y00re*y01im + 
                        c_ab*s2*lnv_rates*y[10]*y[4]*y00re*y01im + c_aa*g2*y[8]*y[5]*y00re*y01im + 
                        c_aa*s2*lnv_rates*y[8]*y[5]*y00re*y01im + c_ab*g2*y[9]*y[5]*y00re*y01im + 
                        c_ab*s2*lnv_rates*y[9]*y[5]*y00re*y01im + c_ab*g2*y[10]*y[5]*y00re*y01im + 
                        c_ab*s2*lnv_rates*y[10]*y[5]*y00re*y01im - c_aa*g2*y[8]*y[7]*pow(y01im,2.) + 
                        c_aa*s2*lnv_rates*y[8]*y[7]*pow(y01im,2.) - c_ab*g2*y[9]*y[7]*pow(y01im,2.) + 
                        c_ab*s2*lnv_rates*y[9]*y[7]*pow(y01im,2.) - c_ab*g2*y[10]*y[7]*pow(y01im,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[7]*pow(y01im,2.) + 2*c_aa*g1_r_to_mu*g1*y[8]*y00im*y01re + 
                        2*c_aa*s1*lnv_rates*y[8]*y00im*y01re + 2*c_ab*g1_r_to_mu*g1*y[9]*y00im*y01re + 
                        2*c_ab*s1*lnv_rates*y[9]*y00im*y01re + 2*c_ab*g1_r_to_mu*g1*y[10]*y00im*y01re + 
                        2*c_ab*s1*lnv_rates*y[10]*y00im*y01re - c_aa*g2*y[8]*y[4]*y00im*y01re - 
                        c_aa*s2*lnv_rates*y[8]*y[4]*y00im*y01re - c_ab*g2*y[9]*y[4]*y00im*y01re - 
                        c_ab*s2*lnv_rates*y[9]*y[4]*y00im*y01re - c_ab*g2*y[10]*y[4]*y00im*y01re - 
                        c_ab*s2*lnv_rates*y[10]*y[4]*y00im*y01re - c_aa*g2*y[8]*y[5]*y00im*y01re - 
                        c_aa*s2*lnv_rates*y[8]*y[5]*y00im*y01re - c_ab*g2*y[9]*y[5]*y00im*y01re - 
                        c_ab*s2*lnv_rates*y[9]*y[5]*y00im*y01re - c_ab*g2*y[10]*y[5]*y00im*y01re - 
                        c_ab*s2*lnv_rates*y[10]*y[5]*y00im*y01re - c_aa*g2*y[8]*y[7]*pow(y01re,2.) + 
                        c_aa*s2*lnv_rates*y[8]*y[7]*pow(y01re,2.) - c_ab*g2*y[9]*y[7]*pow(y01re,2.) + 
                        c_ab*s2*lnv_rates*y[9]*y[7]*pow(y01re,2.) - c_ab*g2*y[10]*y[7]*pow(y01re,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[7]*pow(y01re,2.) - c_ab*g2*y[8]*y[7]*pow(y10im,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[7]*pow(y10im,2.) - c_aa*g2*y[9]*y[7]*pow(y10im,2.) + 
                        c_aa*s2*lnv_rates*y[9]*y[7]*pow(y10im,2.) - c_ab*g2*y[10]*y[7]*pow(y10im,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[7]*pow(y10im,2.) - c_ab*g2*y[8]*y[7]*pow(y10re,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[7]*pow(y10re,2.) - c_aa*g2*y[9]*y[7]*pow(y10re,2.) + 
                        c_aa*s2*lnv_rates*y[9]*y[7]*pow(y10re,2.) - c_ab*g2*y[10]*y[7]*pow(y10re,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[7]*pow(y10re,2.) - 2*c_ab*g1_r_to_mu*g1*y[8]*y10re*y11im - 
                        2*c_ab*s1*lnv_rates*y[8]*y10re*y11im - 2*c_aa*g1_r_to_mu*g1*y[9]*y10re*y11im - 
                        2*c_aa*s1*lnv_rates*y[9]*y10re*y11im - 2*c_ab*g1_r_to_mu*g1*y[10]*y10re*y11im - 
                        2*c_ab*s1*lnv_rates*y[10]*y10re*y11im + c_ab*g2*y[8]*y[4]*y10re*y11im + 
                        c_ab*s2*lnv_rates*y[8]*y[4]*y10re*y11im + c_aa*g2*y[9]*y[4]*y10re*y11im + 
                        c_aa*s2*lnv_rates*y[9]*y[4]*y10re*y11im + c_ab*g2*y[10]*y[4]*y10re*y11im + 
                        c_ab*s2*lnv_rates*y[10]*y[4]*y10re*y11im + c_ab*g2*y[8]*y[5]*y10re*y11im + 
                        c_ab*s2*lnv_rates*y[8]*y[5]*y10re*y11im + c_aa*g2*y[9]*y[5]*y10re*y11im + 
                        c_aa*s2*lnv_rates*y[9]*y[5]*y10re*y11im + c_ab*g2*y[10]*y[5]*y10re*y11im + 
                        c_ab*s2*lnv_rates*y[10]*y[5]*y10re*y11im - c_ab*g2*y[8]*y[7]*pow(y11im,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[7]*pow(y11im,2.) - c_aa*g2*y[9]*y[7]*pow(y11im,2.) + 
                        c_aa*s2*lnv_rates*y[9]*y[7]*pow(y11im,2.) - c_ab*g2*y[10]*y[7]*pow(y11im,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[7]*pow(y11im,2.) + 2*c_ab*g1_r_to_mu*g1*y[8]*y10im*y11re + 
                        2*c_ab*s1*lnv_rates*y[8]*y10im*y11re + 2*c_aa*g1_r_to_mu*g1*y[9]*y10im*y11re + 
                        2*c_aa*s1*lnv_rates*y[9]*y10im*y11re + 2*c_ab*g1_r_to_mu*g1*y[10]*y10im*y11re + 
                        2*c_ab*s1*lnv_rates*y[10]*y10im*y11re - c_ab*g2*y[8]*y[4]*y10im*y11re - 
                        c_ab*s2*lnv_rates*y[8]*y[4]*y10im*y11re - c_aa*g2*y[9]*y[4]*y10im*y11re - 
                        c_aa*s2*lnv_rates*y[9]*y[4]*y10im*y11re - c_ab*g2*y[10]*y[4]*y10im*y11re - 
                        c_ab*s2*lnv_rates*y[10]*y[4]*y10im*y11re - c_ab*g2*y[8]*y[5]*y10im*y11re - 
                        c_ab*s2*lnv_rates*y[8]*y[5]*y10im*y11re - c_aa*g2*y[9]*y[5]*y10im*y11re - 
                        c_aa*s2*lnv_rates*y[9]*y[5]*y10im*y11re - c_ab*g2*y[10]*y[5]*y10im*y11re - 
                        c_ab*s2*lnv_rates*y[10]*y[5]*y10im*y11re - c_ab*g2*y[8]*y[7]*pow(y11re,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[7]*pow(y11re,2.) - c_aa*g2*y[9]*y[7]*pow(y11re,2.) + 
                        c_aa*s2*lnv_rates*y[9]*y[7]*pow(y11re,2.) - c_ab*g2*y[10]*y[7]*pow(y11re,2.) + 
                        c_ab*s2*lnv_rates*y[10]*y[7]*pow(y11re,2.) - c_ab*g2*y[8]*y[7]*pow(y20im,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[7]*pow(y20im,2.) - c_ab*g2*y[9]*y[7]*pow(y20im,2.) + 
                        c_ab*s2*lnv_rates*y[9]*y[7]*pow(y20im,2.) - c_aa*g2*y[10]*y[7]*pow(y20im,2.) + 
                        c_aa*s2*lnv_rates*y[10]*y[7]*pow(y20im,2.) - c_ab*g2*y[8]*y[7]*pow(y20re,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[7]*pow(y20re,2.) - c_ab*g2*y[9]*y[7]*pow(y20re,2.) + 
                        c_ab*s2*lnv_rates*y[9]*y[7]*pow(y20re,2.) - c_aa*g2*y[10]*y[7]*pow(y20re,2.) + 
                        c_aa*s2*lnv_rates*y[10]*y[7]*pow(y20re,2.) - 2*c_ab*g1_r_to_mu*g1*y[8]*y20re*y21im - 
                        2*c_ab*s1*lnv_rates*y[8]*y20re*y21im - 2*c_ab*g1_r_to_mu*g1*y[9]*y20re*y21im - 
                        2*c_ab*s1*lnv_rates*y[9]*y20re*y21im - 2*c_aa*g1_r_to_mu*g1*y[10]*y20re*y21im - 
                        2*c_aa*s1*lnv_rates*y[10]*y20re*y21im + c_ab*g2*y[8]*y[4]*y20re*y21im + 
                        c_ab*s2*lnv_rates*y[8]*y[4]*y20re*y21im + c_ab*g2*y[9]*y[4]*y20re*y21im + 
                        c_ab*s2*lnv_rates*y[9]*y[4]*y20re*y21im + c_aa*g2*y[10]*y[4]*y20re*y21im + 
                        c_aa*s2*lnv_rates*y[10]*y[4]*y20re*y21im + c_ab*g2*y[8]*y[5]*y20re*y21im + 
                        c_ab*s2*lnv_rates*y[8]*y[5]*y20re*y21im + c_ab*g2*y[9]*y[5]*y20re*y21im + 
                        c_ab*s2*lnv_rates*y[9]*y[5]*y20re*y21im + c_aa*g2*y[10]*y[5]*y20re*y21im + 
                        c_aa*s2*lnv_rates*y[10]*y[5]*y20re*y21im - c_ab*g2*y[8]*y[7]*pow(y21im,2.) + 
                        c_ab*s2*lnv_rates*y[8]*y[7]*pow(y21im,2.) - c_ab*g2*y[9]*y[7]*pow(y21im,2.) + 
                        c_ab*s2*lnv_rates*y[9]*y[7]*pow(y21im,2.) - c_aa*g2*y[10]*y[7]*pow(y21im,2.) + 
                        c_aa*s2*lnv_rates*y[10]*y[7]*pow(y21im,2.) + 
                        (c_ab*(y[8] + y[9]) + c_aa*y[10])*(2*g1_r_to_mu*g1 + 2*s1*lnv_rates - (g2 + s2*lnv_rates)*(y[4] + y[5]))*y20im*
                        y21re - (g2 - s2*lnv_rates)*(c_ab*(y[8] + y[9]) + c_aa*y[10])*y[7]*pow(y21re,2.) + 
                        g0*((-2 + y[4] + y[5])*y00re*y01im - 
                        (-2 + y[4] + y[5])*(y00im*y01re - y10re*y11im + y10im*y11re - y20re*y21im + 
                        y20im*y21re) - y[7]*(pow(y01im,2.) + pow(y01re,2.) + pow(y10im,2.) + 
                        pow(y10re,2.) + pow(y11im,2.) + pow(y11re,2.) + pow(y20im,2.) + pow(y20re,2.) + 
                        pow(y21im,2.) + pow(y21re,2.))) - 
                        s0*lnv_rates*((-2 + y[4] + y[5])*y00re*y01im - 
                        (-2 + y[4] + y[5])*(y00im*y01re - y10re*y11im + y10im*y11re - y20re*y21im + 
                        y20im*y21re) + y[7]*(pow(y01im,2.) + pow(y01re,2.) + pow(y10im,2.) + 
                        pow(y10re,2.) + pow(y11im,2.) + pow(y11re,2.) + pow(y20im,2.) + pow(y20re,2.) + 
                        pow(y21im,2.) + pow(y21re,2.))))))/8. +  //stef, new mthM term
                        _MPLTSP_*thmlnv*VNM*((y[4] - y[5])*(y00im*y01im + y00re*y01re + y10im*y11im + y10re*y11re + 
                        y20im*y21im + y20re*y21re) + y[6]*
                        (-pow(y00im,2.) - pow(y00re,2.) + pow(y01im,2.) + pow(y01re,2.) - pow(y10im,2.) - 
                        pow(y10re,2.) + pow(y11im,2.) + pow(y11re,2.) - pow(y20im,2.) - pow(y20re,2.) + 
                        pow(y21im,2.) + pow(y21re,2.)));



        yaux[8]         = (_FD_*_MPLTSP_*(-(g0*y[0]*pow(y00im,2.)) + s0*lnv_rates*y[0]*pow(y00im,2.) + g0*y[4]*pow(y00im,2.) - 
                        s0*lnv_rates*y[4]*pow(y00im,2.) - g0*y[0]*pow(y00re,2.) + s0*lnv_rates*y[0]*pow(y00re,2.) + 
                        g0*y[4]*pow(y00re,2.) - s0*lnv_rates*y[4]*pow(y00re,2.) - 2*g0*y[2]*y00im*y01im + 
                        2*s0*lnv_rates*y[2]*y00im*y01im + 2*g0*y[6]*y00im*y01im - 2*s0*lnv_rates*y[6]*y00im*y01im - 
                        2*g0*y[3]*y00re*y01im - 2*s0*lnv_rates*y[3]*y00re*y01im - 2*g0*y[7]*y00re*y01im - 
                        2*s0*lnv_rates*y[7]*y00re*y01im - g0*y[1]*pow(y01im,2.) + s0*lnv_rates*y[1]*pow(y01im,2.) + 
                        g0*y[5]*pow(y01im,2.) - s0*lnv_rates*y[5]*pow(y01im,2.) + 
                        2*(g0 + s0*lnv_rates)*(y[3] + y[7])*y00im*y01re - 2*(g0 - s0*lnv_rates)*(y[2] - y[6])*y00re*y01re - 
                        (g0 - s0*lnv_rates)*(y[1] - y[5])*pow(y01re,2.) + 
                        c_aa*y[8]*(g2*y[0]*pow(y00im,2.) + s2*lnv_rates*y[0]*pow(y00im,2.) + g2*y[4]*pow(y00im,2.) + 
                        s2*lnv_rates*y[4]*pow(y00im,2.) + g2*y[0]*pow(y00re,2.) + s2*lnv_rates*y[0]*pow(y00re,2.) + 
                        g2*y[4]*pow(y00re,2.) + s2*lnv_rates*y[4]*pow(y00re,2.) + 2*g2*y[2]*y00im*y01im + 
                        2*s2*lnv_rates*y[2]*y00im*y01im + 2*g2*y[6]*y00im*y01im + 2*s2*lnv_rates*y[6]*y00im*y01im + 
                        2*g2*y[3]*y00re*y01im - 2*s2*lnv_rates*y[3]*y00re*y01im - 2*g2*y[7]*y00re*y01im + 
                        2*s2*lnv_rates*y[7]*y00re*y01im + g2*y[1]*pow(y01im,2.) + s2*lnv_rates*y[1]*pow(y01im,2.) + 
                        g2*y[5]*pow(y01im,2.) + s2*lnv_rates*y[5]*pow(y01im,2.) - 
                        2*(g2 - s2*lnv_rates)*(y[3] - y[7])*y00im*y01re + 
                        2*(g2 + s2*lnv_rates)*(y[2] + y[6])*y00re*y01re + (g2 + s2*lnv_rates)*(y[1] + y[5])*pow(y01re,2.) - 
                        2*g1*(pow(y00im,2.) + pow(y00re,2.) + pow(y01im,2.) + pow(y01re,2.)) - 
                        2*s1*lnv_rates*(pow(y00im,2.) + pow(y00re,2.) + pow(y01im,2.) + pow(y01re,2.))) - 
                        c_ab*(y[9] + y[10])*(-(g2*y[0]*pow(y00im,2.)) - s2*lnv_rates*y[0]*pow(y00im,2.) - 
                        g2*y[4]*pow(y00im,2.) - s2*lnv_rates*y[4]*pow(y00im,2.) - g2*y[0]*pow(y00re,2.) - 
                        s2*lnv_rates*y[0]*pow(y00re,2.) - g2*y[4]*pow(y00re,2.) - s2*lnv_rates*y[4]*pow(y00re,2.) - 
                        2*g2*y[2]*y00im*y01im - 2*s2*lnv_rates*y[2]*y00im*y01im - 2*g2*y[6]*y00im*y01im - 
                        2*s2*lnv_rates*y[6]*y00im*y01im - 2*g2*y[3]*y00re*y01im + 2*s2*lnv_rates*y[3]*y00re*y01im + 
                        2*g2*y[7]*y00re*y01im - 2*s2*lnv_rates*y[7]*y00re*y01im - g2*y[1]*pow(y01im,2.) - 
                        s2*lnv_rates*y[1]*pow(y01im,2.) - g2*y[5]*pow(y01im,2.) - s2*lnv_rates*y[5]*pow(y01im,2.) + 
                        2*(g2 - s2*lnv_rates)*(y[3] - y[7])*y00im*y01re - 
                        2*(g2 + s2*lnv_rates)*(y[2] + y[6])*y00re*y01re - (g2 + s2*lnv_rates)*(y[1] + y[5])*pow(y01re,2.) + 
                        2*g1*(pow(y00im,2.) + pow(y00re,2.) + pow(y01im,2.) + pow(y01re,2.)) + 
                        2*s1*lnv_rates*(pow(y00im,2.) + pow(y00re,2.) + pow(y01im,2.) + pow(y01re,2.)))))/2.;

        yaux[9]         = (_FD_*_MPLTSP_*(-(g0*y[0]*pow(y10im,2.)) + s0*lnv_rates*y[0]*pow(y10im,2.) + g0*y[4]*pow(y10im,2.) - 
                        s0*lnv_rates*y[4]*pow(y10im,2.) - g0*y[0]*pow(y10re,2.) + s0*lnv_rates*y[0]*pow(y10re,2.) + 
                        g0*y[4]*pow(y10re,2.) - s0*lnv_rates*y[4]*pow(y10re,2.) - 2*g0*y[2]*y10im*y11im + 
                        2*s0*lnv_rates*y[2]*y10im*y11im + 2*g0*y[6]*y10im*y11im - 2*s0*lnv_rates*y[6]*y10im*y11im - 
                        2*g0*y[3]*y10re*y11im - 2*s0*lnv_rates*y[3]*y10re*y11im - 2*g0*y[7]*y10re*y11im - 
                        2*s0*lnv_rates*y[7]*y10re*y11im - g0*y[1]*pow(y11im,2.) + s0*lnv_rates*y[1]*pow(y11im,2.) + 
                        g0*y[5]*pow(y11im,2.) - s0*lnv_rates*y[5]*pow(y11im,2.) + 
                        2*(g0 + s0*lnv_rates)*(y[3] + y[7])*y10im*y11re - 2*(g0 - s0*lnv_rates)*(y[2] - y[6])*y10re*y11re - 
                        (g0 - s0*lnv_rates)*(y[1] - y[5])*pow(y11re,2.) + 
                        c_aa*y[9]*(g2*y[0]*pow(y10im,2.) + s2*lnv_rates*y[0]*pow(y10im,2.) + g2*y[4]*pow(y10im,2.) + 
                        s2*lnv_rates*y[4]*pow(y10im,2.) + g2*y[0]*pow(y10re,2.) + s2*lnv_rates*y[0]*pow(y10re,2.) + 
                        g2*y[4]*pow(y10re,2.) + s2*lnv_rates*y[4]*pow(y10re,2.) + 2*g2*y[2]*y10im*y11im + 
                        2*s2*lnv_rates*y[2]*y10im*y11im + 2*g2*y[6]*y10im*y11im + 2*s2*lnv_rates*y[6]*y10im*y11im + 
                        2*g2*y[3]*y10re*y11im - 2*s2*lnv_rates*y[3]*y10re*y11im - 2*g2*y[7]*y10re*y11im + 
                        2*s2*lnv_rates*y[7]*y10re*y11im + g2*y[1]*pow(y11im,2.) + s2*lnv_rates*y[1]*pow(y11im,2.) + 
                        g2*y[5]*pow(y11im,2.) + s2*lnv_rates*y[5]*pow(y11im,2.) - 
                        2*(g2 - s2*lnv_rates)*(y[3] - y[7])*y10im*y11re + 
                        2*(g2 + s2*lnv_rates)*(y[2] + y[6])*y10re*y11re + (g2 + s2*lnv_rates)*(y[1] + y[5])*pow(y11re,2.) - 
                        2*g1*(pow(y10im,2.) + pow(y10re,2.) + pow(y11im,2.) + pow(y11re,2.)) - 
                        2*s1*lnv_rates*(pow(y10im,2.) + pow(y10re,2.) + pow(y11im,2.) + pow(y11re,2.))) - 
                        c_ab*(y[8] + y[10])*(-(g2*y[0]*pow(y10im,2.)) - s2*lnv_rates*y[0]*pow(y10im,2.) - 
                        g2*y[4]*pow(y10im,2.) - s2*lnv_rates*y[4]*pow(y10im,2.) - g2*y[0]*pow(y10re,2.) - 
                        s2*lnv_rates*y[0]*pow(y10re,2.) - g2*y[4]*pow(y10re,2.) - s2*lnv_rates*y[4]*pow(y10re,2.) - 
                        2*g2*y[2]*y10im*y11im - 2*s2*lnv_rates*y[2]*y10im*y11im - 2*g2*y[6]*y10im*y11im - 
                        2*s2*lnv_rates*y[6]*y10im*y11im - 2*g2*y[3]*y10re*y11im + 2*s2*lnv_rates*y[3]*y10re*y11im + 
                        2*g2*y[7]*y10re*y11im - 2*s2*lnv_rates*y[7]*y10re*y11im - g2*y[1]*pow(y11im,2.) - 
                        s2*lnv_rates*y[1]*pow(y11im,2.) - g2*y[5]*pow(y11im,2.) - s2*lnv_rates*y[5]*pow(y11im,2.) + 
                        2*(g2 - s2*lnv_rates)*(y[3] - y[7])*y10im*y11re - 
                        2*(g2 + s2*lnv_rates)*(y[2] + y[6])*y10re*y11re - (g2 + s2*lnv_rates)*(y[1] + y[5])*pow(y11re,2.) + 
                        2*g1*(pow(y10im,2.) + pow(y10re,2.) + pow(y11im,2.) + pow(y11re,2.)) + 
                        2*s1*lnv_rates*(pow(y10im,2.) + pow(y10re,2.) + pow(y11im,2.) + pow(y11re,2.)))))/2.;

        yaux[10]        = (_FD_*_MPLTSP_*(-(g0*y[0]*pow(y20im,2.)) + s0*lnv_rates*y[0]*pow(y20im,2.) + g0*y[4]*pow(y20im,2.) - 
                        s0*lnv_rates*y[4]*pow(y20im,2.) - g0*y[0]*pow(y20re,2.) + s0*lnv_rates*y[0]*pow(y20re,2.) + 
                        g0*y[4]*pow(y20re,2.) - s0*lnv_rates*y[4]*pow(y20re,2.) - 2*g0*y[2]*y20im*y21im + 
                        2*s0*lnv_rates*y[2]*y20im*y21im + 2*g0*y[6]*y20im*y21im - 2*s0*lnv_rates*y[6]*y20im*y21im - 
                        2*g0*y[3]*y20re*y21im - 2*s0*lnv_rates*y[3]*y20re*y21im - 2*g0*y[7]*y20re*y21im - 
                        2*s0*lnv_rates*y[7]*y20re*y21im - g0*y[1]*pow(y21im,2.) + s0*lnv_rates*y[1]*pow(y21im,2.) + 
                        g0*y[5]*pow(y21im,2.) - s0*lnv_rates*y[5]*pow(y21im,2.) + 
                        2*(g0 + s0*lnv_rates)*(y[3] + y[7])*y20im*y21re - 2*(g0 - s0*lnv_rates)*(y[2] - y[6])*y20re*y21re - 
                        (g0 - s0*lnv_rates)*(y[1] - y[5])*pow(y21re,2.) + 
                        c_aa*y[10]*(g2*y[0]*pow(y20im,2.) + s2*lnv_rates*y[0]*pow(y20im,2.) + g2*y[4]*pow(y20im,2.) + 
                        s2*lnv_rates*y[4]*pow(y20im,2.) + g2*y[0]*pow(y20re,2.) + s2*lnv_rates*y[0]*pow(y20re,2.) + 
                        g2*y[4]*pow(y20re,2.) + s2*lnv_rates*y[4]*pow(y20re,2.) + 2*g2*y[2]*y20im*y21im + 
                        2*s2*lnv_rates*y[2]*y20im*y21im + 2*g2*y[6]*y20im*y21im + 2*s2*lnv_rates*y[6]*y20im*y21im + 
                        2*g2*y[3]*y20re*y21im - 2*s2*lnv_rates*y[3]*y20re*y21im - 2*g2*y[7]*y20re*y21im + 
                        2*s2*lnv_rates*y[7]*y20re*y21im + g2*y[1]*pow(y21im,2.) + s2*lnv_rates*y[1]*pow(y21im,2.) + 
                        g2*y[5]*pow(y21im,2.) + s2*lnv_rates*y[5]*pow(y21im,2.) - 
                        2*(g2 - s2*lnv_rates)*(y[3] - y[7])*y20im*y21re + 
                        2*(g2 + s2*lnv_rates)*(y[2] + y[6])*y20re*y21re + (g2 + s2*lnv_rates)*(y[1] + y[5])*pow(y21re,2.) - 
                        2*g1*(pow(y20im,2.) + pow(y20re,2.) + pow(y21im,2.) + pow(y21re,2.)) - 
                        2*s1*lnv_rates*(pow(y20im,2.) + pow(y20re,2.) + pow(y21im,2.) + pow(y21re,2.))) - 
                        c_ab*(y[8] + y[9])*(-(g2*y[0]*pow(y20im,2.)) - s2*lnv_rates*y[0]*pow(y20im,2.) - 
                        g2*y[4]*pow(y20im,2.) - s2*lnv_rates*y[4]*pow(y20im,2.) - g2*y[0]*pow(y20re,2.) - 
                        s2*lnv_rates*y[0]*pow(y20re,2.) - g2*y[4]*pow(y20re,2.) - s2*lnv_rates*y[4]*pow(y20re,2.) - 
                        2*g2*y[2]*y20im*y21im - 2*s2*lnv_rates*y[2]*y20im*y21im - 2*g2*y[6]*y20im*y21im - 
                        2*s2*lnv_rates*y[6]*y20im*y21im - 2*g2*y[3]*y20re*y21im + 2*s2*lnv_rates*y[3]*y20re*y21im + 
                        2*g2*y[7]*y20re*y21im - 2*s2*lnv_rates*y[7]*y20re*y21im - g2*y[1]*pow(y21im,2.) - 
                        s2*lnv_rates*y[1]*pow(y21im,2.) - g2*y[5]*pow(y21im,2.) - s2*lnv_rates*y[5]*pow(y21im,2.) + 
                        2*(g2 - s2*lnv_rates)*(y[3] - y[7])*y20im*y21re - 
                        2*(g2 + s2*lnv_rates)*(y[2] + y[6])*y20re*y21re - (g2 + s2*lnv_rates)*(y[1] + y[5])*pow(y21re,2.) + 
                        2*g1*(pow(y20im,2.) + pow(y20re,2.) + pow(y21im,2.) + pow(y21re,2.)) + 
                        2*s1*lnv_rates*(pow(y20im,2.) + pow(y20re,2.) + pow(y21im,2.) + pow(y21re,2.)))))/2.;
}



/*
smooth Sphaleron freeze out
*/
void spahleron(double t, double* y, double (&sph)){
        double vev_T;
        if (t < _XC_){
                vev_T                   = 0;
        }
        else{
                vev_T                   = sqrt(pow(_VEV_,2.)*(1.-_XC_/t))*sqrt(2.)*t/_TSP_;
        }
        double gamma_diff_broken        = exp(-147.7 + 0.83*_TSP_/t)*pow(_TSP_/t,4.);
        double gamma_diff_sym           = 18.0*pow(_ALPHA_W_,5.)*pow(_TSP_/t,4.);
        double gamma_diff               = min(gamma_diff_sym,gamma_diff_broken);
        double gamma_b                  = gamma_diff*pow(3.,2.)*(869.+333.*pow(vev_T,2.))/(792.+306.*pow(vev_T,2.))*pow(t/_TSP_,3.);
        double chi                      = (4.*(27.*pow(vev_T,2.) + 77.))/(333. * pow(vev_T,2.) + 869.);
        double yb_eq                    = 79./28.*_ASYMFAC_*(y[8]+y[9]+y[10])*chi;
        sph                             = -gamma_b*(y[11] - yb_eq)*_MPL_/t * pow(t/_TSP_,2.);
}
/*
set all vector elements to 0
*/
void vector_reset(double vec[11]){
    for (int i = 0; i < 11; ++i){
            vec[i] = 0;
        }
}





/*
Set the System of differential equation
*/
void fun(double t, double* y, double* ydot, void* data){
        vector<double> ini = *(vector<double> *)data;
        double x_ave;
        x_ave= ini[2];

        if (t > 1.1*x_ave){
                y[2] = 0; y[3] = 0; y[6] = 0; y[7] = 0;
        }

        if (t < 1){
                vector_reset(yaux);
                RHS(t, y, data, yaux);
                for (int i = 0; i< 11; i++){
                        ydot[i] = yaux[i];
                }
                if (t > _XC_){
                        spahleron(t,y,sph);
                        ydot[11]        = sph;
                }
                else{
                        ydot[11]        = 0.;
                }

                sph     = 0.;
                vector_reset(yaux);
        }
        else{
                for (int i = 0; i< 12; i++){
                        ydot[i] = 0.0;
                }  
        }


}