typedef complex<double> complx;

/*
Class to calculate the Yukawa matrix in the LN basis 
and the rotation to the (real) diagonal majorana mass basis
Functions to calculate the hermitian conjugate, transpose and conjugate are also implemented 
*/

complx YLN[3][2];
complx Y[3][2];
complx YHC[2][3];
complx YT[2][3];
complx YC[3][2];
class Yukawas{
private:
    double theta12;
    double theta13;
    double theta23;
    double dmsolar;
    double dmatm; 
    double rho; 

    double c12;
    double s12;
    double c13;
    double s13;
    double c23;
    double s23;

    complx Upmns[3][3];
    complx W[2][2];
    complx Wangl = polar(1., -_PI_/2.);
    complx Yaux[3][2];
public:
    /*
    Constructor of the class; rotation matrix
    */
    Yukawas(){
        /*
        set the roation matrix
        */
        W[0][0] = 1./sqrt(2.)*Wangl; 
        W[0][1] = 1./sqrt(2.);
        W[1][0] = -1./sqrt(2.)*Wangl;
        W[1][1] = 1./sqrt(2.);
    }
    
    void Yukawa_anaLN(void* param_ini, complx (&YLN) [3][2] ){
        vector<double> param_var;
        param_var = *(vector<double> *)param_ini;
        double ye1, ymu1, ytau1, ye2, ymu2, ytau2, phie2, phimu2, phitau2; 

        ye1= param_var[2]; ymu1= param_var[3]; ytau1= param_var[4]; 
        ye2= param_var[5]; ymu2= param_var[6]; ytau2= param_var[7]; 
        phie2= param_var[8]; phimu2= param_var[9]; phitau2= param_var[10];

        YLN[0][0] = pow(10,ye1); YLN[1][0] = pow(10,ymu1); YLN[2][0] = pow(10,ytau1); 
        YLN[0][1] = pow(10,ye2)*polar(1.,phie2); YLN[1][1] = pow(10,ymu2)*polar(1.,phimu2); YLN[2][1] = pow(10,ytau2)*polar(1.,phitau2);
    }

    void Yukawa_ana(complx YLN[3][2], complx (&Y) [3][2]){
        /* 
        Multiplying matrix YLN and W to get the full Yukawa matrix in the basis
        where the right handed neutrino mass matrix is real and diagonal
        */
        for(int i = 0; i < 3; ++i){
            for(int j = 0; j < 2; ++j){
                for(int k = 0; k < 2; ++k){
                        Yaux[i][j] += YLN[i][k] * W[k][j];
                }
                Y[i][j] = Yaux[i][j];
                // cout << "Y[" << i << "][" << j << "] = " << Y[i][j] << endl; 
                // cout << "Abs(Y[" << i << "][" << j << "]) = " << abs(Y[i][j]) << endl;  
            }
        }
        Yukawa_reset(Yaux);
    }
    
    void YukawaLN(void* param_ini, complx (&YLN) [3][2]){
        vector<double> param_var;
        param_var = *(vector<double> *)param_ini;
        double hierarchy, M1, DM_M, yukawa, theta, delta, phi;
        hierarchy= param_var[0]; M1 = param_var[1]; DM_M = param_var[2]; yukawa = param_var[3]; 
        theta = param_var[4]; delta = param_var[5]; phi = param_var[6];
        
        /*
        get the complex angles
        */
        complx th     = polar(1., theta/2.);
        complx dirac  = polar(1., delta);
        complx maj    = polar(1., phi);

        if (hierarchy == 0){
            theta12     = 33.45 * _PI_/180;
            theta13     = 8.62  * _PI_/180;
            theta23     = 42.1  * _PI_/180;
            dmsolar     = 7.42  * pow(10,-5) * pow(10,-18);
            dmatm       = 2.510 * pow(10,-3)  * pow(10,-18);
            rho         = (sqrt(dmatm) - sqrt(dmsolar))/(sqrt(dmatm) + sqrt(dmsolar));
        }
        else if (hierarchy == 1){
            theta12     = 33.45 * _PI_/180;
            theta13     = 8.61  * _PI_/180;
            theta23     = 49.0  * _PI_/180;
            dmsolar     = 7.42  * pow(10,-5) * pow(10,-18);
            dmatm       = 2.490 * pow(10,-3)  * pow(10,-18);
            rho         = (sqrt(dmatm)-sqrt(dmatm-dmsolar))/(sqrt(dmatm) + sqrt(dmatm-dmsolar));
        }
        else{
            cout << RED << "[AMIQS]" << RESET << " Hierarchy not correctly passed" << endl;
            cout << "=> pass either 0 (NH) or 1 (IH)" << endl;
            cout << "failing in reading the light neutrino params" << endl;
            exit(0);
        }
        c12 = cos(theta12); s12 = sin(theta12);
        c13 = cos(theta13); s13 = sin(theta13);
        c23 = cos(theta23); s23 = sin(theta23);

        /*
        change to linear scale and calculate related parameters
        */
        M1             = pow(10.,M1);
        DM_M           = pow(10.,DM_M);
        yukawa         = pow(10.,yukawa);
        double M2      = M1 * (DM_M + 1.);
        double Lambda  = 1./2. * (M1 + M2);
        double eta     = 1./2. * (M2 - M1);


        /*
        set up the Upmn matrix
        */
        Upmns[0][0] =  c12*c13;
        Upmns[0][1] =  s12*c13*maj;
        Upmns[0][2] =  s13*conj(dirac);
        Upmns[1][0] = (-s12*c23 - c12*s13*s23*dirac);
        Upmns[1][1] =  (c12*c23 - s12*s13*s23*dirac)*maj;
        Upmns[1][2] =  c13*s23;
        Upmns[2][0] =  (s12*s23 - c12*s13*c23*dirac);
        Upmns[2][1] = (-c12*s23 - s12*s13*c23*dirac)*maj;
        Upmns[2][2] =  c13*c23;

        //  for(int i = 0; i < 3; ++i){
        //    for(int j = 0; j < 3; ++j){
        //      cout << "Upmns[" << i << "][" << j << "] = " << Upmns[i][j] << endl; 
        //    }
        //  }


        /*
        the y' scale
        */
       double yp;
        if (hierarchy == 0){
            yp  = 1/yukawa * Lambda/2./pow(_VEV_,2) * (sqrt(dmatm) + sqrt(dmsolar));
        }
        else if (hierarchy == 1){
            yp  = 1/yukawa * Lambda/2./pow(_VEV_,2) * (sqrt(dmatm) + sqrt(dmatm - dmsolar));
        } 
        else{
            cout << RED << "[AMIQS]" << RESET << " Hierarchy not correctly passed" << endl;
            cout << "=> pass either 0 (NH) or 1 (IH)" << endl;
            cout << "failing in determining the y' scale" << endl;
            exit(0);
        }

        /*
        if the LN constraint is not fulfilled we just return a 0 matrix
        */
        if (yp > yukawa/10.){
            // cout << "[AMIQS] LN constraint not fulfilled with ratio " << scientific << yp/yukawa << endl;
            // cout << "[AMIQS] LN constraint not fulfilled" << endl;
            YLN[0][0]= 0; YLN[1][0]= 0; YLN[2][0]=0; YLN[0][1]=0; YLN[1][1]=0; YLN[2][1]=0;
        }
        else{
            /*
            set up the Yukawa matrix in the LN parametrization
            */
            if (hierarchy == 0){
                YLN[0][0]   = yukawa/sqrt(2)*conj(th)*(sqrt(1+rho)*conj(Upmns[0][2]) + sqrt(1-rho)*conj(Upmns[0][1]));
                YLN[1][0]   = yukawa/sqrt(2)*conj(th)*(sqrt(1+rho)*conj(Upmns[1][2]) + sqrt(1-rho)*conj(Upmns[1][1]));
                YLN[2][0]   = yukawa/sqrt(2)*conj(th)*(sqrt(1+rho)*conj(Upmns[2][2]) + sqrt(1-rho)*conj(Upmns[2][1]));
                YLN[0][1]   = yp/sqrt(2)*th*(sqrt(1+rho)*conj(Upmns[0][2]) - sqrt(1-rho)*conj(Upmns[0][1])) + yukawa/sqrt(2)*eta/2/Lambda*conj(th)*(sqrt(1+rho)*conj(Upmns[0][2]) + sqrt(1-rho)*conj(Upmns[0][1]));
                YLN[1][1]   = yp/sqrt(2)*th*(sqrt(1+rho)*conj(Upmns[1][2]) - sqrt(1-rho)*conj(Upmns[1][1])) + yukawa/sqrt(2)*eta/2/Lambda*conj(th)*(sqrt(1+rho)*conj(Upmns[1][2]) + sqrt(1-rho)*conj(Upmns[1][1]));
                YLN[2][1]   = yp/sqrt(2)*th*(sqrt(1+rho)*conj(Upmns[2][2]) - sqrt(1-rho)*conj(Upmns[2][1])) + yukawa/sqrt(2)*eta/2/Lambda*conj(th)*(sqrt(1+rho)*conj(Upmns[2][2]) + sqrt(1-rho)*conj(Upmns[2][1]));
            }
            else if (hierarchy == 1){
                YLN[0][0]   = yukawa/sqrt(2)*conj(th)*(sqrt(1+rho)*conj(Upmns[0][1]) + sqrt(1-rho)*conj(Upmns[0][0]));
                YLN[1][0]   = yukawa/sqrt(2)*conj(th)*(sqrt(1+rho)*conj(Upmns[1][1]) + sqrt(1-rho)*conj(Upmns[1][0]));
                YLN[2][0]   = yukawa/sqrt(2)*conj(th)*(sqrt(1+rho)*conj(Upmns[2][1]) + sqrt(1-rho)*conj(Upmns[2][0]));
                YLN[0][1]   = yp/sqrt(2)*th*(sqrt(1+rho)*conj(Upmns[0][1]) - sqrt(1-rho)*conj(Upmns[0][0])) + yukawa/sqrt(2)*eta/2/Lambda*conj(th)*(sqrt(1+rho)*conj(Upmns[0][1]) + sqrt(1-rho)*conj(Upmns[0][0]));
                YLN[1][1]   = yp/sqrt(2)*th*(sqrt(1+rho)*conj(Upmns[1][1]) - sqrt(1-rho)*conj(Upmns[1][0])) + yukawa/sqrt(2)*eta/2/Lambda*conj(th)*(sqrt(1+rho)*conj(Upmns[1][1]) + sqrt(1-rho)*conj(Upmns[1][0]));
                YLN[2][1]   = yp/sqrt(2)*th*(sqrt(1+rho)*conj(Upmns[2][1]) - sqrt(1-rho)*conj(Upmns[2][0])) + yukawa/sqrt(2)*eta/2/Lambda*conj(th)*(sqrt(1+rho)*conj(Upmns[2][1]) + sqrt(1-rho)*conj(Upmns[2][0]));
            }
            else{
               cout << RED << "[AMIQS]" << RESET << " Yukawa matrix not calculated due to unknown hierarchy" << endl;
               cout << "=> pass either 0 (NH) or 1 (IH)" << endl;
               exit(0);
            }

        //    for(int i = 0; i < 3; ++i){
        //        for(int j = 0; j < 2; ++j){
        //            cout << "YLN[" << i << "][" << j << "] = " << YLN[i][j] << endl; 
        //           //  cout << "Abs(YLN[" << i << "][" << j << "]) = " << abs(YLN[i][j]) << endl;  
        //        }
        //    }
        }
    }

    void Yukawa(complx YLN[3][2], complx (&Y) [3][2]){
        /* 
        Multiplying matrix YLN and W to get the full Yukawa matrix in the basis
        where the right handed neutrino mass matrix is real and diagonal
        */
        for(int i = 0; i < 3; ++i){
            for(int j = 0; j < 2; ++j){
                for(int k = 0; k < 2; ++k){
                        Yaux[i][j] += YLN[i][k] * W[k][j];
                }
                Y[i][j] = Yaux[i][j];
                // cout << "Y[" << i << "][" << j << "] = " << Y[i][j] << endl; 
                // cout << "Abs(Y[" << i << "][" << j << "]) = " << abs(Y[i][j]) << endl;  
            }
        }
        Yukawa_reset(Yaux);
    }

    void YukawaHC(complx Y[3][2],  complx (&YHC) [2][3]){
        for(int i = 0; i < 3; ++i){
            for(int j = 0; j < 2; ++j){
                YHC[i][j] = conj(Y[j][i]);
            }
        }
    }

    void YukawaT(complx Y[3][2],  complx (&YT) [2][3]){
        for(int i = 0; i < 3; ++i){
            for(int j = 0; j < 2; ++j){
                YT[i][j] = Y[j][i];
            }
        }
    }

    void YukawaC(complx Y[3][2],  complx (&YC) [3][2]){
        for(int i = 0; i < 3; ++i){
            for(int j = 0; j < 2; ++j){
                YC[i][j] = conj(Y[j][j]);
            }
        }
    }

    void Yukawa_reset(complx mat[3][2]){
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 2; ++j){
                mat[i][j] = 0;
            }
        }
    }
};

/*
initialize Yukawa matrix globally
*/
Yukawas Yukawa;

