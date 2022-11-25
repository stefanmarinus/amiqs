/*
Set the LNV contribution to the thermal mass
*/

class Thermal_Mass_LNV{
private:
    /*
    Set variables used within the class
    */
    vector<double> myval_M;
    vector<double> myval_x;
    vector<double> myval_mthM;

    vector <double> vec_M;
    vector<double> vec_x;
    vector<double> vec_mthM;

    const double *double_M;
    const double *double_x;
    const double *double_mthM;

    int num_M;
    int num_x;
    void *za_alloc;
    double* za;

    gsl_spline2d *spline_mthM;

    const gsl_interp2d_type *InterType;
    gsl_interp_accel *xacc;
    gsl_interp_accel *yacc;




public:
    /*
    Constructor of the class; read files and set the shared interpolation settings
    */
    Thermal_Mass_LNV(void){
        ifstream myfile_M;
        myfile_M.open("../../thermal_mass_LNV_cpp/M_V_N_LNV_average.dat");
        if ( myfile_M.is_open() ) {
            double tmp = 0.0;
            while (myfile_M >> tmp) {
                myval_M.push_back(tmp);
            }
        }
        myfile_M.close();
        for (int i = 0; i < int(myval_M.size()); i++){
            vec_M.push_back(myval_M[i]);
        }

        ifstream myfile_x;
        myfile_x.open("../../thermal_mass_LNV_cpp/x_V_N_LNV_average.dat");
        if ( myfile_x.is_open() ) {
            double tmp = 0.0;
            while (myfile_x >> tmp) {
                myval_x.push_back(tmp);
            }
        }
        myfile_x.close();
        for (int i = 0; i < int(myval_x.size()); i++){
            vec_x.push_back(myval_x[i]);
        }

        ifstream myfile_mthM;
        myfile_mthM.open("../../thermal_mass_LNV_cpp/V_N_LNV_average.dat");
        if ( myfile_mthM.is_open() ) {
            double tmp = 0.0;
            while (myfile_mthM >> tmp) {
                myval_mthM.push_back(tmp);
            }
        }
        myfile_mthM.close();
        for (int i = 0; i < int(myval_mthM.size()); i++){
            vec_mthM.push_back(myval_mthM[i]);
        }



        double_M        = &vec_M[0];
        double_x        = &vec_x[0];
        num_M           = vec_M.size()*sizeof(double_M) / sizeof(double); /* x grid points */
        num_x           = vec_x.size()*sizeof(double_x) / sizeof(double); /* y grid points */
        za_alloc        = malloc(num_M * num_x * sizeof(double));
        za              = (double *)za_alloc;

        InterType       = gsl_interp2d_bilinear;
        xacc            = gsl_interp_accel_alloc();
        yacc            = gsl_interp_accel_alloc();


        spline_mthM     = gsl_spline2d_alloc(InterType, num_M, num_x);
    }



    /**************************************************
    Set interpolation grid
    ***************************************************/

    void set_mthM(void){
        double_mthM    = &vec_mthM[0];
        for (int i = 0; i < int(myval_M.size()); i++){
            for (int j = 0; j < int(myval_x.size()); j++){
                gsl_spline2d_set(spline_mthM, za, i, j, double_mthM[int(myval_x.size())*i+j]);
            }
        }
        /* initialize interpolation */
        gsl_spline2d_init(spline_mthM, double_M, double_x, za, num_M, num_x);
    }



    /**************************************************
    Declare methods returning the interpolated value
    ***************************************************/
    double mthM(double M, double x){
        return gsl_spline2d_eval(spline_mthM, M , x , xacc, yacc);
    }

};


/*
set rates globally, accesible by various modules
*/
Thermal_Mass_LNV thermal_mass_lnv;
/*
call this in the int main() (amiqs.cpp) to set the interpolation grid
*/
void set_thermal_mass_lnv(void){
    thermal_mass_lnv.set_mthM();
}
