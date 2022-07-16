/*
Set the interactions rates
*/

class Rates{
private:
    /*
    Set variables used within the class
    */
    vector<double> myval_M;
    vector<double> myval_x;
    vector<double> myval_hc0;
    vector<double> myval_hc1;
    vector<double> myval_hc2;
    vector<double> myval_hv0;
    vector<double> myval_hv1;
    vector<double> myval_hv2;

    vector <double> vec_M;
    vector<double> vec_x;
    vector<double> vec_hc0;
    vector<double> vec_hc1;
    vector<double> vec_hc2;
    vector<double> vec_hv0;
    vector<double> vec_hv1;
    vector<double> vec_hv2;

    const double *double_M;
    const double *double_x;
    const double *double_hc0;
    const double *double_hc1;
    const double *double_hc2;
    const double *double_hv0;
    const double *double_hv1;
    const double *double_hv2;

    int num_M;
    int num_x;
    void *za_alloc;
    double* za;

    gsl_spline2d *spline_hc0;
    gsl_spline2d *spline_hc1;
    gsl_spline2d *spline_hc2;
    gsl_spline2d *spline_hv0;
    gsl_spline2d *spline_hv1;
    gsl_spline2d *spline_hv2;

    const gsl_interp2d_type *InterType;
    gsl_interp_accel *xacc;
    gsl_interp_accel *yacc;





public:
    /*
    Constructor of the class; read files of the rates and set the shared interpolation settings
    */
    Rates(void){
        ifstream myfile_M;
        myfile_M.open("../../rates_cpp/M_mikko_FD_massive_cpp.dat");
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
        myfile_x.open("../../rates_cpp/x_mikko_FD_massive_cpp.dat");
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

        ifstream myfile_hc0;
        myfile_hc0.open("../../rates_cpp/HC0_mikko_FD_massive_cpp.dat");
        if ( myfile_hc0.is_open() ) {
            double tmp = 0.0;
            while (myfile_hc0 >> tmp) {
                myval_hc0.push_back(tmp);
            }
        }
        myfile_hc0.close();
        for (int i = 0; i < int(myval_hc0.size()); i++){
            vec_hc0.push_back(myval_hc0[i]);
        }

        ifstream myfile_hc1;
        myfile_hc1.open("../../rates_cpp/HC1_mikko_FD_massive_cpp.dat");
        if ( myfile_hc1.is_open() ) {
            double tmp = 0.0;
            while (myfile_hc1 >> tmp) {
                myval_hc1.push_back(tmp);
            }
        }
        myfile_hc1.close();
        for (int i = 0; i < int(myval_hc1.size()); i++){
            vec_hc1.push_back(myval_hc1[i]);
        }

        ifstream myfile_hc2;
        myfile_hc2.open("../../rates_cpp/HC2_mikko_FD_massive_cpp.dat");
        if ( myfile_hc2.is_open() ) {
            double tmp = 0.0;
            while (myfile_hc2 >> tmp) {
                myval_hc2.push_back(tmp);
            }
        }
        myfile_hc2.close();
        for (int i = 0; i < int(myval_hc2.size()); i++){
            vec_hc2.push_back(myval_hc2[i]);
        }

        ifstream myfile_hv0;
        myfile_hv0.open("../../rates_cpp/HV0_mikko_FD_massive_cpp.dat");
        if ( myfile_hv0.is_open() ) {
            double tmp = 0.0;
            while (myfile_hv0 >> tmp) {
                myval_hv0.push_back(tmp);
            }
        }
        myfile_hv0.close();
        for (int i = 0; i < int(myval_hv0.size()); i++){
            vec_hv0.push_back(myval_hv0[i]);
        }

        ifstream myfile_hv1;
        myfile_hv1.open("../../rates_cpp/HV1_mikko_FD_massive_cpp.dat");
        if ( myfile_hv1.is_open() ) {
            double tmp = 0.0;
            while (myfile_hv1 >> tmp) {
                myval_hv1.push_back(tmp);
            }
        }
        myfile_hv1.close();
        for (int i = 0; i < int(myval_hv1.size()); i++){
            vec_hv1.push_back(myval_hv1[i]);
        }

        ifstream myfile_hv2;
        myfile_hv2.open("../../rates_cpp/HV2_mikko_FD_massive_cpp.dat");
        if ( myfile_hv2.is_open() ) {
            double tmp = 0.0;
            while (myfile_hv2 >> tmp) {
                myval_hv2.push_back(tmp);
            }
        }
        myfile_hv2.close();
        for (int i = 0; i < int(myval_hv2.size()); i++){
            vec_hv2.push_back(myval_hv2[i]);
        }


        double_M        = &vec_M[0];
        double_x        = &vec_x[0];
        num_M           = vec_M.size()*sizeof(double_M) / sizeof(double); /* x grid points */
        num_x           =  vec_x.size()*sizeof(double_x) / sizeof(double); /* y grid points */
        za_alloc        = malloc(num_M * num_x * sizeof(double));
        za              = (double *)za_alloc;

        InterType       = gsl_interp2d_bilinear;
        xacc            = gsl_interp_accel_alloc();
        yacc            = gsl_interp_accel_alloc();


        spline_hc0      = gsl_spline2d_alloc(InterType, num_M, num_x);
        spline_hc1      = gsl_spline2d_alloc(InterType, num_M, num_x);
        spline_hc2      = gsl_spline2d_alloc(InterType, num_M, num_x);
        spline_hv0      = gsl_spline2d_alloc(InterType, num_M, num_x);
        spline_hv1      = gsl_spline2d_alloc(InterType, num_M, num_x);
        spline_hv2      = gsl_spline2d_alloc(InterType, num_M, num_x);
    }

    /*Destructor of the class*/
    // ~Rates(){
    //     gsl_spline2d_free(spline_hc0);
    //     gsl_spline2d_free(spline_hc1);
    //     gsl_spline2d_free(spline_hc2);
    //     gsl_spline2d_free(spline_hv0);
    //     gsl_spline2d_free(spline_hv1);
    //     gsl_spline2d_free(spline_hv2);
    //     gsl_interp_accel_free(xacc);
    //     gsl_interp_accel_free(yacc);
    //     free(za);
    //     cout << "Momory succesfully freed" << endl;
    // }



    /**************************************************
    Set interpolation grid
    ***************************************************/
    /*
    set the interpolation grid for the HC0 rate
    */
    void set_hc0(void){
        double_hc0    = &vec_hc0[0];
        for (int i = 0; i < int(myval_M.size()); i++){
            for (int j = 0; j < int(myval_x.size()); j++){
                gsl_spline2d_set(spline_hc0, za, i, j, double_hc0[int(myval_x.size())*i+j]);
            }
        }
        /* initialize interpolation */
        gsl_spline2d_init(spline_hc0, double_M, double_x, za, num_M, num_x);
    }

    /*
    set the interpolation grid for the HC1 rate
    */
    void set_hc1(void){
        double_hc1    = &vec_hc1[0];
        for (int i = 0; i < int(myval_M.size()); i++){
            for (int j = 0; j < int(myval_x.size()); j++){
                gsl_spline2d_set(spline_hc1, za, i, j, double_hc1[int(myval_x.size())*i+j]);
            }
        }
        /* initialize interpolation */
        gsl_spline2d_init(spline_hc1, double_M, double_x, za, num_M, num_x);
    }

    /*
    set the interpolation grid for the HC2 rate
    */
    void set_hc2(void){
        double_hc2    = &vec_hc2[0];
        for (int i = 0; i < int(myval_M.size()); i++){
            for (int j = 0; j < int(myval_x.size()); j++){
                gsl_spline2d_set(spline_hc2, za, i, j, double_hc2[int(myval_x.size())*i+j]);
            }
        }
        /* initialize interpolation */
        gsl_spline2d_init(spline_hc2, double_M, double_x, za, num_M, num_x);
    }

    /*
    set the interpolation grid for the HC2 rate
    */
    void set_hv0(void){
        double_hv0    = &vec_hv0[0];
        for (int i = 0; i < int(myval_M.size()); i++){
            for (int j = 0; j < int(myval_x.size()); j++){
                gsl_spline2d_set(spline_hv0, za, i, j, double_hv0[int(myval_x.size())*i+j]);
            }
        }
        /* initialize interpolation */
        gsl_spline2d_init(spline_hv0, double_M, double_x, za, num_M, num_x);
    }

    /*
    set the interpolation grid for the HC2 rate
    */
    void set_hv1(void){
        double_hv1    = &vec_hv1[0];
        for (int i = 0; i < int(myval_M.size()); i++){
            for (int j = 0; j < int(myval_x.size()); j++){
                gsl_spline2d_set(spline_hv1, za, i, j, double_hv1[int(myval_x.size())*i+j]);
            }
        }
        /* initialize interpolation */
        gsl_spline2d_init(spline_hv1, double_M, double_x, za, num_M, num_x);
    }

    /*
    set the interpolation grid for the HC2 rate
    */
    void set_hv2(void){
        double_hv2    = &vec_hv2[0];
        for (int i = 0; i < int(myval_M.size()); i++){
            for (int j = 0; j < int(myval_x.size()); j++){
                gsl_spline2d_set(spline_hv2, za, i, j, double_hv2[int(myval_x.size())*i+j]);
            }
        }
        /* initialize interpolation */
        gsl_spline2d_init(spline_hv2, double_M, double_x, za, num_M, num_x);
    }




    /**************************************************
    Declare methods returning the interpolated value
    ***************************************************/
    /*
    Function returning the interpolated HC0 rate
    */
    double hc0(double M, double x){
        return gsl_spline2d_eval(spline_hc0, M , x , xacc, yacc);
    }

    /*
    Function returning the interpolated HC1 rate
    */
    double hc1(double M, double x){
        return gsl_spline2d_eval(spline_hc1, M , x , xacc, yacc);
    }

    /*
    Function returning the interpolated HC2 rate
    */
    double hc2(double M, double x){
        return gsl_spline2d_eval(spline_hc2, M , x , xacc, yacc);
    }

    /*
    Function returning the interpolated HV0 rate
    */
    double hv0(double M, double x){
        return gsl_spline2d_eval(spline_hv0, M , x , xacc, yacc);
    }

    /*
    Function returning the interpolated HV1 rate
    */
    double hv1(double M, double x){
        return gsl_spline2d_eval(spline_hv1, M , x , xacc, yacc);
    }

    /*
    Function returning the interpolated HV2 rate
    */
    double hv2(double M, double x){
        return gsl_spline2d_eval(spline_hv2, M , x , xacc, yacc);
    }
};


/*
set rates globally, accesible by various modules
*/
Rates rates;
/*
call this in the int main() (amiqs.cpp) to set the interpolation grid
*/
void set_rates(void){
    rates.set_hc0();
    rates.set_hc1();
    rates.set_hc2();
    rates.set_hv0();
    rates.set_hv1();
    rates.set_hv2();
}
