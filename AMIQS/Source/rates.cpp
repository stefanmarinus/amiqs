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
    vector<double> myval_lnc_0;
    vector<double> myval_lnc_1;
    vector<double> myval_lnc_2;
    vector<double> myval_lnv_0;
    vector<double> myval_lnv_1;
    vector<double> myval_lnv_2;

    vector <double> vec_M;
    vector<double> vec_x;
    vector<double> vec_lnc_0;
    vector<double> vec_lnc_1;
    vector<double> vec_lnc_2;
    vector<double> vec_lnv_0;
    vector<double> vec_lnv_1;
    vector<double> vec_lnv_2;

    const double *double_M;
    const double *double_x;
    const double *double_lnc_0;
    const double *double_lnc_1;
    const double *double_lnc_2;
    const double *double_lnv_0;
    const double *double_lnv_1;
    const double *double_lnv_2;

    int num_M;
    int num_x;
    void *za_alloc;
    double* za;

    gsl_spline2d *spline_lnc_0;
    gsl_spline2d *spline_lnc_1;
    gsl_spline2d *spline_lnc_2;
    gsl_spline2d *spline_lnv_0;
    gsl_spline2d *spline_lnv_1;
    gsl_spline2d *spline_lnv_2;

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

        ifstream myfile_lnc_0;
        myfile_lnc_0.open("../../rates_cpp/LNC_0_mikko_FD_massive_cpp.dat");
        if ( myfile_lnc_0.is_open() ) {
            double tmp = 0.0;
            while (myfile_lnc_0 >> tmp) {
                myval_lnc_0.push_back(tmp);
            }
        }
        myfile_lnc_0.close();
        for (int i = 0; i < int(myval_lnc_0.size()); i++){
            vec_lnc_0.push_back(myval_lnc_0[i]);
        }

        ifstream myfile_lnc_1;
        myfile_lnc_1.open("../../rates_cpp/LNC_1_mikko_FD_massive_cpp.dat");
        if ( myfile_lnc_1.is_open() ) {
            double tmp = 0.0;
            while (myfile_lnc_1 >> tmp) {
                myval_lnc_1.push_back(tmp);
            }
        }
        myfile_lnc_1.close();
        for (int i = 0; i < int(myval_lnc_1.size()); i++){
            vec_lnc_1.push_back(myval_lnc_1[i]);
        }

        ifstream myfile_lnc_2;
        myfile_lnc_2.open("../../rates_cpp/LNC_2_mikko_FD_massive_cpp.dat");
        if ( myfile_lnc_2.is_open() ) {
            double tmp = 0.0;
            while (myfile_lnc_2 >> tmp) {
                myval_lnc_2.push_back(tmp);
            }
        }
        myfile_lnc_2.close();
        for (int i = 0; i < int(myval_lnc_2.size()); i++){
            vec_lnc_2.push_back(myval_lnc_2[i]);
        }

        ifstream myfile_lnv_0;
        myfile_lnv_0.open("../../rates_cpp/LNV_0_mikko_FD_massive_cpp.dat");
        if ( myfile_lnv_0.is_open() ) {
            double tmp = 0.0;
            while (myfile_lnv_0 >> tmp) {
                myval_lnv_0.push_back(tmp);
            }
        }
        myfile_lnv_0.close();
        for (int i = 0; i < int(myval_lnv_0.size()); i++){
            vec_lnv_0.push_back(myval_lnv_0[i]);
        }

        ifstream myfile_lnv_1;
        myfile_lnv_1.open("../../rates_cpp/LNV_1_mikko_FD_massive_cpp.dat");
        if ( myfile_lnv_1.is_open() ) {
            double tmp = 0.0;
            while (myfile_lnv_1 >> tmp) {
                myval_lnv_1.push_back(tmp);
            }
        }
        myfile_lnv_1.close();
        for (int i = 0; i < int(myval_lnv_1.size()); i++){
            vec_lnv_1.push_back(myval_lnv_1[i]);
        }

        ifstream myfile_lnv_2;
        myfile_lnv_2.open("../../rates_cpp/LNV_2_mikko_FD_massive_cpp.dat");
        if ( myfile_lnv_2.is_open() ) {
            double tmp = 0.0;
            while (myfile_lnv_2 >> tmp) {
                myval_lnv_2.push_back(tmp);
            }
        }
        myfile_lnv_2.close();
        for (int i = 0; i < int(myval_lnv_2.size()); i++){
            vec_lnv_2.push_back(myval_lnv_2[i]);
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


        spline_lnc_0      = gsl_spline2d_alloc(InterType, num_M, num_x);
        spline_lnc_1      = gsl_spline2d_alloc(InterType, num_M, num_x);
        spline_lnc_2      = gsl_spline2d_alloc(InterType, num_M, num_x);
        spline_lnv_0      = gsl_spline2d_alloc(InterType, num_M, num_x);
        spline_lnv_1      = gsl_spline2d_alloc(InterType, num_M, num_x);
        spline_lnv_2      = gsl_spline2d_alloc(InterType, num_M, num_x);
    }

    /*Destructor of the class*/
    // ~Rates(){
    //     gsl_spline2d_free(spline_lnc_0);
    //     gsl_spline2d_free(spline_lnc_1);
    //     gsl_spline2d_free(spline_lnc_2);
    //     gsl_spline2d_free(spline_lnv_0);
    //     gsl_spline2d_free(spline_lnv_1);
    //     gsl_spline2d_free(spline_lnv_2);
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
    void set_lnc_0(void){
        double_lnc_0    = &vec_lnc_0[0];
        for (int i = 0; i < int(myval_M.size()); i++){
            for (int j = 0; j < int(myval_x.size()); j++){
                gsl_spline2d_set(spline_lnc_0, za, i, j, double_lnc_0[int(myval_x.size())*i+j]);
            }
        }
        /* initialize interpolation */
        gsl_spline2d_init(spline_lnc_0, double_M, double_x, za, num_M, num_x);
    }

    /*
    set the interpolation grid for the HC1 rate
    */
    void set_lnc_1(void){
        double_lnc_1    = &vec_lnc_1[0];
        for (int i = 0; i < int(myval_M.size()); i++){
            for (int j = 0; j < int(myval_x.size()); j++){
                gsl_spline2d_set(spline_lnc_1, za, i, j, double_lnc_1[int(myval_x.size())*i+j]);
            }
        }
        /* initialize interpolation */
        gsl_spline2d_init(spline_lnc_1, double_M, double_x, za, num_M, num_x);
    }

    /*
    set the interpolation grid for the HC2 rate
    */
    void set_lnc_2(void){
        double_lnc_2    = &vec_lnc_2[0];
        for (int i = 0; i < int(myval_M.size()); i++){
            for (int j = 0; j < int(myval_x.size()); j++){
                gsl_spline2d_set(spline_lnc_2, za, i, j, double_lnc_2[int(myval_x.size())*i+j]);
            }
        }
        /* initialize interpolation */
        gsl_spline2d_init(spline_lnc_2, double_M, double_x, za, num_M, num_x);
    }

    /*
    set the interpolation grid for the HC2 rate
    */
    void set_lnv_0(void){
        double_lnv_0    = &vec_lnv_0[0];
        for (int i = 0; i < int(myval_M.size()); i++){
            for (int j = 0; j < int(myval_x.size()); j++){
                gsl_spline2d_set(spline_lnv_0, za, i, j, double_lnv_0[int(myval_x.size())*i+j]);
            }
        }
        /* initialize interpolation */
        gsl_spline2d_init(spline_lnv_0, double_M, double_x, za, num_M, num_x);
    }

    /*
    set the interpolation grid for the HC2 rate
    */
    void set_lnv_1(void){
        double_lnv_1    = &vec_lnv_1[0];
        for (int i = 0; i < int(myval_M.size()); i++){
            for (int j = 0; j < int(myval_x.size()); j++){
                gsl_spline2d_set(spline_lnv_1, za, i, j, double_lnv_1[int(myval_x.size())*i+j]);
            }
        }
        /* initialize interpolation */
        gsl_spline2d_init(spline_lnv_1, double_M, double_x, za, num_M, num_x);
    }

    /*
    set the interpolation grid for the HC2 rate
    */
    void set_lnv_2(void){
        double_lnv_2    = &vec_lnv_2[0];
        for (int i = 0; i < int(myval_M.size()); i++){
            for (int j = 0; j < int(myval_x.size()); j++){
                gsl_spline2d_set(spline_lnv_2, za, i, j, double_lnv_2[int(myval_x.size())*i+j]);
            }
        }
        /* initialize interpolation */
        gsl_spline2d_init(spline_lnv_2, double_M, double_x, za, num_M, num_x);
    }




    /**************************************************
    Declare methods returning the interpolated value
    ***************************************************/
    /*
    Function returning the interpolated LNC 1 rate
    */
    double lnc_0(double M, double x){
        return gsl_spline2d_eval(spline_lnc_0, M , x , xacc, yacc);
    }

    /*
    Function returning the interpolated LNC 1 rate
    */
    double lnc_1(double M, double x){
        return gsl_spline2d_eval(spline_lnc_1, M , x , xacc, yacc);
    }

    /*
    Function returning the interpolated LNC 2 rate
    */
    double lnc_2(double M, double x){
        return gsl_spline2d_eval(spline_lnc_2, M , x , xacc, yacc);
    }

    /*
    Function returning the interpolated LNV 0 rate
    */
    double lnv_0(double M, double x){
        return gsl_spline2d_eval(spline_lnv_0, M , x , xacc, yacc);
    }

    /*
    Function returning the interpolated LNV 1 rate
    */
    double lnv_1(double M, double x){
        return gsl_spline2d_eval(spline_lnv_1, M , x , xacc, yacc);
    }

    /*
    Function returning the interpolated LNV 2 rate
    */
    double lnv_2(double M, double x){
        return gsl_spline2d_eval(spline_lnv_2, M , x , xacc, yacc);
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
    rates.set_lnc_0();
    rates.set_lnc_1();
    rates.set_lnc_2();
    rates.set_lnv_0();
    rates.set_lnv_1();
    rates.set_lnv_2();
}
