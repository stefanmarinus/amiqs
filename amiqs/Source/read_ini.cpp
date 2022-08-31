#include "../include/inicpp.h"

using namespace std;


vector<double> read_ini(string path){
    static int ini_section;
    vector<double> vars(27);
    vector<string> names;
       
    names= {"param", 
            "hierarchy", "M1", "DM_M", "yukawa", "theta", "delta", "phi", 
            "Lambda", "mu", "ye1", "ymu1", "ytau1", "ye2", "ymu2", "ytau2", "phie2", "phimu2", "phitau2",
            "lnv_rates", "nl_approx", "rates_approx", "c_matrix_approx", "g1_approx", "sph_approx", 
            "info", "safe",
            "zre", "zim"};

    ini::IniFile inif;
    inif.load(path);

    ini_section= 0;
    for(const auto &sectionPair : inif){
        const string &sectionName       = sectionPair.first;
        if (sectionName == "Ini"){
            ini_section++;
            for(const auto &fieldPair : sectionPair.second){
                const string &fieldName     = fieldPair.first;
                const ini::IniField &field  = fieldPair.second;
                for (int i=0; i<names.size(); i++){
                    if (fieldName == names[i]){
                        vars[i] = field.as<double>(); 
                    }
                }
            }
        }
    }


    if (ini_section != 0){
        cout << "[AMIQS] Succesfully read .ini file" << endl;
        return vars;
    }
    else{
        cout << "The .ini file needs to contain the section [Ini]" << endl;
        cout << "It is the only section which will be read" << endl;
        exit(0);
    }
    
}

