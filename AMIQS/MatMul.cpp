/*
Class for some convinient matrix multiplications
In particular we implement multiplications of:
    2x2 . 2x2
    3x2 . 2x2
    3x2 . 2x3
    2x3 . 3x2
    2x2 . 2x3
    3x3 . 3x2
    3x3 . 3x3
And the (anti-) commutator for 2x2 matrices
*/
class MatMul{
private:
    complx matres_22_aux[2][2];
    complx matres_32_aux[3][2];
    complx matres_33_aux[3][3];
    complx matres_23_aux[2][3];
public:
    void matmul_22_22(complx mat1[2][2], complx mat2[2][2], complx mat_res [2][2]){
        for (int i = 0; i < 2; ++i){
            for (int j = 0; j < 2; ++j){
                matres_22_aux[i][j] = 0;
            }
        }
        for(int i = 0; i < 2; ++i){
            for(int j = 0; j < 2; ++j){
                for(int k = 0; k < 2; ++k){
                        matres_22_aux[i][j] += mat1[i][k] * mat2[k][j];
                }
                mat_res[i][j] = matres_22_aux[i][j];
            }
        }
        for (int i = 0; i < 2; ++i){
            for (int j = 0; j < 2; ++j){
                matres_22_aux[i][j] = 0;
            }
        }
    }

    void matmul_32_22(complx mat1[3][2], complx mat2[2][2], complx mat_res [3][2]){
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 2; ++j){
                matres_32_aux[i][j] = 0;
            }
        }
        for(int i = 0; i < 3; ++i){
            for(int j = 0; j < 2; ++j){
                for(int k = 0; k < 2; ++k){
                        matres_32_aux[i][j] += mat1[i][k] * mat2[k][j];
                }
                mat_res[i][j] = matres_32_aux[i][j];
            }
        }
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 2; ++j){
                matres_32_aux[i][j] = 0;
            }
        }
    }

    void matmul_32_23(complx mat1[3][2], complx mat2[2][3], complx mat_res [3][3]){
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                matres_33_aux[i][j] = 0;
            }
        }
        for(int i = 0; i < 3; ++i){
            for(int j = 0; j < 3; ++j){
                for(int k = 0; k < 2; ++k){
                        matres_33_aux[i][j] += mat1[i][k] * mat2[k][j];
                }
                mat_res[i][j] = matres_33_aux[i][j];
            }
        }
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                matres_33_aux[i][j] = 0;
            }
        }
    }

    void matmul_23_32(complx mat1[2][3], complx mat2[3][2], complx mat_res [2][2]){
        for (int i = 0; i < 2; ++i){
            for (int j = 0; j < 2; ++j){
                matres_22_aux[i][j] = 0;
            }
        }
        for(int i = 0; i < 2; ++i){
            for(int j = 0; j < 2; ++j){
                for(int k = 0; k < 3; ++k){
                        matres_22_aux[i][j] += mat1[i][k] * mat2[k][j];
                }
                mat_res[i][j] = matres_22_aux[i][j];
            }
        }
        for (int i = 0; i < 2; ++i){
            for (int j = 0; j < 2; ++j){
                matres_22_aux[i][j] = 0;
            }
        }
    }

    void matmul_22_23(complx mat1[2][2], complx mat2[2][3], complx mat_res [2][3]){
        for (int i = 0; i < 2; ++i){
            for (int j = 0; j < 3; ++j){
                matres_23_aux[i][j] = 0;
            }
        }
        for(int i = 0; i < 2; ++i){
            for(int j = 0; j < 3; ++j){
                for(int k = 0; k < 2; ++k){
                        matres_23_aux[i][j] += mat1[i][k] * mat2[k][j];
                }
                mat_res[i][j] = matres_23_aux[i][j];
            }
        }
        for (int i = 0; i < 2; ++i){
            for (int j = 0; j < 3; ++j){
                matres_23_aux[i][j] = 0;
            }
        }
    }

    void matmul_33_32(complx mat1[3][3], complx mat2[3][2], complx mat_res [3][2]){
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 2; ++j){
                matres_32_aux[i][j] = 0;
            }
        }
        for(int i = 0; i < 3; ++i){
            for(int j = 0; j < 2; ++j){
                for(int k = 0; k < 3; ++k){
                        matres_32_aux[i][j] += mat1[i][k] * mat2[k][j];
                }
                mat_res[i][j] = matres_32_aux[i][j];
            }
        }
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 2; ++j){
                matres_32_aux[i][j] = 0;
            }
        }
    }

    void matmul_33_33(complx mat1[3][3], complx mat2[3][3], complx mat_res [3][3]){
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                matres_33_aux[i][j] = 0;
            }
        }
        for(int i = 0; i < 3; ++i){
            for(int j = 0; j < 3; ++j){
                for(int k = 0; k < 3; ++k){
                        matres_33_aux[i][j] += mat1[i][k] * mat2[k][j];
                }
                mat_res[i][j] = matres_33_aux[i][j];
            }
        }
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                matres_33_aux[i][j] = 0;
            }
        }
    }


    void anticom_22_22(complx mat1[2][2], complx mat2[2][2], complx mat_res [2][2]){
        for (int i = 0; i < 2; ++i){
            for (int j = 0; j < 2; ++j){
                matres_22_aux[i][j] = 0;
            }
        }
        for(int i = 0; i < 2; ++i){
            for(int j = 0; j < 2; ++j){
                for(int k = 0; k < 2; ++k){
                        matres_22_aux[i][j] += mat1[i][k] * mat2[k][j] + mat2[i][k] * mat1[k][j];
                }
                mat_res[i][j] = matres_22_aux[i][j];
            }
        }
        for (int i = 0; i < 2; ++i){
            for (int j = 0; j < 2; ++j){
                matres_22_aux[i][j] = 0;
            }
        }
    }

    void com_22_22(complx mat1[2][2], complx mat2[2][2], complx mat_res [2][2]){
        for (int i = 0; i < 2; ++i){
            for (int j = 0; j < 2; ++j){
                matres_22_aux[i][j] = 0;
            }
        }
        for(int i = 0; i < 2; ++i){
            for(int j = 0; j < 2; ++j){
                for(int k = 0; k < 2; ++k){
                        matres_22_aux[i][j] += mat1[i][k] * mat2[k][j] - mat2[i][k] * mat1[k][j];
                }
                mat_res[i][j] = matres_22_aux[i][j];
            }
        }
        for (int i = 0; i < 2; ++i){
            for (int j = 0; j < 2; ++j){
                matres_22_aux[i][j] = 0;
            }
        }
    }


    void mat22_free(complx mat[2][2]){
        for (int i = 0; i < 2; ++i){
            for (int j = 0; j < 2; ++j){
                mat[i][j] = 0;
            }
        } 
    }

    void mat32_free(complx mat[3][2]){
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 2; ++j){
                mat[i][j] = 0;
            }
        } 
    }

    void mat23_free(complx mat[2][3]){
        for (int i = 0; i < 2; ++i){
            for (int j = 0; j < 3; ++j){
                mat[i][j] = 0;
            }
        } 
    }

    void mat33_free(complx mat[3][3]){
        for (int i = 0; i < 3; ++i){
            for (int j = 0; j < 3; ++j){
                mat[i][j] = 0;
            }
        } 
    }
};


/*
initialize globally matrix multiplication class
*/
MatMul matmul;




