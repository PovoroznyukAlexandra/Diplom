#include "diplom_func.h"

///////   Комплексная экспонента   ///////
complex<double> Comp_Exp(double fi) {
    complex<double> z(cos(fi), sin(fi));
    return z;
}


///////   Модуль комплексного вектора   ///////
double Comp_abs(complex<double>* h, int n = 3) {
    double comp_abs = 0;
    for (int i = 0; i < n; i++) {
        comp_abs += abs(h[i]) * abs(h[i]);
    }
    return sqrt(comp_abs);
}


///////   Скалярное произведение 1   ///////
complex<double> Scal_Prod(complex<double>* x, complex<double>* y, int n = 3) {
    complex<double> res = 0;
    for (int i = 0; i < n; i++) {
        res += x[i]*y[i];
    }
    return res;
}


///////   Скалярное произведение 2   ///////
complex<double> Scal_Prod(complex<double>* x, double* y, int n = 3) {
    complex<double> res = 0;
    for (int i = 0; i < n; i++) {
        res += x[i]*y[i];
    }
    return res;
}


///////   Скалярное произведение 3   ///////
double Scal_Prod(double* x, double* y, int n = 3) {
    double res = 0;
    for (int i = 0; i < n; i++) {
        res += x[i]*y[i];
    }
    return res;
}

//---------------------------------------- Структура -----------------------------------------------------//

Data::Data(int* N, double* L, double* A, int n) {
        this->N[0] = N[0]; this->N[1] = N[1]; this->N[2] = N[2]; 
        this->L[0] = L[0]; this->L[1] = L[1]; this->L[2] = L[2];
        this->A[0] = A[0]; this->A[1] = A[1]; this->A[2] = A[2];
        this->n = n;
        this->V_cell = (L[0] * L[1] * L[2])/(N[0] * N[1] * N[2]);
}

//---------------------------------------- Заданная функция -----------------------------------------------------//

int Chi(double *x, Data &data){
    double y[3];
    for (int i = 0; i < 3; i++){
        y[i] = x[i] - data.O_sphere[i];
    }
    if (Scal_Prod(y, y) <= data.R*data.R)
        return 1;
    else
        return 0;
}

///////   Функция Эпсилон   ///////
double Epsilon(double *x) {
    return 4;
}

/*
double Eps_on_cell(double** V_x, Data &data) {
    //cout << "Check_4_Eps_on_cell:" << endl;
    // V_x[0] = Ax, V_x[1] = Bx, V_x[2] = Cx, V_x[3] = Dx
    double x_centre[3];
    //x_centre = A + AB/2 + AC/2 + AD/2
    x_centre[0] = V_x[0][0] + (V_x[1][0] - V_x[0][0])/2; //x_centre[0] = A[0] + (B[0] - A[0])/2
    x_centre[1] = V_x[0][1] + (V_x[2][1] - V_x[0][1])/2; //x_centre[1] = A[1] + (C[1] - A[1])/2
    x_centre[2] = V_x[0][2] + (V_x[3][2] - V_x[0][2])/2; //x_centre[2] = A[2] + (D[2] - A[2])/2
    
    double eps = (Epsilon(x_centre)-1) * Chi(x_centre, data) + 1;
    return eps;
}
*/


void Read_Eps_on_cell(double *E_cell, Data &data) {
    cout << endl;
    cout << "Check_3_Read_Eps_on_cell:" << endl;
    ifstream in("/Users/deti/Desktop/ВТМ/Курсач/Dielectric_scattering/Epsilon.txt"); // окрываем файл для чтения
    double ch;
    int cnt = 0, ind = 0;
    while (in >> ch) {
        if (ind < data.N[0]*data.N[1]*data.N[2]){
            E_cell[ind] = ch;
        }
        //cout << E_cell[ind] << endl;
        cnt++;
        ind++;
    }
    cout << "Чисел в файле Epsilon.txt: " << cnt << endl;
    in.close();
    
    if (cnt == (data.N[0]*data.N[1]*data.N[2])){
        cout << "File Epsilon.txt has been read correctly" << endl;
    }
    else {
        cout << "File Epsilon.txt has been read with errors" << endl;
        return;
    }
    cout << endl;
}


double Theta(double x) {
    double res;
    if (x < 1) {
        res = 3*pow(x, 2) - 2*pow(x, 3);
    }
    else {
        res = 1.0;
    }
    return res;
}

complex<double> Func(Data data, int h, double* y) {
    double v[3];
    complex<double> res;
    complex<double> zero(0, 0);
    for (int i = 0; i < 3; i++) {
        v[i] = data.x[i] - y[i];
    }
    double r = sqrt(Scal_Prod(v, v));
    if (r <= (1e-12)) {
        res = zero;
        return res;
    }
    data.eps = 2*h/data.n; //?

    res = Comp_Exp(K*r)/(4*M_PI*r);
    res *= Theta(r/data.eps);
    return res;
}

//---------------------------------------- Поверхностный интеграл -----------------------------------------------------//

// Sigma_y[0] = Ay, Sigma_y[1] = By, Sigma_y[2] = Dy
complex<double> I_inner_Sigma(double** Sigma_y, Data &data) {
    double p[3];
    double q[3];
    for (int i = 0; i < 3; i++) {
        p[i] = Sigma_y[1][i] - Sigma_y[0][i]; // p = By-Ay
        q[i] = Sigma_y[2][i] - Sigma_y[0][i]; // q = Dy-Ay
    }
    double p_len = sqrt(Scal_Prod(p, p));
    double q_len = sqrt(Scal_Prod(q, q));
    double h = max(p_len/data.n, q_len/data.n);
    
    complex<double> Int(0, 0);
    double y[3];
    for (int i = 0; i < data.n; i++) {
        for (int j = 0; j < data.n; j++) {
            y[0] = Sigma_y[0][0] + (i + 0.5) * p[0]/data.n + (j + 0.5) * q[0]/data.n;
            y[1] = Sigma_y[0][1] + (i + 0.5) * p[1]/data.n + (j + 0.5) * q[1]/data.n;
            y[2] = Sigma_y[0][2] + (i + 0.5) * p[2]/data.n + (j + 0.5) * q[2]/data.n;
            Int += Func(data, h, y);
        }
    }
    double mu = p_len * q_len;
    Int *= mu/(data.n * data.n);
    return Int; // Int * mu(S) / n^2
}

// typedef complex<double> (*func)(Data, int, double*);
// S = ABD - прямоугольник 
// Sigma_x[0] = Ax, Sigma_x[1] = Bx, Sigma_x[2] = Dx
// Sigma_y[0] = Ay, Sigma_y[1] = By, Sigma_y[2] = Dy
complex<double> Integral_Sigma(double** Sigma_x, double** Sigma_y, Data &data) {
    double p[3];
    double q[3];
    for (int i = 0; i < 3; i++) {
        p[i] = Sigma_x[1][i] - Sigma_x[0][i]; //p[i] = Bx[i] - Ax[i];
        q[i] = Sigma_x[2][i] - Sigma_x[0][i]; //q[i] = Dx[i] - Ax[i];
    }
    double p_len = sqrt(Scal_Prod(p, p));
    double q_len = sqrt(Scal_Prod(q, q));
    
    complex<double> Int(0, 0);
    double y[3];
    for (int i = 0; i < data.n; i++) {
        for (int j = 0; j < data.n; j++) {
            data.x[0] = Sigma_x[0][0] + (i + 0.5) * p[0]/data.n + (j + 0.5) * q[0]/data.n; 
            data.x[1] = Sigma_x[0][1] + (i + 0.5) * p[1]/data.n + (j + 0.5) * q[1]/data.n; 
            data.x[2] = Sigma_x[0][2] + (i + 0.5) * p[2]/data.n + (j + 0.5) * q[2]/data.n;
            Int += I_inner_Sigma(Sigma_y, data);
        }
    }
    double mu = p_len * q_len;
    Int *= mu/(data.n * data.n);
    return Int; //Int * mu(S) / n^2
}

//---------------------------------------- Объемный интеграл -----------------------------------------------------//

// V_y[0] = Ay, V_y[1] = By, V_y[2] = Cy, V_y[3] = Dy
complex<double> I_inner_V(double** V_y, Data &data) {
    double p[3];
    double q[3];
    double r[3];
    for (int i = 0; i < 3; i++) {
        p[i] = V_y[1][i] - V_y[0][i]; // p = By-Ay
        q[i] = V_y[2][i] - V_y[0][i]; // q = Cy-Ay
        r[i] = V_y[3][i] - V_y[0][i]; // r = Dy-Ay
    }
    double p_len = sqrt(Scal_Prod(p, p));
    double q_len = sqrt(Scal_Prod(q, q));
    double r_len = sqrt(Scal_Prod(r, r));
    double h = max(max(p_len/data.n, q_len/data.n), r_len/data.n); // max(p_len/data.n, q_len/data.n, r_len/data.n);
    
    complex<double> Int(0, 0);
    double y[3];
    for (int i = 0; i < data.n; i++) {
        for (int j = 0; j < data.n; j++) {
            for (int k = 0; k < data.n; k++) {
                y[0] = V_y[0][0] + (i + 0.5) * p[0]/data.n + (j + 0.5) * q[0]/data.n + (k + 0.5) * r[0]/data.n;
                y[1] = V_y[0][1] + (i + 0.5) * p[1]/data.n + (j + 0.5) * q[1]/data.n + (k + 0.5) * r[1]/data.n;
                y[2] = V_y[0][2] + (i + 0.5) * p[2]/data.n + (j + 0.5) * q[2]/data.n + (k + 0.5) * r[2]/data.n;
                Int += Func(data, h, y);
            }
        }
    }
    double mu = p_len * q_len * r_len;
    Int *= mu/pow(data.n, 3);
    return Int; // Int * mu(V) / n^3
}

// V = ABCD - параллелепипед
// V_x[0] = Ax, V_x[1] = Bx, V_x[2] = Cx, V_x[3] = Dx
// V_y[0] = Ay, V_y[1] = By, V_y[2] = Cy, V_y[3] = Dy
complex<double> Integral_V(double** V_x, double** V_y, Data &data) {
    double p[3];
    double q[3];
    double r[3];
    for (int i = 0; i < 3; i++) {
        p[i] = V_x[1][i] - V_x[0][i]; // p = Bx-Ax
        q[i] = V_x[2][i] - V_x[0][i]; // q = Cx-Ax
        r[i] = V_x[3][i] - V_x[0][i]; // r = Dx-Ax
    }
    double p_len = sqrt(Scal_Prod(p, p));
    double q_len = sqrt(Scal_Prod(q, q));
    double r_len = sqrt(Scal_Prod(r, r));
    
    complex<double> Int(0, 0);
    double y[3];
    for (int i = 0; i < data.n; i++) {
        for (int j = 0; j < data.n; j++) {
            for (int k = 0; k < data.n; k++) {
                data.x[0] = V_x[0][0] + (i + 0.5) * p[0]/data.n + (j + 0.5) * q[0]/data.n + (k + 0.5) * r[0]/data.n;
                data.x[1] = V_x[0][1] + (i + 0.5) * p[1]/data.n + (j + 0.5) * q[1]/data.n + (k + 0.5) * r[1]/data.n;
                data.x[2] = V_x[0][2] + (i + 0.5) * p[2]/data.n + (j + 0.5) * q[2]/data.n + (k + 0.5) * r[2]/data.n;
                Int += I_inner_V(V_y, data);
                //cout << "I_inner_V" << I_inner_V(V_y, data) << endl;
            }
        }
    }
    double mu = p_len * q_len * r_len;
    Int *= mu/pow(data.n, 3);
    return Int; // Int * mu(V) / n^3
}

//---------------------------------------- Набор блоков матрицы размером 3х3-----------------------------------------------------//

// V = ABCD - параллелепипед
// V_x[0] = Ax, V_x[1] = Bx, V_x[2] = Cx, V_x[3] = Dx
double** Sigma_minus(int ind, double** V_x) {
    double **Sigma_x_minus = new double* [3];
    for (int i = 0; i < 3; i++){
        Sigma_x_minus[i] = new double [3]; 
    }

    if (ind == 0) {
    // Sigma_x_minus = Ax Cx Dx
        for (int j = 0; j < 3; j++) {
            Sigma_x_minus[0][j] = V_x[0][j]; // Ax
            Sigma_x_minus[1][j] = V_x[2][j]; // Cx
            Sigma_x_minus[2][j] = V_x[3][j]; // Dx
        }
    }
    else if (ind == 1) {
    // Sigma_x_minus = Ax Bx Dx
        for (int j = 0; j < 3; j++) {
            Sigma_x_minus[0][j] = V_x[0][j]; // Ax
            Sigma_x_minus[1][j] = V_x[1][j]; // Bx
            Sigma_x_minus[2][j] = V_x[3][j]; // Dx
        }
    }
    else {
    // Sigma_x_minus = Ax Bx Cx
        for (int j = 0; j < 3; j++) {
            Sigma_x_minus[0][j] = V_x[0][j]; // Ax
            Sigma_x_minus[1][j] = V_x[1][j]; // Bx
            Sigma_x_minus[2][j] = V_x[2][j]; // Cx
        }
    }
    return Sigma_x_minus;
}

// V = ABCD - параллелепипед
// V_x[0] = Ax, V_x[1] = Bx, V_x[2] = Cx, V_x[3] = Dx
double** Sigma_plus(int ind, double** V_x) {
    double **Sigma_x_plus = new double* [3];
    for (int i = 0; i < 3; i++){
        Sigma_x_plus[i] = new double [3]; 
    }

    if (ind == 0) {
    // Sigma_x_plus = Bx Cx+(Bx-Ax) Dx+(Bx-Ax)
        for (int j = 0; j < 3; j++) {
            Sigma_x_plus[0][j] = V_x[1][j]; // Bx
            Sigma_x_plus[1][j] = V_x[2][j] + (V_x[1][j] - V_x[0][j]); // Cx+(Bx-Ax)
            Sigma_x_plus[2][j] = V_x[3][j] + (V_x[1][j] - V_x[0][j]); // Dx+(Bx-Ax)
        }
    }
    else if (ind == 1) {
    // Sigma_x_plus = Cx Bx+(Cx-Ax) Dx+(Cx-Ax)
        for (int j = 0; j < 3; j++) {
            Sigma_x_plus[0][j] = V_x[2][j]; // Cx
            Sigma_x_plus[1][j] = V_x[1][j] + (V_x[2][j] - V_x[0][j]); // Bx+(Cx-Ax)
            Sigma_x_plus[2][j] = V_x[3][j] + (V_x[2][j] - V_x[0][j]); // Dx+(Cx-Ax)
        }
    }
    else {
    // Sigma_x_plus = Dx Bx+(Dx-Ax) Cx+(Dx-Ax)
        for (int j = 0; j < 3; j++) {
            Sigma_x_plus[0][j] = V_x[3][j]; // Dx
            Sigma_x_plus[1][j] = V_x[1][j] + (V_x[3][j] - V_x[0][j]); // Bx+(Dx-Ax)
            Sigma_x_plus[2][j] = V_x[2][j] + (V_x[3][j] - V_x[0][j]); // Cx+(Dx-Ax)
        }
    }
    return Sigma_x_plus;
}

complex<double> I_1(int l, int k, double** V_x, double** V_y, Data &data) {
    complex<double> Int(0, 0);
    complex<double> Int_1 = Integral_Sigma(Sigma_plus(l, V_x),  Sigma_plus(k, V_y), data);
    complex<double> Int_2 = Integral_Sigma(Sigma_plus(l, V_x),  Sigma_minus(k, V_y), data);
    complex<double> Int_3 = Integral_Sigma(Sigma_minus(l, V_x),  Sigma_plus(k, V_y), data);
    complex<double> Int_4 = Integral_Sigma(Sigma_minus(l, V_x),  Sigma_minus(k, V_y), data);
    Int = Int_2 + Int_3 - Int_1 - Int_4;
    return Int;
}

complex<double> I_2(double** V_x, double** V_y, Data &data) {
    return K*K*Integral_V(V_x, V_y, data);
}

// Задаёт блок A_ij размера 3х3, находящийся на позициях I, J матрицы Matrix размера (3Nx3N)
void Count_Matrix_Block(complex<double>** A_ij, int I, int J, double** V_x, double** V_y, Data &data, double E_cell) {
    complex<double>  Int_2 = I_2(V_x, V_y, data);
    //cout << "I_2 = " << Int_2 << endl;
    //double E = 4; // E(x)
    for (int l = 0; l < 3; l++) {
        for (int k = 0; k < 3; k++) {
            A_ij[l][k] = -(E_cell-1);  // E_cell???
            if (l == k){ // I_1 + I_2
                A_ij[l][k] *= (I_1(l, k, V_x, V_y, data) + Int_2);
                if (I == J) {
                    A_ij[l][k] += data.V_cell;
                    //cout << data.V_cell << endl;
                }
                //cout << "I_1 = " << I_1(l, k, V_x, V_y, data) << endl;
            }
            else { // I_1
                A_ij[l][k] *= I_1(l, k, V_x, V_y, data);
                //cout << "I_1 = " << I_1(l, k, V_x, V_y, data) << endl;
            }
        }
    }
}

//---------------------------------------- Набор матрицы -----------------------------------------------------//

// задаём параллелепипед V_x точками Ax Bx Cx Dx
void Get_V(double** V, double* ind, Data &data) {
    double h[3]; // шаг сетки
    // длины сторон L[0] = AB, L[1] = AC, L[2] = AD параллелепипеда V = AB AC AD
    h[0] = data.L[0]/data.N[0]; h[1] = data.L[1]/data.N[1]; h[2] = data.L[2]/data.N[2];
    
    for (int j = 0; j < 3; j++) {
        V[0][j] = data.A[j] + ind[j]*h[j]; // V[0] = Ax = A + i_1*a + i_2*c + i_3*d
        V[1][j] = V[0][j]; V[2][j] = V[0][j]; V[3][j] = V[0][j];
    }
    V[1][0] += h[0]; // V[1] = Bx = Ax + (B-A)/N[0]
    V[2][1] += h[1]; // V[2] = Cx = Ax + (C-A)/N[1]
    V[3][2] += h[2]; // V[3] = Dx = Ax + (D-A)/N[2]
}

// A[i] = V[0][i], B[i] = V[1][i], C[i] = V[2][i], D[i] = V[3][i]
void Get_Matrix(complex<double>** Matrix, Data &data) {
    cout << "Check_2_Get_Matrix:" << endl;
    complex<double> **A_ij = new complex<double>* [3];    
    for (int i = 0; i < 3; i++){
        A_ij[i] = new complex<double> [3];   
        for (int j = 0; j < 3; j++) {
            A_ij[i][j] = 0;
        }
    }

    double **V_x = new double* [4];
    double **V_y = new double* [4];
    for (int i = 0; i < 4; i++) {
        V_x[i] = new double [3]; 
        V_y[i] = new double [3]; 
    }
    // Набор матрицы
    int I = 0;
    //cout << "I = " << I << endl;
    int J = 0;
    //cout << "J = " << J << endl;
    double *ind_x = new double [3];
    double *ind_y = new double [3];
    int N_123 = data.N[0]*data.N[1]*data.N[2];
    cout << "N_123 = " << N_123 << endl;
    double *E_cell = new double [N_123];
    Read_Eps_on_cell(E_cell, data); //считаем E на ячейке
    //for (int i = 0; i < N_123; i++) {
    //    E_cell[i] = 4;
    //}
    
    // Задаем V_x - I
    for (int i_1 = 0; i_1 < data.N[0]; i_1++) {
        for (int i_2 = 0; i_2 < data.N[1]; i_2++) {
            for (int i_3 = 0; i_3 < data.N[2]; i_3++) {
                J = 0;
                ind_x[0] = i_1; ind_x[1] = i_2; ind_x[2] = i_3; 

                // Задаем V_y - J
                for (int j_1 = 0; j_1 < data.N[0]; j_1++) {
                    for (int j_2 = 0; j_2 < data.N[1]; j_2++) {
                        for (int j_3 = 0; j_3 < data.N[2]; j_3++) {
                            ind_y[0] = j_1; ind_y[1] = j_2; ind_y[2] = j_3;
                            if (I == J) {
                                data.n = 10;
                            }
                            else {
                                data.n = 2;
                            }
                            Get_V(V_x, ind_x, data); // задаём параллелепипед V_x точками Ax Bx Cx Dx
                            Get_V(V_y, ind_y, data); // задаём параллелепипед V_y точками Ay By Cy Dy
                            Count_Matrix_Block(A_ij, I, J, V_x, V_y, data, E_cell[J]); // подсчёт блока 3х3
                            // Заполнение матрицы
                            /*
                            cout << "A_ij:" << endl;
                            for (int q = 0; q < 3; q++) {
                                for (int p = 0; p < 3; p++) {
                                    cout << A_ij[q][p] << ' ';
                                }
                                cout << endl;
                            }
                            cout << endl;

                            cout << "V_x:" << endl;
                            for (int q = 0; q < 4; q++) {
                                for (int p = 0; p < 3; p++) {
                                    cout << V_x[q][p] << ' ';
                                }
                                cout << endl;
                            }
                            cout << endl;

                            cout << "V_y:" << endl;
                            for (int q = 0; q < 4; q++) {
                                for (int p = 0; p < 3; p++) {
                                    cout << V_y[q][p] << ' ';
                                }
                                cout << endl;
                            }
                            cout << endl;
                            return;
                            */

                            for (int q = 0; q < 3; q++) {
                                for (int p = 0; p < 3; p++) {
                                    Matrix[I*3 + q][J*3 + p] = A_ij[q][p];
                                }
                            }
                            J++;
                            //cout << "J = " << J << endl;
                        }
                    }
                }
                if (J != N_123) {
                    cout << "J != N1*N2*N3, J = " << J << endl;
                    return;
                }
                I++;
            }
        }
    }
    if (I != N_123) {
        cout << "I != N1*N2*N3, I = " << I << endl;
        return;
    }
    cout << "Got_Matrix!" << endl;
}

//---------------------------------------- Правая часть -----------------------------------------------------//

complex<double> *Read_Znach_Func(Data &data){
    int N123 = data.N[0] * data.N[1] * data.N[2];
    complex<double> *E_inc = new complex<double> [3*N123];

    double *E_inc_re = new double [3*N123];
    double *E_inc_imag = new double [3*N123];

    cout << endl;
    cout << "Check_4_Read_Znach_Func:" << endl;

    //Right_Part_re.txt
    ifstream in_re("/Users/deti/Desktop/ВТМ/Курсач/Dielectric_scattering/Right_Part_re.txt"); // окрываем файл для чтения
    double ch;
    int cnt = 0, ind = 0;
    while (in_re >> ch) {
        if (ind < 3*N123){
            E_inc_re[ind] = ch;
        }
        cnt++;
        ind++;
    }
    cout << "Чисел в файле Right_Part_re.txt: " << cnt << endl;
    in_re.close();
    
    if (cnt == 3*N123){
        cout << "File Right_Part_re.txt has been read correctly" << endl;
    }
    else {
        cout << "File Right_Part_re.txt has been read with errors" << endl;
        cout << "cnt_re = " << cnt << endl;
    }
    //cout << endl;

    //Right_Part_imag.txt
    ifstream in_imag("/Users/deti/Desktop/ВТМ/Курсач/Dielectric_scattering/Right_Part_imag.txt"); // окрываем файл для чтения
    cnt = 0, ind = 0;
    while (in_imag >> ch) {
        if (ind < 3*N123){
            E_inc_imag[ind] = ch;
        }
        cnt++;
        ind++;
    }
    cout << "Чисел в файле Right_Part_imag.txt: " << cnt << endl;
    in_imag.close();
    
    if (cnt == 3*N123){
        cout << "File Right_Part_imag.txt has been read correctly" << endl;
    }
    else {
        cout << "File Right_Part_imag.txt has been read with errors" << endl;
        cout << "cnt_imag = " << cnt << endl;
    }
    //cout << endl;

    for (int i = 0; i < 3*N123; i++) {
        complex<double> comp(E_inc_re[i], E_inc_imag[i]);
        E_inc[i] = comp;
        //cout << comp << ' ';
    }
    cout << endl;

    return E_inc;
}

//---------------------------------------- Решение СЛАУ -----------------------------------------------------//

void Write_g(complex<double>** g, int n){
    ofstream out;          // поток для записи
    out.open("Gal_Re_g.gv");      // открываем файл для записи
    if (out.is_open()){
        out << "2\t" << n/2 << endl;
        for (int i = 0; i < n; i++){
            out << g[i][0].real() << " \t " << g[i][1].real() << " \t " << g[i][2].real() << endl;
        }
    }
    out.close();
    cout << "File Gal_Re_g.gv has been written" << endl;

    out.open("Gal_Im_g.gv");      // открываем файл для записи
    if (out.is_open()){
        out << "2\t" << n/2 << endl;
        for (int i = 0; i < n; i++){
            out << g[i][0].imag() << " \t " << g[i][1].imag() << " \t " << g[i][2].imag() << endl;
        }
    }
    out.close();
    cout << "File Gal_Im_g.gv has been written" << endl;
}

//---------------------------------------- Подсчет ЭПР -----------------------------------------------------//

///////   Перевод вектора в матрицу  ///////
complex<double>** Vec2Matr(complex<double> *x, int n){
    complex<double> **A = new complex<double> *[n];
    for (int i = 0; i < n; i++){
        A[i] = new complex<double> [3];
        for (int p = 0; p < 3; p++){
            A[i][p] = x[3*i + p];
        }
    }
    return A;
}

void Read_Colloc_Points(double** colloc_points, int n){
    cout << endl;
    cout << "Check_6_Read_Colloc_Points:" << endl;

    ifstream in("/Users/deti/Desktop/ВТМ/Курсач/Dielectric_scattering/Colloc_Points.txt"); // окрываем файл для чтения
    int size;
    double V_cell;
    in >> size;
    cout << "Чисел в файле Colloc_Points.txt: " << size << endl;
    if (size == n){
        cout << "File Colloc_Points.txt has been read correctly" << endl;
    }
    else {
        cout << "File Colloc_Points.txt has been read with errors" << endl;
    }
    in >> V_cell;
    cout << "V_cell = " << V_cell << endl;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < 3; j++) {
            in >> colloc_points[i][j];
        }
    }
    in.close();
}

///////   Подсчет ЭПР   ///////
void Effective_Scattering_Area(double* sigma, complex<double>** g, Data &data){
    double tau[181][3];
    int N123 = data.N[0] * data.N[1] * data.N[2];
    for (int alpha = 0; alpha < 181; alpha++){
        tau [alpha][0] = cos(alpha*M_PI/180);
        tau [alpha][1] = sin(alpha*M_PI/180);
        tau [alpha][2] = 0;
    }

    //complex<double> sigma[181];
    complex<double> h[181][3];
    double E_inc = 1;
    for (int i = 0; i < 181; i++){
        for (int j = 0; j < 3; j++){
            h[i][j] = 0;
        }
    }

    double **colloc_points = new double* [N123];    
    for (int i = 0; i < N123; i++){
        colloc_points[i] = new double [3];  
    }
    Read_Colloc_Points(colloc_points, N123);
    //for (int i = 0; i < N123; i++) {
    //    cout << colloc_points[i][0] << ' ' << colloc_points[i][1] << ' ' << colloc_points[i][2] << endl;
    //}

    double *E_cell = new double [N123];
    Read_Eps_on_cell(E_cell, data); //считаем E на ячейке

    complex<double> exp, coeff;
    for (int alpha = 0 ; alpha < 181; alpha++){
        cout << "alpha = " << alpha << endl;
        for (int i = 0; i < N123; i++){
            cout << '1' << endl;
            cout << "ind = " << i << ' ';
            exp = Comp_Exp(-K * Scal_Prod(tau[alpha], colloc_points[i]));
            coeff = exp * (E_cell[i] - 1) * K * K * data.V_cell;
            for (int p = 0; p < 3; p++){
                h[alpha][p] += coeff * (g[i][p] - Scal_Prod(g[i], tau[alpha]) * tau[alpha][p]);
            }
        }
        cout << endl;
        //sigma[alpha] = (4 * M_PI) * Comp_abs(h[alpha]) * Comp_abs(h[alpha]);
        sigma[alpha] = Comp_abs(h[alpha]) * Comp_abs(h[alpha]) / (E_inc * 4 * M_PI);
    }
}

void Write_ESA(double *e, int n){
    ofstream out;          // поток для записи
    out.open("Gal_ESA.txt");      // открываем файл для записи
    if (out.is_open()){
        for (int i = 0; i < n; i++){
            out << e[i] << ", ";
        }
    }
    out.close();
    cout << "File Gal_ESA.txt has been written" << endl;
}