#ifndef DIPLOM_DIPLOM_FUNC_H
#define DIPLOM_DIPLOM_FUNC_H

#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>
#include <string> 
using namespace std;

///////   Константы   //////
const int C = 299792458 ;    // скорость света
const double Omega = 0.15 * (1e+9) * 2 * M_PI; // 3
const double K = Omega/C; // волновое число

//---------------------------------------- Вспомогательные функции -----------------------------------------------------//

///////   Комплексная экспонента   ///////
complex<double> Comp_Exp(double fi);


///////   Модуль комплексного вектора   ///////
double Comp_abs(complex<double>* h, int n);


///////   Скалярное произведение 1   ///////
complex<double> Scal_Prod(complex<double>* x, complex<double>* y, int n);


///////   Скалярное произведение 2   ///////
complex<double> Scal_Prod(complex<double>* x, double* y, int n);


///////   Скалярное произведение 3   ///////
double Scal_Prod(double* x, double* y, int n);

//---------------------------------------- Структура -----------------------------------------------------//

struct Data {
    int n;
    double eps;
    int N[3]; 
    double L[3];
    double A[3];
    double x[3];
    double V_cell;

    double O_sphere[3];
    double R;

    Data(int* N, double* L, double* A, int n);
};

//---------------------------------------- Заданная функция -----------------------------------------------------//

int Chi(double *x, Data &data);

///////   Функция Эпсилон   ///////
double Epsilon(double *x);

//double Eps_on_cell(double** V_x, Data &data);

void Read_Eps_on_cell(double *E_cell, Data &data);

double Theta(double x);

complex<double> Func(Data data, int h, double* y);

//---------------------------------------- Поверхностный интеграл -----------------------------------------------------//

// Sigma_y[0] = Ay, Sigma_y[1] = By, Sigma_y[2] = Dy
complex<double> I_inner_Sigma(double** Sigma_y, Data &data);

// typedef complex<double> (*func)(Data, int, double*);
// S = ABD - прямоугольник 
// Sigma_x[0] = Ax, Sigma_x[1] = Bx, Sigma_x[2] = Dx
// Sigma_y[0] = Ay, Sigma_y[1] = By, Sigma_y[2] = Dy
complex<double> Integral_Sigma(double** Sigma_x, double** Sigma_y, Data &data);

//---------------------------------------- Объемный интеграл -----------------------------------------------------//

// V_y[0] = Ay, V_y[1] = By, V_y[2] = Cy, V_y[3] = Dy
complex<double> I_inner_V(double** V_y, Data &data);

// V = ABCD - параллелепипед
// V_x[0] = Ax, V_x[1] = Bx, V_x[2] = Cx, V_x[3] = Dx
// V_y[0] = Ay, V_y[1] = By, V_y[2] = Cy, V_y[3] = Dy
complex<double> Integral_V(double** V_x, double** V_y, Data &data);

//---------------------------------------- Набор блоков матрицы размером 3х3-----------------------------------------------------//

// V = ABCD - параллелепипед
// V_x[0] = Ax, V_x[1] = Bx, V_x[2] = Cx, V_x[3] = Dx
double** Sigma_minus(int ind, double** V_x);

// V = ABCD - параллелепипед
// V_x[0] = Ax, V_x[1] = Bx, V_x[2] = Cx, V_x[3] = Dx
double** Sigma_plus(int ind, double** V_x);

complex<double> I_1(int l, int k, double** V_x, double** V_y, Data &data); 

complex<double> I_2(double** V_x, double** V_y, Data &data); 

// Задаёт блок A_ij размера 3х3, находящийся на позициях I, J матрицы Matrix размера (3Nx3N)
void Count_Matrix_Block(complex<double>** A_ij, int I, int J, double** V_x, double** V_y, Data &data, double E_cell); 

//---------------------------------------- Набор матрицы -----------------------------------------------------//

// задаём параллелепипед V_x точками Ax Bx Cx Dx
void Get_V(double** V, double* ind, Data &data);

// A[i] = V[0][i], B[i] = V[1][i], C[i] = V[2][i], D[i] = V[3][i]
void Get_Matrix(complex<double>** Matrix, Data &data);

//---------------------------------------- Правая часть -----------------------------------------------------//

complex<double>* Read_Znach_Func(Data &data);

//---------------------------------------- Решение СЛАУ -----------------------------------------------------//

extern "C" void zgesv(int *n, int *nrhs, complex<double>  *A, int *lda, int* ipiv, complex<double> *b, int *ldb, int *info);

void Write_g(complex<double>** g, int n);

//---------------------------------------- Подсчет ЭПР -----------------------------------------------------//

void Read_Colloc_Points(double **colloc_points, int n);

complex<double>** Vec2Matr(complex<double> *x, int n);

void Effective_Scattering_Area(double* sigma, complex<double>** g, Data &data);

void Write_ESA(double *e, int n);

#endif //DIPLOM_DIPLOM_FUNC_H