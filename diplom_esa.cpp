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


void Read_g(complex<double> **g,int n){
    double **g_re = new double* [n];    
    for (int i = 0; i < n; i++){
        g_re[i] = new double [3];  
    }
    double **g_imag = new double* [n];    
    for (int i = 0; i < n; i++){
        g_imag[i] = new double [3];  
    }
    cout << "n = N^3 = " << n << endl;

    //Gal_Re_g.gv
    ifstream in_re("/Users/deti/Desktop/ВТМ/Diplom/Gal_Re_g.gv"); // окрываем файл для чтения
    double ch;
    in_re >> ch;
    cout << "points in Gal_Re_g.gv: " << ch << " ";
    in_re >> ch;
    cout << ch << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 3; j++) {
            in_re >> g_re[i][j];
        }
    }
    in_re.close();
    
    //Gal_Im_g.gv
    ifstream in_imag("/Users/deti/Desktop/ВТМ/Diplom/Gal_Im_g.gv"); // окрываем файл для чтения
    in_imag >> ch;
    cout << "points in Gal_Im_g.gv: " << ch << " ";
    in_imag >> ch;
    cout << ch << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 3; j++) {
            in_imag >> g_imag[i][j];
        }
    }
    in_imag.close();

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 3; j++) {
            complex<double> comp(g_re[i][j], g_imag[i][j]);
            g[i][j] = comp;
        }
    }
    cout << "Got g" << endl;
    cout << endl;
}

double Read_Colloc_Points(double** colloc_points, int n){
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
    cout << endl;
    return V_cell;
}

void Read_Eps_on_cell(double *E_cell, int n) {
    ifstream in("/Users/deti/Desktop/ВТМ/Курсач/Dielectric_scattering/Epsilon.txt"); // окрываем файл для чтения
    double ch;
    int cnt = 0, ind = 0;
    while (in >> ch) {
        if (ind < n){
            E_cell[ind] = ch;
        }
        //cout << E_cell[ind] << endl;
        cnt++;
        ind++;
    }
    cout << "Чисел в файле Epsilon.txt: " << cnt << endl;
    in.close();
    
    if (cnt == (n)){
        cout << "File Epsilon.txt has been read correctly" << endl;
    }
    else {
        cout << "File Epsilon.txt has been read with errors" << endl;
        return;
    }
    cout << endl;
}

///////   Подсчет ЭПР   ///////
void Effective_Scattering_Area(double* sigma, complex<double>** g, int n){
    double tau[181][3];
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

    double **colloc_points = new double* [n];    
    for (int i = 0; i < n; i++){
        colloc_points[i] = new double [3];  
    }
    double V_cell = Read_Colloc_Points(colloc_points, n);
    //for (int i = 0; i < n; i++) {
    //    cout << colloc_points[i][0] << ' ' << colloc_points[i][1] << ' ' << colloc_points[i][2] << endl;
    //}

    double *E_cell = new double [n];
    Read_Eps_on_cell(E_cell, n); //считаем E на ячейке

    complex<double> exp, coeff;
    for (int alpha = 0 ; alpha < 181; alpha++){
        //cout << "alpha = " << alpha << endl;
        for (int i = 0; i < n; i++){
            exp = Comp_Exp(-K * Scal_Prod(tau[alpha], colloc_points[i]));
            coeff = exp * (E_cell[i] - 1) * K * K * V_cell;
            for (int p = 0; p < 3; p++){
                h[alpha][p] += coeff * (g[i][p] - Scal_Prod(g[i], tau[alpha]) * tau[alpha][p]);
            }
        }
        //sigma[alpha] = (4 * M_PI) * Comp_abs(h[alpha]) * Comp_abs(h[alpha]);
        sigma[alpha] = Comp_abs(h[alpha]) * Comp_abs(h[alpha]) / (E_inc * 4 * M_PI);
    }
}

void Write_Gal_ESA(double *e, int n){
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

int main() {
    int N = 10;
    int n = N*N*N;
    complex<double> **g = new complex<double>* [n];    
    for (int i = 0; i < n; i++){
        g[i] = new complex<double> [3];  
    }
    Read_g(g, n);

    double sigma[181];
    Effective_Scattering_Area(sigma, g, n);

    ///////   Перевод ЭПР в Децибелы (dB) ///////
    //cout << "Check_7_Effective_Scattering_Area" << endl;
    cout << "ЭПР в dB:" << endl;
    double sigma_dB[181];
    for (int i = 0; i < 181; i++){
        sigma_dB[i] = log10(sigma[i] / (M_PI)) * 10.; // ??????/ как будто надо разделить на 128 и ответ сойдется
        //cout << sigma_dB[i] << ", ";
    }
    //cout << endl << endl;
    Write_Gal_ESA(sigma_dB, 181);
}