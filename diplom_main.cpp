#include "diplom_func.h"

//---------------------------------------- main -----------------------------------------------------//

int main() {
    double O_sphere[3] = {0, 0 ,0}; // центр сферы
    double R = 1; //радиус сферы

    int n = 2; // число разбиений по стороне прямоугольной ячейки
    int N[3] = {10, 10, 10}; // разбиение сторон AB, AC, AD соответственно 

    double L[3]; // длины сторон L1 = AB, L2 = AC, L3 = AD параллелепипеда ABCD
    double A[3]; // начальная точка
    for (int i = 0; i < 3; i++) {
        A[i] = O_sphere[i] - R;
        L[i] = 2*R;
    }

    Data data(N, L, A, n);
    cout << "Check_1_data:" << endl;
    data.R = R;
    for(int i = 0; i < 3; i++) {
       data.O_sphere[i] = O_sphere[i]; 
    }
    cout << "N = " << data.N[0] << ' ' << data.N[0] << ' ' << data.N[0] << endl;
    cout << "L = " << data.L[0] << ' ' << data.L[0] << ' ' << data.L[0] << endl;
    cout << "A = " << data.A[0] << ' ' << data.A[0] << ' ' << data.A[0] << endl;
    cout << "O_sphere = " << data.O_sphere[0] << ' ' << data.O_sphere[0] << ' ' << data.O_sphere[0] << endl;
    cout << "R = " << data.R << endl;
    cout << "n = " << data.n << endl;
    cout << "V_cell = " << data.V_cell << endl;
    cout << endl;

    /*
    int N123 = N[0] * N[1] * N[2];
    double **colloc_points = new double* [N123];    
    for (int i = 0; i < N123; i++){
        colloc_points[i] = new double [3];  
    }
    Read_Colloc_Points(colloc_points, N123);
    for (int i = 0; i < N123; i++) {
        cout << colloc_points[i][0] << ' ' << colloc_points[i][1] << ' ' << colloc_points[i][2] << endl;
    }
    return 0;
    */
    
    int N123 = N[0] * N[1] * N[2];
    complex<double> **Matrix = new complex<double>* [3*N123];    // Создаем массив указателей Matrix
    for (int i = 0; i < 3*N123; i++){
        Matrix[i] = new complex<double> [3*N123];    // Создаем элементы Matrix
        for (int j = 0; j < 3*N123; j++){
            Matrix[i][j] = 0;
        }
    }
    Get_Matrix(Matrix, data);
    
    /*
    cout << "Matrix (3n)x(3n):" << endl;
    for (int i = 0; i < 3*N123; i ++) {
        for (int j = 0; j < 3*N123; j++) {
            cout << Matrix[i][j] << ' ';
        }
        cout << endl; 
    }
    cout << endl;
    */

    complex<double> *Matrix_vect = new complex<double> [3*N123 * 3*N123];
    for (int i = 0; i < 3*N123; i++){
        for (int j = 0; j < 3*N123; j++){
            Matrix_vect[i*3*N123 + j] = Matrix[j][i];
        }
    }

    ///////   Решение СЛАУ   ///////
    complex<double> *g = new complex<double> [3*N123]; 
    for (int i = 0; i < 3*N123; i++){
        g[i] = 0;
    }
    g = Read_Znach_Func(data);

    int dim = 3*N123;
    int nrhs = 1;
    int LDA = dim;
    int *ipiv = new int[dim];
    int LDB = dim;
    int info; 
    zgesv(&dim, &nrhs, Matrix_vect, &LDA, ipiv, g, &LDB, &info); // ???

    //cout << "Check _5_Decision:" << endl;
    //cout << "g =" << endl;
    //for (int i = 0; i < dim; i++){
    //    cout << g[i] << ' ';
    //}
    cout << "Info = " << info << endl;
    cout << endl;
    
    cout << "Check_5_Write_g" << endl;
    Write_g(Vec2Matr(g, N123), N123);
    cout << endl;

    /*
    ///////   Подсчет ЭПР  ///////
    double sigma[181];
    Effective_Scattering_Area(sigma, Vec2Matr(g, n), data);
    
    cout << "ЭПР:" << endl;
    for (int i = 0; i < 181; i++){
        cout << sigma[i] << ' ';
    }
    cout << endl;


    ///////   Перевод ЭПР в Децибелы (dB) ///////
    cout << "Check_7_Effective_Scattering_Area" << endl;
    cout << "ЭПР в dB:" << endl;
    double sigma_dB[181];
    for (int i = 0; i < 181; i++){
        sigma_dB[i] = log10(sigma[i] / (M_PI)) * 10.; // ??????/ как будто надо разделить на 128 и ответ сойдется
        //cout << sigma_dB[i] << ", ";
    }
    //cout << endl << endl;
    Write_ESA(sigma_dB, 181);
    */
    
    ///////   Освобождаем память   ///////
    for (int i = 0; i < 3*N123; i++){
        delete[] Matrix[i]; 
    }
    delete [] Matrix;
    delete [] Matrix_vect;
    delete [] g;
    delete [] ipiv;

    return 0;
}
