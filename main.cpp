#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cassert>

#include "typedefs.h"
#include "Givens.h"
#include "Householder.h"
#include "FactorizadorQR.h"
#include "CML.h"

using namespace std;

void imprimirOctave(const matrizDouble& A) {
    unsigned int m = A.size(), n = A[0].size();
    cerr << "[ ";
    for (unsigned int i = 0; i < m; i++) {
        for (unsigned int j = 0; j < n; j++) {
            cerr << A[i][j] << ", ";
        }
        cerr << "; ";
    }
    cerr << " ];" << endl << endl;

}

void imprimirOctave(const vectorDouble& v) {
    unsigned int m = v.size();
    cerr << "[ ";
    for (unsigned int i = 0; i < m; i++) {
        cerr << v[i] << ", ";
    }
    cerr << " ]';" << endl << endl;
}

void imprimirMatriz(const matrizDouble& A) {
    unsigned int m = A.size(), n = A[0].size();
    cerr.precision(5);
    cerr << endl;
    for (unsigned int i = 0; i < m; i++) {
        for (unsigned int j = 0; j < n; j++) {
            cerr << A[i][j] << "\t";
        }
        cerr << endl;
    }
    cerr << endl;

}

void prueba1() {
    unsigned int m = 10, n = 3;
    matrizDouble A(m, vectorDouble(n, 0));
    srand(0);
    for (unsigned int i = 0; i < m; i++) {
        for (unsigned int j = 0; j < n; j++) {
            A[i][j] = round(rand() % 10000) / 10000;
        }
    }
    Givens g1(A, 0, 1);
    imprimirOctave(A);
    g1.aplicarM(A);
    imprimirOctave(A);
    exit(0);
}

void prueba2() {
    unsigned int m = 5, n = 3;
    matrizDouble A(m, vectorDouble(n, 0));
    srand(0);
    for (unsigned int i = 0; i < m; i++) {
        for (unsigned int j = 0; j < n; j++) {
            A[i][j] = round(rand() % 10000) / 10000;
        }
    }
    Givens g01(A, 0, 1);
    imprimirOctave(A);
    g01.aplicarM(A);
    Givens g02(A, 0, 2);
    g02.aplicarM(A);
    Givens g03(A, 0, 3);
    g03.aplicarM(A);
    Givens g04(A, 0, 4);
    g04.aplicarM(A);
    Givens g12(A, 1, 2);
    g12.aplicarM(A);
    Givens g13(A, 1, 3);
    g13.aplicarM(A);
    Givens g14(A, 1, 4);
    g14.aplicarM(A);
    imprimirOctave(A);
    Householder h(A, 2);
    h.aplicarM(A);
    imprimirOctave(A);
    exit(0);
}

void prueba3() {
    unsigned int m = 15, n = 7;
    matrizDouble A(m, vectorDouble(n, 0));
    srand(0);
    for (unsigned int i = 0; i < m; i++) {
        for (unsigned int j = 0; j < n; j++) {
            A[i][j] = round(rand() % 10000) / 10000;
        }
    }
    //imprimirMatriz(A);
    for (unsigned int i = 0; i < n - 1; i++) {
        for (unsigned int j = i + 1; j < m; j++) {
            Givens g(A, i, j);
            g.aplicarM(A);
        }
    }
    //imprimirMatriz(A);
    //imprimirOctave(A);
    Householder h(A, n - 1);
    h.aplicarM(A);
    imprimirOctave(A);
    cerr << endl;
    exit(0);
}

void prueba4() {
    unsigned int m = 15, n = 7;
    matrizDouble A(m, vectorDouble(n, 0));
    matrizDouble B(m, vectorDouble(n, 0));
    srand(0);
    double val;
    for (unsigned int i = 0; i < m; i++) {
        for (unsigned int j = 0; j < n; j++) {
            val= round(rand() % 10000) / 10000;
            A[i][j] = val;
            B[i][j] = val;
        }
    }
    FactorizadorQR qr(A);
    imprimirOctave(A);
    imprimirOctave(B);
    qr.aplicarM(B);
    imprimirOctave(B);
    cerr << endl;
    exit(0);
}

void prueba5(){
    unsigned int m = 15, n = 7;
    matrizDouble A(m, vectorDouble(n, 0));
    vectorDouble b(m,0);
    vectorDouble res(n,0);
    srand(0);
    double val;
    for (unsigned int i = 0; i < m; i++) {
        for (unsigned int j = 0; j < n; j++) {
            val= round(rand() % 10000) / 10000;
            A[i][j] = val;
            
        }
        val= round(rand() % 10000) / 10000;
        b[i]=val;
    }
    cout << "A = ";imprimirOctave(A);
    cout << "b = ";imprimirOctave(b);
    CML cml(A);
    cout << "R= ";imprimirOctave(A);
    cml.resolver(b,res);
    cout << "Qtb = ";imprimirOctave(b);
    cout << "x = ";imprimirOctave(res);
}

void prueba7(){
    unsigned int m = 15, n = 7;
    matrizDouble A(m, vectorDouble(n, 0));
    vectorDouble b(m,0);
    vectorDouble res(n,0);
    srand(0);
    double val;
    for (unsigned int i = 0; i < m; i++) {
        for (unsigned int j = i; j < n; j++) {
            val= round(rand() % 10000) / 10000;
            A[i][j] = val;
            
        }
        val= round(rand() % 10000) / 10000;
        b[i]=val;
    }
    cout << "A=";imprimirOctave(A);
    //FactorizadorQR q(A);
    //cout << "R=";imprimirOctave(A);
    cout << "b=";imprimirOctave(b);
    //q.aplicarV(b);
    //cout << "Qtb=";imprimirOctave(b);
    unsigned int i;
    double acum, a_ii, a_ij;
    assert(n==res.size());
    for (unsigned int _i = 0; _i < n; _i++) {
        i = (n - 1) - _i;
        a_ii = A[i][i];
        acum = 0;
        for (unsigned int j = i + 1; j < n; j++) {
            a_ij =  A[i][j];
            acum += res[j] * a_ij/a_ii;
        }
        res[i] = b[i]/a_ii - acum;
    }
    cout << "x=";imprimirOctave(res);
}

void prueba8(){
    unsigned int m = 15, n = 7;
    matrizDouble A(m, vectorDouble(n, 0));
    vectorDouble b(m,0);
    vectorDouble res(n,0);
    unsigned int i;
    double acum, a_ii, a_ij;
    srand(0);
    double val;
    for (unsigned int i = 0; i < m; i++) {
        for (unsigned int j = 0; j < n; j++) {
            val= round(rand() % 10000) / 10000;
            A[i][j] = val;
            
        }
        val= round(rand() % 10000) / 10000;
        b[i]=val;
    }
    cout << "A=";imprimirOctave(A);
    FactorizadorQR q(A);
    cout << "R=";imprimirOctave(A);
    cout << "b1=";imprimirOctave(b);
    for (unsigned int _i = 0; _i < n; _i++) {
        i = (n - 1) - _i;
        a_ii = A[i][i];
        acum = 0;
        for (unsigned int j = i + 1; j < n; j++) {
            a_ij =  A[i][j];
            acum += res[j] * a_ij/a_ii;
        }
        res[i] = b[i]/a_ii - acum;
    }
    cout << "x1=";imprimirOctave(res);
    q.aplicarV(b);
    cout << "b2=";imprimirOctave(b);
    for (unsigned int _i = 0; _i < n; _i++) {
        i = (n - 1) - _i;
        a_ii = A[i][i];
        acum = 0;
        for (unsigned int j = i + 1; j < n; j++) {
            a_ij =  A[i][j];
            acum += res[j] * a_ij/a_ii;
        }
        res[i] = b[i]/a_ii - acum;
    }
    cout << "x2=";imprimirOctave(res);
}


int main(int argc, char** argv) {
    //prueba5();
    
    //prueba7();
    prueba8();
    return 0;
}

