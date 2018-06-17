#include <cassert>
#include <math.h>
#include <iostream>

#include "Householder.h"
#include "constantes.h"

Householder::Householder(vectorDouble& u) : u(u) {
}

double auxSubNorma2CuadradoColumna(matrizDouble& A, unsigned int _i) {
    double normaCuadrado = 0.0;
    for (unsigned int i = _i; i < A.size(); i++) {
        normaCuadrado += pow(A[i][_i], 2);
    }
    return normaCuadrado;
}

Householder::Householder(matrizDouble& A, unsigned int columna) : u(vectorDouble(A.size() - columna, 0)), columna(columna) {
    double normaCuadrado = auxSubNorma2CuadradoColumna(A, columna);
    double norma = sqrt(normaCuadrado);
    double norma_u = sqrt(2 * normaCuadrado - 2 * A[columna][columna] * norma);
    this->u[0] = (A[columna][columna] - norma) / norma_u;
    for (unsigned int i = 1; i < this->u.size(); i++) {
        this->u[i] = (A[columna + i][columna]) / norma_u;
    }
}

Householder::~Householder() {
}

void Householder::ut_por_A(matrizDouble& A, vectorDouble& res) {
    unsigned int m = A.size() - this->columna;
    unsigned int n = A[0].size() - this->columna;
    assert(n == res.size());
    for (unsigned int i = 0; i < n; i++) {
        res[i] = 0;
        for (unsigned int j = 0; j < m; j++) {
            res[i] += this->u[j] * A[this->columna + j][this->columna + i];
        }
    }
}

double Householder::ut_por_v(vectorDouble& v) {
    unsigned int m = v.size() - this->columna;
    double ret = 0;
    for (unsigned int i = 0; i < m; i++) {
        ret += this->u[i] * v[this->columna + i];
    }
    return ret;
}

void Householder::aplicarM(matrizDouble& matriz) {
    unsigned int m = matriz.size() - this->columna;
    unsigned int n = matriz[0].size() - this->columna;
    vectorDouble ut_x_A = vectorDouble(n, 0);
    ut_por_A(matriz, ut_x_A);
    double calculo;
    for (unsigned int i = 0; i < m; i++) {
        for (unsigned int j = 0; j < n; j++) {
            calculo = matriz[this->columna + i][this->columna + j] - 2 * this->u[i] * ut_x_A[j];
            if (fabs(calculo) < epsilon) matriz[this->columna + i][this->columna + j] = 0;
            else matriz[this->columna + i][this->columna + j] = calculo;
        }
    }
}

void Householder::aplicarV(vectorDouble& vec) {
    unsigned int m = vec.size() - this->columna;
    assert(m == this->u.size());
    double ut_x_v = ut_por_v(vec);
    for (unsigned int i = 0; i < m; i++) {
        vec[this->columna + i] = this->u[i] * ut_x_v;
    }
}

void Householder::transponer() {
}