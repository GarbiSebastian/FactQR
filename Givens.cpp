#include <math.h>
#include <complex>
#include <cassert>
#include <iostream>
#include "Givens.h"
#include "constantes.h"

Givens::Givens(const matrizDouble& A, unsigned int i, unsigned int j) {
    assert(i < j);
    assert(j < A.size());
    double x1 = A[i][i];
    double x2 = A[j][i];
    this->i = i;
    this->j = j;
    this->norma = sqrt(pow(x1, 2) + pow(x2, 2));
    this->x1 = x1 / this->norma;
    this->x2 = x2 / this->norma;
}

Givens::~Givens() {

}

void Givens::aplicarM(matrizDouble& matriz) {
    unsigned int m = matriz.size();
    assert(m > 0);
    unsigned int n = matriz[0].size();
    /*matriz[this->i][this->i]=this->norma;
    matriz[this->j][this->i]=0;*/
    double aik, ajk, calculo;
    //    for (unsigned int k = this->i+1; k < n; k++) {
    for (unsigned int k = this->i; k < n; k++) {
        aik = matriz[this->i][k];
        ajk = matriz[this->j][k];
        calculo = this->x1 * aik + this->x2*ajk;
        if (fabs(calculo) < epsilon) matriz[this->i][k] = 0;
        else matriz[this->i][k] = calculo;
        calculo = this->x1 * ajk - this->x2*aik;
        if (fabs(calculo) < epsilon) matriz[this->j][k] = 0;
        else matriz[this->j][k] = calculo;
    }
}

void Givens::aplicarV(vectorDouble& vec) {
    unsigned int m = vec.size();
    assert(m > 0);
    double vi, vj;
    vi = vec[this->i];
    vj = vec[this->j];
    vec[this->i] = this->x1 * vi + this->x2*vj;
    vec[this->j] = this->x1 * vj - this->x2*vi;
}

void Givens::transponer() {
    this->x2 = -this->x2;
}


