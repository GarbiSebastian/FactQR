#include <cmath>
#include <iostream>

#include "FactorizadorQR.h"
#include "constantes.h"
#include "Givens.h"
#include "Householder.h"

FactorizadorQR::FactorizadorQR(matrizDouble& A) {
    this->transpuesta = false;
    unsigned int m = A.size();
    unsigned int n = A[0].size();
    for (unsigned int i = 0; i < n - 1; i++) {
        for (unsigned int j = i + 1; j < m; j++) {
            if (fabs(A[j][i]) > epsilon) {
                Givens g(A, i, j);
                g.aplicarM(A);
                this->gs.push_back(g);
            }
        }
    }
    if (m > n) {
        Householder h(A, n - 1);
        h.aplicarM(A);
        this->hs.push_back(h);
    }
}

FactorizadorQR::~FactorizadorQR() {
}

void FactorizadorQR::aplicarM(matrizDouble& matriz) {
    unsigned int hs_size = this->hs.size();
    unsigned int gs_size = this->gs.size();
    if (this->transpuesta) {
        for (unsigned int i = 0; i < hs_size; i++) {//debe haber a lo sumo una para la última columna que no se puede resolver por Givens
            Householder h = this->hs[hs_size-1 - i];
            h.aplicarM(matriz);
        }
        for (unsigned int i = 0; i < gs_size; i++) {
            Givens g = this->gs[gs_size -1 - i];
            g.aplicarM(matriz);
        }
    } else {
        for (unsigned int i = 0; i < gs_size; i++) {
            Givens g = this->gs[i];
            g.aplicarM(matriz);
        }
        for (unsigned int i = 0; i < hs_size; i++) {
            Householder h = this->hs[i];
            h.aplicarM(matriz);
        }
    }
}

void FactorizadorQR::aplicarV(vectorDouble& v) {
    unsigned int hs_size = this->hs.size();
    unsigned int gs_size = this->gs.size();
    if (this->transpuesta) {
        for (unsigned int i = 0; i < hs_size; i++) {//debe haber a lo sumo una para la última columna que no se puede resolver por Givens
            this->hs[hs_size -1 - i].aplicarV(v);
        }
        for (unsigned int i = 0; i < gs_size; i++) {
            this->gs[gs_size -1 - i].aplicarV(v);
        }
    } else {
        for (unsigned int i = 0; i < gs_size; i++) {
            this->gs[i].aplicarV(v);
        }
        for (unsigned int i = 0; i < hs_size; i++) {
            this->hs[i].aplicarV(v);
        }
    }
}

void FactorizadorQR::transponer() {
    this->transpuesta = !this->transpuesta;
    unsigned int hs_size = this->hs.size();
    unsigned int gs_size = this->gs.size();
    for (unsigned int i = 0; i < gs_size; i++) {
        this->gs[i].transponer();
    }
    for (unsigned int i = 0; i < hs_size; i++) {
        this->hs[i].transponer();
    }
}