#include <cassert>
#include <iostream>

#include "CML.h"

CML::CML(matrizDouble& A) : A(&A), q(FactorizadorQR(A)) {
    this->q.transponer();
}

void CML::resolver(vectorDouble& b, vectorDouble & respuesta) {
    unsigned int m = this->A->size();
    unsigned int n = (*(this->A))[0].size();
    unsigned int i;
    double acum, a_ii, a_ij;
    assert(n==respuesta.size());
    this->q.aplicarV(b);
    for (unsigned int _i = 0; _i < n; _i++) {
        i = n - 1 - _i;
        a_ii = (*(this->A))[i][i];
        acum = 0;
        for (unsigned int j = i + 1; j < n; j++) {
            a_ij =  (*(this->A))[i][j];
            acum += b[j] * a_ij/a_ii;
        }
        respuesta[i] = b[i]/a_ii - acum;
    }
}

CML::~CML() {
}

