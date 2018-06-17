#ifndef HOUSEHOLDER_H
#define HOUSEHOLDER_H

#include "typedefs.h"
using namespace std;

class Householder {
public:
    Householder(vectorDouble& u);
    Householder(matrizDouble& A, unsigned int columna);
    virtual ~Householder();
    void aplicarM(matrizDouble& matriz);
    void aplicarV(vectorDouble& vec);
    void transponer();
protected:
    void ut_por_A(matrizDouble& A, vectorDouble& res);
    double ut_por_v(vectorDouble& v);
private:
    vectorDouble u;
    unsigned int columna;
};

#endif /* HOUSEHOLDER_H */

