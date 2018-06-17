#ifndef GIVENS_H
#define GIVENS_H

#include "typedefs.h"
using namespace std;

class Givens {
public:
    Givens(const matrizDouble& A, unsigned int i, unsigned int j);
    virtual ~Givens();
    void aplicarM(matrizDouble& matriz);
    void aplicarV(vectorDouble& vec);
    void transponer();
private:
    double x1;
    double x2;
    double norma;
    unsigned int i;
    unsigned int j;

};

#endif /* GIVENS_H */

