#ifndef FACTORIZADORQR_H
#define FACTORIZADORQR_H

#include "typedefs.h"
#include "Givens.h"
#include "Householder.h"

using namespace std;
typedef vector<Givens> vectorGivens;
typedef vector<Householder> vectorHouseholder;

class FactorizadorQR {
public:
    FactorizadorQR(matrizDouble& A);
    virtual ~FactorizadorQR();
    void aplicarM(matrizDouble& matriz);
    void aplicarV(vectorDouble& vec);
    void transponer();
private:
    vectorGivens gs;
    vectorHouseholder hs;
    bool transpuesta;
};

#endif /* FACTORIZADORQR_H */

