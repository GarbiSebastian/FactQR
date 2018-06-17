#ifndef CML_H
#define CML_H

#include "typedefs.h"
#include "FactorizadorQR.h"


class CML {
public:
    CML(matrizDouble& A);
    virtual ~CML();
    void resolver(vectorDouble& b, vectorDouble& res);
private:
    FactorizadorQR q;
    matrizDouble* A;
};

#endif /* CML_H */

