#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <vector>

#ifdef __cplusplus
extern "C" {
#endif
    typedef std::vector<double> vectorDouble;
    typedef std::vector< vectorDouble > matrizDouble;

#ifdef __cplusplus
}
#endif

#endif /* TYPEDEFS_H */

