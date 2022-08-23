//
// Created by Tomas Duchac
//

#ifndef ELLIPTIC_COMMON_H
#define ELLIPTIC_COMMON_H


#include <givaro/modular-integer.h>

using namespace Givaro;

// Z_p, modular integer
typedef Modular<Integer> ZP;
// element of Z_p
typedef ZP::Element Element;

enum CurveType
{
    SHORT_WEIERSTRASS,
    MONTGOMERY,
    EDWARD
};


#endif //ELLIPTIC_COMMON_H
