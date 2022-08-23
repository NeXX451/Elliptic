//
// Created by Tomas Duchac.
//

#ifndef ELLIPTIC_ELLIPTICCURVE_H
#define ELLIPTIC_ELLIPTICCURVE_H

#include "Common.h"
#include "Curve.h"
#include "RationalPoint.h"

class EllipticCurve: public Curve{

public:
    EllipticCurve(const Integer& primeField, const Integer& A, const Integer& B, CurveType type);
    EllipticCurve(const Integer& primeField, const Integer& d, CurveType type);
    ~EllipticCurve();

    const RationalPoint& _inv(RationalPoint& R, RationalPoint& P);
    RationalPoint& _double(RationalPoint& R,  RationalPoint& P);
    RationalPoint& _add(RationalPoint& R, RationalPoint &P,  RationalPoint& Q);
    RationalPoint& _scalar(RationalPoint& R,  RationalPoint& P,Integer k);

    bool verifyPoint(const RationalPoint& P) const;
    void print();

private:
    RationalPoint identity;
    void convertPointsMontgomery_Weierstrass(RationalPoint& Q);

};
#endif //ELLIPTIC_ELLIPTICCURVE_H
