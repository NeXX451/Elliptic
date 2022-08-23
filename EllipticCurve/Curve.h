//
// Created by Tomas Duchac.
//

#ifndef ELLIPTIC_CURVE_H
#define ELLIPTIC_CURVE_H


#include "Common.h"

class Curve {

public:

    Curve();

    /**
     * Create a curve
     * @param primeField prime for the field, Z_primeField
     * @param A
     * @param B
     */
    Curve(const Integer& primeField, const Integer& A, const Integer& B, CurveType type);
    // Edwards Curve
    Curve(const Integer& primeField, const Integer& d, CurveType type);

    ~Curve();

    bool isZeroDiscriminant();

    ZP *getField()
    {
        return this->FField;
    }

    Element getA() const
    {
        return this->A;
    }

    Element getB() const
    {
        return this->B;
    }


    const Element &getDiscriminant() const;

    const Element &getJInvariant() const;

    CurveType getCurveType() const;

    void printField();
    void printDiscriminant();
    void printJInvariant();

protected:

    Element A,B;
    Element Discriminant{};
    Element jInvariant{};

    void computeDiscriminant();
    void compute_jInvariant();

    CurveType curveType;

    ZP* FField{};

    Element MontgomeryA, MontgomeryB;
    Element EdwardsD;

};

#endif //ELLIPTIC_CURVE_H
