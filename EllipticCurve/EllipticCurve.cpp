//
// Created by Tomas Duchac.
//

#include "EllipticCurve.h"
/**
 * Montgomery or Weierstrass Curve
 * @param primeField
 * @param A
 * @param B
 * @param type
 */
EllipticCurve::EllipticCurve(const Integer& primeField, const Integer& A, const Integer& B, CurveType type) :
    Curve(primeField,A,B,type)
{
    this->identity.setIdentity(true);
}

/**
 * Edward's Curve
 * @param primeField
 * @param d
 * @param type
 */
EllipticCurve::EllipticCurve(const Integer &primeField, const Integer &d, CurveType type) :
    Curve(primeField, d, type)
{
    this->identity.setIdentity(true);
}


EllipticCurve::~EllipticCurve() = default;


/**
 * Map Montgomery curve point(s) to Weierstrass
 * @param Q point to map Q(x,y) |-> Q((x+A/3)/B, y/B).
 */
void EllipticCurve::convertPointsMontgomery_Weierstrass(RationalPoint &Q) {
    // Set x
    Element A3, Xnew;
    this->FField->div(A3,this->MontgomeryA, Integer("3"));
    this->FField->add(Xnew, Q.getX(),A3);
    this->FField->divin(Xnew, this->MontgomeryB);
    Q.setX(Xnew);

    // Set y
    Element Ynew;
    this->FField->div(Ynew,Q.getY(),this->MontgomeryB);
    Q.setY(Ynew);
}


const RationalPoint& EllipticCurve::_inv(RationalPoint& R, RationalPoint& P)
{
    if (this->curveType == MONTGOMERY || this->curveType == SHORT_WEIERSTRASS) {
        if (this->curveType == MONTGOMERY) {
            convertPointsMontgomery_Weierstrass(P);
        }

        if (P.isIdentity()) {
            R.setIdentity(true);
            return R;
        } else {
            R.setIdentity(false);
            R.setX(P.getX());
            Element tmp;
            this->FField->sub(tmp, this->FField->zero, P.getY()); // 0 - Point = - Point --> (R,P) -> (R, -P)

            R.setY(tmp);
            return R;
        }
    } else {
        /**
         * Edwards Curve
         * −(x1, y1) = (−x1, y1),
         */
        if (P.isIdentity()) {
            R.setIdentity(true);
            return R;
        }

        R.setIdentity(false);
        R.setY(P.getY());
        Element tmp;
        this->FField->sub(tmp, this->FField->zero, P.getX());
        R.setX(tmp);
        return R;

    }
}


/**
 * For Addition if P == Q
 * @param R
 * @param P
 * @return
 */
RationalPoint& EllipticCurve::_double(RationalPoint& R,  RationalPoint& P)
{
    if (this->curveType == MONTGOMERY || this->curveType == SHORT_WEIERSTRASS) {
        if (this->curveType == MONTGOMERY) {
            convertPointsMontgomery_Weierstrass(P);
        }

        RationalPoint tmp;
        this->_inv(tmp, P);
        if (P.isIdentity() || (tmp == P)) {
            R.setIdentity(true);
            return R;
        } else {

            R.setIdentity(false);
            Element tmp;
            Element xs; // 3*x^2
            this->FField->mul(xs, P.getX(), P.getX());
            this->FField->init(tmp, (uint64_t) 3);
            this->FField->mulin(xs, tmp);
            Element ty; // 2*y
            this->FField->init(tmp, (uint64_t) 2);
            this->FField->mul(ty, P.getY(), tmp);
            Element slope; // m
            this->FField->add(tmp, xs, getA());
            this->FField->div(slope, tmp, ty);
            Element slope2; // m^2
            this->FField->mul(slope2, slope, slope);
            Element tx; // 2x
            this->FField->add(tx, P.getX(), P.getX());
            Element x3; // x_3
            this->FField->sub(x3, slope2, tx);
            Element y3; // y_3
            this->FField->sub(tmp, P.getX(), x3);
            this->FField->sub(y3, this->FField->mulin(tmp, slope), P.getY());

            R.setX(x3);
            R.setY(y3);
            return R;
        }
    } else {
        /**
         * Edward's Curve
         * 2(x1, y1) = ( (2xy) / (1 + dx^2y^2) , (y^2 − x^2) / (1 − dx^2y^2) )
         */
        RationalPoint tmp;
        this->_inv(tmp, P);
        if (P.isIdentity() || (tmp == P)) {
            R.setIdentity(true);
            return R;
        }
        Element xy, twoxy, xsq, ysq, xsqysq, dxsqysq, oneplusdxsqysq, ysqminusxsq, oneminusdxsqysq;
        Element x = P.getX();
        Element y = P.getY();
        this->FField->mul(xy,x,y);
        this->FField->mul(twoxy,Integer("2"), xy);
        this->FField->mul(xsq, x, x);
        this->FField->mul(ysq, y, y);
        this->FField->mul(xsqysq,xsq,ysq);
        this->FField->mul(dxsqysq, this->EdwardsD,xsqysq);
        this->FField->add(oneplusdxsqysq,this->FField->one,dxsqysq);
        this->FField->sub(oneminusdxsqysq,this->FField->one,dxsqysq);
        this->FField->sub(ysqminusxsq,ysq,xsq);

        Element X, Y;
        this->FField->div(X,twoxy,oneplusdxsqysq);
        this->FField->div(Y,ysqminusxsq,oneminusdxsqysq);

        R.setX(X);
        R.setY(Y);
        return R;
    }
}


RationalPoint& EllipticCurve::_add(RationalPoint& R, RationalPoint& P, RationalPoint& Q)
{
    if (this->curveType== MONTGOMERY || this->curveType == SHORT_WEIERSTRASS) {
        if (this->curveType == MONTGOMERY) {
            convertPointsMontgomery_Weierstrass(P);
            convertPointsMontgomery_Weierstrass(Q);
        }

        RationalPoint tmp1, tmp2;
        this->_inv(tmp1, P);
        this->_inv(tmp2, Q);
        if ((P.isIdentity() && Q.isIdentity()) || (tmp1 == Q) || (tmp2 == P)) {
            R.setIdentity(true);
            return R;
        }
        if (P.isIdentity()) {
            R = Q;
            return R;
        }
        if (Q.isIdentity()) {
            R = P;
            return R;
        }
        if (P == Q) {
            return this->_double(R, P);
        }

        R.setIdentity(false);
        Element tmp;
        Element num; // y2 - y1
        this->FField->sub(num, Q.getY(), P.getY());
        Element den; // x2 - x1
        this->FField->sub(den, Q.getX(), P.getX());
        Element slope; // m
        this->FField->div(slope, num, den);
        Element slope2; // m^2
        this->FField->mul(slope2, slope, slope);
        Element x3; // x_3
        this->FField->sub(x3, slope2, this->FField->add(tmp, Q.getX(), P.getX()));
        Element diffx3; // x_1 - x_3
        this->FField->sub(diffx3, P.getX(), x3);
        Element y3; // y_3
        this->FField->mul(tmp, slope, diffx3);
        this->FField->sub(y3, tmp, P.getY());

        R.setX(x3);
        R.setY(y3);
        return R;
    } else {

        /**
         * Edwards Curve
         * (x1, y1) + (x2, y2) = (x3, y3) = ( (x1y2 + x2y1) / (1 + dx1x2y1y2) , (y1y2 − x1x2) /( 1 − dx1x2y1y2 ))
         */

        if (P.isIdentity()) {
            R = Q;
            return R;
        }
        if (Q.isIdentity()) {
            R = P;
            return R;
        }
        if (P == Q) {
            return this->_double(R, P);
        }


        Element x1y2, x2y1;
        this->FField->mul(x1y2, P.getX(),Q.getY());
        this->FField->mul(x2y1, Q.getX(),P.getY());
        Element x1y2plusx2y1;
        this->FField->add(x1y2plusx2y1,x1y2,x2y1);
        Element x1x2, y1y2;
        this->FField->mul(x1x2, P.getX(), Q.getX());
        this->FField->mul(y1y2, P.getY(), Q.getY());
        Element x1x2y1y2;
        this->FField->mul(x1x2y1y2,x1x2,y1y2);
        Element dx1x2y1y2;
        this->FField->mul(dx1x2y1y2, this->EdwardsD, x1x2y1y2);
        Element oneplusdx1x2y1y2;
        this->FField->add(oneplusdx1x2y1y2, this->FField->one, dx1x2y1y2);
        Element y1y2minusx1x2;
        this->FField->sub(y1y2minusx1x2, y1y2, x1x2);
        Element oneminusdx1x2y1y2;
        this->FField->sub(oneminusdx1x2y1y2, this->FField->one, dx1x2y1y2);

        Element X3, Y3;
        this->FField->div(X3, x1y2plusx2y1, oneplusdx1x2y1y2);
        this->FField->div(Y3, y1y2minusx1x2, oneminusdx1x2y1y2);

        R.setX(X3);
        R.setY(Y3);

        return R;
    }
}


RationalPoint& EllipticCurve::_scalar(RationalPoint& R, RationalPoint& P, Integer k)
{

    if(curveType == MONTGOMERY){
        convertPointsMontgomery_Weierstrass(P);
    }

    if(P.isIdentity())
    {
        R.setIdentity(true);
        return R;
    }
    RationalPoint tmp1, tmp2;
    R.setIdentity(true);
    RationalPoint PP = P;
    while(k > 0){
        if (k % 2 == 1){
            this->_add(tmp1,R,PP);
            R = tmp1;
        }
        tmp2 = this->_double(tmp2,PP);
        PP = tmp2;
        k >>= 1;
    }
    return R;
}


bool EllipticCurve::verifyPoint(const RationalPoint& P) const
{
    // y^2 = x^3 + ax + b
    if(P.isIdentity()){
        return true;
    }
    Element x3,y2;
    Element Ax,rhs;
    this->FField->mul(x3,P.getX(),P.getX());
    this->FField->mulin(x3,P.getX()); // x^3
    this->FField->mul(Ax,P.getX(),getA()); // a*x

    this->FField->add(rhs,x3,Ax); // x^3 + (a * x)
    this->FField->addin(rhs,getB()); // x^3 + (a * x) + b


    this->FField->mul(y2,P.getY(),P.getY()); // y^2
    return y2==rhs;
}


void EllipticCurve::print() {
    std::cout<<"Elliptic Curve Defined by ";
    std::cout<<"y^2 = x^3 + ";
    this->FField->write(std::cout,getA());
    std::cout<<"x + ";
    this->FField->write(std::cout,getB());
    std::cout<<std::endl;
}


