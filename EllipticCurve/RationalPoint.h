//
// Created by Tomas Duchac.
//

#ifndef ELLIPTIC_RATIONALPOINT_H
#define ELLIPTIC_RATIONALPOINT_H


#include "Common.h"

class RationalPoint {
public:

    RationalPoint() : identity(true) {}
    explicit RationalPoint(bool b);
    RationalPoint(const Element& x, const Element& y);

    Element getX() const
    {
        return this->x;
    }

    void setX(const Element& _x)
    {
        this->x = _x;
    }

    Element getY() const
    {
        return this->y;
    }

    void setY(const Element& _y)
    {
        this->y = _y;
    }

    bool isIdentity() const
    {
        return this->identity;
    }

    void setIdentity(bool _identity)
    {
        this->identity = _identity;
    }

    RationalPoint& operator=(const RationalPoint& P);

    bool operator==(const RationalPoint& P) const;

    void print();

private:

    bool identity;
    Element x,y;

};


#endif //ELLIPTIC_RATIONALPOINT_H
