//
// Created by Tomas Duchac.
//

#include "RationalPoint.h"
#include "common.h"


RationalPoint::RationalPoint(bool b) : identity(b) {}


RationalPoint::RationalPoint(const Element& x, const Element& y) :
    identity(false), x(x), y(y)
    {}


RationalPoint& RationalPoint::operator =(const RationalPoint& P)
{
    if(P.identity)
    {
        this->identity = true;
        return *this;
    }

    this->identity = P.identity;
    this->x = P.x;
    this->y = P.y;
    return *this;
}


bool RationalPoint::operator ==(const RationalPoint& P) const
{
    return (this->identity && P.identity) ||
           (!this->identity && !P.identity && this->x == P.x && this->y == P.y);
}


void RationalPoint::print()
{
    this->identity ? std::cout << "infinity "<< std::endl :
                     std::cout << "[" << this->x << "," << this->y << "]" << std::endl;
}