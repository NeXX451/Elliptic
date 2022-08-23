//
// Created by Tomas Duchac.
//

#include "Curve.h"


Curve::Curve() = default;


Curve::~Curve()
{
    delete FField;
}

Curve::Curve(const Integer& primeField, const Integer& A, const Integer& B, CurveType type) : curveType(type)
{
    if (type == SHORT_WEIERSTRASS){
        this->FField = new ZP(primeField);
        this->FField->init(this->A,A);
        this->FField->init(this->B,B);
    } else if (type == MONTGOMERY){
        this->MontgomeryA = A;
        this->MontgomeryB = B;
        /**
         * Convert Mongtgomery Curve into short Weierstrass form
         */
        // A = (3 - A^2)/(3 * B^2)
        // B = (2 * A^3 - 9 * A)/(27 * B^3)
        this->FField = new ZP(primeField);

        Element Asqr, Bsqr;
        this->FField->mul(Asqr,A,A);
        this->FField->sub(Asqr,Integer("3"),Asqr);
        this->FField->mul(Bsqr,B,B);
        this->FField->mul(Bsqr,Integer("3"),Bsqr);
        Element AA = this->FField->div(this->A,Asqr,Bsqr);

        Element Acb, A9, Bcb;
        //2*A^3
        this->FField->mul(Acb,A,A);
        this->FField->mulin(Acb,A);
        this->FField->mul(Acb, Integer("2"), Acb);
        // 2*A^3 - 9A
        this->FField->mul(A9,Integer("9"),A);
        this->FField->sub(Acb,Acb,A9);
        //27*B^3
        this->FField->mul(Bcb,B,B);
        this->FField->mulin(Bcb,B);
        this->FField->mulin(Bcb,Integer("27"));
        Element BB = this->FField->div(this->B,Acb,Bcb);

        this->FField->init(this->A,AA);
        this->FField->init(this->B,BB);

    }


    computeDiscriminant();
    compute_jInvariant();
    if (this->isZeroDiscriminant()){
        std::cerr << "[!] Curve not defined, disriminant is Zero" << std::endl;
        std::cout << "[+ INFO ] Field : Z/" << primeField << "Z" << std::endl;
        std::cout << "[+ INFO ] A : " << A << " ; B : " << B << std::endl;
        abort();
    }
}

/**
 * A = 2(1+d)/(1−d)
 * B = 4/(1−d)
 * @param primeField
 * @param d
 * @param type
 */
Curve::Curve(const Integer &primeField, const Integer &d, CurveType type) : curveType(type)
{
    this->FField = new ZP(primeField);

    Element oneplusd, oneminusd;
    this->FField->add(oneplusd,Integer("1"), d);
    this->FField->sub(oneminusd,Integer("1"),d);

    Element oneplusd2;
    this->FField->mul(oneplusd2, Integer("2"), oneplusd);

    this->FField->div(this->A,oneplusd2,oneminusd);
    this->FField->div(this->B, Integer("4"), oneminusd);

    this->FField->init(this->A,this->A);
    this->FField->init(this->B,this->B);

    computeDiscriminant();
    compute_jInvariant();
    if (this->isZeroDiscriminant()){
        std::cerr << "[!] Curve not defined, disriminant is Zero" << std::endl;
        std::cout << "[+ INFO ] Field : Z/" << primeField << "Z" << std::endl;
        std::cout << "[+ INFO ] A : " << A << " ; B : " << B << std::endl;
        abort();
    }
}



void Curve::printField()
{
    std::cout << " Field : Z/" << this->FField->residu() << "Z" << std::endl;
    std::cout << " A :" ;
    this->FField->write(std::cout,this->A);
    std::cout << ", B :";
    this->FField->write(std::cout,this->B);
    std::cout<< std::endl;
}

void Curve::computeDiscriminant()
{

    //// D = −16(4A^3 + 27B^2)
    Integer a = this->A,
            b = this->B;

    Integer Acb, Bsq, Acb4, Bsq27;

    this->FField->mul(Acb, a,a);
    this->FField->mulin(Acb, a);
    this->FField->mul(Acb4, Acb, Integer("4"));

    this->FField->mul(Bsq,b,b);
    this->FField->mul(Bsq27,Bsq, Integer("27"));

    Element r;
    this->FField->add(r,Acb4,Bsq27);

    Element neg16;
    this->FField->sub(neg16,FField->zero, Integer("16"));
    this->FField->mul(this->Discriminant,neg16, r);
}

void Curve::compute_jInvariant()
{
    Integer a = this->A,
            b = this->B;

    Element C, Ccube;
    this->FField->sub(C,this->FField->zero,a); // -A
    this->FField->mulin(C, Integer("48")); // -A*48
    this->FField->mul(Ccube,C, C);
    this->FField->mulin(Ccube,C); // C^3

   this->FField->div(this->jInvariant,Ccube,this->Discriminant);
}


bool Curve::isZeroDiscriminant()
{
    return this->FField->isZero(this->Discriminant);
}

void Curve::printDiscriminant()
{
    std::cout << "Elliptic Discriminant: " << this->Discriminant << std::endl;
}

void Curve::printJInvariant()
{
    /**
     * The j-Invariant for an EC in the Weierstrass form is defined as
     * j = c_4^3 / D
     * where D is the Elliptic Discriminant and c_4^3 = b_2^2 - 24b_4
     */
    std::cout << "j-Invariant: " << this->jInvariant << std::endl;
}

const Element &Curve::getDiscriminant() const {
    return Discriminant;
}

const Element &Curve::getJInvariant() const {
    return jInvariant;
}

CurveType Curve::getCurveType() const {
    return curveType;
}




