//
// Created by Tomas Duchac.
//


#include <iostream>
#include "EllipticCurve.h"

using namespace std;
using namespace Givaro;
int main() {

    /**
     * secure
     * M-221 : y^2 = x^3+117050x^2+x  modulo p = 2^221 - 3 Montgomery curve. https://neuromancer.sk/std/other/M-221
     * E-222 : x^2+y^2 = 1+160102x^2y^2  modulo p = 2^222 - 117 Edwards curve. https://neuromancer.sk/std/other/E-222
     * Curve-1174 : x^2+y^2 = 1-1174x^2y^2  modulo p = 2^251 - 9 Weierstrass curve. https://neuromancer.sk/std/other/Curve1174
     * Curve25519 : y^2 = x^3+486662x^2+x modulo p = 2^255 - 19 Montgomery curve. https://neuromancer.sk/std/other/Curve25519
     * E-382 : x^2+y^2 = 1-67254x^2y^2  modulo p = 2^382 - 105
     * M-383 : y^2 = x^3+2065150x^2+x  modulo p = 2^383 - 187
     * Curve383187 : y^2 = x^3+229969x^2+x  modulo p = 2^383 - 187
     * Ed448-Goldilocks : x^2+y^2 = 1-39081x^2y^2  modulo p = 2^448 - 2^224 - 1
     * M-511 : y^2 = x^3+530438x^2+x  modulo p = 2^511 - 187
     * E-521 : x^2+y^2 = 1-376014x^2y^2  modulo p = 2^521 - 1
     */

    /*// https://neuromancer.sk/std/nist/P-384#
    Integer Mod = "39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319",
            A = "39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112316",
            B = "27580193559959705877849011840389048093056905856361568521428707301988689241309860865136260764883745107765439761230575";
    //Discriminant
    Integer Discriminant("38275261264050278989862136034342276004573039492779555073863190029182890449044186682105480613137214197175883602718257");
    //j-Invariant
    Integer jInvariant("12550029517991417762405079599420518784762671286028430215113399824456742172589190955698027499893480133182923443701083");*/



    //Curve25519 : y^2 = x^3+486662x^2+x modulo p = 2^255 - 19 Montgomery curve.
    /*Integer Mod = "57896044618658097711785492504343953926634992332820282019728792003956564819949",
            A = "486662",
            B = "1";
    //Discriminant
    Integer Discriminant("3789438435840");
    //j-Invariant
    Integer jInvariant("39240375672115510010799456308813573486606784421612167109713554819120306934551");

    EllipticCurve ecc(Mod, A, B, MONTGOMERY);

    ecc.print();
    ecc.printField();
    ecc.printDiscriminant();
    ecc.printJInvariant();
    if (ecc.getDiscriminant() == Discriminant && ecc.getJInvariant() == jInvariant)
        cout << "Discriminant and j-Invariant computed correctly!" << endl;
    else
        cout << "Something went wrong!" << endl;*/

    // E-222 : x^2+y^2 = 1+160102x^2y^2  modulo p = 2^222 - 117 Edwards curve. https://neuromancer.sk/std/other/E-222
    Integer Mod = "6739986666787659948666753771754907668409286105635143120275902562187",
            d = "160102";
    //Discriminant
    Integer Discriminant("3789438435840");
    //j-Invariant
    Integer jInvariant("39240375672115510010799456308813573486606784421612167109713554819120306934551");

    EllipticCurve ecc(Mod, d, EDWARD);

    ecc.print();
    ecc.printField();
    ecc.printDiscriminant();
    ecc.printJInvariant();
    if (ecc.getDiscriminant() == Discriminant && ecc.getJInvariant() == jInvariant)
        cout << "Discriminant and j-Invariant computed correctly!" << endl;
    else
        cout << "Something went wrong!" << endl;



    /*EllipticCurve eccTest(Integer("115792089237316195423570985008687907853269984665640564039457584007908834671663"),
                          Integer("0"),Integer("7"), SHORT_WEIERSTRASS);

    eccTest.print();

    RationalPoint P(Integer("65485170049033141755572552932091555440395722142984193594072873483228125624049"),
                    Integer("73163377763031141032501259779738441094247887834941211187427503803434828368457"));

    if (eccTest.verifyPoint(P)){
        P.print();
        std::cout << "Point : verified!!!!" << std::endl;
    };

    RationalPoint P2;
    eccTest._inv(P2,P);
    if (eccTest.verifyPoint(P2)){
        P2.print();
        std::cout << "Inverse : verified" << std::endl;
    };

    RationalPoint R;
    R = eccTest._double(R, P);
    if (eccTest.verifyPoint(R)){
        R.print();
        std::cout << "double : verified" << std::endl;
    };

    RationalPoint Rp ;
    Rp = eccTest._add(Rp, P , P);
    if (eccTest.verifyPoint(Rp)){
        Rp.print();
        std::cout << "add : verified" << std::endl;
    };

    assert(Rp == R);

    RationalPoint R5;
    R5 = eccTest._scalar(R5,P,Integer("13257"));
    if (eccTest.verifyPoint(R5)){
        R5.print();
        std::cout << "scalar : verified" << std::endl;
    };*/

    return 0;
}