cmake_minimum_required(VERSION 3.22)
project(Elliptic)

set(CMAKE_CXX_STANDARD 11)




add_executable(Elliptic
        EllipticCurve/Common.h EllipticCurve/main.cpp EllipticCurve/RationalPoint.cpp EllipticCurve/RationalPoint.h EllipticCurve/Curve.cpp EllipticCurve/Curve.h EllipticCurve/EllipticCurve.cpp EllipticCurve/EllipticCurve.h)
target_link_libraries(Elliptic PUBLIC givaro gmp)


