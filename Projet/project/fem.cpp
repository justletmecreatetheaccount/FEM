#include "project/fem.hpp"
#include <array>

/*
static const std::array<double, 3> _gaussTri3Xsi = {
    0.166666666666667, 0.666666666666667, 0.166666666666667};
static const std::array<double, 3> _gaussTri3Eta = {
    0.166666666666667, 0.166666666666667, 0.666666666666667};
static const std::array<double, 3> _gaussTri3Weight = {
    0.166666666666667, 0.166666666666667, 0.166666666666667};
static const std::array<double, 2> _gaussEdge2Xsi = {0.577350269189626,
                                                     -0.577350269189626};
static const std::array<double, 2> _gaussEdge2Weight = {1.000000000000000,
                                                        1.000000000000000};
*/
fem::Node::Node() : x(0), y(0) {};
fem::Node::Node(double x_n, double y_n) : x(x_n), y(y_n) {};
