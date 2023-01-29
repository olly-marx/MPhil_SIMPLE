//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk

#ifndef VECTOR_HEADER
#define VECTOR_HEADER

#include <iostream>
#include <array>
#include <string>
#include <vector>

std::array<double,3> cross(std::array<double,3>, std::array<double,3>);
std::array<double,3> sum(std::array<double,3>, std::array<double,3>);
std::array<double,3> diff(std::array<double,3>, std::array<double,3>);
std::array<double,3> scalarMult(double, std::array<double,3>);
double               dot(std::array<double,3>, std::array<double,3>);
double               mod(std::array<double,3>);

#endif
