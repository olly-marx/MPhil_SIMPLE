//AUTHOR:       Oliver Marx - ojm40@cam.ac.uk

#ifndef VECTOR_HEADER
#define VECTOR_HEADER

#include <iostream>
#include <array>
#include <string>
#include <vector>

std::array<double,3> cross(const std::array<double,3>&, const std::array<double,3>&);
std::array<double,3> operator+(const std::array<double,3>&, const std::array<double,3>&);
std::array<double,3> operator-(const std::array<double,3>&, const std::array<double,3>&);
std::array<double,3> operator*(const double&, const std::array<double,3>&);
std::array<double,3> operator*(const std::array<double,3>&, const double&);
double               dot(std::array<double,3>, std::array<double,3>);
double               mod(std::array<double,3>);
std::string          printV(std::array<double,3>);

#endif
