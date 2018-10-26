/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: 
 * @Date  : 03.03.2018
 */

#include <cmath>
#include <limits>
#include <functional>
#include <iostream>
#include <algorithm>

#include "LUDoolittle.hpp"
#include "Matrix.hpp"

#ifndef ANPI_LU_HPP
#define ANPI_LU_HPP

namespace anpi
{


/***
 *  Uses the doolittle algorithm every time because the graphic
 * shows that it is the more efficient method for any size
 * 
 * @param[in] A = Matriz de entrada 
 * @param[out] LU = Matriz de salida con LU
 * @param[out] p = Vector de permutaciones 
 * */
template <typename T>
inline void lu(const anpi::Matrix<T> &A,
               anpi::Matrix<T> &LU,
               std::vector<size_t> &p)
{
    anpi::luDoolittle(A, LU, p);
    
}

} //anpi

#endif
