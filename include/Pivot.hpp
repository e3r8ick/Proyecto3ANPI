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

#include "Exception.hpp"
#include "Matrix.hpp"

#ifndef ANPI_PIVOT_HPP
#define ANPI_PIVOT_HPP

namespace anpi
{


//prints a matrix
/**
 * Imprime una matriz
 * @param[in] M = Matriz de entrada para imprimir
 * 
 * */
template <typename T>
void printM(const Matrix<T> &M){
  int rows = M.rows();
  int cols = M.cols();
  for (int i = 0 ; i < rows; ++i){
    for (int j = 0; j < cols; ++j){
      std::cout << M[i][j] << "\t";
    }
    std::cout << std::endl;
  }
}

/**
 * 
 * Realiza el pivoteo para el cálculo de LU
 * 
 * @param[in] A = Matriz de entrada 
 * @param[out] LU = Matriz de salida con LU
 * @param[out] p = Vector de permutaciones 
 * */
template <typename T>
void pivot(const Matrix<T> &A,
           Matrix<T> &LU,
           std::vector<size_t> &permut)
{
  LU = A;
  int rows = A.rows();
  int columns = A.cols();
  if (rows != columns){
    throw anpi::Exception("La matriz debe ser cuadrada.");
  }
  T mayor;
  int index = 0;
  bool setMayor;
  bool inPermutation;
  for (int j = 0; j < columns; ++j)
  {
    inPermutation = false;
    setMayor = false;
    for (int i = 0; i < rows; ++i)
    {

      for (int k = 0; k < j; ++k)
      {
        if ((int)permut[k] == (int)i)
        {
          inPermutation = true;
          break;
        }
      } //end iteración permutación
      if (!inPermutation)
      {
        if (!setMayor)
        {
          index = i;
          mayor = A[i][j];
          setMayor = true;
        }
        else if (std::abs(mayor) < std::abs(A[i][j]))
        {
          index = i;
          setMayor = true;
          mayor = A[i][j];
        }
      }
      else
      {
        inPermutation = false;
      }
      //not in permutation
    } //end rows A
    permut.push_back(index);
    for (int n = 0; n < columns; ++n)
    {
      LU[j][n] = A[index][n];
    } //end for pivoting
  }   //end columns A
}


} //anpi

#endif
