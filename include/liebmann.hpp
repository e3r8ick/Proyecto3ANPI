/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: 
 * @Date  : 
 */

#include <cmath>
#include <limits>
#include <functional>
#include <iostream>
#include <algorithm>
#include <omp.h>
#include <iostream>
#include <fstream>
#include "Exception.hpp"
#include "Matrix.hpp"
#include "SolveLU.hpp"

#ifndef ANPI_LIEBMANN
#define ANPI_LIEBMANN

namespace anpi
{


/**
   *Calcula el método de Liebmann
   *
   * @param[in] A = Matriz de entrada 
   * @param[in] L = Matriz de salida
   * @param[out] b = vector con los resultados
   */
template <typename T>
void liebmann(const Matrix<T> &A,
              Matrix<T> &L, std::vector<T> b)
{
  size_t m1 = (A.rows())/2;
  size_t n1 = (A.cols())/2;
  size_t n2 = n1*2;
  size_t m2 = m1*2;
  T Aiplus1j, Aiminus1j, Aijminus1, Aijplus1;
  L = anpi::Matrix<T>(m2, n2, 0.0);
  anpi::Matrix<T> Lprev = A;
  T eps = std::numeric_limits<T>::epsilon();
  T maxi = std::numeric_limits<T>::max();
  size_t iter = 0;
  size_t i, j, k;
  k = 0;

  time_t time0;   // create timers.
  time_t time1;

  time(&time0);   // get current time.

  while (iter < maxi)
  {
    /**
     * Se hace de esta manera para paralilizar el cálculo, haciendo dos hilos que 
     * corran en paralelo que cada uno calcule la mita de la matriz
     * 
     * */
    # pragma omp parallel \
      shared ( L, b, Aiplus1j,Aiminus1j, Aijplus1, Aijminus1) \
      private ( i, j, k )

    # pragma omp for
    for (i = 0; i < m1; ++i)
    {
      for (j = 0; j < n1; ++j)
      {
        Aiplus1j = Aiminus1j = Aijminus1 = Aijplus1 = 0.0;
        if (i != 0)
        {
          Aiminus1j = L[i - 1][j];
        }
        if (j != 0)
        {
          Aijminus1 = L[i][j - 1];
        }
        if (i != m1)
        {
          Aiplus1j = L[i + 1][j];
        }
        if (j != n1)
        {
          Aijplus1 = L[i][j + 1];
        }
        L[i][j] = (Aiplus1j + Aiminus1j + Aijplus1 + Aijminus1 - b[k]) / 4;
      }
       ++k;
    }
    # pragma omp parallel \
      shared ( L, b, Aiplus1j,Aiminus1j, Aijplus1, Aijminus1) \
      private ( i, j, k )

    # pragma omp for
    for (i = m2; i < m2; ++i)
    {
      for (j = n2; j < n2; ++j)
      {
        Aiplus1j = Aiminus1j = Aijminus1 = Aijplus1 = 0.0;
        if (i != 0)
        {
          Aiminus1j = L[i - 1][j];
        }
        if (j != 0)
        {
          Aijminus1 = L[i][j - 1];
        }
        if (i != m2)
        {
          Aiplus1j = L[i + 1][j];
        }
        if (j != n2)
        {
          Aijplus1 = L[i][j + 1];
        }
        L[i][j] = (Aiplus1j + Aiminus1j + Aijplus1 + Aijminus1 - b[k]) / 4;
      }
       ++k;
    }
    if (abs(Lprev(2, 3) - L(2, 3)) <= eps)
    {
      break;
    }
    Lprev = L;
    ++iter;

  }

  time(&time1);   // get current time after time pass
  double seconds = time1 - time0;
  std::cout << "seconds since start: " << seconds <<'\n';
}

/**
   *Calcula el método de Liebmann una sola vez
   * para la optimización piramidal de valores inciales
   *
   * @param[in] A = Matriz de entrada 
   * @param[in] L = Matriz de salida
   * @param[out] b = vector con los resultados
   */
template <typename T>
void liebmannOnce(const Matrix<T> &A,
              Matrix<T> &L, std::vector<T> b)
{
  size_t m1 = (A.rows())/2;
  size_t n1 = (A.cols())/2;
  size_t n2 = n1*2;
  size_t m2 = m1*2;
  T Aiplus1j, Aiminus1j, Aijminus1, Aijplus1;
  L = anpi::Matrix<T>(m2, n2, 0.0);
  anpi::Matrix<T> Lprev = A;
  T eps = std::numeric_limits<T>::epsilon();
  size_t iter = 0;
  size_t i, j, k;

    /**
     * Se hace de esta manera para paralilizar el cálculo, haciendo dos hilos que 
     * corran en paralelo que cada uno calcule la mita de la matriz
     * 
     * */
    # pragma omp parallel \
      shared ( L, b, Aiplus1j,Aiminus1j, Aijplus1, Aijminus1) \
      private ( i, j )

    # pragma omp for
    for (i = 0; i < m1; ++i)
    {
      for (j = 0; j < n1; ++j)
      {
        Aiplus1j = Aiminus1j = Aijminus1 = Aijplus1 = 0.0;
        if (i != 0)
        {
          Aiminus1j = L[i - 1][j];
        }
        if (j != 0)
        {
          Aijminus1 = L[i][j - 1];
        }
        if (i != m1)
        {
          Aiplus1j = L[i + 1][j];
        }
        if (j != n1)
        {
          Aijplus1 = L[i][j + 1];
        }
        L[i][j] = (Aiplus1j + Aiminus1j + Aijplus1 + Aijminus1) / 4;
      }
    }
    # pragma omp parallel \
      shared ( L, b, Aiplus1j,Aiminus1j, Aijplus1, Aijminus1) \
      private ( i, j )

    # pragma omp for
    for (i = m2; i < m2; ++i)
    {
      for (j = n2; j < n2; ++j)
      {
        Aiplus1j = Aiminus1j = Aijminus1 = Aijplus1 = 0.0;
        if (i != 0)
        {
          Aiminus1j = L[i - 1][j];
        }
        if (j != 0)
        {
          Aijminus1 = L[i][j - 1];
        }
        if (i != m2)
        {
          Aiplus1j = L[i + 1][j];
        }
        if (j != n2)
        {
          Aijplus1 = L[i][j + 1];
        }
        L[i][j] = (Aiplus1j + Aiminus1j + Aijplus1 + Aijminus1) / 4;
      }
    }
  }

} // namespace anpi

#endif
