/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Josafat Vargas
 * @Date  : 28.10.2018
 */

#include <cmath>
#include <limits>
#include <functional>
#include <iostream>
#include <algorithm>
#include <vector>

#include "Exception.hpp"
#include "Matrix.hpp"
#include "liebmann.hpp"

namespace anpi
{

/**
   *  Expands a known matrix from a 3x3 to it's desired dimensions
   * 
   * @param[in] A the input matrix with desired dimensions and initial values
   * @param[out] expanded matrix
   */
template <typename T>
bool expandMatrix(const anpi::Matrix<T> &A,
						anpi::Matrix<T> &O)
{
	
	int Arows = A.rows();
	int Acols = A.cols();
	anpi::Matrix<T> Old;
	
	for (int n = 0; (2^n+2<Arows || 2^n+2<Acols); n++)
	{
		//Set initial internal values
		int size = 2^n+2;
		if (n==0) {
			New = Matrix<T>(size, size, 0.0); // initialize 3x3 Matrix in 0s
		}
		else {
			New = Matrix<T>(size, size);
			for (i = 1; i<size; i+2) {
				for (j = 1; j<size; j+2) {
					// Each old value becomes 4 new ones
					New[i][j] = Old[i][j];
					New[i+1][j] = Old[i][j];
					New[i][j+1] = Old[i][j];
					New[i+1][j+1] = Old[i][j];
				}
			}
		}
			
		int ini, fin = 0;
		T Top, Bot, Izq, Der = 0;
		
		// Iterates on New Matrix to get edge values
		for (int k = 0; k < size; k++) {
				
			Top = Bot = Izq = Der = 0; //reset avgs	
			
		///Optimizable
		// Get average values of known temperature values Top and Bottom
			ini = (k*Acols/size)+1;
			fin = ((k+1)*Acols/size)+1;
			
			for (int i = ini; i<fin; i++) {
				Top += A[0][i];
				Bot += A[cols][i];
			}
			
			New[0][k] = Top/(fin-ini);
			New[size-1][k] = Bot/(fin-ini);
		// Get average values of known temperature values Left and Right
			ini = (k*Arows/size)+1;
			fin = ((k+1)*Arows/size)+1;
			
			for (int i = ini; i<fin; i++) {
				Izq += A[i][0];
				Der += A[i][rows];
			}
			
			New[k][0] = Izq/(fin-ini);
			New[k][size-1] = Der/(fin-ini);
			
		}
		
		//anpi::liebmannOnce(New, New, V); // calcula los nuevos valores internos
		Old = New;	
	}
	
	// Expand final marix to match desired dimmensions
	
	// Edges
	for (int i = 0; i<Arows; i++) {
		O[0][i] = A[0][i];
		O[Arows][i] = A[Arows][i];
	}
	for (int j = 0; j<Acols; j++) {
		O[0][j] = A[0][j];
		O[Acols][i] = A[Acols][i];
	}
	
	//Set internal value
	
	//O[i][j] = Old[f(i)][f(j)]
	

	
} //anpi
