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
   *  Returns size of matrix
   * 
   * @param[in] n: Iteration value
   * @param[out] gpS: size of col and row
   */
int getpSize(int n){
	int gpS = 1;
	
	for (int i = 0; i<n; i++){
		gpS *= 2;
	}

	return gpS +=2; //temperature rows/cols
}

/**
   *  Fills the border of New matrix by averaging out the 
	 * values of the original
	 *
   * @param[in] A: the original matrix
   * @param[out] New: the matrix we are intializing 
   * @param[in] psize: edge size of New square matrix
   * @param[in] Acols: column size of A matrix
   * @param[in] Arows: rows size of A matrix
   */
template <typename T>
void fillInitialBorders(const anpi::Matrix<T> &A,
															anpi::Matrix<T> &New,
												int psize, int Acols, int Arows) 
	{

	int ini, fin = 0;
	T Top, Bot, Izq, Der = 0;
		
	// Iterates on New Matrix to get edge values
	for (int k = 0; k < psize; k++) {
				
		Top = Bot = Izq = Der = 0; //reset avgs	
			
		///Optimizable
		// Get average values of known temperature values Top and Bottom
			ini = (k*Acols/psize)+1;
			fin = ((k+1)*Acols/psize)+1;
			
			for (int i = ini; i<fin; i++) {
				Top += A[0][i];
				Bot += A[Acols-1][i];
			}
			
			New[0][k] = Top/(fin-ini);
			New[psize-1][k] = Bot/(fin-ini);
		// Get average values of known temperature values Left and Right
			ini = (k*Arows/psize)+1;
			fin = ((k+1)*Arows/psize)+1;
			
			for (int i = ini; i<fin; i++) {
				Izq += A[i][0];
				Der += A[i][Arows-1];
			}
			
			New[k][0] = Izq/(fin-ini);
			New[k][psize-1] = Der/(fin-ini);
		}
}

/**
   *  Fills the inner cells of New matrix by copynig and expanding the
	 * values of the old one
	 *
   * @param[in] Old: the matrix we calculated last itration
   * @param[out] New: the matrix we are intializing 
   * @param[in] psize: edge size of New square matrix
   */
template <typename T>
void fillInitialContents(const anpi::Matrix<T> &Old,
															 anpi::Matrix<T> &New,
													 int psize) 
	{
	for (int i = 1; i<psize; i+2) {
		for (int j = 1; j<psize; j+2) {

			// Each old value becomes 4 new ones
			New[i][j] = Old[i][j];
			New[i+1][j] = Old[i][j];
			New[i][j+1] = Old[i][j];
			New[i+1][j+1] = Old[i][j];
		}
	}

}


/**
   *  Creates intial values for the cells in a Matrix based 
	 * on known values of the edges in order to optimize the Liebmann algorithm
   * 
   * @param[in] A the input matrix with desired dimensions and initial values
   * @param[out] L is A with internal values initialized
   */
template <typename T>
void Piramidal(const anpi::Matrix<T> &A, anpi::Matrix<T> &L)
{
	int Arows = A.rows();
	int Acols = A.cols();
	int psize, n = 0;
	anpi::Matrix<T> Old;

	if (Acols < 3 || Arows < 3) {
		anpi::Exception ("Matrix is too small");
	}

//exception matrix col or row is < 3
	
	while (psize<=Arows && psize<=Acols) // stop before matrix size exceeds desired matrx size
		{		
		psize = anpi::getpSize(n);	// get size for new martix 2+2^n
		n++;												// increase iterator

		anpi::Matrix<T> New = Matrix<T>(psize, psize, 0.0); // initialize New matrix in 0s
		//std::cout << New;

			// Inititalize New matrix
		anpi::fillInitialBorders(&A, &New, Arows, Acols, psize);

		if (n-1 != 0) { // if it's 0 -> Old is empty
			anpi::fillInitialContents(&Old, &New, psize);
		}

			// Calculate new internal values
		//anpi::liebmannOnce(const Matrix<T> &A,
		//            Matrix<T> &L, std::vector<T> b);

		Old = New;	// everything that's new becomes old

		psize = anpi::getpSize(n);	// update stop condition
		}

		//createFinalMatrix

	}


} //anpi
