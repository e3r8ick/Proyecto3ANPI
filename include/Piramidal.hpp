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

		//Optimizable -> separate threads for cols and rows
		//Optimizable -> all at once if A is square 
		//Optimizable -> don't consider corners

	int k = 0; // A col iterator
	
	// average out top and bottom line
	for (int i = 0; i<psize; i++) {
		int iupper = (int) Acols*(i+1)/psize; 
		int usedValues = 0; //divider
		float top, bot = 0; //accumulators
																						//std::cout << "i=" << i << ": ";
		for (; k<iupper; k++) {
			usedValues++;
			top += A[0][k];
			bot += A[Arows-1][k];
																						//std::cout << A[Arows-1][k] << ", ";
		}

																						//std::cout << "acum: " << bot << " | ";
																						//std::cout << "used: " << usedValues << " | ";

		//get averages and avoid div by 0					
		if (usedValues != 0) {
				top /= usedValues;
				bot /= usedValues;
		}	else {
				top /= 1;
				bot /= 1;
		}
																						//std::cout << " average: " << bot << std::endl;

		//store averages
		New[0][i] = top;
		New[psize-1][i] = bot;

		// hard-reset acumulators
		top -= top;
		bot -= bot;
	}

	int m = 0; // A row iterator
	
	// average out left and right lines
	for (int j = 0; j<psize; j++) {
		int jupper = (int) Arows*(j+1)/psize; 
		int usedValues = 0; //divider
		int izq, der = 0; //accumulators
		for (; m<jupper; m++) {
			usedValues++;
			izq += A[m][0];
			der += A[m][Acols-1];
		}
		//get averages
		if (usedValues != 0) {
			izq /= usedValues;
			der /= usedValues;
		}	else {
				izq /= 1;
				der /= 1;
		}

		//store averages
		New[j][0] = izq;
		New[j][psize-1] = der;

		// hard-reset acumulators
		izq -= izq;
		der -= der;
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

	// Iterator for Old matrix
	int l, k = 1;

	for (int i = 1; i<=psize-2; i+=2) {
		l = 1;
		for (int j = 1; j<=psize-2; j+=2) {
			// Each old value becomes 4 new ones
			New[i][j] = Old[k][l];
			New[i+1][j] = Old[k][l];
			New[i][j+1] = Old[k][l];
			New[i+1][j+1] = Old[k][l];
			l++;
		}
	k++;
	}

}

template <typename T>
void initializeA (anpi::Matrix<T> &L,
								  anpi::Matrix<T> &New,
							int psize, int Acols, int Arows) {

	//ignore borders
	psize -= 2;
	Acols -= 2;
	Arows -= 2;
	
	// i, j are coordinates from New
	// aI,aJ are the target coordinate from A
	for (int j = 0; j<psize; j++){
		int Ari = Arows*j/psize+1;
		int Arf = Arows*(j+1)/psize+1;

		for (int aJ = Ari; aJ < Arf; aJ++) {
																						//std::cout << "j=" << j << ": ";
			for (int i = 0; i<psize; i++){
				int Aci = Acols*i/psize+1;
				int Acf = Acols*(i+1)/psize+1;
																						//std::cout << "i=" << i << ": ";
				for (int aI = Aci; aI < Acf; aI++) {
					L[aI][aJ] = New[i+1][j+1];
																						//std::cout << aI << ", ";
				}
																						//std::cout << "| " ;		
			}
																						//std::cout << std::endl;	
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
	int psize, n, lsize = 0;
	anpi::Matrix<T> Old;

	if (Acols < 3 || Arows < 3) {
		anpi::Exception ("Matrix is too small");
	}

//exception matrix col or row is < 3
	
	while (psize<=Arows && psize<=Acols) // stop before matrix size exceeds desired matrx size
		{		
		psize = anpi::getpSize(n);	// get size for new matix 2+2^n
		lsize = anpi::getpSize(n);	// get size last size actually used
		n++;												// increase iterator

		anpi::Matrix<T> New = Matrix<T>(psize, psize, 0.0); // initialize New matrix in 0s

			// Inititalize New matrix
		anpi::fillInitialBorders(&A, &New, Arows, Acols, psize);

		if (n-1 != 0) { // if it's 0 -> Old is empty
			anpi::fillInitialContents(&Old, &New, psize);
		}

			// Calculate new internal values using onw iteratio of the liebmann algorithm
		//std::vector<T> b
		//anpi::liebmannOnce(A, New, b);

		Old = New;	// everything that's new becomes old

		psize = anpi::getpSize(n);	// update stop condition
		}

		initializeA (L, Old, lsize, Acols, Arows); 

}


} //anpi
