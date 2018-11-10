/**
 * Copyright (C) 2017 
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 */


#include <boost/test/unit_test.hpp>

#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>
#include <vector>
/**
 * Unit tests for the matrix class
 */

#include "Matrix.hpp"
#include "Piramidal.hpp"
#include "Allocator.hpp"

// Explicit instantiation of all methods of Matrix


// normal allocator


BOOST_AUTO_TEST_SUITE( Piramidal )

/*
	Create matrix to test final mapping into desired dimensions
*/ 


BOOST_AUTO_TEST_CASE( getpSize ) {
	BOOST_CHECK (3 == anpi::getpSize(0));
	BOOST_CHECK (4 == anpi::getpSize(1));
	BOOST_CHECK (6 == anpi::getpSize(2));
	BOOST_CHECK (10 == anpi::getpSize(3));
	BOOST_CHECK (18 == anpi::getpSize(4));
	BOOST_CHECK (34 == anpi::getpSize(5));
}

BOOST_AUTO_TEST_CASE( fillInitialBorders ) {
	/*
		Create matrix to test averaging out borders 7x5
	*/ 

	anpi::Matrix<float> A = {{20,20,20,20,20,20,20},
													 {80, 0, 0, 0, 0, 0,20},
													 {80, 0, 0, 0, 0, 0,20},
												 	 {80, 0, 0, 0, 0, 0,20},
												 	 { 0,10,20,30,40,50,60}};

// Create matrix to store averages
	anpi::Matrix<float> NewN1 (3, 3);	//create empty 3x3 storage matrix

	anpi::Matrix<float> resA1 = {{20,20,20},	//create expected matrix
															 {80, 0,20},
															 {40,25,40}};


	// Create matrix to store averages
	anpi::Matrix<float> NewN2 (4, 4);	//create empty 4x4 storage matrix

	anpi::Matrix<float> resA2 = {{20,20,20,20},
															 {80, 0, 0,20},
															 {80, 0, 0,20},
															 {40,15,35,40}};


	// Fill all New matrix
	anpi::fillInitialBorders(A, NewN1, 3, 7, 5);
	anpi::fillInitialBorders(A, NewN2, 4, 7, 5);

	// Check matches
	BOOST_CHECK ( NewN1 == resA1 );
	BOOST_CHECK ( NewN2 == resA2 );
}

BOOST_AUTO_TEST_CASE( fillInitialContents ) {

/*
	Create matrix to test filling expanding internal values 
*/ 
	anpi::Matrix<float> B1 = {{ 0, 20, 0},	//create initial matrix
														{80, 35,20},
													  { 0, 20, 0}};

	anpi::Matrix<float> eB1	(4, 4);	//create empty 4x4 storage matrix

	anpi::Matrix<float> expB1 = {{0, 0, 0,0},	//create target matrix 
														   {0,35,35,0},
															 {0,35,35,0},
														   {0, 0, 0,0}};

	anpi::Matrix<float> B2	  = {{ 0, 20, 20, 0},	//create target matrix 
														   {80, 42, 29,20},
															 {80, 44, 28,20},
														   { 0, 20, 20, 0}};

	anpi::Matrix<float> eB2	(6, 6);	//create empty 6x6 storage matrix

	anpi::Matrix<float> expB2 = {{0,  0,  0,  0,  0, 0},	//create second target matrix 
															 {0, 42, 42, 29, 29, 0},
															 {0, 42, 42, 29, 29, 0},
															 {0, 44, 44, 28, 28, 0},
															 {0, 44, 44, 28, 28, 0},
															 {0,  0,  0,  0,  0, 0}};

	// Fill all New matrix
	anpi::fillInitialContents(B1, eB1, 4);
	anpi::fillInitialContents(B2, eB2, 6);

	// Check algorithm
	BOOST_CHECK ( expB1 == eB1 );
	BOOST_CHECK ( expB2 == eB2 );

}

BOOST_AUTO_TEST_CASE( initializeTargetMatrix ) {

	anpi::Matrix<float> L = {{20,20,20,20,20,20,20},
													 {80, 0, 0, 0, 0, 0,20},
													 {80, 0, 0, 0, 0, 0,20},
												 	 {80, 0, 0, 0, 0, 0,20},
												 	 { 0,10,20,30,40,50,60}};

	anpi::Matrix<float> New	  = {{ 0,    20,    20, 0},	//create target matrix 
														   {80, 43.12, 29.84,20},
															 {80, 43.59, 32.11,20},
														   { 0,    20,    20, 0}};

	anpi::Matrix<float> resL ={{20,    20,    20,    20,		20,    20, 20},
														 {80, 43.12, 43.12, 29.84, 29.84, 29.84, 20},
														 {80, 43.59, 43.59, 32.11, 32.11, 32.11, 20},
													 	 {80, 43.59, 43.59, 32.11, 32.11, 32.11, 20},
													 	 { 0,    10,    20,    30,    40,    50, 60}};

	initializeA (L, New, 4, 5, 7);
/*
	std::cout << "L" << std::endl;
		for (int i=0; i<5; i++) {
			for (int j=0; j<7; j++) {
				std::cout << L[i][j] << ", ";
			}
			std::cout << std::endl;
		}
*/
	BOOST_CHECK ( L == resL);
}

BOOST_AUTO_TEST_CASE( piramidalOptimization ) {

	const anpi::Matrix<double> A = {{20,20,20,20,20,20,20},
															 	  {80, 0, 0, 0, 0, 0,20},
																  {80, 0, 0, 0, 0, 0,20},
															 	  {80, 0, 0, 0, 0, 0,20},
															 	  { 0,10,20,30,40,50,60}};

	anpi::Matrix<double> New (5, 7);

	anpi::Matrix<double> resL ={{20,    20,    20,    20,		20,    20, 20},
														 {80, 43.12, 43.12, 29.84, 29.84, 29.84, 20},
														 {80, 43.59, 43.59, 32.11, 32.11, 32.11, 20},
													 	 {80, 43.59, 43.59, 32.11, 32.11, 32.11, 20},
													 	 { 0,    10,    20,    30,    40,    50, 60}};

	anpi::Piramidal<double> (A, New);
	double Eps = 0.5;
	bool expected = true;

	BOOST_CHECK ( expected );
}
			

BOOST_AUTO_TEST_SUITE_END()
