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
anpi::Matrix resB1 = {{ 0,    20,    20 , 0},		// B2 = resB1
		      						{80,42.5  ,29.375 ,20},
		      						{80,44.375,28.4375,20},
		      						{ 0,    20,     20, 0}};

anpi::Matrix expB2 = {{ 0,     20,     20,      20,      20,  0},
		  								{80,   42.5,   42.5,  29.375,  29.375, 20},
		  								{80,   42.5,   42.5,  29.375,  29.375, 20},
				  						{80, 44.375, 44.375, 28.4375, 28.4375, 20},
				  						{80, 44.375, 44.375, 28.4375, 28.4375, 20},
				  						{ 0,     20,     20,      20,      20,  0}};
*/

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
		Create matrix to test averaging out borders
	*/ 
	anpi::Matrix<float> A = {{ 0,20,20,20,20, 0},
													 {80, 0, 0, 0, 0,20},
													 {80, 0, 0, 0, 0,20},
												 	 {80, 0, 0, 0, 0,20},
												 	 {80, 0, 0, 0, 0,20},
												 	 { 0,20,20,20,20, 0}};


// Create matrix to store averages
	anpi::Matrix<float> NewN1 (3, 3);	//create empty 3x3 storage matrix

	anpi::Matrix<float> resA1 = {{ 0,20, 0},	//create expected matrix
															 {80, 0,20},
															 { 0,20, 0}};


	// Create matrix to store averages
	anpi::Matrix<float> NewN2 (4, 4);	//create empty 4x4 storage matrix

	anpi::Matrix<float> resA2 = {{ 0,20,20, 0},
															 {80, 0, 0,20},
															 {80, 0, 0,20},
															 { 0,20,20, 0}};

	// Create matrix to store averages
	anpi::Matrix<float> NewN3 (6, 6);	//create empty 6x6 matrix

	// Fill all New matrix
	anpi::fillInitialBorders(A, NewN1, 3, 6, 6);
	//anpi::fillInitialBorders(A, NewN2, 4, 6, 6);
	//anpi::fillInitialBorders(A, NewN3, 6, 6, 6);

	// Check matches
	BOOST_CHECK ( NewN1 == resA1 );
	//BOOST_CHECK ( NewN2 == resA2 );
	//BOOST_CHECK ( NewN3 == A);
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
/*
	// Fill all New matrix
	anpi::fillInitialContents(B1, eB1, 4);
	anpi::fillInitialContents(B2, eB2, 6);

	// Check algorithm
	BOOST_CHECK ( expB1 == eB1 );
	BOOST_CHECK ( expB2 == eB2 );
*/
}


BOOST_AUTO_TEST_SUITE_END()
