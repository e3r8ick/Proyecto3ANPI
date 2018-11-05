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

/**
 * Instantiate and test the methods of the Matrix class
 */
BOOST_AUTO_TEST_CASE( getpSize ) {
	BOOST_CHECK (3 == anpi::getpSize(0));
	BOOST_CHECK (4 == anpi::getpSize(1));
	BOOST_CHECK (6 == anpi::getpSize(2));
	BOOST_CHECK (10 == anpi::getpSize(3));
	BOOST_CHECK (18 == anpi::getpSize(4));
	BOOST_CHECK (34 == anpi::getpSize(5));
}


BOOST_AUTO_TEST_SUITE_END()
