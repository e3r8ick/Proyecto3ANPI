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
#include "Allocator.hpp"

// Explicit instantiation of all methods of Matrix


// normal allocator


typedef std::complex<double> dcomplex;

typedef std::allocator<float> alloc;

template class anpi::Matrix<dcomplex,alloc>;
template class anpi::Matrix<double  ,alloc>;
template class anpi::Matrix<float   ,alloc>;
template class anpi::Matrix<int     ,alloc>;

typedef anpi::Matrix<dcomplex,alloc> cmatrix;
typedef anpi::Matrix<double  ,alloc> dmatrix;
typedef anpi::Matrix<float   ,alloc> fmatrix;
typedef anpi::Matrix<int     ,alloc> imatrix;

// aligned allocator
typedef anpi::aligned_allocator<float> aalloc;

template class anpi::Matrix<dcomplex,aalloc>;
template class anpi::Matrix<double  ,aalloc>;
template class anpi::Matrix<float   ,aalloc>;
template class anpi::Matrix<int     ,aalloc>;

typedef anpi::Matrix<dcomplex,aalloc> acmatrix;
typedef anpi::Matrix<double  ,aalloc> admatrix;
typedef anpi::Matrix<float   ,aalloc> afmatrix;
typedef anpi::Matrix<int     ,aalloc> aimatrix;

// row aligned allocator
typedef anpi::aligned_row_allocator<float> aralloc;

template class anpi::Matrix<dcomplex,aralloc>;
template class anpi::Matrix<double  ,aralloc>;
template class anpi::Matrix<float   ,aralloc>;
template class anpi::Matrix<int     ,aralloc>;

typedef anpi::Matrix<dcomplex,aralloc> arcmatrix;
typedef anpi::Matrix<double  ,aralloc> ardmatrix;
typedef anpi::Matrix<float   ,aralloc> arfmatrix;
typedef anpi::Matrix<int     ,aralloc> arimatrix;

#if 1
# define dispatchTest(func) \
  func<cmatrix>();          \
  func<dmatrix>();          \
  func<fmatrix>();          \
  func<imatrix>();          \
                            \
  func<acmatrix>();         \
  func<admatrix>();         \
  func<afmatrix>();         \
  func<aimatrix>();         \
                            \
  func<arcmatrix>();        \
  func<ardmatrix>();        \
  func<arfmatrix>();        \
  func<arimatrix>();

#else
# define dispatchTest(func) func<arfmatrix>(); 
#endif




BOOST_AUTO_TEST_SUITE( Matrix )

template<class M>
void testConstructors() {
  // Constructors
  { // default
    M a;
    BOOST_CHECK( a.rows() == 0);
    BOOST_CHECK( a.cols() == 0);
    BOOST_CHECK( a.dcols() == 0);
  }
  { // unitilialized
    M a(2,3,anpi::DoNotInitialize);
    BOOST_CHECK( a.rows() == 2);
    BOOST_CHECK( a.cols() == 3);
    BOOST_CHECK( a.dcols() >= 3);
  }
  { // default initialized
    M a(3,2);
    BOOST_CHECK( a.rows() == 3);
    BOOST_CHECK( a.cols() == 2);
    BOOST_CHECK( a(0,0) == typename M::value_type(0));
  }
  { // default initialized
    M a(3,2,typename M::value_type(4));
    BOOST_CHECK( a.rows() == 3);
    BOOST_CHECK( a.cols() == 2);
    BOOST_CHECK( a(0,0) == typename M::value_type(4));
  }
  { // initializer_list
    M a = { {1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15} };
    BOOST_CHECK( a.rows() == 3);
    BOOST_CHECK( a.cols() == 5);
    
    BOOST_CHECK( a(0,0) == typename M::value_type(1));
    BOOST_CHECK( a(1,2) == typename M::value_type(8));
    BOOST_CHECK( a(2,3) == typename M::value_type(14));
  }
  { // Copy constructor
    M a = { {1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15} };
    M b(a);

    BOOST_CHECK( a==b );
    BOOST_CHECK( b.rows() == 3 );
    BOOST_CHECK( b.cols() == 5 );
    BOOST_CHECK( b.data() != a.data());
  }

  { // Move constructor
    M a = { {1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15} };
    M b(std::move(a));

    BOOST_CHECK( b.rows() == 3 );
    BOOST_CHECK( b.cols() == 5 );

    BOOST_CHECK( a.empty() );
  }
  { // Mem constructor
    M a = { {1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15} };
    M b(a.rows(),a.cols(),a.data());

    BOOST_CHECK( a==b );
    BOOST_CHECK( b.rows() == 3 );
    BOOST_CHECK( b.cols() == 5 );
    BOOST_CHECK( b.data() != a.data() );
  }
}
/**
 * Instantiate and test the methods of the Matrix class
 */
BOOST_AUTO_TEST_CASE( Constructors ) {
  dispatchTest(testConstructors);
}

template<class M>
void testComparison() {
  // == and !=
  M a = { {1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15} };
  M b = { {1,2,3,4,5},{6,7,9,9,10},{11,12,13,14,15} };

  BOOST_CHECK( (a!=b) );
  
  b(1,2)=typename M::value_type(8);
  
  BOOST_CHECK( (a==b) );
}  

BOOST_AUTO_TEST_CASE(Comparison) 
{
  dispatchTest(testComparison);
}

template<class M>
void testAssignment() {
  { // Move assignment
    M a = { {1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15} };
    M c(a);
    M b;
    b=std::move(a);
    BOOST_CHECK(a.empty() );
    BOOST_CHECK(!b.empty() );
    BOOST_CHECK(b.rows()==3 );
    BOOST_CHECK(b.cols()==5 );
    BOOST_CHECK(b==c );
  }
  { // assignment
    M a = { {1,2,3,4,5},{5,6,7,8,9},{9,10,11,12,13} };
    M b;
    b=a;
    BOOST_CHECK(a==b );
  }  
  { // swap
    M a = { {1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15} };
    M b = { {13,14},{15,16} };

    M c(a);
    M d(b);

    BOOST_CHECK( a==c );
    BOOST_CHECK( d==b );
    
    c.swap(d);
    BOOST_CHECK( a==d );
    BOOST_CHECK( b==c );
  }
  { // column
    M a = { {1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15} };
    std::vector<typename M::value_type> col = a.column(1);
    std::vector<typename M::value_type> ref = {2,7,12};
    BOOST_CHECK( col == ref );
  }
}

BOOST_AUTO_TEST_CASE(Assignment)
{
  dispatchTest(testAssignment);
}

template<class M>
void testArithmetic() {
  
  {
    M a = { {1,2,3},{ 4, 5, 6} };
    M b = { {7,8,9},{10,11,12} };
    M r = { {8,10,12},{14,16,18} };
    
    M c(a);
    c+=b;
    BOOST_CHECK(c==r );
    c=a+b;
    BOOST_CHECK(c==r );


    c=M{ {1,2,3},{ 4, 5, 6} } + b;
    BOOST_CHECK(c==r );

    c=a+M{ {7,8,9},{10,11,12} };
    BOOST_CHECK(c==r );
  }

  {
    M a = { {1,2,3},{ 4, 5, 6} };
    M b = { {7,8,9},{10,11,12} };
    M r = { {-6,-6,-6},{-6,-6,-6} };
    
    M c(a);
    c-=b;
    BOOST_CHECK( c==r );
    c=a-b;
    BOOST_CHECK( c==r );


    c=M{ {1,2,3},{ 4, 5, 6} } - b;
    BOOST_CHECK( c==r );

    c=a-M{ {7,8,9},{10,11,12} };
    BOOST_CHECK( c==r );
  } 
}

BOOST_AUTO_TEST_CASE(Arithmetic) {
  dispatchTest(testArithmetic);  
}



template <class M, typename T>

void mxmtest(){
  T eps = std::numeric_limits<T>::epsilon();
  M a = {{1, 2, 3}, 
        {-1, 0, 3}}; //2x3;
  M b = {{2, 3, 4, 5}, 
         {1, 3, 0 ,1}}; //2x4;

  //matrix * matrix dimensions
  try{
    M c = a*b;
    BOOST_CHECK_MESSAGE(false, "Not propertly checking dimensions");
  }catch(anpi::Exception& exc){
    BOOST_CHECK_MESSAGE(true,"Properly checking dimensions");
  }

  //check multi by zero
  M zero = {{0 , 0, 0},{0, 0, 0}, {0,0,0}};
  a = {{1, 2, 3}, {2, 3, 4}, {3, 4, 5}};
  M p = a * zero;
  for (int i = 0; i < (int) a.rows(); ++i){
    for (int j = 0; j < (int) a.cols(); ++j){
      BOOST_CHECK(abs(p[i][j])<eps);
    }
  }

  //check multi by identity
  M ident = {{1,0,0},{0,1,0},{0,0,1}};
  p = a * ident;
  for (int i = 0; i < 3; ++i){
    for (int j = 0; j < 3; ++j){
      BOOST_CHECK(abs(p[i][j] - a[i][j]) < eps);
    }
  }

  //checks multi by inverse
  a = {{1,2,3},{0,1,4}, {5,6,0}};
  M inv = {{-24,18,5},{20,-15,-4},{-5,4,1}};
  p = a * inv;
  for (int i = 0; i < (int) a.rows(); ++i){
    for (int j = 0; j < (int) a.cols(); ++j){
      if (i == j){
        BOOST_CHECK(abs(p[i][j] - 1 ) < eps);
      }else{
        BOOST_CHECK(abs(p[i][j]) < eps);
      }
    }
  }
}

template <class M, typename T>
void mxvtest(){

  //checks matrix * vector dimensions

  M a = {{1, 2, 3},
         {-1, 0, 3}}; //2x3;
  std::vector<T> b = {2, 3, 4, 5};//4x1
  try{
    std::vector<T> c = a*b;
    BOOST_CHECK_MESSAGE(false, "Not propertly checking the dimensions");
  }catch(anpi::Exception& exc){
    BOOST_CHECK_MESSAGE(true,"Properly checking dimensions");
  }

  //checks vector * matrix dimensions
  try{
    std::vector<T> c = b*a;
    BOOST_CHECK_MESSAGE(false, "Not propertly checking the dimensions");
  }catch(anpi::Exception& exc){
    BOOST_CHECK_MESSAGE(true,"Properly checking dimensions");
  }



}

BOOST_AUTO_TEST_CASE(ProductMxM){
  mxmtest<fmatrix, float>();
  mxmtest<dmatrix, double>();
}

BOOST_AUTO_TEST_CASE(ProductMxV){
  mxvtest<fmatrix, float>();
  mxvtest<dmatrix, double>();
}
  
BOOST_AUTO_TEST_SUITE_END()
