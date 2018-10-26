#include <cmath>
#include <limits>
#include <functional>
#include <iostream>
#include <algorithm>
#include "Exception.hpp"
#include "Matrix.hpp"
#include "SolveLU.hpp"

#ifndef ANPI_SPLINE
#define ANPI_SPLINE

namespace anpi
{


/**
 * Calcula la interpolación cúbica
 * 
 * @param[in] v = Vector de entrada de temperaturas
 * @param[out] x = Vector de resultado con las constantes de las 
 *                  interpolaciones
 * */

template <typename T>
void spline(const std::vector<T> &v, std::vector<T> &x)
{
    int n = v.size();
    std::vector<T> a = v;
    std::vector<T> b, d, u, h, alpha, c, l ,z;
    h = std::vector<T> (n-1, 1); 
    c = l = z = std::vector<T> (n, 0);
    u = d = b = std::vector<T> (n-1,0);
    alpha = std::vector<T> (n-2,0);
    for (int i = 1; i <= n-2; ++i){
        alpha[i] = 3*(a[i+1]-a[i])/h[i] - 3*(a[i]-a[i-1])/h[i-1];
    }
    for (int i = 1; i <= n-2; ++i){
        l[i] = 2*( (i+1) - (i-1) ) - h[i-1]*u[i-1];
        u[i] = h[i]/l[i];
        z[i] = (alpha[i] - h[i-1]*z[i-1])/l[i];
    }

    l[n-1] = 1;
    z[n-1] = 0;
    c[n-1] = 0;

    for (int j = n-2; j >= 0; --j ){
        c[j] = z[j] - u[j]*c[j+1];
        b[j] = (a[j+1] - a[j])/h[j] - h[j]*(c[j+1]+ 2*c[j])/3;
        d[j] = (c[j+1]-c[j])/(3*h[j]);
    }

    for (int i = 0; i < n-1; ++i){
        x.push_back(a[i]);
        x.push_back(b[i]);
        x.push_back(c[i]);
        x.push_back(d[i]);
    }
}

} // namespace anpi

#endif