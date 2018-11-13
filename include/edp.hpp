#include <cmath>
#include <limits>
#include <functional>
#include <iostream>
#include <algorithm>

#include "Exception.hpp"
#include "Matrix.hpp"
#include "MapeoTemperatura.hpp"
#include "liebmann.hpp"
#ifndef ANPI_EDP
#define ANPI_EDP

namespace anpi
{


/**
   *Mapea los valores lineales y de interpolación cúbica
   *
   * @param[out] A = Matriz de salida con EDP
   * @param[in] h = tamaño horizontal de la matriz
   * @param[in] v = tamaño vertical de la matriz
   * @param[in] Top = Masa térmica en la parte superior de la placa
   * @param[in] Bot = Masa térmica en la parte inferior de la placa
   * @param[in] Top = Masa térmica en la parte izquierda de la placa
   * @param[in] Top = Masa térmica en la parte derecha de la placa
   * @param[in] aislados = Lados aislados de la placa
   * @param[out] valores = Suma de la tempreraturas 
   */
template <typename T>
void formEDP(anpi::Matrix<T> &A, int h, int v, const std::vector<T> &tTop, const std::vector<T> &tBot,
                        const std::vector<T> &tLeft, const std::vector<T> &tRight, bool aislados[], std::vector<T> &bs)
{ 
    std::vector<T> topTemp, botTemp, leftTemp, rightTemp;
    anpi::mapeo(tTop, v, topTemp);
    anpi::mapeo(tBot, v, botTemp);
    anpi::mapeo(tLeft, h, leftTemp);
    anpi::mapeo(tRight, h, rightTemp);

    size_t lenA = v* h;
    A = anpi::Matrix<T>(lenA, lenA, 0.0);
    std::string solu;
    int k = 0; // filas
    
    k = 0;
    T b = 0;
    size_t bordeDer, bordeInf;
    bordeDer = v-1;
    bordeInf = h-1;
    for (size_t i = 0; i < (size_t)v; ++i)
    {
        b = 0;
        for (size_t j = 0; j <(size_t) h; ++j, ++k)
        {
            A(k, k) = -4;
            if (i == 0)
            {   

                if (!aislados[0]){
                    b-=topTemp[j];
                }
                
            }
            else
            {   
                A(k, k-v) = 1;
            }
            if (i == bordeDer)
            {

                if (!aislados[1]){
                    b-=botTemp[j];
                }
            }
            else
            {       
                A(k, k+v) = 1;
            }
            if (j == 0)
            {

                if (!aislados[2]){
                    b-=leftTemp[j];
                }
            }
            else
            {
                A(k, k-1) = 1;
            }
            if (j == bordeInf)
            {

                if (!aislados[3]){
                    b-=rightTemp[j];
                }
            }
            else
            {
                A(k, k + 1) = 1;            
            }
            bs.push_back(b);
            b=0;
        } //for j
    }     //for i
    
} //end edp

} // namespace anpi
#endif