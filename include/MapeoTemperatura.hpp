/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: 
 * @Date  : 03.03.2018
 */

#include <cmath>
#include <limits>
#include <functional>
#include <iostream>
#include <algorithm>

#include "Exception.hpp"
#include "Matrix.hpp"
#include "spline.hpp"

#ifndef ANPI_MAPEO
#define ANPI_MAPEO

namespace anpi
{

/**
   *Mapea los valores lineales y de interpolación cúbica
   *
   * @param[in] val = val; valores de entrada 
   * @param[in] tamano = tamaño del vector;
   * @param[out] valores = valores de salida mapeados de acuerdo al tamaño
   */
template <typename T>
void mapeo(const std::vector<T>& val, const size_t tamano, std::vector<T>& valores)
{
    //asignacion de valores iniciales
    valores = std::vector<T>(0,tamano);

    size_t Tval = val.size();
    if(Tval == 1){
        for(size_t i = 0; i < tamano; ++i){
            valores.push_back(val[0]);
        }
    }
    else if (Tval == 2){
        T valPrev;
        for(size_t i = 0; i < tamano; ++i){
            T tmp = (val[1]-val[0])/(tamano-1);
            valPrev = val[0] +  tmp*i;
            valores.push_back(valPrev);

        }
    }
    else{
        std::vector <T> sol, x;
        anpi::spline(val, x);
        std::vector<T> va, vb, vc, vd;
        for (size_t i = 0; i < x.size(); i+=4){
            va.push_back(x[i]);
            vb.push_back(x[i+1]);
            vc.push_back(x[i+2]);
            vd.push_back(x[i+3]);
        }
        T iter = (T)  (val.size()-1) / (tamano-1);
        size_t index = 0;
        T i = 0;
        T counts = 1;

        while (counts <= tamano){
            if (index == va.size()){
                i = 1;
                index--;
                valores.push_back((vd[index] * i*i*i) + (vc[index] * i *i) + (vb[index] * i) + va[index]);
                break;

            }
            valores.push_back((vd[index] * i*i*i) + (vc[index] * i *i) + (vb[index] * i) + va[index]);
            i+=iter;
            if (i >= 1){
                index++;
                i -= 1;
            }
            counts++;            
        }
   }
    
}

} // namespace anpi

#endif