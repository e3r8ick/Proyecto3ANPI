/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: 
 * @Date  : 24.02.2018
 */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <boost/program_options.hpp>
#include <string>
#include "Matrix.hpp"
#include "edp.hpp"
#include "liebmann.hpp"
#include "SolveLU.hpp"


namespace po = boost::program_options;
std::string findEdge(std::string line)
{
  if (line == "top")
  {
    return "t";
  }
  else if (line == "bottom")
  {
    return "b";
  }
  else if (line == "left")
  {
    return "l";
  }
  else if (line == "right")
  {
    return "d";
  }
  else
  {
    return "";
  }
}

int main(int ac, char *av[])
{
  /**************DECLARA E INICIALIZA LAS VARIABLES *****************/
  std::vector<double> tempsTop, tempBot, tempLeft, tempRight;
  bool visualize, flow;
  bool i[4] = {true, true, true, true}; //aislar tblr;
  std::string fileName;
  double value;
  int hori, vert, grid;
  hori = vert = 20;
  grid = 5;
  flow = false;
  visualize = true;

  std::string line, borde;
  borde == "";
  std::vector<double> values;
  bool notEOF = true;
  std::ifstream inFile;

  /***************DEFINE LAS OPCIONES *********************************/
  po::options_description desc("Optiones");
  desc.add_options()("help", "produce help message")
                    ("t,t", po::value<double>(), "Indica temperatura en borde superior")
                    ("b,b", po::value<double>(), "Indica temperatura en borde inferior")
                    ("l,l", po::value<double>(), "Indica temperatura en borde izquierdo")
                    ("d,d", po::value<double>(), "Indica temperatura en borde izquierdo")
                    ("i,i", po::value<std::string>(), "Borde a aislar")
                    ("p,p", po::value<std::string>(), "Nombre del archivo con el perfil término")
                    ("h,h", po::value<int>(), "Pixeles horizontales")
                    ("v,v", po::value<int>(), "Pixeles verticales")
                    ("q,q", "Desactivar visualización")
                    ("f,f", "Activa el cálculo de flujo de calor")
                    ("g,g", po::value<int>(), "Tamaño de rejilla");
  po::variables_map map;
  try
  {
    po::store(po::parse_command_line(ac, av, desc), map);
    po::notify(map);
  }
  catch (std::exception e)
  {
    std::cerr << "Opción inválida o repetida" << std::endl;
    return -1;
  }

  /******************COMIENZA A INTERPRETAR LOS COMANDOS *************************/
  if (map.count("help"))
  {
    std::cout << desc << "\n";
    return 0;
  }
  if (map.count("h"))
  {
    hori = map["h"].as<int>();
  }
  if (map.count("v"))
  {
    vert = map["v"].as<int>();
  }
  if (map.count("g"))
  {
    grid = map["g"].as<int>();
  }
  if (map.count("p"))
  {
    fileName = map["p"].as<std::string>();
  }
  if (map.count("q"))
  {
    visualize = false;
  }
  if (map.count("f"))
  {
    flow = true;
  }

  /************LEE EL ARCHIVO**********************/
  inFile.open("../" + fileName);
  if (!inFile)
  {
    std::cerr << "No se puede abrir el archivo" << std::endl;
    return -1;
  }

  /****ENCUENTRA EL PRIMER BORDE *****/
  if (inFile >> line)
  {
    borde = findEdge(line);
    if (borde == "")
    {
      std::cerr << "Debe indicar un borde válido" << std::endl;
      return -1;
    }
  }
  else
  {
    std::cerr << "El archivo está vacío" << std::endl;
    return -1;
  }

  /****ENCUENTRA LOS VALORES DE TEMPERATURA PARA TODOS LOS BORDES */
  do
  {
    if (!(inFile >> line))
    {
      notEOF = false;
      line = "EOF";
    }
    try
    {
      value = std::stod(line);
      if (borde == "")
      {
        std::cerr << "Debe indicar un borde válido" << std::endl;
        return -1;
      }
      values.push_back(value);
    }
    catch (std::exception e)
    {
      if (values.size() == 0)
      { //no se indico temp para el borde (error)
        std::cerr << "Debe indicar al menos 1 valor de temperatura para el borde" << std::endl;
        return -1;
      }
      /*setear los valores a los bordes */
      if (borde == "t")
      {
        tempsTop = values;
        i[0] = false;
      }
      else if (borde == "b")
      {
        tempBot = values;
        i[1] = false;
      }
      else if (borde == "l")
      {
        tempLeft = values;
        i[2] = false;
      }
      else if (borde == "d")
      {
        tempRight = values;
        i[3] = false;
      }
      else
      {
        std::cerr << "ocurrio un error inesperado" << std::endl;
        return -1;
      }
      values.clear();
      if (!notEOF)
      {
        break;
      }
      borde = findEdge(line);
    }               //end catch
  } while (notEOF); //end while
  inFile.close();

  /***************TERMINA DE INTERPRETAR LOS COMANDOS QUE DEPENDEN DE TEMPERATURA***********/

  if (map.count("t"))
  {
    tempsTop.clear();
    tempsTop.push_back(map["t"].as<double>());
    i[0] = false;
  }
  if (map.count("b"))
  {

    tempBot.clear();
    tempBot.push_back(map["b"].as<double>());
    i[1] = false;
  }
  if (map.count("l"))
  {
    tempLeft.clear();
    tempLeft.push_back(map["l"].as<double>());
    i[2] = false;
  }
  if (map.count("d"))
  {

    tempRight.clear();
    tempRight.push_back(map["d"].as<double>());
    i[3] = false;
  }
  if (map.count("i"))
  {
    std::string aislados = map["i"].as<std::string>();
    size_t len = aislados.length();
    if (len == 0)
    {
      std::cerr << "Error: si indica la option -i debe indicar al menos 1 borde para aislar" << std::endl;
    }
    else
    {
      for (size_t j = 0; j < len; ++j)
      {
        if (aislados[j] != 't' && aislados[j] != 'b' && aislados[j] != 'l' && aislados[j] != 'd')
        {
          std::cerr << "Error: Borde invalido para aislar, debe ser 't', 'b', 'l' o 'd'" << std::endl;
          return -1;
        }
      } // end for
      if (aislados.find('t') != std::string::npos)
      {
        i[0] = true;
      }
      if (aislados.find('b') != std::string::npos)
      {
        i[1] = true;
      }
      if (aislados.find('l') != std::string::npos)
      {
        i[2] = true;
      }
      if (aislados.find('d') != std::string::npos)
      {
        i[3] = true;
      }
    } // end if else len == 0
  }   // end maps i

  /********************HACE LAS LLAMADAS**************************/
  anpi::Matrix<double> edp;
  std::vector<double> sol;
  std::vector<double> x;
  anpi::formEDP(edp, hori, vert, tempsTop, tempBot, tempLeft, tempRight, i, sol);
  anpi::Matrix<double> L;
  anpi::liebmann(edp, L, sol);
  anpi::solveLU(edp, x, sol);

  /* guarda la matriz en un archivo */
  std::ofstream matrixFile;  
  int count = 0;
  matrixFile.open("../Matrix.txt");
  for (size_t i = 0; i < x.size(); ++i){
    if (count == hori){
      count = 0;
      matrixFile << "\n";
    }
     matrixFile << x[i] << "\t";
     count++;
  }
  matrixFile.close();

  matrixFile.open("../flags.txt");
  if (visualize){
    matrixFile << "True" << "\t";
  }else{
    matrixFile << "False" << "\t";
  }
  if (flow){
    matrixFile << "True" << "\t";
  }else{
    matrixFile << "False" << "\t";
  }
  matrixFile << grid << "\t";
  matrixFile << "\n";
  matrixFile.close();
  return 0;
}