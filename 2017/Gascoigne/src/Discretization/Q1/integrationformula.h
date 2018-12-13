/**
*
* Copyright (C) 2004, 2006, 2010 by the Gascoigne 3D authors
*
* This file is part of Gascoigne 3D
*
* Gascoigne 3D is free software: you can redistribute it
* and/or modify it under the terms of the GNU General Public
* License as published by the Free Software Foundation, either
* version 3 of the License, or (at your option) any later
* version.
*
* Gascoigne 3D is distributed in the hope that it will be
* useful, but WITHOUT ANY WARRANTY; without even the implied
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE.  See the GNU General Public License for more
* details.
*
* Please refer to the file LICENSE.TXT for further information
* on this license.
*
**/


#ifndef __integrationformula_h
#define __integrationformula_h

#include  "integrationformulabase.h"

/*------------------------------------------------------------*/

namespace Gascoigne
{
typedef IntegrationFormulaBase<1>  IntegrationFormula1d;
typedef IntegrationFormulaBase<2>  IntegrationFormula2d;
typedef IntegrationFormulaBase<3>  IntegrationFormula3d;

/*------------------------------------------------------------*/

class LineMidPoint : public IntegrationFormula1d{
public:  LineMidPoint();};

/*------------------------------------------------------------*/

class LineTrapez : public IntegrationFormula1d{
public:  LineTrapez();};

/*------------------------------------------------------------*/

class LineSimpson : public IntegrationFormula1d{
public:  LineSimpson();};

/*------------------------------------------------------------*/

class LineGauss1 : public IntegrationFormula1d{
public:  LineGauss1();};

/*------------------------------------------------------------*/

class LineGauss2 : public IntegrationFormula1d{
public:  LineGauss2();};

/*------------------------------------------------------------*/

class LineGauss3 : public IntegrationFormula1d{
public:  LineGauss3();};

/*------------------------------------------------------------*/

class LineGauss4 : public IntegrationFormula1d{
public:  LineGauss4();};

/*------------------------------------------------------------*/

class LineGauss5 : public IntegrationFormula1d{
public:  LineGauss5();};

/*------------------------------------------------------------*/

class LineGauss6 : public IntegrationFormula1d{
public:  LineGauss6();};

/*------------------------------------------------------------*/

class LineGauss7 : public IntegrationFormula1d{
public:  LineGauss7();};

/*------------------------------------------------------------*/

class LineGauss8 : public IntegrationFormula1d{
public:  LineGauss8();};

/*------------------------------------------------------------*/

class LineGauss9 : public IntegrationFormula1d{
public:  LineGauss9();};

/*------------------------------------------------------------*/

class LineGauss10 : public IntegrationFormula1d{
public:  LineGauss10();};

/*------------------------------------------------------------*/

class LineGauss11 : public IntegrationFormula1d{
public:  LineGauss11();};

/*------------------------------------------------------------*/

class LineGauss12 : public IntegrationFormula1d{
public:  LineGauss12();};

/*------------------------------------------------------------*/

class LineGauss13 : public IntegrationFormula1d{
public:  LineGauss13();};

/*------------------------------------------------------------*/

class LineGauss14 : public IntegrationFormula1d{
public:  LineGauss14();};

/*------------------------------------------------------------*/

template<int N, class LineFormula>
class TensorFormula2d : public IntegrationFormula2d
{
public: TensorFormula2d();};

/*------------------------------------------------------------*/

template<int N, class Line>
class TensorFormula3d : public IntegrationFormula3d{
public: TensorFormula3d();};

/*------------------------------------------------------------*/

typedef TensorFormula2d<1,LineMidPoint> QuadMidPoint;
typedef TensorFormula2d<2,LineTrapez>   QuadTrapez;
typedef TensorFormula2d<3,LineSimpson>  QuadSimpson;

typedef TensorFormula2d<1,LineGauss1>   QuadGauss1; 
typedef TensorFormula2d<2,LineGauss2>   QuadGauss4; 
typedef TensorFormula2d<3,LineGauss3>   QuadGauss9;
typedef TensorFormula2d<4,LineGauss4>   QuadGauss16; 
typedef TensorFormula2d<5,LineGauss5>   QuadGauss25; 
typedef TensorFormula2d<6,LineGauss6>   QuadGauss36; 
typedef TensorFormula2d<7,LineGauss7>   QuadGauss49; 
typedef TensorFormula2d<8,LineGauss8>   QuadGauss64; 
typedef TensorFormula2d<9,LineGauss9>   QuadGauss81; 
typedef TensorFormula2d<10,LineGauss10> QuadGauss100; 

/*------------------------------------------------------------*/

typedef TensorFormula3d<2,LineTrapez> HexTrapez;
typedef TensorFormula3d<2,LineGauss2> HexGauss8; 
typedef TensorFormula3d<3,LineGauss3> HexGauss27; 
typedef TensorFormula3d<4,LineGauss4> HexGauss64; 
typedef TensorFormula3d<5,LineGauss5> HexGauss125; 
}

/*------------------------------------------------------------*/

#endif
