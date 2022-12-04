
// MIT License

// Copyright (c) 2022 unmanned

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// Rabinovich-fabrikant & Chen Lee: https://marmphco.com/dynamical/writeup.pdf
// Original MATLAB Code for Chua oscillator: http://www.chuacircuits.com/matlabsim.php


#pragma once
#include <math.h>

////////////////////////////////////////////////////////////////////////////////////////
// Roessler ////////////////////////////////////////////////////////////////////////////
typedef struct
{
	float x;
	float y;
	float z;
	float a;
	float b;
	float c;
	float delta;

} roessler_t;

void roessler_init(roessler_t* o);
void roessler_iterate(roessler_t* o);

////////////////////////////////////////////////////////////////////////////////////////
// Hopf ////////////////////////////////////////////////////////////////////////////////
typedef struct
{
	float x;
	float y;
	float p;
    float t;

} hopf_t;

void hopf_init(hopf_t* o);
void hopf_iterate(hopf_t* o);

////////////////////////////////////////////////////////////////////////////////////////
// Helmholz ////////////////////////////////////////////////////////////////////////////
typedef struct
{
	float x;
	float y;
    float z;
	float theta;
	float delta;
    float t;

} helmholz_t;

void helmholz_init(helmholz_t* o);
void helmholz_iterate(helmholz_t* o);

////////////////////////////////////////////////////////////////////////////////////////
// Sprott //////////////////////////////////////////////////////////////////////////////
typedef struct
{
	float x;
	float y;
    float z;
    float t;

} sprott_t;

void sprott_init(sprott_t* o);
void sprott_iterate(sprott_t* o);

////////////////////////////////////////////////////////////////////////////////////////
// Sprott-Linz /////////////////////////////////////////////////////////////////////////
typedef struct
{
	float x;
	float y;
    float z;
    float a;
    float t;

} linz_t;

void linz_init(linz_t* o);
void linz_iterate(linz_t* o);

////////////////////////////////////////////////////////////////////////////////////////
// Sprott-Linz D ///////////////////////////////////////////////////////////////////////
typedef struct
{
	float x;
	float y;
    float z;
   	float a;
   	float t;
   
} linzd_t;

void linzd_init(linzd_t* o);
void linzd_iterate(linzd_t* o);

////////////////////////////////////////////////////////////////////////////////////////
// Sprott 6-term ///////////////////////////////////////////////////////////////////////
typedef struct
{
	float x;
	float y;
    float z;
	float a;
	float b;
    float c;
    float d;
    float t;
    
} sprottst_d;

void sprottst_init(sprottst_d* o);
void sprottst_iterate(sprottst_d* o);

////////////////////////////////////////////////////////////////////////////////////////
// Rayleigh-Benard /////////////////////////////////////////////////////////////////////
typedef struct
{
	double x;
	double y;
    double z;
	float a;
	float r;
    float b;
    float t;
    
} rayleigh_t;

void rayleigh_init(rayleigh_t* o);
void rayleigh_iterate(rayleigh_t* o);

////////////////////////////////////////////////////////////////////////////////////////
// Wang ////////////////////////////////////////////////////////////////////////////////
typedef struct
{
	float x;
	float y;
    float z;
    float w;
	float a;
	float b;
   	float c;
   	float d;
   	float h;
   	float t;
    
} wang_t;

void wang_init(wang_t* o);
void wang_iterate(wang_t* o);

////////////////////////////////////////////////////////////////////////////////////////
// Yu-Wang /////////////////////////////////////////////////////////////////////////////
typedef struct
{
	float x;
	float y;
    float z;
	float a;
	float b;
    float c;
	float d;
    float t;

} yu_wang_t;

void yu_wang_init(yu_wang_t* o);
void yu_wang_iterate(yu_wang_t* o);

////////////////////////////////////////////////////////////////////////////////////////
// Three-Scroll Unified Chaotic System (TSUCS) /////////////////////////////////////////
typedef struct
{
	float x;
	float y;
    float z;
	float a;
	float b;
    float c;
    float d;
    float e;
   	float t;
   
} tsucs_t;

void tsucs_init(tsucs_t* o);
void tsucs_iterate(tsucs_t* o);

////////////////////////////////////////////////////////////////////////////////////////
// Three-Scroll Unified Chaotic System (TSUCS2) ////////////////////////////////////////
typedef struct
{
	float x;
	float y;
    float z;
	float a;
	float b;
    float c;
    float d;
    float e;
    float f;
    float t;
    
} tsucs2_t;

void tsucs2_init(tsucs2_t* o);
void tsucs2_iterate(tsucs2_t* o);

////////////////////////////////////////////////////////////////////////////////////////
// Lorenz //////////////////////////////////////////////////////////////////////////////
typedef struct
{
    float a;
    float b;
    float c;
    float t; 
    float x; 
    float y;
    float z; 
    float delta;

} lorenz_t;

void lorenz_init(lorenz_t* o);
void lorenz_iterate(lorenz_t* o);

////////////////////////////////////////////////////////////////////////////////////////
// Aizawa //////////////////////////////////////////////////////////////////////////////
typedef struct
{
    float a;
    float b;
    float c;
    float d;
    float e;
    float f;
    float t; 
    float x; 
    float y;
    float z; 

} aizawa_t;

void aizawa_init(aizawa_t* o);
void aizawa_iterate(aizawa_t* o);

////////////////////////////////////////////////////////////////////////////////////////
// Halvorsen ///////////////////////////////////////////////////////////////////////////
typedef struct
{
    float a;
    float t; 
    float x; 
    float y;
    float z; 

} halvorsen_t;

void halvorsen_init(halvorsen_t* o);
void halvorsen_iterate(halvorsen_t* o);

////////////////////////////////////////////////////////////////////////////////////////
// Ikeda ///////////////////////////////////////////////////////////////////////////////
typedef struct
{
	float u;
	float x;
	float y;
	float t;

} ikeda_t;

void ikeda_init(ikeda_t* o);
void ikeda_iterate(ikeda_t* o);

////////////////////////////////////////////////////////////////////////////////////////
// Duffing /////////////////////////////////////////////////////////////////////////////
typedef struct
{
	float x;
	float y;
	float a;
	float b;

} duffing_t;

void duffing_init(duffing_t* o);
void duffing_iterate(duffing_t* o);

////////////////////////////////////////////////////////////////////////////////////////
// Henon ///////////////////////////////////////////////////////////////////////////////
typedef struct
{
	float x;
	float y;
    float dy;
    float dx;
	float a;
	float b;
    float t;

} henon_t;

void henon_init(henon_t* o);
void henon_iterate(henon_t* o);

////////////////////////////////////////////////////////////////////////////////////////
// Gingerbreadman //////////////////////////////////////////////////////////////////////
typedef struct
{
	float x;
	float y;
    float t;

} gingerbreadman_t;

void gingerbreadman_init(gingerbreadman_t* o);
void gingerbreadman_iterate(gingerbreadman_t* o);

////////////////////////////////////////////////////////////////////////////////////////
// Van Der Pol /////////////////////////////////////////////////////////////////////////
typedef struct
{
	float x;
	float y;
	float f;
    float t;
	float m;

} vanderpol_t;

void vanderpol_init(vanderpol_t* o);
void vanderpol_iterate(vanderpol_t* o);

////////////////////////////////////////////////////////////////////////////////////////
// Kaplan-Yorke ////////////////////////////////////////////////////////////////////////
typedef struct
{
    long a;
    long b;
    double x;
	double y;
    double t;
	double alpha;

} kaplan_yorke_t;

void kaplan_yorke_init(kaplan_yorke_t* o);
void kaplan_yorke_iterate(kaplan_yorke_t* o);

////////////////////////////////////////////////////////////////////////////////////////
// Rabinovich-Fabrikant ////////////////////////////////////////////////////////////////
typedef struct 
{
    float theta;
    float alpha;
    float x;
    float y;
    float z;
    float t;

} rabinovich_fabrikant_t;

void rabinovich_fabrikant_init(rabinovich_fabrikant_t* o);
void rabinovich_fabrikant_iterate(rabinovich_fabrikant_t* o);

////////////////////////////////////////////////////////////////////////////////////////
// Chen-Lee ////////////////////////////////////////////////////////////////////////////
typedef struct
{
    float a;
    float b;
    float c;
    float x;
    float y;
    float z;
    float t;

} chen_lee_t;

void chen_lee_init(chen_lee_t* o);
void chen_lee_iterate(chen_lee_t* o);

// ////////////////////////////////////////////////////////////////////////////////////////
// // Chua ////////////////////////////////////////////////////////////////////////////////
// struct chua
// {
//     float x = 1, y = 1, z = 1;  // Ins

//     float alpha  =  15.6;
//     float beta   =  28; 
//     float ma     = -1.143;
//     float mb     = -0.714;
//     float h;
//     float t = 0.1f;
    
//     float dx;
//     float dy;
//     float dz;
   
//     void iterate();
// };

// void chua::iterate()
// {
//     h  = ma * x + 0.5f * (ma - mb) * (abs(x + 1.0f) - abs(x - 1.0f));
//     dx = t * (alpha * (y - x - h));
//     dy = t * (x - y + z);
//     dz = t * (- beta * y);

//     x=dx;
//     y=dy;
//     z=dz;
// }

// ////////////////////////////////////////////////////////////////////////////////////////
// // Chua circuit emulation //////////////////////////////////////////////////////////////
// struct realchua
// {
//     float x, y, z;         // In
//     float C1    = 10e-9;   // 10nF
//     float C2    = 100e-9;  // 100nF
//     float R     = 1800;    // 1.8k Ohms
//     float G     = 1/R;

//     // Chua Diode //////////////////////
//     float R1    = 220;
//     float R2    = 220;
//     float R3    = 2200;
//     float R4    = 22000;
//     float R5    = 22000;
//     float R6    = 3300;

//     float Esat  = 9; // 9V phasecvbatteries
//     float E1    = R3/(R2+R3)*Esat;
//     float E2    = R6/(R5+R6)*Esat;

//     float m12   = -1/R6;
//     float m02   = 1/R4;
//     float m01   = 1/R1;
//     float m11   = -1/R3;

//     float m0;
//     float m1    = m12+m11;

//     // Gyrator ////////////////////
//     float R7    = 100;  // 100 Ohms
//     float R8    = 1000; // phasecv 1k Ohms
//     float R9    = 1000; // 1k Ohms
//     float R10   = 1800;
//     float C     = 100*10^(-9);    // 100nF
//     float L     = R7*R9*C*R10/R8; // 18mH 

//     float xdot;
//     float ydot;
//     float zdot;
    
//     void iterate();
// };


// void realchua::iterate()
// {
//     if(E1>E2)
//     m0 = m11 + m02;
//     else
//     m0 = m12 + m01;   

//     float mm1 = m01 + m02;
//     float Emax = std::max(E1, E2);
//     float Emin = std::min(E1, E2);

//     float g;
//     if (abs(x) < Emin) g = x*m1;     
//     else if(abs(x) < Emax )
//     {
//         g = x*m0;
//         if (x > 0) g = g + Emin*(m1-m0);    
//         else g = g + Emin*(m0-m1);  
//     }

//     else if (abs(x) >= Emax)
//     {
//         g = x*mm1;    
//         if (x > 0)
//             g = g + Emax*(m0-mm1) + Emin*(m1-m0);
//         else
//             g = g + Emax*(mm1-m0) +  Emin*(m0-m1);
//     }
//     xdot = (1/C1)*(G*(y-x)-g);
//     ydot = (1/C2)*(G*(x-y)+z);
//     zdot  = -(1/L)*y;
// }

typedef struct 
{
    float zx;
    float zy;
    float cx;
    float cy;
    float t;

} julia_t;


void julia_init(julia_t* o);
void julia_iterate(julia_t* o);