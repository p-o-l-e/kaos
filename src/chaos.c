
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


#include "chaos.h"

////////////////////////////////////////////////////////////////////////////////////////
// Roessler ////////////////////////////////////////////////////////////////////////////
void roessler_init(roessler_t* o)
{
    o->x = 1.0f;
	o->y = 1.0f;
	o->z = 1.0f;
	o->a = 0.2f;
	o->b = 0.2f;
	o->c = 5.7f;
	o->delta = 0.01f;
}

void roessler_iterate(roessler_t* o)
{
	o->x += (-o->y - o->z) * o->delta;
	o->y += (o->x + o->a * o->y) * o->delta;
	o->z += (o->b + o->z * (o->x - o->c)) * o->delta;
}

////////////////////////////////////////////////////////////////////////////////////////
// Hopf ////////////////////////////////////////////////////////////////////////////////
void hopf_init(hopf_t* o)
{
    o->x = 0.01f;
	o->y = 0.01f;
	o->p = 0.11f;
    o->t = 0.01f;
}

void hopf_iterate(hopf_t* o)
{
	o->x += o->t * ( -o->y + o->x * (o->p - (o->x*o->x + o->y*o->y)));
	o->y += o->t * (  o->x + o->y * (o->p - (o->x*o->x + o->y*o->y)));
}

////////////////////////////////////////////////////////////////////////////////////////
// Helmholz ////////////////////////////////////////////////////////////////////////////
void helmholz_init(helmholz_t* o)
{
	o->x = 0.1f;
	o->y = 0.1f;
    o->z = 0.1f;
	o->theta = 5.11f;
	o->delta = 0.55f;
    o->t = 0.01f;
}

void helmholz_iterate(helmholz_t* o)
{
	o->x += o->t * o->y;
	o->y += o->t * o->theta * o->z;
	o->z += o->t * ( -o->z - o->delta * o->y - o->x - o->x * o->x );
}


////////////////////////////////////////////////////////////////////////////////////////
// Sprott //////////////////////////////////////////////////////////////////////////////
void sprott_init(sprott_t* o)
{
	o->x = 0.1f;
	o->y = 0.1f;
    o->z = 0.1f;
    o->t = 0.1f;
}

void sprott_iterate(sprott_t* o)
{
	o->x += o->t * o->y;
	o->y += o->t * (o->y * o->z - o->x);
    o->z += o->t * (1.0f - o->y * o->y);
}

////////////////////////////////////////////////////////////////////////////////////////
// Sprott-Linz /////////////////////////////////////////////////////////////////////////
void linz_init(linz_t* o)
{
	o->x = 0.1f;
	o->y = 0.1f;
    o->z = 0.1f;
    o->a = 0.5f;
    o->t = 0.1f;
}

void linz_iterate(linz_t* o)
{
	o->x += o->t * (o->y + o->z);
	o->y += o->t * (o->y * o->a - o->x);
    o->z += o->t * (o->x * o->x - o->z);
}

////////////////////////////////////////////////////////////////////////////////////////
// Sprott-Linz D ///////////////////////////////////////////////////////////////////////
void linzd_init(linzd_t* o)
{
	o->x = 0.1f;
	o->y = 0.1f;
    o->z = 0.1f;
   	o->a = 3.0f;
   	o->t = 0.01f;
}

inline void linzd_iterate(linzd_t* o)
{
	o->x += o->t * (-o->y);
	o->y += o->t * (o->x + o->z);
    o->z += o->t * (o->x * o->z + o->a*o->y*o->y);
}

////////////////////////////////////////////////////////////////////////////////////////
// Sprott 6-term ///////////////////////////////////////////////////////////////////////
void sprottst_init(sprottst_d* o)
{
	o->x = 0.1f;
	o->y = 0.1f;
    o->z = 0.1f;
	o->a = 0.8f;
	o->b = 0.5f;
    o->c = 0.1f;
    o->d = 1.0f;
    o->t = 0.01f;
}

inline void sprottst_iterate(sprottst_d* o)
{
	o->x += o->t * o->y * o->a;
	o->y += o->t * (-o->y*o->z - o->x);
    o->z += o->t * (o->b * o->y * o->y - o->c * o->x - o->d);
}


////////////////////////////////////////////////////////////////////////////////////////
// Rayleigh-Benard /////////////////////////////////////////////////////////////////////
void rayleigh_init(rayleigh_t* o)
{
    o->x = 0.01f;
	o->y = 0.0f;
    o->z = 0.0f;
	o->a = 9.00f;
	o->r = 12.0f;
    o->b = 5.00f;
    o->t = 0.19f;
}

inline void rayleigh_iterate(rayleigh_t* o)
{
	o->x = o->t * (- o->a*o->x + o->a*o->y);
	o->y = o->t * (o->r*o->x - o->y - o->x*o->z);
    o->z = o->t * (o->x*o->y - o->b*o->z);
}


////////////////////////////////////////////////////////////////////////////////////////
// Wang ////////////////////////////////////////////////////////////////////////////////
void wang_init(wang_t* o)
{
	o->x = 0.1f;
	o->y = 0.1f;
    o->z = 0.1f;
    o->w = 0.0f;
	o->a = 27.5f;
	o->b = 3.0f;
   	o->c = 19.3f;
   	o->d = 3.3f;
   	o->h = 2.9f;
   	o->t = 0.001f;
}

void wang_iterate(wang_t* o)
{
	o->x += o->t * o->a * (o->y - o->x);
	o->y += o->t * (o->b * o->x + o->c * o->y - o->x * o->z + o->w);
    o->z += o->t * (o->y * o->y - o->h * o->z);
    o->w  = o->d * -o->y;
}


////////////////////////////////////////////////////////////////////////////////////////
// Yu-Wang /////////////////////////////////////////////////////////////////////////////
void yu_wang_init(yu_wang_t* o)
{
    o->x = 0.1f;
	o->y = 0.1f;
    o->z = 0.1f;
	o->a = 10.0f;
	o->b = 40.0f;
    o->c = 2.0f;
	o->d = 2.5f;
    o->t = 0.001f;
}

void yu_wang_iterate(yu_wang_t* o)
{
	o->x += o->t * o->a * (o->y - o->x);
	o->y += o->t * (o->b * o->x - o->c * o->x * o->z);
    o->z += o->t * (powf(M_E, o->x*o->y) - o->d * o->z);
}


////////////////////////////////////////////////////////////////////////////////////////
// Three-Scroll Unified Chaotic System (TSUCS) /////////////////////////////////////////
void tsucs_init(tsucs_t* o)
{
	o->x = 1.0f;
	o->y = 1.0f;
    o->z = 1.0f;
	o->a = 40.0f;
	o->b = 0.5f;
    o->c = 20.0f;
    o->d = 0.833f;
    o->e = 0.65f;
   	o->t = 0.001;
}

void tsucs_iterate(tsucs_t* o)
{
	o->x += o->t * (o->a*(o->y-o->x) + o->b*o->x*o->z);
	o->y += o->t * (o->c*o->y - o->x*o->z);
    o->z += o->t * (o->d*o->z + o->x*o->y - o->e*o->x*o->x);
}


////////////////////////////////////////////////////////////////////////////////////////
// Three-Scroll Unified Chaotic System (TSUCS2) ////////////////////////////////////////
void tsucs2_init(tsucs2_t* o)
{
	o->x = 1.0f;
	o->y = 1.0f;
    o->z = 1.0f;
	o->a = 40.0f;
	o->b = 0.16f;
    o->c = 55.0f;
    o->d = 20.0f;
    o->e = 1.833f;
    o->f = 0.65f;
    o->t = 0.001f;
}

void tsucs2_iterate(tsucs2_t* o)
{
	o->x += o->t * (o->a*(o->y-o->x) + o->b*o->x*o->z);
	o->y += o->t * (o->c*o->y - o->x*o->z + o->d*o->y);
    o->z += o->t * (o->e*o->z + o->x*o->y - o->f*o->x*o->x);
}

////////////////////////////////////////////////////////////////////////////////////////
// Lorenz //////////////////////////////////////////////////////////////////////////////
void lorenz_init(lorenz_t* o)
{
    o->a = 10.0f;
    o->b = 28.0f;
    o->c = 8.0f / 3.0f;
    o->t = 0.01f; 
    o->x = 0.1f; 
    o->y = 0.0f;
    o->z = 0.0f; 
    o->delta = 1.0f;
}

void lorenz_iterate(lorenz_t* o)
{
    o->x += o->t * o->a * (o->y - o->x) * o->delta;
    o->y += o->t * (o->x * (o->b - o->z) - o->y) * o->delta;
    o->z += o->t * (o->x * o->y - o->c * o->z) * o->delta;
}

////////////////////////////////////////////////////////////////////////////////////////
// Aizawa //////////////////////////////////////////////////////////////////////////////
void aizawa_init(aizawa_t* o)
{
    o->a = 0.95f;
    o->b = 0.7f;
    o->c = 0.6f;
    o->d = 3.5f;
    o->e = 0.25f;
    o->f = 0.1f;
    o->t = 0.01f; 
    o->x = 0.1f; 
    o->y = 0.0f;
    o->z = 0.0f; 
}

void aizawa_iterate(aizawa_t* o)
{
    o->x += o->t * ((o->z - o->b) * o->x - o->d * o->y);
    o->y += o->t * ((o->z - o->b) * o->y + o->d * o->x);
    o->z += o->t * (o->c + o->a*o->z - o->z*o->z*o->z/3.0f 
    - (o->x*o->x+o->y*o->y)*(1.0f+o->e*o->z) 
    + o->f*o->z*o->x*o->x*o->x);
}

////////////////////////////////////////////////////////////////////////////////////////
// Halvorsen ///////////////////////////////////////////////////////////////////////////
void halvorsen_init(halvorsen_t* o)
{
    o->a = 1.4f;
    o->t = 0.01f; 
    o->x = 0.1f; 
    o->y = 0.0f;
    o->z = 0.0f; 
}

void halvorsen_iterate(halvorsen_t* o)
{
    o->x += o->t * (-o->a * o->x - 4.0f * o->y - 4.0f * o->z - o->y * o->y);
    o->y += o->t * (-o->a * o->y - 4.0f * o->z - 4.0f * o->x - o->z * o->z);
    o->z += o->t * (-o->a * o->z - 4.0f * o->x - 4.0f * o->y - o->x * o->x);
}

////////////////////////////////////////////////////////////////////////////////////////
// Ikeda ///////////////////////////////////////////////////////////////////////////////
void ikeda_init(ikeda_t* o)
{
    o->u = 0.918f;
	o->x = 0.8f;
	o->y = 0.7f;
	o->t = 0.0f;
}

void ikeda_iterate(ikeda_t* o)
{ 
    o->t  = 0.4f - 6.0f / (1.0f + o->x * o->x + o->y * o->y);
    o->x  = 1.0f + o->u * (o->x * cos(o->t) - o->y * sin(o->t));
    o->y  = o->u * (o->x * sin(o->t) + o->y * cos(o->t));
}

////////////////////////////////////////////////////////////////////////////////////////
// Duffing /////////////////////////////////////////////////////////////////////////////
void duffing_init(duffing_t* o)
{
    o->x = 0.1f;
	o->y = 0.1f;
	o->a = 2.75f;
	o->b = 0.2f;
}

void duffing_iterate(duffing_t* o)
{
	o->x = o->y;
	o->y = (-o->b*o->x + o->a*o->y - o->y*o->y*o->y);
}

////////////////////////////////////////////////////////////////////////////////////////
// Henon ///////////////////////////////////////////////////////////////////////////////
void henon_init(henon_t* o)
{
    o->x = 0.0f;
	o->y = 0.0f;
    o->dy = 0.0f;
    o->dx = 0.0f;
	o->a = 1.4f;
	o->b = 0.3f;
    o->t = 1.0f;
}

void henon_iterate(henon_t* o)
{
	o->dx = o->t * (1.0f - o->a * o->x * o->x + o->y);
	o->dy = o->t * o->b * o->x;
    o->x = o->dx;
    o->y = o->dy;
}

////////////////////////////////////////////////////////////////////////////////////////
// Gingerbreadman //////////////////////////////////////////////////////////////////////
void gingerbreadman_init(gingerbreadman_t* o)
{
    o->x = 0.0f;
    o->y = 0.0f;
    o->t = 0.001f;
}

void gingerbreadman_iterate(gingerbreadman_t* o)
{
	o->x += o->t*(1.0f - o->y + fabsf(o->x));
	o->y += o->t*o->x;
}

////////////////////////////////////////////////////////////////////////////////////////
// Van Der Pol /////////////////////////////////////////////////////////////////////////
void vanderpol_init(vanderpol_t* o)
{
	o->x = 0.1f;
	o->y = 0.1f;
	o->f = 1.2f;
    o->t = 0.1f;
	o->m = 1.0f;
}

void vanderpol_iterate(vanderpol_t* o)
{
	o->x += o->t * o->y;
	o->y += o->t * (o->m * (o->f - o->x * o->x) * o->y - o->x);
}

////////////////////////////////////////////////////////////////////////////////////////
// Kaplan-Yorke ////////////////////////////////////////////////////////////////////////
void kaplan_yorke_init(kaplan_yorke_t* o)
{
    o->a = 0xFFFFFFFF;
    o->b = 2147483647;
    o->x = 0.0;
	o->y = 0.0;
    o->t = 0.1f;
	o->alpha = 0.0;
}

void kaplan_yorke_iterate(kaplan_yorke_t* o)
{
    long f = 2 * o->a % o->b;
	o->x += o->t * ((double)o->a / (double)(o->b));
	o->y += o->t * (o->alpha*o->y + cos(4.0f * M_PI * o->x));
    o->a = f;
}

////////////////////////////////////////////////////////////////////////////////////////
// Rabinovich-Fabrikant ////////////////////////////////////////////////////////////////
void rabinovich_fabrikant_init(rabinovich_fabrikant_t* o)
{
    o->theta = 0.87f;
    o->alpha = 1.1f;
    o->x = 0.1f;
    o->y = 0.1f;
    o->z = 0.1f;
    o->t = 0.01f;
}

void rabinovich_fabrikant_iterate(rabinovich_fabrikant_t* o) 
{
    o->x += o->t * (o->y*(o->x-1.0f+o->x*o->x)+o->theta*o->x);
    o->y += o->t * (o->x*(3.0f*o->z+1.0f-o->x*o->x)+o->theta*o->y);
    o->z += o->t * (-2.0f*o->z*(o->alpha+o->x*o->y));
}

////////////////////////////////////////////////////////////////////////////////////////
// Chen-Lee ////////////////////////////////////////////////////////////////////////////
void chen_lee_init(chen_lee_t* o)
{
    o->a = 45.0f;
    o->b = 3.0f;
    o->c = 28.0f;
    o->x = 1.0f;
    o->y = 1.0f;
    o->z = 1.0f;
    o->t = 0.0035f;
}

void chen_lee_iterate(chen_lee_t* o)
{
    o->x += o->t * o->a * (o->y - o->x);
    o->y += o->t * ((o->c - o->a) * o->x - o->x*o->z + o->c*o->y);
    o->z += o->t * (o->x*o->y - o->b*o->z);
}

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

void julia_init(julia_t* o)
{
    o->zx = 1.2f;
    o->zy = 0.8f;
    o->cx = 0.2f;
    o->cy = 0.3f;
    o->t = 0.001f;
}

void julia_iterate(julia_t* o)
{
    float xm = o->t * (o->zx * o->zx - o->zy * o->zy);
    o->zy += o->t * (2.0f * o->zx * o->zy  + o->cy) ;
    o->zx += o->t * (xm + o->cx);
}



