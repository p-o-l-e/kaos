
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
#include <cmath>

////////////////////////////////////////////////////////////////////////////////////////
// Roessler ////////////////////////////////////////////////////////////////////////////
struct roessler
{
	float	x = 1.0f;
	float	y = 1.0f;
	float	z = 1.0f;

	float	a = 0.2f;
	float	b = 0.2f;
	float	c = 5.7f;

	float	delta = 0.01f;
    
    void    iterate();
};

void roessler::iterate()
{
	x += (-y - z)*delta;
	y += (x + a * y)*delta;
	z += (b + z * (x - c))*delta;
}

////////////////////////////////////////////////////////////////////////////////////////
// Hopf ////////////////////////////////////////////////////////////////////////////
struct hopf
{
	float	x = 0.01f;
	float	y = 0.01f;

	float	p = 0.11f;
    float   t = 0.01f;

    void    iterate();
};

void hopf::iterate()
{
	x += t * ( -y + x * (p - (x*x + y*y)));
	y += t * (  x + y * (p - (x*x + y*y)));
}

////////////////////////////////////////////////////////////////////////////////////////
// Helmholz ////////////////////////////////////////////////////////////////////////////
struct helmholz
{
	float	x = 0.1f;
	float	y = 0.1f;
    float   z = 0.1f;

	float	gamma = 5.11f;
	float	delta = 0.55f;
    float   t     = 0.01;
    
    void    iterate();
};

void helmholz::iterate()
{
	x += t * y;
	y += t * gamma * z;
    z += t * ( -z - delta * y - x - x * x );
}


////////////////////////////////////////////////////////////////////////////////////////
// Sprott ////////////////////////////////////////////////////////////////////////////
struct sprott
{
	float	x = 0.1f;
	float	y = 0.1f;
    float   z = 0.1f;

    float   t = 0.1;
    
    void    iterate();
};

void sprott::iterate()
{
	x += t * y;
	y += t * (y * z - x);
    z += t * (1.0f - y * y);
}

////////////////////////////////////////////////////////////////////////////////////////
// Sprott-Linz ////////////////////////////////////////////////////////////////////////////
struct linz
{
	float	x = 0.1f;
	float	y = 0.1f;
    float   z = 0.1f;

    float   a = 0.5f;
    float   t = 0.1;
    
    void    iterate();
};

void linz::iterate()
{
	x += t * (y + z);
	y += t * (y * a - x);
    z += t * (x * x - z);
}

////////////////////////////////////////////////////////////////////////////////////////
// Sprott-Linz D ////////////////////////////////////////////////////////////////////////////
struct linz_d
{
	float	x = 0.1f;
	float	y = 0.1f;
    float   z = 0.1f;

    float   a = 3.0f;
    float   t = 0.01;
    
    void    iterate();
};

void linz_d::iterate()
{
	x += t * (-y);
	y += t * (x + z);
    z += t * (x * z + a*y*y);
}

////////////////////////////////////////////////////////////////////////////////////////
// Sprott 6-term ////////////////////////////////////////////////////////////////////////////
struct sprott_st
{
	float	x = 0.1f;
	float	y = 0.1f;
    float   z = 0.1f;

	float	a = 0.8f;
	float	b = 0.5f;
    float   c = 0.1f;
    float   d = 1.0f;
    float   t = 0.01f;
    
    void    iterate();
};

void sprott_st::iterate()
{
	x += t * y * a;
	y += t * (- y*z - x);
    z += t * (b * y * y - c * x - d);
}


////////////////////////////////////////////////////////////////////////////////////////
// Rayleigh-Benard ////////////////////////////////////////////////////////////////////////////
struct rayleigh
{
	double	x = 0.01f;
	double	y = 0.0f;
    double  z = 0.0f;

	float	a = 9.00f;
	float	r = 12.0f;
    float   b = 5.00f;

    float   t = 0.19f;
    
    void    iterate();
};

void rayleigh::iterate()
{
	x = t * (- a*x + a*y);
	y = t * (r*x - y - x*z);
    z = t * (x*y - b*z);

}


////////////////////////////////////////////////////////////////////////////////////////
// Wang ////////////////////////////////////////////////////////////////////////////
struct wang
{
	float	x = 0.1f;
	float	y = 0.1f;
    float   z = 0.1f;
    float   w;

	float	a = 27.5f;
	float	b = 3.0f;
    float   c = 19.3f;
    float   d = 3.3f;
    float   h = 2.9f;
    float   t = 0.001;
    
    void    iterate();
};

void wang::iterate()
{
	x += t * a * (y - x);
	y += t * (b * x + c * y - x * z + w);
    z += t * (y * y - h * z);
    w  = d * -y;
}


////////////////////////////////////////////////////////////////////////////////////////
// Yu-Wang ////////////////////////////////////////////////////////////////////////////
struct yu_wang
{
	float	x = 0.1f;
	float	y = 0.1f;
    float   z = 0.1f;

	float	a = 10.0f;
	float	b = 40.0f;
    float   c = 2.0f;
    float   d = 2.5f;

    float   t = 0.001;
    
    void    iterate();
};

void yu_wang::iterate()
{
	x += t * a * (y - x);
	y += t * (b * x - c * x * z);
    z += t * (pow(M_E, x*y) - d * z);

}


////////////////////////////////////////////////////////////////////////////////////////
// Three-Scroll Unified Chaotic System (TSUCS) ////////////////////////////////////////////////////////////////////////////
struct tsucs
{
	float	x = 1.0f;
	float	y = 1.0f;
    float   z = 1.0f;

	float	a = 40.0f;
	float	b = 0.5f;
    float   c = 20.0f;
    float   d = 0.833f;
    float   e = 0.65f;

    float   t = 0.001;
    
    void    iterate();
};

void tsucs::iterate()
{
	x += t * (a*(y-x) + b*x*z);
	y += t * (c*y - x*z);
    z += t * (d*z + x*y - e*x*x);
}


////////////////////////////////////////////////////////////////////////////////////////
// Three-Scroll Unified Chaotic System (TSUCS2) ////////////////////////////////////////////////////////////////////////////
struct tsucs2
{
	float	x = 1.0f;
	float	y = 1.0f;
    float   z = 1.0f;

	float	a = 40.0f;
	float	b = 0.16f;
    float   c = 55.0f;
    float   d = 20.0f;
    float   e = 1.833f;
    float   f = 0.65f;

    float   t = 0.001;
    
    void    iterate();
};

void tsucs2::iterate()
{
	x += t * (a*(y-x) + b*x*z);
	y += t * (c*y - x*z + d*y);
    z += t * (e*z + x*y - f*x*x);
}

////////////////////////////////////////////////////////////////////////////////////////
// Lorenz //////////////////////////////////////////////////////////////////////////////
struct lorenz
{
        float a = 10.0;
        float b = 28.0;
        float c = 8.0 / 3.0;
        float t = 0.01; 

        float x = 0.1; 
        float y = 0;
        float z = 0; 

        float delta = 1.0f;

        void iterate()
        {
                x += t * a * (y - x) * delta;
                y += t * (x * (b - z) - y) * delta;
                z += t * (x * y - c * z) * delta;
        }
};

////////////////////////////////////////////////////////////////////////////////////////
// Aizawa //////////////////////////////////////////////////////////////////////////////
struct aizawa
{
        float a = 0.95f;
        float b = 0.7f;
        float c = 0.6f;
        float d = 3.5f;
        float e = 0.25f;
        float f = 0.1f;

        float t = 0.01f; 

        float x = 0.1; 
        float y = 0;
        float z = 0; 

        void iterate()
        {
                x += t * ((z - b) * x - d * y);
                y += t * ((z - b) * y + d * x);
                z += t * (c + a*z - z*z*z/3.0f - (x*x+y*y)*(1.0f+e*z) + f*z*x*x*x);
        }
};


////////////////////////////////////////////////////////////////////////////////////////
// Halvorsen //////////////////////////////////////////////////////////////////////////////
struct halvorsen
{
        float a = 1.4;
        float t = 0.01; 

        float x = 0.1; 
        float y = 0;
        float z = 0; 



        void iterate()
        {
                x += t * (-a * x - 4.0f * y - 4.0f * z - y * y);
                y += t * (-a * y - 4.0f * z - 4.0f * x - z * z);
                z += t * (-a * z - 4.0f * x - 4.0f * y - x * x);
        }
};


////////////////////////////////////////////////////////////////////////////////////////
// Ikeda ///////////////////////////////////////////////////////////////////////////////
struct ikeda
{

	float u = 0.918;
	float x = 0.8;
	float y = 0.7;
	float t;

    void iterate()
	{ 
        t  = 0.4f - 6.0f / (1.0f + x * x + y * y);
        x  = 1.0f + u * (x * cos(t) - y * sin(t));
        y  = u * (x * sin(t) + y * cos(t));
	}
};

////////////////////////////////////////////////////////////////////////////////////////
// Duffing /////////////////////////////////////////////////////////////////////////////
struct duffing
{
	float x = 0.1f;
	float y = 0.1f;

	float a = 2.75f;
	float b = 0.2f;

	void iterate()
	{
		x = y;
		y = (-b*x + a*y - y*y*y);
	}
};

////////////////////////////////////////////////////////////////////////////////////////
// Henon ///////////////////////////////////////////////////////////////////////////////
struct henon
{
	float x;
	float y;
    float dy, dx;

	float a = 1.4f;
	float b = 0.3f;

    float t = 1.0f;

	void iterate()
	{
		dx = t * (1.0f - a * x * x + y);
		dy = t * b * x;

        x = dx;
        y = dy;
	}
};


////////////////////////////////////////////////////////////////////////////////////////
// Gingerbreadman //////////////////////////////////////////////////////////////////////
struct gingerbreadman
{
	float x;
	float y;

	void iterate()
	{
		x = 1.0f - y + abs(x);
		y = x;
	}
};

////////////////////////////////////////////////////////////////////////////////////////
// Van Der Pol /////////////////////////////////////////////////////////////////////////
struct vanderpol
{
	float x = 0.1f;
	float y = 0.1f;
	float f = 1.2f;
    float t = 0.1f;
	float m = 1.0f;

	void iterate();
};

void vanderpol::iterate()
{
	x += t * y;
	y += t * (m * (f - x * x) * y - x);
}

////////////////////////////////////////////////////////////////////////////////////////
// Kaplan-Yorke ////////////////////////////////////////////////////////////////////////
struct kaplan_yorke
{
	float x;
	float y;
    int a   = 0xFFFFFFFF;
    int b   = 2147483647;
    float t = 0.1f;
	float alpha;

	void iterate();
};

void kaplan_yorke::iterate()
{
    int aa = 2 * a % b;
	x += t * (float(a) / float(b));
	y += t * (alpha*y + cos(4.0f * M_PI * x));
    a = aa;
}

////////////////////////////////////////////////////////////////////////////////////////
// Rabinovich-Fabrikant ////////////////////////////////////////////////////////////////
struct rabinovich_fabrikant
{
    float gamma = 0.87f;
    float alpha = 1.1f;
    float x = 0.1f, y = 0.1f, z = 0.1f;
    float t = 0.01f;

    void iterate();

};

void rabinovich_fabrikant::iterate() 
{
    x += t * (y*(x-1+x*x)+gamma*x);
    y += t * (x*(3*z+1-x*x)+gamma*y);
    z += t * (-2*z*(alpha+x*y));
}

////////////////////////////////////////////////////////////////////////////////////////
// Chen-Lee ////////////////////////////////////////////////////////////////////////////
struct chen_lee
{
    float a = 45.0f, b = 3.0f, c = 28.0f;
    float x = 1, y = 1, z = 1;
    float t = 0.0035f;

    void iterate();
};

void chen_lee::iterate()
{
    x += t * a * (y - x);
    y += t * ((c - a) * x - x*z + c*y);
    z += t * (x*y - b*z);
}

////////////////////////////////////////////////////////////////////////////////////////
// Chua ////////////////////////////////////////////////////////////////////////////////
struct chua
{
    float x = 1, y = 1, z = 1;  // Ins

    float alpha  =  15.6;
    float beta   =  28; 
    float ma     = -1.143;
    float mb     = -0.714;
    float h;
    float t = 0.1f;
    
    float dx;
    float dy;
    float dz;
   
    void iterate();
};

void chua::iterate()
{
    h  = ma * x + 0.5f * (ma - mb) * (abs(x + 1.0f) - abs(x - 1.0f));
    dx = t * (alpha * (y - x - h));
    dy = t * (x - y + z);
    dz = t * (- beta * y);

    x=dx;
    y=dy;
    z=dz;
}

////////////////////////////////////////////////////////////////////////////////////////
// Chua circuit emulation //////////////////////////////////////////////////////////////
struct realchua
{
    float x, y, z;         // In
    float C1    = 10e-9;   // 10nF
    float C2    = 100e-9;  // 100nF
    float R     = 1800;    // 1.8k Ohms
    float G     = 1/R;

    // Chua Diode //////////////////////
    float R1    = 220;
    float R2    = 220;
    float R3    = 2200;
    float R4    = 22000;
    float R5    = 22000;
    float R6    = 3300;

    float Esat  = 9; // 9V phasecvbatteries
    float E1    = R3/(R2+R3)*Esat;
    float E2    = R6/(R5+R6)*Esat;

    float m12   = -1/R6;
    float m02   = 1/R4;
    float m01   = 1/R1;
    float m11   = -1/R3;

    float m0;
    float m1    = m12+m11;

    // Gyrator ////////////////////
    float R7    = 100;  // 100 Ohms
    float R8    = 1000; // phasecv 1k Ohms
    float R9    = 1000; // 1k Ohms
    float R10   = 1800;
    float C     = 100*10^(-9);    // 100nF
    float L     = R7*R9*C*R10/R8; // 18mH 

    float xdot;
    float ydot;
    float zdot;
    
    void iterate();
};


void realchua::iterate()
{
    if(E1>E2)
    m0 = m11 + m02;
    else
    m0 = m12 + m01;   

    float mm1 = m01 + m02;
    float Emax = std::max(E1, E2);
    float Emin = std::min(E1, E2);

    float g;
    if (abs(x) < Emin) g = x*m1;     
    else if(abs(x) < Emax )
    {
        g = x*m0;
        if (x > 0) g = g + Emin*(m1-m0);    
        else g = g + Emin*(m0-m1);  
    }

    else if (abs(x) >= Emax)
    {
        g = x*mm1;    
        if (x > 0)
            g = g + Emax*(m0-mm1) + Emin*(m1-m0);
        else
            g = g + Emax*(mm1-m0) +  Emin*(m0-m1);
    }
    xdot = (1/C1)*(G*(y-x)-g);
    ydot = (1/C2)*(G*(x-y)+z);
    zdot  = -(1/L)*y;
}

struct julia
{
    float zx = 1.2, zy = 0.8;
    float cx = 0.2, cy = 0.3;

    float t = 0.1f;

    void iterate()
    {
        float   xm = t * (zx * zx - zy * zy);
                zy = t * (2.0f * zx * zy  + cy) ;
                zx = t * (xm + cx);
    }

};




