/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*	MODULE
*	    INTEGRAL.C
*	Module to compute integrals of functions
*
*	CONTENTS:
*	    real IgTrapecium(f, xinf, xsup, ni)
*	    real IgSimpson(f, xinf, xsup, ni)
*
*	DESIGNED BY:
*	    Jos� Ram�n Valverde Carrillo
*
*	ULTIMA MODIFICACION:
*	    14 - dic - 1988 (Implemented on IBM PC-AT)
*	    16 - dic - 1988 (Ported to Mac)
*	    17 - dic - 1988 (Library split into modules)
*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

#include <stdio.h>

#include <portable.h>

#include <math.h>

typedef double real;	    /* to control precision */

typedef unsigned short iter_no;

typedef unsigned Sys_ord;

typedef unsigned count;

typedef unsigned counter;

typedef real *matrix;

typedef real *vector;

#define odd(x)		((x) & 1)
	/* zeroes all bits but the least significant. If the number
	is even, then the result is 0 and returns FALSE, but if it
	was odd, then the result is 1 and returns TRUE */

#define even(x)		(! odd(x))


/*--------------------------*
|			    |
|	INTEGRATION	    |
|			    |
*---------------------------*/

/* Integration by the trapezoidal rule */

public real IgTrapecium(f, xinf, xsup, ni)
real (*f)();	    /* tunction to integrate */
real xinf, xsup;    /* upper and lower limits of integration interval */
count ni;	    /* number of integration subintervals */
  {
    real y1,	    /* value of y sub i */
         y2,	    /* value of y sub i + 1 */
         sum,	    /* value ofl summatorium */
         h,	    /* interval */
         x;	    /* value to compute */
    counter i;
    
    y1 = (*f)(xinf);
    h = (xsup - xinf) / ni;
    sum = 0.0;
    x = xinf;
    for (i = 1; i <= ni; i++)
      {
        x += h;
        y2 = (*f)(x);
        sum += h * (y1 + y2) / 2.0;
        y1 = y2;
      }
    return sum;
  }

/* integration by Simpson's method */

public real IgSimpson(f, xinf, xsup, ni)
real (*f)();	    /* funtion to integrate */
real xinf, xsup;    /* upper and lower limits of integration interval */
count ni;	    /* number of integration subintervals */
  {
    counter i, nt;
    real x, y, h, sum;
    
    x = xinf;
    if (odd(ni))
      /* incorrect number of tracts. Must be even */
      nt = ni + 1;
    else
      nt = ni;
    
    /* Proceso */
    sum = (*f)(x);
    h = (xsup - xinf) / nt;
    for (i = 2; i <= nt; i++)
      {
        x += h;
        y = (*f)(x);
        sum += (odd(i))? (4.0 * y) : (2.0 * y);
      }
    sum += (*f)(x + h) * h / 3.0;
    
    return sum;
  }

