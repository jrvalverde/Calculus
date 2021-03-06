/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*	MODULE
*	    ODE.C
*	Functions for the solution of Ordinary Differential Equations
*
*	CONTENTS:
*	    ODEEuler()
*	    ODEEuler_Gauss()
*	    ODERunge_Kutta()
*
*	USES:
*	    STDIO.H
*	    PORTABLE.H
*	    MATH.H
*
*	DESIGNED BY:
*	    J. R. Valverde Carrillo
*
*	LAST MODIFIED:
*	    25 - feb - 1989
*	    08 - sep - 1993
*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/


#include <stdio.h>

#include "portable.h"

#include <math.h>

typedef double real;		/* To control precission */

typedef unsigned short iter_no;

typedef unsigned Sys_ord;

typedef unsigned count;

typedef unsigned counter;

typedef real *matrix;

typedef real *vector;

#define odd(x)		((x) & 1)

#define even(x)		(! odd(x))

/*----------------------------------------------*
|						|
|	ORDINARY DIFFERENTIAL EQUATIONS 	|
|						|
*-----------------------------------------------*/

/*
    Euler's method for the approximate solution
    of ordinary differential equations.
*/

public void ODEEuler(x0, y0, uvx, np, vals, func)
real x0;	    /* initial value of x */
real y0;	    /* initial value of y */
real uvx;	    /* last value of x */
count np;	    /* number of points */
real *vals;	    /* computed (0..np) values of y */
real (*func)();	    /* function to integrate */
  {
    real inter;	    /* interval */
    real x, y;
    counter i;
    
    x = x0;
    y = y0;
    inter = (uvx - x) / np;
    vals[0] = y0;
    for (i = 1; i < np; i++)
      {
        y += inter * (*func)(x, y);
        x += inter;
        vals[i] = y;
      }
  }

public void ODEEuler_Gauss(x, uvx, y, np, cerr, nmc, vals, func)
real x;		/* initial value of x */
real uvx;	/* last value of x */
real y;		/* initial value of y */
int np;		/* number of points */
real cerr;	/* allowed error */
count nmc;	/* maximum number of corrections */
real *vals;	/* set of computed y values */
real (*func)();
  {
    real x1,	/* value of x at point i in each iteration */
    	 x2,	/* value of x in point i + 1 ... */
    	 y1,	/* value of y in point i in each iteration */
    	 yp2,	/* predictor value of y in each iteration */
    	 yc2,	/* corrector value of y in each iteration */
    	 inter,	/* interval */
    	 erit,	/* relative error of the iteration */
    	 f1, f2;
    counter nc,	/* number of corrections by iteration */
    	i;	/* loop counter */
    
    x1 = x;
    y1 = y;
    inter = (uvx - x1);
    vals[0] = y1;
    for (i = 1; i < np; i++)
      {
        x2 = x1 + inter;
        yp2 = y1 + inter * (*func)(x1, y1);
        nc = 0;
        do
          {
            f1 = (*func)(x1, y1);
            f2 = (*func)(x2, yp2);
            yc2 = y1 + inter * (f1 + f2) / 2.0;
            nc++;
            yp2 = yc2;
          }
        while (nc < nmc);
        vals[i] = yp2;
        x1 = x2;
        y1 = yp2;
      }
  }

/*
    Runge - Kutta method
*/

public void ODERunge_Kutta(x, y, uvx, np, vals, func)
real x,		    /* x value in each iteration */
     y,		    /* y value in each iteration */
     uvx;	    /* last value of x */
count np;	    /* no. of points */
real vals[];	    /* Matrix with values got for y */
real (*func)();	    /* Function to be integrated */
  {
    real inter,		    /* interval */
	 t1, t2, t3, t4;    /* temporal variables */
    counter i;			/* loop counter */
    
    inter = (uvx - x) / np;
    vals[0] = y;
    for (i = 1; i < np; i++)
      {
        t1 = (*func)(x, y);
        t2 = (*func)(x + inter / 2, y + inter * t1 / 2);
        t3 = (*func)(x + inter / 2, y + inter * t2 / 2);
        t4 = (*func)(x + inter, y + inter * t3);
        y += inter * (t1 + 2 * t2 + 2 * t3 + t4) / 6;
        x += inter;
        vals[i] = y;
      }
  }

