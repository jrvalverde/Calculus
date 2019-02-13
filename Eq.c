/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*	MODULE
*	    EQUATION.C
*	Contains functions for solving equations
*
*	CONTENTS:
*	    int EqSuccApprox(sol, f, x, cerr, nmi)
*	    int EqNewton_Raphson(sol, func, derivf, x0, cerr, nmi)
*	    int EqApprNewtonRaphson(sol, func, x0, cerr, nmi)
*	    int EqRegulaFalsi(sol, f, a, b, cerr, nmi)
*	    int EqSecant(sol, f, a, b, cerr, nmi)
*	    int EqBipartition(sol, f, a, b, cerr, nmi)
*
*	DESIGNED BY:
*	    J. R. Valverde Carrillo
*
*	ULTIMA MODIFICACION:
*	    14 - dic - 1988	(implemented in IBM PC-AT)
*	    16 - dic - 1988	(port to Mac)
*	    17 - dic - 1988	(code fragmentation)
*	    08 - sep - 1993	(ported to UNIX)
*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

#include <stdio.h>

#include "portable.h"

#include <math.h>

typedef double real;	    /* To control precission */

typedef unsigned short iter_no;

typedef unsigned Sys_ord;

typedef unsigned counter;

typedef real *matrix;

typedef real *vector;

#define odd(x)		((x) & 1)

#define even(x)		(! odd(x))

/*--------------------------------------*
|					|
|	SOLUTION OF ONE EQUATION	|
|					|
*---------------------------------------*/

/* Solution of an equation by the method of successive
approximations							*/

public status EqSuccApprox(sol, f, x, cerr, nmi)
real *sol;	    /* Solution */
real (*f)();	    /* Equation to solve */
real x;		    /* Initial value */
iter_no nmi;	    /* Maximum number of iterations */
real cerr;	    /* Maximum admisible error */
  {
    real x1,	    /* value of x in the iteration i */
         x2,	    /* value of x at the iteration i + 1 */
         erit;	    /* relative error at each iteration */
    iter_no ni;	    /* number of iterations */
    
    ni = 0;
    x1 = x;
    do
      {
        x2 = (*f)(x1);
        erit = fabs((x2 - x1) * 100.0 / x2);
        x1 = x2;
        ni++;
      }
    while ((erit > cerr) && (ni < nmi));
    *sol = x2;
    return (ni >= nmi)? FAIL : SUCCESS;
  }

/* Newton - Raphson method for the solution of one equation */

public status EqNewton_Raphson(sol, func, derivf, x0, cerr, nmi)
real *sol;	    		/* solution */
real (*func)(), (*derivf)();	/* ecuation and its derivative */
real x0;	    		/* initial value */
iter_no nmi;			/* maximum number of iterations */
real cerr;	    		/* maximum error allowed */
  {
    iter_no ni;		/* no. of iterations */
    real x,		/* x value at iteration i */
         x_1,		/* x value at iteration i - 1 */
         erit;		/* relative error at each iteration */
    
    x_1 = x0;
    ni = 0;
    do
      {
        x = x_1 - (*func)(x_1) / (*derivf)(x_1);
        erit = fabs((x - x_1) * 100.0 / x);
        x_1 = x;
        ni++;
      }
    while ((erit > cerr) && (ni < nmi));
    *sol = x;
    return (ni >= nmi)? FAIL : SUCCESS;
  }

/* Newton - Rpahson 2:
    We can substitute the computation of the derivative of f, f'(xi-1) by
    
	f'(xi-1) approx= (f(xi-1) - f(xi-2)) / (xi-1 - xi-2)
	
    with which we have
    
	xi = xi-1 - (((xi-1 - xi-2) / (f(xi-1) - f(xi-2)) . f(xi-1))
*/

public status EqApprNewtonRaphson(sol, func, x0, cerr, nmi)
real *sol;	    /* solution */
real (*func)();	    /* ecuation */
real x0;	    /* initial value */
iter_no nmi;	    /* maximum number of iterations */
real cerr;	    /* allowed error */
  {
    iter_no ni;	    /* no of iterations */
    real x,	    /* x value at iteration i */
    	 x_1,	    /* x value at iteration i - 1 */
    	 x_2,	    /* x value at iteration i - 2 */
    	 fx_1,	    /* value of f(x_1) */
    	 fx_2,	    /* value of f(x_2) */
    	 erit;	    /* relative error at each iteration */
    
    x_2 = x0;
    x_1 = fx_2 = (*func)(x_2);
    ni = 0;
    do
      {
        fx_1 = (*func)(x_1);
        x = x_1 - ((x_1 - x_2) / (fx_1 - fx_2)) * fx_1;
        erit = fabs((x - x_1) * 100.0 / x);
        x_2 = x_1;
        x_1 = x;
        fx_2 = fx_1;
        ni++;
      }
    while ((erit > cerr) && (ni < nmi));
    *sol = x;
    return (ni >= nmi)? FAIL : SUCCESS;
  }

/*  Regula Falsi Method

    Let f be an equation supposed continuous in the interval [a, b]
    and with only one radix f(x) = 0 in [a, b] => f(a) * f(b) < 0.
     from {a, f(a)} and {b, f(b)}
	x0 = b - (((b - a) f(b)) / (f(b) - f(a)))
    We verify if f(x0) . f(a) < 0
    If so, we have b = x0 and repeat the computation, otherwise a = x0

    Convergence is slower than with Newton-Raphson, but it doesn't
    require computation of the derivative of f.
*/

public status EqRegulaFalsi(sol, f, a, b, cerr, nmi)
real *sol;
real (*f)();
real a, b, cerr;
iter_no nmi;
  {
    real inf, sup, oldmed, med, finf, fsup, fmed, erit;
    iter_no ni;
    
    if (fsign((*f)(a)) == fsign((*f)(b)))
      return FAIL;
      /* If not, then it has no solution or more than one */
    
    ni = 0;
    inf = a;
    sup = b;
    oldmed = a;
    do
      {
        finf = (*f)(inf);
        fsup = (*f)(sup);
        med = sup - (((sup - inf) * fsup) / (fsup - finf));
        fmed = (*f)(med);
        erit = fabs((med - oldmed) * 100.0 / med);
        oldmed = med;
        if ((fmed > 0.0) && (finf < 0.0))
          sup = med;
        else
          inf = med;
        ni++;
      }
    while ((erit > cerr) && (ni < nmi));
    *sol = med;
    return (ni >= nmi)? FAIL : SUCCESS;
  }

/* Secant method */

public status EqSecant(sol, f, a, b, cerr, nmi)
real *sol;
real (*f)();
real a, b, cerr;
iter_no nmi;
  {
    real xn_1,	    /* x(n - 1) */
    	 xn,	    /* x(n) */
    	 xn1,	    /* x(n + 1) */
    	 fxn_1, fxn, erit;
    iter_no ni;
    
    if (fsign((*f)(a)) == fsign((*f)(b)))
      return FAIL;
    
    ni = 0;
    xn_1 = a;
    xn = b;
    do
      {
        fxn = (*f)(xn);
        fxn_1 = (*f)(xn_1);
        xn1 = ((fxn * xn_1) - (fxn_1 * xn)) / (fxn - fxn_1);
        erit = fabs((xn1 - xn) * 100.0 / xn1);
        xn_1 = xn;
        xn = xn1;
        ni++;
      }
    while ((erit > cerr) && (ni < nmi));
    *sol = xn1;
    return (ni >= nmi)? FAIL : SUCCESS;
  }

/* Bipartition method
    If f(x) is such that it is monotonous and continuous in
    [a, b] and that f(a) and f(b) have different signs then
    there is a unique solution X* in [a, b].
*/

public status EqBipartition(sol, f, a, b, cerr, nmi)
real *sol;
real (*f)();
real a, b, cerr;
iter_no nmi;
  {
    real inf, sup, med, old,
    	 finf, fmed, erit;
    iter_no ni;
    
    if (fsign((*f)(a)) == fsign((*f)(b)))
      return FAIL;
    ni = 0;
    inf = a;
    sup = b;
    old = a;
    do
      {
        med = (inf + sup) / 2.0;
        finf = (*f)(inf);
        fmed = (*f)(med);
        erit = fabs((med - old) * 100.0 / med);
        if (fsign(finf) == fsign(fmed))
          old = inf = med;
        else
          old = sup = med;
        ni++;
      }
    while ((erit > cerr) && (ni < nmi));
    *sol = med;
    return (ni >= nmi)? FAIL : SUCCESS;
  }

