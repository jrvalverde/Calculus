/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*	COMPLEX.C
*
*	Procedures for portable complex number manipulation.
*
*	Contents:
*	    c_set(c, real, imag)
*	    p_set(z, modulus, argument)
*	    c_comp(a, b)
*	    c_topolar(r, z)
*	    c_copy(a, b)
*	    c_add(r, a, b)
*	    c_subs(r, a, b)
*	    c_mult(r, a, b)
*	    c_conj(r, z)
*	    c_inv(r, z)
*	    c_div(r, a, b)
*	    c_mod(r, a)
*	    c_exp(r, z)
*	    c_log(r, z)
*	    c_logx(r, z, n)
*	    c_pow(r, z, w)
*	    c_radix(r, z, w)
*	    c_rpow(r, z, n)
*	    c_rradix(r, z, n, order)
*	    c_sin(r, z)
*	    c_cos(r, z)
*	    c_tan(r, z)
*	    p_tocomplex(r, z)
*	    p_copy(dest, src)
*	    p_mult(r, a, b)
*	    p_div(r, a, b)
*	    p_rpow(r, a, n)
*	    p_rradix(r, a, n, order)
*	    p_mod(r, z)
*	    ComplError()
*	    getComplError()
*	    resetComplError()
*	    setComplError()
*
*	Uses:
*	    PORTABLE.H
*	    MATH.H
*
*	Caveats:
*	    NOT DEBUGGED.
*
*	Designed by: J. R. Valverde Carrillo
*
*	Last modified:
*	    17 - oct - 1988  (written as funcs)
*	    29 - oct - 1988  (compiles without errors)
*	    19 - nov - 1988  (transformed to procs)
*	    14 - dec - 1988  (ported from IBM-AT to Mac)
*	    16 - dec - 1988  (error checking)
*	    17 - dec - 1988  (various embellishments)
*	    07 - sep - 1993  (ported to UNIX)
*
*	Copyright:
*	    © YoEgo.This is PUBLIC DOMAIN SOFTWARE
*	Using it for commercial purposes would be a nasty thing
*	Don't even dream of doing that without prior permission.
*	Never relay on any oral consent, consents are only valid
*	if written.
*						    YoEgo.
*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

#include "portable.h"

#include <math.h>

typedef double real;	/* to control precission */

typedef struct
    {
	real r, i;
    } complex;
    
typedef struct
    {
	real mod,	/* mod = sqrt(r * r + i * i)	*/
	     arg;	/* arg = atan(i / r)		*/
    } polar;
    
typedef union
    {
	complex c;
	polar	p;
    } Complex;

static Complex aux, aux2;
    
static int _cError = 0;

#define NOCERROR	 0
#define DIVBYZERO	-1
#define NEGSQROOT	-2
#define LOGOFNEG	-3
#define NEGRPOWER	-4
#define NEGRRADIX	-5


public void c_set(comp, p_real, p_imag)
/* initializes a complex number */
complex *comp;
real p_real, p_imag;
  {
    comp->r = p_real;
    comp->i = p_imag;
  }

public void p_set(pol, modulus, angle)
/* initializes a complex number in polar form */
polar *pol;
real modulus, angle;
  {
    pol->mod = modulus;
    pol->arg = angle;
  }

public void c_topolar(pol, z)
/* converts a complex number to polar form */
polar *pol;
complex *z;
  {
    extern Complex aux;    /* just in case pol == z */
    extern int _cError;
    real tmp;
    
    if (z->r == 0.0)
      {
        _cError = DIVBYZERO;
        return;
      }
    tmp = (z->r * z->r) + (z->i * z->i);
    if (tmp < 0.0)
      {
        _cError = NEGSQROOT;
        return;
      }
    aux.p.mod = sqrt(tmp);
    aux.p.arg = atan(z->i / z->r);
    
    pol->mod = aux.p.mod;
    pol->arg = aux.p.arg;
  }

public bool c_comp(a, b)
/* verifies if a and b are equal */
complex *a, *b;
  {
    if ((a->r == b->r) && (a->i == b->i))
      return TRUE;
    else
      return FALSE;
  }

public void c_copy(a, b)
/* copies number b over number a */
complex *a, *b;
  {
    a->r = b->r;
    a->i = b->i;
  }

public void c_add(sum, a, b)
/* adds a to b and stores result in sum */
complex *sum, *a, *b;
  {
    sum->r = a->r + b->r;
    sum->i = a->i + b->i;
  }
  
public void c_subs(res, min, sus)
/* res = minuend - substraend */
complex *res, *min, *sus;
  {
    res->r = min->r - sus->r;
    res->i = min->i - sus->i;
  }

public void c_mult(prod, a, b)
/* multiplies a by b and stores result in prod */
complex *prod, *a, *b;
  {
    extern Complex aux;
    
    aux.c.r = (a->r * b->r) - (a->i * b->i);
    aux.c.i = (a->r * b->i) + (a->i * b->r);
    prod->r = aux.c.r;
    prod->i = aux.c.i;
  }

public void c_conj(cnj, a)
/* finds the conjugate of complex number a */
complex *cnj, *a;
  {
    cnj->r = a->r;
    cnj->i = -(a->i);
 }
  
public void c_inv(in, a)
/* computes the invert of a */
complex *in, *a;
  {
    extern Complex aux;
    extern int _cError;
    real tmp;
    
    tmp = (a->r * a->r) + (a->i * a->i);
    if (tmp == 0)
      {
        _cError = DIVBYZERO;
        return;
      }  
    aux.c.r = a->r /tmp;
    aux.c.i = - (a->i) /tmp;
    in->r = aux.c.r;
    in->i = aux.c.i;
  }

public void c_div(coc, a, b)
complex *coc, *a, *b;
/* computes the quotient a / b */
  {
    extern Complex aux;
    extern int _cError;
    
    c_inv(&aux, b);
    if (_cError)
      {
        /* _cError is already set */
        return;
      }
    c_mult(&aux, a, &aux);
    coc->r = aux.c.r;
    coc->i = aux.c.i;
  }

public void c_mod(mod, a)
real *mod;
complex *a;
/* computes the modulus of complex number a */
  {
    real tmp;
    extern int _cError;
    
    tmp = (a->r * a->r) + (a->i * a->i);
    if (tmp < 0.0)
      {
        _cError = NEGSQROOT;
        return;
      }
    *mod = sqrt(tmp);
  }

public void c_exp(e_to_z, z)
complex *e_to_z, *z;
/* computes exp(z) being z a complex number:
	 z
	e				    */
  {
    real tmp;
    
    tmp = exp(z->r);
    e_to_z->r = tmp * cos(z->i);
    e_to_z->i = tmp * sin(z->i);
  }

public void c_log(lo, z)
complex *lo, *z;
/* Computes the neperian logarithm of z */
  {
    extern Complex aux;
    extern int _cError;
    
    c_topolar(&aux, z);
    if (aux.p.mod < 0.0)
      {
        _cError = LOGOFNEG;
        return;
      }
    lo->r = log(aux.p.mod);
    lo->i = aux.p.arg;
  }

public void c_logx(lo, base, z)
real base;
complex *lo, *z;
/* finds logarithm lo in any base of z */
  {
    extern Complex aux;
    extern Complex aux2;
    extern int _cError;
    
    if (base < 0.0)
      {
        _cError = LOGOFNEG;
        return;
      }
    aux2.c.r = log(base);
    aux2.c.i = 0.0;
    c_log(&aux, z);
    if (_cError)
      {
        return;
        /* _cError is already set */
      }
    c_div(&aux2, &aux, &aux2);
    if (_cError)
      {
        return;
      }
    lo->r = aux2.c.r;
    lo->i = aux2.c.i;
  }

public void c_pow(p, z, w)
complex *p, *z, *w;
/* 	 w
    p = z	*/
  {
    extern Complex aux;
    extern int _cError;
    
    c_log(&aux, z);
    if (_cError)
      {
        return;
      }
    c_mult(&aux, w, &aux);
    c_exp(&aux);
    if (_cError)
      {
        return;
      }
    p->r = aux.c.r;
    p->i = aux.c.i;
  }

public void c_radix(ra, z, w)
complex *ra, *z, *w;
/* finds the w complex radix of z */
  {
    extern Complex aux;
    extern int _cError;
    
    c_log(&aux, z);
    if (_cError)
      return;
    c_div(&aux, &aux, w);
    if (_cError)
      return;
    c_exp(&aux, &aux);
    if (_cError)
      return;
    ra->r = aux.c.r;
    ra->i = aux.c.i;
  }

public void c_rpow(rp, z, n)
complex *rp, *z;
real n;
/* raises z to p being z complex and p real */
  {
    extern Complex aux;
    extern int _cError;
    extern p_rpow();
    extern p_tocomplex();
    
    c_topolar(&aux, z);
    if (_cError) return;
      
    p_rpow(&aux, &aux, n);
    if (_cError) return;
    
    p_tocomplex(&aux, &aux);
    
    rp->r = aux.c.r;
    rp->i = aux.c.i;
  }

public void c_rradix(ra, z, n, order)
complex *ra, *z;
real n;
unsigned int order;
/* finds ra, the real radix of order order of a complex z */
  {
    extern Complex aux;
    extern int _cError;
    extern p_rradix();
    extern p_tocomplex();
    
    c_topolar(&aux, z);
    if (_cError) return;
    
    p_rradix(&aux, z, n, order);
    if (_cError) return;
    
    p_tocomplex(&aux, &aux);
    
    ra->r = aux.c.r;
    ra->i = aux.c.i;
  }

public void c_sin(s, z)
complex *s, *z;
/* finds the sinus of a complex number z */
  {
    extern Complex aux;
    
    aux.c.r = sin(z->r) * cosh(z->i);
    aux.c.i = cos(z->r) * sinh(z->i);
    s->r = aux.c.r;
    s->i = aux.c.i;
  }

public void c_cos(co, z)
complex *co, *z;
/* cosinus of z */
  {
    extern Complex aux;
    
    aux.c.r = cos(z->r) * cosh(z->i);
    aux.c.i = - (sin(z->r) * sinh(z->i));
    co->r = aux.c.r;
    co->i = aux.c.i;
  }

public void c_tan(tan, z)
real *tan;
complex *z;
/* computes the tangent of z; angles are in radians */
  {
    real aux;

    *tan = (sin(z->r + z->r) + sinh(z->i + z->i));
    aux = (cos(z->r + z->r) + cosh(z->i + z->i));
    if (aux == 0) {
	_cError = DIVBYZERO;
	return;
    }
    *tan /= aux;
  }

public void p_tocomplex(c, z)
complex *c;
polar *z;
/* converts a polar number to complex form */
  {
    extern Complex aux;
    
    aux.c.r = z->mod * cos(z->arg);
    aux.c.i = z->mod * sin(z->arg);
    c->r = aux.c.r;
    c->i = aux.c.i;
  }

public void p_copy(p1, p2)
polar *p1, *p2;
/* copies a polar number over other */
  {
    p1->mod = p2->mod;
    p1->arg = p2->arg;
  }

public void p_mult(prod, z1, z2)
polar *prod, *z1, *z2;
/* multiplies two complex numbers in polar form */
  {
    prod->mod = z1->mod * z2->mod;
    prod->arg = z1->arg + z2->arg;
  }

public void p_div(coc, z1, z2)
polar *coc, *z1, *z2;
/* divides two complex numbers in polar form */
  {
    if (z2->mod == 0.0)
      {
        _cError = DIVBYZERO;
        return;
      }
    coc->mod = z1->mod / z2->mod;
    coc->arg = z1->arg + z2->arg;
  }

public void p_mod(modulus, z)
real *modulus;
polar *z;
/* returns the modulus of a complex number in polar form */
  {
    *modulus = z->mod;
  }

public void p_rpow(p, z, n)
polar *p, *z;
real n;
/* real power of a complex number in polar form */
  {
    extern int _cError;
    
    if (n < 0.0)
      {
        _cError = NEGRPOWER;
        return;
      }
    p->mod = pow(z->mod, n);
    p->arg = z->arg * n;
    
  }

public void p_rradix(pr, z, n, order)
polar *pr, *z;
real n;
unsigned int order;
/* returns radix of order order of the n-th radix of z:
	 (1/n)
	z	    (order = 0, ... , n - 1)		*/
  {
    extern int _cError;
    
    if (n <= 0.0)
      {
        _cError = NEGRRADIX;
        return;
      }
    pr->mod = pow(z->mod, 1 / n);
    pr->arg = (z->arg / n) + ((360 * order) / n);
    return;
  }

public int ComplError()
  {
    return _cError;
  }

public int getComplError()
/* returns the last error occurred in complex number
   calculations and resets the error switch */
  {
    extern int _cError;
    register int theError;
    
    theError = _cError;
    _cError = NOCERROR;
    return theError;
    
  }

public void resetComplError()
/* resets the complex error calculus switch */
  {
    extern int _cError;
    
    _cError = NOCERROR;
  }

public void setComplError(err)
int err;
/* sets the error switch to the specified value */
  {
    extern int _cError;
    
    _cError = err;
  }
