/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*	MODULE
*	    EIGENVAL_VECT.C
*	Calculation of eigenvalues and eigenvectors of a matrix
*
*	CONTENTS:
*	    void EVLeverrierFaddev(n, a, p)
*	    int EVSuccApprox(n, a, x, avmax, cerr, nmi)
*
*	DESIGNED BY:
*	    J. R. Valverde Carrillo.
*
*	LAST MODIFICATION:
*	    14 - dic - 1988	(IBM PC-AT)
*	    16 - dic - 1988	(Mac)
*	    17 - dic - 1988	(Subdivision in small modules)
*	    07 - sep - 1993	(UNIX)
*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

#include <stdio.h>
#include <stdlib.h>

#include "portable.h"

#include <math.h>

typedef double real;	    /* To be able to control precision */

typedef unsigned short iter_no;

typedef unsigned Sys_ord;

typedef unsigned counter;

typedef real *matrix;

typedef real *vector;

#define odd(x)		((x) & 1)
	/* Sets all bits to zero but the last. If it is even it
	returns 0 (FALSE) and if not, 1 (TRUE) */

#define even(x)		(! odd(x))

/*--------------------------------------------------*
|						    |
|	CALCULO DE AUTOVALORES Y AUTOVECTORES	    |
|						    |
*---------------------------------------------------*/

/*
	Leverrier - Faddeev Method
    to get the coefficients of the characteristic
    polynomion.
*/

public status EVVLeverrierFaddeev(n, a, p, inv)
Sys_ord n;	/* Order of the System */
matrix a;	/* Given matrix */
vector p;	/* System vector keeping the coefficients */
matrix inv;   	/* a matrix to store the inverse of A (or NULL) */
  {
    matrix b;	/* matrix = ([A] - tr([A]) [I]) [A] */
    vector c;	/* vector used in the computation of [A] x [B] */
    counter i, j, k, l;
    
    b = (matrix) calloc(n * n, sizeof(real));
    c = (vector) calloc(n, sizeof(real));
    
    if ( (b == NULL) || (c == NULL) )
      return FAIL;
    
    /*
     *	The algorithm runs like this:
     *
     *	B_1 = A     	    	    and 	p_1 = Tr[B_1]
     *
     *  B_2 = A (B_1 - p_1 I)	    and     	p_2 = 1 / 2 Tr[B_2]
     *
     *	. . .
     *
     *  B_k = A (B_k-1 - p_k-1 I)   and     	p_k = 1 / k Tr[B_k]
     *
     *	. . .
     *
     *	B_n = A (B_n-1 - p_n-1 I)   and     	p_n = 1 / n Tr[B_n]
     *
     * In addition the inverse matrix is given by
     *
     *  A^-1 = (1 / p_n) (B_n-1 - p_n-1 I)
     */
    
    /* Process */
    /*	Start with B = A */
    j = n * n;
    for (i = 0; i < j; i++)
      b[i] = a[i];
      
    for (i = 0; i < n; i++)
      {
#ifdef EVV_DEBUG
	printf("\nB%d\n", i+1);
	for (j = 0; j < n; j++) {
	    for (k = 0; k < n; k++) 
	    	printf("%+3.3g ", b[j * n + k]);
	    printf("\n");
	}
#endif	
      	/* set p = (1/i) Trace(B)  --where trace(B) = b_1,1 + b_2,2 + ... + b_n,n */
	p[i] = 0.0;
        for (j = 0; j < n; j++)
          p[i] += b[(j * n) + j];
        p[i] /= (i + 1);
        /* resulting p[i] */
#ifdef EVV_DEBUG
    	printf("\np%d\n%+3.3g", i+1, p[i]);
#endif


	/* compute B_i = B_i - (p_i * I) */
	for (j = 0; j < n; j++)
          b[(j * n) + j] -= p[i];

#ifdef EVV_DEBUG
	printf("\nB%d - p%d I\n", i+1, i+1);
	for (j = 0; j < n; j++) {
	    for (k = 0; k < n; k++) 
	    	printf("%+3.3g ", b[j * n + k]);
	    printf("\n");
	}
#endif
    	/* if inverse matrix requested, save but-last B - p I */
	if ((i == n - 2) && (inv != NULL)) {
	    k = n * n;
	    for (j = 0; j < k; j++)
	    	inv[j] = b[j];
#ifdef EVV_DEBUG
	printf("\nA^-1%d - p%d I\n", i+1, i+1);
	for (j = 0; j < n; j++) {
	    for (k = 0; k < n; k++) 
	    	printf("%+3.3g ", inv[j * n + k]);
	    printf("\n");
	}
#endif
	}
	
    	/* we don't need to compute the last B */
    	if (i < (n - 1)) {
    	/* compute B_i+1 = A * (B_i - p_i * I) */
    	/* using square matrix multiplication */
        for (j = 0; j < n; j++)
          {
            for (k = 0; k < n; k++)
              {
                c[k] = b[(k * n) + j];
                b[(k * n) + j] = 0.0;
              }
            for (k = 0; k < n; k++)
              for (l = 0; l < n; l++)
                b[(k * n) + j] += a[(k * n) + l] * c[l];
          }
	}

      }
    free(b);
    free(c);
#ifdef EVV_INVERT_LF_SIGN
    for (i = 0; i < n ; i++)
    	p[i] = (-p[i]);
#endif

    if (inv != NULL) {
    	/* Using Leverrier-Faddev it is very easy to further
	compute the inverse matrix. If the user gave us some
	space to store it, we'll do so */
	/*
	 * The inverse of the matrix can easily be computed
	 * as
	 *
	 *  A^-1 = (1 / p_n) * (B_n-1 - p_n-1 * I)
	 */
	 j = n * n;
	 /* B already contains (B_n-1 - p_n-1 * I) */
	 /* So we only need to divide by p_n --which corresponds 
	  * to p[n-1] since we use zero-offset arrays */
	 for (i = 0; i < j; i++)
	    inv[i] = inv[i] / p[n-1];
    }

    return SUCCESS;
  }

/*
 * Find maximum eigenvalue and associated eigenvector of a square matrix
 *
 * Takes
 *  	Matrix size (matrix is square n x n)
 *  	Matrix
 *  	Initial eigenvector guess
 *  	pointer to store maximum eigenvalue
 *  	maximum allowed error
 *  	maximum number of iterations
 */
public status EVVSuccApprox(n, a, x, avmax, cerr, nmi)
Sys_ord n;	/* System order */
matrix a;	/* Data matrix  (n x n)*/
vector x;	/* Eigenvector. MUST HAVE AN INITIAL VALUE */
real *avmax;	/* Maximum eigenvalue */
real cerr;	/* Max acceptable error */
iter_no nmi;	/* Max number of iterations */
  {
    counter i, j;
    iter_no iter;   /* iteration counter */
    vector c;	    /* vector product of a . x */
    real temp,	    /* temporal variable */
    	 maxi;	    /* max value of c in each iteration */
    
    c = (vector) calloc(n, sizeof(real));
    if (c == NULL) return FAIL;
    
    iter = 0;
    maxi = *avmax;
    do
      {
        iter++;
        for (i = 0; i < n; i++)
          for (j = 0; j < n; j++)
            c[i] += a[i * n + j] * x[j];
        temp = maxi;
        maxi = 0.0;
        for (i = 0; i < n; i++)
          if (fabs(c[i]) > fabs(maxi))
            maxi = c[i];
        for (i = 0; i < n; i++)
          {
            x[i] = c[i] / maxi;
            c[i] = 0.0;
          }
      }
    while ((fabs(temp - maxi) > cerr) && (iter != nmi));
    
    *avmax = maxi;
    free(c);
    return (iter == nmi)? FAIL : SUCCESS;
  }


#ifdef NOTDEF
EVJacobi(n, a, v)
SysOrder_t n;	    /* Order of matrices */
matrix a;   	    /* matrix */
matrix v;   	    /* eigenvectors */
{
    int count, i, j, k, limit, skipped;
    real c, p, q, s, t;
    char ch;
    boolean oki, okj, rotn;
    
    count = 0;
    limit = 30;     /* arbitrarily chosen after Eberlein */
    skipped = 0;    /* no rotations skipped so far */
#define EPSILON 0.00000001
    
    while ((count <= limit) && (skipped < (n * (n -1)) / 2)) {
    	/* this is a check to avoid infinite execution */
	count++;    	/* sweep count */
	skipped = 0;	/* rotations skipped during the sweep */
	for (i = 0; i < n-1; i++) {
	    for (j = i+1; j < n; j++) {
	    	rotn = TRUE;	/* carry out a rotation unless unneeded */
		P = 0.5 * a[i * n + j] + a[j * n + i];
		q = a[i * n + i] * a [j * n + j];
		t = sqrt(4.0 * p * p + q * q);
		if (t < EPSILON) {
		    /* if t == 0 no rotation needed */
		    rotn = FALSE;
		}
		else {
		    if (q >= 0.0) {
		    	/* rotation for eigenvalue approximations in order */
		    	if (rotn) {
			    c = sqrt((t + q) / (2.0 * t));
			    s = p / (t * c);
			}
		    } else {
		    	/* q < 0.0 always rotate */
			rotn = TRUE;
			s = sqrt((t - q) / (2.0 * t);
			if ( p < 0.0) s = -s;
			c = p / (t * s);
		    }
		    if (fabs(s) < EPSILON)
		    	rotn = FALSE;
		}
	    }
	    if (rotn) {
	    	for (k = 0; k < n; k++) {
		    q = a[i *n + k];
		    a[i * n + k] = c * q + s * a[j * n + k];
		    a[j * n + k] = -s * q + c * a[j * n + k];
		}
		for (k = 0; k < n; k++) {
		    a = a[k * n + i];
		    a[k * n + i] = c * q + s * a[k * n + j];
		    a[k * n + j] = -s * q + c * a[k * n + j];
		    /* eigenvectors */
		    q = v[k * n + i];
		    v[k * n + i] = c * q + s * v[k * n + j];
		    v[k * n + j] = -s * q + c * v[k * n + j];

		}
	    }
	    else skipped++;
	}
    }
}
#endif
