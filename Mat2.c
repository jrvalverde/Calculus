#include <stdio.h>
#include "portable.h"
#include <math.h>
#include "Mat.h"

nop()
  {
  return;
  }

GaussShort(x, a, b, n)	    /* NOT VERIFIED	*/
matrix	x,  /* solutions				*/
	a,  /* matrix of equation coefficients		*/
	b;  /* results vector				*/
int	n;  /* dimensions: a[n x n], x[n], b[n]		*/
/*
    Computes the solution of the system of equations A . X = B
    Aborts if |A| = 0 since then the solution is not unique
    and in the elimination process we get a divide by zero.
    Also if in any step one of the diagonals is too close to zero
    due to overload of the accumulator: too big an exponent.
    Sometimes it can yield erroneous results:
	an unrefined algorithm for valiants.
*/
  {
    double m, d;
    int i, j, k;
    
    /* Make n - 1 elimination steps.
    	Current step is number k.			*/
    for (k = 0; k < (n - 1); k++)
      {
        /* Compute multiplier Mik to eliminate Aik	*/
        for (i = k + 1; i < n; i++)
          {
            m = -a[i * n + k] / a[k * n + k];
            /* make operations over the rows of A	*/
            for (j = k + 1; j < n; j++)
              a[i * n + j] += (m * a[k * n + j]);
            /* make operations over the rows of B	*/
            b[i] = b[i] + m * b[k];
          }
      }
    /* Compute x(n)..x(1) by backward substitution	*/
    x[n] = b[n-1] / a[n * n - 1];
    for (i = n - 2; i >= 0; i--)
      {
        d = b[i];
        for (j = i + 1; j < n; j++)
          d -= a[i * n + j] * x[j];
        x[i] = d / a[i * n + j];
      }
  }

/*------------------------------------------------------------------*/

Gauss(x, r, det, a, b, n)	/* NOT VERIFIED	    */
matrix	x,	    /* solutions			    */
	r;	    /* residual vector			    */
double	det;	    /* Value of determinant of A	    */
matrix	a,	    /* coefficients of the equations	    */
	b;	    /* vector of results		    */
int	n;	    /* Dimension: A[n x n] B[n] X[n] R[n]   */
/*
    Solution of systems of equations by Gauss' method.
    Computes the residual vector to verify the result
	r = b - ax
    If r has small values then the equation is satisfied.
    To compute the conditioning of A we also compute its
    determinant |A| in the Gauss elimination: if it is
    small it is badly conditioned.
	If r[i] -> 0	=> well solved
	If |A|  -> 0	=> bad precision
    The former caveats are still valid	!
*/
  {
    matrix  a2,	    /* dimension n x n 	*/
    	    b2;	    /* dimension n	*/
    double d, m, w;
    int i, j, k;
    
    /* Get copies of the original A, B in A2, B2   */
    a2 = mat2alloc(n, n);
    b2 = matalloc(n);
    massign(a2, a, n, n);
    massign(b2, b, n, 1);
    
    /* As in GaussShort but accumulating the determinant of A in det  */
    det = 1.0;
    for (k = 0; k < (n - 1); k++)
      for(i = k + 1; i < n; i++)
        {
          m = - a2[i * n + k] / a2[k * n + k];
          for (j = k + 1; j < n; j++)
            a2[i * n + j] += (m * a2[k * n + j]);
          b[i] += (m * b2[k]);
        }
    /* Compute det = det(A) */
    for (i = 0; i < n; i++)
      det *= a2[i * n + i];
      
    /* Compute x(n)..x(1)  */
    x[n] = b2[n-1] / a2[n * n - 1];
    for (i = n - 2; i >= 0; i--)
      {
        d = b2[i];
        for (j = i + 1; j < n; j++)
          d -= a2[i * n + j] * x[j];
        x[i] = d / a2[i * n + i];
      }
    /* Solution is now in X  */
    
    /* Compute residual vector	r = b - ax  */
    for (i = 0; i < n; i++)
      {
        w = b[i];   /* the original value    */
        for (j = 0; j < n; j++)
          w -= a[i * n + j] * x[j];
        r[i] = w;
      }
  }
  
GaussPlus(x, r, det, a, b, n)
matrix	x,	    /* Solutions			    */
	r;	    /* Residual vector			    */
double	det;	    /* Value of determinant of A	    */
matrix	a,	    /* Coefficients of the equations	    */
	b;	    /* Vector of results		    */
int	n;	    /* Dimension: A[n x n] B[n] X[n] R[n]   */
		    /* N = nunber of equations.		    */
/*
    Enhanced solution of equations by GAUSS method
    To avoid roundoff errors and divide by zero abort
    we change the pivot: if all the multipliers have
    magnitud < 1.0 we control roundoff.
    Computes the residual value B - AX and the determinant of A.
    Solves the same system with new independent terminis without
    repeating all the calculi.
    Stores the multipliers in the lower triangle of A.
*/
  {
    matrix  a2,
    	    b2,
    	    l;
    int i, j, k, l1;
    double c, d, m, w;
    bool otro;
    
    a2 = mat2alloc(n, n);
    b2 = mat2alloc(n, 1);
    l = mat2alloc(n, 1);
    
    massign(a2, a, n, n);
    massign(b2, b, n, 1);
    
    /*	Make N-1 elimination steps.
    	Accumulate det(a) in det.	    */
    /*	Current step is k		    */
    for (k = 0; k < (n - 1); k++)
      {
        /*  Look for bigger element in column k.
            record the no. of the row l(k) of that element. */
        d = 0;
        for (i = k; i < n; i++)
          {
            c = fabs(a2[i * n + k]);
            if (d < c)
              {
                d = c;
                l[k] = i;
              }
          }
        l1 = l[k];
        if (l1 != k)
          {
            /* change sign of the determinant	*/
            det = -det;
            /* swap rows k and l[k]		*/
            for (j = k; j < n; j++)
              {
                d = a2[k * n + j];
                a2[k * n + j] = a2[l1 * n + j];
                a2[l1 * n + j] = d;
              }
           } 
         /* compute multiplier m[i,k] to eliminate a[i,k]	*/
         for (i = k + 1; i < n; i++)
           {
             m = -a2[i * n + k] / a2[k * n + k];
             /* store m in a[i,k]		*/
             a2[i * n + k] = m;
             /* operate on rows of a		*/
             for (j = k + 1; j < n; j++)
               a2[i * n + j] += m * a2[k * n + j];
           }
       }
       
     /* compute determinant det of a	*/
     for (i = 0; i < n; i++)
       det *= a2[i * n + i];
     
     /* swap rows k and l[k] of b	*/
__Solve__:
     for (k = 0; k < (n - 1); k++)
       {
         l1 = l[k];
         if (l1 != k)
           {
             d = b[k];
             b2[k] = b2[l1];
             b2[l1] = d;
           }
         /* operate over rows of b */
         for (i = k + 1; i < n; i++)
           {
             m = a2[i * n + k];
             b2[i] += m * b2[k];
           }
       }
     
     /* compute x[n]..x[1] by backward substitution	*/
     x[n] = b2[n - 1] / a2[(n * n) - 1];
     for (i = n - 2; i == 0; i--)
       {
         d = b[i];
         for (j = i + 1; j < n; j++)
           d -= a2[i * n + j] * x[j];
         x[i] = d / a[i * n + i];
       }
     
     /* compute residual vector r = b - ax;
     for (i = 0; i < n; i++)
       {
         w = b[i];
         for (j = 0; j < n; j++)
           w -= a[i * n + j] * x[j];
         r[i] = w;
       }
     
    /* solve the system for a new b if needed    */
    /* !!!!!!!!!!!!!!!!    MAKE    !!!!!!!!!!!!!!!!!    */
    otro = FALSE;
    if (otro)
      nop();
    goto __Solve__;
  }

/*--------------------------------------------------------------*/
