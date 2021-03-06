/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*	MODULE
*	    PDE.C
*
*	Functions to solve partial differential equations.
*
*	CONTENTS:
*
*	USES:
*	    STDIO.H
*	    PORTABLE.H
*	    MATH.H
*
*	DESIGNED BY:
*	    Jos� Ram�n Valverde Carrillo
*
*	LAST MODIFICATION:
*	    19 - mar - 1989
*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/


#include <stdio.h>

#include <portable.h>

#include <math.h>

typedef double real;	    /* This is to be able to control precision */

typedef unsigned short iter_no;

typedef unsigned Sys_ord;

typedef unsigned count;

typedef unsigned counter;

typedef real *matrix;

typedef real *vector;

#define odd(x)		((x) & 1)
	/* sets all bits but last to zero. If the number was even
	we get zero and return FALSE, otherwise we get 1 and TRUE */

#define even(x)		(! odd(x))

/*--------------------------------------------------*
|						    |
|	    SOLUTION OF PARTIAL DIFFERENTIAL	    |
|		    EQUATIONS			    |
|						    |
*--------------------------------------------------*/

/*
    Finite differences method for eliptic equations
*/

public status PDEElipFinDif(vals, nf, nc, nmi, cerr)
matrix vals;	    /* Matrix of values of f in each point of the mesh */
Sys_ord nf;	    /* Number of rows */
Sys_ord nc;	    /* Number of columns */
iter_no nmi;	    /* Maximum number of iterations */
real cerr;	    /* Allowed error */
/*  Upon exit Sys contains the computed values */
  {
    iter_no iter;   /* Number of this iteration */
    flag done;	    /* Relative error in this iteration */
    real temp;	    /* Temporal variable */
    counter i, j;   /* Counters */
    
    iter = 0;
    do
      {
        iter++;
        done = TRUE;
        for (i = 0; i < (nf - 1); i++)
          for (j = 1; j < (nc - 1); j++)
            {
              temp = vals[((i - 1) * nf) + j] +
              	     vals[((i + 1) * nf) + j] +
              	     vals[(i * nf) + (j - 1)] +
              	     vals[(i * nf) + (j + 1)];
              temp /= 4.0;
              if (fabs(vals[(i * nf) + j] - temp) > cerr)
                done = !done;
              vals[(i * nf) + j] = temp;
            }
      }
    /* New iteration */
    while ((!done) && (iter != nmi));
    
    return (done)? SUCCESS : FAIL;
  }

/*
    Finite Differences in parabolic equations
*/

public status PDEParFinDif(vals, nf, nc, tiv, tm, cerr)
matrix vals;	    /* Matrix of vslues of f */
Sys_ord nf;	    /* Number of rows */
Sys_ord nc;	    /* Number of columns */
real tiv;	    /* Initial temperature of bar */
real tm;	    /* Medium temperature */
real cerr;	    /* Uniformity limit */
  {
    vector aux;	    /* Vector of values of f at time i */
    flag done;	    /* Error flag */
    counter i, j, k;
    
    /* Allocate space for aux */
    aux = (vector) calloc(nc + 1, sizeof(real));
    if (aux == NULL)
      return FAIL;
    
    i = 0;
    vals[0] = aux[0] = tm;
    for (j = 1; j < nc; j++)
      vals[j] = aux[j] = tiv;
    /* We have assigned values to the first row */
    
    do
      {
        for (k = 1; k < nc; k++)
          aux[k] = vals[(i * nf) + k];
        i++;
        aux[nc + 1] = aux[nc - 1];
        done = TRUE;
        vals[i * nf] = aux[0];
        for (j = 1; j < nc; j++)
          {
            vals[(i * nf) + j] = 0.5 * (aux[j - 1] + aux[j + 1]);
            if ( fabs(vals[(i * nf) + j] - vals[(i * nf) + j - 1]) > cerr)
              done = !done;
          }
      }
    while ( (!done) && (i < nf) );
    
    return (done)? SUCCESS : FAIL;
  }


/*
    Finite Differences for hiperbolic equations
*/

public status PDEHipFinDif(vals, ne, ti)
matrix vals;	/* Matrix of values */
counter ne,	/* Number of elements */
	ti;	/* Time to lapse */
 {
   /* vals[0][] contains the initial values */
   /* vals[i][j] will have the values of each point at time t = i + 1 */
   counter i, j;
      
   for (i = 1; i < (ne - 1); i++)
     {
       vals[ne + i] = (vals[i - 1] + vals[i + 1]) / 2;
     }
   /* We have got the second vector now */
   
   for (i = 2; i < (ti - 1); i++)
     {
       for (j = 1; j < (ne - 1); j++)
         {
           vals[(ne * i) + j] = vals[(ne * (i - 1)) + j - 1]
           		      + vals[(ne * (i - 1)) + j + 1]
           		      - vals[(ne * (i - 2)) + j];
         }
       /* We have got the vector position for t = i + 1 */
       
      }
  
    return SUCCESS;
  }

