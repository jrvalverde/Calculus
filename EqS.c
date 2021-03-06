/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*	MODULE
*	    EqSys.C
*	Procedures to solve equation systems, both through direct
*   and iterative methods.
*
*	CONTENTS:
*	    int EqSGauss(n, a, x)
*	    int EqS_Gauss(n, system)
*	    int EqSGauss_Jordan(n, a, x)
*	    int EqSCholesky(n, a, b, x)
*	    int EqSJacobi(n, a, x, cerr, nmi)
*	    int EqSGauss_Seidel(n, eqns, solns, precis, nmi)
*	    int EqSGaussSeidel(n, a, x, cerr, nmi)
*
*	DESIGNED BY:
*	    Jos� Ram�n Valverde Carrillo
*
*	LAST MODIFICATION:
*	    14 - dic - 1988 (implementation on IBM PC-AT)
*	    16 - dic - 1988 (ported to Mac)
*	    17 - dic - 1988 (library splitted and corrected)
*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

#include <stdio.h>

#include <portable.h>

#include <math.h>

typedef double real;	    /* to control precision */

typedef unsigned short iter_no;

typedef unsigned Sys_ord;

typedef unsigned counter;

typedef real *matrix;

typedef real *vector;

#define odd(x)		((x) & 1)
#define even(x)		(! odd(x))

/*--------------------------------------------------*
|						    |
|	SOLUTION OF SYSTEMS OF EQUATIONS	    |
|						    |
*---------------------------------------------------*/

/* Gauss' method */

public status EqSGauss(n, a, x)
Sys_ord n;	/* order of the system */
matrix a;	/* amplied coefficient matrix */
vector x;	/* unknowns vector */
  {
    counter up, i, j, k;
    real maxi, tmp;
    
    /* Triangulation module */
    
    /* Look for pivot row */
    for (i = 0; i < (n - 1); i++)
      {
        maxi = 0.0;
        for (j = 0; j < n; j++)
          {
            tmp = fabs(a[j * n + i]);
            if (tmp > maxi)
              {
                maxi = tmp;
                up = j;
              }
          }
        /* Test singular matrix */
        if (maxi == 0.0)
          return FAIL;
        /* swap rows */
        for (j = 0; j <= n; j++)
          {
            maxi = a[up * n + j];
            a[up * n + j] = a[i * n + j];
            a[i * n + j] = maxi;
          }
        for (j = i + 1; j < n; j++)
          for (k = i + 1; k < (n + 1); k++)
            a[j * n + k] -= a[i * n + k] * a[j * n  + i] / a[i * n + i];
      }
    
    /* Backwards solution/substitution module */
    
    for (i = n - 1; i >= 0; i--)
      {
        maxi = 0.0;
        if (i != n)
          for (j = 0; j < n; j++)
            maxi += a[i * n + j] * x[j];
        x[i] = (a[(i * n) + (n + 1)] - maxi) / a[i * n + i];
      }

    return SUCCESS;
  }

public EqS_Gauss(n, system)
Sys_ord n;
matrix system;		/* system[n, n + 1] */
/*
Given the no. of equation in the system, and the system, an
amplied matrix in the form
	| a11, a12, ... , a1n | b1 |
	| a21, a22, ... , a2n | b2 |
	| ...  ...  ...   ... | .. |
	| an1, an2, ... , ann | bn |
this procedure finds the solution matrix corresponding
	|  1 ,  0 , ... ,  0  | x1 |
	|  0 ,  1 , ... ,  0  | x2 |
	| ...  ...  ...   ... | .. |
	|  0,   0 , ... ,  0  | xn |
where ax = b		(dim(a) = n * (n + 1))
*/
  {
    real pivot;	    /* pivot element in a row */
    counter m,	    /* number of columns = n + 1  */
        i, j, k;    /* iteration and subindex counters */
        
    /* Gaussian elimination */
    /*	    Forward elimination */
    m = n + 1;
    for (i = 0; i < n; i++)
      {
        /* divide each element in the row by pivot */
        pivot = system[i * n + i];
        for (j = 0; j < m; j++)
          system[i * n + j] /= pivot;
        
        /* rests a multiple of the row to each row */
        for (j = i + 1; j < n; j++)
          {
            pivot = system[j * n + i];
            for (k = 0; k < m; k++)
              system[j * n + k] -= pivot * system[i * n + k];
            system[j * n + i] = 0.0;
          }
      }
    
    /* Backwards substitution */
    for (i = n - 1; i > 0; i--)
      {
        for (j = 0; j < (i - 1); j++)
          {
            system[j * n + m] -= system[j * n + i] * system[i * n + m];
            system[j * n + i] = 0.0;
          }
      }
  }

public status EqSGauss_Jordan(n, a, x)
Sys_ord n;  /* order of the system */
matriz a;   /* amplied coefficient matrix */
vector x;   /* unknowns vector */
  {
    flag *ind;	/* vector of flags of use of the row as a pivot
    		ind[i] = 1 ==> row i has been used
    		ind[i] = 0 ==> row i has not been used	    */
    counter *ord;   /* vector storing the order of use of row i as a pivot */
    real temp;	/* temporal variable */
    counter np,	/* number of pivot row */
    	i, j, k;    /* loop control variables */
    
    /* allocate space for auxiliary vectors */
    ind = (flag *) calloc(n, sizeof(flag));
    ord = (counter *) calloc(n, sizeof(counter));
    
    /* inicialization */
    for (i = 0; i < n; i++)
      {
        ind[i] = 0;
        ord[i] = 0;
      }
    
    /* Diagonalization module */
    /* Look for pivot row */
    for (i = 0; i < n; i++)
      {
        temp = 0.0;
        for (j = 0; j < n; j++)
          if (ind[j] != 1)
            if (fabs(a[j * n + i]) > fabs(temp))
              {
                temp = a[j * n + i];
                np = j;
              }
        
        /* test singular matrix */
        if (temp = 0.0)
          {
            free(ind);
            free(ord);
            return FAIL;
          }
        ind[np] = 1;
        ord[i] = np;
        for (j = 0; j < n; j++)
          if (j != np)
            for (k = i + 1; k < (n + 1); k++)
              a[j * n + k] -= a[np * n + k] * a[j * n + i] / temp;
      }
    
    /* Solution module */
    for (i = 0; i < n; i++)
      {
        np = ord[i];
        x[i] = a[(np * n) + (np + 1)] / a[np * n + i];
      }
    
    free(ind);
    free(ord);
    return SUCCESS;
  }

/* Cholesky's method */

public EqSCholesky(n, a, b, x)
Sys_ord n;  /* system order */
matrix a;   /* coefficient matrix of order n x n */
vector b;   /* vector of independent termini */
vector x;   /* vector of unknowns */
  {
    real sum;	/* auxiliary variable (summatories) */
    real temp;	/* temporal variable */
    counter i, j, k;
    
    /* inicialization */
    for (i = 0; i < n; i++)
      x[i] = 0.0;
    sum = temp = 0.0;
    
    /* Decomposition module */
    for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
        {
          sum = a[i * n + j];
          for (k = 0; k < (i - 1); k++)
            sum -= a[k * n + i] * a[k * n + j];
          if (i != j)
            a[i * n + j] + sum * temp;
          else
            {
              temp = 1 / sqrt(sum);
              a[i * n + j] = temp;
            }
        }
      
    /* Solution module. Get x* */
    for (i = 0; i < n; i++)
      {
        sum = b[i];
        for (k = 0; k < (i - 1); k++)
          sum -= a[k * n + i] * x[k];
        x[i] = sum * a[i * n + i];
      }
    
    /* Solution module. Get x */
    for (i = n - 1; i >= 0; i--)
      {
        sum = x[i];
        for (k = i + 1; k < n; k++)
          sum -= a[i * n + k] * x[k];
        x[i] = sum * a[i * n + i];
      }
  }

/* Jacobi's iterative method */

public status EqSJacobi(n, a, x, cerr, nmi)
Sys_ord n;	/* system order */
matrix a;	/* amplied matrix of the system */
vector x;	/* solution. Iteration k of vector of unknowns */
real cerr;	/* allowed error */
iter_no nmi;	/* maximum number of iterations */
  {
    vector y;	    /* iteration k + 1 of unknowns vector */
    real temp,	    /* temporal variable */
    	 sum;	    /* summatorium of the algorithm */
    iter_no iter;   /* iteration counter */
    counter i, j;   /* loop indexes and counters */
    flag erit;	    /* error flag of the iteration */
    
    y = (vector) calloc(n, sizeof(real));
    
    /* Iteration module */
    iter = 0;
    for (;;)
      {
        iter++;
        erit = FALSE;
        for (i = 0; i < n; i++)
          {
            y[i] = a[(i * n) + n + 1];
            for (j = 1; j < n; j++)
              if (j != i)
                y[i] -= a[i * n + j] * x[j];
            y[i] /= a[i * n + i];
            if (fabs(y[i] - x[i]) > cerr)
              erit = TRUE;
          }
        if (erit == FALSE)
          break;
        else
          if (iter == nmi)
            {
              free(y);
              return FAIL;
            }
          else
            {
              for (i = 0; i < n; i++)
                x[i] = y[i];
            }
      }
    free(y);
    return SUCCESS;
  }

public status EqSGauss_Seidel(n, eqns, solns, precis, nmi)
Sys_ord n;
matrix eqns;
vector solns;
iter_no nmi;
/*
  Given the parameters n, eqns, solns and precis as described
  	n - no. of equations
  	eqns - matrix of equations, sized n x n
  	solns - solution vector of size n
  	precis - desired precision of the result
  	nmi - maximum number of iterations
  this procedure solves the set of linear simultaneous
  equations by the iterative method of Gauss - Seidel.
    Other used variables are
	sum - the sum of squares of the differences
	old - a previous approximation
	done - indicates iteration convergence
	i, j, k - iteration indexes and counters 
*/
  {
    real sum, old;
    counter i, j, k;
    bool done;
    
    /* Gauss - Seidel iteration */
    done = FALSE;
    i = 1;
    /* iterate up to maximum of nmi */
    while ((!done) && (i <= nmi))
      {
        /* compute a new approximation set */
        sum = 0.0;
        for (j = 0; j < n; j++)
          {
            old = solns[j];
            solns[j] = eqns[j * n + 0];
            for (k = 1; k < n; k++)
              {
                if (k != j)
                  solns[j] += eqns[j * n + k] * solns[k - 1];
                else
                  solns[j] += eqns[j * n + k] * solns[k];
              }
            sum += (solns[j] - old) * (solns[j] - old);
          }
        
        /* test precision */
        if (sum < precis)
          done = TRUE;
        else
          i++;
      }
    
    /* converges ? */
    if (done)
      return SUCCESS;
    else
      return FAIL;
  }

public status EqSGaussSeidel(n, a, x, cerr, nmi)
Sys_ord n;	    /* system order */
matrix a;	    /* amplied coefficient matrix */
vector x;	    /* unknowns vector */
real cerr;	    /* allowed error */
iter_no nmi;	    /* number of max. iterations */
  {
    real temp,	    /* temporal variable */
    	 sum;	    /* algorithm's sumatorium */
    bool erit;	    /* error flag */
    iter_no iter;   /* iteration counter */
    counter i, j;
    
    /* get an initial solution */
    for (i = 0; i < n; i++)
      x[i] = a[(i * n) + n + 1] / a[i * n + i];
    
    /* iterate */
    iter = 0;
    do
      {
        erit = TRUE;
        for (i = 0; i < n; i++)
          {
            sum = 0.0;
            for (j = 0; j < n; j++)
              if (j != i)
                sum += a[i * n + j] * x[j];
            temp = (a[(i * n) + n + 1] - sum) / a[i * n + i];
            if ((fabs(temp) - x[i]) > cerr)
              erit = FALSE;
            x[i] = temp;
          }
        if (erit == FALSE)
          break;
        else
          iter++;
      }
    while (iter <= nmi);
    
    if (erit == TRUE)
      return FAIL;
    else
      return SUCCESS;
  }

