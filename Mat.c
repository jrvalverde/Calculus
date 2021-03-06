/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*	MODULE MATRIX.C
*
*	DESCRIPTION:
*	Generic implementation of bidimensional matrix
*   calculus. It contains functions to allocate memory
*   space and, by this time, elementary matrix calculus.
*
*	CONTENTS:
*	typedef double *matrix;
*
*	double *matalloc(size)
*	    int size;
*	matrix mat2alloc(m, n)
*	    int, m, n;
*	minit(mat, m, n, value)
*	    matrix mat;
*	    int m, n;
*	    double value;
*	massign(mat1, mat2, m, n)
*	    matrix mat1, *mat2;
*	    int, m, n;
*	mtransp(trn, mat, m, n)
*	    matrix trn, mat;
*	    int m, n;
*	mident(mat, n)
*	    matrix mat;
*	    int n;
*	msum(sum, mat2, mat2, m, n)
*	    matrix sum, mat1, mat2;
*	    int m, n;
*	msubst(res, min, sust, m, n)
*	    matrix res, min, sus;
*	    int m, n;
*	mmult(prod, fact1, fact2, m, n, p)
*	    matrix prod, fact1, fact2;
*	mscprod(prod, mat, m, n, lambda)
*	    matrix prod, mat;
*	    int m, n;
*	    double lambda;
*	int minv(inv, mat, n)
*	    matrix inv, mat;
*	    int n;
*
*	RETURNS:
*	    The two allocation functions, matalloc()
*	and mat2alloc() return a pointer to the allocated
*	memory (or null on failure). See matrix type.
*	    The other functions do not return anything,
*	except for minv(), that returns FAIL if it can't
*	allocate the memory needed for the auxiliary
*	matrix, M_SINGULAR if the matrix to invert is
*	singular and SUCCESS if all goes well.
*
*	CAVEATS:
*	    You must be cautious in sending all the required
*	parameters, or else, you can get a pretty hanging of
*	the computer. (I discovered this the hard way when
*	debugging).
*	    A matrix of dimensions M x N must not be declared
*		double matrix[m-1][n-1];
*	but
*		double matrix[m][n];
*	although you handle it as  a normal array, referencing
*	its elements from 0 to m-1 and n-1.
*	    If the matrixes are declared as pointers to double,
*	allocating space with the provided functions all goes
*	softly.
*	    You must remember to check the return values to
*	avoid gross errors, specially to check that the
*	matrix passed to minv() could be singular.
*
*	FILES:
*		STDIO.H
*		PORTABLE.H
*		MATH.H
*		UNIX.H
*		STORAGE.H
*
*	NOTES:
*	    Estoy hasta los cojones.
*	    Today (1 - oct - 1988) is my father's birthday
*	and here I am like a dumb loser.
*	    I intend to add functions for the solution
*	of systems of equations, manipulation of eigenvalues
*	and eigenvectors.
*
*	SEE ALSO:
*	    MatTest.C:	    Is an auxiliary module with
*	the test code of these functions used to debug this
*	module (last time everything went OK).
*
*	DESIGNED BY:
*	    Jos� Ram�n Valverde Carrillo.
*
*	LAST MODIFICATION:
*	    1 - october - 1988.
*
*	COPYRIGHT:
*	    � YoEgo.	Since I have no cash, I can't
*	register this (nor do I believe I should). So
*	this module is left in the PUBLIC DOMAIN.
*	    It is furthermore forbidden its use for
*	commercial purposes unless I get a share on
*	the profits.
*	    I say.
*						YoEgo.
*
*
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "portable.h"

typedef double *matrix;

#define M_SINGULAR  (-2)

#define matfree(x)	free(x)

/*----------------------------------------------------------------*/

public double *matalloc(size)
int size;	    /* Size of the matrix: m x n in 2D	*/

/*
    Allocates memory space for a general matrix of size
    elements
*/
  {
    auto double *theMatrix;

    theMatrix = (double *) calloc(size, sizeof(double));
    
    return (theMatrix);
    
  }

/*----------------------------------------------------------------*/

public double *mat2alloc(m, n)
int m, n;	    /* dimensions of the matrix to allocate  */
/*
    Allocates memory space for a bidimensional matrix
    of M x N elements
*/
  {
    auto double *theMatrix;

    theMatrix = (double *) calloc(m * n, sizeof(double));
    
    return (theMatrix);
    
  }


/*----------------------------------------------------------------*/

public minit(theMatrix, m, n, value)
double *theMatrix;		/* matrix to initialize	*/
int m, n;			/* dimensions		*/
double value;			/* fill value		*/
/*
    Sets the [m] x [n] specified elements to value
*/
  {
    register int i, j;
    
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
        theMatrix[(i * n) + j] = value;
  }

/*----------------------------------------------------------------*/

public massign(mat1, mat2, m, n)
double *mat1, *mat2;	    /* receptor and donor matrices    */
int m, n;		    /* dimensions			    */
/*
  MAT mat1 = mat2; transfer ([m] x [n]) elements from 2 to 1.
*/
  {
    register int i, j;
    
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
        mat1[(i * n) + j] = mat2[(i * n) + j];
  }
  
/*----------------------------------------------------------------*/

public mtransp(trn, theMatrix, m, n)
double	*trn, 		/* transposed matrix	*/
	*theMatrix;	/* matrix to transpose	*/
int m, n;		/* dimensions		*/
/*
    MAT trn = TRN(theMatrix); computes the transposed matrix 
    of theMatrix[m] x [n] and stores it in en trn[n] x [m]
*/
  {
    register int i, j;
    
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
        trn[(j * m) + i] = theMatrix[(i * n) + j];
  }

/*----------------------------------------------------------------*/

public mident(theMatrix, n)
double *theMatrix;	    /* matrix to initialize as an identity	*/
int n;			    /* dimension (is a square one)		*/
/*
    Compute the identity matrix of dimension [n] x [n]:
    Every element ij / i <> j is 0, every ij / i = j is 1.
*/
  {
    register int i, j;
    
    for (i = 0; i < n; i++)
      {
        for (j = 0; j < n; j++)
          theMatrix[(i * n) + j] = 0.0;
        theMatrix[(i * n) + i] = 1.0;
      }
  }
  
/*----------------------------------------------------------------*/

public msum(sum, mat1, mat2, m, n)
double *sum, *mat1, *mat2;	/* sum y sumands	*/
int m, n;			/* dimensions		*/
/*
    MAT sum = mat1 + mat2;
    computes the sum of matrices mat1 and mat2, of dimensions m x n,
    and stores it in sum
*/
  {
    register int i, j;
    
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
        sum[(i* n) + j] = (mat1[(i * n) + j]) + (mat2[(i * n) + j]);
  }

/*----------------------------------------------------------------*/

public msubst(res, mat1, mat2, m, n)
double *res, *mat1, *mat2;	/* rest, minuend y sustraend	    */
int m, n;			/* their dimensions		    */
/*
    MAT res = mat1 - mat2;
    Substracts mat2 from mat1 and stores result in res.
*/
  {
    register int i, j;
    
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
        res[(i * n) +j] = mat1[(i * n) +j] - mat2[(i * n) + j];
  }     

/*----------------------------------------------------------------*/

public mmult(prod, fact1, fact2, m, n, p)
double *prod, *fact1, *fact2;	    /* product and factors		*/
int m, n, p;			    /* abbreviature of their dimensions	*/
/*
    MAT prod[m][p] = fact1[m][n] * fact2[n][p]
*/
  {
    register int i, j, k;
    auto double temp;
    
    for (i = 0; i < m; i++)
      for (j = 0; j < p; j++)
        {
          temp = 0.0;
          for (k = 0; k < n; k++)
            temp += (fact1[(i * n) + k] * fact2[(k * p) + j]);
          prod[(i * p) + j] = temp;
        }
  }

/*----------------------------------------------------------------*/

public mscprod(scprod, theMatrix, m, n, lambda)
double *scprod, *theMatrix;	/* producto y padre de la criatura  */
int m, n;			/* dimensiones			    */
double lambda;			/* el tan ansiado escalar	    */
  {
    register int i, j;
    
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
        scprod[(i * n) + j] = theMatrix[(i * n) + j] * lambda;
  }
  
/*----------------------------------------------------------------*/

public int minv(b, c, n)
double *b, *c;		/* matrices inverse and to invert	*/
int n;			/* both are square matrices n x n	*/
/*
    MAT B = INV ( C );
    Computes the inverted matrix of C, which is a square matrix n x n
    and stores the result on B.
    The algorithm employed is a modification of Gauss' methos.
    To avoid changing matrix C we get first a working copy that
    we store in matrix A, whose space is reserved ex-profeso and
    is freed upon termination.
*/
  {
    int i, j, k, l;		/* indexes		    */
    double temp;		/* auxiliary variable	    */
    double *a;			/* backup matrix	    */
    bool   singular = FALSE;	/* TRUE if c is singular    */
   
    a = mat2alloc(n, n);    /* make up a */
    if (a == NULL)	    /* (if we can) */
      return FAIL;
    massign(a, c, n, n);    /* copy c into a to avoid changing it */
    
    /* make b equal the identity matrix */
    mident(b, n);
    
    for (j = 0; j < n; j++)
      {
        for (i = j; i < n; i++)
          if (a[i * n + j] != 0)    /* is singular ?	*/
            goto __Compute__;	    /* NO		*/
        singular = TRUE;	    /* SI		*/
        goto __end__;
__Compute__:
        for (k = 0; k < n; k++)
          {
            /* Permutar a[j][k] por a[i][k]*/
            
            temp = a[j * n + k];
            a[j * n + k] = a[i * n + k];
            a[i * n + k] = temp;
            
            /* Permutar inv[j][k] por inv[i][k] */
            
            temp = b[j * n + k];
            b[j * n + k] = b[i * n + k];
            b[i * n + k] = temp;
          }	/* K */
           
        temp = 1 / a[j * n + j];
        for (k = 0; k < n; k++)
          {
            a[j * n + k] *= temp;
            b[j * n + k] *= temp;
          }	/* K */
       
        for (l = 0; l < n; l++)
          {
            if (l != j)
              {
                temp = -a[l * n + j];
                for (k = 0; k < n; k++)
                  {
                    a[l * n + k] += 
                            (temp * a[j * n + k]);
                    
                    b[l * n + k] += (temp * b[j * n + k]);
                  } 	/* K */
              }	    /* IF */
          }	/* L */  
      }	    /* J */
__end__:     
    /* b = inverse of c 		    */
    /* a is no longer needed. kiss good bye.  */
    matfree(a);	    
    if (singular == TRUE)
      return (M_SINGULAR);	/* so sorry	*/
    else
      return (SUCCESS);		/* happy end	*/
  }

/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*							*
*		END OF MODULE MATRIX.C			*
*		    (T'was time!)			*
*							*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/
