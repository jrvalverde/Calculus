/**
 *  @file mat_ops.c
 *
 *  @brief 2D-matrix and elements basic operations.
 *
 *  This module implements basic operations on 2D matrices: assignment
 *  of one matrix contents to another, computing the transpose and the
 *  inverse.
 *
 *  @pre    portable.h
 *
 *  @pre    matrix.h
 *
 *  @see    matrix.h for a general introduction to matrices
 *
 *  @see    mat_housekeep.c to learn how to create/destroy matrices
 *	    before using this module.
 *
 *  @see    mat_init.c to learn how you can assign initial values to
 *	    a matrix or its elements.
 *
 *  @see    mat_arith.c to learn more about how to perform basic
 *	    matrix arithmetic.
 *
 *  @author Jos� Ram�n Valverde Carrillo    (jrvalverde@acm.org)
 *
 *  @version	3.0
 *
 *  @date   23 - february - 2004    v3.0
 *
 *  @date   11 - february - 2004    v2.0
 *
 *  @date    1 - october - 1988     Last modification of v1.0
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
 * $Id: mat_ops.c,v 1.1 2004/03/05 18:50:14 jr Exp $
 * $Log: mat_ops.c,v $
 * Revision 1.1  2004/03/05 18:50:14  jr
 * Initial revision
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <bits/nan.h>

#include "portable.h"
#include "matrix.h"



/*----------------------------------------------------------------*/
/*			MATRIX OPERATIONS			  */
/*----------------------------------------------------------------*/

/** @defgroup matrix_operations
 *  @{
 */

/**
 * @brief MAT dest = orig
 *
 * Transfer ([m] x [n]) elements from 2 to 1.
 *
 * Assigns the values of one matrix to another. Both matrices must have
 * the same dimensions, otherwise an error result will occur.
 *
 * If the origin matrix is NULL, then the destination matrix will be
 * set to the null matrix (i.e. all of its elements initialized to zero).
 *
 * @param dest  the destination matrix
 * @param orig  the origin matrix
 *
 * @return SUCCESS if all went well, an error code otherwise
 *
 */
public status mat_assign(matrix dest, matrix orig)
{
#ifndef MAT_OPTIMIZE
    register int i, j, maxrow, maxcol;
    real **val1, **val2;
#else
    register int i, size;
    real *element1;
    real *element2;

#endif

    if (orig == NULL)
	return mat_init(dest, 0.0);

#ifdef MAT_PARANOIA
    if ((dest->rows != orig->rows) || (dest->cols != orig->cols))
	return MAT_BOUNDSCHECK;
#endif

#ifndef MAT_OPTIMIZE
    maxrow = dest->rows;
    maxcol = dest->cols;
    val1 = dest->values;
    val2 = orig->values;
    for (i = 1; i <= maxrow; i++)
	for (j = 1; j <= maxcol; j++)
	    val1[i][j] = val2[i][j];
#else
    /* we know all values are stored contiguously after mat->values[1][1] */
    size = dest->rows * dest->cols;
    element1 = dest->values[1];
    element2 = orig->values[1];
    for (i = 1; i <= size; i++) {
	/* element1[i] = element2[i]; */
	*element1 = *element2;
	element1++;
	element2++;
    }
#endif
    return SUCCESS;
}

/**
 * @brief MAT trn = TRN(mat)
 *	  Compute the transpose of a matrix
 *
 * The transpose A<sup>T</sup> of a matrix A is
 * A' = A<sup>T</sup> / row A' = col A =>		 
 *
 *  (A)<sub>ij</sub> = (A<sup>T</sup>)<sub>ji</sub>	 
 *
 * Thus, if A<sub>m�n</sub>, then A'<sub>n.m</sub> =>	 
 *
 *  If A � B = C -> C' = A' � B'			 
 *
 * We may therefore describe a column vector as a row vector
 *
 @f[
 \pmatrix {
    x_1 \cr
    x_2 \cr
    \cdots \cr
    x_{n-1} \cr
    x_n \cr
    } = \pmatrix {x_1 & x_2 & \cdots & x_n \cr }
 @f]
 *
 * @param trn	a matrix where we will store the transpose
 * @param mat	the matrix whose transpose we want to take
 */
public status mat_transpose(matrix trn, matrix mat)
{

    register int i, j;

#   ifdef MAT_PARANOIA
    /* Check dimensions */
    if ((trn->rows != mat->cols) || (trn->cols != mat->rows))
	return MAT_BOUNDSCHECK;
#   endif

#   ifndef MAT_OPTIMIZE
    for (i = 1; i <= mat->rows; i++)
	for (j = 1; j <= mat->cols; j++)
	    trn->values[j][i] = mat->values[i][j];
#   else
    /* we know all elements are stored contiguously */
    {
	int i, size;
	real *trn_elem;
	real *mat_elem;

	size = mat->rows * mat->cols;
	trn_elem = trn->values[1] + size;
	mat_elem = mat->values[1];
	for (i = 1; i <= size; i++) {
	    /* trn_elem[size - i + 1] = mat_elem[size] */
	    *trn_elem = *mat_elem;
	    trn_elem--;
	    mat_elem++;
	}
    }
#   endif
    return SUCCESS;
}


/**
 * @brief MAT B = INV ( C );
 *
 *   Computes the inverted matrix of C, which is a square matrix n x n
 *   and stores the result on B.
 *   The algorithm employed is a modification of Gauss' method. This
 *   is probably not the best (nor is it the nicest implementation),
 *   but for the moment will do.
 *
 *   Given the identity matrix I<sub>nn</sub>, we may think of another
 *   such that
 *
 *  A B = B A = I => B = A<sup>-1</sup>
 *
 *  A A<sup>-1</sup> = I = A<sup>-1</sup> A
 *
 *   Using the cofactors of A<sub>ij</sub> we may compute
 *
 *  A<sup>-1</sup> = 1 / |A| � adj(A)
 *
 *   (calculating through the adjunct matrix). If |A| = 0 then it is not
 *  defined and A is SINGULAR.
 *   
 *   The inverse of a matrix A, if it exists, is unique and may be
 *   found:
 *   We first form for A<sub>nn</sub> the matrix A<sub>n x 2n</sub>
 *   @f[ (A, I) = \pmatrix {
	a_{11} & a_{12} & \cdots & a_{1n} & 1 & 0 & \cdots & 0 \cr
	a_{21} & a_{22} & \cdots & a_{2n} & 0 & 1 & \cdots & 0  \cr
	\vdots & \vdots & \ddots & \vdots & \vdots & \vdots & \ddots & 0  \cr
	a_{n1} & a_{n2} & \cdots & a_{nn} & 0 & 0 & \cdots & 1  \cr }
     @f]
 *   That is, the left half is A and the right half I, the identity
 *   matrix. Using a modified Gauss method we transform the former
 *   matrix into
 *   @f[ (I, B) =
  \pmatrix {
    1	   & 0      & \cdots & 0    &	b_{11} & b_{12} & \cdots & b_{1n}\cr
    0	   & 1      & \cdots & 0    &	b_{21} & b_{22} & \cdots & b_{2n}\cr
    \vdots & \vdots & \ddots & 0    &	\vdots & \vdots & \ddots & \vdots\cr
    0	   & 0      & \cdots & 1    &	b_{n1} & b_{n2} & \cdots & b_{nn}\cr }
     @f]
 *   Now the left half is I and the right half, B is the inverse of A.
 *
 *   To avoid changing matrix C we get first a working copy that
 *   we store in matrix A, whose space is reserved ex-profeso and
 *   is freed upon termination.
 *
 * @callgraph
 */

public status mat_invert(matrix B, matrix C)
				/* matrices inverse and to invert	*/
				/* both are square matrices n x n	*/
{
    int i, j, k, l;		/* indexes		    */
    int n;			/* dimension		    */
    real temp;  		/* auxiliary variable	    */
    matrix A;			/* backup matrix	    */
    boolean singular = FALSE;	/* TRUE if c is singular    */
    real **a, **b;

    n = B->rows;
    mat_alloc(&A, n, n);		/* create A */
    if (A == NULL)		/* (if we can) */
	return MAT_NOMEMORY;
    mat_assign(A, C);		/* copy C into A to avoid changing it */

    /* make B equal the identity matrix */
    mat_identity(B);

    a = A->values;
    b = B->values;
    
    for (j = 1; j <= n; j++) {
	for (i = j; i <= n; i++)
	    if (a[i][j] != 0)		/* is it singular ?	   */
		goto __Compute__;	/* NO			   */
	    /* The matrix is singular, EXIT */
	    mat_free(A);
	    return MAT_SINGULAR;

      __Compute__:
	for (k = 1; k <= n; k++) {
	    /* swap a[j][k] and a[i][k] */
	    temp = a[j][k];
	    a[j][k] = a[i][k];
	    a[i][k] = temp;

	    /* swap inv[j][k] and inv[i][k] */
	    temp = b[j][k];
	    b[j][k] = b[i][k];
	    b[i][k] = temp;
	}			/* K */

	temp = 1 / a[j][j];
	for (k = 1; k <= n; k++) {
	    a[j][k] *= temp;
	    b[j][k] *= temp;
	}			/* K */

	for (l = 1; l <= n; l++) {
	    if (l != j) {
		temp = -a[l][j];
		for (k = 1; k <= n; k++) {
		    a[l][k] += (temp * a[j][k]);

		    b[l][k] += (temp * b[j][k]);
		}		/* K */
	    }			/* IF */
	}			/* L */
    }				/* J */
    /* B = inverse of C 		    */
    /* A is no longer needed, kiss it good bye.  */
    mat_free(A);
    if (singular == TRUE)
	return (MAT_SINGULAR);  /* so sorry	*/
    else
	return (SUCCESS);	/* happy end	*/
}

/** @} */
