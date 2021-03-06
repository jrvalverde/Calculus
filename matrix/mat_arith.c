/**
 *
 *  @file   mat_arith.c
 *
 *  @brief  Basic arithmetic operations with matrices.
 *
 *  An implementation of marix arithmetics. This module includes functions
 *  to perform matrix summation, substraction, multiplication and computing
 *  the scalar product. These operations are base to more advanced matrix
 *  calculus manipulations.
 *
 *  @pre    portable.h
 *
 *  @pre    matrix.h
 *
 *  @see    matrix.h for a general introduction to matrices
 *
 *  @see    mat_init.c to learn how you can assign initial values to
 *	    a matrix or its elements.
 *
 *  @see    mat_ops.c	to learn more about how to perform basic
 *	    matrix operations.
 *
 *  @see    mat_test.c      Is an auxiliary module with
 *	    the test code of these functions used to debug this
 *	    module (last time everything went OK).
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
 * $Id: mat_arith.c,v 1.1 2004/03/05 18:50:14 jr Exp $
 * $Log: mat_arith.c,v $
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
/*			MATRIX ARITHMETIC			  */
/*----------------------------------------------------------------*/

/** \defgroup matrix_arithmetic
 * @{
 */
 

/** 
 * \brief MAT C = A + B;
 *
 *  sum two matrices and store result in another.
 *
 *  function mat_sum adds two matrices and stores the result in a
 *  separate one. Both matrices must have the same number of rows
 *  and columns:
 *
 *  \f[
    \sum A_{mn} + B_{mn} = C_{mn} / c_{ij} = a_{ij} + b_{ij} 
    (i = 1 \ldots m, j = 1 \ldots n)
    \f]
 *
 * \param result    a matrix to store the sum
 * \param mat1      sumand
 * \param mat2      sumand
 *
 *  \return SUCCESS if all went well, an error code otherwise
 *
 */
public status mat_sum(matrix result, matrix mat1, matrix mat2)
{
#   ifdef MAT_PARANOIA
    /* check bounds */
    if ((result->rows != mat1->rows) || (mat1->rows != mat2->rows) ||
	(result->cols != mat1->cols) || (mat1->cols != mat2->cols))
	return MAT_BOUNDSCHECK;
#   endif

#ifndef MAT_OPTIMIZE
    {
	register int i, j;
	real **sum, **val1, **val2;

	sum = result->values;
	val1 = mat1->values;
	val2 = mat2->values;
	for (i = 1; i <= result->rows; i++)
	    for (j = 1; j <= result->cols; j++)
		sum[i][j] = val1[i][j] + val2[i][j];
    }
#else
    {
	/* we know elements are stored contiguously, hence we can
	   treat is as unidimensional as well, saving on indirections */
	register int i, size;
	real *res, *elem1, *elem2;

	size = result->rows * result->cols;
	res = result->values[1];
	elem1 = mat1->values[1];
	elem2 = mat2->values[1];
	    for (i = 1; i <= size; i++)
	    res[i] = elem1[i] + elem2[i];
    }
#endif
    return SUCCESS;
}


/**
 * \brief MAT res = mat1 - mat2
 *
 * Substracts mat2 from mat1 and stores result in res.
 *
 *  Subtract one matrix from another and store the result in another one.
 *  All three matrices must have the same dimensions:
 *
 *  \f[ A_{mn} - B_{mn} = E_{mn} / e_{ij} = a_{ij} - b_{ij}
 *  (i = 1 \ldots m; j = 1 \ldots n) \f]
 *
 * \param result	The resulting substraction matrix
 * \param mat1  Minuend matrix
 * \param mat2  Substraend matrix
 *
 * \return SUCCESS if all went well, an error code otherwise
 *
 */
public status mat_substract(matrix result, matrix mat1, matrix mat2)
{
#   ifdef MAT_PARANOIA
    /* check bounds */
    if ((result->rows != mat1->rows) || (mat1->rows != mat2->rows) ||
	(result->cols != mat1->cols) || (mat1->cols != mat2->cols))
	return MAT_BOUNDSCHECK;
#   endif

#ifndef MAT_OPTIMIZE
    {
	register int i, j;
	real **sub, **val1, **val2;

	sub = result->values;
	val1 = mat1->values;
	val2 = mat2->values;
	for (i = 1; i <= result->rows; i++)
	    for (j = 1; j <= result->cols; j++)
		sub[i][j] = val1[i][j] - val2[i][j];
    }
#else
    {
	/* we know elements are stored contiguously, hence we can
	   treat is as unidimensional as well, saving on indirections */
	register int i, size;
	real *res, *elem1, *elem2;

	size = result->rows * result->cols;
	res = result->values[1];
	elem1 = mat1->values[1];
	elem2 = mat2->values[1];
	    for (i = 1; i <= size; i++)
	    res[i] = elem1[i] - elem2[i];
    }
#endif
    return SUCCESS;
}


/**
 * \brief   MAT prod[m][p] = fact1[m][n] * fact2[n][p];
 *  Compute the product of two matrices
 *
 * Compute the product of two matrices and store result in another one.
 * The number of columns in the first factor matrix must be equal to
 * the number of rows in the second factor matrix. The resulting product 
 * matrix will have same number of rows as first factor and same number
 * of columns as second factor.
 *
 * Let A<sub>mn</sub> and B<sub>np</sub>:
 * \f[
    A � B = C / A_{mn} B_{np} = C_{mp} /
    c_{ij} = \sum_{k=1}^n a_{ik} b_{kj}
    (i = 1 \ldots m; j = 1 \ldots p)
   \f]
 *
 * If p == 1 (i. e. B is a vector with a single column) then
 * A � x for A<sub>mn</sub> and x<sub>n1</sub> is
 * \f[
    Y_{m1} / y_i = \sum_{k=1}^n a_{ik} x_k  (i = 1 \ldots m)
 \f]
 *
 * \param prod  a matrix to store the resulting product
 * \param fact1 first factor to multiply
 * \param fact2 second factor matrix to multiply
 *
 *  \return SUCCESS if all went well, an error code otherwise
 *
 */
public status mat_multiply(matrix prod, matrix fact1, matrix fact2)
{
#   ifdef MAT_PARANOIA
    if ((prod->rows != fact1->rows) ||
	(prod->cols != fact2->cols) || (fact1->cols != fact2->rows))
	return MAT_BOUNDSCHECK;
#   endif
    {
	auto int i, j, k;
	auto int m, n, p;
	auto real temp;
	auto real **p_vals, **f1_vals, **f2_vals;
	m = prod->rows;
	n = fact1->cols;
	p = fact2->cols;
	p_vals = prod->values;
	f1_vals = fact1->values;
	f2_vals = fact2->values;

	for (i = 1; i <= m; i++)
	    for (j = 1; j <= p; j++) {
		temp = 0.0;
		for (k = 1; k <= n; k++)
		    temp += (f1_vals[i][k] * f2_vals[k][j]);
		p_vals[i][j] = temp;
	    }
    }
    return SUCCESS;
}


/**
 *  \brief compute the scalar product 
 *
 *  Compute the product of a matrix by a scalar quantity and store
 * the result in another matrix: \f[
    B = \lambda \cdot A (\lambda = constant)
    \f]\f[
    b_{ij} = \lambda \cdot a_{ij} (i = 1 \ldots m, j = 1 \ldots n)
 \f]
 *
 * \param scprod    the matrix to store the result
 * \param mat	    the matrix to multiply
 * \param lambda    the scalar operand
 *
 *  \return SUCCESS if all went well, an error code otherwise
 *
 */
public status mat_scalar_product(matrix scprod, matrix mat, real lambda)
{
#ifndef MAT_OPTIMIZE
    register int i, j;

    for (i = 1; i <= mat->rows; i++)
	for (j = 1; j <= mat->cols; j++)
	    scprod->values[i][j] = mat->values[i][j] * lambda;
#else
    register int i, size;
    register real *result, *val;

    /* we know all rows were allocated contiguously starting
     * at mat->values[1] hence we may save on pointer indirections */
    size = mat->rows * mat->cols;
    result = scprod->values[1];
    val = mat->values[1];
    for (i = 1; i <= size; i++)
	result[i] = val[i] * lambda;

#endif
}


/** @} */
