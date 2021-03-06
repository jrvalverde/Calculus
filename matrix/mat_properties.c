/**
 *  @file mat_properties.c
 *
 *  @brief  functions to query 2D matrix properties
 *
 *  When using 2D matrices it is often useful to find about the properties
 *  of a matrix before performing some operations on it, either to be able
 *  to select a better algorithm, to avoid problems or for any other reason.
 *
 *  The functions in this module allow you to find out if a given matrix
 *  has the queried characteristic.
 *
 *  @note   The results found might be cached in specific fields of the
 *	    matrix structure, with a .dirty field to signal when a given
 *	    operation may have changed the properties of a matrix. The
 *	    general usage of matrices should be investigated to decide if
 *	    this approach is worth the trouble.
 *
 *  @pre    portable.h
 *
 *  @pre    matrix.h
 *
 *  @see    matrix.h for a general introduction to matrices
 *
 *  @see    mat_test.c      Is an auxiliary module with
 *	the test code of these functions used to debug this
 *	module (last time everything went OK).
 *
 *  @author    Jos� Ram�n Valverde Carrillo.
 *
 *  @version	2.0
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
 * $Id: mat_properties.c,v 1.1 2004/03/05 18:50:14 jr Exp $
 * $Log: mat_properties.c,v $
 * Revision 1.1  2004/03/05 18:50:14  jr
 * Initial revision
 *
 * Revision 1.1  2004/02/22 21:13:12  jr
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
/*			MATRIX PROPERTIES			  */
/*----------------------------------------------------------------*/

/** \defgroup matrix_properties
 * @{
 */

/**
 *  \brief query if a matrix is square
 *
 *  Verify if both dimensions of a matrix are the same (i.e. it is square)
 *
 *  \param mat the matrix to verify
 *
 *  \return TRUE if the matrix is square, false otherwise
 */
public boolean mat_is_square(matrix mat)
{
    return (mat->rows == mat->cols);
}

/**
 * \brief check if a matrix is symmetric 
 *
 * Check if a matrix is symmetric. A symmetric matric has all its
 * elements symmetric with respect to the main diagonal (a<sub>ij</sub>
 * = a<sub>ji</sub>): A=A'
 *
 * This is useful for handling quadratic forms: Given A<sub>mn</sub>
 * and X<sub>n1</sub> then \f[
     X^T A X = \sum_{i=1}^n \sum{j=1}^n a_{ij} x{i} x{j}
 \f]
 * and conversely, given a quadratic form, it may be converted to matrix
 * form, but there are many ways to choose a matrix A. One of them is to
 * select A a symmetric matrix: it is then easier to select its values.
 *
 * \param mat the matrix to verify
 *
 * \return TRUE if the matrix is symmetric, FALSE otherwise.
 */
public boolean mat_is_symmetric(matrix mat)
{
#ifndef MAT_OPTIMIZE
    boolean symmetric = TRUE;
    int i, j;
    
    if (mat->rows != mat->cols)
	/* the matrix must be square to start with */
	return FALSE;

    for (i = 1; i <= mat->rows - 1; i++)
	for (j = i + 1; j <= mat->cols; j++)
	    if (mat->values[i][j] != mat->values[j][i])
		return FALSE;
#else
    /* trivial optimizations */
    boolean symmetric = TRUE;
    int i, j, dim;
    real **element;
    
    if ((dim = mat->rows) != mat->cols)
	return FALSE;

    element = mat->values;
    for (i = 1; i <= dim; i++)
	for (j = 1; j <= dim; j++)
	    if (element[i][j] != element[j][i])
		return FALSE;
#endif
}

/**
 * \brief query if this is an identity matrix
 *
 *  Return if the matrix is an identity matrix. The matrix must be a
 * square matrix and have all elements in the main diagonal set to 
 * one, and all off-diagonal elements set to zero.
 *
 * An identity matrix of order n, I<sub>nn</sub> verifies
 * 
 * A<sub>mn</sub> I<sub>nn</sub> = A<sub>mn</sub>
 *
 * I<sub>nn</sub> B<sub>np</sub> = B<sub>np</sub>
 *
 * \param mat the matrix to verify
 * 
 * \return TRUE if the matrix is an identity matrix.
 */
public boolean mat_is_identity(matrix mat)
{
    int i, j, dim;
    
    if ((dim = mat->rows) != mat->cols)
	return FALSE;
    
    for (i = 1; i <= dim; i++)
	for (j = 1; j <= dim; j++)
	    if (i == j) {
		if (mat->values[i][i] != 1.0)
		    return FALSE;
	    }
	    else
		if (mat->values[i][j] != 0.0)
		    return FALSE;
    return TRUE;
}

/**
 * \brief query if a matrix is null
 *
 * Find out if we are dealing with the null matrix. A null matrix has
 * all of its elements set to zero.
 *
 * \param mat the matrix to verify
 * 
 * \return TRUE if the matrix is a null matrix.
 */
public boolean mat_is_null(matrix mat) {
#   ifndef MAT_OPTIMIZE
    int i, j;
    
    for (i = 1; i <= mat->rows; i++)
	for (j = 1; j <= mat->cols; j++)
	    if (mat->values[i][j] != 0.0)
		return FALSE;
#   else
    int i, size;
    real *elem;
    
    size = mat->rows * mat->cols;
    elem = mat->values[i];
    for (i = 1; i <= size; i++)
	if (elem[i] != 0.0)
	    return FALSE;
#   endif
    return TRUE;
}

/** @} */
