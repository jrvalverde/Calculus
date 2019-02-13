/**
 *  @file mat_init.c
 *
 *  @brief 2D-matrix and elements initialization routines.
 *
 *  This module implements functions to initialize matrix element values.
 *  In addition to allowing assignment to any matrix element, functions
 *  are provided to initialize a full matrix to a set of predefined
 *  common values (e.g. the identity matrix) or to set all elements to the
 *  same value (e.g. 0.0 for the null matrix).
 *
 *  There is a second way of initializing/accessing matrix elements, which
 *  consists in querying the matrix for its full value set, or a whole row
 *  and handling it like a normal C array (albeit a one-offset one). This
 *  may be faster under some conditions for performing many asignments
 *  or operations. However, in general, you are advised to use these
 *  functions for initialization.
 *
 *  @note   Should you decide to access the matrix values directly as a
 *	    C array, you must always keep in mind that they are all one-offset
 *	    (i.e. subindexes start at one) instead of zero-offset like in C.
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
 *  @see    mat_ops.c	to learn more about how to perform basic
 *	    matrix operations.
 *
 *  @author José Ramón Valverde Carrillo    (jrvalverde@acm.org)
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
 *	    © YoEgo.	Since I have no cash, I can't
 *	register this (nor do I believe I should). So
 *	this module is left in the PUBLIC DOMAIN.
 *	    It is furthermore forbidden its use for
 *	commercial purposes unless I get a share on
 *	the profits.
 *	    I say.
 *						YoEgo.
 *
 * $Id: mat_init.c,v 1.1 2004/03/05 18:50:14 jr Exp $
 * $Log: mat_init.c,v $
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



/*--------------------------------------------------------------*/
/*			MATRIX INITIALIZATION			*/
/*--------------------------------------------------------------*/

/** @defgroup matrix_initialization
 *  @{
 */

/**
 *  @fn mat_set(matrix mat, int row, int col, real value)
 *
 *  @brief set matrix element mat[row][col] to value
 *
 *  @param mat  the matrix whose value we want to set
 *  @param row  value row offset
 *  @param col  value column offset
 *  @param value    value to assign to mat[row][col]
 *
 *  @return SUCCESS if all went well, an error code otherwise
 *
 *  @note if MAT_PARANOID is defined, extra bound checking is
 *	performed at a high execution penalty.
 */
public status mat_set(matrix mat, int row, int col, real value)
{
#   ifdef MAT_PARANOIA
    /* check bounds */
    if ((row <= 0) || (row > mat->rows) || (col <= 0) || (col > mat->cols))
	return MAT_BOUNDSCHECK;
#   endif

    mat->values[row][col] = value;
    return SUCCESS;
}


/**
 *  @brief Set all matrix elements to the specified value
 *
 *  @param mat  a matrix allocated by mat_alloc() whose values will be all
 *		set to the real value specified
 *  @param value the value to assign to all matrix elements
 *
 *  @return SUCCESS if all went well, an error code otherwise
 *
 */
public status mat_init(matrix mat, real value)
{
#ifndef MAT_OPTIMIZE
    register int i, j, maxrow, maxcol;
    real **val;

    val = mat->values;
    maxrow = mat->rows;
    maxcol = mat->cols;
    for (i = 1; i <= maxrow; i++)
	for (j = 1; j <= maxcol; j++)
	    val[i][j] = value;
#else
    register int i, size;
    register real *values;

    /* we know all rows were allocated contiguously starting
     * at mat->values[1] hence we may save on pointer indirections */
    size = mat->rows * mat->cols;
    values = mat->values[1];
    for (i = 1; i <= size; i++)
	values[i] = value;
#endif
    return SUCCESS;
}


/**
 *  @brief Compute the identity matrix of dimension [n] x [n]:
 *	    Every element ij / i <> j is 0, every ij / i = j is 1.
 *
 *  @param mat  	a matrix allocated by mat_alloc() that will be set to
 *			the identity matrix
 *
 *  @return SUCCESS if all went well, an error code otherwise
 *
 */
public status mat_identity(matrix mat)
{
    register int i, j;
    real **val;

#   ifdef MAT_PARANOIA
    if (mat->rows != mat->cols)
	/* must be a square matrix */
	return MAT_NOTSQUARE;
#   endif

    val = mat->values;
    for (i = 1; i <= mat->rows; i++) {
	for (j = 1; j <= mat->cols; j++)
	    val[i][j] = 0.0;
	val[i][i] = 1.0;
    }
    return SUCCESS;
}

/*  @} */
