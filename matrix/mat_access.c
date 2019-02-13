/**
 *  @file   mat_access.c
 *
 *  @brief  Routines to access matrix elements and various bits of data
 *
 *  This module hids the matrix implementation by providing various ways
 *  to query and obtain data from a matrix data object. While you may
 *  access some data structure fields directly, it is strongly advised 
 *  that you don't and that all accesses to internal or descriptive data
 *  be done through the interfaces included in this module. This way,
 *  should the implementation change, your program won't break.
 *
 *  For all purposes, a matrix is a bi-dimensional object that is composed
 *  of 'real' numbers arranged in rows and columns. A matrix M<sub>m·n</sub>
 *  has m rows and n columns. You may query de data stored either individually
 *  (indexed by row and column), arranged as rows (a uni-dimensional C vector
 *  or array with offset at one --1), or as a whole bi-dimensional C array
 *  with offsets at one --1)
 *
 *  @note   It is important to notice that offsets in matrices all start at 
 *	    one unlike conventional C arrays. This is intentionally done to
 *	    comply with the standard mathematical convention and ease use
 *	    and translation of mathematical methods.
 *
 *  @note   In previous (v1) implementations, matrix elements were stored
 *	    at zero-offset fr maximal compatibility with C, but this proved
 *	    conceptually more cumbersome to extend and more difficult to 
 *	    implement and debug.
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
 * $Id: mat_access.c,v 1.1 2004/03/05 18:50:14 jr Exp $
 * $Log: mat_access.c,v $
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


/*-------------------------------------------------------------------*/
/*			    MATRIX ACCESS			     */
/*-------------------------------------------------------------------*/

/** @defgroup matrix_access
 *  @{
 */

/** 
 * @brief return number of rows of a matrix
 *
 * @param mat a matrix allocated by mat_alloc() whose row order we want to know
 * @return  the number of rows of the matrix
 */
public int mat_rows(matrix mat)
{
    return mat->rows;
}

/**
 * @brief return number of columns of a matrix
 *
 * @param mat a matrix allocated by mat_alloc() whose column dimension we want to know
 * @return  the number of columns of the matrix
 */
public int mat_cols(matrix mat)
{
    return mat->cols;
}


/**
 *  @brief return a (handy) pointer to the table of values of a matrix
 *
 *  @param mat the matrix whose value table we want to access
 *  @return a 2D table with the values of the matrix arranged by rows
 *	    and columns.
 */
public real **mat_values(matrix mat)
{
    return mat->values;
}

/**
 *  @brief   return a row of a matrix
 *
 *  This function returns a pointer to the values of the specified row
 *  R<sub>i</sub>(mat)
 *
 *  @param mat  	a matrix allocated by mat_alloc
 *  @param row_number	the index of the row to return
 *  @return		a pointer to an array [1..dim] containing the row values
 */
public real *mat_row(matrix mat, int row_number)
{
#   ifdef MAT_PARANOIA
    if (row_number > mat->rows)

#   endif
    return mat->values[row_number];
}


/**
 *  @fn public real mat_element(matrix mat, int row, int col)
 *
 *  @brief return value of matrix element mat[row][col]
 *
 *  @param mat  the matrix whose value we want to set
 *  @param row  value row offset
 *  @param col  value column offset
 *
 *  @return the value of the element.
 *
 *  @note if MAT_PARANOIA is defined, extra bound checking is
 *	performed and IEEE NaN returned on error.
 */
public real mat_element(matrix mat, int row, int col)
{
#   ifdef MAT_PARANOIA
    /* check bounds */
    if ((row <= 0) || (row > mat->rows) || (col <= 0) || (col > mat->cols))
	/* IEEE Not A Number */
	return NAN;
#   endif

    return mat->values[row][col];
}

/** @} */
