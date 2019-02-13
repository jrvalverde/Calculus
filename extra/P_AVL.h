/*
 *  P_AVL.H
 *
 *	File to include to provide access to the routines for
 *  special manipulation of AVL trees in module P_AVL.C. The
 *  remaining functions must be provided by the module P_Stree.C
 *  (binary search tree subroutines).
 *
 *	Designed by:
 *	    J. R. Valverde	8  - Apr - 1990
 */


typedef void *item_t;

struct avl_node {
    item_t	    key;
    struct avl_node *left;
    struct avl_node *right;
    flag	    equi;
};

typedef struct avl_node *avl_tree;

extern void avl_insert(avl_tree *, item_t *, int (*)(), bool *);

extern void avl_delete(avl_tree *, item_t *, int(*)(), bool *);

/*
 *
 */
