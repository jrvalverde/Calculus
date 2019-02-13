/*
 *  P_DLIST.H
 *
 *	Include file containing the definitions for the module
 *  P_DList.C
 *
 *	Requires PORTABLE.H and P_TYPES.H
 *
 *	Designed by:
 *	    J. R. Valverde	8  - apr - 1990
 */

extern dlist_t new_dlist();

extern dl_node_t *new_dnode(info_t *, dl_node_t *, dl_node_t *);

extern boolean dlist_put(info_t *, dlist_t );

extern boolean dlist_append(info_t *, dlist_t );

extern info_t *dlist_del_first(dlist_t);

extern bool dlist_insert(info_t *, dl_node_t *, dlist_t);

extern info_t *dlist_del_last(dlist_t);

extern info_t *dlist_del_next(info_t *, dlist_t);

extern info_t *dlist_del_node(dl_node_t *, dlist_t);

extern int dlist_length( dlist_t );

extern info_t *dlist_first( dlist_t );

extern info_t *dlist_last( dlist_t );

extern info_t *dlist_next( dlist_t );

extern info_t *dlist_prev( dlist_t );

extern info_t *dlist_begset( dlist_t );

extern info_t *dlist_endset( dlist_t );

extern bool dlist_ffind(dlist_t, info_t *, int (*)());

extern bool dlist_bkfind(dlist_t, info_t *, int (*)());

extern void dlist_ftraverse(dlist_t, void (*)());

extern void dlist_bktraverse(dlist_t, void (*)());

extern info_t *dlist_del_first(dlist_t);

extern info_t *dlist_del_last(dlist_t);

extern info_t *dlist_fdelete(dlist_t, info_t *, int (*)());

extern info_t *dlist_bkdelete(dlist_t, info_t *, int (*)());

/*
 *  Nothing compares 2 u!
 */

