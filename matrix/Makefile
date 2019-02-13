# Makefile for matrix calculus module
#
#   Note: should be automatically generated.
#
#   $Id: Makefile,v 1.2 2004/03/05 18:49:57 jr Exp $
#   $Log: Makefile,v $
#   Revision 1.2  2004/03/05 18:49:57  jr
#   Makefile for matrix calculus module [j]
#
#

CC = cc
CFLAGS = -g -I. -I..
CDEFS = -DMAT_OPTIMIZE
DIFF = xxdiff

INC = matrix.h ../portable.h
LIBSRC = mat_housekeep.c mat_init.c mat_access.c mat_ops.c \
      mat_arith.c mat_properties.c
LIBOBJ = mat_housekeep.o mat_init.o mat_access.o mat_ops.o \
      mat_arith.o mat_properties.o

.c.o:
	$(CC) -c $(CFLAGS) $*.c -o $*.o

libmatrix.a: $(LIBOBJ)
	ar r libmatrix.a $?
	ranlib libmatrix.a

debug: mat_test.c libmatrix.a $(INC)
	$(CC) $(CFLAGS) -o mt $(CDEFS) -I. -I.. mat_test.c -L. -lmatrix -lm
	-./mt > mat_test.out
	$(DIFF) mat_test.out mat_test.out.ok

clean: 
	-rm *.o
	-rm core core.*
	-rm mon.out
	-rm *.Addrs
	-rm *.Counts
	-rm *.pixie
