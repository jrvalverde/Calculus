/*
 *
 *  P_sort.c
 *
 *    Library of generic array sorting procedures.
 *
 *  CONTENTS:
 *    DirSort(arr, nelm, comp)
 *    BinSort(arr, nelm, comp)
 *    SelSort(arr, nelm, comp)
 *    BubSort(arr, nelm, comp)
 *    ShkSort(arr, nelm, comp)
 *    ShellSort(arr, nelms, comp)
 *    HeapSort(arr, nelms, comp)
 *    RQkSort(arr, beg, end, comp)
 *    RQuickSort(arr, nelms, comp)
 *    IQuickSort(arr, nelms, comp)
 *
 *  SEE ALSO:
 *    Algorithms + data structures = programs
 *    Niklaus Wirth.
 *
 *  DESIGNED BY:
 *    Jos� Ram�n Valverde Carrillo.
 *
 *  LAST MODIFICATION:
 *     3 - Dic - 1988 J. R. Valverde - Mac version
 *    10 - Dic - 1988 J. R. Valverde - Port to IBM-PC
 *    12 - Jul - 1989 J. R. Valverde - Port to uVAX-VMS
 *     1 - Aug - 1991 J. R. Valverde - Port to Ultrix
 *    18 - May - 1997 J. R. Valverde - Port to Linux
 *    24 - May - 1997 J. R. Valverde - Port to FreeBSD
 *
 * $Log: P_sort.c,v $
 * Revision 1.2  1997/05/24 17:13:16  jr
 * Ported to FreeBSD
 *
 * Revision 1.1  1997/05/18 15:12:36  jr
 * Initial revision
 *
 */

#define MOD_DEBUG   0

#define MINOR	-1
#define EQUAL	 0
#define GREATER	 1

#define TRUE 	1
#define FALSE	0

typedef char flag;

typedef int bool;

typedef void *pitem;

typedef int index_t;

/*********************************************************************
*
*   Direct methods.
*
*    2
*   n   comparisons. Fasters for n small.
*   
*   C = N� of key comparisons
*   M = N� of item movements
*
*   J. R. VALVERDE
*********************************************************************/

/*--------------------------------------------------------------------
    DIRSORT
	Direct insertion without sentinel.

    J. R. VALVERDE
--------------------------------------------------------------------*/

void DirSort(arr, nelm, comp)
pitem arr[];	    /* array of pointers to items */
index_t nelm;
int (*comp)();
  {
    pitem x;
    index_t i, j;
    
    for (i = 0; i < nelm; i++)
      {
        x = arr[i];
        j = i;
	/* while there is still an element below and x's value
	is less than it's */
        while ((j) && ((*comp)(x, arr[j - 1]) == MINOR))
          {
            /* move element below one place upwards */
            arr[j] = arr[j - 1];
            /* try next */
            j--;
          }
        /* once found its place, set it */
        arr[j] = x;
      }
  }
  
/*--------------------------------------------------------------------
    BINSORT
	Binary insertion.
    Preferable if destination sequence is already sorted.
	
    C = n (log n - log e +- 0.5)
         2
    M = n

    J. R. VALVERDE
--------------------------------------------------------------------*/

void BinSort(arr, nelm, comp)
pitem arr[];
index_t nelm;
int (*comp)();
  {
    index_t i, j, left, right, m;
    pitem x;
    
    for (i = 1; i < nelm; i++) {
        x = arr[i];
        left = 0;
        right = i - 1;
        while (left <= right) {
            m = (left + right) >> 1;
            if ((*comp)(x, arr[m]) == MINOR)
                right = m - 1;
            else
              left = m + 1;
        }
        for (j = i - 1; j >= left; j--)
          arr[j + 1] = arr[j];
        arr[left] = x;
    }
}

/*--------------------------------------------------------------------
    SELSORT
	Direct selection sort.
		Selects minimum key element
		Exchanges it with the first
		Follows for n = 1 to n
                2
    C = 1 / 2 (n  - n)
    
    Mmin = 3 (n - 1)
                  2
    Mmax = trunc(n  / 4) + 3 (n - 1)
    
    Mmed = n (ln n + E)

(E = Euler's constant = 0.577216...)

    J. R. VALVERDE
--------------------------------------------------------------------*/

void SelSort(arr, nelm, comp)
pitem arr[];
index_t nelm;
int (*comp)();
  {
    index_t i, j, k;
    pitem x;
    
    for (i = 0; i < (nelm - 1); i++)
      {
        k = i;
        x = arr[i];
        for (j = i + 1; j < nelm; j++)
          if ((*comp)(arr[j], x) == MINOR)
            {
              k = j;
              x = arr[j];
            }
        arr[k] = arr[i];
        arr[i] = x;
      }
  }

/*-------------------------------------------------------------------
    BUBSORT
	Direct exchange sort.
	Bubble method.
                  2
      C = 1 / 2 (n  - n)
      
      Mmin = 0;
                     2
      Mmax = 3 / 2 (n  - n)
                     2
      Mmed = 3 / 4 (n  - n)

    J. R. VALVERDE
--------------------------------------------------------------------*/

void BubSort(arr, nelm, comp)
pitem arr[];
index_t nelm;
int (*comp)();
  {
    index_t i, j, k;
    pitem x;
    
    for (i = 0; i < nelm; i++)
    {
        k = i;
        for (j = nelm - 1; j > i; j--)
            if ((*comp)(arr[j - 1], arr[j]) == GREATER) {
            /* if (arr[j - 1] > arr[j]) then swap */
            	x = arr[j - 1];
              	arr[j - 1] = arr[j];
              	arr[j] = x;
              	/* store last change position */
              	k = j;
            }
	/* if in this pass we have not made any change
	then it is already sorted. */
        if (k == i) break;
        else
      	/* below k it is all already sorted. Skip it */
            i = k - 1;
        
      }
  }

/*--------------------------------------------------------------------
    SHKSORT
	Shake method sort.
    
    Does not decrease M with respect to BubSort, but it does with C.
    Advantageous when items are almost sorted.

    J. R. VALVERDE
--------------------------------------------------------------------*/

void ShkSort(arr, nelm, comp)
pitem arr[];
index_t nelm;
int (*comp)();
  {
    index_t j, k, iz, de;
    pitem x;
    
    iz = 1;
    de = k = nelm - 1;
    do
      {
        for (j = de; j >= iz; j--)
          if ((*comp)(arr[j - 1], arr[j]) == GREATER)
            {
              x = arr[j - 1];
              arr[j - 1] = arr[j];
              arr[j] = x;
              k = j;
            }
        iz = k + 1;
        for (j = iz; j <= de; j++)
          if ((*comp)(arr[j - 1], arr[j]) == GREATER)
            {
              x = arr[j - 1];
              arr[j - 1] = arr[j];
              arr[j] = x;
              k = j;
            }
        de = k - 1;
      }
    while (iz <= de);
  }

/*--------------------------------------------------------------------

    SHELLSORT
	D. L. Shell's method sorting.

    Mathematical effort required is proportional to
         1.2
	n
    It sorts by jumps, from h to h. Some common increment sequences are
    
      h_i = 3 * h_i-1 + 1
      h_t = 1
      t = [log3 n] - 1
      	    ( h = 1, 4, 13, 40, 121, ...)
      
    or as well
      h_i = 2 * h_i-1 + 1
      h_t = 1
      t = [log2 n] - 1
      	    ( h = 1, 3, 7, 15, 31, ...)

    It is not clear which increment sequence gives best results, but
    it is desirable not to be ones multiples of the others. If a sequence
    is sorted from k to k and it is sorted from i to i, it still remains
    sorted from k to k.
    The algorithm runs acceptably well for moderately sized arrays (say
    N < 5000).

    J. R. VALVERDE
--------------------------------------------------------------------*/

void ShellSort(arr, nelms, comp)
pitem *arr;
index_t nelms;
int (*comp)();
{
    index_t i, j, h;
    pitem aux;

    /* *** jr *** HACK!!!
     * I took this from a previous Pascal implementation 
     * using 1-offset instead of 0-offset arrays.Next is a
     * trick to change the array offest.
     */
    arr--;
    /* first find biggest h value */
    for (h = 1; h <= nelms / 9; h = 3 * h + 1) ;
    /* and now go down from bigger to smaller step sizes */
    for ( ; h > 0; h /= 3)
    	for ( i = h + 1; i <= nelms; i +=h) {
    	    aux = arr[i];
    	    j = i;
    	    while ((j > h) && ((*comp)(arr[j - h], aux) == GREATER)) {
    	    	arr[j] = arr[j - h];
    	    	j -= h;
    	    }
    	arr[j] = aux;
    }
}


/*********************************************************************
*
*   Tree sorting methods.
*   
*   C = n * log2 n
*   
*********************************************************************/

/*--------------------------------------------------------------------
    HEAPSORT
	Heap method.
    Not advisable for small numbers of items, but it is good for
    big N. Its efficiency is greater as N grows.
    It usually works better with original sequences in which items
    are more or less in inverse order.

    	A heap is a sequence of keys {h_l2, h_l2+1, ... hr} /
    h_i <= h_2i and h_i <= h_2i+1 for every i = left .. right/2
    Element h_1 is the minimum.

    J. R. VALVERDE
--------------------------------------------------------------------*/

static void screen(x, a, iz, de, comp)
pitem x, a[];
index_t iz, de;
int (*comp)();
  {
    index_t i, j;
    
    i = iz;
    j = i << 1;
    x = a[i];
    while (j <= de)
      {
        if (j < de)
            if ((*comp)(a[j], a[j + 1]) == MINOR)
              j++;
        if ((*comp)(x, a[j]) == MINOR) {
            a[i] = a[j];
            i = j;
            j = i << 1;
         }
         else
           break;
       }
     a[i] = x;
  }

void HeapSort(arr, nelms, comp)
pitem *arr;
index_t nelms;
int (*comp)();
{
    void criba();
    index_t iz, de;
    pitem x;
    
    /* translated from pascal: we need a 1-offset array */
    arr--;
    iz = (nelms >> 1) + 1;
    de = nelms;
    while (iz > 1) {
        iz--;
        screen(&x, arr, iz, de, comp);
    }
    while (de > 1) {
        x = arr[1];
        arr[1] = arr[de];
        arr[de] = x;
        de--;
        screen(&x, arr, iz, de, comp);
    }
}

/********************************************************************
*
*   Partition sorting: Quick Sort.
*
********************************************************************/

/*-------------------------------------------------------------------
    RQUICKSORT
	Quick method. Recursive version.

    J. R. VALVERDE
--------------------------------------------------------------------*/

void RQkSort(ar, iz, de, comp)
pitem ar[];
index_t iz, de;
int (*comp)();
  {
    index_t i, j;
    pitem x;
    
    i = iz;
    j = de;
    x = ar[(iz + de) >> 1];
    do
      {
        while ( (*comp)(ar[i], x) == MINOR ) i++;
        while ( (*comp)(x, ar[j]) == MINOR ) j--;
        if (i <= j)
          {
            pitem w;
            
            w = ar[i];
            ar[i] = ar[j];
            ar[j] = w;
            i++;
            j--;
          }
      }
    while (i <= j);
    if (iz < j)
      RQkSort(ar, iz, j, comp);
    if (i < de)
      RQkSort(ar, i, de, comp);
  }

void RQuickSort(arr, nelms, comp)
pitem arr[];
index_t nelms;
int (*comp)();
  {
    void RQkSort();

    RQkSort(arr, 0, nelms - 1, comp);
  }

/*--------------------------------------------------------------------
    IQUICKSORT
	Quick method. Iterative version.
    
    C = n;
    
    M = (n / 6) - (1 /(6 * n)) aprox= n / 6
    
    Caverage = n log n
    Mmin = log2 n

    J. R. VALVERDE
--------------------------------------------------------------------*/

#include <stdio.h>
#include <signal.h>

void IQuickSort(arr, nelms, comp)
pitem arr[];
index_t nelms;
int (*comp)();
  {

/* If I'm not mistaken, this should work for up to 2**256 items
 * Maybe it should rather be (sizeof(index_t) << 3) the number of
 * bits used for indexing the array.
 */
#define MAX_STACK_DEPTH	    256
    
    index_t i, j, iz, de;
    pitem x, w;
    short p;	    /* index in the stack */
    struct
      {
        index_t iz;
        index_t de;
      }
    stack[MAX_STACK_DEPTH];
    
    p = 0;
    stack[0].iz = 0;
    stack[0].de = nelms - 1;
    do
      {
	/* take upper request in the stack */
        iz = stack[p].iz;
        de = stack[p].de;
        p--;
        do
          {
            /* subdivision arr[iz] ... arr[de] */
            i = iz;
            j = de;
            x = arr[(iz + de) >> 1];
            do
              {
                while ((*comp)(arr[i], x) == MINOR) i++;
                while ((*comp)(x, arr[j]) == MINOR) j--;
                if (i <= j)
                  {
                    w = arr[i];
                    arr[i] = arr[j];
                    arr[j] = w;
                    i++;
                    j--;
                  }
              }
            while (i <= j);
            if (i < de)
              {
                /* store request in the stack */
                p++;
                if (p >= MAX_STACK_DEPTH) {
                    /* do something useful */
                    fprintf(stderr,"\nERROR: stack overrun on IQuickSort\n");
                    raise(SIGSEGV);
                }
                stack[p].iz = i;
                stack[p].de = de;
              }
            de = j;
          }
        while (iz < de);
      }
    while (p != -1);
  }


/*---------------------------------------------------------------------------*/


#if MOD_DEBUG != 0

#include <stdio.h>
#include <stdlib.h>

#define ARRAYSIZE   100000

int icmp(pi, pj)
int *pi, *pj;
{
    return ( (*pi > *pj)? GREATER : ((*pi == *pj)? EQUAL : MINOR) );
}

long int iary[ARRAYSIZE];
long int *pary[ARRAYSIZE];

main()
{
    long int i, j, k;

    /* initialize array */
    for (i = 0; i < ARRAYSIZE; i++) {
    	iary[i] = random();
    	pary[i] = &iary[i];
    }

    /* sort it */
#if MOD_DEBUG == 1
    DirSort(pary, ARRAYSIZE, icmp);
#endif
#if MOD_DEBUG == 2
    BinSort(pary, ARRAYSIZE, icmp);
#endif
#if MOD_DEBUG == 3
    SelSort(pary, ARRAYSIZE, icmp);
#endif
#if MOD_DEBUG == 4
    BubSort(pary, ARRAYSIZE, icmp);
#endif
#if MOD_DEBUG == 5
    ShkSort(pary, ARRAYSIZE, icmp);
#endif
#if MOD_DEBUG == 6
    ShellSort(pary, ARRAYSIZE, icmp);
#endif
#if MOD_DEBUG == 7
    HeapSort(pary, ARRAYSIZE, icmp);
#endif
#if MOD_DEBUG == 8
    RQuickSort(pary, ARRAYSIZE, icmp);
#endif
#if MOD_DEBUG == 9
    IQuickSort(pary, ARRAYSIZE, icmp);
#endif
#if MOD_DEBUG == 13
    qsort(iary, ARRAYSIZE, sizeof (int *), icmp);
#endif

    /* verify sort */
    for (j = 0, i = 1; i < ARRAYSIZE; i++, j++)
    	if (*pary[j] > *pary[i])
    	    fprintf(stderr, 
    	    	"ERROR: sorted array out of order [%d]=%d > [%d]=%d\n",
    	    	j, *pary[j], i, *pary[i]);
/*
    for (i = 0; i < ARRAYSIZE; i++)
    	fprintf(stdout, "%d ", iary[i]);
    fprintf(stdout, "\n\n");
    for (i = 0; i < ARRAYSIZE; i++)
    	fprintf(stdout, "%d ", *pary[i]);
    fprintf(stdout, "\n\n");
*/
    fprintf(stderr, " \nOK\n");
}
#endif  
