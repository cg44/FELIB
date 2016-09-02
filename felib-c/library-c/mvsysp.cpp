/* mvsysp.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"

/* *********************************************************************** */
/* Subroutine */ int mvsys_(a, ia, row, irow, diag, idiag, v, iv, w, iw, n, 
	itest)
doublereal *a;
integer *ia, *row, *irow, *diag, *idiag;
doublereal *v;
integer *iv;
doublereal *w;
integer *iw, *n, *itest;
{

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer jrow, i, j, id, in;

    /* Parameter adjustments */
    --w;
    --v;
    --diag;
    --row;
    --a;

    /* Function Body */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	w[i] = (float)0.;
/* L1: */
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	id = diag[i];
	in = diag[i + 1] - 1;
	i__2 = in;
	for (j = id; j <= i__2; ++j) {
	    jrow = row[j];
	    w[jrow] += a[j] * v[i];
/* L50: */
	}
	++id;
	if (in < id) {
	    goto L100;
	}
	i__2 = in;
	for (j = id; j <= i__2; ++j) {
	    jrow = row[j];
	    w[i] += a[j] * v[jrow];
/* L60: */
	}
L100:
	;
    }
    return 0;
} /* mvsys_ */

