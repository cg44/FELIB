/* csysub.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int csysub_(kb, ikb, jkb, loads, iloads, n, hband, itest)
doublereal *kb;
integer *ikb, *jkb;
doublereal *loads;
integer *iloads, *n, *hband, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "CSYSUB";

    /* System generated locals */
    integer kb_dim2, kb_offset, i__1, i__2;

    /* Local variables */
    static integer i, j, k, l, m, w;
    static doublereal x, y, ai;
    static integer ij;
    static doublereal ar, xi, xr;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      CSYSUB performs forward and backward substitution on a matrix */
/*      reduced by CSYRDN */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 3.0  29 Oct 1984 (CRIE) */
/*      Commented    10 Oct 1985 (CG) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      KB      array of dimension (IKB, JKB).  contains the */
/*              elements of the lower half of the complex symmetric */
/*              band matrix of order N and semi-bandwidth HBAND. KB */
/*             should previously have been reduced using CSYSOL or CSYRDN.
*/
/*      IKB     first dimension of KB (.GE. N) */
/*      JKB     second dimension of KB (.GE. HBAND) */
/*      LOADS   on entry, contains the vector of rhs's */
/*      ILOADS  dimension of LOADS (.GE. N) */
/*      N       order of matrix KB */
/*      HBAND   semi-bandwidth of KB */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      LOADS   on exit, contains the solution vector */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE CSYSUB(KB,IKB,JKB,LOADS,ILOADS,N,HBAND,ITEST) */
/* ********************************************************************** 
*/

    /* Parameter adjustments */
    loads -= 3;
    kb_dim2 = *ikb;
    kb_offset = (kb_dim2 + 1 << 1) + 1;
    kb -= kb_offset;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*iloads < *n) {
	    ierror = 3;
	}
	if (*ikb < *n || *jkb < *hband) {
	    ierror = 2;
	}
	if (*n <= 0 || *hband <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    w = *hband - 1;
    x = loads[3];
    y = loads[4];
    ar = kb[((w + 1) * kb_dim2 + 1 << 1) + 1];
    ai = kb[((w + 1) * kb_dim2 + 1 << 1) + 2];
    if (abs(ar) > abs(ai)) {
	loads[3] = (x + ai / ar * y) / (ai / ar * ai + ar);
	loads[4] = (y - ai / ar * x) / (ai / ar * ai + ar);
    } else {
	loads[3] = (ar / ai * x + y) / (ar / ai * ar + ai);
	loads[4] = (ar / ai * y - x) / (ar / ai * ar + ai);
    }
    i__1 = *n;
    for (i = 2; i <= i__1; ++i) {
	xr = 0.;
	xi = 0.;
	k = 1;
	if (i <= w + 1) {
	    k = w - i + 2;
	}
	i__2 = w;
	for (j = k; j <= i__2; ++j) {
	    ij = i + j - w - 1;
	    x = kb[(i + j * kb_dim2 << 1) + 1];
	    y = kb[(i + j * kb_dim2 << 1) + 2];
	    ar = loads[(ij << 1) + 1];
	    ai = loads[(ij << 1) + 2];
	    xr = xr + x * ar - y * ai;
	    xi = xi + x * ai + y * ar;
/* L1000: */
	}
	x = loads[(i << 1) + 1] - xr;
	y = loads[(i << 1) + 2] - xi;
	ar = kb[(i + (w + 1) * kb_dim2 << 1) + 1];
	ai = kb[(i + (w + 1) * kb_dim2 << 1) + 2];
	if (abs(ar) > abs(ai)) {
	    loads[(i << 1) + 1] = (x + ai / ar * y) / (ai / ar * ai + ar);
	    loads[(i << 1) + 2] = (y - ai / ar * x) / (ai / ar * ai + ar);
	} else {
	    loads[(i << 1) + 1] = (ar / ai * x + y) / (ar / ai * ar + ai);
	    loads[(i << 1) + 2] = (ar / ai * y - x) / (ar / ai * ar + ai);
	}
/* L1010: */
    }
    x = loads[(*n << 1) + 1];
    y = loads[(*n << 1) + 2];
    ar = kb[(*n + (w + 1) * kb_dim2 << 1) + 1];
    ai = kb[(*n + (w + 1) * kb_dim2 << 1) + 2];
    if (abs(ar) > abs(ai)) {
	loads[(*n << 1) + 1] = (x + ai / ar * y) / (ai / ar * ai + ar);
	loads[(*n << 1) + 2] = (y - ai / ar * x) / (ai / ar * ai + ar);
    } else {
	loads[(*n << 1) + 1] = (ar / ai * x + y) / (ar / ai * ar + ai);
	loads[(*n << 1) + 2] = (ar / ai * y - x) / (ar / ai * ar + ai);
    }
    i = *n - 1;
L1020:
    xr = 0.;
    xi = 0.;
    l = i + w;
    if (i > *n - w) {
	l = *n;
    }
    m = i + 1;
    i__1 = l;
    for (j = m; j <= i__1; ++j) {
	ij = w + i - j + 1;
	x = kb[(j + ij * kb_dim2 << 1) + 1];
	y = kb[(j + ij * kb_dim2 << 1) + 2];
	ar = loads[(j << 1) + 1];
	ai = loads[(j << 1) + 2];
	xr = xr + x * ar - y * ai;
	xi = xi + x * ai + y * ar;
/* L1030: */
    }
    x = loads[(i << 1) + 1] - xr;
    y = loads[(i << 1) + 2] - xi;
    ar = kb[(i + (w + 1) * kb_dim2 << 1) + 1];
    ai = kb[(i + (w + 1) * kb_dim2 << 1) + 2];
    if (abs(ar) > abs(ai)) {
	loads[(i << 1) + 1] = (x + ai / ar * y) / (ai / ar * ai + ar);
	loads[(i << 1) + 2] = (y - ai / ar * x) / (ai / ar * ai + ar);
    } else {
	loads[(i << 1) + 1] = (ar / ai * x + y) / (ar / ai * ar + ai);
	loads[(i << 1) + 2] = (ar / ai * y - x) / (ar / ai * ar + ai);
    }
    --i;
    if (i != 0) {
	goto L1020;
    }

} /* csysub_ */

