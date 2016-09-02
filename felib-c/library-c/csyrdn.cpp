/* csyrdn.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b4 = 1.;


/* Subroutine */ int csyrdn_(kb, ikb, jkb, n, hband, itest)
doublereal *kb;
integer *ikb, *jkb, *n, *hband, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "CSYRDN";

    /* System generated locals */
    integer kb_dim2, kb_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(), d_sign();

    /* Local variables */
    static integer a, b, i, j, k, l, w;
    static doublereal x, y, ai, ar;
    static integer ik, lk;
    static doublereal xi, xr;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      CSYRDN performs symmetric reduction on a complex symmetric banded 
*/
/*      banded matrix.  Only the lower band and diagonal are stored */
/*      in a rectangular array KB. */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 3.0  29 Oct 1984 (CRIE) */
/*      Commented     1 Nov 1985 (CG) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      KB      on entry, contains lower band and diagonal of */
/*              complex symmetric matrix */
/*      IKB     first dimension of KB (.GE. N) */
/*      JKB     second dimension of KB (.GE. HBAND) */
/*      N       order of matrix KB */
/*      HBAND   semi-bandwidth of KB  (including diagonal) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      KB      reduced matrix L, where the input matrix KB  has */
/*              been reduced to triangular matrices L and LT */
/*              where KB =L*LT */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE CSYRDN(KB,IKB,JKB,N,HBAND,ITEST) */
/* ********************************************************************** 
*/

    /* Parameter adjustments */
    kb_dim2 = *ikb;
    kb_offset = (kb_dim2 + 1 << 1) + 1;
    kb -= kb_offset;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
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

/*     Main loops */

    w = *hband - 1;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	xr = 0.;
	xi = 0.;
	i__2 = w;
	for (j = 1; j <= i__2; ++j) {
	    x = kb[(i + j * kb_dim2 << 1) + 1];
	    y = kb[(i + j * kb_dim2 << 1) + 2];
	    xr = xr + x * x - y * y;
	    xi += x * y * 2.;
/* L1000: */
	}

	if (*itest != -1) {
	    ierror = 0;
	    x = kb[(i + (w + 1) * kb_dim2 << 1) + 1] - xr;
	    y = kb[(i + (w + 1) * kb_dim2 << 1) + 2] - xi;
	    if (sqrt(x * x + y * y) <= 0.) {
		ierror = 3;
	    }
	    *itest = errmes_(itest, &ierror, srname, 6L);
	    if (*itest != 0) {
		return 0;
	    }
	}

	if (kb[(i + (w + 1) * kb_dim2 << 1) + 1] - xr >= 0.) {
	    ar = kb[(i + (w + 1) * kb_dim2 << 1) + 1] - xr;
	    ai = kb[(i + (w + 1) * kb_dim2 << 1) + 2] - xi;
	    x = sqrt((ar + sqrt(ar * ar + ai * ai)) / 2.);
	    y = 0.;
	    if (x > 0.) {
		y = ai / (x * 2.);
	    }
	    kb[(i + (w + 1) * kb_dim2 << 1) + 1] = x;
	    kb[(i + (w + 1) * kb_dim2 << 1) + 2] = y;
	} else {
	    ar = kb[(i + (w + 1) * kb_dim2 << 1) + 1] - xr;
	    ai = kb[(i + (w + 1) * kb_dim2 << 1) + 2] - xi;
	    y = d_sign(&c_b4, &ai) * sqrt((abs(ar) + sqrt(ar * ar + ai * ai)) 
		    / 2.);
	    x = 0.;
	    if (y > 0.) {
		x = ai / (y * 2.);
	    }
	    kb[(i + (w + 1) * kb_dim2 << 1) + 1] = x;
	    kb[(i + (w + 1) * kb_dim2 << 1) + 2] = y;
	}
	i__2 = w;
	for (k = 1; k <= i__2; ++k) {
	    xr = 0.;
	    xi = 0.;
	    if (i + k <= *n) {
		if (k != w) {
		    l = w - k;
L1010:
		    ik = i + k;
		    lk = l + k;
		    x = kb[(ik + l * kb_dim2 << 1) + 1];
		    y = kb[(ik + l * kb_dim2 << 1) + 2];
		    ar = kb[(i + lk * kb_dim2 << 1) + 1];
		    ai = kb[(i + lk * kb_dim2 << 1) + 2];
		    xr = xr + x * ar - y * ai;
		    xi = xi + x * ai + y * ar;
		    --l;
		    if (l != 0) {
			goto L1010;
		    }
		}
		a = i + k;
		b = w - k + 1;
		x = kb[(a + b * kb_dim2 << 1) + 1] - xr;
		y = kb[(a + b * kb_dim2 << 1) + 2] - xi;
		ar = kb[(i + (w + 1) * kb_dim2 << 1) + 1];
		ai = kb[(i + (w + 1) * kb_dim2 << 1) + 2];
		if (abs(ar) > abs(ai)) {
		    kb[(a + b * kb_dim2 << 1) + 1] = (x + ai / ar * y) / (ai /
			     ar * ai + ar);
		    kb[(a + b * kb_dim2 << 1) + 2] = (y - ai / ar * x) / (ai /
			     ar * ai + ar);
		} else {
		    kb[(a + b * kb_dim2 << 1) + 1] = (ar / ai * x + y) / (ar /
			     ai * ar + ai);
		    kb[(a + b * kb_dim2 << 1) + 2] = (ar / ai * y - x) / (ar /
			     ai * ar + ai);
		}
	    }
/* L1020: */
	}
/* L1030: */
    }

} /* csyrdn_ */

