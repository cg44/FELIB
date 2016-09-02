/* select.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int select_(v, iv, steer, isteer, n, w, iw, itest)
doublereal *v;
integer *iv, *steer, *isteer, *n;
doublereal *w;
integer *iw, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "SELECT";

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i, j, jtest;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      SELECT constructs an element value vector from a full system */
/*      vector using the steering vector */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 3.0  29 Oct 1985 (CG) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      V       vector of system values */
/*      IV      first dimension of vector V */
/*      STEER   element steering vector - contains freedom nos. */
/*      ISTEER  first dimension of vector STEER (.GE. N) */
/*      N       number of freedoms to be assembled */
/*      IW      first dimension of vector W (.GE. N) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      W       vector of element values */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE SELECT(V,IV,STEER,ISTEER,N,W,IW,ITEST) */
/* ***********************************************************************
 */


    /* Parameter adjustments */
    --w;
    --steer;
    --v;

    /* Function Body */

/*     Parameter checking */

    jtest = *itest;
    if (jtest != -1) {
	ierror = 0;
	if (*isteer < *n || *iw < *n) {
	    ierror = 2;
	}
	if (*n <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main loops */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	w[i] = 0.;
	j = (i__2 = steer[i], abs(i__2));
	if (j != 0) {

	    if (jtest != -1) {
		ierror = 0;
		if (j > *iv) {
		    ierror = 3;
		}
		*itest = errmes_(&jtest, &ierror, srname, 6L);
		if (*itest != 0) {
		    return 0;
		}
	    }

	    w[i] = v[j];
	}
/* L1000: */
    }

} /* select_ */

