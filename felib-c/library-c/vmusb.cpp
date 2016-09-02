/* vmusb.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int vmusb_(v, iv, a, ia, ja, w, iw, n, hband, itest)
doublereal *v;
integer *iv;
doublereal *a;
integer *ia, *ja;
doublereal *w;
integer *iw, *n, *hband, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "VMUSB ";

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer band, i, j;
    static doublereal x;
    static integer ij, ji;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      VMUSB pre-multiplies a real unsymmetric banded matrix stored as */

/*      a rectangular array by a vector */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    14 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      V       vector of dimension IV */
/*      IV      dimension of vector V (.GE. N) */
/*      A       array of dimension (IA, JA). Contains the */
/*              elements of the real unsymmetric BAND matrix */
/*              of order N and semi-bandwidth HBAND */
/*      IA      first dimension of A (.GE. N) */
/*      JA      second dimension of A (.GE. 2*HBAND-1) */
/*      IW      dimension of vector W (.GE. N) */
/*      N       order of the real unsymmetric BAND matrix */
/*      HBAND   semi-bandwidth of the real unsymmetric BAND matrix */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      W       vector of dimension IW. Contains the result of */
/*              the operation W=V*A */

/* ROUTINES called */
/*      ERRMES */

/*     SUBROUTINE VMUSB(V,IV,A,IA,JA,W,IW,N,HBAND,ITEST) */
/* ***********************************************************************
 */


    /* Parameter adjustments */
    --w;
    a_dim1 = *ia;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --v;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*iw < *n) {
	    ierror = 4;
	}
	if (*iv < *n) {
	    ierror = 3;
	}
	if (*ia < *n || *ja < (*hband << 1) - 1) {
	    ierror = 2;
	}
	if (*n <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }
    band = (*hband << 1) - 1;

/*     Main body */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	x = 0.;
	i__2 = band;
	for (j = 1; j <= i__2; ++j) {
	    ij = j - *hband + i;
	    ji = band + 1 - j;
	    if (ij >= 1 && ij <= *n) {
		x += a[ij + ji * a_dim1] * v[ij];
	    }
/* L1000: */
	}
	w[i] = x;
/* L1010: */
    }

} /* vmusb_ */

