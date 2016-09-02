/* mvusb.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int mvusb_(a, ia, ja, v, iv, w, iw, n, hband, itest)
doublereal *a;
integer *ia, *ja;
doublereal *v;
integer *iv;
doublereal *w;
integer *iw, *n, *hband, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "MVUSB ";

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i, j, k;
    static doublereal x;
    static integer if_, is;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      MVUSB post-multiplies a real unsymmetric banded matrix stored as 
*/
/*      a rectangular array by A vector */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    14 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      A       array of dimension (IA, JA). Contains the */
/*              elements of the real unsymmetric band matrix */
/*              of order N and semi-bandwidth HBAND */
/*      IA      first dimension of A (.GE. N) */
/*      JA      second dimension of A (.GE. 2*HBAND-1) */
/*      V       vector of dimension IV */
/*      IV      dimension of vector V (.GE. N) */
/*      IW      dimension of vector W (.GE. N) */
/*      N       order of the real unsymmetric band matrix */
/*      HBAND   semi-bandwidth of the real unsymmetric band matrix */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      W       vector of dimension IW.  contains the result of */
/*              the operation W=A*V */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE MVUSB(A,IA,JA,V,IV,W,IW,N,HBAND,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    --w;
    --v;
    a_dim1 = *ia;
    a_offset = a_dim1 + 1;
    a -= a_offset;

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

/*     Main body */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	x = 0.;
/* Computing MAX */
	i__2 = *hband + 1 - i;
	is = max(i__2,1);
/* Computing MIN */
	i__2 = *n + *hband - i, i__3 = (*hband << 1) - 1;
	if_ = min(i__2,i__3);
	i__2 = if_;
	for (j = is; j <= i__2; ++j) {
	    k = i + j - *hband;
	    x += a[i + j * a_dim1] * v[k];
/* L1000: */
	}
	w[i] = x;
/* L1010: */
    }

} /* mvusb_ */

