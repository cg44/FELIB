/* b2c2.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int b2c2_(b, ib, jb, der, ider, jder, nodel, itest)
doublereal *b;
integer *ib, *jb;
doublereal *der;
integer *ider, *jder, *nodel, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "B2C2  ";

    /* System generated locals */
    integer b_dim1, b_offset, der_dim1, der_offset, i__1;

    /* Local variables */
    static integer k, l, m;
    extern integer errmes_();
    extern /* Subroutine */ int matnul_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      B2C2 forms the strain-displacement matrix for 2d plane */
/*      elasticity. */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (IMS) */
/*      Commented    12 Feb 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IB      first dimension of array B (.GE. 3) */
/*      JB      second dimension of B (.GE. 2*NODEL) */
/*      DER     DER(I,J) contains the derivative of the j'th */
/*              shape function with respect to the i'th global */
/*              coordinate */
/*      IDER    first dimension of DER (.GE. 2) */
/*      JDER    second dimension of DER (.GE. NODEL) */
/*      NODEL   number of nodes on the element */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      B       contains the values of the strain-displacement */
/*              matrix */

/* ROUTINES called */
/*      MATNUL  ERRMES */

/*      SUBROUTINE B2C2(B,IB,JB,DER,IDER,JDER,NODEL,ITEST) */
/* ***********************************************************************
 */



    /* Parameter adjustments */
    der_dim1 = *ider;
    der_offset = der_dim1 + 1;
    der -= der_offset;
    b_dim1 = *ib;
    b_offset = b_dim1 + 1;
    b -= b_offset;

    /* Function Body */

/*     Intialisation */

    k = 3;
    l = *nodel << 1;

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*ider < 2 || *jder < *nodel) {
	    ierror = 3;
	}
	if (*ib < 3 || *jb < l) {
	    ierror = 2;
	}
	if (*nodel <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Initialise matrix */

    matnul_(&b[b_offset], ib, jb, &k, &l, itest);

    i__1 = *nodel;
    for (m = 1; m <= i__1; ++m) {
	k = m << 1;
	l = k - 1;
	b[l * b_dim1 + 1] = der[m * der_dim1 + 1];
	b[k * b_dim1 + 3] = der[m * der_dim1 + 1];
	b[k * b_dim1 + 2] = der[m * der_dim1 + 2];
	b[l * b_dim1 + 3] = der[m * der_dim1 + 2];
/* L1000: */
    }

} /* b2c2_ */

