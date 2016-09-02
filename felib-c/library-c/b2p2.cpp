/* b2p2.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int b2p2_(b, ib, jb, der, ider, jder, fun, ifun, coord, 
	icoord, jcoord, nodel, itest)
doublereal *b;
integer *ib, *jb;
doublereal *der;
integer *ider, *jder;
doublereal *fun;
integer *ifun;
doublereal *coord;
integer *icoord, *jcoord, *nodel, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "B2P2  ";

    /* System generated locals */
    integer b_dim1, b_offset, coord_dim1, coord_offset, der_dim1, der_offset, 
	    i__1;

    /* Local variables */
    extern doublereal vtol_();
    static integer k, l, m;
    static doublereal x;
    static integer jtest;
    extern integer errmes_();
    extern /* Subroutine */ int matnul_();
    static integer ierror;
    static doublereal tol, sum;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      B2P2 forms the strain-displacement matrix for axisymmetric */
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
/*      FUN     FUN(I) conatins the value of the i'th shape */
/*              function at the point under consideration */
/*      IFUN    first dimension of FUN (.GE. 4) */
/*      COORD   COORD(I,J) contains the j'th global coordinate */
/*              of the i'th node */
/*      ICOORD  first dimension of COORD (.GE. number of nodes */
/*              in the mesh) */
/*      JCOORD  second dimension of COORD (.GE. dimensionality */
/*              of problem) */
/*      NODEL   number of nodes on element */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      B       conatins values of strain-displacement matrix */

/* ROUTINES called */
/*      MATNUL  ERRMES  VTOL */

/*      SUBROUTINE B2P2(B,IB,JB,DER,IDER,JDER,FUN,IFUN,COORD,ICOORD, */
/*     *                JCOORD,NODEL,ITEST) */
/* ***********************************************************************
 */



    /* Parameter adjustments */
    coord_dim1 = *icoord;
    coord_offset = coord_dim1 + 1;
    coord -= coord_offset;
    --fun;
    der_dim1 = *ider;
    der_offset = der_dim1 + 1;
    der -= der_offset;
    b_dim1 = *ib;
    b_offset = b_dim1 + 1;
    b -= b_offset;

    /* Function Body */

/*     Intialisation */

    k = 4;
    l = *nodel << 1;
    tol = vtol_(&x);

/*     Parameter checking */

    jtest = *itest;
    if (jtest != -1) {
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
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main loops */

    matnul_(&b[b_offset], ib, jb, &k, &l, itest);

/*     Calculate average radius */

    sum = 0.;
    i__1 = *nodel;
    for (k = 1; k <= i__1; ++k) {
	sum += fun[k] * coord[k + coord_dim1];
/* L1000: */
    }

/*     Check value of SUM rbar */

    if (jtest != -1) {
	ierror = 0;
	if (abs(sum) < tol) {
	    ierror = 4;
	}
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

    i__1 = *nodel;
    for (m = 1; m <= i__1; ++m) {
	k = m << 1;
	l = k - 1;
	b[l * b_dim1 + 1] = der[m * der_dim1 + 1];
	b[k * b_dim1 + 3] = der[m * der_dim1 + 1];
	b[k * b_dim1 + 2] = der[m * der_dim1 + 2];
	b[l * b_dim1 + 3] = der[m * der_dim1 + 2];
	b[l * b_dim1 + 4] = fun[m] / sum;
/* L1010: */
    }

} /* b2p2_ */

