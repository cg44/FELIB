/* linqua.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int linqua_(xi, eta, geom, igeom, jgeom, nodel, sidnum, alen,
	 itest)
doublereal *xi, *eta, *geom;
integer *igeom, *jgeom, *nodel, *sidnum;
doublereal *alen;
integer *itest;
{
    /* Initialized data */

    static integer dimen = 2;
    static integer ider = 2;
    static integer ifun = 12;
    static integer ijac = 2;
    static integer jder = 12;
    static integer jjac = 2;
    static char srname[6+1] = "LINQUA";

    /* System generated locals */
    integer geom_dim1, geom_offset;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    extern /* Subroutine */ int quam4_(), quam8_();
    static integer n;
    extern /* Subroutine */ int quam12_();
    static integer jtest;
    extern integer errmes_();
    extern /* Subroutine */ int matmul_();
    static integer ierror;
    static doublereal jac[4]	/* was [2][2] */, der[24]	/* was [2][12]
	     */, fun[12];

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      LINQUA calculates a unit of length along the side of a */
/*      rectangular element (4, 8 or 12 noded) */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 2.0   1 Feb 1981 (CG) */
/*      Commented    24 Jul 1981 (CG) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      XI      local coordinate of point at which length is */
/*              required */
/*      ETA     local coordinate of point at which length is */
/*              required */
/*      GEOM    local coordinate array containing the global */
/*              coordiantes of each node on an element in the */
/*              local order */
/*      IGEOM   first dimension of array GEOM (.GE. NODEL) */
/*      JGEOM   second dimension of array GEOM (.GE. 2) */
/*      NODEL   number of nodes on the element */
/*      SIDNUM  the side number of the side of the element to */
/*              be used in calculating the length (.LE. 3) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      ALEN    the unit of length at the specified point */

/* ROUTINES called */
/*      ERRMES  QUAM4  QUAM8  QUAM12  MATMUL */


/*      SUBROUTINE LINQUA(XI,ETA,GEOM,IGEOM,JGEOM,NODEL,SIDNUM,ALEN,ITEST)
 */
/* ***********************************************************************
 */



    /* Parameter adjustments */
    geom_dim1 = *igeom;
    geom_offset = geom_dim1 + 1;
    geom -= geom_offset;

    /* Function Body */

/*     Parameter checking */

    jtest = *itest;
    if (jtest != -1) {
	ierror = 0;
	if (*nodel <= 0 || *sidnum <= 0) {
	    ierror = 1;
	}
	if (*igeom < *nodel || *jgeom < 2) {
	    ierror = 2;
	}
	if (*sidnum < 1 || *sidnum > 4) {
	    ierror = 3;
	}
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    n = *nodel / 4;

/*     Range checking on N (should 1, 2 or 3) */

    if (jtest != -1) {
	ierror = 0;
	if (n <= 0 || n >= 4) {
	    ierror = 4;
	}
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

    switch ((int)n) {
	case 1:  goto L1000;
	case 2:  goto L1010;
	case 3:  goto L1020;
    }

L1000:
    quam4_(fun, &ifun, der, &ider, &jder, xi, eta, itest);
    goto L1030;
L1010:
    quam8_(fun, &ifun, der, &ider, &jder, xi, eta, itest);
    goto L1030;
L1020:
    quam12_(fun, &ifun, der, &ider, &jder, xi, eta, itest);

L1030:
    matmul_(der, &ider, &jder, &geom[geom_offset], igeom, jgeom, jac, &ijac, &
	    jjac, &dimen, nodel, &dimen, itest);

/*     SELECT correct side of element */

    switch ((int)*sidnum) {
	case 1:  goto L1040;
	case 2:  goto L1050;
	case 3:  goto L1040;
	case 4:  goto L1050;
    }
L1040:
/* Computing 2nd power */
    d__1 = jac[1];
/* Computing 2nd power */
    d__2 = jac[3];
    *alen = sqrt(d__1 * d__1 + d__2 * d__2);
    return 0;
L1050:
/* Computing 2nd power */
    d__1 = jac[0];
/* Computing 2nd power */
    d__2 = jac[2];
    *alen = sqrt(d__1 * d__1 + d__2 * d__2);

} /* linqua_ */

