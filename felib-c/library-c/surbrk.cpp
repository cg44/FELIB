/* surbrk.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int surbrk_(xi, eta, zeta, geom, igeom, jgeom, nodel, facnum,
	 uarea, itest)
doublereal *xi, *eta, *zeta, *geom;
integer *igeom, *jgeom, *nodel, *facnum;
doublereal *uarea;
integer *itest;
{
    /* Initialized data */

    static integer dimen = 3;
    static integer ider = 3;
    static integer ifun = 32;
    static integer ijac = 3;
    static integer jder = 32;
    static integer jjac = 3;
    static char srname[6+1] = "SURBRK";

    /* System generated locals */
    integer geom_dim1, geom_offset;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    extern /* Subroutine */ int brk20_(), brk32_();
    static integer i, j, n, jtest;
    static doublereal g1, g2, g3;
    extern integer errmes_();
    extern /* Subroutine */ int matmul_();
    static integer ierror;
    static doublereal jac[9]	/* was [3][3] */, der[96]	/* was [3][32]
	     */, fun[32];
    extern /* Subroutine */ int brk8_();

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      SURBRK calculates a unit of area on the face of a */
/*      brick element (8, 20 or 32 noded) */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 2.0  1  Feb 1981 (CG) */
/*      Commented    24 Jul 1981 (CG) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      XI      local coordinate of point at which area is */
/*              required */
/*      ETA     local coordinate of point at which area is */
/*              required */
/*      ZETA    local coordinate of point at which area is */
/*              required */
/*      GEOM    local coordinate array containing the global */
/*              coordiantes of each node on an element in the */
/*              local order */
/*      IGEOM   first dimension of array GEOM (.GE. NODEL) */
/*      JGEOM   second dimension of array GEOM (.GE. 3) */
/*      NODEL   number of nodes on the element */
/*      FACNUM  the face number of the face of the element to */
/*              be used in calculating the area (.LE. 6) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      UAREA    the unit of area at the specified point */

/* ROUTINES called */
/*      ERRMES  BRK8  BRK20  BRK32  MATMUL */

/*     SUBROUTINE SURBRK(XI,ETA,ZETA,GEOM,IGEOM,JGEOM,NODEL,FACNUM,UAREA, 
*/
/*     *                  ITEST) */
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
	if (*nodel <= 0 || *facnum <= 0) {
	    ierror = 1;
	}
	if (*igeom < *nodel || *jgeom < 3) {
	    ierror = 2;
	}
	if (*facnum < 1 || *facnum > 6) {
	    ierror = 3;
	}
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    n = *nodel / 8;

/*     Range checking on N */

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

/*     Get shape function derivatives */

    switch ((int)n) {
	case 1:  goto L1020;
	case 2:  goto L1000;
	case 3:  goto L1010;
    }
    if (jtest == -1) {
	goto L1020;
    } else {
	ierror = 4;
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	return 0;
    }
L1000:
    brk20_(fun, &ifun, der, &ider, &jder, xi, eta, zeta, itest);
    goto L1030;
L1010:
    brk32_(fun, &ifun, der, &ider, &jder, xi, eta, zeta, itest);
    goto L1030;

L1020:
    brk8_(fun, &ifun, der, &ider, &jder, xi, eta, zeta, itest);

/*     Calculate global derivatives */

L1030:
    matmul_(der, &ider, &jder, &geom[geom_offset], igeom, jgeom, jac, &ijac, &
	    jjac, &dimen, nodel, &dimen, itest);

    switch ((int)*facnum) {
	case 1:  goto L1040;
	case 2:  goto L1040;
	case 3:  goto L1050;
	case 4:  goto L1060;
	case 5:  goto L1050;
	case 6:  goto L1060;
    }
    goto L1050;

/*     ZETA = constant */

L1040:
    i = 1;
    j = 2;
    goto L1070;

/*     XI = constant */

L1050:
    i = 2;
    j = 3;
    goto L1070;

/*     ETA = constant */

L1060:
    i = 1;
    j = 3;

L1070:
    g1 = jac[i + 2] * jac[j + 5] - jac[i + 5] * jac[j + 2];
    g2 = jac[i + 5] * jac[j - 1] - jac[i - 1] * jac[j + 5];
    g3 = jac[i - 1] * jac[j + 2] - jac[i + 2] * jac[j - 1];
/* Computing 2nd power */
    d__1 = g1;
/* Computing 2nd power */
    d__2 = g2;
/* Computing 2nd power */
    d__3 = g3;
    *uarea = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);

} /* surbrk_ */

