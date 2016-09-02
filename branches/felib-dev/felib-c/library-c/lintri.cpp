/* lintri.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int lintri_(xi, eta, geom, igeom, jgeom, nodel, sidnum, alen,
	 itest)
doublereal *xi, *eta, *geom;
integer *igeom, *jgeom, *nodel, *sidnum;
doublereal *alen;
integer *itest;
{
    /* Initialized data */

    static char srname[6+1] = "LINTRI";

    /* System generated locals */
    integer geom_dim1, geom_offset;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal dxdl, dydl;
    static integer n, jtest;
    static doublereal l1, l2, l3, t1, t2, t3, t4;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      LINTRI calculates a unit of length along the side of a */
/*      triangular element (3, 6 or 10 noded) */

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
/*              be used in calculating the length */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      ALEN    the unit of length at the specified point */

/* ROUTINES called */
/*      ERRMES */


/*      SUBROUTINE LINTRI(XI,ETA,GEOM,IGEOM,JGEOM,NODEL,SIDNUM,ALEN,ITEST)
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
	if (*sidnum < 1 || *sidnum > 3) {
	    ierror = 3;
	}
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    n = *nodel / 3;

/*     Range checking on N (should be 1, 2 or 3) */

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
	case 2:  goto L1040;
	case 3:  goto L1080;
    }
L1000:
    switch ((int)*sidnum) {
	case 1:  goto L1010;
	case 2:  goto L1020;
	case 3:  goto L1030;
    }

/*     Code for 3-node triangle */

L1010:
    dxdl = geom[geom_dim1 + 1] - geom[geom_dim1 + 2];
    dydl = geom[(geom_dim1 << 1) + 1] - geom[(geom_dim1 << 1) + 2];
    goto L1120;
L1020:
    dxdl = geom[geom_dim1 + 2] - geom[geom_dim1 + 3];
    dydl = geom[(geom_dim1 << 1) + 2] - geom[(geom_dim1 << 1) + 3];
    goto L1120;
L1030:
    dxdl = geom[geom_dim1 + 3] - geom[geom_dim1 + 1];
    dydl = geom[(geom_dim1 << 1) + 3] - geom[(geom_dim1 << 1) + 1];
    goto L1120;

/*     Code for 6-node triangle */

L1040:
    switch ((int)*sidnum) {
	case 1:  goto L1050;
	case 2:  goto L1060;
	case 3:  goto L1070;
    }
L1050:
    l1 = (*xi * 2. + 1.) * .33333333333333331;
    dxdl = (l1 * 4. - 1.) * geom[geom_dim1 + 1] + (1. - l1 * 2.) * 4. * geom[
	    geom_dim1 + 2] + (l1 * 4. - 3.) * geom[geom_dim1 + 3];
    dydl = (l1 * 4. - 1.) * geom[(geom_dim1 << 1) + 1] + (1. - l1 * 2.) * 4. *
	     geom[(geom_dim1 << 1) + 2] + (l1 * 4. - 3.) * geom[(geom_dim1 << 
	    1) + 3];
    goto L1120;
L1060:
    l2 = (1. - *xi - sqrt(3.) * *eta) * .33333333333333331;
    dxdl = (l2 * 4. - 1.) * geom[geom_dim1 + 3] + (1. - l2 * 2.) * 4. * geom[
	    geom_dim1 + 4] + (l2 * 4. - 3.) * geom[geom_dim1 + 5];
    dydl = (l2 * 4. - 1.) * geom[(geom_dim1 << 1) + 3] + (1. - l2 * 2.) * 4. *
	     geom[(geom_dim1 << 1) + 4] + (l2 * 4. - 3.) * geom[(geom_dim1 << 
	    1) + 5];
    goto L1120;
L1070:
    l3 = (1. - *xi + sqrt(3.) * *eta) * .33333333333333331;
    dxdl = (l3 * 4. - 3.) * geom[geom_dim1 + 1] + (l3 * 4. - 1.) * geom[
	    geom_dim1 + 5] + (1. - l3 * 2.) * 4. * geom[geom_dim1 + 6];
    dydl = (l3 * 4. - 3.) * geom[(geom_dim1 << 1) + 1] + (l3 * 4. - 1.) * 
	    geom[(geom_dim1 << 1) + 5] + (1. - l3 * 2.) * 4. * geom[(
	    geom_dim1 << 1) + 6];
    goto L1120;

/*     Code for 10-node triangle */

L1080:
    switch ((int)*sidnum) {
	case 1:  goto L1090;
	case 2:  goto L1100;
	case 3:  goto L1110;
    }
L1090:
    l1 = (*xi * 2. + 1.) * .33333333333333331;
    l2 = 1. - l1;
    t1 = (l1 * 27. * l1 - l2 * 14. + 2.) * .5;
    t2 = (l1 * 3. * l1 + l1 * 6. * l2 - l1 - l2) * 4.5;
    t3 = (l2 * 3. * l2 + l1 * 6. * l2 - l1 - l2) * 4.5;
    t4 = (l2 * 27 * l2 - l2 * 14. + 2.) * .5;
    dxdl = t1 * geom[geom_dim1 + 1] + t2 * geom[geom_dim1 + 2] + t3 * geom[
	    geom_dim1 + 4] + t4 * geom[geom_dim1 + 4];
    dydl = t1 * geom[(geom_dim1 << 1) + 1] + t2 * geom[(geom_dim1 << 1) + 2] 
	    + t3 * geom[(geom_dim1 << 1) + 4] + t4 * geom[(geom_dim1 << 1) + 
	    4];
    goto L1120;
L1100:
    l2 = (1. - *xi - sqrt(3.) * *eta) * .33333333333333331;
    l3 = 1. - l2;
    t1 = (l2 * 27. * l2 - l2 * 14. + 2.) * .5;
    t2 = (l2 * 3. * l2 + l2 * 6. * l3 - l2 - l3) * 4.5;
    t3 = (l3 * 3. * l3 + l2 * 6. * l3 - l2 - l3) * 4.5;
    t4 = (l3 * 27. * l3 - l3 * 14. + 2.) * .5;
    dxdl = t1 * geom[geom_dim1 + 4] + t2 * geom[geom_dim1 + 5] + t3 * geom[
	    geom_dim1 + 6] + t4 * geom[geom_dim1 + 7];
    dxdl = t1 * geom[(geom_dim1 << 1) + 4] + t2 * geom[(geom_dim1 << 1) + 5] 
	    + t3 * geom[(geom_dim1 << 1) + 6] + t4 * geom[(geom_dim1 << 1) + 
	    7];
    goto L1120;
L1110:
    l3 = (1. - *xi + sqrt(3.) * *eta) * .33333333333333331;
    l1 = 1. - l3;
    t1 = (l1 * 27. * l1 - l1 * 14. + 2.) * .5;
    t2 = (l3 * 27. * l3 - l3 * 14. + 2.) * .5;
    t3 = (l3 * 3. * l3 + l1 * 6. * l3 - l1 - l3) * 4.5;
    t4 = (l1 * 3. * l1 + l1 * 6. * l3 - l1 - l3) * 4.5;
    dxdl = t1 * geom[geom_dim1 + 1] + t2 * geom[geom_dim1 + 7] + t3 * geom[
	    geom_dim1 + 8] + t4 * geom[geom_dim1 + 9];
    dxdl = t1 * geom[(geom_dim1 << 1) + 1] + t2 * geom[(geom_dim1 << 1) + 7] 
	    + t3 * geom[(geom_dim1 << 1) + 8] + t4 * geom[(geom_dim1 << 1) + 
	    9];

/*     Do the calculation */

L1120:
/* Computing 2nd power */
    d__1 = dxdl;
/* Computing 2nd power */
    d__2 = dydl;
    *alen = sqrt(d__1 * d__1 + d__2 * d__2);

} /* lintri_ */

