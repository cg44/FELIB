/* dpsn.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int dpsn_(d, id, jd, e, nu, numss, itest)
doublereal *d;
integer *id, *jd;
doublereal *e, *nu;
integer *numss, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "DPSN  ";

    /* System generated locals */
    integer d_dim1, d_offset;

    /* Local variables */
    static doublereal nunu;
    static integer i, j;
    extern integer errmes_();
    extern /* Subroutine */ int matnul_();
    static integer ierror;
    static doublereal nu1;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      DPSN forms the 3 by 3 stress-strain matrix for use in plane */
/*      strain problems */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1980 (IMS) */
/*      Commented    14 Feb 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      ID      first dimension of array D (.GE. 3) */
/*      JD      second dimension of D (.GE. 3) */
/*      E       young's modulus */
/*      NU      poisson's ratio */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      D       stress-strain matrix for plane strain */
/*      NUMSS   order of stress-strain matrix (3) */

/* ROUTINES called */
/*      MATNUL  ERRMES */

/*      SUBROUTINE DPSN(D,ID,JD,E,NU,NUMSS,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    d_dim1 = *id;
    d_offset = d_dim1 + 1;
    d -= d_offset;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*nu < 0. || *nu >= .5) {
	    ierror = 2;
	}
	if (*id < *numss || *jd < *numss) {
	    ierror = 1;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    *numss = 3;

/*     Initialise matrix D */

    matnul_(&d[d_offset], id, jd, numss, numss, itest);

    nu1 = *nu / (1. - *nu);
    nunu = (.5 - *nu) / (1. - *nu);

    d[d_dim1 + 1] = 1.;
    d[(d_dim1 << 1) + 2] = 1.;
    d[d_dim1 * 3 + 3] = nunu;
    d[(d_dim1 << 1) + 1] = nu1;
    d[d_dim1 + 2] = nu1;

    for (i = 1; i <= 3; ++i) {
	for (j = 1; j <= 3; ++j) {
	    d[i + j * d_dim1] = d[i + j * d_dim1] * *e / ((*nu + 1.) * 2. * 
		    nunu);
/* L1000: */
	}
/* L1010: */
    }

} /* dpsn_ */

