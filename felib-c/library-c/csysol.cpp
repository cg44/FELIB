/* csysol.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int csysol_(kb, ikb, jkb, loads, iloads, n, hband, itest)
doublereal *kb;
integer *ikb, *jkb;
doublereal *loads;
integer *iloads, *n, *hband, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "CSYSOL";

    /* System generated locals */
    integer kb_dim2, kb_offset;

    /* Local variables */
    static integer jtest;
    extern integer errmes_();
    static integer ierror;
    extern /* Subroutine */ int csyrdn_(), csysub_();

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      CSYSOL solves a set of complex symmetric banded equations with a 
*/
/*      single right hand side using a symmetric decomposition. Only */
/*      the lower band and diagonal are stored in a rectangular array KB. 
*/

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 3.0  29 Oct 1979 (CRIE) */
/*      Commented    12 Feb 1985 (CG) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      KB      on entry contains lower half of complex symmetric */
/*              band matrix stored as a rectangular array */
/*      IKB     first dimension of KB (.GE. N) */
/*      JKB     second dimension of KB (.GE. HBAND) */
/*      LOADS   contains elements of right hand side */
/*      ILOADS  dimension of LOADS (.GE. N) */
/*      N       order of matrix KB */
/*      HBAND   semi-bandwidth of KB (includes diagonal) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      KB      on exit, contains lower triangular reduced */
/*      LOADS   matrix solution vector */

/* ROUTINES called */
/*      ERRMES CSYRDN CSYSUB */

/*      SUBROUTINE CSYSOL(KB,IKB,JKB,LOADS,ILOADS,N,HBAND,ITEST) */
/************************************************************************
***/



    /* Parameter adjustments */
    loads -= 3;
    kb_dim2 = *ikb;
    kb_offset = (kb_dim2 + 1 << 1) + 1;
    kb -= kb_offset;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*iloads < *n) {
	    ierror = 3;
	}
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

/*     Main body */

    jtest = 1;
    csyrdn_(&kb[kb_offset], ikb, jkb, n, hband, &jtest);

/*     Parameter checking */

    if (*itest != -1) {
	ierror = jtest;
	if (jtest == 3) {
	    ierror = 4;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

    csysub_(&kb[kb_offset], ikb, jkb, &loads[3], iloads, n, hband, &jtest);

/*     Parameter checking */

    if (*itest != -1) {
	ierror = jtest;
	if (jtest == 3) {
	    ierror = 4;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest == 0) {
	    return 0;
	}
    }


} /* csysol_ */

