/* norm.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


doublereal norm_(rhs, irhs, totdof, itest)
doublereal *rhs;
integer *irhs, *totdof, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "NORM  ";

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer i;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      NORM computes an L2 norm of a vector for use */
/*      in terminating non-linear iterations */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 3.0  29 Jun 1986 (CJH) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      RHS     vector containing values */
/*      IRHS    dimension of vector RHS (IRHS .GE. TOTDOF) */
/*      TOTDOF  number of entries in RHS */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      NORM    the value of the NORM (function value) */

/* ROUTINES called */
/*      ERRMES */


/*     DOUBLE PRECISION FUNCTION NORM(RHS,IRHS,TOTDOF,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    --rhs;

    /* Function Body */

/*     Parameter checking */

    ret_val = 0.;

    if (*itest != -1) {
	ierror = 0;
	if (*totdof > *irhs) {
	    ierror = 1;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return ret_val;
	}
    }

/*     Compute NORM */

    i__1 = *totdof;
    for (i = 1; i <= i__1; ++i) {
/* Computing 2nd power */
	d__1 = rhs[i];
	ret_val += d__1 * d__1;
/* L1000: */
    }
    ret_val = sqrt(ret_val);

    return ret_val;
} /* norm_ */

