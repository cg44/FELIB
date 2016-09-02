/* rasrhs.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int rasrhs_(rhs, irhs, value, ivalue, steer, isteer, dofel, 
	itest)
doublereal *rhs;
integer *irhs;
doublereal *value;
integer *ivalue, *steer, *isteer, *dofel, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "RASRHS";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k, l, jtest;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      RASRHS the routine adds into the right-hand side of a system */
/*      the values contianed in an element vector, thus forming */
/*      the right-hand side VALUE is real, RHS is complex. */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 3.0   1 Jul 1985 (CRIE) */
/*      Commented     1 Nov 1985 (CG) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      RHS     the right-hand side of the system */
/*              complex ( ordered pairs ) */
/*      IRHS    dimension of array RHS */
/*      VALUE   the element vector of the current element to */
/*              be added into the right-hand side */
/*      ivale   dimension of array VALUE */
/*      STEER   the steering vector containing the freedom */
/*              numbers of the freedoms associated with the */
/*              current element in the local order */
/*      ISTEER  dimension of array STEER */
/*      DOFEL   the maximum number of degrees of freedom on */
/*              an element of the current type */
/*      ITEST   error checking option */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE RASRHS(RHS,IRHS,VALUE,IVALUE,STEER,ISTEER,DOFEL,ITEST) 
*/
/* ***********************************************************************
 */

    /* Parameter adjustments */
    --steer;
    --value;
    rhs -= 3;

    /* Function Body */

/*     Parameter checking */

    jtest = *itest;
    if (*itest != -1) {
	ierror = 0;
	if (*isteer < *dofel) {
	    ierror = 3;
	}
	if (*ivalue < *dofel) {
	    ierror = 2;
	}
	if (*dofel <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main loops */

    i__1 = *dofel;
    for (k = 1; k <= i__1; ++k) {
	if (steer[k] != 0) {
	    l = steer[k];
	    if (jtest != -1) {
		ierror = 0;
		if (l > *irhs) {
		    ierror = 4;
		}
		*itest = errmes_(itest, &ierror, srname, 6L);
		if (*itest != 0) {
		    return 0;
		}
	    }
	    rhs[(l << 1) + 1] += value[k];
	}
/* L1000: */
    }

} /* rasrhs_ */

