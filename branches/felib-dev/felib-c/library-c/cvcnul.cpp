/* cvcnul.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int cvcnul_(vec, ivec, n, itest)
doublereal *vec;
integer *ivec, *n, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "CVCNUL";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      CVCNUL sets the first N elements of a complex vector to zero */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 3.0  29 Oct 1984 (CRIE) */
/*      Commented    23 Oct 1985 (CG) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IVEC    length of vector VEC (.GE.N) */
/*      N       number of elements to be set to zero */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      VEC     vector of length IVEC. VEC(1,I)=0.0D0 and VEC(2,I)=0.0D0 
*/
/*              for I=1(1)N */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE CVCNUL(VEC,IVEC,N,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    vec -= 3;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*n > *ivec) {
	    ierror = 2;
	}
	if (*n <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);

/*     Main loops */

	if (ierror != 0) {
	    return 0;
	}
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	vec[(i << 1) + 1] = 0.;
	vec[(i << 1) + 2] = 0.;
/* L1000: */
    }

} /* cvcnul_ */

