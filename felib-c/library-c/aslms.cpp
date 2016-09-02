/* aslms.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int aslms_(sysm, isysm, elm, ielm, jelm, steer, isteer, 
	dofel, dofnod, size, itest)
doublereal *sysm;
integer *isysm;
doublereal *elm;
integer *ielm, *jelm, *steer, *isteer, *dofel, *dofnod;
doublereal *size;
integer *itest;
{
    /* Initialized data */

    static doublereal zero = 0.;
    static char srname[6+1] = "ASLMS ";

    /* System generated locals */
    integer elm_dim1, elm_offset, i__1, i__2;

    /* Local variables */
    static integer i, j;
    static doublereal x;
    static integer jtest;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      ASLMS assembles the contribution from an element to the */
/*      (diagonal) system matrix which is stored as a vector. */
/*      the 'lumped mass' approximation is assumed, the diagonal */
/*      elements of the element 'consistent mass' matrix being */
/*      used, suitably biassed to conserve the relevant quantity */
/*      (eg mass) on the element. */

/* HISTORY */

/*      Copyright (C) 1999 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    31 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      ISYSM   dimension of vector SYSM (.GE. total number of */
/*              unconstrained freedoms) */
/*      ELM     element mass matrix of dimension (IELEM, JELEM) */
/*              on entry, ELM(I,J) should contain the consistent */
/*              mass approximations for the element for */
/*              I=1(1)DOFEL */
/*      IELM    first dimension of ELM (.GE. DOFEL) */
/*      JELM    second dimension of ELM (.GE. DOFEL) */
/*      STEER   integer vector of length ISTEER containing */
/*              freedom numbers associating element matrix */
/*              contributions to system freedom numbers */
/*      ISTEER  length of vector STEER (.GE. DOFEL) */
/*      DOFEL   maximum number of degrees of freedom associated */
/*              with the element type */
/*      DOFNOD  number of degrees of freedom per node on the */
/*              element */
/*      SIZE    in two dimensions, area of the element */
/*              in three dimensions, volume of the element */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      SYSM    vector of length ISYSM containing the diagonal */
/*              elements of the (diagonal) system matrix */
/*      ELM     element 'mass' matrix, of dimension (IELM,JELM). */
/*              on exit, ELM(I,J) contains ZERO if I .NE. J, and */
/*              the calculated 'lumped mass' values if I=J */

/* ROUTINES called */
/*      ERRMES */

/*     SUBROUTINE ASLMS(SYSM,ISYSM,ELM,IELM,JELM,STEER,ISTEER,DOFEL, */
/*    *                 DOFNOD,SIZE,ITEST) */
/* ***********************************************************************
 */




    /* Parameter adjustments */
    --steer;
    elm_dim1 = *ielm;
    elm_offset = elm_dim1 + 1;
    elm -= elm_offset;
    --sysm;

    /* Function Body */

/*     Parameter checking */

    jtest = *itest;
    if (jtest != -1) {
	ierror = 0;
	if (*ielm < *dofel || *jelm < *dofel) {
	    ierror = 2;
	}
	if (*isteer < *dofel) {
	    ierror = 3;
	}
	if (*dofel == 0 || *dofnod == 0) {
	    ierror = 1;
	}
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main loops */

    x = 0.;
    i__1 = *dofel;
    for (i = 1; i <= i__1; ++i) {
	x += elm[i + i * elm_dim1];
/* L1000: */
    }
    x = *size / x * (doublereal) ((real) (*dofnod));
    i__1 = *dofel;
    for (i = 1; i <= i__1; ++i) {
	i__2 = *dofel;
	for (j = 1; j <= i__2; ++j) {
	    if (i == j) {
		elm[i + j * elm_dim1] *= x;
	    } else {
		elm[i + j * elm_dim1] = zero;
	    }
/* L1010: */
	}
/* L1020: */
    }
    i__1 = *dofel;
    for (i = 1; i <= i__1; ++i) {
	j = steer[i];
	if (j != 0) {

/*     Range checking on I and J */

	    if (jtest != -1) {
		ierror = 0;
		if (*isysm < j) {
		    ierror = 4;
		}
		*itest = errmes_(&jtest, &ierror, srname, 6L);
		if (*itest != 0) {
		    return 0;
		}
	    }

	    sysm[j] += elm[i + i * elm_dim1];
	}
/* L1030: */
    }

} /* aslms_ */

