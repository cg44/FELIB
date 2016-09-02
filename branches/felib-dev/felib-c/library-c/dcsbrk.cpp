/* dcsbrk.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int dcsbrk_(jacin, ijacin, jjacin, facnum, cosin, icosin, 
	itest)
doublereal *jacin;
integer *ijacin, *jjacin, *facnum;
doublereal *cosin;
integer *icosin, *itest;
{
    /* Initialized data */

    static integer dimen = 3;
    static integer iwork = 3;
    static char srname[6+1] = "DCSBRK";

    /* System generated locals */
    integer jacin_dim1, jacin_offset;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal amod, work[3];
    extern /* Subroutine */ int matvec_();
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      DCSBRK calculates the direction cosines of the outward normal */
/*      to the specified face of a hexahedral element given the */
/*      jacobian inverse */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 2.0   1 Feb 1981 (CG) */
/*      Commented     1 Feb 1981 (CG) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      JACIN   array of dimension (IJACIN, JJACIN). contains the */
/*              inverse of the transformation jacobian. */
/*      IJACIN  first dimension of JACIN (.GE. 3) */
/*      JJACIN  second dimension of JACIN (.GE. 3) */
/*      FACNUM  the side number for which the outward normal */
/*              is required (.LE. 6) */
/*      ICOSIN  dimension of vector COSIN (.GE. 3) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      COSIN   vector of dimension ICOSIN. contains the dir- */
/*              ection cosines of the outward normal */

/* ROUTINES called */
/*      ERRMES  MATVEC */

/*      SUBROUTINE DCSBRK(JACIN,IJACIN,JJACIN,FACNUM,COSIN,ICOSIN,ITEST) 
*/
/* ***********************************************************************
 */

    /* Parameter adjustments */
    --cosin;
    jacin_dim1 = *ijacin;
    jacin_offset = jacin_dim1 + 1;
    jacin -= jacin_offset;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*ijacin < 3 || *jjacin < 3) {
	    ierror = 3;
	}
	if (*icosin < 3) {
	    ierror = 2;
	}
	if (*facnum <= 0 || *facnum > 6) {
	    ierror = 1;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    switch ((int)*facnum) {
	case 1:  goto L1000;
	case 2:  goto L1010;
	case 3:  goto L1020;
	case 4:  goto L1030;
	case 5:  goto L1040;
	case 6:  goto L1050;
    }
    goto L1020;

/*     Face number 1 : zeta = -1 */

L1000:
    work[0] = 0.;
    work[1] = 0.;
    work[2] = -1.;
    goto L1060;

/*     Face number 2 : zeta = +1 */

L1010:
    work[0] = 0.;
    work[1] = 0.;
    work[2] = 1.;
    goto L1060;

/*     Face number 3 : xi = -1 */

L1020:
    work[0] = -1.;
    work[1] = 0.;
    work[2] = 0.;
    goto L1060;

/*     Face number 4 : eta = +1 */

L1030:
    work[0] = 0.;
    work[1] = 1.;
    work[2] = 0.;
    goto L1060;

/*     Face number 5 : xi = +1 */

L1040:
    work[0] = 1.;
    work[1] = 0.;
    work[2] = 0.;
    goto L1060;

/*     Face number 6 : eta = -1 */

L1050:
    work[0] = 0.;
    work[1] = -1.;
    work[2] = 0.;

/*     Calculate direction cosines */

L1060:
    matvec_(&jacin[jacin_offset], ijacin, jjacin, work, &iwork, &dimen, &
	    dimen, &cosin[1], icosin, itest);
    amod = sqrt(cosin[1] * cosin[1] + cosin[2] * cosin[2] + cosin[3] * cosin[
	    3]);
    cosin[1] /= amod;
    cosin[2] /= amod;
    cosin[3] /= amod;

} /* dcsbrk_ */

