/* dcsqua.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int dcsqua_(jacin, ijacin, jjacin, sidnum, cosin, icosin, 
	itest)
doublereal *jacin;
integer *ijacin, *jjacin, *sidnum;
doublereal *cosin;
integer *icosin, *itest;
{
    /* Initialized data */

    static integer dimen = 2;
    static integer iwork = 2;
    static char srname[6+1] = "DCSQUA";

    /* System generated locals */
    integer jacin_dim1, jacin_offset;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal amod, work[2];
    extern /* Subroutine */ int matvec_();
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      DCSQUA calculates the direction cosines of the outward normal */
/*      to the specified side of a quadrilateral element given the */
/*      jacobian inverse */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 2.0   1 Feb 1981 (CG) */
/*      Commented     1 Feb 1981 (CG) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      JACIN   array of dimension (IJACIN,JJACIN). contains the */
/*              inverse of the transformation jacobian. */
/*      IJACIN  first dimension of JACIN (.GE.2) */
/*      JJACIN  second dimension of JACIN (.GE.2) */
/*      SIDNUM  the side number for which the outward normal */
/*              is required (.LE.4) */
/*      ICOSIN  dimension of vector COSIN (.GE.2) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      COSIN   vector of dimension ICOSIN. contains the dir- */
/*              ection cosines of the outward normal */

/* ROUTINES called */
/*      ERRMES  MATVEC */

/*      SUBROUTINE DCSQUA(JACIN,IJACIN,JJACIN,SIDNUM,COSIN,ICOSIN,ITEST) 
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
	if (*ijacin < 2 || *jjacin < 2) {
	    ierror = 3;
	}
	if (*icosin < 2) {
	    ierror = 2;
	}
	if (*sidnum <= 0 || *sidnum > 4) {
	    ierror = 1;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    switch ((int)*sidnum) {
	case 1:  goto L1000;
	case 2:  goto L1010;
	case 3:  goto L1020;
	case 4:  goto L1030;
    }

/*     Side number 1 : xi = -1 */

L1000:
    work[0] = -1.;
    work[1] = 0.;
    goto L1040;

/*     Side number 2 : eta = +1 */

L1010:
    work[0] = 0.;
    work[1] = 1.;
    goto L1040;

/*     Side number 3 : xi = +1 */

L1020:
    work[0] = 1.;
    work[1] = 0.;
    goto L1040;

/*     Side number 4 : eta = -1 */

L1030:
    work[0] = 0.;
    work[1] = -1.;

/*     Calculate direction cosines */

L1040:
    matvec_(&jacin[jacin_offset], ijacin, jjacin, work, &iwork, &dimen, &
	    dimen, &cosin[1], icosin, itest);
    amod = sqrt(cosin[1] * cosin[1] + cosin[2] * cosin[2]);
    cosin[1] /= amod;
    cosin[2] /= amod;

} /* dcsqua_ */

