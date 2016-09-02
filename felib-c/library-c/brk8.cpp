/* brk8.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int brk8_(fun, ifun, der, ider, jder, xi, eta, zeta, itest)
doublereal *fun;
integer *ifun;
doublereal *der;
integer *ider, *jder;
doublereal *xi, *eta, *zeta;
integer *itest;
{
    /* Initialized data */

    static char srname[6+1] = "BRK8  ";
    static doublereal half = .5;

    /* System generated locals */
    integer der_dim1, der_offset;

    /* Local variables */
    static doublereal etam, etap;
    extern doublereal veps_();
    static doublereal zetam, zetap, dummy;
    extern integer errmes_();
    static integer ierror;
    static doublereal val, xim, xip;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      BRK8 returns the values of shape functions and derivatives */
/*      at a specified point for an 8-noded brick element. The shape */
/*      function is continuous across element boundaries. */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    10 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IFUN    dimension of vector FUN (.GE.8) */
/*      IDER    first dimension of array DER (.GE.3) */
/*      JDER    second dimension of array DER (.GE.8) */
/*      XI      value of local coordinate at which function and */
/*              derivative values required */
/*      ETA     value of second local coordinate */
/*      ZETA    value of third local coordinate */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      FUN     real vector of dimension IFUN. FUN(I) contains */
/*              value of i'th shape function at (XI, ETA, ZETA) */
/*      DER     real array of dimensions (IDER, JDER). DER(I, J) */
/*              contains the derivative of the j'th shape */
/*              function with respect to the i'th coordinate at */
/*              the point (XI, ETA, ZETA) */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE BRK8(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    der_dim1 = *ider;
    der_offset = der_dim1 + 1;
    der -= der_offset;
    --fun;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*ifun < 8) {
	    ierror = 1;
	}
	if (*ider < 3 || *jder < 8) {
	    ierror = 2;
	}
	val = veps_(&dummy) + 1.;
	if (abs(*xi) > val || abs(*eta) > val || abs(*zeta) > val) {
	    ierror = 3;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main code */

    etam = (1. - *eta) * half;
    etap = (*eta + 1.) * half;
    xim = (1. - *xi) * half;
    xip = (*xi + 1.) * half;
    zetam = (1. - *zeta) * half;
    zetap = (*zeta + 1.) * half;

/*     Set shape functions */

    fun[1] = xim * etam * zetam;
    fun[2] = xim * etap * zetam;
    fun[3] = xip * etap * zetam;
    fun[4] = xip * etam * zetam;
    fun[5] = xim * etam * zetap;
    fun[6] = xim * etap * zetap;
    fun[7] = xip * etap * zetap;
    fun[8] = xip * etam * zetap;

/*     Set derivatves of shape functions */

    der[der_dim1 + 1] = -etam * zetam * half;
    der[(der_dim1 << 1) + 1] = -etap * zetam * half;
    der[der_dim1 * 3 + 1] = etap * zetam * half;
    der[(der_dim1 << 2) + 1] = etam * zetam * half;
    der[der_dim1 * 5 + 1] = -etam * zetap * half;
    der[der_dim1 * 6 + 1] = -etap * zetap * half;
    der[der_dim1 * 7 + 1] = etap * zetap * half;
    der[(der_dim1 << 3) + 1] = etam * zetap * half;
    der[der_dim1 + 2] = -xim * zetam * half;
    der[(der_dim1 << 1) + 2] = xim * zetam * half;
    der[der_dim1 * 3 + 2] = xip * zetam * half;
    der[(der_dim1 << 2) + 2] = -xip * zetam * half;
    der[der_dim1 * 5 + 2] = -xim * zetap * half;
    der[der_dim1 * 6 + 2] = xim * zetap * half;
    der[der_dim1 * 7 + 2] = xip * zetap * half;
    der[(der_dim1 << 3) + 2] = -xip * zetap * half;
    der[der_dim1 + 3] = -xim * etam * half;
    der[(der_dim1 << 1) + 3] = -xim * etap * half;
    der[der_dim1 * 3 + 3] = -xip * etap * half;
    der[(der_dim1 << 2) + 3] = -xip * etam * half;
    der[der_dim1 * 5 + 3] = xim * etam * half;
    der[der_dim1 * 6 + 3] = xim * etap * half;
    der[der_dim1 * 7 + 3] = xip * etap * half;
    der[(der_dim1 << 3) + 3] = xip * etam * half;

} /* brk8_ */

