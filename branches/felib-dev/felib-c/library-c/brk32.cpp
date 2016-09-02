/* brk32.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int brk32_(fun, ifun, der, ider, jder, xi, eta, zeta, itest)
doublereal *fun;
integer *ifun;
doublereal *der;
integer *ider, *jder;
doublereal *xi, *eta, *zeta;
integer *itest;
{
    /* Initialized data */

    static char srname[6+1] = "BRK32 ";

    /* System generated locals */
    integer der_dim1, der_offset;

    /* Local variables */
    extern doublereal veps_();
    static doublereal mthrd, pthrd, dummy, m1, p1;
    extern integer errmes_();
    static integer ierror;
    static doublereal val;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      BRK32 calculates the values of the shape functions and their */
/*      derivatives at a point for a 32-noded brick element. The shape */
/*      function is continuous across element boundaries */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    10 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IFUN    dimension of vector FUN (.GE.32) */
/*      IDER    first dimension of array DER (.GE.3) */
/*      JDER    second dimension of array DER (.GE.32) */
/*      XI      first local coordinate */
/*      ETA     second local coordinate */
/*      ZETA    third local coordinate */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      FUN     vector containing shape functions.  FUN(I) */
/*              contains the value of the i'th shape function at */
/*              (XI, ETA, ZETA) */
/*      DER     array containing the derivatives of the shape */
/*              functions.  DER(I, J) contains the value of the */
/*              derivative of the j'th shape function with */
/*              respect to the i'th coordinate at (XI, ETA, ZETA) */

/* ROUTINES called */
/*      ERRMES */


/*      SUBROUTINE BRK32(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    der_dim1 = *ider;
    der_offset = der_dim1 + 1;
    der -= der_offset;
    --fun;

    /* Function Body */

/*     Statement functions */


/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*ifun < 32) {
	    ierror = 1;
	}
	if (*ider < 3 || *jder < 32) {
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

/*     Body of code */

    p1 = 1.;
    m1 = -1.;
    pthrd = .33333333333333331;
    mthrd = -.33333333333333331;

/*     Set shape functions */

    fun[1] = (*xi * m1 + 1.) * .015625 * (*eta * m1 + 1.) * (*zeta * m1 + 1.) 
	    * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.);
    fun[2] = (1. - *eta * *eta) * .140625 * (*eta * 9. * mthrd + 1.) * (*xi * 
	    m1 + 1.) * (*zeta * m1 + 1.);
    fun[3] = (1. - *eta * *eta) * .140625 * (*eta * 9. * pthrd + 1.) * (*xi * 
	    m1 + 1.) * (*zeta * m1 + 1.);
    fun[4] = (*xi * m1 + 1.) * .015625 * (*eta * p1 + 1.) * (*zeta * m1 + 1.) 
	    * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.);
    fun[5] = (1. - *xi * *xi) * .140625 * (*xi * 9. * mthrd + 1.) * (*eta * 
	    p1 + 1.) * (*zeta * m1 + 1.);
    fun[6] = (1. - *xi * *xi) * .140625 * (*xi * 9. * pthrd + 1.) * (*eta * 
	    p1 + 1.) * (*zeta * m1 + 1.);
    fun[7] = (*xi * p1 + 1.) * .015625 * (*eta * p1 + 1.) * (*zeta * m1 + 1.) 
	    * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.);
    fun[8] = (1. - *eta * *eta) * .140625 * (*eta * 9. * pthrd + 1.) * (*xi * 
	    p1 + 1.) * (*zeta * m1 + 1.);
    fun[9] = (1. - *eta * *eta) * .140625 * (*eta * 9. * mthrd + 1.) * (*xi * 
	    p1 + 1.) * (*zeta * m1 + 1.);
    fun[10] = (*xi * p1 + 1.) * .015625 * (*eta * m1 + 1.) * (*zeta * m1 + 1.)
	     * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.);
    fun[11] = (1. - *xi * *xi) * .140625 * (*xi * 9. * pthrd + 1.) * (*eta * 
	    m1 + 1.) * (*zeta * m1 + 1.);
    fun[12] = (1. - *xi * *xi) * .140625 * (*xi * 9. * mthrd + 1.) * (*eta * 
	    m1 + 1.) * (*zeta * m1 + 1.);
    fun[13] = (1. - *zeta * *zeta) * .140625 * (*zeta * 9. * mthrd + 1.) * (*
	    xi * m1 + 1.) * (*eta * m1 + 1.);
    fun[14] = (1. - *zeta * *zeta) * .140625 * (*zeta * 9. * mthrd + 1.) * (*
	    xi * m1 + 1.) * (*eta * p1 + 1.);
    fun[15] = (1. - *zeta * *zeta) * .140625 * (*zeta * 9. * mthrd + 1.) * (*
	    xi * p1 + 1.) * (*eta * p1 + 1.);
    fun[16] = (1. - *zeta * *zeta) * .140625 * (*zeta * 9. * mthrd + 1.) * (*
	    xi * p1 + 1.) * (*eta * m1 + 1.);
    fun[17] = (1. - *zeta * *zeta) * .140625 * (*zeta * 9. * pthrd + 1.) * (*
	    xi * m1 + 1.) * (*eta * m1 + 1.);
    fun[18] = (1. - *zeta * *zeta) * .140625 * (*zeta * 9. * pthrd + 1.) * (*
	    xi * m1 + 1.) * (*eta * p1 + 1.);
    fun[19] = (1. - *zeta * *zeta) * .140625 * (*zeta * 9. * pthrd + 1.) * (*
	    xi * p1 + 1.) * (*eta * p1 + 1.);
    fun[20] = (1. - *zeta * *zeta) * .140625 * (*zeta * 9. * pthrd + 1.) * (*
	    xi * p1 + 1.) * (*eta * m1 + 1.);
    fun[21] = (*xi * m1 + 1.) * .015625 * (*eta * m1 + 1.) * (*zeta * p1 + 1.)
	     * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.);
    fun[22] = (1. - *eta * *eta) * .140625 * (*eta * 9. * mthrd + 1.) * (*xi *
	     m1 + 1.) * (*zeta * p1 + 1.);
    fun[23] = (1. - *eta * *eta) * .140625 * (*eta * 9. * pthrd + 1.) * (*xi *
	     m1 + 1.) * (*zeta * p1 + 1.);
    fun[24] = (*xi * m1 + 1.) * .015625 * (*eta * p1 + 1.) * (*zeta * p1 + 1.)
	     * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.);
    fun[25] = (1. - *xi * *xi) * .140625 * (*xi * 9. * mthrd + 1.) * (*eta * 
	    p1 + 1.) * (*zeta * p1 + 1.);
    fun[26] = (1. - *xi * *xi) * .140625 * (*xi * 9. * pthrd + 1.) * (*eta * 
	    p1 + 1.) * (*zeta * p1 + 1.);
    fun[27] = (*xi * p1 + 1.) * .015625 * (*eta * p1 + 1.) * (*zeta * p1 + 1.)
	     * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.);
    fun[28] = (1. - *eta * *eta) * .140625 * (*eta * 9. * pthrd + 1.) * (*xi *
	     p1 + 1.) * (*zeta * p1 + 1.);
    fun[29] = (1. - *eta * *eta) * .140625 * (*eta * 9. * mthrd + 1.) * (*xi *
	     p1 + 1.) * (*zeta * p1 + 1.);
    fun[30] = (*xi * p1 + 1.) * .015625 * (*eta * m1 + 1.) * (*zeta * p1 + 1.)
	     * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.);
    fun[31] = (1. - *xi * *xi) * .140625 * (*xi * 9. * pthrd + 1.) * (*eta * 
	    m1 + 1.) * (*zeta * p1 + 1.);
    fun[32] = (1. - *xi * *xi) * .140625 * (*xi * 9. * mthrd + 1.) * (*eta * 
	    m1 + 1.) * (*zeta * p1 + 1.);

/*     Set derivatives */

    der[der_dim1 + 1] = (*eta * m1 + 1.) * .015625 * (*zeta * m1 + 1.) * (m1 *
	     ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *xi * 
	    18. * (*xi * m1 + 1.));
    der[der_dim1 + 2] = (*xi * m1 + 1.) * .015625 * (*zeta * m1 + 1.) * (m1 * 
	    ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *eta * 
	    18. * (*eta * m1 + 1.));
    der[der_dim1 + 3] = (*xi * m1 + 1.) * .015625 * (*eta * m1 + 1.) * (m1 * (
	    (*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *zeta * 
	    18. * (*zeta * m1 + 1.));
    der[(der_dim1 << 1) + 1] = m1 * .140625 * (1. - *eta * *eta) * (*eta * 9. 
	    * mthrd + 1.) * (*zeta * m1 + 1.);
    der[(der_dim1 << 1) + 2] = (*xi * m1 + 1.) * .140625 * (*zeta * m1 + 1.) *
	     (mthrd * 9. * (1. - *eta * *eta) - *eta * 2. * (*eta * 9. * 
	    mthrd + 1.));
    der[(der_dim1 << 1) + 3] = m1 * .140625 * (1. - *eta * *eta) * (*eta * 9. 
	    * mthrd + 1.) * (*xi * m1 + 1.);
    der[der_dim1 * 3 + 1] = m1 * .140625 * (1. - *eta * *eta) * (*eta * 9. * 
	    pthrd + 1.) * (*zeta * m1 + 1.);
    der[der_dim1 * 3 + 2] = (*xi * m1 + 1.) * .140625 * (*zeta * m1 + 1.) * (
	    pthrd * 9. * (1. - *eta * *eta) - *eta * 2. * (*eta * 9. * pthrd 
	    + 1.));
    der[der_dim1 * 3 + 3] = m1 * .140625 * (1. - *eta * *eta) * (*eta * 9. * 
	    pthrd + 1.) * (*xi * m1 + 1.);
    der[(der_dim1 << 2) + 1] = (*eta * p1 + 1.) * .015625 * (*zeta * m1 + 1.) 
	    * (m1 * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *
	    xi * 18. * (*xi * m1 + 1.));
    der[(der_dim1 << 2) + 2] = (*xi * m1 + 1.) * .015625 * (*zeta * m1 + 1.) *
	     (p1 * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *
	    eta * 18. * (*eta * p1 + 1.));
    der[(der_dim1 << 2) + 3] = (*xi * m1 + 1.) * .015625 * (*eta * p1 + 1.) * 
	    (m1 * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *
	    zeta * 18. * (*zeta * m1 + 1.));
    der[der_dim1 * 5 + 1] = (*eta * p1 + 1.) * .140625 * (*zeta * m1 + 1.) * (
	    mthrd * 9. * (1. - *xi * *xi) - *xi * 2. * (*xi * 9. * mthrd + 1.)
	    );
    der[der_dim1 * 5 + 2] = p1 * .140625 * (1. - *xi * *xi) * (*xi * 9. * 
	    mthrd + 1.) * (*zeta * m1 + 1.);
    der[der_dim1 * 5 + 3] = m1 * .140625 * (1. - *xi * *xi) * (*xi * 9. * 
	    mthrd + 1.) * (*eta * p1 + 1.);
    der[der_dim1 * 6 + 1] = (*eta * p1 + 1.) * .140625 * (*zeta * m1 + 1.) * (
	    pthrd * 9. * (1. - *xi * *xi) - *xi * 2. * (*xi * 9. * pthrd + 1.)
	    );
    der[der_dim1 * 6 + 2] = p1 * .140625 * (1. - *xi * *xi) * (*xi * 9. * 
	    pthrd + 1.) * (*zeta * m1 + 1.);
    der[der_dim1 * 6 + 3] = m1 * .140625 * (1. - *xi * *xi) * (*xi * 9. * 
	    pthrd + 1.) * (*eta * p1 + 1.);
    der[der_dim1 * 7 + 1] = (*eta * p1 + 1.) * .015625 * (*zeta * m1 + 1.) * (
	    p1 * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *xi 
	    * 18. * (*xi * p1 + 1.));
    der[der_dim1 * 7 + 2] = (*xi * p1 + 1.) * .015625 * (*zeta * m1 + 1.) * (
	    p1 * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *
	    eta * 18. * (*eta * p1 + 1.));
    der[der_dim1 * 7 + 3] = (*xi * p1 + 1.) * .015625 * (*eta * p1 + 1.) * (
	    m1 * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *
	    zeta * 18. * (*zeta * m1 + 1.));
    der[(der_dim1 << 3) + 1] = p1 * .140625 * (1. - *eta * *eta) * (*eta * 9. 
	    * pthrd + 1.) * (*zeta * m1 + 1.);
    der[(der_dim1 << 3) + 2] = (*xi * p1 + 1.) * .140625 * (*zeta * m1 + 1.) *
	     (pthrd * 9. * (1. - *eta * *eta) - *eta * 2. * (*eta * 9. * 
	    pthrd + 1.));
    der[(der_dim1 << 3) + 3] = m1 * .140625 * (1. - *eta * *eta) * (*eta * 9. 
	    * pthrd + 1.) * (*xi * p1 + 1.);
    der[der_dim1 * 9 + 1] = p1 * .140625 * (1. - *eta * *eta) * (*eta * 9. * 
	    mthrd + 1.) * (*zeta * m1 + 1.);
    der[der_dim1 * 9 + 2] = (*xi * p1 + 1.) * .140625 * (*zeta * m1 + 1.) * (
	    mthrd * 9. * (1. - *eta * *eta) - *eta * 2. * (*eta * 9. * mthrd 
	    + 1.));
    der[der_dim1 * 9 + 3] = m1 * .140625 * (1. - *eta * *eta) * (*eta * 9. * 
	    mthrd + 1.) * (*xi * p1 + 1.);
    der[der_dim1 * 10 + 1] = (*eta * m1 + 1.) * .015625 * (*zeta * m1 + 1.) * 
	    (p1 * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *
	    xi * 18. * (*xi * p1 + 1.));
    der[der_dim1 * 10 + 2] = (*xi * p1 + 1.) * .015625 * (*zeta * m1 + 1.) * (
	    m1 * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *
	    eta * 18. * (*eta * m1 + 1.));
    der[der_dim1 * 10 + 3] = (*xi * p1 + 1.) * .015625 * (*eta * m1 + 1.) * (
	    m1 * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *
	    zeta * 18. * (*zeta * m1 + 1.));
    der[der_dim1 * 11 + 1] = (*eta * m1 + 1.) * .140625 * (*zeta * m1 + 1.) * 
	    (pthrd * 9. * (1. - *xi * *xi) - *xi * 2. * (*xi * 9. * pthrd + 
	    1.));
    der[der_dim1 * 11 + 2] = m1 * .140625 * (1. - *xi * *xi) * (*xi * 9. * 
	    pthrd + 1.) * (*zeta * m1 + 1.);
    der[der_dim1 * 11 + 3] = m1 * .140625 * (1. - *xi * *xi) * (*xi * 9. * 
	    pthrd + 1.) * (*eta * m1 + 1.);
    der[der_dim1 * 12 + 1] = (*eta * m1 + 1.) * .140625 * (*zeta * m1 + 1.) * 
	    (mthrd * 9. * (1. - *xi * *xi) - *xi * 2. * (*xi * 9. * mthrd + 
	    1.));
    der[der_dim1 * 12 + 2] = m1 * .140625 * (1. - *xi * *xi) * (*xi * 9. * 
	    mthrd + 1.) * (*zeta * m1 + 1.);
    der[der_dim1 * 12 + 3] = m1 * .140625 * (1. - *xi * *xi) * (*xi * 9. * 
	    mthrd + 1.) * (*eta * m1 + 1.);
    der[der_dim1 * 13 + 1] = m1 * .140625 * (1. - *zeta * *zeta) * (*zeta * 
	    9. * mthrd + 1.) * (*eta * m1 + 1.);
    der[der_dim1 * 13 + 2] = m1 * .140625 * (1. - *zeta * *zeta) * (*zeta * 
	    9. * mthrd + 1.) * (*xi * m1 + 1.);
    der[der_dim1 * 13 + 3] = (*xi * m1 + 1.) * .140625 * (*eta * m1 + 1.) * (
	    mthrd * 9. * (1. - *zeta * *zeta) - *zeta * 2. * (*zeta * 9. * 
	    mthrd + 1.));
    der[der_dim1 * 14 + 1] = m1 * .140625 * (1. - *zeta * *zeta) * (*zeta * 
	    9. * mthrd + 1.) * (*eta * p1 + 1.);
    der[der_dim1 * 14 + 2] = p1 * .140625 * (1. - *zeta * *zeta) * (*zeta * 
	    9. * mthrd + 1.) * (*xi * m1 + 1.);
    der[der_dim1 * 14 + 3] = (*xi * m1 + 1.) * .140625 * (*eta * p1 + 1.) * (
	    mthrd * 9. * (1. - *zeta * *zeta) - *zeta * 2. * (*zeta * 9. * 
	    mthrd + 1.));
    der[der_dim1 * 15 + 1] = p1 * .140625 * (1. - *zeta * *zeta) * (*zeta * 
	    9. * mthrd + 1.) * (*eta * p1 + 1.);
    der[der_dim1 * 15 + 2] = p1 * .140625 * (1. - *zeta * *zeta) * (*zeta * 
	    9. * mthrd + 1.) * (*xi * p1 + 1.);
    der[der_dim1 * 15 + 3] = (*xi * p1 + 1.) * .140625 * (*eta * p1 + 1.) * (
	    mthrd * 9. * (1. - *zeta * *zeta) - *zeta * 2. * (*zeta * 9. * 
	    mthrd + 1.));
    der[(der_dim1 << 4) + 1] = p1 * .140625 * (1. - *zeta * *zeta) * (*zeta * 
	    9. * mthrd + 1.) * (*eta * m1 + 1.);
    der[(der_dim1 << 4) + 2] = m1 * .140625 * (1. - *zeta * *zeta) * (*zeta * 
	    9. * mthrd + 1.) * (*xi * p1 + 1.);
    der[(der_dim1 << 4) + 3] = (*xi * p1 + 1.) * .140625 * (*eta * m1 + 1.) * 
	    (mthrd * 9. * (1. - *zeta * *zeta) - *zeta * 2. * (*zeta * 9. * 
	    mthrd + 1.));
    der[der_dim1 * 17 + 1] = m1 * .140625 * (1. - *zeta * *zeta) * (*zeta * 
	    9. * pthrd + 1.) * (*eta * m1 + 1.);
    der[der_dim1 * 17 + 2] = m1 * .140625 * (1. - *zeta * *zeta) * (*zeta * 
	    9. * pthrd + 1.) * (*xi * m1 + 1.);
    der[der_dim1 * 17 + 3] = (*xi * m1 + 1.) * .140625 * (*eta * m1 + 1.) * (
	    pthrd * 9. * (1. - *zeta * *zeta) - *zeta * 2. * (*zeta * 9. * 
	    pthrd + 1.));
    der[der_dim1 * 18 + 1] = m1 * .140625 * (1. - *zeta * *zeta) * (*zeta * 
	    9. * pthrd + 1.) * (*eta * p1 + 1.);
    der[der_dim1 * 18 + 2] = p1 * .140625 * (1. - *zeta * *zeta) * (*zeta * 
	    9. * pthrd + 1.) * (*xi * m1 + 1.);
    der[der_dim1 * 18 + 3] = (*xi * m1 + 1.) * .140625 * (*eta * p1 + 1.) * (
	    pthrd * 9. * (1. - *zeta * *zeta) - *zeta * 2. * (*zeta * 9. * 
	    pthrd + 1.));
    der[der_dim1 * 19 + 1] = p1 * .140625 * (1. - *zeta * *zeta) * (*zeta * 
	    9. * pthrd + 1.) * (*eta * p1 + 1.);
    der[der_dim1 * 19 + 2] = p1 * .140625 * (1. - *zeta * *zeta) * (*zeta * 
	    9. * pthrd + 1.) * (*xi * p1 + 1.);
    der[der_dim1 * 19 + 3] = (*xi * p1 + 1.) * .140625 * (*eta * p1 + 1.) * (
	    pthrd * 9. * (1. - *zeta * *zeta) - *zeta * 2. * (*zeta * 9. * 
	    pthrd + 1.));
    der[der_dim1 * 20 + 1] = p1 * .140625 * (1. - *zeta * *zeta) * (*zeta * 
	    9. * pthrd + 1.) * (*eta * m1 + 1.);
    der[der_dim1 * 20 + 2] = m1 * .140625 * (1. - *zeta * *zeta) * (*zeta * 
	    9. * pthrd + 1.) * (*xi * p1 + 1.);
    der[der_dim1 * 20 + 3] = (*xi * p1 + 1.) * .140625 * (*eta * m1 + 1.) * (
	    pthrd * 9. * (1. - *zeta * *zeta) - *zeta * 2. * (*zeta * 9. * 
	    pthrd + 1.));
    der[der_dim1 * 21 + 1] = (*eta * m1 + 1.) * .015625 * (*zeta * p1 + 1.) * 
	    (m1 * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *
	    xi * 18. * (*xi * m1 + 1.));
    der[der_dim1 * 21 + 2] = (*xi * m1 + 1.) * .015625 * (*zeta * p1 + 1.) * (
	    m1 * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *
	    eta * 18. * (*eta * m1 + 1.));
    der[der_dim1 * 21 + 3] = (*xi * m1 + 1.) * .015625 * (*eta * m1 + 1.) * (
	    p1 * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *
	    zeta * 18. * (*zeta * p1 + 1.));
    der[der_dim1 * 22 + 1] = m1 * .140625 * (1. - *eta * *eta) * (*eta * 9. * 
	    mthrd + 1.) * (*zeta * p1 + 1.);
    der[der_dim1 * 22 + 2] = (*xi * m1 + 1.) * .140625 * (*zeta * p1 + 1.) * (
	    mthrd * 9. * (1. - *eta * *eta) - *eta * 2. * (*eta * 9. * mthrd 
	    + 1.));
    der[der_dim1 * 22 + 3] = p1 * .140625 * (1. - *eta * *eta) * (*eta * 9. * 
	    mthrd + 1.) * (*xi * m1 + 1.);
    der[der_dim1 * 23 + 1] = m1 * .140625 * (1. - *eta * *eta) * (*eta * 9. * 
	    pthrd + 1.) * (*zeta * p1 + 1.);
    der[der_dim1 * 23 + 2] = (*xi * m1 + 1.) * .140625 * (*zeta * p1 + 1.) * (
	    pthrd * 9. * (1. - *eta * *eta) - *eta * 2. * (*eta * 9. * pthrd 
	    + 1.));
    der[der_dim1 * 23 + 3] = p1 * .140625 * (1. - *eta * *eta) * (*eta * 9. * 
	    pthrd + 1.) * (*xi * m1 + 1.);
    der[der_dim1 * 24 + 1] = (*eta * p1 + 1.) * .015625 * (*zeta * p1 + 1.) * 
	    (m1 * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *
	    xi * 18. * (*xi * m1 + 1.));
    der[der_dim1 * 24 + 2] = (*xi * m1 + 1.) * .015625 * (*zeta * p1 + 1.) * (
	    p1 * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *
	    eta * 18. * (*eta * p1 + 1.));
    der[der_dim1 * 24 + 3] = (*xi * m1 + 1.) * .015625 * (*eta * p1 + 1.) * (
	    p1 * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *
	    zeta * 18. * (*zeta * p1 + 1.));
    der[der_dim1 * 25 + 1] = (*eta * p1 + 1.) * .140625 * (*zeta * p1 + 1.) * 
	    (mthrd * 9. * (1. - *xi * *xi) - *xi * 2. * (*xi * 9. * mthrd + 
	    1.));
    der[der_dim1 * 25 + 2] = p1 * .140625 * (1. - *xi * *xi) * (*xi * 9. * 
	    mthrd + 1.) * (*zeta * p1 + 1.);
    der[der_dim1 * 25 + 3] = p1 * .140625 * (1. - *xi * *xi) * (*xi * 9. * 
	    mthrd + 1.) * (*eta * p1 + 1.);
    der[der_dim1 * 26 + 1] = (*eta * p1 + 1.) * .140625 * (*zeta * p1 + 1.) * 
	    (pthrd * 9. * (1. - *xi * *xi) - *xi * 2. * (*xi * 9. * pthrd + 
	    1.));
    der[der_dim1 * 26 + 2] = p1 * .140625 * (1. - *xi * *xi) * (*xi * 9. * 
	    pthrd + 1.) * (*zeta * p1 + 1.);
    der[der_dim1 * 26 + 3] = p1 * .140625 * (1. - *xi * *xi) * (*xi * 9. * 
	    pthrd + 1.) * (*eta * p1 + 1.);
    der[der_dim1 * 27 + 1] = (*eta * p1 + 1.) * .015625 * (*zeta * p1 + 1.) * 
	    (p1 * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *
	    xi * 18. * (*xi * p1 + 1.));
    der[der_dim1 * 27 + 2] = (*xi * p1 + 1.) * .015625 * (*zeta * p1 + 1.) * (
	    p1 * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *
	    eta * 18. * (*eta * p1 + 1.));
    der[der_dim1 * 27 + 3] = (*xi * p1 + 1.) * .015625 * (*eta * p1 + 1.) * (
	    p1 * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *
	    zeta * 18. * (*zeta * p1 + 1.));
    der[der_dim1 * 28 + 1] = p1 * .140625 * (1. - *eta * *eta) * (*eta * 9. * 
	    pthrd + 1.) * (*zeta * p1 + 1.);
    der[der_dim1 * 28 + 2] = (*xi * p1 + 1.) * .140625 * (*zeta * p1 + 1.) * (
	    pthrd * 9. * (1. - *eta * *eta) - *eta * 2. * (*eta * 9. * pthrd 
	    + 1.));
    der[der_dim1 * 28 + 3] = p1 * .140625 * (1. - *eta * *eta) * (*eta * 9. * 
	    pthrd + 1.) * (*xi * p1 + 1.);
    der[der_dim1 * 29 + 1] = p1 * .140625 * (1. - *eta * *eta) * (*eta * 9. * 
	    mthrd + 1.) * (*zeta * p1 + 1.);
    der[der_dim1 * 29 + 2] = (*xi * p1 + 1.) * .140625 * (*zeta * p1 + 1.) * (
	    mthrd * 9. * (1. - *eta * *eta) - *eta * 2. * (*eta * 9. * mthrd 
	    + 1.));
    der[der_dim1 * 29 + 3] = p1 * .140625 * (1. - *eta * *eta) * (*eta * 9. * 
	    mthrd + 1.) * (*xi * p1 + 1.);
    der[der_dim1 * 30 + 1] = (*eta * m1 + 1.) * .015625 * (*zeta * p1 + 1.) * 
	    (p1 * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *
	    xi * 18. * (*xi * p1 + 1.));
    der[der_dim1 * 30 + 2] = (*xi * p1 + 1.) * .015625 * (*zeta * p1 + 1.) * (
	    m1 * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *
	    eta * 18. * (*eta * m1 + 1.));
    der[der_dim1 * 30 + 3] = (*xi * p1 + 1.) * .015625 * (*eta * m1 + 1.) * (
	    p1 * ((*xi * *xi + *eta * *eta + *zeta * *zeta) * 9. - 19.) + *
	    zeta * 18. * (*zeta * p1 + 1.));
    der[der_dim1 * 31 + 1] = (*eta * m1 + 1.) * .140625 * (*zeta * p1 + 1.) * 
	    (pthrd * 9. * (1. - *xi * *xi) - *xi * 2. * (*xi * 9. * pthrd + 
	    1.));
    der[der_dim1 * 31 + 2] = m1 * .140625 * (1. - *xi * *xi) * (*xi * 9. * 
	    pthrd + 1.) * (*zeta * p1 + 1.);
    der[der_dim1 * 31 + 3] = p1 * .140625 * (1. - *xi * *xi) * (*xi * 9. * 
	    pthrd + 1.) * (*eta * m1 + 1.);
    der[(der_dim1 << 5) + 1] = (*eta * m1 + 1.) * .140625 * (*zeta * p1 + 1.) 
	    * (mthrd * 9. * (1. - *xi * *xi) - *xi * 2. * (*xi * 9. * mthrd + 
	    1.));
    der[(der_dim1 << 5) + 2] = m1 * .140625 * (1. - *xi * *xi) * (*xi * 9. * 
	    mthrd + 1.) * (*zeta * p1 + 1.);
    der[(der_dim1 << 5) + 3] = p1 * .140625 * (1. - *xi * *xi) * (*xi * 9. * 
	    mthrd + 1.) * (*eta * m1 + 1.);

} /* brk32_ */

