/* quam12.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int quam12_(fun, ifun, der, ider, jder, xi, eta, itest)
doublereal *fun;
integer *ifun;
doublereal *der;
integer *ider, *jder;
doublereal *xi, *eta;
integer *itest;
{
    /* Initialized data */

    static char srname[6+1] = "QUAM12";

    /* System generated locals */
    integer der_dim1, der_offset;

    /* Local variables */
    extern doublereal veps_();
    static doublereal dummy;
    extern integer errmes_();
    static integer ierror;
    static doublereal val;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      QUAM12 returns the values of shape functions and their */
/*      derivatives at a specified point for an 12-noded c0 */
/*      continuous quadrilateral element.  The approximated */
/*      function will be continuous across element boundaries */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    21 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IFUN    length of vector FUN (.GE. 12) */
/*      IDER    first dimension of array DER (.GE. 2) */
/*      JDER    second dimension of array DER (.GE. 12) */
/*      XI      first local coordinate */
/*      ETA     second local coordinate */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      FUN     vector of length IFUN. FUN(I) contains the */
/*              value of the i'th shape function at (XI, ETA) */
/*              for i=1(1)12 */
/*      DER     array of dimension (IDER, JDER). DER(I, J) */
/*              contains the value of the derivative of the j'th */
/*              shape function with respect to the i'th */
/*              coordinate, for i=1(1)2 and j=1(1)12 */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE QUAM12(FUN,IFUN,DER,IDER,JDER,XI,ETA,ITEST) */
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
	if (*ifun < 12) {
	    ierror = 1;
	}
	if (*ider < 2 || *jder < 12) {
	    ierror = 2;
	}
	val = veps_(&dummy) + 1.;
	if (abs(*xi) > val || abs(*eta) > val) {
	    ierror = 3;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    fun[1] = (1. - *xi) * .03125 * (1. - *eta) * ((*xi * *xi + *eta * *eta) * 
	    9. - 10.);
    fun[2] = (1. - *xi) * .28125 * (1. - *eta * *eta) * (1. - *eta * 3.);
    fun[3] = (1. - *xi) * .28125 * (1. - *eta * *eta) * (*eta * 3. + 1.);
    fun[4] = (1. - *xi) * .03125 * (*eta + 1.) * ((*xi * *xi + *eta * *eta) * 
	    9. - 10.);
    fun[5] = (*eta + 1.) * .28125 * (1. - *xi * *xi) * (1. - *xi * 3.);
    fun[6] = (*eta + 1.) * .28125 * (1. - *xi * *xi) * (*xi * 3. + 1.);
    fun[7] = (*xi + 1.) * .03125 * (*eta + 1.) * ((*xi * *xi + *eta * *eta) * 
	    9. - 10.);
    fun[8] = (*xi + 1.) * .28125 * (1. - *eta * *eta) * (*eta * 3. + 1.);
    fun[9] = (*xi + 1.) * .28125 * (1. - *eta * *eta) * (1. - *eta * 3.);
    fun[10] = (*xi + 1.) * .03125 * (1. - *eta) * ((*xi * *xi + *eta * *eta) *
	     9. - 10.);
    fun[11] = (1. - *eta) * .28125 * (1. - *xi * *xi) * (*xi * 3. + 1.);
    fun[12] = (1. - *eta) * .28125 * (1. - *xi * *xi) * (1. - *xi * 3.);

    der[der_dim1 + 1] = (1. - *eta) * .03125 * (10. - *xi * 27. * *xi + *xi * 
	    18. - *eta * 9. * *eta);
    der[der_dim1 + 2] = (1. - *xi) * .03125 * (10. - *eta * 27. * *eta + *eta 
	    * 18. - *xi * 9. * *xi);
    der[(der_dim1 << 1) + 1] = (1. - *eta * *eta) * -.28125 * (1. - *eta * 3.)
	    ;
    der[(der_dim1 << 1) + 2] = (1. - *xi) * .28125 * (*eta * 9. * *eta - *eta 
	    * 2. - 3.);
    der[der_dim1 * 3 + 1] = (1. - *eta * *eta) * -.28125 * (*eta * 3. + 1.);
    der[der_dim1 * 3 + 2] = (1. - *xi) * .28125 * (*eta * -9. * *eta - *eta * 
	    2. + 3.);
    der[(der_dim1 << 2) + 1] = (*eta + 1.) * .03125 * (10. - *xi * 27. * *xi 
	    + *xi * 18. - *eta * 9. * *eta);
    der[(der_dim1 << 2) + 2] = (1. - *xi) * .03125 * (*eta * 27. * *eta - 10. 
	    + *eta * 18. + *xi * 9. * *xi);
    der[der_dim1 * 5 + 1] = (*eta + 1.) * .28125 * (*xi * 9. * *xi - *xi * 2. 
	    - 3.);
    der[der_dim1 * 5 + 2] = (1. - *xi * *xi) * .28125 * (1. - *xi * 3.);
    der[der_dim1 * 6 + 1] = (*eta + 1.) * .28125 * (3. - *xi * 2. - *xi * 9. *
	     *xi);
    der[der_dim1 * 6 + 2] = (1. - *xi * *xi) * .28125 * (*xi * 3. + 1.);
    der[der_dim1 * 7 + 1] = (*eta + 1.) * .03125 * (*xi * 27. * *xi - 10. + *
	    xi * 18. + *eta * 9. * *eta);
    der[der_dim1 * 7 + 2] = (*xi + 1.) * .03125 * (*eta * 27. * *eta - 10. + *
	    eta * 18. + *xi * 9. * *xi);
    der[(der_dim1 << 3) + 1] = (1. - *eta * *eta) * .28125 * (*eta * 3. + 1.);

    der[(der_dim1 << 3) + 2] = (*xi + 1.) * .28125 * (3. - *eta * 2. - *eta * 
	    9. * *eta);
    der[der_dim1 * 9 + 1] = (1. - *eta * *eta) * .28125 * (1. - *eta * 3.);
    der[der_dim1 * 9 + 2] = (*xi + 1.) * .28125 * (*eta * 9. * *eta - *eta * 
	    2. - 3.);
    der[der_dim1 * 10 + 1] = (1. - *eta) * .03125 * (*xi * 27. * *xi - 10. + *
	    xi * 18. + *eta * 9. * *eta);
    der[der_dim1 * 10 + 2] = (*xi + 1.) * .03125 * (10. - *eta * 27. * *eta + 
	    *eta * 18. - *xi * 9. * *xi);
    der[der_dim1 * 11 + 1] = (1. - *eta) * .28125 * (3. - *xi * 2. - *xi * 9. 
	    * *xi);
    der[der_dim1 * 11 + 2] = (1. - *xi * *xi) * -.28125 * (*xi * 3. + 1.);
    der[der_dim1 * 12 + 1] = (1. - *eta) * .28125 * (*xi * 9. * *xi - *xi * 
	    2. - 3.);
    der[der_dim1 * 12 + 2] = (1. - *xi * *xi) * -.28125 * (1. - *xi * 3.);

} /* quam12_ */

