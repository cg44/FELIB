/* tet4.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int tet4_(fun, ifun, der, ider, jder, xi, eta, zeta, itest)
doublereal *fun;
integer *ifun;
doublereal *der;
integer *ider, *jder;
doublereal *xi, *eta, *zeta;
integer *itest;
{
    /* Initialized data */

    static char srname[6+1] = "TET4  ";

    /* System generated locals */
    integer der_dim1, der_offset;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      TET4 returns the values of shape functions and derivatives */
/*      at a specified point for an 10-noded tetrahedral element. */
/*      The function is continuous across element boundaries. */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    10 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IFUN    dimension of vector FUN (.GE. 4) */
/*      IDER    first dimension of array DER (.GE. 3) */
/*      JDER    second dimension of array DER (.GE. 4) */
/*      XI      value of local coordinate at which function and */
/*              derivative values required */
/*      ETA     value of second local coordinate */
/*      ZETA    value of third local coordinate */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      FUN     real vector of dimension IFUN.  FUN(I) contains */
/*              value of i'th shape function at (XI, ETA, ZETA) */
/*      DER     real array of dimensions (IDER, JDER).  DER(I, J) */
/*              contains the derivative of the j'th shape */
/*              function with respect to the i'th coordinate at */
/*              the point (XI, ETA, ZETA) */

/* ROUTINES called */
/*      ERRMES */

/*     SUBROUTINE TET4(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST) */
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
	if (*ifun < 4) {
	    ierror = 1;
	}
	if (*ider < 3 || *jder < 4) {
	    ierror = 2;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    fun[1] = (*xi * 8. + 3. - sqrt(2.) * 2. * *zeta) * .083333333333333329;
    fun[2] = (3. - (*xi + sqrt(3.) * *eta) * 4. - sqrt(2.) * 2. * *zeta) * 
	    .083333333333333329;
    fun[3] = (3. - (*xi - sqrt(3.) * *eta) * 4. - sqrt(2.) * 2. * *zeta) * 
	    .083333333333333329;
    fun[4] = (sqrt(2.) * 2. * *zeta + 1.) * .25;

    der[der_dim1 + 1] = .66666666666666663;
    der[der_dim1 + 2] = 0.;
    der[der_dim1 + 3] = -sqrt(2.) / 6.;
    der[(der_dim1 << 1) + 1] = -.33333333333333331;
    der[(der_dim1 << 1) + 2] = -sqrt(3.) / 3.;
    der[(der_dim1 << 1) + 3] = -sqrt(2.) / 6.;
    der[der_dim1 * 3 + 1] = -.33333333333333331;
    der[der_dim1 * 3 + 2] = sqrt(3.) / 3.;
    der[der_dim1 * 3 + 3] = -sqrt(2.) / 6.;
    der[(der_dim1 << 2) + 1] = 0.;
    der[(der_dim1 << 2) + 2] = 0.;
    der[(der_dim1 << 2) + 3] = sqrt(2.) / 2.;

} /* tet4_ */

