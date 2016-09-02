/* tet10.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int tet10_(fun, ifun, der, ider, jder, xi, eta, zeta, itest)
doublereal *fun;
integer *ifun;
doublereal *der;
integer *ider, *jder;
doublereal *xi, *eta, *zeta;
integer *itest;
{
    /* Initialized data */

    static char srname[6+1] = "TET10 ";

    /* System generated locals */
    integer der_dim1, der_offset;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal l1, l2, l3, l4;
    extern integer errmes_();
    static integer ierror;
    static doublereal dl1x, dl1y, dl1z, dl2x, dl2y, dl2z, dl3x, dl3y, dl3z, 
	    dl4x, dl4y, dl4z;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      TET10 returns the values of shape functions and derivatives */
/*      at a specified point for an 10-noded tetrahedral element. */
/*      The function is continuous across element boundaries. */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    10 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IFUN    dimension of vector FUN (.GE. 10) */
/*      IDER    first dimension of array DER (.GE. 3) */
/*      JDER    second dimension of array DER (.GE. 10) */
/*      XI      value of local coordinate at which function and */
/*              derivative values required */
/*      ETA     value of second local coordinate */
/*      ZETA    value of third local coordinate */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      FUN     real vector of dimension IFUN. FUN(I) contains */
/*              value of i'th shape function at (XI, ETA, ZETA) */
/*      DER     real array of dimensions (IDER, JDER). DER(I,J) */
/*              contains the derivative of the j'th shape */
/*              function with respect to the i'th coordinate at */
/*              the point (XI, ETA, ZETA) */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE TET10(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST) */
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
	if (*ifun < 10) {
	    ierror = 1;
	}
	if (*ider < 3 || *jder < 10) {
	    ierror = 2;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    l1 = (*xi * 8. + 3. - sqrt(2.) * 2. * *zeta) * .083333333333333329;
    l2 = (3. - (*xi + sqrt(3.) * *eta) * 4. - sqrt(2.) * 2. * *zeta) * 
	    .083333333333333329;
    l3 = (3. - (*xi - sqrt(3.) * *eta) * 4. - sqrt(2.) * 2. * *zeta) * 
	    .083333333333333329;
    l4 = (sqrt(2.) * 2. * *zeta + 1.) * .25;

    dl1x = .66666666666666663;
    dl1y = 0.;
    dl1z = -1. / (sqrt(2.) * 3.);
    dl2x = -.33333333333333331;
    dl2y = -1. / sqrt(3.);
    dl2z = -1. / (sqrt(2.) * 3.);
    dl3x = -.33333333333333331;
    dl3y = 1. / sqrt(3.);
    dl3z = -1. / (sqrt(2.) * 3.);
    dl4x = 0.;
    dl4y = 0.;
    dl4z = 1. / sqrt(2.);

/*     Shape functions */

    fun[1] = (l1 * 2. - 1.) * l1;
    fun[2] = l1 * 4. * l2;
    fun[3] = (l2 * 2. - 1.) * l2;
    fun[4] = l2 * 4. * l3;
    fun[5] = (l3 * 2. - 1.) * l3;
    fun[6] = l3 * 4. * l1;
    fun[7] = l1 * 4. * l4;
    fun[8] = l2 * 4. * l4;
    fun[9] = l3 * 4. * l4;
    fun[10] = (l4 * 2. - 1.) * l4;

/*     Derivatives */

    der[der_dim1 + 1] = (l1 * 4. - 1.) * dl1x;
    der[der_dim1 + 2] = (l1 * 4. - 1.) * dl1y;
    der[der_dim1 + 3] = (l1 * 4. - 1.) * dl1z;
    der[(der_dim1 << 1) + 1] = l2 * 4. * dl1x + l1 * 4. * dl2x;
    der[(der_dim1 << 1) + 2] = l2 * 4. * dl1y + l1 * 4. * dl2y;
    der[(der_dim1 << 1) + 3] = l2 * 4. * dl1z + l1 * 4. * dl2z;
    der[der_dim1 * 3 + 1] = (l2 * 4. - 1.) * dl2x;
    der[der_dim1 * 3 + 2] = (l2 * 4. - 1.) * dl2y;
    der[der_dim1 * 3 + 3] = (l2 * 4. - 1.) * dl2z;
    der[(der_dim1 << 2) + 1] = l3 * 4. * dl2x + l2 * 4. * dl3x;
    der[(der_dim1 << 2) + 2] = l3 * 4. * dl2y + l2 * 4. * dl3y;
    der[(der_dim1 << 2) + 3] = l3 * 4. * dl2z + l2 * 4. * dl3z;
    der[der_dim1 * 5 + 1] = (l3 * 4. - 1.) * dl3x;
    der[der_dim1 * 5 + 2] = (l3 * 4. - 1.) * dl3y;
    der[der_dim1 * 5 + 3] = (l3 * 4. - 1.) * dl3z;
    der[der_dim1 * 6 + 1] = l1 * 4. * dl3x + l3 * 4. * dl1x;
    der[der_dim1 * 6 + 2] = l1 * 4. * dl3y + l3 * 4. * dl1y;
    der[der_dim1 * 6 + 3] = l1 * 4. * dl3z + l3 * 4. * dl1z;
    der[der_dim1 * 7 + 1] = l4 * 4. * dl1x + l1 * 4. * dl4x;
    der[der_dim1 * 7 + 2] = l4 * 4. * dl1y + l1 * 4. * dl4y;
    der[der_dim1 * 7 + 3] = l4 * 4. * dl1z + l1 * 4. * dl4z;
    der[(der_dim1 << 3) + 1] = l4 * 4. * dl2x + l2 * 4. * dl4x;
    der[(der_dim1 << 3) + 2] = l4 * 4. * dl2y + l2 * 4. * dl4y;
    der[(der_dim1 << 3) + 3] = l4 * 4. * dl2z + l2 * 4. * dl4z;
    der[der_dim1 * 9 + 1] = l4 * 4. * dl3x + l3 * 4. * dl4x;
    der[der_dim1 * 9 + 2] = l4 * 4. * dl3y + l3 * 4. * dl4y;
    der[der_dim1 * 9 + 3] = l4 * 4. * dl3z + l3 * 4. * dl4z;
    der[der_dim1 * 10 + 1] = (l4 * 4. - 1.) * dl4x;
    der[der_dim1 * 10 + 2] = (l4 * 4. - 1.) * dl4y;
    der[der_dim1 * 10 + 3] = (l4 * 4. - 1.) * dl4z;

} /* tet10_ */

