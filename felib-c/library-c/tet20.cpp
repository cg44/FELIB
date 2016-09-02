/* tet20.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int tet20_(fun, ifun, der, ider, jder, xi, eta, zeta, itest)
doublereal *fun;
integer *ifun;
doublereal *der;
integer *ider, *jder;
doublereal *xi, *eta, *zeta;
integer *itest;
{
    /* Initialized data */

    static char srname[6+1] = "TET20 ";

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
/*      TET20 returns the values of shape functions and derivatives */
/*      at a specified point for an 20-noded tetrahedral element. */
/*      The function is continuous across element boundaries. */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    10 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IFUN    dimension of vector FUN (.GE. 20) */
/*      IDER    first dimension of array DER (.GE. 3) */
/*      JDER    second dimension of array DER (.GE. 20) */
/*      XI      value of local coordinate at which function and */
/*              derivative values required */
/*      ETA     value of second local coordinate */
/*      ZETA    value of third local coordinate */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      FUN     real vector of dimension IFUN.  FUN(I) contains */
/*              value of i'th shape function at (XI, ETA, ZETA) */
/*      DER     real array of dimensions (IDER, JDER).  DER(I ,J) */
/*              contains the derivative of the j'th shape */
/*              function with respect to the i'th coordinate at */
/*              the point (XI, ETA, ZETA) */

/* ROUTINES called */
/*      ERRMES */

/*     SUBROUTINE TET20(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST) */
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
	if (*ifun < 20) {
	    ierror = 1;
	}
	if (*ider < 3 || *jder < 20) {
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

    fun[1] = (l1 * 3. - 1.) * .5 * (l1 * 3. - 2.) * l1;
    fun[2] = l1 * 4.5 * l2 * (l1 * 3. - 1.);
    fun[3] = l2 * 4.5 * l1 * (l2 * 3. - 1.);
    fun[4] = (l2 * 3. - 1.) * .5 * (l2 * 3. - 2.) * l2;
    fun[5] = l2 * 4.5 * l3 * (l2 * 3. - 1.);
    fun[6] = l3 * 4.5 * l2 * (l3 * 3. - 1.);
    fun[7] = (l3 * 3. - 1.) * .5 * (l3 * 3. - 2.) * l3;
    fun[8] = l3 * 4.5 * l1 * (l3 * 3. - 1.);
    fun[9] = l1 * 4.5 * l3 * (l1 * 3. - 1.);
    fun[10] = l1 * 27. * l2 * l3;
    fun[11] = l1 * 4.5 * l4 * (l1 * 3. - 1.);
    fun[12] = l4 * 27. * l1 * l2;
    fun[13] = l2 * 4.5 * l4 * (l2 * 3. - 1.);
    fun[14] = l2 * 27. * l3 * l4;
    fun[15] = l3 * 4.5 * l4 * (l3 * 3. - 1.);
    fun[16] = l3 * 27. * l4 * l1;
    fun[17] = l4 * 4.5 * l1 * (l4 * 3. - 1.);
    fun[18] = l4 * 4.5 * l2 * (l4 * 3. - 1.);
    fun[19] = l4 * 4.5 * l3 * (l4 * 3. - 1.);
    fun[20] = (l4 * 3. - 1.) * .5 * (l4 * 3. - 2.) * l4;

/*     Derivatives */

    der[der_dim1 + 1] = (l1 * 27. * l1 - l1 * 18. + 2.) * .5 * dl1x;
    der[der_dim1 + 2] = (l1 * 27. * l1 - l1 * 18. + 2.) * .5 * dl1y;
    der[der_dim1 + 3] = (l1 * 27. * l1 - l1 * 18. + 2.) * .5 * dl1z;
    der[(der_dim1 << 1) + 1] = l2 * 4.5 * (l1 * 6. - 1.) * dl1x + l1 * 4.5 * (
	    l1 * 3. - 1.) * dl2x;
    der[(der_dim1 << 1) + 2] = l2 * 4.5 * (l1 * 6. - 1.) * dl1y + l1 * 4.5 * (
	    l1 * 3. - 1.) * dl2y;
    der[(der_dim1 << 1) + 3] = l2 * 4.5 * (l1 * 6. - 1.) * dl1z + l1 * 4.5 * (
	    l1 * 3. - 1.) * dl2z;
    der[der_dim1 * 3 + 1] = l2 * 4.5 * (l2 * 3. - 1.) * dl1x + l1 * 4.5 * (l2 
	    * 6. - 1.) * dl2x;
    der[der_dim1 * 3 + 2] = l2 * 4.5 * (l2 * 3. - 1.) * dl1y + l1 * 4.5 * (l2 
	    * 6. - 1.) * dl2y;
    der[der_dim1 * 3 + 3] = l2 * 4.5 * (l2 * 3. - 1.) * dl1z + l1 * 4.5 * (l2 
	    * 6. - 1.) * dl2z;
    der[(der_dim1 << 2) + 1] = (l2 * 27. * l2 - l2 * 18. + 2.) * .5 * dl2x;
    der[(der_dim1 << 2) + 2] = (l2 * 27. * l2 - l2 * 18. + 2.) * .5 * dl2y;
    der[(der_dim1 << 2) + 3] = (l2 * 27. * l2 - l2 * 18. + 2.) * .5 * dl2z;
    der[der_dim1 * 5 + 1] = l3 * 4.5 * (l2 * 6. - 1.) * dl2x + l2 * 4.5 * (l2 
	    * 3. - 1.) * dl3x;
    der[der_dim1 * 5 + 2] = l3 * 4.5 * (l2 * 6. - 1.) * dl2y + l2 * 4.5 * (l2 
	    * 3. - 1.) * dl3y;
    der[der_dim1 * 5 + 3] = l3 * 4.5 * (l2 * 6. - 1.) * dl2z + l2 * 4.5 * (l2 
	    * 3. - 1.) * dl3z;
    der[der_dim1 * 6 + 1] = l3 * 4.5 * (l3 * 3. - 1.) * dl2x + l2 * 4.5 * (l3 
	    * 6. - 1.) * dl3x;
    der[der_dim1 * 6 + 2] = l3 * 4.5 * (l3 * 3. - 1.) * dl2y + l2 * 4.5 * (l3 
	    * 6. - 1.) * dl3y;
    der[der_dim1 * 6 + 3] = l3 * 4.5 * (l3 * 3. - 1.) * dl2z + l2 * 4.5 * (l3 
	    * 6. - 1.) * dl3z;
    der[der_dim1 * 7 + 1] = (l3 * 27. * l3 - l3 * 18. + 2.) * .5 * dl3x;
    der[der_dim1 * 7 + 2] = (l3 * 27. * l3 - l3 * 18. + 2.) * .5 * dl3y;
    der[der_dim1 * 7 + 3] = (l3 * 27. * l3 - l3 * 18. + 2.) * .5 * dl3z;
    der[(der_dim1 << 3) + 1] = l1 * 4.5 * (l3 * 6. - 1.) * dl3x + l3 * 4.5 * (
	    l3 * 3. - 1.) * dl1x;
    der[(der_dim1 << 3) + 2] = l1 * 4.5 * (l3 * 6. - 1.) * dl3y + l3 * 4.5 * (
	    l3 * 3. - 1.) * dl1y;
    der[(der_dim1 << 3) + 3] = l1 * 4.5 * (l3 * 6. - 1.) * dl3z + l3 * 4.5 * (
	    l3 * 3. - 1.) * dl1z;
    der[der_dim1 * 9 + 1] = l1 * 4.5 * (l1 * 3. - 1.) * dl3x + l3 * 4.5 * (l1 
	    * 6. - 1.) * dl1x;
    der[der_dim1 * 9 + 2] = l1 * 4.5 * (l1 * 3. - 1.) * dl3y + l3 * 4.5 * (l1 
	    * 6. - 1.) * dl1y;
    der[der_dim1 * 9 + 3] = l1 * 4.5 * (l1 * 3. - 1.) * dl3z + l3 * 4.5 * (l1 
	    * 6. - 1.) * dl1z;
    der[der_dim1 * 10 + 1] = l2 * 27. * l3 * dl1x + l3 * 27. * l1 * dl2x + l1 
	    * 27. * l2 * dl3x;
    der[der_dim1 * 10 + 2] = l2 * 27. * l3 * dl1y + l3 * 27. * l1 * dl2y + l1 
	    * 27. * l2 * dl3y;
    der[der_dim1 * 10 + 3] = l2 * 27. * l3 * dl1z + l3 * 27. * l1 * dl2z + l1 
	    * 27. * l2 * dl3z;
    der[der_dim1 * 11 + 1] = l4 * 4.5 * (l1 * 6. - 1.) * dl1x + l1 * 4.5 * (
	    l1 * 3. - 1.) * dl4x;
    der[der_dim1 * 11 + 2] = l4 * 4.5 * (l1 * 6. - 1.) * dl1y + l1 * 4.5 * (
	    l1 * 3. - 1.) * dl4y;
    der[der_dim1 * 11 + 3] = l4 * 4.5 * (l1 * 6. - 1.) * dl1z + l1 * 4.5 * (
	    l1 * 3. - 1.) * dl4z;
    der[der_dim1 * 12 + 1] = l2 * 27. * l4 * dl1x + l1 * 27. * l4 * dl2x + l1 
	    * 27. * l2 * dl4x;
    der[der_dim1 * 12 + 2] = l2 * 27. * l4 * dl1y + l1 * 27. * l4 * dl2y + l1 
	    * 27. * l2 * dl4y;
    der[der_dim1 * 12 + 3] = l2 * 27. * l4 * dl1z + l1 * 27. * l4 * dl2z + l1 
	    * 27. * l2 * dl4z;
    der[der_dim1 * 13 + 1] = l4 * 4.5 * (l2 * 6. - 1.) * dl2x + l2 * 4.5 * (
	    l2 * 3. - 1.) * dl4x;
    der[der_dim1 * 13 + 2] = l4 * 4.5 * (l2 * 6. - 1.) * dl2y + l2 * 4.5 * (
	    l2 * 3. - 1.) * dl4y;
    der[der_dim1 * 13 + 3] = l4 * 4.5 * (l2 * 6. - 1.) * dl2z + l2 * 4.5 * (
	    l2 * 3. - 1.) * dl4z;
    der[der_dim1 * 14 + 1] = l3 * 27. * l4 * dl2x + l2 * 27. * l4 * dl3x + l2 
	    * 27. * l3 * dl4x;
    der[der_dim1 * 14 + 2] = l3 * 27. * l4 * dl2y + l2 * 27. * l4 * dl3y + l2 
	    * 27. * l3 * dl4y;
    der[der_dim1 * 14 + 3] = l3 * 27. * l4 * dl2z + l2 * 27. * l4 * dl3z + l2 
	    * 27. * l3 * dl4z;
    der[der_dim1 * 15 + 1] = l4 * 4.5 * (l3 * 6. - 1.) * dl3x + l3 * 4.5 * (
	    l3 * 3. - 1.) * dl4x;
    der[der_dim1 * 15 + 2] = l4 * 4.5 * (l3 * 6. - 1.) * dl3y + l3 * 4.5 * (
	    l3 * 3. - 1.) * dl4y;
    der[der_dim1 * 15 + 3] = l4 * 4.5 * (l3 * 6. - 1.) * dl3z + l3 * 4.5 * (
	    l3 * 3. - 1.) * dl4z;
    der[(der_dim1 << 4) + 1] = l3 * 27. * l4 * dl1x + l1 * 27. * l4 * dl3x + 
	    l1 * 27. * l3 * dl4x;
    der[(der_dim1 << 4) + 2] = l3 * 27. * l4 * dl1y + l1 * 27. * l4 * dl3y + 
	    l1 * 27. * l3 * dl4y;
    der[(der_dim1 << 4) + 3] = l3 * 27. * l4 * dl1z + l1 * 27. * l4 * dl3z + 
	    l1 * 27. * l3 * dl4z;
    der[der_dim1 * 17 + 1] = l4 * 4.5 * (l4 * 3. - 1.) * dl1x + l1 * 4.5 * (
	    l4 * 6. - 1.) * dl4x;
    der[der_dim1 * 17 + 2] = l4 * 4.5 * (l4 * 3. - 1.) * dl1y + l1 * 4.5 * (
	    l4 * 6. - 1.) * dl4y;
    der[der_dim1 * 17 + 3] = l4 * 4.5 * (l4 * 3. - 1.) * dl1z + l1 * 4.5 * (
	    l4 * 6. - 1.) * dl4z;
    der[der_dim1 * 18 + 1] = l4 * 4.5 * (l4 * 3. - 1.) * dl2x + l2 * 4.5 * (
	    l4 * 6. - 1.) * dl4x;
    der[der_dim1 * 18 + 2] = l4 * 4.5 * (l4 * 3. - 1.) * dl2y + l2 * 4.5 * (
	    l4 * 6. - 1.) * dl4y;
    der[der_dim1 * 18 + 3] = l4 * 4.5 * (l4 * 3. - 1.) * dl2z + l2 * 4.5 * (
	    l4 * 6. - 1.) * dl4z;
    der[der_dim1 * 19 + 1] = l4 * 4.5 * (l4 * 3. - 1.) * dl3x + l3 * 4.5 * (
	    l4 * 6. - 1.) * dl4x;
    der[der_dim1 * 19 + 2] = l4 * 4.5 * (l4 * 3. - 1.) * dl3y + l3 * 4.5 * (
	    l4 * 6. - 1.) * dl4y;
    der[der_dim1 * 19 + 3] = l4 * 4.5 * (l4 * 3. - 1.) * dl3z + l3 * 4.5 * (
	    l4 * 6. - 1.) * dl4z;
    der[der_dim1 * 20 + 1] = (l4 * 27. * l4 - l4 * 18. + 2.) * .5 * dl4x;
    der[der_dim1 * 20 + 2] = (l4 * 27. * l4 - l4 * 18. + 2.) * .5 * dl4y;
    der[der_dim1 * 20 + 3] = (l4 * 27. * l4 - l4 * 18. + 2.) * .5 * dl4z;

} /* tet20_ */

