/* update.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int update_(phi, iphi, rhs, irhs, totnod, dofnod, totdof, nf,
	 inf, jnf, itest)
doublereal *phi;
integer *iphi;
doublereal *rhs;
integer *irhs, *totnod, *dofnod, *totdof, *nf, *inf, *jnf, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "UPDATE";

    /* System generated locals */
    integer nf_dim1, nf_offset, i__1, i__2;

    /* Local variables */
    static integer i, j, k, l;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      UPDATE takes a full solution vector and a set of */
/*      updates and updates the solution vector */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 3.0  29 Jun 1986 (CJH, CG) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      PHI     solution vector */
/*      IPHI    dimension of vector PHI (IPHI .GE. TOTNOD*DOFNOD) */
/*      RHS     vector of updates */
/*      IRHS    dimension of vector RHS (IRHS .GE. TOTDOF) */
/*      TOTNOD  the number of nodes in the problem */
/*      DOFNOD  the maximum number of nodes per node */
/*      TOTDOF  the number of freedoms in RHS */
/*      NF      the nodal freedom array */
/*      INF     first dimension of NF (INF .GE. TOTNOD) */
/*      JNF     second dimension of NF (JNF .GE. DOFNOD) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      PHI     the updated solution vector */

/* ROUTINES called */
/*      ERRMES */

/*     SUBROUTINE UPDATE(PHI,IPHI,RHS,IRHS,TOTNOD,DOFNOD,TOTDOF,NF,INF, */

/*    *                  JNF,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    nf_dim1 = *inf;
    nf_offset = nf_dim1 + 1;
    nf -= nf_offset;
    --rhs;
    --phi;

    /* Function Body */

/*     Check array bounds */

    if (*itest != -1) {
	ierror = 0;
	if (*iphi < *totnod * *dofnod) {
	    ierror = 1;
	}
	if (*irhs < *totdof) {
	    ierror = 2;
	}
	if (*inf < *totnod) {
	    ierror = 3;
	}
	if (*jnf < *dofnod) {
	    ierror = 4;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

    i__1 = *totnod;
    for (i = 1; i <= i__1; ++i) {

	i__2 = *dofnod;
	for (j = 1; j <= i__2; ++j) {

	    k = nf[i + j * nf_dim1];

/*     Don't UPDATE restrained variables */

	    if (k != 0) {

		l = i + j - 1;

/*     UPDATE solution vector */

		phi[l] += rhs[k];
	    }
/* L1000: */
	}
/* L1010: */
    }


} /* update_ */

