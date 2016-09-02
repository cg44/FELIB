/* elgeom.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int elgeom_(nele, eltop, ieltop, jeltop, coord, icoord, 
	jcoord, geom, igeom, jgeom, dimen, itest)
integer *nele, *eltop, *ieltop, *jeltop;
doublereal *coord;
integer *icoord, *jcoord;
doublereal *geom;
integer *igeom, *jgeom, *dimen, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "ELGEOM";

    /* System generated locals */
    integer coord_dim1, coord_offset, eltop_dim1, eltop_offset, geom_dim1, 
	    geom_offset, i__1, i__2;

    /* Local variables */
    static integer idim, inod, jnod, nodel, jtest;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      ELGEOM constructs the element geometry array for the specified */
/*      element in the local node order */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    14 Feb 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      NELE    element number for which the nodal coordinate */
/*              array is required */
/*      ELTOP   contains element type, number of nodes, and */
/*              element topologies for all the elements */
/*      IELTOP  first dimension of ELTOP (.GE. maximum node */
/*              number on region) */
/*      JELTOP  second dimension of ELTOP (.GE. number of nodes */
/*              on element + 2) */
/*      COORD   COORD(I, J) contains the j'th global coordinate */
/*              of the i'th node */
/*      ICOORD  first dimension of COORD (.GE. maximum number of */
/*              nodes in the mesh) */
/*      JCOORD  second dimension of COORD (.GE. DIMEN) */
/*      IGEOM   first dimension of array GEOM (.GE. number of */
/*              nodes on the element) */
/*      JGEOM   second dimension of GEOM (.GE. DIMEN) */
/*      DIMEN   dimensionality of the mesh */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      GEOM    contains global coordinates of the nodes */
/*              associated with element NELE in the local node */
/*              numbering order */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE ELGEOM(NELE,ELTOP,IELTOP,JELTOP,COORD,ICOORD,JCOORD, */

/*     *                  GEOM,IGEOM,JGEOM,DIMEN,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    geom_dim1 = *igeom;
    geom_offset = geom_dim1 + 1;
    geom -= geom_offset;
    coord_dim1 = *icoord;
    coord_offset = coord_dim1 + 1;
    coord -= coord_offset;
    eltop_dim1 = *ieltop;
    eltop_offset = eltop_dim1 + 1;
    eltop -= eltop_offset;

    /* Function Body */

/*     Parameter checking */

    jtest = *itest;
    if (jtest != -1) {
	ierror = 0;
	if (*jgeom < *dimen) {
	    ierror = 4;
	}
	if (*jcoord < *dimen) {
	    ierror = 3;
	}
	if (*ieltop < *nele) {
	    ierror = 2;
	}
	if (*nele <= 0 || *dimen <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    nodel = eltop[*nele + (eltop_dim1 << 1)];

/*     Range checking on NODEL : NODEL+2 <= JELTOP and NODEL <= IGEOM */

    if (jtest != -1) {
	ierror = 0;
	if (*jeltop < nodel + 2) {
	    ierror = 5;
	}
	if (*igeom < nodel) {
	    ierror = 7;
	}
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

    i__1 = nodel;
    for (jnod = 1; jnod <= i__1; ++jnod) {
	inod = eltop[*nele + (jnod + 2) * eltop_dim1];

/*     Range checking on INOD : INOD <= ICOORD */

	if (jtest != -1) {
	    ierror = 0;
	    if (*icoord < inod) {
		ierror = 6;
	    }
	    *itest = errmes_(&jtest, &ierror, srname, 6L);
	    if (*itest != 0) {
		return 0;
	    }
	}

/*     Construct geometry array GEOM */

	i__2 = *dimen;
	for (idim = 1; idim <= i__2; ++idim) {
	    geom[jnod + idim * geom_dim1] = coord[inod + idim * coord_dim1];
/* L1000: */
	}
/* L1010: */
    }

} /* elgeom_ */

