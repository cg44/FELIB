/* errmes.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;


integer errmes_(itest, ierror, srname, srname_len)
integer *itest, *ierror;
char *srname;
ftnlen srname_len;
{
    /* Initialized data */

    static integer get = 0;
    static integer jtest = 1;

    /* Format strings */
    static char fmt_9980[] = "(\002 RELEASE 3.0  -  1 JAN 87\002)";
    static char fmt_9990[] = "(\002 ERROR DETECTED BY LEVEL 0 LIBRARY ROUTIN\
E \002,a,\002 - ITEST = \002,i5,//)";

    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    integer s_wsfe(), e_wsfe(), do_fio();
    /* Subroutine */ int s_stop();

    /* Local variables */
    static integer unit;
    extern /* Subroutine */ int adunit_(), erunit_();

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 0, 0, fmt_9980, 0 };
    static cilist io___5 = { 0, 0, 0, fmt_9990, 0 };
    static cilist io___6 = { 0, 0, 0, fmt_9990, 0 };


/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      ERRMES returns the value of IERROR or terminates the program, */
/*      printing a failure message */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    14 Feb 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      ITEST   contains either 0 (hard fail) or 1 (soft fail). */
/*              any other entry gives hard fail. */
/*      IERROR  contains the number of the detected error */
/*      SRNAME  contains up to 8 characters - usually a library */
/*              routine name */

/* ARGUMENTS out */
/*      ERRMES  routine name, contains the value of IERROR */

/* ROUTINES called */
/*      can call auxiliary routine in some versions of library */
/*      ERUNIT and ADUNIT */


/*     INTEGER FUNCTION ERRMES(ITEST,IERROR,SRNAME) */
/* ***********************************************************************
 */



/*     Hard failure */

    if (*itest == -99) {

/*     To return Release message */

	adunit_(&unit, &get, &jtest);
	io___4.ciunit = unit;
	s_wsfe(&io___4);
	e_wsfe();

    } else if (*itest == 1 || *ierror == 0) {

/*     Soft failure */

	ret_val = *ierror;
	if (*itest != 0 && *ierror != 0) {
	    adunit_(&unit, &get, &jtest);
	    io___5.ciunit = unit;
	    s_wsfe(&io___5);
	    do_fio(&c__1, srname, 6L);
	    do_fio(&c__1, (char *)&(*ierror), (ftnlen)sizeof(integer));
	    e_wsfe();
	}

    } else {

	erunit_(&unit, &get, &jtest);
	io___6.ciunit = unit;
	s_wsfe(&io___6);
	do_fio(&c__1, srname, 6L);
	do_fio(&c__1, (char *)&(*ierror), (ftnlen)sizeof(integer));
	e_wsfe();
	s_stop("", 0L);

    }
    return ret_val;
} /* errmes_ */

