/*  -- translated by f2c (version of 26 February 1990  17:38:00).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* *********************************************************************** */
/* *********************************************************************** */

/*    COPYRIGHT (C) 1987 : SERC, RUTHERFORD APPLETON LABORATORY */
/*                         CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX */

/* Main program */ MAIN__()
{
    /* Initialized data */

    static integer iabss = 3;
    static integer ijacin = 3;
    static integer ilder = 3;
    static integer ip = 3;
    static integer ipd = 3;
    static integer iscvec = 8;
    static integer isteer = 8;
    static integer iwght = 9;
    static integer jabss = 9;
    static integer jcoord = 3;
    static integer jdtpd = 8;
    static integer idtpd = 8;
    static integer jelk = 8;
    static integer jgder = 8;
    static integer jgdert = 3;
    static integer jgeom = 3;
    static integer jjac = 3;
    static integer jjacin = 3;
    static integer jlder = 8;
    static integer jnf = 1;
    static integer jp = 3;
    static integer jpd = 8;
    static integer ielk = 8;
    static doublereal scale = 1e10;
    static integer icoord = 100;
    static integer ieltop = 100;
    static integer inf = 100;
    static integer irhs = 100;
    static integer isysk = 100;
    static integer jeltop = 10;
    static integer jsysk = 25;
    static integer nin = 5;
    static integer nout = 6;
    static integer ielq = 8;
    static integer ifun = 8;
    static integer igder = 3;
    static integer igdert = 8;
    static integer igeom = 8;
    static integer ijac = 3;

    /* Format strings */
    static char fmt_9010[] = "(//\002 **** NODAL GEOMETRY ****\002//\002 \
\002)";
    static char fmt_8010[] = "(16i5)";
    static char fmt_9020[] = "(\002 \002,16i5)";
    static char fmt_8020[] = "(i5,6f10.0)";
    static char fmt_9030[] = "(\002 \002,i5,6f10.5)";
    static char fmt_9040[] = "(//\002 **** ELEMENT TOPOLOGY ****\002//\002\
 \002)";
    static char fmt_9050[] = "(//\002 **** PERMEABILITIES ****\002//\002 \
\002)";
    static char fmt_8030[] = "(2f10.0)";
    static char fmt_9060[] = "(\002 \002,2f10.5)";
    static char fmt_9070[] = "(//\002 **** SOURCE STRENGTH ****\002//\002\
 \002)";
    static char fmt_9080[] = "(//\002 **** BOUNDARY CONDITIONS ****\002//\
\002 \002)";
    static char fmt_8040[] = "(i5,6f10.0)";
    static char fmt_9090[] = "(//\002 **** NODAL POTENTIALS ****\002//\002\
 \002)";

    /* System generated locals */
    integer i_1, i_2, i_3, i_4;

    /* Builtin functions */
    integer s_wsfe(), e_wsfe(), s_rsfe(), do_fio(), e_rsfe();
    /* Subroutine */ int s_stop();

    /* Local variables */
    static doublereal gder[24]	/* was [3][8] */;
    static integer nele;
    static doublereal bval[30], geom[24]	/* was [8][3] */, abss[27]	
	    /* was [3][9] */, lder[24]	/* was [3][8] */, dtpd[64]	/* 
	    was [8][8] */, wght[9], quot, sysk[2500]	/* was [100][25] */;
    extern /* Subroutine */ int quam4_(), qqua4_();
    static integer i, j, hband;
    static doublereal p[9]	/* was [3][3] */, jacin[9]	/* was [3][3] 
	    */;
    static integer bnode[30], dofel, dimen;
    static doublereal x, y;
    static integer iquad, nodel;
    static doublereal coord[300]	/* was [100][3] */, gdert[24]	/* 
	    was [8][3] */, scvec[8];
    static integer elnum;
    extern /* Subroutine */ int asrhs_();
    static integer eltop[1000]	/* was [100][10] */, steer[8];
    static logical first;
    static integer itest;
    extern /* Subroutine */ int assym_();
    static integer eltyp, nf[100]	/* was [100][1] */;
    static doublereal pd[24]	/* was [3][8] */;
    extern /* Subroutine */ int vecadd_(), matadd_(), fredif_();
    static doublereal xi;
    static integer bndnod;
    extern /* Subroutine */ int elgeom_();
    static integer dofnod;
    extern /* Subroutine */ int direct_(), scaprd_(), matran_(), chosol_(), 
	    vecnul_();
    static integer totdof, nodnum;
    extern /* Subroutine */ int matnul_(), matmul_(), matinv_();
    static integer totnod;
    extern /* Subroutine */ int prtval_();
    static integer totels;
    static doublereal strgth, jac[9]	/* was [3][3] */;
    static integer dif;
    static doublereal eta, elk[64]	/* was [8][8] */, det, elq[8];
    extern doublereal src_();
    static doublereal fun[8], rhs[100];
    static integer nqp;

    /* Fortran I/O blocks */
    static cilist io__41 = { 0, 0, 0, fmt_9010, 0 };
    static cilist io__42 = { 0, 0, 0, fmt_8010, 0 };
    static cilist io__45 = { 0, 0, 0, fmt_9020, 0 };
    static cilist io__47 = { 0, 0, 0, fmt_8020, 0 };
    static cilist io__51 = { 0, 0, 0, fmt_9030, 0 };
    static cilist io__52 = { 0, 0, 0, fmt_9040, 0 };
    static cilist io__53 = { 0, 0, 0, fmt_8010, 0 };
    static cilist io__57 = { 0, 0, 0, fmt_9020, 0 };
    static cilist io__58 = { 0, 0, 0, fmt_8010, 0 };
    static cilist io__61 = { 0, 0, 0, fmt_9020, 0 };
    static cilist io__62 = { 0, 0, 0, fmt_9050, 0 };
    static cilist io__64 = { 0, 0, 0, fmt_8030, 0 };
    static cilist io__65 = { 0, 0, 0, fmt_9060, 0 };
    static cilist io__66 = { 0, 0, 0, fmt_9070, 0 };
    static cilist io__67 = { 0, 0, 0, fmt_8030, 0 };
    static cilist io__69 = { 0, 0, 0, fmt_9060, 0 };
    static cilist io__70 = { 0, 0, 0, fmt_9080, 0 };
    static cilist io__71 = { 0, 0, 0, fmt_8010, 0 };
    static cilist io__73 = { 0, 0, 0, fmt_9020, 0 };
    static cilist io__74 = { 0, 0, 0, fmt_8010, 0 };
    static cilist io__76 = { 0, 0, 0, fmt_9020, 0 };
    static cilist io__77 = { 0, 0, 0, fmt_8040, 0 };
    static cilist io__80 = { 0, 0, 0, fmt_9030, 0 };
    static cilist io__113 = { 0, 0, 0, fmt_9090, 0 };



/*                            PROBLEM SIZE DEPENDENT ARRAYS */



/*                            PROBLEM SIZE DEPENDENT DATA STATEMENTS */



/*                            SET ITEST FOR FULL CHECKING */

    itest = 0;

/* *                           ********************** */
/* *                           *                    * */
/* *                           * INPUT DATA SECTION * */
/* *                           *                    * */
/* *                           ********************** */

/*                            INPUT OF NODAL GEOMETRY */

    io__41.ciunit = nout;
    s_wsfe(&io__41);
    e_wsfe();
    io__42.ciunit = nin;
    s_rsfe(&io__42);
    do_fio(&c__1, (char *)&totnod, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&dimen, (ftnlen)sizeof(integer));
    e_rsfe();
    io__45.ciunit = nout;
    s_wsfe(&io__45);
    do_fio(&c__1, (char *)&totnod, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&dimen, (ftnlen)sizeof(integer));
    e_wsfe();
    i_1 = totnod;
    for (i = 1; i <= i_1; ++i) {
	io__47.ciunit = nin;
	s_rsfe(&io__47);
	do_fio(&c__1, (char *)&nodnum, (ftnlen)sizeof(integer));
	i_2 = dimen;
	for (j = 1; j <= i_2; ++j) {
	    do_fio(&c__1, (char *)&coord[nodnum + j * 100 - 101], (ftnlen)
		    sizeof(doublereal));
	}
	e_rsfe();
	io__51.ciunit = nout;
	s_wsfe(&io__51);
	do_fio(&c__1, (char *)&nodnum, (ftnlen)sizeof(integer));
	i_2 = dimen;
	for (j = 1; j <= i_2; ++j) {
	    do_fio(&c__1, (char *)&coord[nodnum + j * 100 - 101], (ftnlen)
		    sizeof(doublereal));
	}
	e_wsfe();
/* L1010: */
    }

/*                            INPUT OF ELEMENT TOPOLOGY */

    io__52.ciunit = nout;
    s_wsfe(&io__52);
    e_wsfe();
    io__53.ciunit = nin;
    s_rsfe(&io__53);
    do_fio(&c__1, (char *)&eltyp, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&totels, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nodel, (ftnlen)sizeof(integer));
    e_rsfe();
    io__57.ciunit = nout;
    s_wsfe(&io__57);
    do_fio(&c__1, (char *)&eltyp, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&totels, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nodel, (ftnlen)sizeof(integer));
    e_wsfe();
    i_1 = totels;
    for (i = 1; i <= i_1; ++i) {
	io__58.ciunit = nin;
	s_rsfe(&io__58);
	do_fio(&c__1, (char *)&elnum, (ftnlen)sizeof(integer));
	i_2 = nodel;
	for (j = 1; j <= i_2; ++j) {
	    do_fio(&c__1, (char *)&eltop[elnum + (j + 2) * 100 - 101], (
		    ftnlen)sizeof(integer));
	}
	e_rsfe();
	io__61.ciunit = nout;
	s_wsfe(&io__61);
	do_fio(&c__1, (char *)&elnum, (ftnlen)sizeof(integer));
	i_2 = nodel;
	for (j = 1; j <= i_2; ++j) {
	    do_fio(&c__1, (char *)&eltop[elnum + (j + 2) * 100 - 101], (
		    ftnlen)sizeof(integer));
	}
	e_wsfe();
	eltop[elnum - 1] = eltyp;
	eltop[elnum + 99] = nodel;
/* L1020: */
    }

/*                            INPUT OF PERMEABILITIES, CON- */
/*                            STRUCTION OF PERMEABILITY MATRIX P AND */
/*                            SOURCE STRENGTH */

    io__62.ciunit = nout;
    s_wsfe(&io__62);
    e_wsfe();
    matnul_(p, &ip, &jp, &dimen, &dimen, &itest);
    io__64.ciunit = nin;
    s_rsfe(&io__64);
    i_1 = dimen;
    for (i = 1; i <= i_1; ++i) {
	do_fio(&c__1, (char *)&p[i + i * 3 - 4], (ftnlen)sizeof(doublereal));
    }
    e_rsfe();
    io__65.ciunit = nout;
    s_wsfe(&io__65);
    i_1 = dimen;
    for (i = 1; i <= i_1; ++i) {
	do_fio(&c__1, (char *)&p[i + i * 3 - 4], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();

    io__66.ciunit = nout;
    s_wsfe(&io__66);
    e_wsfe();
    io__67.ciunit = nin;
    s_rsfe(&io__67);
    do_fio(&c__1, (char *)&strgth, (ftnlen)sizeof(doublereal));
    e_rsfe();
    io__69.ciunit = nout;
    s_wsfe(&io__69);
    do_fio(&c__1, (char *)&strgth, (ftnlen)sizeof(doublereal));
    e_wsfe();

/*                            INPUT OF NUMBER OF DEGREES OF FREEDOM */
/*                            PER NODE, INPUT OF BOUNDARY CON- */
/*                            DITIONS AND CONSTRUCTION OF NODAL */
/*                            FREEDOM ARRAY NF */

    io__70.ciunit = nout;
    s_wsfe(&io__70);
    e_wsfe();
    io__71.ciunit = nin;
    s_rsfe(&io__71);
    do_fio(&c__1, (char *)&dofnod, (ftnlen)sizeof(integer));
    e_rsfe();
    io__73.ciunit = nout;
    s_wsfe(&io__73);
    do_fio(&c__1, (char *)&dofnod, (ftnlen)sizeof(integer));
    e_wsfe();
    io__74.ciunit = nin;
    s_rsfe(&io__74);
    do_fio(&c__1, (char *)&bndnod, (ftnlen)sizeof(integer));
    e_rsfe();
    io__76.ciunit = nout;
    s_wsfe(&io__76);
    do_fio(&c__1, (char *)&bndnod, (ftnlen)sizeof(integer));
    e_wsfe();
    i_1 = bndnod;
    for (i = 1; i <= i_1; ++i) {
	io__77.ciunit = nin;
	s_rsfe(&io__77);
	do_fio(&c__1, (char *)&bnode[i - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&bval[i - 1], (ftnlen)sizeof(doublereal));
	e_rsfe();
	io__80.ciunit = nout;
	s_wsfe(&io__80);
	do_fio(&c__1, (char *)&bnode[i - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&bval[i - 1], (ftnlen)sizeof(doublereal));
	e_wsfe();
/* L1030: */
    }
    totdof = 0;
    i_1 = totnod;
    for (i = 1; i <= i_1; ++i) {
	i_2 = dofnod;
	for (j = 1; j <= i_2; ++j) {
	    ++totdof;
	    nf[i + j * 100 - 101] = totdof;
/* L1040: */
	}
/* L1050: */
    }

/*                            CALCULATION OF SEMI-BANDWIDTH */

    first = TRUE_;
    i_1 = totels;
    for (nele = 1; nele <= i_1; ++nele) {
	fredif_(&nele, eltop, &ieltop, &jeltop, nf, &inf, &jnf, &dofnod, &
		first, &dif, &itest);
/* L1060: */
    }
    hband = dif + 1;

/* *                           ************************************ */
/* *                           *                                  * */
/* *                           * SYSTEM STIFFNESS MATRIX ASSEMBLY * */
/* *                           *                                  * */
/* *                           ************************************ */

    vecnul_(rhs, &irhs, &totdof, &itest);
    matnul_(sysk, &isysk, &jsysk, &totdof, &hband, &itest);
    dofel = nodel * dofnod;
    qqua4_(wght, &iwght, abss, &iabss, &jabss, &nqp, &itest);
    i_1 = totels;
    for (nele = 1; nele <= i_1; ++nele) {
	elgeom_(&nele, eltop, &ieltop, &jeltop, coord, &icoord, &jcoord, geom,
		 &igeom, &jgeom, &dimen, &itest);

/*                            INTEGRATION LOOP FOR ELEMENT STIFFNESS 
*/
/*                            USING NQP QUADRATURE POINTS */

	matnul_(elk, &ielk, &jelk, &dofel, &dofel, &itest);
	vecnul_(elq, &ielq, &dofel, &itest);
	i_2 = nqp;
	for (iquad = 1; iquad <= i_2; ++iquad) {

/*                            FORM LINEAR SHAPE FUNCTION AND 
SPACE */
/*                            DERIVATIVES IN THE LOCAL 
CORRDINATES. */
/*                            TRANSFORM LOCAL DERIVATIVES TO 
GLOBAL */
/*                            COORDINATE SYSTEM */

	    xi = abss[iquad * 3 - 3];
	    eta = abss[iquad * 3 - 2];
	    quam4_(fun, &ifun, lder, &ilder, &jlder, &xi, &eta, &itest);

	    scaprd_(geom, &igeom, fun, &ifun, &nodel, &x, &itest);
	    scaprd_(&geom[8], &igeom, fun, &ifun, &nodel, &y, &itest);

	    matmul_(lder, &ilder, &jlder, geom, &igeom, &jgeom, jac, &ijac, &
		    jjac, &dimen, &nodel, &dimen, &itest);
	    matinv_(jac, &ijac, &jjac, jacin, &ijacin, &jjacin, &dimen, &det, 
		    &itest);
	    matmul_(jacin, &ijacin, &jjacin, lder, &ilder, &jlder, gder, &
		    igder, &jgder, &dimen, &dimen, &nodel, &itest);

/*                            FORMATION OF ELEMENT STIFFNESS ELK 
*/

	    matmul_(p, &ip, &jp, gder, &igder, &jgder, pd, &ipd, &jpd, &dimen,
		     &dimen, &dofel, &itest);
	    matran_(gder, &igder, &jgder, gdert, &igdert, &jgdert, &dimen, &
		    dofel, &itest);
	    matmul_(gdert, &igdert, &jgdert, pd, &ipd, &jpd, dtpd, &idtpd, &
		    jdtpd, &dofel, &dimen, &dofel, &itest);
	    quot = abs(det) * wght[iquad - 1];
	    i_3 = dofel;
	    for (i = 1; i <= i_3; ++i) {
		i_4 = dofel;
		for (j = 1; j <= i_4; ++j) {
		    dtpd[i + (j << 3) - 9] *= quot;
/* L1070: */
		}
		scvec[i - 1] = fun[i - 1] * src_(&x, &y, &strgth) * quot;
/* L1080: */
	    }
	    matadd_(elk, &ielk, &jelk, dtpd, &idtpd, &jdtpd, &dofel, &dofel, &
		    itest);
	    vecadd_(elq, &ielq, scvec, &iscvec, &dofel, &itest);
/* L1090: */
	}

/*                            ASSEMBLY OF SYSTEM STIFFNESS MATRIX */

	direct_(&nele, eltop, &ieltop, &jeltop, nf, &inf, &jnf, &dofnod, 
		steer, &isteer, &itest);
	assym_(sysk, &isysk, &jsysk, elk, &ielk, &jelk, steer, &isteer, &
		hband, &dofel, &itest);
	asrhs_(rhs, &irhs, elq, &ielq, steer, &isteer, &dofel, &itest);
/* L1100: */
    }

/* *                           ********************* */
/* *                           *                   * */
/* *                           * EQUATION SOLUTION * */
/* *                           *                   * */
/* *                           ********************* */

/*                            MODIFICATION OF STIFFNESS MATRIX AND */
/*                            RIGHT-HAND SIDE TO IMPLEMENT BOUNDARY */
/*                            CONDITIONS */

    i_1 = bndnod;
    for (i = 1; i <= i_1; ++i) {
	j = bnode[i - 1];
	sysk[j + hband * 100 - 101] *= scale;
	rhs[j - 1] = sysk[j + hband * 100 - 101] * bval[i - 1];
/* L1110: */
    }

/*                            SOLUTION OF SYSTEM MATRIX FOR THE */
/*                            NODAL VALUES OF THE POTENTIAL */

    chosol_(sysk, &isysk, &jsysk, rhs, &irhs, &totdof, &hband, &itest);
    io__113.ciunit = nout;
    s_wsfe(&io__113);
    e_wsfe();
    prtval_(rhs, &irhs, nf, &inf, &jnf, &dofnod, &totnod, &nout, &itest);
    s_stop("", 0L);
} /* MAIN__ */


doublereal src_(x, y, strgth)
doublereal *x, *y, *strgth;
{
    /* System generated locals */
    doublereal ret_val;

    ret_val = 0.;
    if (*x > 1. && *x < 2. && *y > 1. && *y < 2.) {
	ret_val = *strgth;
    }
    return ret_val;
} /* src_ */

