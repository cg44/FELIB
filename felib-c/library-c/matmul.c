# 6 "matmul.cpp"
# 1 "f2c.h" 1
# 10 "f2c.h"
typedef long int integer;
typedef char *address;
typedef short int shortint;
typedef float real;
typedef double doublereal;
typedef struct { real r, i; } complex;
typedef struct { doublereal r, i; } doublecomplex;
typedef long int logical;
typedef short int shortlogical;
# 36 "f2c.h"
typedef long flag;
typedef long ftnlen;
typedef long ftnint;



typedef struct
{ flag cierr;
        ftnint ciunit;
        flag ciend;
        char *cifmt;
        ftnint cirec;
} cilist;


typedef struct
{ flag icierr;
        char *iciunit;
        flag iciend;
        char *icifmt;
        ftnint icirlen;
        ftnint icirnum;
} icilist;


typedef struct
{ flag oerr;
        ftnint ounit;
        char *ofnm;
        ftnlen ofnmlen;
        char *osta;
        char *oacc;
        char *ofm;
        ftnint orl;
        char *oblnk;
} olist;


typedef struct
{ flag cerr;
        ftnint cunit;
        char *csta;
} cllist;


typedef struct
{ flag aerr;
        ftnint aunit;
} alist;


typedef struct
{ flag inerr;
        ftnint inunit;
        char *infile;
        ftnlen infilen;
        ftnint *inex;
        ftnint *inopen;
        ftnint *innum;
        ftnint *innamed;
        char *inname;
        ftnlen innamlen;
        char *inacc;
        ftnlen inacclen;
        char *inseq;
        ftnlen inseqlen;
        char *indir;
        ftnlen indirlen;
        char *infmt;
        ftnlen infmtlen;
        char *inform;
        ftnint informlen;
        char *inunf;
        ftnlen inunflen;
        ftnint *inrecl;
        ftnint *innrec;
        char *inblank;
        ftnlen inblanklen;
} inlist;



union Multitype {
        shortint h;
        integer i;
        real r;
        doublereal d;
        complex c;
        doublecomplex z;
        };

typedef union Multitype Multitype;

typedef long Long;

struct Vardesc {
        char *name;
        char *addr;
        ftnlen *dims;
        int type;
        };
typedef struct Vardesc Vardesc;

struct Namelist {
        char *name;
        Vardesc **vars;
        int nvars;
        };
typedef struct Namelist Namelist;
# 157 "f2c.h"
typedef int (*U_fp)(...);
typedef shortint (*J_fp)(...);
typedef integer (*I_fp)(...);
typedef real (*R_fp)(...);
typedef doublereal (*D_fp)(...), (*E_fp)(...);
typedef void (*C_fp)(...);
typedef void (*Z_fp)(...);
typedef logical (*L_fp)(...);
typedef shortlogical (*K_fp)(...);
typedef void (*H_fp)(...);
typedef int (*S_fp)(...);
# 182 "f2c.h"
typedef void C_f;
typedef void H_f;
typedef void Z_f;
typedef doublereal E_f;
# 7 "matmul.cpp" 2


                 int matmul_(a, ia, ja, b, ib, jb, c, ic, jc, l, m, n, itest)
doublereal *a;
integer *ia, *ja;
doublereal *b;
integer *ib, *jb;
doublereal *c;
integer *ic, *jc, *l, *m, *n, *itest;
{


    static char srname[6+1] = "MATMUL";


    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2,
            i__3;


    static integer i, j, k;
    static doublereal x;
    extern integer errmes_();
    static integer ierror;
# 73 "matmul.cpp"
    c_dim1 = *ic;
    c_offset = c_dim1 + 1;
    c -= c_offset;
    b_dim1 = *ib;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    a_dim1 = *ia;
    a_offset = a_dim1 + 1;
    a -= a_offset;





    if (*itest != -1) {
        ierror = 0;
        if (*l > *ic || *n > *jc) {
            ierror = 4;
        }
        if (*m > *ib || *n > *jb) {
            ierror = 3;
        }
        if (*l > *ia || *m > *ja) {
            ierror = 2;
        }
        if (*l <= 0 || *m <= 0 || *n <= 0) {
            ierror = 1;
        }
        *itest = errmes_(itest, &ierror, srname, 6L);
        if (*itest != 0) {
            return 0;
        }
    }



    i__1 = *l;
    for (i = 1; i <= i__1; ++i) {
        i__2 = *n;
        for (j = 1; j <= i__2; ++j) {
            x = 0.;
            i__3 = *m;
            for (k = 1; k <= i__3; ++k) {
                x += a[i + k * a_dim1] * b[k + j * b_dim1];

            }
            c[i + j * c_dim1] = x;

        }

    }

}
