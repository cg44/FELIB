
MODULE mod_prtval

USE FELIB_GLOBALS,only : wp, stdout

PRIVATE
PUBLIC prtval

contains

SUBROUTINE prtval(val,nf,nout,itest)

!-----------------------------------------------------------------------
! PURPOSE
!      PRTVAL prints out the nodal values of the solution in a
!      standard format
!
! HISTORY
!
!      Copyright (C) 2000 : CCLRC, Rutherford Appleton Laboratory
!                           Chilton, Didcot, Oxfordshire OX11 0QX
!
!      Release 1.0   2 Jul 2000 (CG)
!
! ARGUMENTS in
!      VAL     vector of dimension IVAL containing solution
!              values,and prescribed boundary values
!      NF      integer array of dimension (INF, JNF) containing
!              freedom numbers associated with each NODE
!      DOFNOD  number of degrees of freedom at each NODE
!      TOTNOD  total number of nodes in mesh
!      NOUT    fortran unit number
!
! ROUTINES called
!
!***********************************************************************

IMPLICIT NONE

! Dummy arguments

INTEGER :: dofnod , totnod
INTEGER , OPTIONAL :: nout, itest
INTEGER , POINTER :: nf(:,:)
REAL (wp), POINTER :: val(:)


! Local variables

INTEGER :: i , ierror , j , jtest , k , n, lnout, ktest
INTEGER :: node(5)
REAL (wp)  :: work(5)
CHARACTER(6) :: srname

INTRINSIC abs
!
DATA srname/'PRTVAL'/

!
!     Main body
!
if(present(itest))then
ktest=itest
else
ktest=0
endif
if(present(nout))then
lnout=nout
else
lnout=stdout
endif
dofnod=size(nf,2)
totnod=size(nf,1)

IF ( dofnod==1 ) THEN
!
   WRITE (lnout,FMT=99005)
   n = 0
   DO i = 1 , totnod
       n = n + 1
       node(n) = i
       k = nf(i,1)
       work(n) = 0.0D0
       IF ( k/=0 ) THEN
           IF ( k<=0 ) k = abs(k)
           work(n) = val(k)
       END IF
       IF ( n==4 ) THEN
           WRITE (lnout,FMT=99007) (node(j),work(j),j=1,4)
           n = 0
       END IF
   END DO
!
   IF ( n/=0 ) WRITE (nout,FMT=99007) (node(j),work(j),j=1,n)
ELSE
   IF ( dofnod==2 ) WRITE (lnout,FMT=99001)
   IF ( dofnod==3 ) WRITE (lnout,FMT=99002)
   IF ( dofnod==4 ) WRITE (lnout,FMT=99003)
   IF ( dofnod>=5 ) WRITE (lnout,FMT=99004)
   n = 0
   DO i = 1 , totnod
       DO j = 1 , dofnod
           n = n + 1
           node(n) = i
           k = nf(i,j)
!
           work(n) = 0.0D0
           IF ( k/=0 ) THEN
               IF ( k<=0 ) k = abs(k)
               work(n) = val(k)
           END IF
           IF ( ((dofnod/=2) .OR. (n==4)) .AND. (dofnod==2) ) THEN
               WRITE (lnout,FMT=99008) node(1) , work(1) , work(2) ,     &
     &                                node(3) , work(3) , work(4)
               n = 0
           ELSE IF ( (dofnod/=2) .AND. (n==dofnod) ) THEN
               WRITE (nout,FMT=99006) i , (work(k),k=1,dofnod)
               n = 0
           END IF
       END DO
   END DO
!
   IF ( n/=0 ) THEN
       IF ( dofnod==2 ) WRITE (lnout,FMT=99008) node(1) , work(1) ,      &
     &                         work(2)
       IF ( dofnod/=2 ) WRITE (lnout,FMT=99006) i , (work(k),k=1,n)
   END IF
END IF
!
99001 FORMAT (/' ',2('NODE',6X,2('VALUE',9X)),/' ')
99002 FORMAT (/' NODE',6X,3('VALUE',9X),/' ')
99003 FORMAT (/' NODE',6X,4('VALUE',9X),/' ')
99004 FORMAT (/' NODE',6X,5('VALUE',9X),/' ')
99005 FORMAT (/' ',4('NODE',5X,'VALUE',5X),/' ')
99006 FORMAT (' ',I4,1X,5(D12.5,2X),/' ',5X,5(D12.5,2X),/' ',5X,        &
     &        5(D12.5,2X))
99007 FORMAT (' ',4(I4,1X,D12.5,2X))
99008 FORMAT (' ',2(I4,1X,2(D12.5,2X),5X))

END SUBROUTINE PRTVAL
END MODULE mod_prtval
