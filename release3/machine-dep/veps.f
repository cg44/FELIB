C***********************************************************************
      DOUBLE PRECISION FUNCTION VEPS(X)
C-----------------------------------------------------------------------
C PURPOSE
C      RETURNS THE VALUE OF THE SMALLEST NUMBER VEPS SUCH THAT
C      1.D0+VEPS > 1.DO ON THIS COMPUTER
C
C      *********************************************
C      **********MACHINE DEPENDENT ROUTINE**********
C      *********************************************
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 1.1  29 OCT 1979 (CG)
C      COMMENTED    23 OCT 1980 (KR)
C      EXTENDED FOR SUN AND STELLAR 29/6/1989, PYRAMID 11/10/1989
C
C ARGUMENTS IN
C      X       DUMMY ARGUMENT
C
C
C     DOUBLE PRECISION FUNCTION VEPS(X)
C***********************************************************************
C
      DOUBLE PRECISION X, VALUE
C
C     DATA VALUE /1.D13/
C
C+PRIME EPS = 2.0D0**(-44)
C      INTEGER*2  L(4)
C      EQUIVALENCE (VALUE,L(1))
C      DATA L(1),L(2),L(3),L(4)/:040000,:000000,:000000,:000125/
C+GEC+IBM
C+GEC DATA VALUE /Z3410000000000000/
C+SUN+STELLAR EPS=2.0D0**(-52)
C+SUN ON SUN THIS IS SPECIFIED AS
C     DATA VALUE /Z'3CB0000000000000'/
C+STELLAR ON STELLAR IT IS SPECIFIED AS
C     DATA VALUE /'3CB0000000000000'X/
C     PYRAMID DOESN'T RECOGNISE EITHER OF THESE HEXADECIMAL FORMS, USE:
C
      VALUE = 2.0D-17
C
      VEPS = VALUE
      RETURN
      END
C***********************************************************************
