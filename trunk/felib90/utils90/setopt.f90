 MODULE mod_setopt

! MODULE setopt provides a simple way of setting the initial
! value of an OPTIONAL argument in a procedure call. The real
! for SETOPT is to aviod multiple repititions of the same
! sequence of tests at the head of a routine.

! System

    USE felib_globals, ONLY : wp

    PRIVATE
    PUBLIC setopt

! Generic interfaces

    INTERFACE setopt
       MODULE PROCEDURE complex_setopt
       MODULE PROCEDURE real_setopt
       MODULE PROCEDURE integer_setopt
       MODULE PROCEDURE character_setopt
    END INTERFACE

 CONTAINS

    SUBROUTINE complex_setopt(dummy,value,actual)

       IMPLICIT NONE
       INTRINSIC present
       COMPLEX (wp) :: dummy, value
       COMPLEX (wp), OPTIONAL :: actual

       IF (present(actual)) THEN
          dummy = actual
       ELSE
          dummy = value
       END IF


    END SUBROUTINE complex_setopt

    SUBROUTINE real_setopt(dummy,value,actual)

       IMPLICIT NONE
       INTRINSIC present
       REAL (wp) :: dummy, value
       REAL (wp), OPTIONAL :: actual

       IF (present(actual)) THEN
          dummy = actual
       ELSE
          dummy = value
       END IF


    END SUBROUTINE real_setopt

    SUBROUTINE integer_setopt(dummy,value,actual)

       IMPLICIT NONE
       INTRINSIC present
       INTEGER :: dummy, value
       INTEGER, OPTIONAL :: actual


       IF (present(actual)) THEN
          dummy = actual
       ELSE
          dummy = value
       END IF


    END SUBROUTINE integer_setopt

    SUBROUTINE character_setopt(dummy,value,actual)

       IMPLICIT NONE
       CHARACTER (*) :: dummy, value
       CHARACTER (*), OPTIONAL :: actual
       INTRINSIC present

       IF (present(actual)) THEN
          dummy = actual
       ELSE
          dummy = value
       END IF

    END SUBROUTINE character_setopt

 END MODULE mod_setopt

