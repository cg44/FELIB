      SUBROUTINE FORMXI(FSOIL,FMAX,RF,RM,R0,XI)
C
C      THIS SUBROUTINE FORMS PART OF THE SPRING STIFFNESS TERM
C      FOR PROGRAM 12.2
C
      PHI=FSOIL*R0*RF/FMAX
      XI=LOG((RM-PHI)/(R0-PHI))+PHI*(RM-R0)/((RM-PHI)*(R0-PHI))
      RETURN
      END