      SUBROUTINE BEAMKP(KP,ELL)
C
C      THIS SUBROUTINE FORMS THE TERMS OF THE BEAM STIFFNESS
C      MATRIX DUE TO AXIAL LOADING
C
      REAL KP(4,4)
      KP(1,1)=1.2/ELL
      KP(1,2)=0.1
      KP(2,1)=0.1
      KP(1,3)=-1.2/ELL
      KP(3,1)=-1.2/ELL
      KP(1,4)=0.1
      KP(4,1)=0.1
      KP(2,2)=2.0*ELL/15.0
      KP(2,3)=-0.1
      KP(3,2)=-0.1
      KP(2,4)=-ELL/30.0
      KP(4,2)=-ELL/30.0
      KP(3,3)=1.2/ELL
      KP(3,4)=-0.1
      KP(4,3)=-0.1
      KP(4,4)=2.0*ELL/15.0
      RETURN
      END
