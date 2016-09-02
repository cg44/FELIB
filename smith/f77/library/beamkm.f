      SUBROUTINE BEAMKM(KM,EI,ELL)
C
C      THIS SUBROUTINE FORMS THE STIFFNESS MATRIX OF A
C      HORIZONTAL BEAM ELEMENT(BENDING ONLY)
C
      REAL KM(4,4)
      KM(1,1)=12.*EI/(ELL*ELL*ELL)
      KM(3,3)=KM(1,1)
      KM(1,2)=6.*EI/(ELL*ELL)
      KM(2,1)=KM(1,2)
      KM(1,4)=KM(1,2)
      KM(4,1)=KM(1,4)
      KM(1,3)=-KM(1,1)
      KM(3,1)=KM(1,3)
      KM(3,4)=-KM(1,2)
      KM(4,3)=KM(3,4)
      KM(2,3)=KM(3,4)
      KM(3,2)=KM(2,3)
      KM(2,2)=4.*EI/ELL
      KM(4,4)=KM(2,2)
      KM(2,4)=2.*EI/ELL
      KM(4,2)=KM(2,4)
      RETURN
      END
