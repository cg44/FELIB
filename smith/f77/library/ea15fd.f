      SUBROUTINE EA15FD(N,ITAPE,NLAN,IO,V)
C STORE OR RECOVER A LANCZOS VECTOR.
      REAL             V(N)
C N IS THE ORDER OF THE MATRIX (ALSO THE ORDER OF THE VECTOR BEING
C     STORED OR RECOVERED). IT MUST NOT BE ALTERED.
C ITAPE IS THE UNIT NUMBER OF THE TAPE UNIT. IT IS NOT ALTERED.
C NLAN IS THE LANCZOS ITERATION NUMBER. THE VECTORS ARE ALWAYS SENT FOR
C     STORAGE IN ORDER (NLAN=1,2,...) AND ARE RECOVERED IN ORDER
C     (NLAN=1,2,...). IF THE FIRST VECTOR (NLAN=1) IS SENT FOR STORAGE
C     THEN IT CAN BE ASSUMED THAT A NEW PROBLEM IS IN HAND AND THE OLD
C     VECTORS MAY BE OVERWRITTEN. NLAN MUST NOT BE ALTERED.
C IO IS 1 IF A VECTOR IS TO BE STORED AND TO 2 IF ONE IS TO BE
C     RECOVERED. IT MUST NOT BE ALTERED.
C V HOLDS THE VECTOR TO BE STORED OR IS SET TO THE VECTOR RECOVERED.
      COMMON/EA15HD/LAST
C LAST IS THE VALUE OF NLAN ON THE LAST CALL.
      IF(ITAPE.LE.0)GO TO 50
      IF(NLAN.NE.1)GO TO 10
      REWIND ITAPE
      LAST=0
C IF NECESSARY, SKIP OVER THE VECTORS WRITTEN EARLIER
   10 NGAP=NLAN-LAST-1
      IF(NGAP.LE.0)GO TO 30
      DO 20 I=1,NGAP
      READ(ITAPE)
   20 CONTINUE
C PERFORM ACTUAL READ OR WRITE
   30 IF(IO.EQ.1)WRITE(ITAPE)V
      IF(IO.EQ.2)READ(ITAPE)V
   40 LAST=NLAN
   50 RETURN
      END