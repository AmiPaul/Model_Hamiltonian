      PROGRAM PWR_SPEC 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER NSTEP, I, K
      PARAMETER (NSTEP=50000)

      DOUBLE PRECISION DT
      DOUBLE PRECISION QSAVE(NSTEP)
      DOUBLE PRECISION POWER(NSTEP/2)
      DOUBLE PRECISION FREQ(NSTEP/2)

      DT = 0.0002D0

      DO I = 1,NSTEP
      READ(5,*)A,B,C,QSAVE(I)
      ENDDO

      CALL POWER_SPEC(QSAVE,NSTEP,DT,FREQ,POWER)

      OPEN(20,FILE='Q_power.dat',STATUS='REPLACE')

      DO 200 K=1,NSTEP/2
         WRITE(20,'(2E20.10)') FREQ(K), POWER(K)
 200  CONTINUE

      CLOSE(20)
      END


      SUBROUTINE POWER_SPEC(Q,N,DT,FREQ,POWER)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER N,K,J
      DOUBLE PRECISION Q(N),FREQ(N/2),POWER(N/2)

      DOUBLE PRECISION PI,QAVG,QR
      DOUBLE PRECISION RE,IM,ANGLE

      PI = 4.0D0*DATAN(1.0D0)

C----- Remove DC component
      QAVG = 0.0D0
      DO 10 J=1,N
         QAVG = QAVG + Q(J)
 10   CONTINUE
      QAVG = QAVG/DFLOAT(N)

C----- DFT
      DO 100 K=1,N/2

         RE = 0.0D0
         IM = 0.0D0

         DO 50 J=1,N

            QR = Q(J) - QAVG

            ANGLE = 2.0D0*PI*DFLOAT(K-1)*DFLOAT(J-1)
     &             /DFLOAT(N)

            RE = RE + QR*DCOS(ANGLE)
            IM = IM - QR*DSIN(ANGLE)

 50      CONTINUE

         POWER(K) = RE*RE + IM*IM

C----- Frequency in 1/ps
         FREQ(K) = DFLOAT(K-1)/(DFLOAT(N)*DT)

 100  CONTINUE

      RETURN
      END      
