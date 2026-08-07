      PROGRAM HBOND_REAL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER NSTEP, I
      PARAMETER (NT = 200,NS=15000)
      PARAMETER (C1=418.4)
      DOUBLE PRECISION DT

C --- Coordinates (Å) and momenta (amu*Å/ps)
      DOUBLE PRECISION qD, qA, Q
      DOUBLE PRECISION pD, pA, pQ
      DOUBLE PRECISION A, phi
      DOUBLE PRECISION omega, rnd, pi

C --- Forces (kcal/mol/Å)
      DOUBLE PRECISION fD, fA, fQ

C --- Masses (amu)
      DOUBLE PRECISION mD, mA, mQ

C --- Morse params (kcal/mol, Å^-1)
      DOUBLE PRECISION DD, DA, aD, aA

C --- HB force constant (kcal/mol/Å^2)
      DOUBLE PRECISION DHB, aHB

C --- Couplings (kcal/mol/Å^2)
      DOUBLE PRECISION kappa, delta, alphaD, alphaA, betaA, betaD

C ---  (energy currents
      DOUBLE PRECISION qqA(NT,NS),qqD(NT,NS),QQ(NT,NS)
      DOUBLE PRECISION EED(NT,NS),EEA(NT,NS),EEHB(NT,NS)
      DOUBLE PRECISION HH(NT,NS),VV(NT,NS),TT(NT,NS)

C --- Energies
      DOUBLE PRECISION ED_IVR(NT,NS),EA_IVR(NT,NS)
      DOUBLE PRECISION EA_HB(NT,NS),ED_HB(NT,NS)
      DOUBLE PRECISION JHB(NT,NS),JD(NT,NS),JA(NT,NS)
      DOUBLE PRECISION JJHB(NT),JJD(NT),JJA(NT)

C-------------------------------
      DOUBLE PRECISION H,V,T
      DOUBLE PRECISION JHBENS, JDENS, JAENS      
      DOUBLE PRECISION qAA, qDD, QQQ, EDD, EAA, EHHB
      DOUBLE PRECISION EDD_IVR, EAA_IVR, EDD_HB, EAA_HB
      DOUBLE PRECISION JHHB, JDD, JAA
      DOUBLE PRECISION HHT, VVT, TTT

      DT = 0.0002D0
      PI = 4.0D0*DATAN(1.0D0)

C --- Masses
      mD = 1.0D0
      mA = 1.0D0
      mQ = 10.0D0

C --- Morse
      DHB = 4.0D0*C1
      DD = 110.0D0*C1
      DA = 110.0D0*C1
      aHB = 1.0D0
      aD = 2.2D0
      aA = 2.2D0
      

C --- Couplings
      kappa  = 1.0D0*C1
      delta  = 3.5D0*C1
      alphaD = 0.5D0*C1
      alphaA = 2.5D0*C1
      betaA = 0.0D0*C1
      betaD = 0.0D0*C1



      A = 0.1205d0  ! Angstrom (for v = 1 of D)
      omega = 668 !ps-1

      JHBENS = 0.0D0
      JDENS = 0.0D0
      JAENS = 0.0D0

      OPEN(10, FILE='traj_real_D.dat',status="replace")
      OPEN(11, FILE='traj_JD.dat',status="replace")
      OPEN(12, FILE='traj_ED_Comp.dat',status="replace")
      OPEN(13, FILE='traj_ET_D.dat',status="replace")

      DO 200 ITRAJ = 1, NT
      JJHB(ITRAJ)= 0.0d0
      JJD(ITRAJ)= 0.0d0
      JJA(ITRAJ)= 0.0d0

C --- Initial conditions
      qD = 0.0D0
      qA = 0.0D0
      Q  = 0.0D0

      pD = 0.0D0
      pA = 0.0D0   
      pQ = 0.0D0

      call random_number(rnd)
      phi = 2.0d0*pi*rnd

c     qA = A*cos(phi)
c     pA = mA*omega*A*sin(phi)
      qD = A*cos(phi)
      pD = mD*omega*A*sin(phi)
c     qD = 0.0d0
c     pD = 8.0d0



         CALL FORCES(qD,qA,Q,fD,fA,fQ,DHB,DD,DA,aHB,aD,aA,
     &         kappa,delta,alphaD,alphaA,betaA,betaD)


      DO 100 I=1,NS

C --- Half step momenta
         pD = pD + 0.5D0*DT*fD
         pA = pA + 0.5D0*DT*fA
         pQ = pQ + 0.5D0*DT*fQ

C --- Positions
         qD = qD + DT*pD/mD
         qA = qA + DT*pA/mA
         Q  = Q  + DT*pQ/mQ


C --- New forces
         CALL FORCES(qD,qA,Q,fD,fA,fQ,DHB,DD,DA,aHB,aD,aA,
     &         kappa,delta,alphaD,alphaA,betaA,betaD)

C --- Complete momenta
         pD = pD + 0.5D0*DT*fD
         pA = pA + 0.5D0*DT*fA
         pQ = pQ + 0.5D0*DT*fQ

C --- Positions
c        qD = qD + DT*pD/mD
c        qA = qA + DT*pA/mA
c        Q  = Q  + DT*pQ/mQ

C --- Energies (kcal/mol)
         ED  = 0.5D0*pD*pD/mD
         EA  = 0.5D0*pA*pA/mA
         EHB = 0.5D0*pQ*pQ/mQ

C----------------ENERGIES
      CALL ENERGY (pD,pA,pQ,qD,qA,Q,DHB,DD,DA,aHB,aD,aA,
     & mA,mD,mQ,kappa,delta,alphaD,alphaA,betaA,betaD,H,V,T)    
         HH(ITRAJ,I)=H
         VV(ITRAJ,I)=V
         TT(ITRAJ,I)=T
c        write(7,*)I*DT,HH(itraj,i),vv(itraj,i),tt(itraj,i)

C----SUMMING UP--------------
         qqA(ITRAJ,I)=qA
         qqD(ITRAJ,I)=qD
         QQ(ITRAJ,I)=Q
         EED(ITRAJ,I)=ED         
         EEA(ITRAJ,I)=EA
         EEHB(ITRAJ,I)=EHB
C---- IVR COMPONENT OF ENERGY CURRENT---         
      ED_IVR(ITRAJ,I)=kappa*qA+2.0*delta*qA*qD
      EA_IVR(ITRAJ,I)=kappa*qD+delta*qD*qD

C---- HB COMPONENT OF ENERGY CURRENT---         
      ED_HB(ITRAJ,I)=alphaD*Q+2.0*betaD*qD*Q
      EA_HB(ITRAJ,I)=alphaA*Q+2.0*betaA*qA*Q

C-----Energy current(J)-----
      JHB(ITRAJ,I)=-pQ*(alphaA*qA+alphaD*qD+betaA*qA*qA+betaD*qD*qD)/mQ   
      JD(ITRAJ,I)=-pD*(kappa*qA+2.0*delta*qA*qD+alphaD*Q+
     &             2.0*betaD*qD*Q)/mD 
      JA(ITRAJ,I)=-pA*(kappa*qD+delta*qD*qD+alphaA*Q+2.0*betaA*qA*Q)/mA      
      JJHB(ITRAJ)=JJHB(ITRAJ)+JHB(ITRAJ,I)
      JJD(ITRAJ)=JJD(ITRAJ)+JD(ITRAJ,I)
      JJA(ITRAJ)=JJA(ITRAJ)+JA(ITRAJ,I)
c     write(6,*)qA,qD,Q,ED,ED_IVR,ED_HB,JJD
c     WRITE(12+ITRAJ,*)qA,qD,Q,ED,EA,EHB

 100  CONTINUE
c     write(7,*)

C-----AVERAGING----------------      
      JJHB(ITRAJ)=JJHB(ITRAJ)/NS
      JJD(ITRAJ)=JJD(ITRAJ)/NS
      JJA(ITRAJ)=JJA(ITRAJ)/NS
      JHBENS = JHBENS+JJHB(ITRAJ)
      JDENS = JDENS+JJD(ITRAJ)
      JAENS = JAENS+JJA(ITRAJ)

 200  CONTINUE
      JHBENS=JHBENS/NT
      JDENS = JDENS/NT
      JAENS = JAENS/NT
      WRITE(6,*)"Averaged JHB, JD, and JA  ",JHBENS/C1,JDENS/C1,JAENS/C1

      DO 300 ISTEP=1,NS
        qAA=0.0d0
        qDD=0.0d0
        QQQ=0.0d0
        EDD=0.0d0
        EAA=0.0d0
        EHHB=0.0d0
        EDD_IVR=0.0d0
        EAA_IVR=0.0d0
        EDD_HB=0.0d0
        EAA_HB=0.0d0
        JHHB=0.0d0
        JDD=0.0d0
        JAA=0.0d0
        HHT=0.0d0
        VVT=0.0d0
        TTT=0.0d0
      DO 400 ITRAJ=1,NT
        qAA=qAA+qqA(ITRAJ,ISTEP)
        qDD=qDD+qqD(ITRAJ,ISTEP)
        QQQ=QQQ+QQ(ITRAJ,ISTEP)
        EDD=EDD+EED(ITRAJ,ISTEP)
        EAA=EAA+EEA(ITRAJ,ISTEP)
        EHHB=EHHB+EEHB(ITRAJ,ISTEP)
        EDD_IVR=EDD_IVR+ED_IVR(ITRAJ,ISTEP)
        EAA_IVR=EAA_IVR+EA_IVR(ITRAJ,ISTEP)
        EDD_HB=EDD_HB+ED_HB(ITRAJ,ISTEP)
        EAA_HB=EAA_HB+EA_HB(ITRAJ,ISTEP)
        JHHB=JHHB+JHB(ITRAJ,ISTEP)
        JDD=JDD+JD(ITRAJ,ISTEP)
        JAA=JAA+JA(ITRAJ,ISTEP)
        HHT=HHT+HH(ITRAJ,ISTEP)
        VVT=VVT+VV(ITRAJ,ISTEP)
        TTT=TTT+TT(ITRAJ,ISTEP)
 400  CONTINUE
        qAA=qAA/NT
        qDD=qDD/NT
        QQQ=QQQ/NT  
        EDD=EDD/NT
        EAA=EAA/NT
        EHHB=EHHB/NT
        EDD_IVR=EDD_IVR/NT
        EAA_IVR=EAA_IVR/NT
        EDD_HB=EDD_HB/NT
        EAA_HB=EAA_HB/NT
        JHHB=JHHB/NT
        JDD=JDD/NT
        JAA=JAA/NT
        HHT=HHT/NT
        VVT=VVT/NT
        TTT=TTT/NT
         WRITE(10,*) ISTEP*DT, qDD, qAA, QQQ, EDD/C1, EAA/C1, EHHB/C1
         WRITE(11,*) ISTEP*DT, JAA/C1, JDD/C1, JHHB/C1
         WRITE(12,*) ISTEP*DT,EDD_IVR/C1,EAA_IVR/C1,EDD_HB/C1,EAA_HB/C1 
         WRITE(13,*) ISTEP*DT,HHT/C1,VVT/C1,TTT/C1
 300  CONTINUE
      CLOSE(10)
      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      END


C =====================================================
      SUBROUTINE FORCES(qD,qA,Q,fD,fA,fQ,DHB,DD,DA,aHB,aD,aA,
     &         kappa,delta,alphaD,alphaA,betaA,betaD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION eD, eA, eHB
      DOUBLE PRECISION kappa, dVA, dVD, dVQ

      eD = EXP(-aD*qD)
      eA = EXP(-aA*qA)
      eHB = EXP(-aHB*Q)

C --- Morse forces
      dVD = 2.0D0*DD*(1.0D0 - eD)*(aD*eD)
      dVA = 2.0D0*DA*(1.0D0 - eA)*(aA*eA)

C --- HB harmonic
      dVQ = 2.0D0*DHB*(1.0D0 - eHB)*(aHB*eHB)

C --- Forces
      fD = -dVD - kappa*qA - 2.0*delta*qA*qD - alphaD*Q
      fA = -dVA - kappa*qD - delta*qD*qD - alphaA*Q
      fQ = -dVQ - alphaD*qD - alphaA*qA
      fQ = fQ - betaA*qA*qA - betaD*qD*qD
      fA = fA - 2.0D0*betaA*qA*Q
      fD = fD - 2.0D0*betaD*qD*Q

      RETURN
      END
 
C --------------------------------------------------------
      SUBROUTINE ENERGY (pD,pA,pQ,qD,qA,Q,DHB,DD,DA,aHB,aD,aA,
     & mA,mD,mQ,kappa,delta,alphaD,alphaA,betaA,betaD,H,V,T)    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION mA, mD, mQ, kappa

      eD = EXP(-aD*qD)
      eA = EXP(-aA*qA)
      eHB = EXP(-aHB*Q)
      HD = 0.5D0*pD*pD/mD
      HA = 0.5D0*pA*pA/mA
      HQ = 0.5D0*pQ*pQ/mQ
      H=HD+HA+HQ
      VMD = DD*(1.0D0 - eD)**2
      VMA = DA*(1.0D0 - eA)**2
      VMHB = DHB*(1.0D0 - eHB)**2
      V=VMD+VMA+VMHB
      V=V+kappa*qA*qD+alphaD*qD*Q+alphaA*qA*Q
      V=V+delta*qD*qD*qA+betaA*qA*qA*Q+betaD*qD*qD*Q
      T=H+V

      RETURN
      END



