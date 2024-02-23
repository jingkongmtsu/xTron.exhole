      SUBROUTINE b05NDOpp1P(THRESH,VAL_P,
     $     RA,RB,GAA,GBB,TA,TB,LA,LB,UA,UB,HIRWT,
     $     F0,DRhoA,DGAA,DTA,DLA,DUA,DRhoB,DGBB,DTB,DLB,DUB)

      Implicit none

      !Input variables.

      REAL*8 THRESH,VAL_P,RA,RB,GAA,GBB,TA,TB,LA,LB,UA,UB,HIRWT

      !Return variables.
      REAL*8 F0,DRhoA,DGAA,DTA,DLA,DUA,DRhoB,DGBB,DTB,DLB,DUB
      !Internal variables
      REAL*8 ZERO,ONE,F12,TEMP
      REAL*8 FAC
      REAL*8 VA
      REAL*8 VB
      REAL*8 DFDRA
      REAL*8 DFDGAA
      REAL*8 DFDTA
      REAL*8 DFDLA
      REAL*8 DFDUA
      REAL*8 DFDRB
      REAL*8 DFDGBB
      REAL*8 DFDTB
      REAL*8 DFDLB
      REAL*8 DFDUB

      ZERO = 0.0D0
      ONE  = 1.0D0
      F12  = 0.5D0

      ! now derive the f factor
         FAC    = ZERO
         DFDRA  = ZERO
         DFDGAA = ZERO
         DFDTA  = ZERO
         DFDLA  = ZERO
         DFDUA  = ZERO
         DFDRB  = ZERO
         DFDGBB = ZERO
         DFDTB  = ZERO
         DFDLB  = ZERO
         DFDUB  = ZERO
         CALL becke05_f(THRESH,VAL_P,HIRWT,RA,GAA,TA,LA,UA,
     $RB,GBB,TB,LB,UB,
     $FAC,DFDRA,DFDGAA,DFDTA,DFDLA,DFDUA,DFDRB,DFDGBB,DFDTB,DFDLB,DFDUB)

         ! remember, here we use potential for exchange energy
         ! therefore it should be divided by rho
         VA = ZERO
         VB = ZERO
         IF (RA > THRESH) THEN
            VA = UA/RA
         END IF
         IF (RB > THRESH) THEN
            VB = UB/RB
         END IF

         ! now compose the finaly energy
         F0 = F12*FAC*(RA*VB+RB*VA)
c         write(*,*) 'usopp = ', F(i)  !jk. usopp value.
         ! derivatives
         TEMP = RA*VB+RB*VA
         DGAA = F12*DFDGAA*TEMP
         DGBB = F12*DFDGBB*TEMP
         DTA  = F12*DFDTA*TEMP
         DTB  = F12*DFDTB*TEMP
         DLA  = F12*DFDLA*TEMP
         DLB  = F12*DFDLB*TEMP

         ! derivatives for UA and UB
         DUA  = F12*DFDUA*TEMP 
         DRhoA= F12*DFDRA*TEMP + F12*FAC*VB
         IF (RA > THRESH) THEN
            DUA  = DUA + F12*FAC*(RB/RA)
            IF (RA*RA > THRESH) THEN
               DRhoA= DRhoA - F12*FAC*(RB/RA)*(UA/RA)
            END IF
         END IF
         DUB  = F12*DFDUB*TEMP 
         DRhoB= F12*DFDRB*TEMP + F12*FAC*VA
         IF (RB > THRESH) THEN
            DUB  = DUB + F12*FAC*RA/RB
            IF (RB*RB > THRESH) THEN
               DRhoB= DRhoB - F12*FAC*(RA/RB)*(UB/RB)
            END IF
         END IF
      return
      end
