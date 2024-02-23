      subroutine b05NDPar1P(NA, NB, THRESH, VAL_P,VAL_Q,
     $     RA,RB,GAA,GBB,TA,TB,LA,LB,UA,UB,HIRWT,
     $     F0,DRhoA,DGAA,DTA,DLA,DUA,DRhoB,DGBB,DTB,DLB,DUB)

         implicit none
         !Input
         INTEGER NA, NB
         real*8 THRESH,VAL_P,VAL_Q,RA,RB,GAA,GBB,TA,TB,LA,LB,UA,UB,HIRWT
         !output
         real*8 F0,DRhoA,DGAA,DTA,DLA,DUA,DRhoB,DGBB,DTB,DLB,DUB,DGAB
         !Internal
         INTEGER iSpin, NE
         real*8 A, DADRA, DADGAA, DADTA, DADLA, DADUA
         real*8    DADRB, DADGBB, DADTB, DADLB, DADUB
         real*8 M1, DM1DRA, DM1DGAA, DM1DTA, DM1DLA, DM1DUA
         real*8     DM1DRB, DM1DGBB, DM1DTB, DM1DLB, DM1DUB
         real*8 ZERO, ONE, F12

      ZERO = 0.0D0
      ONE  = 1.0D0
      F12  = 0.5D0
      NE   = NA+NB


         ! now derive the alpha part first
         iSpin  = 1
         A      = ZERO
         DADRA  = ZERO
         DADGAA = ZERO
         DADTA  = ZERO
         DADLA  = ZERO
         DADUA  = ZERO
         DADRB  = ZERO
         DADGBB = ZERO
         DADTB  = ZERO
         DADLB  = ZERO
         DADUB  = ZERO
         CALL becke05_A_sigma(iSpin,NA,NE,THRESH,
     $VAL_P,VAL_Q,HIRWT,
     $RA,GAA,TA,LA,UA,RB,GBB,TB,LB,UB,A,DADRA,DADGAA,DADTA,DADLA,
     $DADUA,DADRB,DADGBB,DADTB,DADLB,DADUB)

         ! now derive M1
         M1      = ZERO
         DM1DRA  = ZERO
         DM1DGAA = ZERO
         DM1DTA  = ZERO
         DM1DLA  = ZERO
         DM1DUA  = ZERO
         CALL becke05_M1(THRESH,RA,GAA,TA,LA,UA,HIRWT,M1,
     $DM1DRA,DM1DGAA,DM1DTA,DM1DLA,DM1DUA)

         ! now compose the finaly energy in terms of alpha
         F0 = -F12*RA*A*M1
         ! derivatives
         DRhoA= -F12*A*M1-F12*RA*DADRA*M1-F12*RA*A*DM1DRA
         DGAA = -F12*RA*DADGAA*M1-F12*RA*A*DM1DGAA
         DTA  = -F12*RA*DADTA*M1-F12*RA*A*DM1DTA
         DLA  = -F12*RA*DADLA*M1-F12*RA*A*DM1DLA
         DUA  = -F12*RA*DADUA*M1-F12*RA*A*DM1DUA
         DRhoB= -F12*RA*DADRB*M1
         DGBB = -F12*RA*DADGBB*M1
         DTB  = -F12*RA*DADTB*M1
         DLB  = -F12*RA*DADLB*M1
         DUB  = -F12*RA*DADUB*M1


         ! now let's do beta part
         iSpin  = 2
         A      = ZERO
         DADRA  = ZERO
         DADGAA = ZERO
         DADTA  = ZERO
         DADLA  = ZERO
         DADUA  = ZERO
         DADRB  = ZERO
         DADGBB = ZERO
         DADTB  = ZERO
         DADLB  = ZERO
         DADUB  = ZERO
         CALL becke05_A_sigma(iSpin,NB,NE,THRESH,
     $VAL_P,VAL_Q,HIRWT,
     $RA,GAA,TA,LA,UA,RB,GBB,TB,LB,UB,A,DADRA,DADGAA,DADTA,DADLA,
     $DADUA,DADRB,DADGBB,DADTB,DADLB,DADUB)

         ! now derive M1
         M1      = ZERO
         DM1DRB  = ZERO
         DM1DGBB = ZERO
         DM1DTB  = ZERO
         DM1DLB  = ZERO
         DM1DUB  = ZERO
         CALL becke05_M1(THRESH,RB,GBB,TB,LB,UB,HIRWT,M1,
     $DM1DRB,DM1DGBB,DM1DTB,DM1DLB,DM1DUB)

         ! now compose the finaly energy in terms of alpha
         F0 = -F12*RB*A*M1 + F0
         ! derivatives
         DRhoA= -F12*RB*DADRA*M1 + DRhoA
         DGAA = -F12*RB*DADGAA*M1 + DGAA
         DTA  = -F12*RB*DADTA*M1 + DTA
         DLA  = -F12*RB*DADLA*M1 + DLA
         DUA  = -F12*RB*DADUA*M1 + DUA
         DRhoB= -F12*A*M1-F12*RB*DADRB*M1-F12*RB*A*DM1DRB + DRhoB
         DGBB = -F12*RB*DADGBB*M1-F12*RB*A*DM1DGBB + DGBB
         DTB  = -F12*RB*DADTB*M1-F12*RB*A*DM1DTB + DTB
         DLB  = -F12*RB*DADLB*M1-F12*RB*A*DM1DLB + DLB
         DUB  = -F12*RB*DADUB*M1-F12*RB*A*DM1DUB + DUB

      return
      end

      subroutine b05NDPar1Pxx(NA, NB, THRESH, VAL_P,VAL_Q,
     $     RA,RB,GAA,GBB,TA,TB,LA,LB,UA,UB,HIRWT,
     $     F0,DRhoA,DGAA,DTA,DLA,DUA,DRhoB,DGBB,DTB,DLB,DUB)

         implicit none
         !Input
         INTEGER NA, NB
         real*8 THRESH,VAL_P,VAL_Q,RA,RB,GAA,GBB,TA,TB,LA,LB,UA,UB,HIRWT
         !output
         real*8 F0,DRhoA,DGAA,DTA,DLA,DUA,DRhoB,DGBB,DTB,DLB,DUB,DGAB
         !Internal
         INTEGER iSpin, NE
         real*8 A, DADRA, DADGAA, DADTA, DADLA, DADUA
         real*8    DADRB, DADGBB, DADTB, DADLB, DADUB
         real*8 M1, DM1DRA, DM1DGAA, DM1DTA, DM1DLA, DM1DUA
         real*8     DM1DRB, DM1DGBB, DM1DTB, DM1DLB, DM1DUB
         real*8 ZERO, ONE, F12

      ZERO = 0.0D0
      ONE  = 1.0D0
      F12  = 0.5D0
      NE   = NA+NB


         ! now derive the alpha part first
         iSpin  = 1
         A      = ZERO
         DADRA  = ZERO
         DADGAA = ZERO
         DADTA  = ZERO
         DADLA  = ZERO
         DADUA  = ZERO
         DADRB  = ZERO
         DADGBB = ZERO
         DADTB  = ZERO
         DADLB  = ZERO
         DADUB  = ZERO
         CALL becke05_A_sigma(iSpin,NA,NE,THRESH,
     $VAL_P,VAL_Q,HIRWT,
     $RA,GAA,TA,LA,UA,RB,GBB,TB,LB,UB,A,DADRA,DADGAA,DADTA,DADLA,
     $DADUA,DADRB,DADGBB,DADTB,DADLB,DADUB)

         ! now derive M1
         M1      = ZERO
         DM1DRA  = ZERO
         DM1DGAA = ZERO
         DM1DTA  = ZERO
         DM1DLA  = ZERO
         DM1DUA  = ZERO
         CALL becke05_M1(THRESH,RA,GAA,TA,LA,UA,HIRWT,M1,
     $DM1DRA,DM1DGAA,DM1DTA,DM1DLA,DM1DUA)

         ! now compose the finaly energy in terms of alpha
         F0 = -F12*RA*A*M1

         ! derivatives
         DRhoA= -F12*A*M1-F12*RA*DADRA*M1-F12*RA*A*DM1DRA
         DGAA = -F12*RA*DADGAA*M1-F12*RA*A*DM1DGAA
         DTA  = -F12*RA*DADTA*M1-F12*RA*A*DM1DTA
         DLA  = -F12*RA*DADLA*M1-F12*RA*A*DM1DLA
         DUA  = -F12*RA*DADUA*M1-F12*RA*A*DM1DUA
         DRhoB= -F12*RA*DADRB*M1
         DGBB = -F12*RA*DADGBB*M1
         DTB  = -F12*RA*DADTB*M1
         DLB  = -F12*RA*DADLB*M1
         DUB  = -F12*RA*DADUB*M1
      return
      end
