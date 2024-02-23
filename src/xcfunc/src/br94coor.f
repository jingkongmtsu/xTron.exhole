
      SUBROUTINE br94coor_par(INFOR,NDEN,NA,NB,NG,THRESH,RhoA,RhoB,
     $DRA,DRB,TauA,TauB,LapA,LapB,F,D1F)
c
c    ******************************************************************
c    *   This is BR94 correlation functional which defined in:        *
c    *   A. D. Becke, Int. J. Quant. Chem (symp), 28, 625 (1994)      *
c    *   this is for the parallel spin component                      *         
c    *                                                                *
c    *   input:                                                       *
c    *   INFOR:      information array about variable position        *
c    *   NDEN:       number of density(close shell: 1, open shell: 2) *
c    *   NA,NB:      number of alpha/beta electrons                   *
c    *   NG:         number grid points                               *         
c    *   THRESH:     threshold to determine the small rho etc.        *         
c    *   RhoA,RhoB:  spin density                                     *
c    *   DRA,DRB:    gradient of spin density                         *
c    *   TauA,TauB:  kinetic density                                  *
c    *   LapA,LapB:  laplacian of density                             *
c    *                                                                *
c    *   output:                                                      *
c    *   F:        functional value                                   *
c    *   D1F:      functional derivatives                             *         
c    *                                                                *
c    *   fenglai liu and E. Proynov                                   *
c    *                                                                *
c    ******************************************************************
c
      IMPLICIT NONE
      INTEGER NDEN,NG,I
      INTEGER NA,NB
      INTEGER INFOR(*)
      REAL*8 THRESH
      REAL*8 RhoA(NG),RhoB(NG),DRA(NG,3),DRB(NG,3)
      REAL*8 TauA(NG),TauB(NG),LapA(NG),LapB(NG)
      REAL*8 F(NG),D1F(NG,*)
      REAL*8 RA,RB
      REAL*8 GAA,GBB 
      REAL*8 TA,TB
      REAL*8 LA,LB
#include "fderiv1.inc"
      INTEGER  D1VARS(N_FUNC_DERIV_1)
      INTEGER  ID_RA_POS, ID_GAA_POS, ID_GAB_POS, ID_TA_POS, ID_LA_POS
      INTEGER  ID_RB_POS, ID_GBB_POS, ID_TB_POS, ID_LB_POS

      ! now let's define the parameter
      REAL*8 CAA,CBB,CAB
      REAL*8 CAA2,CBB2,CAB2

      ! variables for U and it's derivatives
      REAL*8 UA,UA_RA,UA_GAA,UA_TA,UA_LA    
      REAL*8 UB,UB_RB,UB_GBB,UB_TB,UB_LB    

      ! variables for ZAA,ZAB,ZBB and it's derivatives
      REAL*8 ZAA,ZAA_RA,ZAA_GAA,ZAA_TA,ZAA_LA    
      REAL*8 ZBB,ZBB_RB,ZBB_GBB,ZBB_TB,ZBB_LB    
      REAL*8 ZAA2,ZAA3,ZAA4
      REAL*8 ZBB2,ZBB3,ZBB4
      REAL*8 ZAA_U,ZAB_U,ZBB_U
      REAL*8 FZ0,FZ1,FZ,FZ_Z
      REAL*8 ONEDU

      ! D and it's derivatives
      REAL*8 D,DDDR,DDDG,DDDT

      ! variables for EC and it's derivatives
      REAL*8 ECAA,ECAA_RA,ECAA_GAA,ECAA_TA,ECAA_LA  
      REAL*8 ECBB,ECBB_RB,ECBB_GBB,ECBB_TB,ECBB_LB  
      REAL*8 ECAB,ECAB_RA,ECAB_GAA,ECAB_TA,ECAB_LA  
      REAL*8      ECAB_RB,ECAB_GBB,ECAB_TB,ECAB_LB  

C       
C      the tau used in this program is for "big tau", that 
C      is to say, without factor of 1/2
C

       ! constants defined in the paper
       CAA   = -0.88D0
       CBB   = -0.88D0
       CAB   = -0.63D0
       CAA2  = -0.01D0
       CBB2  = -0.01D0
       CAB2  = -0.8D0

       ! firstly initilize variable position information
       CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)
      
       ! loop over NG
      DO I = 1,NG

         ! variables
         RA  = RhoA(i)
         RB  = RhoB(i)
         GAA = DRA(i,1)*DRA(i,1) + DRA(i,2)*DRA(i,2)
     &+ DRA(i,3)*DRA(i,3)
         GBB = DRB(i,1)*DRB(i,1) + DRB(i,2)*DRB(i,2)
     &+ DRB(i,3)*DRB(i,3)
         TA  = TauA(i)
         TB  = TauB(i)
         LA  = LapA(i)
         LB  = LapB(i)
c         write(6,*)"our VAR",RA,RB,GAA,GBB,TA,TB,LA,LB
C         write(6,*)"D1RA(i,1)", DRA(i,1)
C         write(6,*)"D1RA(i,2)", DRA(i,2)
C         write(6,*)"D1RA(i,3)", DRA(i,3)
         
         ! EC_{alpha,alpha}
         UA       = 0.0D0
         UA_RA    = 0.0D0
         UA_GAA   = 0.0D0
         UA_TA    = 0.0D0
         UA_LA    = 0.0D0
         ECAA     = 0.0D0
         ECAA_RA  = 0.0D0
         ECAA_GAA = 0.0D0
         ECAA_TA  = 0.0D0
         ECAA_LA  = 0.0D0
         IF (NA > 1) THEN
            IF (RA > THRESH) THEN
   
               ! BR89 exchange hole
               CALL BR89HOLE(THRESH,RA,GAA,TA,LA,UA,UA_RA,UA_GAA,
     $                       UA_TA,UA_LA)
   
               ! now it's ZAA
               ZAA       = 0.0D0
               ZAA_RA    = 0.0D0
               ZAA_GAA   = 0.0D0
               ZAA_TA    = 0.0D0
               ZAA_LA    = 0.0D0
c               write(6,*)"UA",UA
               IF (ABS(UA) > THRESH) THEN
                  ONEDU   = 1.0D0/UA
                  ZAA     = CAA*2.0D0*ONEDU
                  ZAA_U   = -2.0D0*CAA/(UA*UA)
                  ZAA_RA  = ZAA_U*UA_RA
                  ZAA_GAA = ZAA_U*UA_GAA
                  ZAA_TA  = ZAA_U*UA_TA
                  ZAA_LA  = ZAA_U*UA_LA
               END IF
c               write(6,*)"ZAA",ZAA
c               write(6,*)"ZAA_RA",ZAA_RA
   
               ! now it's function of ZAA
               FZ   = 0.0D0
               FZ_Z = 0.0D0
               IF (ABS(ZAA) > THRESH) THEN
                  ZAA2 = ZAA*ZAA
                  ZAA3 = ZAA2*ZAA
                  ZAA4 = ZAA3*ZAA
                  FZ0  = DLOG(1.0D0+ZAA/2.0D0)
                  FZ1  = 1.0D0-(2.0D0/ZAA)*FZ0
                  FZ   = ZAA4*FZ1
                  FZ_Z = 4.0D0*ZAA3*FZ1 + ZAA4*(2.0D0*FZ0/ZAA2 -
     $                   1.0D0/(ZAA*(1.0D0+ZAA/2.0D0)))
   
               END IF
c               write(6,*)"FZ for ZAA",FZ
   
               ! D and it's derivatives
               D    = 0.0D0
               DDDR = 0.0D0
               DDDG = 0.0D0
               DDDT = 0.0D0
               IF (RA > THRESH) THEN
                  D    = TA - 0.25D0*GAA/RA
                  DDDR = 0.25D0*GAA/(RA*RA)
                  DDDG = -0.25D0/RA
                  DDDT = 1.0D0
               END IF
C               write(6,*)"D",D
c               write(6,*)"D",D
c               write(6,*)"D_RA",DDDR
   
               ! now let's assemble all of pices together
               ECAA      = CAA2*RA*D*FZ
               ECAA_RA   = CAA2*(D*FZ+RA*DDDR*FZ+RA*D*FZ_Z*ZAA_RA)
               ECAA_GAA  = CAA2*RA*(DDDG*FZ+D*FZ_Z*ZAA_GAA)
               ECAA_TA   = CAA2*RA*(DDDT*FZ+D*FZ_Z*ZAA_TA)
               ECAA_LA   = CAA2*RA*D*FZ_Z*ZAA_LA
c               write(6,*)"ECAA",ECAA
            END IF
         END IF

         ! EC_{beta,beta}
         UB       = 0.0D0
         UB_RB    = 0.0D0
         UB_GBB   = 0.0D0
         UB_TB    = 0.0D0
         UB_LB    = 0.0D0
         ECBB     = 0.0D0
         ECBB_RB  = 0.0D0
         ECBB_GBB = 0.0D0
         ECBB_TB  = 0.0D0
         ECBB_LB  = 0.0D0
         IF (NB > 1) THEN
            IF (NDEN == 1) THEN
               UB       = UA
               UB_RB    = UA_RA
               UB_GBB   = UA_GAA
               UB_TB    = UA_TA
               UB_LB    = UA_LA
               ECBB     = ECAA
               ECBB_RB  = ECAA_RA
               ECBB_GBB = ECAA_GAA
               ECBB_TB  = ECAA_TA
               ECBB_LB  = ECAA_LA
            ELSE IF (RB > THRESH) THEN
   
               ! BR89 exchange hole
               CALL BR89HOLE(THRESH,RB,GBB,TB,LB,UB,UB_RB,UB_GBB,
     $                       UB_TB,UB_LB)
   
               ! now it's ZBB
               ZBB      = 0.0D0
               ZBB_RB   = 0.0D0
               ZBB_GBB  = 0.0D0
               ZBB_TB   = 0.0D0
               ZBB_LB   = 0.0D0
c               write(6,*)"UB",UB
               IF (ABS(UB) > THRESH) THEN
                  ONEDU   = 1.0D0/UB
                  ZBB     = CBB*2.0D0*ONEDU
                  ZBB_U   = -2.0D0*CBB/(UB*UB)
                  ZBB_RB  = ZBB_U*UB_RB
                  ZBB_GBB = ZBB_U*UB_GBB
                  ZBB_TB  = ZBB_U*UB_TB
                  ZBB_LB  = ZBB_U*UB_LB
               END IF
c               write(6,*)"ZBB",ZBB
   
               ! now it's function of ZBB
               FZ   = 0.0D0
               FZ_Z = 0.0D0
               IF (ABS(ZBB) > THRESH) THEN
                  ZBB2 = ZBB*ZBB
                  ZBB3 = ZBB2*ZBB
                  ZBB4 = ZBB3*ZBB
                  FZ0  = DLOG(1.0D0+ZBB/2.0D0)
                  FZ1  = 1.0D0-(2.0D0/ZBB)*FZ0
                  FZ   = ZBB4*FZ1
                  FZ_Z = 4.0D0*ZBB3*FZ1 + ZBB4*(2.0D0*FZ0/ZBB2 -
     $                   1.0D0/(ZBB*(1.0D0+ZBB/2.0D0)))
   
               END IF
c               write(6,*)"FZ for ZBB",FZ
   
               ! D and it's derivatives
               D    = TB - 0.25D0*GBB/RB
               DDDR = 0.25D0*GBB/(RB*RB)
               DDDG = -0.25D0/RB
               DDDT = 1.0D0
c              write(6,*)"D",D
   
               ! now let's assemble all of pices together
               ECBB      = CBB2*RB*D*FZ
               ECBB_RB   = CBB2*(D*FZ+RB*DDDR*FZ+RB*D*FZ_Z*ZBB_RB)
               ECBB_GBB  = CBB2*RB*(DDDG*FZ+D*FZ_Z*ZBB_GBB)
               ECBB_TB   = CBB2*RB*(DDDT*FZ+D*FZ_Z*ZBB_TB)
               ECBB_LB   = CBB2*RB*D*FZ_Z*ZBB_LB
c               write(6,*)"ECBB",ECBB
            END IF
         END IF

         ! collect energy
         F(i) = F(i) + ECAA + ECBB

         ! collect all of terms for alpha derivatives
         ID_RA_POS=D1VARS(ID_RA)
         D1F(i, ID_RA_POS)  = D1F(i, ID_RA_POS) + ECAA_RA 
         ID_GAA_POS=D1VARS(ID_GAA)
         D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) + ECAA_GAA 
         ID_GAB_POS=D1VARS(ID_GAB)
         D1F(i, ID_GAB_POS) = 0.0D0
         ID_TA_POS=D1VARS(ID_TA)
         D1F(i, ID_TA_POS)  = D1F(i, ID_TA_POS) + ECAA_TA 
         ID_LA_POS=D1VARS(ID_LA)
         D1F(i, ID_LA_POS)  = D1F(i, ID_LA_POS) + ECAA_LA 
         
         ! collect all of terms for beta derivatives
         ID_RB_POS=D1VARS(ID_RB)
         D1F(i, ID_RB_POS)  = D1F(i, ID_RB_POS) + ECBB_RB 
         ID_GBB_POS=D1VARS(ID_GBB)
         D1F(i, ID_GBB_POS) = D1F(i, ID_GBB_POS) + ECBB_GBB 
         ID_TB_POS=D1VARS(ID_TB)
         D1F(i, ID_TB_POS)  = D1F(i, ID_TB_POS) + ECBB_TB 
         ID_LB_POS=D1VARS(ID_LB)
         D1F(i, ID_LB_POS)  = D1F(i, ID_LB_POS) + ECBB_LB 

      END DO

      RETURN 
      END

#if 0
      SUBROUTINE br94coor_op(INFOR,NDEN,NA,NB,NG,THRESH,RhoA,RhoB,
     $DRA,DRB,TauA,TauB,LapA,LapB,F,D1F)
c
c    ******************************************************************
c    *   This is BR94 correlation functional which defined in:        *
c    *   A. D. Becke, Int. J. Quant. Chem (symp), 28, 625 (1994)      *
c    *   this is for the opposite spin component of the functional    *         
c    *                                                                *
c    *   input:                                                       *
c    *   INFOR:      information array about variable position        *
c    *   NDEN:       number of density(close shell: 1, open shell: 2) *
c    *   NA,NB:      number of alpha/beta electrons                   *
c    *   NG:         number grid points                               *         
c    *   THRESH:     threshold to determine the small rho etc.        *         
c    *   RhoA,RhoB:  spin density                                     *
c    *   DRA,DRB:    gradient of spin density                         *
c    *   TauA,TauB:  kinetic density                                  *
c    *   LapA,LapB:  laplacian of density                             *
c    *                                                                *
c    *   output:                                                      *
c    *   F:        functional value                                   *
c    *   D1F:      functional derivatives                             *         
c    *                                                                *
c    *   fenglai liu and E. Proynov                                   *
c    *                                                                *
c    ******************************************************************
c
      IMPLICIT NONE
      INTEGER NDEN,NG,I
      INTEGER NA,NB
      INTEGER INFOR(*)
      REAL*8 THRESH
      REAL*8 RhoA(NG),RhoB(NG),DRA(NG,3),DRB(NG,3)
      REAL*8 TauA(NG),TauB(NG),LapA(NG),LapB(NG)
      REAL*8 F(NG),D1F(NG,*)
      REAL*8 RA,RB
      REAL*8 GAA,GBB 
      REAL*8 TA,TB
      REAL*8 LA,LB
#include "fderiv1.inc"
      INTEGER  D1VARS(N_FUNC_DERIV_1)
      INTEGER  ID_RA_POS, ID_GAA_POS, ID_GAB_POS, ID_TA_POS, ID_LA_POS
      INTEGER  ID_RB_POS, ID_GBB_POS, ID_TB_POS, ID_LB_POS

      ! now let's define the parameter
      REAL*8 CAA,CBB,CAB
      REAL*8 CAA2,CBB2,CAB2

      ! variables for U and it's derivatives
      REAL*8 UA,UA_RA,UA_GAA,UA_TA,UA_LA    
      REAL*8 UB,UB_RB,UB_GBB,UB_TB,UB_LB    

      ! variables for ZAA,ZAB,ZBB and it's derivatives
      REAL*8 ZAB,ZAB_RA,ZAB_GAA,ZAB_TA,ZAB_LA    
      REAL*8     ZAB_RB,ZAB_GBB,ZAB_TB,ZAB_LB    
      REAL*8 ZAA_U,ZAB_U,ZBB_U
      REAL*8 FZ0,FZ1,FZ,FZ_Z
      REAL*8 ONEDU

      ! D and it's derivatives
      REAL*8 D,DDDR,DDDG,DDDT

      ! variables for EC and it's derivatives
      REAL*8 ECAB,ECAB_RA,ECAB_GAA,ECAB_TA,ECAB_LA  
      REAL*8      ECAB_RB,ECAB_GBB,ECAB_TB,ECAB_LB  

C       
C      the tau used in this program is for "big tau", that 
C      is to say, without factor of 1/2
C

       ! constants defined in the paper
       CAA   = -0.88D0
       CBB   = -0.88D0
       CAB   = -0.63D0
       CAA2  = -0.01D0
       CBB2  = -0.01D0
       CAB2  = -0.8D0

       ! firstly initilize variable position information
       CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)
      
       ! loop over NG
      DO I = 1,NG

         ! variables
         RA  = RhoA(i)
         RB  = RhoB(i)
         GAA = DRA(i,1)*DRA(i,1) + DRA(i,2)*DRA(i,2)
     &+ DRA(i,3)*DRA(i,3)
         GBB = DRB(i,1)*DRB(i,1) + DRB(i,2)*DRB(i,2)
     &+ DRB(i,3)*DRB(i,3)
         TA  = TauA(i)
         TB  = TauB(i)
         LA  = LapA(i)
         LB  = LapB(i)
c         write(6,*)"our VAR",RA,RB,GAA,GBB,TA,TB,LA,LB
C         write(6,*)"D1RA(i,1)", DRA(i,1)
C         write(6,*)"D1RA(i,2)", DRA(i,2)
C         write(6,*)"D1RA(i,3)", DRA(i,3)
         
         ! EC_{alpha,alpha}
         UA       = 0.0D0
         UA_RA    = 0.0D0
         UA_GAA   = 0.0D0
         UA_TA    = 0.0D0
         UA_LA    = 0.0D0
         IF (NA > 0) THEN
            IF (RA > THRESH) THEN
   
               ! BR89 exchange hole
               CALL BR89HOLE(THRESH,RA,GAA,TA,LA,UA,UA_RA,UA_GAA,
     $                       UA_TA,UA_LA)
c               write(6,*)"UA in ZAB", UA
            END IF
         END IF

         ! EC_{beta,beta}
         UB       = 0.0D0
         UB_RB    = 0.0D0
         UB_GBB   = 0.0D0
         UB_TB    = 0.0D0
         UB_LB    = 0.0D0
         IF (NB > 0) THEN
            IF (NDEN == 1) THEN
               UB       = UA
               UB_RB    = UA_RA
               UB_GBB   = UA_GAA
               UB_TB    = UA_TA
               UB_LB    = UA_LA
            ELSE IF (RB > THRESH) THEN
   
               ! BR89 exchange hole
               CALL BR89HOLE(THRESH,RB,GBB,TB,LB,UB,UB_RB,UB_GBB,
     $                       UB_TB,UB_LB)
c               write(6,*)"UB in ZAB", UB
            END IF
         END IF

         ! now it's ECAB
         ECAB     = 0.0D0
         ECAB_RA  = 0.0D0
         ECAB_GAA = 0.0D0
         ECAB_TA  = 0.0D0
         ECAB_LA  = 0.0D0
         ECAB_RB  = 0.0D0
         ECAB_GBB = 0.0D0
         ECAB_TB  = 0.0D0
         ECAB_LB  = 0.0D0
         IF (RA > THRESH .and. RB > THRESH .and. 
     $       NA > 0 .and. NB > 0) THEN

            ! now it's ZAB
            ZAB      = 0.0D0
            ZAB_RA   = 0.0D0
            ZAB_GAA  = 0.0D0
            ZAB_TA   = 0.0D0
            ZAB_LA   = 0.0D0
            ZAB_RB   = 0.0D0
            ZAB_GBB  = 0.0D0
            ZAB_TB   = 0.0D0
            ZAB_LB   = 0.0D0
            IF (ABS(UA) > THRESH) THEN
               ONEDU   = 1.0D0/UA
               ZAB     = CAB*ONEDU
               ZAB_U   = -CAB/(UA*UA)
               ZAB_RA  = ZAB_U*UA_RA
               ZAB_GAA = ZAB_U*UA_GAA
               ZAB_TA  = ZAB_U*UA_TA
               ZAB_LA  = ZAB_U*UA_LA
            END IF
            IF (NDEN == 1) THEN
               ZAB     = 2.0D0*ZAB 
               ZAB_RB  = ZAB_RA
               ZAB_GBB = ZAB_GAA
               ZAB_TB  = ZAB_TA
               ZAB_LB  = ZAB_LA
            ELSE IF (ABS(UB) > THRESH) THEN
               ONEDU   = 1.0D0/UB
               ZAB     = ZAB + CAB*ONEDU
               ZAB_U   = -CAB/(UB*UB)
               ZAB_RB  = ZAB_U*UB_RB
               ZAB_GBB = ZAB_U*UB_GBB
               ZAB_TB  = ZAB_U*UB_TB
               ZAB_LB  = ZAB_U*UB_LB
            END IF
c            write(6,*)"ZAB",ZAB

            ! now it's function of ZAB
            FZ   = 0.0D0
            FZ_Z = 0.0D0
            IF (ABS(ZAB) > THRESH) THEN
               FZ0  = DLOG(1.0D0+ZAB)/ZAB
               FZ   = ZAB*ZAB*(1.0D0-FZ0)
               FZ_Z = 2.0D0*ZAB*(1.0D0-FZ0) + ZAB*ZAB*(FZ0/ZAB -
     $                1.0D0/((1.0D0+ZAB)*ZAB))

            END IF
c            write(6,*)"FZ for ZAB",FZ

            ! now let's assemble all of pices together
            ECAB      = CAB2*RA*RB*FZ
            ECAB_RA   = CAB2*RB*(FZ+RA*FZ_Z*ZAB_RA)
            ECAB_RB   = CAB2*RA*(FZ+RB*FZ_Z*ZAB_RB)
            ECAB_GAA  = CAB2*RA*RB*FZ_Z*ZAB_GAA
            ECAB_GBB  = CAB2*RA*RB*FZ_Z*ZAB_GBB
            ECAB_TA   = CAB2*RA*RB*FZ_Z*ZAB_TA
            ECAB_TB   = CAB2*RA*RB*FZ_Z*ZAB_TB
            ECAB_LA   = CAB2*RA*RB*FZ_Z*ZAB_LA
            ECAB_LB   = CAB2*RA*RB*FZ_Z*ZAB_LB
         END IF

         ! collect energy
         F(i) = F(i) + ECAB 

         ! collect all of terms for alpha derivatives
         ID_RA_POS=D1VARS(ID_RA)
         D1F(i, ID_RA_POS)  = D1F(i, ID_RA_POS) + ECAB_RA
         ID_GAA_POS=D1VARS(ID_GAA)
         D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) + ECAB_GAA
         ID_GAB_POS=D1VARS(ID_GAB)
         D1F(i, ID_GAB_POS) = 0.0D0
         ID_TA_POS=D1VARS(ID_TA)
         D1F(i, ID_TA_POS)  = D1F(i, ID_TA_POS) + ECAB_TA
         ID_LA_POS=D1VARS(ID_LA)
         D1F(i, ID_LA_POS)  = D1F(i, ID_LA_POS) + ECAB_LA
         
         ! collect all of terms for beta derivatives
         ID_RB_POS=D1VARS(ID_RB)
         D1F(i, ID_RB_POS)  = D1F(i, ID_RB_POS) + ECAB_RB
         ID_GBB_POS=D1VARS(ID_GBB)
         D1F(i, ID_GBB_POS) = D1F(i, ID_GBB_POS) + ECAB_GBB
         ID_TB_POS=D1VARS(ID_TB)
         D1F(i, ID_TB_POS)  = D1F(i, ID_TB_POS) + ECAB_TB
         ID_LB_POS=D1VARS(ID_LB)
         D1F(i, ID_LB_POS)  = D1F(i, ID_LB_POS) + ECAB_LB

      END DO

      RETURN 
      END

#else

      SUBROUTINE br94coor_op(INFOR,NDEN,NA,NB,NG,THRESH,RhoA,RhoB,
     $DRA,DRB,TauA,TauB,LapA,LapB,F,D1F)
c
c    ******************************************************************
c    *   This is BR94 correlation functional which defined in:        *
c    *   A. D. Becke, Int. J. Quant. Chem (symp), 28, 625 (1994)      *
c    *   this is for the opposite spin component of the functional    *
c    *   Note: this is the lamda-averaged formula.
c    *                                                                *
c    *   input:                                                       *
c    *   INFOR:      information array about variable position        *
c    *   NDEN:       number of density(close shell: 1, open shell: 2) *
c    *   NA,NB:      number of alpha/beta electrons                   *
c    *   NG:         number grid points                               *
c    *   THRESH:     threshold to determine the small rho etc.        *
c    *   RhoA,RhoB:  spin density                                     *
c    *   DRA,DRB:    gradient of spin density                         *
c    *   TauA,TauB:  kinetic density                                  *
c    *   LapA,LapB:  laplacian of density                             *
c    *                                                                *
c    *   output:                                                      *
c    *   F:        functional value                                   *
c    *   D1F:      functional derivatives                             *
c    *                                                                *
c    *   fenglai liu and E. Proynov                                   *
c    *                                                                *
c    ******************************************************************
c
      IMPLICIT NONE
      INTEGER NDEN,NG,I
      INTEGER NA,NB
      INTEGER INFOR(*)
      REAL*8 THRESH
      REAL*8 RhoA(NG),RhoB(NG),DRA(NG,3),DRB(NG,3)
      REAL*8 TauA(NG),TauB(NG),LapA(NG),LapB(NG)
      REAL*8 F(NG),D1F(NG,*)
      REAL*8 RA,RB
      REAL*8 GAA,GBB
      REAL*8 TA,TB
      REAL*8 LA,LB
#include "fderiv1.inc"
      INTEGER  D1VARS(N_FUNC_DERIV_1)
      INTEGER  ID_RA_POS, ID_GAA_POS, ID_GAB_POS, ID_TA_POS, ID_LA_POS
      INTEGER  ID_RB_POS, ID_GBB_POS, ID_TB_POS, ID_LB_POS

      ! variables for EC and it's derivatives
      REAL*8 ECAB,ECAB_RA,ECAB_GAA,ECAB_TA,ECAB_LA
      REAL*8      ECAB_RB,ECAB_GBB,ECAB_TB,ECAB_LB

       ! firstly initilize variable position information

       CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)

       ! loop over NG
      DO I = 1,NG

         ! variables
         RA  = RhoA(i)
         RB  = RhoB(i)
         GAA = DRA(i,1)*DRA(i,1) + DRA(i,2)*DRA(i,2)
     &+ DRA(i,3)*DRA(i,3)
         GBB = DRB(i,1)*DRB(i,1) + DRB(i,2)*DRB(i,2)
     &+ DRB(i,3)*DRB(i,3)
         TA  = TauA(i)
         TB  = TauB(i)
         LA  = LapA(i)
         LB  = LapB(i)

         call br94CorrOpp1P(NDEN,NA,NB,THRESH,
     $                RA,RB,GAA,GBB,TA,TB,LA,LB,
     $                ECAB,ECAB_RA,ECAB_GAA,ECAB_TA,ECAB_LA,
     $                ECAB_RB,ECAB_GBB,ECAB_TB,ECAB_LB)

         ! collect energy
         F(i) = F(i) + ECAB

         ID_RA_POS=D1VARS(ID_RA)
         D1F(i, ID_RA_POS)  = D1F(i, ID_RA_POS) + ECAB_RA
         ID_GAA_POS=D1VARS(ID_GAA)
         D1F(i, ID_GAA_POS) = D1F(i, ID_GAA_POS) + ECAB_GAA
         ID_GAB_POS=D1VARS(ID_GAB)
         D1F(i, ID_GAB_POS) = 0.0D0
         ID_TA_POS=D1VARS(ID_TA)
         D1F(i, ID_TA_POS)  = D1F(i, ID_TA_POS) + ECAB_TA
         ID_LA_POS=D1VARS(ID_LA)
         D1F(i, ID_LA_POS)  = D1F(i, ID_LA_POS) + ECAB_LA

         ! collect all of terms for beta derivatives
         ID_RB_POS=D1VARS(ID_RB)
         D1F(i, ID_RB_POS)  = D1F(i, ID_RB_POS) + ECAB_RB
         ID_GBB_POS=D1VARS(ID_GBB)
         D1F(i, ID_GBB_POS) = D1F(i, ID_GBB_POS) + ECAB_GBB
         ID_TB_POS=D1VARS(ID_TB)
         D1F(i, ID_TB_POS)  = D1F(i, ID_TB_POS) + ECAB_TB
         ID_LB_POS=D1VARS(ID_LB)
         D1F(i, ID_LB_POS)  = D1F(i, ID_LB_POS) + ECAB_LB

      END DO

      RETURN
      END
#endif

      subroutine br94CorrOpp1P(NDEN,NA,NB,THRESH,
     $                RA,RB,GAA,GBB,TA,TB,LA,LB,
     $                ECAB,ECAB_RA,ECAB_GAA,ECAB_TA,ECAB_LA,
     $                     ECAB_RB,ECAB_GBB,ECAB_TB,ECAB_LB)
      implicit none
      INTEGER NDEN,NA,NB
      REAL*8 THRESH,RA,RB,GAA,GBB,TA,TB,LA,LB

      REAL*8 ECAB,ECAB_RA,ECAB_GAA,ECAB_TA,ECAB_LA
      REAL*8      ECAB_RB,ECAB_GBB,ECAB_TB,ECAB_LB

      ! now let's define the parameter
      REAL*8 CAB2

      ! variables for U and it's derivatives
      REAL*8 UA,UA_RA,UA_GAA,UA_TA,UA_LA    
      REAL*8 UB,UB_RB,UB_GBB,UB_TB,UB_LB    

      ! variables for ZAA,ZAB,ZBB and it's derivatives
      REAL*8 ZAB,ZAB_RA,ZAB_GAA,ZAB_TA,ZAB_LA
      REAL*8     ZAB_RB,ZAB_GBB,ZAB_TB,ZAB_LB

      !For the call of br94fzab
      REAL*8 FZ,FZ_UA,FZ_UB

       ! constants defined in the paper
       CAB2  = -0.8D0

         ! EC_{alpha,alpha}
         UA       = 0.0D0
         UA_RA    = 0.0D0
         UA_GAA   = 0.0D0
         UA_TA    = 0.0D0
         UA_LA    = 0.0D0
         IF (NA > 0) THEN
            IF (RA > THRESH) THEN
   
               ! BR89 exchange hole
               CALL BR89HOLE(THRESH,RA,GAA,TA,LA,UA,UA_RA,UA_GAA,
     $                       UA_TA,UA_LA)
            END IF
         END IF

         ! EC_{beta,beta}
         UB       = 0.0D0
         UB_RB    = 0.0D0
         UB_GBB   = 0.0D0
         UB_TB    = 0.0D0
         UB_LB    = 0.0D0
         IF (NB > 0) THEN
            IF (NDEN == 1) THEN
               UB       = UA
               UB_RB    = UA_RA
               UB_GBB   = UA_GAA
               UB_TB    = UA_TA
               UB_LB    = UA_LA
            ELSE IF (RB > THRESH) THEN
   
               ! BR89 exchange hole
               CALL BR89HOLE(THRESH,RB,GBB,TB,LB,UB,UB_RB,UB_GBB,
     $                       UB_TB,UB_LB)
c               write(6,*)"UB in ZAB", UB
            END IF
         END IF

         ! now it's ECAB
         ECAB     = 0.0D0
         ECAB_RA  = 0.0D0
         ECAB_GAA = 0.0D0
         ECAB_TA  = 0.0D0
         ECAB_LA  = 0.0D0
         ECAB_RB  = 0.0D0
         ECAB_GBB = 0.0D0
         ECAB_TB  = 0.0D0
         ECAB_LB  = 0.0D0
         IF (RA > THRESH .and. RB > THRESH .and. 
     $       NA > 0 .and. NB > 0) THEN

            call br94CorrOppFZab(FZ, FZ_UA, FZ_UB, NDEN, THRESH, 
     $                    UA, UB)

            ! now let's assemble all of pices together
            ECAB      = CAB2*RA*RB*FZ
            ECAB_RA   = CAB2*RB*(FZ+RA*FZ_UA*UA_RA)
            ECAB_GAA  = CAB2*RA*RB*FZ_UA*UA_GAA
            ECAB_TA   = CAB2*RA*RB*FZ_UA*UA_TA
            ECAB_LA   = CAB2*RA*RB*FZ_UA*UA_LA
            ECAB_RB   = CAB2*RA*(FZ+RB*FZ_UB*UB_RB)
            ECAB_GBB  = CAB2*RA*RB*FZ_UB*UB_GBB
            ECAB_TB   = CAB2*RA*RB*FZ_UB*UB_TB
            ECAB_LB   = CAB2*RA*RB*FZ_UB*UB_LB
         END IF
      return
      end

      subroutine br94CorrOppFZab(FZ, FZ_UA, FZ_UB, NDEN, THRESH, 
     $                    UA, UB) 
      implicit none
      !Output
      REAL*8 FZ,FZ_UA,FZ_UB
      !Input
      INTEGER NDEN
      REAL*8 THRESH, UA, UB
      !Internal
      REAL*8 FZ0,FZ_Z,ZAB,ZAB_UA,ZAB_UB,ONEDU
      !Constant
      REAL*8 CAB

      CAB   = -0.63D0

      ZAB = 0D0
      
            IF (ABS(UA) > THRESH) THEN
               ONEDU   = 1.0D0/UA
               ZAB     = CAB*ONEDU
               ZAB_UA  = -CAB/(UA*UA)
            END IF
            IF (NDEN == 1) THEN
               ZAB     = 2.0D0*ZAB 
               ZAB_UB  = ZAB_UA
            ELSE IF (ABS(UB) > THRESH) THEN
               ONEDU   = 1.0D0/UB
               ZAB     = ZAB + CAB*ONEDU
               ZAB_UB  = -CAB/(UB*UB)
            END IF
            ! now it's function of ZAB

            FZ_Z = 0.0D0
            IF (ABS(ZAB) > THRESH) THEN
               FZ0  = DLOG(1.0D0+ZAB)/ZAB
               FZ   = ZAB*ZAB*(1.0D0-FZ0)
               FZ_Z = 2.0D0*ZAB*(1.0D0-FZ0) + ZAB*ZAB*(FZ0/ZAB -
     $                1.0D0/((1.0D0+ZAB)*ZAB))

            END IF
      FZ_UA = FZ_Z*ZAB_UA
      FZ_UB = FZ_Z*ZAB_UB
      return
      end
