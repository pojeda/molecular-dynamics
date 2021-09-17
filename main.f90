MODULE PARAMETERS
! MODULE WITH THE PARAMETERS NEEDED IN THE SIMULATION
! PEDRO OJEDA,  07/MAY/2011

IMPLICIT NONE
INTEGER, PARAMETER  :: NUM_RES=30                       !NUMBER OF AMINOACIDS IN EACH SEQUENCE
INTEGER, PARAMETER  :: K4B=SELECTED_INT_KIND(9)
INTEGER, PARAMETER  :: DP = SELECTED_REAL_KIND(12, 60)
INTEGER, PARAMETER  :: TIME=9000                        !TIME OF SIMULATION
INTEGER, SAVE       :: NSTEP
INTEGER             :: CLASE(NUM_RES)                   !SEQUENCES OF AMINOACIDS
INTEGER             :: I,J,K,L,M
INTEGER, PARAMETER  :: SECU=0                           !THE NUMBER OF SEQUENCE                            
REAL(DP),PARAMETER  :: PI=3.141592653589793D0
REAL(DP),PARAMETER  :: DT=0.001                         !TIME STEP
REAL(DP),PARAMETER  :: BOXL=100.0                       !BOX SIZE=BOXL 
REAL(DP),PARAMETER  :: GAMMA=1.0                        !GAMMA=6.5*6.0*PI FROM STOKES LAW 
REAL(DP),PARAMETER  :: TEMP1=9000.0D0 , TEMP2=500.0D0       !TWO TEMPERATURES OF THE SYSTEM
REAL(DP),PARAMETER  :: KBT_CONST = 0.5D0                  ! KB * T 
REAL(DP)            :: BOXI                             !1.0/BOXL 
REAL(DP)            :: EPS_CONST(4,4)                   !ARRAY FOR EPSILON
REAL(DP)            :: SIGMA_CONST(NUM_RES,NUM_RES)     !ARRAY OF SIGMA
CHARACTER*1         :: AMINO(4)=(/'S','C','P','N'/)     !ALPHABET OF FOUR LETTERS        
REAL(DP)            :: POT_ENER



REAL(DP)            :: XDIS(NUM_RES,NUM_RES),YDIS(NUM_RES,NUM_RES),ZDIS(NUM_RES,NUM_RES)
REAL(DP)            :: DIS(NUM_RES,NUM_RES),RDN(NUM_RES,3)




!PARAMETERS FOR THE BERENDSEN THERMOSTAT
REAL(DP)            :: LAMBDA_BERENDSEN
REAL(DP)            :: TEMP
REAL(DP),PARAMETER  :: N_F = 3.0D0*NUM_RES-3.0D0, TAU = 1.0D0 
REAL(DP), SAVE      :: TEMP0 


TYPE PARTICLE_STRUCTURE
REAL(DP)            :: COOR(3),VEL(3), GRAD(3)
REAL(DP)            :: MASS
REAL(DP)            :: CHARGE
REAL(DP)            :: INERTIA
END TYPE PARTICLE_STRUCTURE


TYPE (PARTICLE_STRUCTURE), SAVE, DIMENSION(NUM_RES) :: PARTICLE
!TYPE (PARTICLE_STRUCTURE), SAVE, DIMENSION(NUM_RES) :: PARTICLE_BACKUP
TYPE (PARTICLE_STRUCTURE), SAVE, DIMENSION(NUM_RES) :: PARTICLE_FOLDED

END MODULE PARAMETERS



!MODULE WL_MODULE
!
!! MODULE FOR THE PARAMETERS AND VARIABLES USED BY THE WANG-LANDAU METHOD
!INTEGER, PARAMETER       :: DP = SELECTED_REAL_KIND(12, 60)
!INTEGER, PARAMETER       :: E_MIN = -20000 , E_MAX = 0
!INTEGER                  :: E_BIN
!INTEGER,SAVE             :: E_BIN_OLD, E_BIN_NEW
!REAL(DP), SAVE           :: T_HIST( E_MIN : E_MAX ), ALPHA( E_MIN : E_MAX ), DOS( E_MIN : E_MAX )
!REAL(DP), SAVE           :: CUM_HIST( E_MIN : E_MAX ), VISITS( E_MIN : E_MAX )
!REAL(DP), SAVE           :: GAUSSIAN( -40 : 40 ) 
!REAL(DP), PARAMETER      :: DELTA_ENER = 0.1D0
!REAL(DP), SAVE           :: FACT_F                                   ! LOG(F) = MODIFICATION FACTOR
!REAL(DP), SAVE           :: DF 
!REAL(DP), SAVE           :: SCALE_WL 
!REAL(DP), SAVE           :: OMEGA
!LOGICAL,SAVE             :: ACCEPT 
!INTEGER, SAVE            :: INDICE
!INTEGER, SAVE            :: INTERV_SIZE
!
!END MODULE WL_MODULE




PROGRAM MOLECULAR_DYNAMICS

! MAIN PROGRAM TO RUN THE MOLECULAR DYNAMICS
! PEDRO OJEDA,  07/MAY/2011

USE PARAMETERS
USE Ziggurat
REAL                :: T1, T2,TT1, TT2
INTEGER             :: IMAC
INTEGER             :: P


REAL(DP)            :: ENER
REAL(DP)            :: HIST(-2000:2000),TRANS(-2000:0,-2000:0),CMX,CMY,CMZ
INTEGER	            :: OBIN,EBIN,NEBIN

CHARACTER*40 NAME
CHARACTER*40 NAME4


	IMAC=-35000
        CALL ZIGSET( IMAC )


!READ THE SEQUENCE
        IF(SECU.LT.10) THEN
	   WRITE(NAME4,"('sec_',i1,'.dat')")SECU
        ELSE
	  WRITE(NAME4,"('sec_',i2,'.dat')")SECU
        ENDIF
        OPEN(55,FILE=NAME4,STATUS='UNKNOWN')
	    DO P=1,NUM_RES
		READ(55,*)	CLASE(P)
	    ENDDO
        CLOSE(55)


!ASSIGN THE VALUES FOR THE EPSILONS AND SIGMAS
        CALL SECUENCIAS

        CALL INITIALIZE

        OPEN(88,FILE='min_0.dat',STATUS='UNKNOWN')
                   DO P=1,NUM_RES
                        READ(88,*) PARTICLE(P)%COOR(1),&
                        PARTICLE(P)%COOR(2),PARTICLE(P)%COOR(3)
                   ENDDO
        CLOSE(88)

        CALL FORCES

        write(6,*) 'Potential Energy', POT_ENER

        OPEN(88,FILE='KOORDINATEN_1T.xyz',STATUS='UNKNOWN')

                   WRITE(88,*) NUM_RES
                   WRITE(88,*) '  '
                   DO P=1,NUM_RES
                        WRITE(88,*) 'C', 0.25*PARTICLE(P)%COOR(1),&
                        0.25*PARTICLE(P)%COOR(2),0.25*PARTICLE(P)%COOR(3)
                   ENDDO


!MAIN PART OF MOLECULAR DINAMICS

!        DO NSTEP=1,6000000
!
!!               CALL WRITE_COOR
!               CALL UPDATE_COOR_LEAP_FROG_BERENDSEN
!
!               IF(MOD(NSTEP,100000).EQ.0) THEN
!
!                   WRITE(6,*) TEMP, POT_ENER
!
!                   WRITE(88,*) NUM_RES
!                   WRITE(88,*) '  '
!                   DO P=1,NUM_RES
!                        WRITE(88,*) 'C', 0.25*PARTICLE_FOLDED(P)%COOR(1), &
!                        0.25*PARTICLE_FOLDED(P)%COOR(2),0.25*PARTICLE_FOLDED(P)%COOR(3)
!                   ENDDO
!               ENDIF
!
!               IF(MOD(NSTEP,10000000) .EQ. 0 ) THEN
!                      CALL WRITE_WANG_LANDAU_MD_HIST
!                      CALL RE_INITIALIZE
!               ENDIF
!
!               IF(FACT_F < 0.0004) THEN
!                     IF(MOD(NSTEP,3000000) .EQ. 0 ) THEN
!                      !CALL WRITE_WANG_LANDAU_MD_HIST
!                         CALL RE_INITIALIZE
!                     ENDIF
!               ENDIF
!
!	 ENDDO
        CLOSE(88)

!         Activate WL sampling
!         CALL WRITE_WANG_LANDAU_MD_HIST

END PROGRAM MOLECULAR_DYNAMICS
