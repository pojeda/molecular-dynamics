MODULE PARAMETERS
! MODULE WITH THE PARAMETERS NEEDED IN THE SIMULATION
! PEDRO OJEDA,  07/MAY/2011

IMPLICIT NONE
   INTEGER, PARAMETER  :: NUM_RES=30                       !NUMBER OF AMINOACIDS IN EACH SEQUENCE
   INTEGER, PARAMETER  :: DP = SELECTED_REAL_KIND(12, 60)
   INTEGER, PARAMETER  :: TIME=90000                       !TIME OF SIMULATION
   INTEGER, SAVE       :: NSTEP
   INTEGER             :: CLASE(NUM_RES)                   !SEQUENCES OF AMINOACIDS
   INTEGER             :: I,J,K,L,M
   INTEGER, PARAMETER  :: SECU=0                           !THE NUMBER OF SEQUENCE                            
   REAL(DP),PARAMETER  :: PI=3.141592653589793D0           !PI CONSTANT
   REAL(DP),PARAMETER  :: DT=0.0001D0                      !TIME STEP
   REAL(DP)            :: DT2=DT*DT                        !TIME STEP SQUARED
   REAL(DP)            :: DTP5=0.5D0*DT                    !TIME STEP HALFED
   REAL(DP),PARAMETER  :: KBT_CONST = 0.5D0                !KB * T 
   REAL(DP)            :: EPS_CONST(4,4)                   !ARRAY FOR EPSILON
   REAL(DP)            :: SIGMA_CONST(NUM_RES,NUM_RES)     !ARRAY OF SIGMA
   REAL(DP),PARAMETER  :: R_CERO=3.8D0                     !SPRING EQUILIBRIUM SEPARATION
   REAL(DP),PARAMETER  :: A_CONST=50.0D0                   !SPRING CONSTANT
   REAL(DP)            :: POT_ENER                         !POTENTIAL ENERGY
   REAL(DP)            :: KIN_ENER                         !KINETIC ENERGY
   CHARACTER*1         :: AMINO(4)=(/'S','C','P','N'/)     !ALPHABET OF FOUR LETTERS        
   
   INTEGER,SAVE        :: O = 1, N = 1                     !INDICES FOR OLD AND NEW CONFIGURATIONS
   
   TYPE PARTICLE_STRUCTURE
   REAL(DP)            :: COORX,COORY,COORZ
   REAL(DP)            :: VELX,VELY,VELZ
   REAL(DP)            :: GRADX(2),GRADY(2),GRADZ(2)
   REAL(DP)            :: MASS
   REAL(DP)            :: CHARGE
   REAL(DP)            :: INERTIA
   END TYPE PARTICLE_STRUCTURE
   
   TYPE (PARTICLE_STRUCTURE), SAVE, DIMENSION(NUM_RES) :: PARTICLE
END MODULE PARAMETERS


PROGRAM MOLECULAR_DYNAMICS
! MAIN PROGRAM TO RUN THE MOLECULAR DYNAMICS
! PEDRO OJEDA,  07/MAY/2011

USE PARAMETERS
USE Ziggurat
INTEGER             :: IMAC
INTEGER             :: P
INTEGER             :: TMP
REAL(DP)            :: ENER

CHARACTER*40 NAME4

IMAC=-35000
CALL ZIGSET( IMAC )


!READ THE SEQUENCE
IF(SECU.LT.10) THEN
   WRITE(NAME4,"('../Data/sec_',i1,'.dat')")SECU
ELSE
   WRITE(NAME4,"('../Data/sec_',i2,'.dat')")SECU
ENDIF
OPEN(55,FILE=NAME4,STATUS='UNKNOWN')
DO P=1,NUM_RES
    READ(55,*)	CLASE(P)
ENDDO
CLOSE(55)


!ASSIGN THE VALUES FOR THE EPSILONS AND SIGMAS
        CALL SECUENCIAS

        CALL INITIALIZE

        OPEN(88,FILE='../Data/ini.dat',STATUS='UNKNOWN')
                   DO P=1,NUM_RES
                        READ(88,*) PARTICLE(P)%COORX,&
                        PARTICLE(P)%COORY,PARTICLE(P)%COORZ
                   ENDDO
        CLOSE(88)

        CALL FORCES


        OPEN(88,FILE='KOORDINATEN_1T.xyz',STATUS='UNKNOWN')



!MAIN PART OF MOLECULAR DINAMICS
        N = 2
        DO NSTEP=1,10000000

               CALL UPDATE_VEL_VERLET_COOR
               CALL FORCES
               CALL UPDATE_VEL_VERLET_VEL
               CALL KINETIC

               TMP = O
               O = N 
               N = TMP

               IF(MOD(NSTEP,10000).EQ.0) THEN 
                   WRITE(88,*) NUM_RES
                   WRITE(88,*) '  '
                   DO P=1,NUM_RES
                        WRITE(88,*) 'C', PARTICLE(P)%COORX, PARTICLE(P)%COORY,PARTICLE(P)%COORZ
                   ENDDO
                   WRITE(6,*) POT_ENER, KIN_ENER, POT_ENER + KIN_ENER
               ENDIF

        ENDDO
        CLOSE(88)

END PROGRAM MOLECULAR_DYNAMICS



SUBROUTINE LENNARD_JONES 
! SUBROUTINE TO COMPUTE THE LENNARD-JONES POTENTIAL AND FORCE
! PEDRO OJEDA,  07/MAY/2011

USE PARAMETERS
REAL(DP)            :: FDUMMY1,FDUMMY2,FDUMMY4,FDUMMY5
REAL(DP)            :: FDUMMYD,FDUMMYE
REAL(DP)            :: FTERM
REAL(DP)            :: ENER1
REAL(DP)            :: DISX,DISY,DISZ,DISD

ENER1=0.0D0

!VAN DER WAALS CONTRIBUTION
DO I=1,NUM_RES-2
     DO J=I+2,NUM_RES
        DISX=PARTICLE(I)%COORX-PARTICLE(J)%COORX                          
        DISY=PARTICLE(I)%COORY-PARTICLE(J)%COORY                          
        DISZ=PARTICLE(I)%COORZ-PARTICLE(J)%COORZ                           
        DISD=DISX*DISX + DISY*DISY + DISZ*DISZ        

        FDUMMY1=SIGMA_CONST(I,J)*SIGMA_CONST(I,J)/DISD
        K=CLASE(I)
        L=CLASE(J)
        FDUMMY2=FDUMMY1*FDUMMY1                             ! POWER 4
        FDUMMY4=FDUMMY2*FDUMMY1                             ! POWER 6
        FDUMMY5=FDUMMY4*FDUMMY4                             ! POWER 12
        FDUMMYD=1.0D0/DISD                                  ! 1/DISTANCE_IJ^2
        FDUMMYE=FDUMMY5-FDUMMY4                             ! (SIGMA/R)^12 - (SIGMA/R)^6
        ENER1=ENER1+EPS_CONST(K,L)*FDUMMYE                  ! ENERGY = EPSILON*[(SIGMA/R)^12 - (SIGMA/R)^6 ]
        FTERM=(12.0D0*FDUMMY5-6.0D0*FDUMMY4)*FDUMMYD        ! FORCE = -DV/DR = [12*(SIGMA/R)^12 - 6*(SIGMA/R)^6 ]*(1/R^2)

        PARTICLE(I)%GRADX(N)=PARTICLE(I)%GRADX(N)+EPS_CONST(K,L)*FTERM*DISX   ! X COMPONENT OF THE GRADIENT
        PARTICLE(I)%GRADY(N)=PARTICLE(I)%GRADY(N)+EPS_CONST(K,L)*FTERM*DISY   ! Y COMPONENT OF THE GRADIENT
        PARTICLE(I)%GRADZ(N)=PARTICLE(I)%GRADZ(N)+EPS_CONST(K,L)*FTERM*DISZ   ! Z COMPONENT OF THE GRADIENT
        PARTICLE(J)%GRADX(N)=PARTICLE(J)%GRADX(N)-EPS_CONST(K,L)*FTERM*DISX   ! X COMPONENT OF THE GRADIENT
        PARTICLE(J)%GRADY(N)=PARTICLE(J)%GRADY(N)-EPS_CONST(K,L)*FTERM*DISY   ! Y COMPONENT OF THE GRADIENT
        PARTICLE(J)%GRADZ(N)=PARTICLE(J)%GRADZ(N)-EPS_CONST(K,L)*FTERM*DISZ   ! Z COMPONENT OF THE GRADIENT
     ENDDO
ENDDO

POT_ENER = POT_ENER + ENER1

END SUBROUTINE LENNARD_JONES


SUBROUTINE SPRING 
! SUBROUTINE TO COMPUTE THE HARMONIC OSCILLATOR POTENTIAL AND FORCE
! PEDRO OJEDA,  07/MAY/2011

USE PARAMETERS

REAL(DP)            :: FDUMMY1,FDUMMY2,FDUMMY3,FDUMMY4
REAL(DP)            :: ENER2
REAL(DP)            :: DISX,DISY,DISZ,DISPX,DISD,DISPY,DISPZ,DISPD

ENER2=0.0D0

!HARMONIC POTENTIAL
FDUMMY3=2.0D0*A_CONST

DO I=2,NUM_RES-1
    DISX=PARTICLE(I)%COORX-PARTICLE(I+1)%COORX                          
    DISY=PARTICLE(I)%COORY-PARTICLE(I+1)%COORY                          
    DISZ=PARTICLE(I)%COORZ-PARTICLE(I+1)%COORZ                           
    DISD=SQRT(DISX*DISX + DISY*DISY + DISZ*DISZ)        

    DISPX=PARTICLE(I)%COORX-PARTICLE(I-1)%COORX                          
    DISPY=PARTICLE(I)%COORY-PARTICLE(I-1)%COORY                          
    DISPZ=PARTICLE(I)%COORZ-PARTICLE(I-1)%COORZ                           
    DISPD=SQRT(DISPX*DISPX + DISPY*DISPY + DISPZ*DISPZ)       

    FDUMMY1=(DISD-R_CERO)
    FDUMMY2=FDUMMY1*FDUMMY1
    ENER2=ENER2+FDUMMY2
    FDUMMY4=(DISPD-R_CERO)
    PARTICLE(I)%GRADX(N)=PARTICLE(I)%GRADX(N)-FDUMMY3*DISX*FDUMMY1/DISD - FDUMMY3*DISPX*FDUMMY4/DISPD
    PARTICLE(I)%GRADY(N)=PARTICLE(I)%GRADY(N)-FDUMMY3*DISY*FDUMMY1/DISD - FDUMMY3*DISPY*FDUMMY4/DISPD
    PARTICLE(I)%GRADZ(N)=PARTICLE(I)%GRADZ(N)-FDUMMY3*DISZ*FDUMMY1/DISD - FDUMMY3*DISPZ*FDUMMY4/DISPD
ENDDO


    DISX=PARTICLE(1)%COORX-PARTICLE(2)%COORX                          
    DISY=PARTICLE(1)%COORY-PARTICLE(2)%COORY                          
    DISZ=PARTICLE(1)%COORZ-PARTICLE(2)%COORZ                           
    DISD=SQRT(DISX*DISX + DISY*DISY + DISZ*DISZ)       
    FDUMMY1=DISD-R_CERO
    FDUMMY2=FDUMMY1*FDUMMY1
    ENER2=ENER2+FDUMMY2
    PARTICLE(1)%GRADX(N)=PARTICLE(1)%GRADX(N)-FDUMMY3*DISX*FDUMMY1/DISD
    PARTICLE(1)%GRADY(N)=PARTICLE(1)%GRADY(N)-FDUMMY3*DISY*FDUMMY1/DISD
    PARTICLE(1)%GRADZ(N)=PARTICLE(1)%GRADZ(N)-FDUMMY3*DISZ*FDUMMY1/DISD


    DISX=PARTICLE(NUM_RES)%COORX-PARTICLE(NUM_RES-1)%COORX                          
    DISY=PARTICLE(NUM_RES)%COORY-PARTICLE(NUM_RES-1)%COORY                          
    DISZ=PARTICLE(NUM_RES)%COORZ-PARTICLE(NUM_RES-1)%COORZ                           
    DISD=SQRT(DISX*DISX + DISY*DISY + DISZ*DISZ)       
    FDUMMY1=DISD-R_CERO
    PARTICLE(NUM_RES)%GRADX(N)=PARTICLE(NUM_RES)%GRADX(N)-FDUMMY3*DISX*FDUMMY1/DISD
    PARTICLE(NUM_RES)%GRADY(N)=PARTICLE(NUM_RES)%GRADY(N)-FDUMMY3*DISY*FDUMMY1/DISD
    PARTICLE(NUM_RES)%GRADZ(N)=PARTICLE(NUM_RES)%GRADZ(N)-FDUMMY3*DISZ*FDUMMY1/DISD

POT_ENER = POT_ENER + A_CONST*ENER2

END SUBROUTINE SPRING


SUBROUTINE KINETIC 
! SUBROUTINE TO COMPUTE THE KINETIC ENERGY
! PEDRO OJEDA,  17/SEP/2021

USE PARAMETERS

REAL(DP)            :: ENER2

ENER2=0.0D0

DO I=1,NUM_RES
      ENER2 = ENER2 + PARTICLE(I)%VELX**2 + PARTICLE(I)%VELY**2 + PARTICLE(I)%VELZ**2 
ENDDO

      KIN_ENER = 0.5D0*ENER2

END SUBROUTINE KINETIC


SUBROUTINE FORCES
! SUBROUTINE TO CALCULATE THE FORCES DEPENDING ON THE INTERACTIONS DESIRED
! PEDRO OJEDA,  07/MAY/2011

USE PARAMETERS 

  CALL INITIALIZE_GRAD
  CALL SPRING
  CALL LENNARD_JONES

END SUBROUTINE FORCES


SUBROUTINE UPDATE_VEL_VERLET_COOR
! SUBROUTINE TO UPDATE THE VELOCITIES BY MEANS OF THE VELOCITY VERLET ALGORITHM
! PEDRO OJEDA,  17/SEP/2021
USE PARAMETERS
INTEGER             :: P

DO P=1,NUM_RES
      PARTICLE(P)%COORX = PARTICLE(P)%COORX + DT*PARTICLE(P)%VELX + 0.5D0*DT2*PARTICLE(P)%GRADX(O) 
      PARTICLE(P)%COORY = PARTICLE(P)%COORY + DT*PARTICLE(P)%VELY + 0.5D0*DT2*PARTICLE(P)%GRADY(O) 
      PARTICLE(P)%COORZ = PARTICLE(P)%COORZ + DT*PARTICLE(P)%VELZ + 0.5D0*DT2*PARTICLE(P)%GRADZ(O)
ENDDO

END SUBROUTINE UPDATE_VEL_VERLET_COOR


SUBROUTINE UPDATE_VEL_VERLET_VEL
! SUBROUTINE TO UPDATE THE VELOCITIES BY MEANS OF THE VELOCITY VERLET ALGORITHM
! PEDRO OJEDA,  17/SEP/2021
USE PARAMETERS
INTEGER             :: P

DO P=1,NUM_RES
    PARTICLE(P)%VELX = PARTICLE(P)%VELX + DTP5 * ( PARTICLE(P)%GRADX(N) + PARTICLE(P)%GRADX(O) )
    PARTICLE(P)%VELY = PARTICLE(P)%VELY + DTP5 * ( PARTICLE(P)%GRADY(N) + PARTICLE(P)%GRADY(O) )
    PARTICLE(P)%VELZ = PARTICLE(P)%VELZ + DTP5 * ( PARTICLE(P)%GRADZ(N) + PARTICLE(P)%GRADZ(O) )
ENDDO

END SUBROUTINE UPDATE_VEL_VERLET_VEL


SUBROUTINE INITIALIZE
! SUBROUTINE TO INITIALIZE THE COORDINATES, VELOCITIES AND GRADIENT
! PEDRO OJEDA,  07/MAY/2011

USE PARAMETERS
USE Ziggurat

INTEGER             :: P

DO P=1,NUM_RES
    PARTICLE(P)%COORX=  (P-15)*4.0D0
    PARTICLE(P)%COORY=  0.0D0
    PARTICLE(P)%COORZ=  0.0D0

    PARTICLE(P)%COORX=  PARTICLE(P)%COORX
    PARTICLE(P)%COORY=  PARTICLE(P)%COORY
    PARTICLE(P)%COORZ=  PARTICLE(P)%COORZ

    PARTICLE(P)%VELX=  0.001D0*RNOR()
    PARTICLE(P)%VELY=  0.001D0*RNOR()
    PARTICLE(P)%VELZ=  0.001D0*RNOR()

    PARTICLE(P)%GRADX(N)=  0.0D0
    PARTICLE(P)%GRADY(N)=  0.0D0
    PARTICLE(P)%GRADZ(N)=  0.0D0

    PARTICLE(P)%GRADX(O)=  0.0D0
    PARTICLE(P)%GRADY(O)=  0.0D0
    PARTICLE(P)%GRADZ(O)=  0.0D0

    PARTICLE(P)%MASS=  1.0D0

ENDDO

POT_ENER = 0.0

END SUBROUTINE INITIALIZE


SUBROUTINE INITIALIZE_GRAD
! SUBROUTINE TO INITIALIZE THE GRADIENT AND POTENTIAL ENERGY
! PEDRO OJEDA,  07/MAY/2011

USE PARAMETERS

INTEGER             :: P

DO P=1,NUM_RES
    PARTICLE(P)%GRADX(N)=  0.0
    PARTICLE(P)%GRADY(N)=  0.0
    PARTICLE(P)%GRADZ(N)=  0.0
ENDDO

POT_ENER = 0.0D0         
KIN_ENER = 0.0D0           

END SUBROUTINE INITIALIZE_GRAD

