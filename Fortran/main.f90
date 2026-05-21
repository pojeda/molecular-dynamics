MODULE PARAMETERS
! MODULE WITH THE PARAMETERS NEEDED IN THE SIMULATION
! PEDRO OJEDA,  07/MAY/2011

IMPLICIT NONE
   INTEGER, PARAMETER  :: NUM_RES=1000                     !NUMBER OF AMINOACIDS IN EACH SEQUENCE
   INTEGER, PARAMETER  :: DP = SELECTED_REAL_KIND(12, 60)
   INTEGER, PARAMETER  :: NSTEPS=100000                       !TIME OF SIMULATION
   INTEGER, PARAMETER  :: WRITE_FREQ = 1000                ! FREQUENCY OF WRITING FILES
   REAL(DP),PARAMETER  :: PI=3.141592653589793_DP           !PI CONSTANT
   REAL(DP),PARAMETER  :: DT=0.001_DP                       !TIME STEP
   REAL(DP)            :: DT2=DT*DT                        !TIME STEP SQUARED
   REAL(DP)            :: DTP5=0.5_DP*DT                    !TIME STEP HALFED
   REAL(DP),PARAMETER  :: KBT_CONST = 0.5_DP                !KB * T 
   REAL(DP),PARAMETER  :: EPS_CONST = 40.0_DP              !ARRAY FOR EPSILON
   REAL(DP),PARAMETER  :: SIGMA_CONST = 6.5_DP             !ARRAY OF SIGMA
   REAL(DP),PARAMETER  :: R_CERO=3.8_DP                    !SPRING EQUILIBRIUM SEPARATION
   REAL(DP),PARAMETER  :: A_CONST=50.0_DP                  !SPRING CONSTANT
   REAL(DP)            :: POT_ENER                         !POTENTIAL ENERGY
   REAL(DP)            :: KIN_ENER                         !KINETIC ENERGY
   CHARACTER*1         :: AMINO(4)=(/'S','C','P','N'/)     !ALPHABET OF FOUR LETTERS        
   
   INTEGER             :: OLD = 1, NEW = 2                     !INDICES FOR OLD AND NEW CONFIGURATIONS
   
   ! PARTICLE DATA (STRUCTURE OF ARRAYS)
   REAL(DP), SAVE :: COOR(3,NUM_RES)
   REAL(DP), SAVE :: VEL(3,NUM_RES)
   REAL(DP), SAVE :: FORCE(3,2,NUM_RES)
   REAL(DP), SAVE :: MASS(NUM_RES)
END MODULE PARAMETERS

MODULE RNG_UTILS
  USE PARAMETERS
  IMPLICIT NONE
CONTAINS

  SUBROUTINE INIT_RANDOM_SEED()
    INTEGER :: i, n, clock
    INTEGER, ALLOCATABLE :: seed(:)

    CALL RANDOM_SEED(SIZE=n)
    ALLOCATE(seed(n))

    CALL SYSTEM_CLOCK(COUNT=clock)

    DO i = 1, n
       seed(i) = clock + 37*i
    END DO

    CALL RANDOM_SEED(PUT=seed)
    DEALLOCATE(seed)
  END SUBROUTINE INIT_RANDOM_SEED


  FUNCTION RANDN() RESULT(z)
    REAL(DP) :: z
    REAL(DP) :: u1, u2

    CALL RANDOM_NUMBER(u1)
    CALL RANDOM_NUMBER(u2)

    IF (u1 <= 1.0D-12) u1 = 1.0D-12

    z = SQRT(-2.0_DP*LOG(u1)) * COS(2.0_DP*PI*u2)
  END FUNCTION RANDN

END MODULE RNG_UTILS

PROGRAM MOLECULAR_DYNAMICS
! MAIN PROGRAM TO RUN THE MOLECULAR DYNAMICS
! PEDRO OJEDA,  07/MAY/2011

USE PARAMETERS
USE RNG_UTILS
IMPLICIT NONE
INTEGER             :: NSTEP
INTEGER             :: P
INTEGER             :: TMP
REAL(DP)            :: ENER


  CALL init_random_seed()
  CALL initialize()
  CALL compute_forces()

        OPEN(88,FILE='KOORDINATEN_1T.xyz',STATUS='UNKNOWN')



!MAIN PART OF MOLECULAR DINAMICS
        DO NSTEP=1,NSTEPS

     CALL update_coordinates()
     CALL compute_forces()
     CALL update_velocities()
     CALL kinetic_energy()

     tmp = old
     old = new
     new = tmp

     IF (MOD(nstep, WRITE_FREQ) == 0) THEN
        WRITE(88,*) NUM_RES
        WRITE(88,*) ' '
        DO p = 1, NUM_RES
           WRITE(88,*) 'C', coor(1,p), coor(2,p), coor(3,p)
        ENDDO
        WRITE(*,*) pot_ener, kin_ener, pot_ener + kin_ener
     ENDIF
        ENDDO
        CLOSE(88)

END PROGRAM MOLECULAR_DYNAMICS





SUBROUTINE compute_forces()
  USE PARAMETERS
  IMPLICIT NONE

  force(:,new,:) = 0.0_DP
  pot_ener = 0.0_DP
  kin_ener = 0.0_DP

  CALL spring_force()
  CALL lennard_jones_force()
END SUBROUTINE compute_forces

SUBROUTINE lennard_jones_force()
  USE PARAMETERS
  IMPLICIT NONE

  INTEGER :: i, j
  REAL(DP) :: rij(3), r2
  REAL(DP) :: sr2, sr4, sr6, sr12
  REAL(DP) :: fterm, eij

  DO i = 1, NUM_RES-2
     DO j = i+2, NUM_RES

        rij = coor(:,i) - coor(:,j)
        r2 = SUM(rij*rij)

        sr2  = SIGMA_CONST*SIGMA_CONST/r2
        sr4  = sr2*sr2
        sr6  = sr4*sr2
        sr12 = sr6*sr6

        eij = EPS_CONST*(sr12 - sr6)
        fterm = EPS_CONST*(12.0_DP*sr12 - 6.0_DP*sr6)/r2

        pot_ener = pot_ener + eij

        force(:,new,i) = force(:,new,i) + fterm*rij
        force(:,new,j) = force(:,new,j) - fterm*rij

     ENDDO
  ENDDO

END SUBROUTINE lennard_jones_force

SUBROUTINE spring_force()
  USE PARAMETERS
  IMPLICIT NONE

  INTEGER :: i
  REAL(DP) :: rij(3), r, dr, fij(3)

  DO i = 1, NUM_RES-1

     rij = coor(:,i) - coor(:,i+1)
     r = SQRT(SUM(rij*rij))
     dr = r - R_CERO

     pot_ener = pot_ener + A_CONST*dr*dr

     fij = -2.0_DP*A_CONST*dr*rij/r

     force(:,new,i)   = force(:,new,i)   + fij
     force(:,new,i+1) = force(:,new,i+1) - fij

  ENDDO

END SUBROUTINE spring_force

SUBROUTINE update_coordinates()
  USE PARAMETERS
  IMPLICIT NONE

  INTEGER :: p

  DO p = 1, NUM_RES
     coor(:,p) = coor(:,p) + DT*vel(:,p) + 0.5_DP*DT2*force(:,old,p)/mass(p)
  ENDDO
END SUBROUTINE update_coordinates

SUBROUTINE update_velocities()
  USE PARAMETERS
  IMPLICIT NONE

  INTEGER :: p

  DO p = 1, NUM_RES
     vel(:,p) = vel(:,p) + DTP5*(force(:,new,p) + force(:,old,p))/mass(p)
  ENDDO
END SUBROUTINE update_velocities

SUBROUTINE kinetic_energy()
  USE PARAMETERS
  IMPLICIT NONE

  INTEGER :: p

  kin_ener = 0.0_DP

  DO p = 1, NUM_RES
     kin_ener = kin_ener + mass(p)*SUM(vel(:,p)*vel(:,p))
  ENDDO

  kin_ener = 0.5_DP*kin_ener
END SUBROUTINE kinetic_energy

SUBROUTINE initialize()
  USE PARAMETERS
  USE RNG_UTILS
  IMPLICIT NONE

  INTEGER :: p

  DO p = 1, NUM_RES
     coor(:,p) = [ (p-15)*4.0_DP, 0.0_DP, 0.0_DP ]
     vel(:,p)  = [ randn(), randn(), randn() ]
     force(:,:,p) = 0.0_DP
     mass(p) = 1.0_DP
  ENDDO

  pot_ener = 0.0_DP
  kin_ener = 0.0_DP
END SUBROUTINE initialize


