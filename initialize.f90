SUBROUTINE INITIALIZE
! SUBROUTINE TO INITIALIZE THE COORDINATES, VELOCITIES AND GRADIENT
! PEDRO OJEDA,  07/MAY/2011

USE PARAMETERS

USE Ziggurat

INTEGER             :: P

             DO P=1,NUM_RES
	         PARTICLE(P)%COORX(1)=  (P-15)*6.5D0
	         PARTICLE(P)%COORY(1)=  0.0D0
	         PARTICLE(P)%COORZ(1)=  0.0D0

	         PARTICLE(P)%COORX(2)=  PARTICLE(P)%COORX(1)
	         PARTICLE(P)%COORY(2)=  PARTICLE(P)%COORY(1)
	         PARTICLE(P)%COORZ(2)=  PARTICLE(P)%COORZ(1)

	         PARTICLE_FOLDED(P)%COORX(1)=  PARTICLE(P)%COORX(1)
	         PARTICLE_FOLDED(P)%COORY(1)=  PARTICLE(P)%COORY(1)
	         PARTICLE_FOLDED(P)%COORZ(1)=  PARTICLE(P)%COORZ(1)


	         PARTICLE(P)%VELX(1)=  2.0D0*UNI()-1.0D0
	         PARTICLE(P)%VELY(1)=  2.0D0*UNI()-1.0D0
	         PARTICLE(P)%VELZ(1)=  2.0D0*UNI()-1.0D0

	         PARTICLE(P)%VELX(2)=  PARTICLE(P)%VELX(1)
	         PARTICLE(P)%VELY(2)=  PARTICLE(P)%VELY(1)
	         PARTICLE(P)%VELZ(2)=  PARTICLE(P)%VELZ(1)

	         PARTICLE(P)%GRADX(1)=  0.0D0
	         PARTICLE(P)%GRADY(1)=  0.0D0
	         PARTICLE(P)%GRADZ(1)=  0.0D0

	         PARTICLE(P)%GRADX(2)=  PARTICLE(P)%GRADX(1)
	         PARTICLE(P)%GRADY(2)=  PARTICLE(P)%GRADY(1)
	         PARTICLE(P)%GRADZ(2)=  PARTICLE(P)%GRADZ(1)

	         PARTICLE(P)%MASS=  1.0D0

             ENDDO

             POT_ENER = 0.0
             BOXI = 1.0 / BOXL

END SUBROUTINE INITIALIZE


!SUBROUTINE RE_INITIALIZE
!
!! SUBROUTINE TO INITIALIZE THE COORDINATES, VELOCITIES AND GRADIENT
!! PEDRO OJEDA,  07/MAY/2011
!
!USE PARAMETERS
!
!USE Ziggurat
!
!INTEGER             :: P
!
!             DO 10 P=1,NUM_RES
!	         PARTICLE(P)%COOR(1)=  (P-15)*6.5
!	         PARTICLE(P)%COOR(2)=  0.0
!	         PARTICLE(P)%COOR(3)=  0.0
!
!
!	         PARTICLE_FOLDED(P)%COOR(1)=  PARTICLE(P)%COOR(1)
!	         PARTICLE_FOLDED(P)%COOR(2)=  PARTICLE(P)%COOR(2)
!	         PARTICLE_FOLDED(P)%COOR(3)=  PARTICLE(P)%COOR(3)
!
!
!	         PARTICLE(P)%VEL(1)=  2.0*UNI()-1.0
!	         PARTICLE(P)%VEL(2)=  2.0*UNI()-1.0
!	         PARTICLE(P)%VEL(3)=  2.0*UNI()-1.0
!
!
!	         PARTICLE(P)%GRAD(1)=  0.0
!	         PARTICLE(P)%GRAD(2)=  0.0
!	         PARTICLE(P)%GRAD(3)=  0.0
!
!	         PARTICLE(P)%MASS=  1.0
!
!10           CONTINUE
!
!
!END SUBROUTINE RE_INITIALIZE



SUBROUTINE INITIALIZE_GRAD

! SUBROUTINE TO INITIALIZE THE GRADIENT AND POTENTIAL ENERGY
! PEDRO OJEDA,  07/MAY/2011

USE PARAMETERS

INTEGER             :: P

             DO P=1,NUM_RES
	         PARTICLE(P)%GRADX(1)=  0.0
	         PARTICLE(P)%GRADY(1)=  0.0
	         PARTICLE(P)%GRADZ(1)=  0.0
             ENDDO

             POT_ENER = 0.0D0         
             KIN_ENER = 0.0D0           

END SUBROUTINE INITIALIZE_GRAD


