SUBROUTINE UPDATE_VEL_LEAP_FROG

! SUBROUTINE TO UPDATE THE VELOCITIES BY MEANS OF THE VELOCITY VERLET ALGORITHM
! PEDRO OJEDA,  17/SEP/2021
USE PARAMETERS
INTEGER             :: P


	    DO P=1,NUM_RES
		    PARTICLE(P)%COORX(N) = PARTICLE(P)%COORX(O) + DT*PARTICLE(P)%VELX(O) + &
		    0.5D0*DT2*PARTICLE(P)%GRADX(O) 
		    PARTICLE(P)%COORY(N) = PARTICLE(P)%COORY(O) + DT*PARTICLE(P)%VELY(O) + &
		    0.5D0*DT2*PARTICLE(P)%GRADY(O) 
		    PARTICLE(P)%COORZ(N) = PARTICLE(P)%COORZ(O) + DT*PARTICLE(P)%VELZ(O) + &
		    0.5D0*DT2*PARTICLE(P)%GRADZ(O) 
            ENDDO

	    DO P=1,NUM_RES
		    PARTICLE(P)%VELX(N) = PARTICLE(P)%VELX(O) + DTP5 * ( PARTICLE(P)%GRADX(N) + PARTICLE(P)%GRADX(O) )
		    PARTICLE(P)%VELY(N) = PARTICLE(P)%VELY(O) + DTP5 * ( PARTICLE(P)%GRADY(N) + PARTICLE(P)%GRADY(O) )
		    PARTICLE(P)%VELZ(N) = PARTICLE(P)%VELZ(O) + DTP5 * ( PARTICLE(P)%GRADZ(N) + PARTICLE(P)%GRADZ(O) )
            ENDDO

END SUBROUTINE UPDATE_VEL_LEAP_FROG


!SUBROUTINE UPDATE_VEL_LEAP_FROG
!
!! SUBROUTINE TO UPDATE THE VELOCITIES BY MEANS OF THE LEAP-FROG ALGORITHM
!! PEDRO OJEDA,  07/MAY/2011
!USE PARAMETERS
!INTEGER             :: P
!
!
!	    DO 10 P=1,NUM_RES
!		    PARTICLE(P)%VEL(1) = PARTICLE(P)%VEL(1) + DT * PARTICLE(P)%GRAD(1)/PARTICLE(P)%MASS
!		    PARTICLE(P)%VEL(2) = PARTICLE(P)%VEL(2) + DT * PARTICLE(P)%GRAD(2)/PARTICLE(P)%MASS
!		    PARTICLE(P)%VEL(3) = PARTICLE(P)%VEL(3) + DT * PARTICLE(P)%GRAD(3)/PARTICLE(P)%MASS
!10          CONTINUE
!
!END SUBROUTINE UPDATE_VEL_LEAP_FROG



!SUBROUTINE UPDATE_VEL_LEAP_FROG_LANGEVIN
!! SUBROUTINE TO UPDATE THE VELOCITIES BY MEANS OF THE LEAP-FROG ALGORITHM
!! PEDRO OJEDA,  07/MAY/2011
!USE PARAMETERS
!USE Ziggurat
!INTEGER             :: P
!REAL(DP)            :: FACTOR_COCIENT, FACTOR_DENOMINATOR, FACTOR_LANGEVIN
!
!FACTOR_DENOMINATOR = 1.0 + 0.5 * DT 
!FACTOR_COCIENT = ( 1.0 - 0.5 * DT ) / FACTOR_DENOMINATOR
!FACTOR_LANGEVIN = SQRT( 2.0 * GAMMA * KBT_CONST / DT )
!
!	    DO 10 P=1,NUM_RES
!
!		    PARTICLE(P)%VEL(1) = FACTOR_COCIENT * PARTICLE(P)%VEL(1)  &
!                   + DT *( FACTOR_LANGEVIN * RNOR() + PARTICLE(P)%GRAD(1) ) /FACTOR_DENOMINATOR
!
!		    PARTICLE(P)%VEL(2) = FACTOR_COCIENT * PARTICLE(P)%VEL(2)  &
!                   + DT *( FACTOR_LANGEVIN * RNOR() + PARTICLE(P)%GRAD(2) ) /FACTOR_DENOMINATOR
!
!		    PARTICLE(P)%VEL(3) = FACTOR_COCIENT * PARTICLE(P)%VEL(3)  &
!                   + DT *( FACTOR_LANGEVIN * RNOR() + PARTICLE(P)%GRAD(3) ) /FACTOR_DENOMINATOR
!		  
!10          CONTINUE
!
!END SUBROUTINE UPDATE_VEL_LEAP_FROG_LANGEVIN



!SUBROUTINE UPDATE_COOR_LEAP_FROG_BERENDSEN
!! SUBROUTINE TO UPDATE COORDINATES BY MEANS OF THE LEAP-FROG ALGORITHM AND
!! COUPLED TO THE BERENDSEN THERMOSTAT,  H. BERENDSEN, J. CHEM. PHYS. 81, 3684,1984 
!! PEDRO OJEDA,  07/MAY/2011
!
!USE PARAMETERS
!INTEGER             :: P
!
!        CALL PERIODIC_BOUNDARIES_FOLDED
!	CALL FORCES
!        CALL WANG_LANDAU_MD                              ! WANG-LANDAU SIMULATION
!        CALL PERIODIC_BOUNDARIES_FOLDED
!!        CALL COPY_VEL
!        CALL BERENDSEN_FACTOR
!	CALL UPDATE_VEL_LEAP_FROG
!        CALL SCALING_VEL_BERENDSEN
!
!	DO 10 P=1,NUM_RES
!
!	      PARTICLE(P)%COOR(1) = PARTICLE(P)%COOR(1) + DT * PARTICLE(P)%VEL(1)
!
!	      PARTICLE(P)%COOR(2) = PARTICLE(P)%COOR(2) + DT * PARTICLE(P)%VEL(2)
!
!	      PARTICLE(P)%COOR(3) = PARTICLE(P)%COOR(3) + DT * PARTICLE(P)%VEL(3)
!	      
!10      CONTINUE
!
!END SUBROUTINE UPDATE_COOR_LEAP_FROG_BERENDSEN




!SUBROUTINE UPDATE_COOR_LEAP_FROG_BERENDSEN_SHAKE
!! SUBROUTINE TO UPDATE COORDINATES BY MEANS OF THE LEAP-FROG ALGORITHM AND
!! COUPLED TO THE BERENDSEN THERMOSTAT,  H. BERENDSEN, J. CHEM. PHYS. 81, 3684,1984 
!! FOR THE IMPLEMENTATION OF SHAKE WITH THE LEAP-FROG ALGORITHM SEE J. COMP. CHEM. 22, 501, (2001) 
!! PEDRO OJEDA,  07/MAY/2011
!
!USE PARAMETERS
!INTEGER             :: P
!
!
!        CALL PERIODIC_BOUNDARIES_FOLDED
!
!	CALL FORCES
!
!!        CALL COPY_VEL
!
!        CALL COPY_COOR
!
!!        CALL WRITE_COOR
!
!!        CALL WRITE_COOR_BACKUP
!
!!         CALL WRITE_VEL
!
!        CALL BERENDSEN_FACTOR
!
!	CALL UPDATE_VEL_LEAP_FROG
!
!        CALL SCALING_VEL_BERENDSEN
!
!
!!       UPDATE POSITIONS WITHOUT ANY CONSTRAINT
!	DO 10 P=1,NUM_RES
!
!	      PARTICLE(P)%COOR(1) = PARTICLE(P)%COOR(1) + DT * PARTICLE(P)%VEL(1)
!
!	      PARTICLE(P)%COOR(2) = PARTICLE(P)%COOR(2) + DT * PARTICLE(P)%VEL(2)
!
!	      PARTICLE(P)%COOR(3) = PARTICLE(P)%COOR(3) + DT * PARTICLE(P)%VEL(3)
!	      
!10      CONTINUE
!
!        CALL SHAKE
!
!END SUBROUTINE UPDATE_COOR_LEAP_FROG_BERENDSEN_SHAKE



!SUBROUTINE SHAKE
!! SUBROUTINE TO UPDATE COORDINATES BY MEANS OF THE LEAP-FROG ALGORITHM AND
!! COUPLED TO THE BERENDSEN THERMOSTAT,  H. BERENDSEN, J. CHEM. PHYS. 81, 3684,1984 
!! FOR THE IMPLEMENTATION OF SHAKE WITH THE LEAP-FROG ALGORITHM SEE J. COMP. CHEM. 22, 501, (2001) 
!! THERE ARE LINES TAKEN FROM THE F.08 PROGRAM OF TILDESLEY'S BOOK
!! PEDRO OJEDA,  12/MAY/2011
!
!USE PARAMETERS
!REAL(DP),PARAMETER     :: TOL = 0.001
!REAL(DP)               :: DTSQ, TOL2, KIN, WC       
!INTEGER,PARAMETER      :: NA = NUM_RES, NB = NA-1
!LOGICAL     MOVING(NA)
!LOGICAL     MOVED(NA)
!LOGICAL     DONE
!INTEGER,PARAMETER      :: MAXIT=1000
!INTEGER                :: IT
!REAL(DP)               :: PXAB, PYAB, PZAB, PABSQ
!REAL(DP)               :: PXI(NA), PYI(NA), PZI(NA)
!REAL(DP)               :: RXAB, RYAB, RZAB, RABSQ, DIFFSQ, RPAB
!REAL(DP)               :: GAB, DX, DY, DZ, VXIA, VYIA, VZIA
!REAL(DP)               :: RPTOL, RMA, RMB
!REAL(DP)               :: DSQ(NA)
!
!PARAMETER ( RPTOL = 1.0E-2 )
!
!           DTSQ   = DT * DT
!
!	   TOL2   = 2.0 * TOL
!
!	   KIN      = 0.0
!
!	   WC     = 0.0
!
!           DO 100 I = 1, NA
!
!              PXI(I) = PARTICLE(I)%COOR(1)
!
!              PYI(I) = PARTICLE(I)%COOR(2)
!
!              PZI(I) = PARTICLE(I)%COOR(3)
!
!              MOVING(I) = .FALSE.
!
!              MOVED(I)  = .TRUE.
!
!              DSQ(I)    =  6.5*6.5
!
!100        CONTINUE
!
!
!           IT = 0
!
!           DONE = .FALSE.
!
!
!1000       IF ( ( .NOT. DONE ) .AND. ( IT .LE. MAXIT ) ) THEN
!
!              DONE = .TRUE.
!
!              DO 300 I = 1, NB
!
!                 J = I + 1
!
!                 IF ( J .GT. NA ) J = 1
!
!                 IF ( MOVED(I) .OR. MOVED(J) ) THEN
!
!
!			PXAB = PXI(I) - PXI(J)
!
!			PXAB = PXAB - ANINT ( PXAB * BOXI ) * BOXL
!
!			PYAB = PYI(I) - PYI(J)
!
!			PYAB = PYAB - ANINT ( PYAB * BOXI ) * BOXL
!
!			PZAB = PZI(I) - PZI(J)
!
!			PZAB = PZAB - ANINT ( PZAB * BOXI ) * BOXL
!
!			PABSQ  = PXAB *PXAB + PYAB * PYAB + PZAB * PZAB
!
!			RABSQ  = DSQ(I)
!
!			DIFFSQ = RABSQ - PABSQ
!
!			    IF ( ABS(DIFFSQ) .GT. ( RABSQ * TOL2 ) ) THEN
!
!				    RXAB = PARTICLE_BACKUP(I)%COOR(1) - PARTICLE_BACKUP(J)%COOR(1)
!
!				    RXAB = RXAB - ANINT ( RXAB * BOXI ) * BOXL
!
!				    RYAB = PARTICLE_BACKUP(I)%COOR(2) - PARTICLE_BACKUP(J)%COOR(2)
!
!				    RYAB = RYAB - ANINT ( RYAB * BOXI ) * BOXL
!
!				    RZAB = PARTICLE_BACKUP(I)%COOR(3) - PARTICLE_BACKUP(J)%COOR(3)
!
!				    RZAB = RZAB - ANINT ( RZAB * BOXI ) * BOXL
!
!				    RPAB = RXAB * PXAB + RYAB * PYAB + RZAB * PZAB
!
!				    IF ( RPAB .LT. ( RABSQ * RPTOL ) ) THEN
!
!					WRITE(6,*) RPAB, RABSQ
!
!					STOP 'CONSTRAINT FAILURE '
!
!				    ENDIF
!
!				    RMA = 1.0 / PARTICLE(I)%MASS
!
!				    RMB = 1.0 / PARTICLE(J)%MASS
!
!				    GAB = DIFFSQ / ( 2.0 * ( RMA + RMB ) * RPAB )
!
!				    WC  = WC + GAB * RABSQ
!
!				    DX  = RXAB * GAB
!
!				    DY  = RYAB * GAB
!
!				    DZ  = RZAB * GAB
!
!				    PXI(I) = PXI(I) + RMA * DX
!
!				    PYI(I) = PYI(I) + RMA * DY
!
!				    PZI(I) = PZI(I) + RMA * DZ
!
!				    PXI(J) = PXI(J) - RMB * DX
!
!				    PYI(J) = PYI(J) - RMB * DY
!
!				    PZI(J) = PZI(J) - RMB * DZ
!
!				    MOVING(I) = .TRUE.
!
!				    MOVING(J) = .TRUE.
!
!				    DONE = .FALSE.
!
!			    ENDIF
!                 ENDIF
!
!300           CONTINUE
!
!              DO 400 I = 1, NA
!
!                 MOVED(I) = MOVING(I)
!
!                 MOVING(I) = .FALSE.
!
!400           CONTINUE
!
!              IT = IT + 1
!
!              GOTO 1000
!           
!           ENDIF
!
!
!           IF ( .NOT. DONE ) THEN
!
!              WRITE(*,'('' TOO MANY CONSTRAINT ITERATIONS '')')
!
!              STOP
!
!           ENDIF
!
!           DO 500 I = 1, NA
!
!              VXIA = 0.5 * ( PXI(I) - PARTICLE_BACKUP(I)%COOR(1) ) / DT
!
!              VYIA = 0.5 * ( PYI(I) - PARTICLE_BACKUP(I)%COOR(2) ) / DT
!
!              VZIA = 0.5 * ( PZI(I) - PARTICLE_BACKUP(I)%COOR(3) ) / DT
!
!              PARTICLE(I)%VEL(1) = VXIA
!
!              PARTICLE(I)%VEL(2) = VYIA
!
!              PARTICLE(I)%VEL(3) = VZIA
!
!              KIN  = KIN + ( VXIA ** 2 + VYIA ** 2 + VZIA ** 2 ) * PARTICLE(I)%MASS
!
!              PARTICLE_BACKUP(I)%COOR(1) = PARTICLE(I)%COOR(1)
!
!              PARTICLE_BACKUP(I)%COOR(2) = PARTICLE(I)%COOR(2)
!
!              PARTICLE_BACKUP(I)%COOR(3) = PARTICLE(I)%COOR(3) 
!
!              PARTICLE(I)%COOR(1) = PXI(I)
!
!              PARTICLE(I)%COOR(2) = PYI(I)
!
!              PARTICLE(I)%COOR(3) = PZI(I)
!
!500        CONTINUE
!
!
!
!!        TEMP = KIN/N_F 
!
!!        WRITE(6,*) 'TEMPERATURE = ', TEMP
!
!        WC = WC / DTSQ / 3.0
!
!        KIN  = 0.5 * KIN
!
!
!END SUBROUTINE SHAKE




!SUBROUTINE UPDATE_COOR_LEAP_FROG_LANGEVIN
!! SUBROUTINE TO UPDATE COORDINATES BY MEANS OF THE LEAP-FROG ALGORITHM AND
!! COUPLED TO THE BERENDSEN THERMOSTAT,  H. BERENDSEN, J. CHEM. PHYS. 81, 3684,1984 
!! PEDRO OJEDA,  07/MAY/2011
!
!USE PARAMETERS
!INTEGER             :: P
!
!        CALL PERIODIC_BOUNDARIES
!
!	CALL FORCES
!
!        CALL COPY_VEL
!
!	CALL UPDATE_VEL_LEAP_FROG_LANGEVIN
!
!
!	DO 10 P=1,NUM_RES
!
!	      PARTICLE(P)%COOR(1) = PARTICLE(P)%COOR(1) + DT * PARTICLE(P)%VEL(1)
!
!	      PARTICLE(P)%COOR(2) = PARTICLE(P)%COOR(2) + DT * PARTICLE(P)%VEL(2)
!
!	      PARTICLE(P)%COOR(3) = PARTICLE(P)%COOR(3) + DT * PARTICLE(P)%VEL(3)
!	      
!10      CONTINUE
!
!
!END SUBROUTINE UPDATE_COOR_LEAP_FROG_LANGEVIN




!SUBROUTINE COPY_VEL
!! SUBROUTINE TO BACKUP THE VELOCITIES 
!! PEDRO OJEDA,  07/MAY/2011
!
!USE PARAMETERS
!INTEGER             :: P
!
!	DO 10 P=1,NUM_RES
!
!	      PARTICLE_BACKUP(P)%VEL(1) =  PARTICLE(P)%VEL(1)
!
!	      PARTICLE_BACKUP(P)%VEL(2) =  PARTICLE(P)%VEL(2)
!
!	      PARTICLE_BACKUP(P)%VEL(3) =  PARTICLE(P)%VEL(3)
!	      
!10      CONTINUE
!
!
!END SUBROUTINE COPY_VEL




!SUBROUTINE COPY_COOR
!! SUBROUTINE TO BACKUP THE COORDINATES 
!! PEDRO OJEDA,  07/MAY/2011
!
!USE PARAMETERS
!INTEGER             :: P
!
!	DO 10 P=1,NUM_RES
!
!	      PARTICLE_BACKUP(P)%COOR(1) =  PARTICLE(P)%COOR(1)
!
!	      PARTICLE_BACKUP(P)%COOR(2) =  PARTICLE(P)%COOR(2)
!
!	      PARTICLE_BACKUP(P)%COOR(3) =  PARTICLE(P)%COOR(3)
!	      
!10      CONTINUE
!
!END SUBROUTINE COPY_COOR




!SUBROUTINE WRITE_COOR
!! SUBROUTINE TO BACKUP THE COORDINATES 
!! PEDRO OJEDA,  07/MAY/2011
!
!USE PARAMETERS 
!INTEGER            :: P
!
!       WRITE(6,*) NUM_RES
!
!       WRITE(6,*) '   '
!
!	DO 10 P=1,NUM_RES
!
!              WRITE(6,*) 'C',0.3*PARTICLE(P)%COOR(1),0.3*PARTICLE(P)%COOR(2),0.3*PARTICLE(P)%COOR(3)
!	      
!10      CONTINUE
!
!END SUBROUTINE WRITE_COOR
!
!
!
!
!SUBROUTINE WRITE_VEL
!! SUBROUTINE TO BACKUP THE COORDINATES 
!! PEDRO OJEDA,  07/MAY/2011
!
!USE PARAMETERS 
!INTEGER            :: P
!
!        WRITE(6,*) 'WRITING VELOCITIES  '
!
!	DO 10 P=1,NUM_RES
!
!              WRITE(6,*) PARTICLE(P)%VEL(1),PARTICLE(P)%VEL(2),PARTICLE(P)%VEL(3)
!	      
!10      CONTINUE
!
!END SUBROUTINE WRITE_VEL




!SUBROUTINE WRITE_COOR_BACKUP
!! SUBROUTINE TO BACKUP THE COORDINATES 
!! PEDRO OJEDA,  07/MAY/2011
!
!USE PARAMETERS 
!INTEGER            :: P
!
!
!	DO 10 P=1,NUM_RES
!
!              WRITE(6,*) 0.3*PARTICLE_BACKUP(P)%COOR(1),0.3*PARTICLE_BACKUP(P)%COOR(2),0.3*PARTICLE_BACKUP(P)%COOR(3)
!	      
!10      CONTINUE
!
!END SUBROUTINE WRITE_COOR_BACKUP



!SUBROUTINE PERIODIC_BOUNDARIES_FOLDED
!! SUBROUTINE TO MAKE THE PERIODIC BOUNDARY CONDITIONS
!! PEDRO OJEDA,  07/MAY/2011
!
!USE PARAMETERS
!INTEGER             :: P
!
!	DO 10 P=1,NUM_RES
!
!	      PARTICLE_FOLDED(P)%COOR(1) =  PARTICLE(P)%COOR(1)
!
!	      PARTICLE_FOLDED(P)%COOR(2) =  PARTICLE(P)%COOR(2)
!
!	      PARTICLE_FOLDED(P)%COOR(3) =  PARTICLE(P)%COOR(3)
!	      
!10      CONTINUE
!
!
!	DO 20 P=1,NUM_RES
!
!	      PARTICLE_FOLDED(P)%COOR(1) =  PARTICLE_FOLDED(P)%COOR(1) - BOXL * ANINT ( PARTICLE_FOLDED(P)%COOR(1) * BOXI )
!
!	      PARTICLE_FOLDED(P)%COOR(2) =  PARTICLE_FOLDED(P)%COOR(2) - BOXL * ANINT ( PARTICLE_FOLDED(P)%COOR(2) * BOXI )
!
!	      PARTICLE_FOLDED(P)%COOR(3) =  PARTICLE_FOLDED(P)%COOR(3) - BOXL * ANINT ( PARTICLE_FOLDED(P)%COOR(3) * BOXI )
!	      
!20      CONTINUE
!
!END SUBROUTINE PERIODIC_BOUNDARIES_FOLDED
!
!
!
!SUBROUTINE PERIODIC_BOUNDARIES
!! SUBROUTINE TO MAKE THE PERIODIC BOUNDARY CONDITIONS
!! PEDRO OJEDA,  07/MAY/2011
!
!USE PARAMETERS
!INTEGER             :: P
!
!
!	DO 10 P=1,NUM_RES
!
!	      PARTICLE(P)%COOR(1) =  PARTICLE(P)%COOR(1) - BOXL * ANINT ( PARTICLE(P)%COOR(1) * BOXI )
!
!	      PARTICLE(P)%COOR(2) =  PARTICLE(P)%COOR(2) - BOXL * ANINT ( PARTICLE(P)%COOR(2) * BOXI )
!
!	      PARTICLE(P)%COOR(3) =  PARTICLE(P)%COOR(3) - BOXL * ANINT ( PARTICLE(P)%COOR(3) * BOXI )
!	      
!10      CONTINUE
!
!END SUBROUTINE PERIODIC_BOUNDARIES
!
!
!
!
!SUBROUTINE BERENDSEN_FACTOR
!! SUBROUTINE TO COMPUTE THE BERENDSEN FACTOR,  H. BERENDSEN, J. CHEM. PHYS. 81, 3684,1984 
!! PEDRO OJEDA,  07/MAY/2011
!
!USE PARAMETERS
!INTEGER             :: P
!
!        TEMP=0.0
!
!	DO 10 P=1,NUM_RES
!
!	      TEMP = TEMP + PARTICLE(P)%VEL(1)*PARTICLE(P)%VEL(1)  & 
!                     + PARTICLE(P)%VEL(2)*PARTICLE(P)%VEL(2) + PARTICLE(P)%VEL(3)*PARTICLE(P)%VEL(3)
!	     
!10      CONTINUE
!
!        TEMP = TEMP/ N_F 
!
!        LAMBDA_BERENDSEN = SQRT(1.0 + TAU*( TEMP0/TEMP - 1.0) )
!
!END SUBROUTINE BERENDSEN_FACTOR
!
!
!
!
!SUBROUTINE SCALING_VEL_BERENDSEN
!! SUBROUTINE TO RESCALE VELOCITIES BY THE BERENDSEN FACTOR,  H. BERENDSEN, J. CHEM. PHYS. 81, 3684,1984 
!! PEDRO OJEDA,  07/MAY/2011
!
!USE PARAMETERS
!INTEGER             :: P
!
!
!	DO 10 P=1,NUM_RES
!
!		    PARTICLE(P)%VEL(1) = PARTICLE(P)%VEL(1) * LAMBDA_BERENDSEN
!
!		    PARTICLE(P)%VEL(2) = PARTICLE(P)%VEL(2) * LAMBDA_BERENDSEN
!
!		    PARTICLE(P)%VEL(3) = PARTICLE(P)%VEL(3) * LAMBDA_BERENDSEN
!
!10      CONTINUE
!
!
!END SUBROUTINE SCALING_VEL_BERENDSEN
