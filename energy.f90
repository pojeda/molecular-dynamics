!SUBROUTINE DELTA_MATRIX 
!
!! SUBROUTINE TO COMPUTE THE DISTANCES AND VECTORS BETWEEN "NUM_RES" PARTICLES
!! PEDRO OJEDA,  07/MAY/2011
!USE PARAMETERS                                               
!REAL(DP)            ::  FDUMMY1,FDUMMY2,FDUMMY3
!
!DO I=1,NUM_RES  
!  DO J=I,NUM_RES
!        DIS(I,J)=0.0D0   
!        DIS(J,I)=0.0D0 
!  ENDDO                
!ENDDO                  
!
!!THE DISTANCES BETWEEN AMINO ACIDS ARE COMPUTED HERE 
!DO 10 I=1,NUM_RES-1                                    
!  DO 20 J=I+1,NUM_RES                                     
!    FDUMMY1=PARTICLE(I)%COORX(1)-PARTICLE(J)%COORX(1)                          
!    FDUMMY2=PARTICLE(I)%COORY(1)-PARTICLE(J)%COORY(1)                          
!    FDUMMY3=PARTICLE(I)%COORZ(1)-PARTICLE(J)%COORZ(1)                           
!            
!
!    XDIS(I,J)=FDUMMY1
!    YDIS(I,J)=FDUMMY2
!    ZDIS(I,J)=FDUMMY3
!    XDIS(J,I)=-FDUMMY1
!    YDIS(J,I)=-FDUMMY2
!    ZDIS(J,I)=-FDUMMY3
!    DIS(I,J)=SQRT(FDUMMY1*FDUMMY1+FDUMMY2*FDUMMY2+FDUMMY3*FDUMMY3)
!    DIS(J,I)=DIS(I,J)
!
!20  CONTINUE
!
!10  CONTINUE
!
!END SUBROUTINE DELTA_MATRIX






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

    DISX=PARTICLE(I)%COORX(1)-PARTICLE(J)%COORX(1)                          
    DISY=PARTICLE(I)%COORY(1)-PARTICLE(J)%COORY(1)                          
    DISZ=PARTICLE(I)%COORZ(1)-PARTICLE(J)%COORZ(1)                           
    DISD=DISX*DISX + DISY*DISY + DISZ*DISZ        


           FDUMMY1=SIGMA_CONST(I,J)*SIGMA_CONST(I,J)/DISD
           K=CLASE(I)
           L=CLASE(J)
	   FDUMMY2=FDUMMY1*FDUMMY1                             ! POWER 4
	   FDUMMY4=FDUMMY2*FDUMMY1                             ! POWER 6
	   FDUMMY5=FDUMMY4*FDUMMY4                             ! POWER 12
	   FDUMMYD=1.0/DISD                                    ! 1/DISTANCE_IJ^2
	   FDUMMYE=FDUMMY5-FDUMMY4                             ! (SIGMA/R)^12 - (SIGMA/R)^6
	   ENER1=ENER1+EPS_CONST(K,L)*FDUMMYE                  ! ENERGY = EPSILON*[(SIGMA/R)^12 - (SIGMA/R)^6 ]
	   FTERM=(12.*FDUMMY5-6.*FDUMMY4)*FDUMMYD              ! FORCE = -DV/DR = [12*(SIGMA/R)^12 - 6*(SIGMA/R)^6 ]*(1/R^2)

           PARTICLE(I)%GRADX(1)=PARTICLE(I)%GRADX(1)+EPS_CONST(K,L)*FTERM*DISX   ! X COMPONENT OF THE GRADIENT
           PARTICLE(I)%GRADY(1)=PARTICLE(I)%GRADY(1)+EPS_CONST(K,L)*FTERM*DISY   ! Y COMPONENT OF THE GRADIENT
           PARTICLE(I)%GRADZ(1)=PARTICLE(I)%GRADZ(1)+EPS_CONST(K,L)*FTERM*DISZ   ! Z COMPONENT OF THE GRADIENT
           PARTICLE(J)%GRADX(1)=PARTICLE(J)%GRADX(1)-EPS_CONST(K,L)*FTERM*DISX   ! X COMPONENT OF THE GRADIENT
           PARTICLE(J)%GRADY(1)=PARTICLE(J)%GRADY(1)-EPS_CONST(K,L)*FTERM*DISY   ! Y COMPONENT OF THE GRADIENT
           PARTICLE(J)%GRADZ(1)=PARTICLE(J)%GRADZ(1)-EPS_CONST(K,L)*FTERM*DISZ   ! Z COMPONENT OF THE GRADIENT

     ENDDO

ENDDO

POT_ENER = POT_ENER + ENER1

END SUBROUTINE LENNARD_JONES


SUBROUTINE SPRING 
! SUBROUTINE TO COMPUTE THE HARMONIC OSCILLATOR POTENTIAL AND FORCE
! PEDRO OJEDA,  07/MAY/2011

USE PARAMETERS

REAL(DP)            :: FDUMMY1,FDUMMY2,FDUMMY3,FDUMMY4,FDUMMY5
REAL(DP)            :: FDUMMY0,FDUMMYD,FDUMMYE
REAL(DP)            :: FTERM
REAL(DP)            :: ENER2,ENER
REAL(DP)            :: R_CERO,A_CONST
REAL(DP)            :: DISX,DISY,DISZ,DISPX,DISD,DISPY,DISPZ,DISPD



	R_CERO=3.8D0
	A_CONST=50.0D0
	ENER2=0.0D0


!HARMONIC POTENTIAL
FDUMMY3=2.0D0*A_CONST

DO I=2,NUM_RES-1

    DISX=PARTICLE(I)%COORX(1)-PARTICLE(I+1)%COORX(1)                          
    DISY=PARTICLE(I)%COORY(1)-PARTICLE(I+1)%COORY(1)                          
    DISZ=PARTICLE(I)%COORZ(1)-PARTICLE(I+1)%COORZ(1)                           
    DISD=SQRT(DISX*DISX + DISY*DISY + DISZ*DISZ)        

    DISPX=PARTICLE(I)%COORX(1)-PARTICLE(I-1)%COORX(1)                          
    DISPY=PARTICLE(I)%COORY(1)-PARTICLE(I-1)%COORY(1)                          
    DISPZ=PARTICLE(I)%COORZ(1)-PARTICLE(I-1)%COORZ(1)                           
    DISPD=SQRT(DISPX*DISPX + DISPY*DISPY + DISPZ*DISPZ)       

      FDUMMY1=(DISD-R_CERO)
      FDUMMY2=FDUMMY1*FDUMMY1
      ENER2=ENER2+FDUMMY2
      FDUMMY4=(DISPD-R_CERO)
      PARTICLE(I)%GRADX(1)=PARTICLE(I)%GRADX(1)-FDUMMY3*DISX*FDUMMY1/DISD - FDUMMY3*DISPX*FDUMMY4/DISPD
      PARTICLE(I)%GRADY(1)=PARTICLE(I)%GRADY(1)-FDUMMY3*DISY*FDUMMY1/DISD - FDUMMY3*DISPY*FDUMMY4/DISPD
      PARTICLE(I)%GRADZ(1)=PARTICLE(I)%GRADZ(1)-FDUMMY3*DISZ*FDUMMY1/DISD - FDUMMY3*DISPZ*FDUMMY4/DISPD
ENDDO


    DISX=PARTICLE(1)%COORX(1)-PARTICLE(2)%COORX(1)                          
    DISY=PARTICLE(1)%COORY(1)-PARTICLE(2)%COORY(1)                          
    DISZ=PARTICLE(1)%COORZ(1)-PARTICLE(2)%COORZ(1)                           
    DISD=SQRT(DISX*DISX + DISY*DISY + DISZ*DISZ)       
      FDUMMY1=DISD-R_CERO
      FDUMMY2=FDUMMY1*FDUMMY1
      ENER2=ENER2+FDUMMY2
      PARTICLE(1)%GRADX(1)=PARTICLE(1)%GRADX(1)-FDUMMY3*DISX*FDUMMY1/DISD
      PARTICLE(1)%GRADY(1)=PARTICLE(1)%GRADY(1)-FDUMMY3*DISY*FDUMMY1/DISD
      PARTICLE(1)%GRADZ(1)=PARTICLE(1)%GRADZ(1)-FDUMMY3*DISZ*FDUMMY1/DISD


    DISX=PARTICLE(NUM_RES)%COORX(1)-PARTICLE(NUM_RES-1)%COORX(1)                          
    DISY=PARTICLE(NUM_RES)%COORY(1)-PARTICLE(NUM_RES-1)%COORY(1)                          
    DISZ=PARTICLE(NUM_RES)%COORZ(1)-PARTICLE(NUM_RES-1)%COORZ(1)                           
    DISD=SQRT(DISX*DISX + DISY*DISY + DISZ*DISZ)       
      FDUMMY1=DISD-R_CERO
      PARTICLE(NUM_RES)%GRADX(1)=PARTICLE(NUM_RES)%GRADX(1)-FDUMMY3*DISX*FDUMMY1/DISD
      PARTICLE(NUM_RES)%GRADY(1)=PARTICLE(NUM_RES)%GRADY(1)-FDUMMY3*DISY*FDUMMY1/DISD
      PARTICLE(NUM_RES)%GRADZ(1)=PARTICLE(NUM_RES)%GRADZ(1)-FDUMMY3*DISZ*FDUMMY1/DISD

      POT_ENER = POT_ENER + A_CONST*ENER2

END SUBROUTINE SPRING


SUBROUTINE KINETIC 
! SUBROUTINE TO COMPUTE THE KINETIC ENERGY
! PEDRO OJEDA,  17/SEP/2021

USE PARAMETERS

REAL(DP)            :: ENER2

	ENER2=0.0D0

DO I=1,NUM_RES
      ENER2 = ENER2 + PARTICLE(I)%VELX(N)**2 + PARTICLE(I)%VELY(N)**2 + PARTICLE(I)%VELZ(N)**2 
ENDDO

      KIN_ENER = 0.5D0*ENER2

END SUBROUTINE KINETIC


SUBROUTINE FORCES
! SUBROUTINE TO CALCULATE THE FORCES DEPENDING ON THE INTERACTIONS DESIRED
! PEDRO OJEDA,  07/MAY/2011

USE PARAMETERS 


	  CALL INITIALIZE_GRAD

          CALL KINETIC

	  !CALL DELTA_MATRIX 

	  CALL SPRING

          CALL LENNARD_JONES


END SUBROUTINE FORCES




!SUBROUTINE CENTER_MASS(CMX,CMY,CMZ)
!
!! SUBROUTINE TO COMPUTE THE CENTER OF MASS FOR "NUM_RES" PARTICLES
!
!! PEDRO OJEDA,  07/MAY/2011
!
!USE PARAMETERS                                  
!             
!REAL(DP)            :: XCM,YCM,ZCM,CMX,CMY,CMZ,DIFF
!
!
!XCM=0.0D0
!
!YCM=0.0D0
!
!ZCM=0.0D0
!
!DO I=1,NUM_RES            
!                                                     
!    XCM= XCM + PARTICLE(I)%COOR(1)               
!                        
!    YCM= YCM + PARTICLE(I)%COOR(2)                
!                       
!    ZCM= ZCM + PARTICLE(I)%COOR(3)                 
!                      
!ENDDO
!
!CMX=XCM/(1.0*NUM_RES)
!CMY=YCM/(1.0*NUM_RES)
!CMZ=ZCM/(1.0*NUM_RES)
!
!IF( ABS(CMX) > 20.0D0 ) THEN
!
!        DO I=1,NUM_RES                  
!                                               
!            PARTICLE(I)%COOR(1)=PARTICLE(I)%COOR(1)-CMX             
!                          
!        ENDDO
!
!        I=CMX
!
!        J=MOD(I,20)
!
!        IF(CMX > 0) THEN
!
!            CMX= -20.0D0+ABS(1.0D0*J)
!
!        ELSE
!
!            CMX=  20.0D0-ABS(1.0D0*J) 
!
!        ENDIF
!
!        DO I=1,NUM_RES                                     
!                            
!            PARTICLE(I)%COOR(1)=PARTICLE(I)%COOR(1)+CMX                       
!                
!        ENDDO
!
!ENDIF
!
!IF( ABS(CMY) > 20.0D0 ) THEN
!
!        DO I=1,NUM_RES                        
!                                         
!            PARTICLE(I)%COOR(2)=PARTICLE(I)%COOR(2)-CMY            
!                           
!        ENDDO
!
!        I=CMY
!
!        J=MOD(I,20)
!
!        IF(CMY > 0) THEN
!
!            CMY= -20.0+ABS(1.0*J)
!
!        ELSE
!
!            CMY=  20.0-ABS(1.0*J) 
!
!        ENDIF
!
!        DO I=1,NUM_RES                                      
!                           
!            PARTICLE(I)%COOR(2)=PARTICLE(I)%COOR(2)+CMY                         
!              
!        ENDDO
!
!ENDIF
!
!IF( ABS(CMZ) > 20.0 ) THEN
!
!        DO I=1,NUM_RES                            
!                                     
!            PARTICLE(I)%COOR(3)=PARTICLE(I)%COOR(3)-CMZ                  
!                     
!        ENDDO
!
!        I=CMZ
!
!        J=MOD(I,20)
!
!        IF(CMZ > 0) THEN
!
!            CMZ= -20.0+ABS(1.0*J)
!
!        ELSE
!
!            CMZ=  20.0-ABS(1.0*J) 
!
!        ENDIF
!
!        DO I=1,NUM_RES                                 
!                                
!            PARTICLE(I)%COOR(3)=PARTICLE(I)%COOR(3)+CMZ                         
!              
!        ENDDO
!
!ENDIF
!
!
!END SUBROUTINE CENTER_MASS
