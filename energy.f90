SUBROUTINE DELTA_MATRIX 

! SUBROUTINE TO COMPUTE THE DISTANCES AND VECTORS BETWEEN "NUM_RES" PARTICLES
! PEDRO OJEDA,  07/MAY/2011
USE PARAMETERS                                               
REAL(DP)            ::  FDUMMY1,FDUMMY2,FDUMMY3

DO I=1,NUM_RES  
  DO J=I,NUM_RES
        DIS(I,J)=0.0D0   
        DIS(J,I)=0.0D0 
  ENDDO                
ENDDO                  

!THE DISTANCES BETWEEN AMINO ACIDS ARE COMPUTED HERE 
DO 10 I=1,NUM_RES-1                                    
  DO 20 J=I+1,NUM_RES                                     
    FDUMMY1=PARTICLE(I)%COORX(1)-PARTICLE(J)%COORX(1)                          
    FDUMMY2=PARTICLE(I)%COORY(1)-PARTICLE(J)%COORY(1)                          
    FDUMMY3=PARTICLE(I)%COORZ(1)-PARTICLE(J)%COORZ(1)                           
            

    XDIS(I,J)=FDUMMY1
    YDIS(I,J)=FDUMMY2
    ZDIS(I,J)=FDUMMY3
    XDIS(J,I)=-FDUMMY1
    YDIS(J,I)=-FDUMMY2
    ZDIS(J,I)=-FDUMMY3
    DIS(I,J)=SQRT(FDUMMY1*FDUMMY1+FDUMMY2*FDUMMY2+FDUMMY3*FDUMMY3)
    DIS(J,I)=DIS(I,J)

20  CONTINUE

10  CONTINUE

END SUBROUTINE DELTA_MATRIX






SUBROUTINE LENNARD_JONES 
! SUBROUTINE TO COMPUTE THE LENNARD-JONES POTENTIAL AND FORCE
! PEDRO OJEDA,  07/MAY/2011

USE PARAMETERS
REAL(DP)            :: FDUMMY1,FDUMMY2,FDUMMY3,FDUMMY4,FDUMMY5
REAL(DP)            :: FDUMMY0,FDUMMYD,FDUMMYE
REAL(DP)            :: FTERM
REAL(DP)            :: ENER1


	ENER1=0.0D0

!VAN DER WAALS CONTRIBUTION
DO I=1,NUM_RES-2
     DO J=I+2,NUM_RES
           FDUMMY1=SIGMA_CONST(I,J)/DIS(I,J)
           K=CLASE(I)
           L=CLASE(J)
	   FDUMMY2=FDUMMY1*FDUMMY1                             ! POWER 2
	   FDUMMY3=FDUMMY2*FDUMMY2                             ! POWER 4
	   FDUMMY4=FDUMMY2*FDUMMY3                             ! POWER 6
	   FDUMMY5=FDUMMY4*FDUMMY4                             ! POWER 12
	   FDUMMY0=1.0/DIS(I,J)                                 ! 1/DISTANCE_IJ
	   FDUMMYD=FDUMMY0*FDUMMY0                             ! 1/(DISTANCE_IJ)^2
	   FDUMMYE=FDUMMY5-FDUMMY4                             ! (SIGMA/R)^12 - (SIGMA/R)^6
	   ENER1=ENER1+EPS_CONST(K,L)*FDUMMYE                  ! ENERGY = EPSILON*[(SIGMA/R)^12 - (SIGMA/R)^6 ]
	   FTERM=(12.*FDUMMY5-6.*FDUMMY4)*FDUMMYD              ! FORCE = -DV/DR = [12*(SIGMA/R)^12 - 6*(SIGMA/R)^6 ]*(1/R^2)

           PARTICLE(I)%GRADX(1)=PARTICLE(I)%GRADX(1)+EPS_CONST(K,L)*FTERM*XDIS(I,J)   ! X COMPONENT OF THE GRADIENT
           PARTICLE(I)%GRADY(1)=PARTICLE(I)%GRADY(1)+EPS_CONST(K,L)*FTERM*YDIS(I,J)   ! Y COMPONENT OF THE GRADIENT
           PARTICLE(I)%GRADZ(1)=PARTICLE(I)%GRADZ(1)+EPS_CONST(K,L)*FTERM*ZDIS(I,J)   ! Z COMPONENT OF THE GRADIENT
           PARTICLE(J)%GRADX(1)=PARTICLE(J)%GRADX(1)-EPS_CONST(K,L)*FTERM*XDIS(I,J)   ! X COMPONENT OF THE GRADIENT
           PARTICLE(J)%GRADY(1)=PARTICLE(J)%GRADY(1)-EPS_CONST(K,L)*FTERM*YDIS(I,J)   ! Y COMPONENT OF THE GRADIENT
           PARTICLE(J)%GRADZ(1)=PARTICLE(J)%GRADZ(1)-EPS_CONST(K,L)*FTERM*ZDIS(I,J)   ! Z COMPONENT OF THE GRADIENT

     ENDDO

ENDDO

DO I=1,NUM_RES
           write(6,*) PARTICLE(I)%GRADX(1),PARTICLE(I)%GRADY(1),PARTICLE(I)%GRADZ(1)
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



	R_CERO=3.8D0
	A_CONST=50.0D0
	ENER2=0.0D0


!HARMONIC POTENTIAL
FDUMMY3=2.*A_CONST

DO I=2,NUM_RES-1
      FDUMMY1=(DIS(I,I+1)-R_CERO)
      FDUMMY2=FDUMMY1*FDUMMY1
      ENER2=ENER2+FDUMMY2
      FDUMMY4=(DIS(I,I-1)-R_CERO)
      PARTICLE(I)%GRADX(1)=PARTICLE(I)%GRADX(1)-FDUMMY3*XDIS(I,I+1)*FDUMMY1/DIS(I,I+1) - FDUMMY3*XDIS(I,I-1)*FDUMMY4/DIS(I-1,I)
      PARTICLE(I)%GRADY(1)=PARTICLE(I)%GRADY(1)-FDUMMY3*YDIS(I,I+1)*FDUMMY1/DIS(I,I+1) - FDUMMY3*YDIS(I,I-1)*FDUMMY4/DIS(I-1,I)
      PARTICLE(I)%GRADZ(1)=PARTICLE(I)%GRADZ(1)-FDUMMY3*ZDIS(I,I+1)*FDUMMY1/DIS(I,I+1) - FDUMMY3*ZDIS(I,I-1)*FDUMMY4/DIS(I-1,I)
ENDDO


      FDUMMY1=DIS(1,2)-R_CERO
      FDUMMY2=FDUMMY1*FDUMMY1
      ENER2=ENER2+FDUMMY2
      PARTICLE(1)%GRADX(1)=PARTICLE(1)%GRADX(1)-FDUMMY3*XDIS(1,2)*FDUMMY1/DIS(1,2)
      PARTICLE(1)%GRADY(1)=PARTICLE(1)%GRADY(1)-FDUMMY3*YDIS(1,2)*FDUMMY1/DIS(1,2)
      PARTICLE(1)%GRADZ(1)=PARTICLE(1)%GRADZ(1)-FDUMMY3*ZDIS(1,2)*FDUMMY1/DIS(1,2)
      FDUMMY1=DIS(NUM_RES,NUM_RES-1)-R_CERO
      PARTICLE(NUM_RES)%GRADX(1)=PARTICLE(NUM_RES)%GRADX(1)-FDUMMY3*XDIS(NUM_RES,NUM_RES-1)*FDUMMY1/DIS(NUM_RES,NUM_RES-1)
      PARTICLE(NUM_RES)%GRADY(1)=PARTICLE(NUM_RES)%GRADY(1)-FDUMMY3*YDIS(NUM_RES,NUM_RES-1)*FDUMMY1/DIS(NUM_RES,NUM_RES-1)
      PARTICLE(NUM_RES)%GRADZ(1)=PARTICLE(NUM_RES)%GRADZ(1)-FDUMMY3*ZDIS(NUM_RES,NUM_RES-1)*FDUMMY1/DIS(NUM_RES,NUM_RES-1)

      POT_ENER = POT_ENER + A_CONST*ENER2

END SUBROUTINE SPRING




SUBROUTINE FORCES
! SUBROUTINE TO CALCULATE THE FORCES DEPENDING ON THE INTERACTIONS DESIRED
! PEDRO OJEDA,  07/MAY/2011

USE PARAMETERS 


	  CALL INITIALIZE_GRAD

	  CALL DELTA_MATRIX 

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
