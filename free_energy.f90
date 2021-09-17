SUBROUTINE INITIALIZE_WANG_LANDAU_MD

! SUBROUTINE TO UPDATE THE MICROCANONICAL TEMPERATURE BY USING THE SCHEME OF

! KIM ET. AL., PRL 97, 050601

! PEDRO OJEDA,  07/MAY/2011

USE PARAMETERS

USE WL_MODULE
  
USE Ziggurat


FACT_F   = 1.0

SCALE_WL = 1.0        

DF = FACT_F / ( 2.0 * DELTA_ENER )

OMEGA = 5.0

ACCEPT = .TRUE.

INDICE = 10000

TEMP0 = 15.0

INTERV_SIZE = 550

DO 10 I = E_MIN , E_MAX 

           ALPHA( I )  = 1.0

           T_HIST( I ) = 1.0  !* ( 1.0*I - 1.0*E_MAX ) / (1.0*E_MAX - 1.0*E_MIN) + 1.0 

           DOS( I ) = 0.0   !( 1.0*I - 1.0*E_MAX ) / (1.0*E_MAX - 1.0*E_MIN) + 1.0 

           CUM_HIST( I ) = 1.0

           VISITS( I ) = 1.0

10  CONTINUE



DO 20 I = -40, 40
           
           GAUSSIAN(I) = EXP( -(0.1*I)*(0.1*I) )

20    CONTINUE



END SUBROUTINE INITIALIZE_WANG_LANDAU_MD





SUBROUTINE WANG_LANDAU_MD

! SUBROUTINE TO UPDATE THE MICROCANONICAL TEMPERATURE BY USING THE SCHEME OF

! KIM ET. AL., PRL 97, 050601

! PEDRO OJEDA,  07/MAY/2011

USE PARAMETERS 

USE Ziggurat

USE WL_MODULE
          
INTEGER                :: P, Q, LIMIT, DIFF

REAL                   :: R, ETA, PENDIENTE, SCALING
         
REAL                   :: XMEAN, YMEAN, XMEANSQ, YMEANSQ, XYMEAN, NUMBER_BINS, ACOEF, BCOEF , DENOMIN

REAL                   :: OLD_B, DELT_DIFF_INT

IF(POT_ENER < -100.0 .AND. ACCEPT .EQV. .TRUE. ) THEN

    E_BIN_OLD = ANINT( POT_ENER / DELTA_ENER )

    ACCEPT =.FALSE.

    WRITE(6,*) 'BIN INICIAL', E_BIN_OLD

ENDIF


IF(MOD(NSTEP,10000000) .EQ. 0 ) THEN

 FACT_F = 0.5 * FACT_F

 WRITE(6,*) 'UPDATE FACTOR F', FACT_F

DO 15 I = E_MIN, E_MAX
          T_HIST( I ) = DOS( I )
15     CONTINUE 


IF(FACT_F > 0.001) THEN
!IF(INTERV_SIZE > 10) THEN
!DO P = 1,33
        OLD_B =  0.0
P = 1
60      CONTINUE
        
        XMEAN = 0.0
        YMEAN = 0.0
        XMEANSQ = 0.0
        YMEANSQ = 0.0 
        XYMEAN = 0.0
        NUMBER_BINS = 0.0
        DO Q = E_MAX - INTERV_SIZE * (P - 1) , E_MAX - INTERV_SIZE * P , -1
             XMEAN = XMEAN + 1.0*Q
             YMEAN = YMEAN + DOS(Q)
             XMEANSQ = XMEANSQ + 1.0*Q*Q
             YMEANSQ = YMEANSQ + DOS(Q)*DOS(Q) 
             XYMEAN = XYMEAN + 1.0*Q*DOS(Q)
             NUMBER_BINS = NUMBER_BINS + 1.0
        ENDDO 



        XMEAN = XMEAN / NUMBER_BINS
        YMEAN = YMEAN / NUMBER_BINS
        DENOMIN = XMEANSQ - NUMBER_BINS * XMEAN * XMEAN

        ACOEF = ( YMEAN * XMEANSQ - XMEAN * XYMEAN ) / DENOMIN
        BCOEF =  ( XYMEAN - NUMBER_BINS * XMEAN * YMEAN ) / DENOMIN

        WRITE(6,*) ACOEF, BCOEF

        IF( ABS(OLD_B) > 0.0 ) THEN
              DELT_DIFF_INT = OLD_B - 1.0*( ACOEF + BCOEF * (1.0 * E_MAX - INTERV_SIZE * (P - 1)) )
        ENDIF


        DO Q = E_MAX - INTERV_SIZE * (P - 1) , E_MAX - INTERV_SIZE * P , -1
              IF( ABS(OLD_B) > 0.0 ) THEN
   
                    DOS(Q) = ACOEF + BCOEF * 1.0 * Q + DELT_DIFF_INT

               ELSE
                    DOS(Q) = ACOEF + BCOEF * 1.0 * Q                  
!             WRITE(6,*) Q, T_HIST(Q)
              ENDIF
        ENDDO 

        OLD_B = DOS( E_MAX - INTERV_SIZE * P ) 


        IF(Q  > -16500 ) THEN
            P = P+1
            GOTO 60 
        ELSE
            GOTO 70
        ENDIF

70      CONTINUE


        IF(INTERV_SIZE >50 ) THEN
             INTERV_SIZE = INTERV_SIZE - 50
!        ELSE
!             INTERV_SIZE = INTERV_SIZE - 5
        ENDIF

 WRITE(6,*) 'INTERVAL SAMPLE', INTERV_SIZE, 'FACT_F', FACT_F

!ENDIF
ENDIF



ENDIF

E_BIN = ANINT( POT_ENER / DELTA_ENER )


!WRITE(6,*) NSTEP, POT_ENER 

IF(E_BIN  .LE. E_MAX .AND. E_BIN .GE. E_MIN) THEN

              DO 10 P = -40, 40
                  IF( E_BIN + P .LE. E_MAX )  THEN   
                     DOS( E_BIN + P ) =  DOS( E_BIN + P ) + FACT_F * GAUSSIAN( P )
                  ENDIF
10         CONTINUE

!        DOS(E_BIN) = DOS(E_BIN) + FACT_F

        IF(NSTEP > 100000) THEN

!               T_HIST(E_BIN) = (DOS(E_BIN + 100) - DOS(E_BIN))/(100.0*DELTA_ENER)

              

!                 WRITE(6,*) NSTEP, E_BIN, T_HIST(E_BIN), 'TRUE'
       ! ELSE

            ! T_HIST(E_BIN) = T_HIST(-100)

             


!              IF( MOD(NSTEP,1000) .EQ. 0) THEN   

          !        WRITE(6,*) 'BIN DE PRUEBA', E_BIN

!                 IF(ABS(T_HIST(E_BIN)) > 0.01) THEN

          !         T_HIST(E_BIN) = (DOS(E_BIN) - DOS(E_BIN_OLD)) / ( (1.0*E_BIN-1.0*E_BIN_OLD)*DELTA_ENER )


               LIMIT = -11000 + 500 * UNI()

               DIFF = AINT( 300.0*FACT_F*UNI() + 3 )

               IF( E_BIN > LIMIT ) THEN

                   PENDIENTE = (DOS(E_BIN) - DOS(E_BIN - DIFF)) / ( (1.0*DIFF)*DELTA_ENER )  ! ESTO FUNCIONA BIEN

!                   PENDIENTE = (DOS(E_BIN) - DOS(E_BIN - 1)) / ( DELTA_ENER )
           
                   SCALE_WL = PENDIENTE
               ELSE

                   PENDIENTE = (DOS(E_BIN + DIFF) - DOS(E_BIN)) / ( (1.0*DIFF)*DELTA_ENER )  ! ESTO FUNCIONA BIEN

!                   PENDIENTE = (DOS(E_BIN) - DOS(E_BIN - 1)) / ( DELTA_ENER )
           
                   SCALE_WL = PENDIENTE

               ENDIF        


          !         SCALE_WL = T_HIST(E_BIN) + PENDIENTE * (POT_ENER - E_BIN * DELTA_ENER)

          !         WRITE(6,*) NSTEP, SCALE_WL, POT_ENER ! E_BIN, T_HIST(E_BIN)
         !          WRITE(6,*) 'BINS',DOS(E_BIN), DOS(E_BIN_OLD) , E_BIN, E_BIN_OLD
!                   E_BIN_OLD = E_BIN

!                 ENDIF
!              ENDIF

        ELSE

                  SCALE_WL = 1.0

        ENDIF 

ELSE

     SCALE_WL = 1.0

ENDIF

               
         CALL RESCALE_FORCES   
     

!IF( MOD(NSTEP,INDICE) .EQ. 0 ) THEN


!        IF(ACCEPT) THEN

!               E_BIN_OLD = ANINT( POT_ENER / DELTA_ENER )

!               CALL COPY_COOR

!               CALL COPY_VEL

!               INDICE = 200

!               ACCEPT = .FALSE.

!        ENDIF

!       WRITE(6,*) E_BIN_OLD



!        E_BIN_NEW = ANINT( POT_ENER / DELTA_ENER )

!        IF(E_BIN_NEW .LE. E_MAX) THEN

!        ETA= DOS(E_BIN_NEW) - DOS(E_BIN_OLD)

!        IF(ETA.LE.0.0) THEN    ! COMPARACION DE R

!          DOS(E_BIN_NEW)=DOS(E_BIN_NEW)+FACT_F

!	  DO 10 P=1,NUM_RES

!	      PARTICLE_BACKUP(P)%COOR(1) = PARTICLE(P)%COOR(1) 

!	      PARTICLE_BACKUP(P)%COOR(2) = PARTICLE(P)%COOR(2) 

!	      PARTICLE_BACKUP(P)%COOR(3) = PARTICLE(P)%COOR(3) 

!	      PARTICLE_BACKUP(P)%VEL(1) = PARTICLE(P)%VEL(1) 

!	      PARTICLE_BACKUP(P)%VEL(2) = PARTICLE(P)%VEL(2) 

!	      PARTICLE_BACKUP(P)%VEL(3) = PARTICLE(P)%VEL(3)

	      
!10       CONTINUE

!              E_BIN_OLD = E_BIN_NEW

!        ELSE                 ! COMPARACION DE R

!        R=UNI()

!        IF(R.LT.EXP(-ETA)) THEN  !SEGUNDO

!          DOS(E_BIN_NEW)=DOS(E_BIN_NEW)+FACT_F

!	  DO 20 P=1,NUM_RES

!	      PARTICLE_BACKUP(P)%COOR(1) = PARTICLE(P)%COOR(1) 

!	      PARTICLE_BACKUP(P)%COOR(2) = PARTICLE(P)%COOR(2) 

!	      PARTICLE_BACKUP(P)%COOR(3) = PARTICLE(P)%COOR(3) 

!	      PARTICLE_BACKUP(P)%VEL(1) = PARTICLE(P)%VEL(1) 

!	      PARTICLE_BACKUP(P)%VEL(2) = PARTICLE(P)%VEL(2) 

!	      PARTICLE_BACKUP(P)%VEL(3) = PARTICLE(P)%VEL(3)
	      
!20       CONTINUE

!              E_BIN_OLD = E_BIN_NEW

!        ELSE                     !SEGUNDO

!          DOS(E_BIN_OLD)=DOS(E_BIN_OLD)+FACT_F

!	  DO 30 P=1,NUM_RES

!	      PARTICLE(P)%COOR(1) = PARTICLE_BACKUP(P)%COOR(1) + UNI()

!	      PARTICLE(P)%COOR(2) = PARTICLE_BACKUP(P)%COOR(2) + UNI()

!	      PARTICLE(P)%COOR(3) = PARTICLE_BACKUP(P)%COOR(3) + UNI()

!	      PARTICLE(P)%VEL(1) = PARTICLE_BACKUP(P)%VEL(1) 

!	      PARTICLE(P)%VEL(2) = PARTICLE_BACKUP(P)%VEL(2) 

!	      PARTICLE(P)%VEL(3) = PARTICLE_BACKUP(P)%VEL(3)
	      
!30       CONTINUE

!        ENDIF                    !SEGUNDO

!        ENDIF                ! COMPARACION DE R


!      ENDIF

!ENDIF


!IF( (E_BIN) .LE. E_MAX ) THEN 


!              DO 10 P = -40, 40

!                     DOS( E_BIN + P ) =  DOS( E_BIN + P ) + GAUSSIAN( P )

!10         CONTINUE


!         DOS( E_BIN ) = DOS( E_BIN ) + FACT_F

!        IF( MOD(NSTEP, 100000) .EQ. 0 ) THEN

!              DO 10 P = E_MIN , E_MAX !- 400

                !   MIDDLE = ANINT ( ( 1.0*E_MAX - 1.0*E_MIN )/2.0 )

               !    IF (E_BIN .LT. MIDDLE ) THEN

               !        T_HIST( P ) = ( DOS( P + 400 ) - DOS( P ) ) / ( 400.0 * DELTA_ENER )

               !    ELSE

              !  P = E_BIN     

              !         T_HIST( P ) = ( DOS( P ) - DOS( P - 1 ) ) / ( DELTA_ENER ) 

              !     ENDIF

!10     CONTINUE                
       
!        ENDIF
!
!	ALPHA( E_BIN + 1 )  = 1.0 / ( 1.0 - DF * T_HIST( E_BIN + 1 ) )
!
!	ALPHA( E_BIN - 1 )  = 1.0 / ( 1.0 + DF * T_HIST( E_BIN - 1 ) )
!
!	T_HIST( E_BIN + 1 ) = ALPHA( E_BIN + 1 ) * T_HIST( E_BIN + 1 )
!
!	T_HIST( E_BIN - 1 ) = ALPHA( E_BIN - 1 ) * T_HIST( E_BIN - 1 )
!

!	T_HIST( E_BIN ) = ( T_HIST( E_BIN + 1 ) +  T_HIST( E_BIN - 1 ) ) / 2.0


!           VISITS( E_BIN ) = VISITS( E_BIN ) + 1.0

!           CUM_HIST( E_BIN ) = CUM_HIST( E_BIN ) + T_HIST( E_BIN ) 


       ! IF(  T_HIST( E_BIN )  .LT. 0.0 ) THEN 

         !     SCALE_WL =  ABS( TEMP * T_HIST ( E_BIN ) )

        !     T_HIST( E_BIN ) = 0.1

         !     SCALE_WL = 1.0

        ! ENDIF

        !IF(  T_HIST( E_BIN )  .GT. 2.0 ) THEN 

         !     SCALE_WL =  ABS( TEMP * T_HIST ( E_BIN ) )

        !     T_HIST( E_BIN ) = 10.0

              

       ! ENDIF
         !    SCALE_WL = TEMP * T_HIST ( E_BIN )


!ENDIF

          

!              DO 20 P = E_MIN, E_MAX-2

!                     T_HIST( P ) = ( (1.0 - OMEGA) *DOS( P ) + (OMEGA/DELTA_ENER**4) * ( -DOS(P-2) &
!                     +4.0*DOS(P-1) + 4.0*DOS(P+1) -DOS(P+2) ) ) / (1.0-OMEGA + 6.0*OMEGA/DELTA_ENER**4) 

!20         CONTINUE


!
!              SCALE_WL = 1.0

!        CALL RESCALE_FORCES    



END SUBROUTINE WANG_LANDAU_MD




SUBROUTINE RESCALE_FORCES

! SUBROUTINE TO UPDATE THE MICROCANONICAL TEMPERATURE BY USING THE SCHEME OF

! KIM ET. AL., PRL 97, 050601

! PEDRO OJEDA,  07/MAY/2011

USE PARAMETERS

USE WL_MODULE
          
INTEGER             :: P



             DO 10 P=1,NUM_RES

	         PARTICLE(P)%GRAD(1) =  SCALE_WL * PARTICLE(P)%GRAD(1)

	         PARTICLE(P)%GRAD(2) =  SCALE_WL * PARTICLE(P)%GRAD(2)

	         PARTICLE(P)%GRAD(3) =  SCALE_WL * PARTICLE(P)%GRAD(3)

10           CONTINUE



END SUBROUTINE RESCALE_FORCES





SUBROUTINE WRITE_WANG_LANDAU_MD_HIST

! SUBROUTINE TO UPDATE THE MICROCANONICAL TEMPERATURE BY USING THE SCHEME OF

! KIM ET. AL., PRL 97, 050601

! PEDRO OJEDA,  07/MAY/2011

USE WL_MODULE

INTEGER            :: I 


DO 10 I = E_MIN, E_MAX
 
         WRITE(89,*) I, DOS( I )

10     CONTINUE

DO 20 I = E_MIN, E_MAX
 
         WRITE(90,*) I,  T_HIST( I )

20     CONTINUE

         WRITE(89,*) '&'

         WRITE(90,*) '&'

END SUBROUTINE WRITE_WANG_LANDAU_MD_HIST


