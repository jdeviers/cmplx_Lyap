MODULE mod_lpck
	USE mod_prec
	USE mod_proc	
	implicit none

	contains

	SUBROUTINE INIT_A(N,A)
		implicit none

!	.. Parameters ..
		COMPLEX(8),INTENT(OUT) :: A(:,:) 
		INTEGER,   INTENT(IN)  :: N

!	.. Local scalars ..
		INTEGER,PARAMETER      :: IDIST = 1
		INTEGER                :: i,ISEED(4)

		DO i=1,N
			CALL ZLARNV(IDIST,ISEED,N,A(:,i))
		END DO

	END SUBROUTINE INIT_A

	SUBROUTINE SCHUR_DECOMP(A,TAU)
		implicit none

!	.. Parameters ..
		COMPLEX(8),ALLOCATABLE,INTENT(INOUT) :: A(:,:)
		COMPLEX(8),ALLOCATABLE,INTENT(OUT)   :: TAU(:)

!	.. Local cmplx ..
		COMPLEX(8),ALLOCATABLE               :: WORK(:,:)

!	.. Local scalars ..
		INTEGER                              :: M,N,LDA,LWORK,INFO


		M=UBOUND(A,1);N=UBOUND(A,2);LDA=M
		ALLOCATE( WORK(1,1), TAU(MIN(M,N)) )

! -- Dry run to determine optimal WORK size
		CALL ZGEQRF(M,N,A,LDA,TAU,WORK,-1,INFO)

! -- Set LWORK and WORK, actual computation of QR
		LWORK=INT(WORK(1,1))
		DEALLOCATE(WORK);ALLOCATE(WORK(1,LWORK))
		CALL ZGEQRF(M,N,A,LDA,TAU,WORK,LWORK,INFO)
		DEALLOCATE(WORK)

		IF (INFO.NE.0) THEN
			WRITE(*,'(A)') ' *** QR FAILED *** '
			STOP
		END IF

	END SUBROUTINE SCHUR_DECOMP

END MODULE mod_lpck