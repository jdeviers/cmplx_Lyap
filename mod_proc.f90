MODULE mod_proc
	USE mod_prec
	implicit none

	contains

	SUBROUTINE PRINT_CMAT(M)
		implicit none

!	.. Parameters ..
		COMPLEX(8),INTENT(IN) :: M(:,:)

!	.. Local scalars ..
		INTEGER               :: i,j

		DO i = 1,UBOUND(M,1)
			WRITE(*,FMT_SCMPLX) (M(i,j), j=1,UBOUND(M,2))
		END DO

	END SUBROUTINE PRINT_CMAT
! ----------
	SUBROUTINE PRINT_CVEC(V)
		implicit none

!	.. Parameters ..
		COMPLEX(8),INTENT(IN) :: V(:)

!	.. Local scalars ..
		INTEGER               :: i

		DO i = 1,UBOUND(V,1)
			WRITE(*,FMT_SCMPLX) V(i)
		END DO

	END SUBROUTINE PRINT_CVEC
! ----------
	SUBROUTINE UNPACK_QR(QR,Q,R)
		implicit none

!	.. Parameters ..
		COMPLEX(8),ALLOCATABLE,INTENT(INOUT) :: QR(:,:)
		COMPLEX(8),ALLOCATABLE,INTENT(OUT)   :: Q(:,:),R(:,:)

!	.. Local scalars ..
		INTEGER                              :: i,j
		INTEGER                              :: M,N,K

		M = UBOUND(QR,1);N = UBOUND(QR,2);K = MIN(M,N)
		ALLOCATE( Q(K,K), R(K,N) )
		Q = (0.,0.); R = (0.,0.)
		DO i = 1,M
			Q(i,i) = (1.,0.)
			DO j = 1,N
				IF (j.GE.i) THEN
					R(i,j) = QR(i,j)
				ELSE
					Q(i,j) = QR(i,j)
				END IF
			END DO
		END DO

	END SUBROUTINE UNPACK_QR
! ----------
	SUBROUTINE BUILD_CID(N,ID)
		implicit none

!	.. Parameters ..
		INTEGER,               INTENT(IN)  :: N
		COMPLEX(8),ALLOCATABLE,INTENT(OUT) :: ID(:,:)

!	.. Local scalars .. 
		INTEGER :: i

		ALLOCATE(ID(N,N))
		ID = (0.,0.)
		DO i = 1,N
			ID(i,i) = (1.,0.)
		END DO

	END SUBROUTINE BUILD_CID
! ----------
	SUBROUTINE REBUILD_Q(Q,TAU)
		implicit none

!	.. Parameters ..
		COMPLEX(8),INTENT(INOUT) :: Q(:,:)
		COMPLEX(8),INTENT(IN)    :: TAU(:)

!	.. Local scalars ..
		INTEGER :: M,i

!	.. Local cmplx ..
		COMPLEX(8),ALLOCATABLE   :: ID(:,:)
		COMPLEX(8),ALLOCATABLE   :: H_i(:,:),H_tmp(:,:)
		COMPLEX(8),ALLOCATABLE   :: v(:,:),v_t(:,:)

		M = UBOUND(Q,1)
		CALL BUILD_CID(M,ID)

		ALLOCATE(H_i(M,M),H_tmp(M,M))
		ALLOCATE(v(M,1),v_t(1,M))

		H_tmp = ID
		DO i=1,M 
			v(:,1)   =       Q(1:M,i)
			v_t(1,:) = CONJG(Q(1:M,i))

			! WRITE(*,'(/,A)') 'v = '
			! WRITE(*,FMT_CMPLX) (v(j,1), j=1,M)
			! WRITE(*,'(/,A)') 'v_t = '
			! WRITE(*,FMT_CMPLX) (v_t(1,j), j=1,M)

			H_i = ID - ( TAU(i) * MATMUL(v,v_t) )
			H_tmp = MATMUL(H_tmp,H_i)
		END DO
		Q = H_tmp

		DEALLOCATE(H_i,H_tmp)
		DEALLOCATE(v,v_t)
	END SUBROUTINE REBUILD_Q
! ----------
	SUBROUTINE CHECK_UNITARITY(Q)
		implicit none

!	.. Parameters ..
		COMPLEX(8),INTENT(IN)  :: Q(:,:)

!	.. Local cmplx ..
		COMPLEX(8),ALLOCATABLE :: Q_t(:,:)

		WRITE(*,'(/,A)') 'Check unitarity: should be an identity matrix:'

		ALLOCATE( Q_t(UBOUND(Q,1),UBOUND(Q,2)) )
		Q_t = CONJG(Q)

		CALL PRINT_CMAT( MATMUL(Q,TRANSPOSE(Q_t)) )
		DEALLOCATE(Q_t)

	END SUBROUTINE CHECK_UNITARITY
! ----------
	SUBROUTINE CHECK_REBUILD_A(Q,R)
		implicit none

!	.. Parameters ..
		COMPLEX(8),ALLOCATABLE,INTENT(IN) :: Q(:,:),R(:,:)

!	.. Local cmplx ..
		COMPLEX(8),ALLOCATABLE            :: Q_t(:,:)

		ALLOCATE( Q_t(UBOUND(Q,1),UBOUND(Q,2)) )
		Q_t = CONJG(Q)

		WRITE(*,'(/,A)') 'Checking QR worked: rebuild A = Q * R * Q^H'
		CALL PRINT_CMAT( MATMUL(Q,MATMUL(R,TRANSPOSE(Q_t))) )
		DEALLOCATE(Q_t)

	END SUBROUTINE CHECK_REBUILD_A

END MODULE mod_proc