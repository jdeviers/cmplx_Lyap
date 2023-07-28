PROGRAM cmplx_lyap
	USE mod_prec
	USE mod_init
	USE mod_proc
	USE mod_lpck
	implicit none

!	.. Work scalars ..
	INTEGER                :: i,j
	REAL(dp)               :: mu
	COMPLEX(8)             :: mu_c
!	.. Param scalars ..
	INTEGER,PARAMETER      :: N = 12
!	.. Cmplx arrays ..
	COMPLEX(8),ALLOCATABLE :: A(:,:),D(:,:),U(:,:),btmp(:,:)
	COMPLEX(8),ALLOCATABLE :: Q(:,:),R(:,:),B(:,:),CID(:,:)
	COMPLEX(8),ALLOCATABLE :: TAU(:),u_vec(:,:),b_rowvec(:,:),b_colvec(:,:)


! ! -- This block to use random initial A matrix
! 	ALLOCATE(A(N,N))
! 	WRITE(11,*) A ! Necessary to initialise A, it seems
! 	CALL INIT_A_RND(N,A)

	CALL INIT_A(A)
	WRITE(*,'(/,A)') 'Initial A matrix:'
	CALL PRINT_CMAT(A)

! -- This is really expensive:
	CALL SCHUR_DECOMP(A,TAU)
	WRITE(*,'(/,A)') 'Compact QR matrix:'
	CALL PRINT_CMAT(A)

	CALL UNPACK_QR(A,Q,R)

	WRITE(*,'(/,A)') 'Q before rebuilding; for now, just a concatenation of columns vectors v.'
	WRITE(*,'(A)')   'H_i = I - TAU(i) * (v * v*), and Q = H_1 * H_2 * ... * H_i * ... * H_m, with m the nbr of columns in Q.'
	CALL PRINT_CMAT(Q)

	CALL REBUILD_Q(Q,TAU)
	WRITE(*,'(/,A)') 'After rebuilding: Q is a unitary matrix (QQ* = Q*Q = I, where Q* is the conjugate transpose):'
	CALL PRINT_CMAT(Q)
	WRITE(*,'(/,A)') 'And R is an upper triangular matrix containing the eigenvalues of A on its diagonal:'
	CALL PRINT_CMAT(R)

	CALL CHECK_UNITARITY(Q)
! -- OK

!	IF (UBOUND(A,1).EQ.UBOUND(A,2)) CALL CHECK_REBUILD_A(Q,R) ! Need square matrix

	CALL INIT_D(D)
	WRITE(*,'(/,A)') 'D matrix, from AX + XA^H = D; must take the negative b/c we solve here AX + XA^H + (-D) = 0:'
	CALL PRINT_CMAT(D)

	D = -D
	CALL CHOLESKY(N,mat_to_packed(N,D),B)
	WRITE(*,'(/,A)') 'Cholesky factorisation B of D such that D = B * B^T:'
	CALL PRINT_CMAT(D)

! -- Modified Hammarling algorithm starts here:
	CALL INIT_U(N,U)
	ALLOCATE(b_rowvec(1,N))
	ALLOCATE(b_colvec(N,1))

	DO j=N,2,-1
! -- (1)
		b_rowvec(1,:) = B(j,:)
		mu            = NORM2(REALPART(b_rowvec))
		mu_c = SQRT( -2 * R(j,j) ) ! has to be complex b/c R(j,j) can be positive
! -- (2)
		ALLOCATE(u_vec(j-1,1))
		IF (mu .GT. 0.d0) THEN
			PRINT*, 'in'
			b_rowvec = b_rowvec / mu
			CALL BUILD_CID(j-1,CID)
! -- btmp construction block --
			ALLOCATE(btmp(j-1,1))
			b_colvec = TRANSPOSE(CONJG(b_rowvec))
			btmp = ( mu_c * MATMUL( B(1:j-1,:),b_colvec ) ) + ( R(1:j-1,j:j) * mu / mu_c) ! btmp is of dim (j-1,1)
! -----------------------------
			CALL SOLVE_FOR_UVEC(R(1:j-1,1:j-1),CID(:,:),btmp(:,:),u_vec(:,:))
			B(1:j-1,:) = B(1:j-1,:) - ( MATMUL( u_vec,b_rowvec ) * mu_c )
			DEALLOCATE(CID,btmp)
		ELSE
			u_vec = (0.d0,0.d0)
		END IF
! -- (3)
		U(j,j) = mu / mu_c
! -- (4)
		U(1:j-1,j:j) = u_vec
		DEALLOCATE(u_vec)
	END DO

	U(1,1) = NORM2(REALPART(B(1,:)),1) / SQRT( -2 * R(1,1) )
	WRITE(*,'(/,A)') 'Cholesky factorisation U of the solution P such that P = U * U^H:'
	CALL PRINT_CMAT(U)

	DEALLOCATE(A,D)   ! Initial matrices
	DEALLOCATE(Q,R)   ! Algo cmplx matrices
	DEALLOCATE(b_rowvec,b_colvec) ! Algo cmplx vectors
	DEALLOCATE(U,B)   ! Algo real matrices

END PROGRAM cmplx_lyap
