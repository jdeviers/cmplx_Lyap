PROGRAM cmplx_lyap
	USE mod_prec
	USE mod_init
	USE mod_proc
	USE mod_lpck
	implicit none

	INTEGER,PARAMETER      :: N = 12
	COMPLEX(8),ALLOCATABLE :: A(:,:),D(:,:)
	COMPLEX(8),ALLOCATABLE :: Q(:,:),R(:,:),B(:,:)
	COMPLEX(8),ALLOCATABLE :: TAU(:)

! ! -- This block to use random initial A matrix
! 	ALLOCATE(A(N,N))
! 	WRITE(11,*) A ! Necessary to initialise A, it seems
! 	CALL INIT_A_RND(N,A)

	CALL INIT_A(A)
	WRITE(*,'(/,A)') 'Initial A matrix:'
	CALL PRINT_CMAT(A)

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

	IF (UBOUND(A,1).EQ.UBOUND(A,2)) CALL CHECK_REBUILD_A(Q,R) ! Need square matrix
! -- This fails.

	CALL INIT_D(D)
	WRITE(*,'(/,A)') 'D matrix, from AX + XA^H = D; must take the negative b/c we solve here AX + XA^H + (-D) = 0:'
	CALL PRINT_CMAT(D)
	CALL CHOLESKY(N,mat_to_packed(N,D),B)


	DEALLOCATE(A,D)
	DEALLOCATE(Q,R,B)

END PROGRAM cmplx_lyap