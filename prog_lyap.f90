PROGRAM cmplx_lyap
	USE mod_prec
	USE mod_proc
	USE mod_lpck
	implicit none

	INTEGER,PARAMETER      :: N = 5
	COMPLEX(8),ALLOCATABLE :: A(:,:)
	COMPLEX(8),ALLOCATABLE :: Q(:,:),R(:,:)
	COMPLEX(8),ALLOCATABLE :: TAU(:)

	COMPLEX(8),ALLOCATABLE :: B(:,:)

	ALLOCATE(A(N,N))
	WRITE(11,*) A ! Necessary to initialise A, it seems
	CALL INIT_A(N,A)

	WRITE(*,'(/,A)') 'Initial A matrix:'
	CALL PRINT_CMAT(A)

!	TEST:
	ALLOCATE(B(3,4))
    B(1,1)=(-1.60,0.10);B(1,2)=( 0.30, 1.70);B(1,3)=( 0.30, 0.20);B(1,4)=(-0.50,-1.80)
	B(2,1)=(-1.20,0.00);B(2,2)=(-0.90,-0.50);B(2,3)=( 1.50, 0.80);B(2,4)=( 1.50,-1.10)
    B(3,1)=(-0.10,1.30);B(3,2)=(-1.10, 0.50);B(3,3)=( 0.40,-1.30);B(3,4)=( 1.60, 0.70)

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

	DEALLOCATE(A,B)

END PROGRAM cmplx_lyap