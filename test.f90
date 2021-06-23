PROGRAM test ! Solve A*X + B = 0, with A upper triangular
    implicit none

	CHARACTER(LEN=*),PARAMETER :: &
	    FMT_SCMPLX = '(*( "(" SP F8.3 X SP F8.3 " i)" 2X ))'  , & ! Complex nbr format as float 
	    FMT_LCMPLX = '(*( "(" SP E10.3 X SP E10.3 " i)" 2X ))'    ! Complex nbr format as power

    COMPLEX(8),PARAMETER :: &
        z  = (0.,0.)

    COMPLEX(8),PARAMETER :: &
        a1 = (1.,0.),       &
        a2 = (2.,0.),       &
        a3 = (3.,0.),       &
        a4 = (4.,0.)

    COMPLEX(8),PARAMETER :: &
        b1 = (5.,0.),       &
        b2 = (6.,0.),       &
        b3 = (7.,0.),       &
        b4 = (8.,0.)

    COMPLEX(8) :: A(4,4) = RESHAPE( &
            (/ a4,a3,a2,a1, &
               z ,a4,a3,a2, &
               z ,z ,a4,a3, &
               z ,z ,z ,a4  &
            /),(/4,4/),ORDER=(/2,1/) )

    COMPLEX(8) :: B(4,4) = RESHAPE( &
            (/ b4,b3,b2,b1, &
               b1,b4,b3,b2, &
               b2,b1,b4,b3, &
               b3,b2,b1,b4  &
            /),(/4,4/),ORDER=(/2,1/) )

    COMPLEX(8),ALLOCATABLE :: X(:,:)
    INTEGER                :: i,j
!
! -- PROGRAM STARTS HERE --
!
    CALL solve(A,B,X)

    CALL CMAT(B,msg='B matrix:')
    CALL CMAT(X,msg='X matrix:')

    B = -MATMUL(A,X)
    CALL CMAT(B,msg='Check X, reform B:')

    DEALLOCATE(X)
!
! -- Internal routines:
!
    contains
    SUBROUTINE solve(A,B,X)
        implicit none

        COMPLEX(8),            INTENT(IN)  :: A(:,:),B(:,:)
        COMPLEX(8),ALLOCATABLE,INTENT(OUT) :: X(:,:)

        INTEGER :: M,N,P
        INTEGER :: i,j,k

! -- Cannot initialise them because A,B shapes are still unknown at compilation time.
! -- But deferred-shape arrays' shapes become known at runtime, so assignment is possible instead.
        M = UBOUND(A,1)
        N = UBOUND(A,2)
        P = UBOUND(B,2)

        IF (UBOUND(B,1) .NE. N) ERROR STOP '*** WRONG A OR B MATRIX DIMENSIONS; EXITING. ***'

        ALLOCATE(X(M,P))
        X = (0.,0.)
! Do systematic recursive solve:
        DO i=N,1,-1
            DO j=1,P
                X(i,j) = ( -B(i,j) - DOT_PRODUCT(A(i,:),X(:,j))) / A(i,i)
            END DO
        END DO

    END SUBROUTINE solve

    SUBROUTINE CMAT(M,msg)
        implicit none

        COMPLEX(8),               INTENT(IN) :: M(:,:)
        CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: msg

        INTEGER :: i,j

        IF (PRESENT(msg)) WRITE(*,'(/,A)') msg
    	DO i = 1,UBOUND(M,1)
	    	WRITE(*,FMT_SCMPLX) (M(i,j), j=1,UBOUND(M,2))
    	END DO
    END SUBROUTINE CMAT

END PROGRAM test