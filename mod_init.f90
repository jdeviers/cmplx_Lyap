MODULE mod_init
	USE mod_prec
	implicit none

! -- AX + XA^H = Q

	contains

	SUBROUTINE init_A(A)
		implicit none

		COMPLEX(8),ALLOCATABLE,INTENT(OUT) :: A(:,:)
		REAL(8)                ::                    &
			a1  =   -1.0d0, a2 = 0.0d0, a3 = -0.5d0, &
			a4  =   -1.7097622964577763e-15,         &
			a5  =   -8.168140899333462     ,         &
			a6  =    2.2371143170757382e-17,         &
			a7  =  154.5663585566178       ,         &
			a8  = -148.34632244319613      ,         &
			a9  =    4.3982297150257095    ,         &
			a10 =    6.220036113421712     ,         &
			a11 =  160.78639467003956


		ALLOCATE(A(12,12))
!		A = (0.,0.)

!		A(1,1) =cmplx(a1,a4),    A(1,5)  =cmplx(a2,a5), A(1,8) =cmplx(a6,a7),                                             &
!		A(2,2) =cmplx(a1,a2),    A(2,6)  =cmplx(a2,-a5),A(2,10)=cmplx(a2,a5),                                             &
!		A(3,3) =cmplx(a1,-a4),   A(3,9)  =cmplx(a6,-a7),A(3,11)=cmplx(a2,-a5),                                            &
!		A(4,4) =cmplx(a3,a8),    A(4,8)  =cmplx(a2,a9),                                                                   &
!		A(5,1) =cmplx(a2,a5),    A(5,5)  =cmplx(a3,a10),A(5,7) =cmplx(a2,a9), A(5,8)  =cmplx(a2,a5),                      &
!		A(6,2) =cmplx(a2,-a5),   A(6,6)  =cmplx(a3,a11),A(6,7) =cmplx(a2,a5), A(6,9)  =cmplx(a2,a9),                      &
!		A(7,5) =cmplx(a2,a9),    A(7,6)  =cmplx(a2,a5), A(7,7) =cmplx(a3,a2), A(7,10) =cmplx(a2,a5),A(7,11)=cmplx(a2,a9), &
!		A(8,1) =cmplx(a6/2.,a7), A(8,4)  =cmplx(a2,a9), A(8,5) =cmplx(a2,a5), A(8,8)  =cmplx(a3,a4),A(8,10)=cmplx(a2,a9), &
!		A(9,3) =cmplx(a6/2.,-a7),A(9,6)  =cmplx(a2,a9), A(9,9) =cmplx(a3,-a4),A(9,11) =cmplx(a2,a5),A(9,12)=cmplx(a2,a9), &
!		A(10,2)=cmplx(a2,a5),    A(10,7) =cmplx(a2,a5), A(10,8)=cmplx(a2,a9), A(10,10)=cmplx(a3,-a8),                     &
!		A(11,3)=cmplx(a2,-a5),   A(11,7) =cmplx(a2,a9), A(11,9)=cmplx(a2,a5), A(11,11)=cmplx(a3,-a10),                    &
!		A(12,9)=cmplx(a2,a9),    A(12,12)=cmplx(a3,-a11)

!               1/7                2/8                3/9                4/10               5/11               6/12
		A = RESHAPE(                                                                                                           &
		 (/ cmplx(a1,a4,8)    ,cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(a2,a5,8)    ,cmplx(0.d0,0.d0,8), & 
			cmplx(0.d0,0.d0,8),cmplx(a6,a7,8)    ,cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8), & ! 1

			cmplx(0.d0,0.d0,8),cmplx(a1,a2,8)    ,cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(a2,-a5,8)   , & 
			cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(a2,a5,8)    ,cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8), & ! 2

			cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(a1,-a4,8)   ,cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8), & 
			cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(a6,-a7,8)   ,cmplx(0.d0,0.d0,8),cmplx(a2,-a5,8)   ,cmplx(0.d0,0.d0,8), & ! 3

			cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(a3,a8,8)    ,cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8), & 
			cmplx(0.d0,0.d0,8),cmplx(a2,a9,8)    ,cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8), & ! 4

			cmplx(a2,a5,8)    ,cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(a3,a10,8)   ,cmplx(0.d0,0.d0,8), & 
			cmplx(a2,a9,8)    ,cmplx(a2,a5,8)    ,cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8), & ! 5

			cmplx(0.d0,0.d0,8),cmplx(a2,-a5,8)   ,cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(a3,a11,8)   , & 
			cmplx(a2,a5,8)    ,cmplx(0.d0,0.d0,8),cmplx(a2,a9,8)    ,cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8), & ! 6

			cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(a2,a9,8)    ,cmplx(a2,a5,8)    , & 
			cmplx(a3,a2,8)    ,cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(a2,a5,8)    ,cmplx(a2,a9,8)    ,cmplx(0.d0,0.d0,8), & ! 7

			cmplx(a6/2.,a7,8) ,cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(a2,a9,8)    ,cmplx(a2,a5,8)    ,cmplx(0.d0,0.d0,8), & 
			cmplx(0.d0,0.d0,8),cmplx(a3,a4,8)    ,cmplx(0.d0,0.d0,8),cmplx(a2,a9,8)    ,cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8), & ! 8

			cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(a6/2.,-a7,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(a2,a9,8)    , & 
			cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(a3,-a4,8)   ,cmplx(0.d0,0.d0,8),cmplx(a2,a5,8)    ,cmplx(a2,a9,8)    , & ! 9

			cmplx(0.d0,0.d0,8),cmplx(a2,a5,8)    ,cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8), & 
			cmplx(a2,a5,8)    ,cmplx(a2,a9,8)    ,cmplx(0.d0,0.d0,8),cmplx(a3,-a8,8)   ,cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8), & ! 10

			cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(a2,-a5,8)   ,cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8), & 
			cmplx(a2,a9,8)    ,cmplx(0.d0,0.d0,8),cmplx(a2,a5,8)    ,cmplx(0.d0,0.d0,8),cmplx(a3,-a10,8)  ,cmplx(0.d0,0.d0,8), & ! 11

			cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8), & 
			cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(a2,a9,8)    ,cmplx(0.d0,0.d0,8),cmplx(0.d0,0.d0,8),cmplx(a3,-a11,8) /),& ! 12

		    (/12,12/),ORDER=(/2,1/) )

	END SUBROUTINE init_A

	SUBROUTINE init_D(D)
		implicit none

		REAL(dp),ALLOCATABLE,INTENT(OUT) :: D(:,:)

		ALLOCATE(D(12,12))

		D = (0.d0,0.d0)

		D(1,1) = (0.33333d0); D(2,2) = (0.33333d0); D(3,3) = (0.33333d0)
		D(4,4) = (0.33333d0); D(5,5) = (0.33333d0); D(6,6) = (0.33333d0)
		D(7,7) = (0.33333d0); D(8,8) = (0.33333d0); D(9,9) = (0.33333d0)
		D(10,10) = (0.33333d0); D(11,11) = (0.33333d0); D(12,12) = (0.33333d0)


	END SUBROUTINE init_D

END MODULE mod_init
