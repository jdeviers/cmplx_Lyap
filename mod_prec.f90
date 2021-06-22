MODULE mod_prec
	implicit none

	CHARACTER(LEN=*),PARAMETER :: &
	FMT_SCMPLX='(*( "(" SP F6.3 X SP F6.3 " i)" 2X ))'  , & ! Complex nbr format as float 
	FMT_LCMPLX='(*( "(" SP E10.3 X SP E10.3 " i)" 2X ))', & ! Complex nbr format as power
	FMT_SR='(*(F6.3))'                                  , & ! Real nbr format as float 
	FMT_LR='(*(E10.3))'                                     ! Real nbr format as power

	INTEGER(KIND=4),PARAMETER  ::    &
	sp = SELECTED_REAL_KIND(6,37)  , & ! 32-bits precision
	dp = SELECTED_REAL_KIND(15,307)    ! 64-bits precision

END MODULE mod_prec