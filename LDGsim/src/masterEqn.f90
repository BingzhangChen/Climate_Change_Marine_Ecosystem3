!A test simulation only simulating drift and dispersal between two neighboring
!communities under two different temperatures based on the master equation approach
SUBROUTINE neutral_2box(Nens, Nini, NBox, Duration, dtdays, theta, beta, delta, delta1, pmig)
USE mGf90
IMPLICIT NONE

integer, intent(in)  :: Nens  !The number of ensemble runs

integer, intent(in)  :: Nini  !Initial total number of individuals 

!Number of boxes
integer, intent(in)  :: NBox

integer, intent(in)  :: Duration    !the number of days for simulation

real,    intent(in)  :: theta(NBox)       !biodiversity variable (product of community size and per capita speciation rate)

!Time step (in days)
real,    intent(in)  :: dtdays     !how much fraction of a day

!Density independent birth rate 
real,    intent(in)  :: beta(NBox) 

!Density independent death rate
real,    intent(in)  :: delta(NBox)

!Density dependent death rate
real,    intent(in)  :: delta1(NBox)

!mixing (dispersal) rate
real,    intent(in)  :: pmig

integer              :: boxInd(NBox)     !The indexes of all boxes (1:NBox)

!Define individuals
TYPE    IND
	real         :: b_time     = 1.  !Time needed for next birth
	real         :: d_time     = 1.
	real         :: clock      = 0.  !Biological clock (in days) indicating the current life cycle
	real         :: Trait      = 0.  !A conserved trait (0~1) independent of fitness
	integer      :: ibox       = 1   !The box that the individual is residing
	integer      :: spname     = 0   !Species name
END TYPE IND

!The community of all the individuals
Type (IND), allocatable  ::  Z(:)

!A scratch community of Z
Type (IND), allocatable  ::  Zc(:)
Type (IND), allocatable  ::  Zd(:)

real                 :: nu(NBox)     !Per capita Speciation rate
real                 :: Kcc(NBox)    !Estimated carrying capacity of each box

integer              :: Abun_Box(NBox)     !Abundances in each box

!Scratch variable
integer              ::  i,j,k,m,jj,im, LL, LiD

!Total number of iterations
integer              ::  NT         = 0

!Total number of individuals
integer              ::  N          = 0

integer              ::  N_birth(NBox)       !Number of births during the simulation
integer              ::  N_death(NBox)       !Number of deaths during the simulation
integer              ::  N_spec(NBox)        !Number of speciations during the simulation
integer              ::  N_ext(NBox)         !Number of local extinction events during the simulation
integer              ::  N_emig(NBox)        !Number of emigrations during the simulation
integer              ::  S_emig(NBox)        !Number of emigrated new species during the simulation

integer              ::  status     = 0        

!Time interval to save
integer              ::  nsave      = 0  !Save results every 1 day

!Current day
real                 ::  curr_day   = 0.

!Scratch random number
real                 :: rnd = 0.

!Variance of trait mutation
real                 :: TraitVAR(1,1) = 0.01
real                 :: rnd2(1)  = [-999.]	
!File name of output
character(len=6)     :: indfile  = 'Ind'
character(len=7)     :: ratfile  = 'Rat'

!unit of filename
integer, parameter   :: indunit  = 10
integer, parameter   :: ratunit  = 11

!the indexes of individuals that die
integer, allocatable :: ideath(:)

integer, allocatable :: add_death(:)

integer              :: NewName = 0   !A scratch new species name

integer, allocatable :: AllSpnames(:)  !Recording all the species names
integer, allocatable :: add_Spnames(:) !A scratch variable to dynamically add a species name

logical              :: match         = .true.
logical              :: conspecific   = .true.
logical              :: within_ideath = .true.
logical              :: new_immigrant = .true.

character(len=10), parameter :: format_string = "(A3,I0)"

!External functions
real,    external    :: rexp
integer, external    :: sample, neighbor
!End of declaration-----------------------------------

NT    = int(dble(Duration)/dtdays)  !Number of iterations

nsave = int(1./dtdays)  !Save results every 1 day

CALL RANDOM_SEED()

!Initialize boxindex
do i = 1, NBox
	 boxInd(i) = i     !The indexes of all boxes (1:NBox)

	 !Calculate carrying capacity
     Kcc(i)    = (beta(i) - delta(i))/delta1(i)
enddo

do i = 1, NBox
     nu(i)     = theta(i)/(sum(Kcc(:))/dble(NBox))  !Calculate per capita speciation rate
enddo

!Start individual simulation
DO im = 1, Nens

	N = Nini

	!Initialize all the individuals (including living and dead)
	allocate(Z(N),  STAT = Status)
	IF (Status /= 0) STOP "*** Problem in allocating Z ***"

	allocate(AllSpnames(N), stat = Status)
	IF (Status /= 0) STOP "*** Problem in allocating AllSpnames ***"

	!Initialize trait
	call random_number(rnd)

	do k = 1, N  !N should be fixed for each iteration
		!Initialize the individuals
		Z(k)%clock   = 0.d0
		Z(k)%ibox    = sample(NBox, boxInd)
		Z(k)%spname  = 1
		Z(k)%Trait   = rnd

		!Fill AllSpnames (start from only one species)
		AllSpnames(k)= 1  

		Z(k)%b_time = rexp(beta(Z(k)%ibox))
		Z(k)%d_time = rexp(delta(Z(k)%ibox))
	enddo

	!Count how many individuals in each box (to initialize d_time)
	Abun_Box(:) = 0
	do k = 1, NBox
		do i = 1, N
			if (Z(i)%ibox .eq. k) Abun_Box(k) = Abun_Box(k) + 1
		enddo
	enddo

	do k = 1, N
		Z(k)%d_time = rexp( delta(Z(k)%ibox) + delta1(Z(k)%ibox) * Abun_Box(Z(k)%ibox) )
	enddo

	!Create file and Save initial condition

	!Create new file name
  	WRITE(indfile, format_string) 'Ind', im

	OPEN(indunit, file=indfile, status='replace')
	WRITE(indunit, 1000) 'Timestep', 'Day', 'IBOX', 'Spname', 'Trait'

	!write into file
	do i = 1, N
		write(indunit, 1001) 0, 0., Z(i)%ibox, Z(i)%spname, Z(i)%Trait
	enddo

	CLOSE(indunit)

	!Initialize the numbers of births, deaths, emigrations, speciations, and extinction rates
	do i = 1, NBox
		N_birth(i) = 0
		N_death(i) = 0
		N_spec(i)  = 0
		N_ext(i)   = 0
		N_emig(i)  = 0
		S_emig(i)  = 0
	enddo

	!Write data of N_birth etc into the rate file
    write(ratfile, format_string) 'Rat', im

	OPEN(ratunit, file=ratfile, status='replace')
	WRITE(ratunit, 1002) 'Timestep', 'Day', 'Box', 'N_birth', 'N_death', 'N_spec', 'N_ext', 'N_emig', 'S_emig'

	CLOSE(ratunit)

	!Start iteration
	DO k = 1, NT

		if (k .gt. 1 .and. N .le. 1) then
		   write(6, *) 'Community size too small! STOP simulation!'
		   exit
		else if (N .gt. 1e8) then
		   write(6, *) 'Community size too large! STOP simulation!'
		   exit
		endif

		!Current day
		curr_day  = dble(k)*dtdays

		DO j = 1, N

		   !Advance internal biological clock
		   Z(j)%clock = Z(j)%clock + dtdays

		   !Migrate if the probability is smaller than pmig
		   call random_number(rnd)

		   If (rnd .lt. pmig) then  !migrate to the nearest box

			!Update N_emig
	       N_emig(Z(j)%iBox) = N_emig(Z(j)%iBox) + 1

			!Find the neighbor box
           LL = neighbor(NBox, Z(j)%iBox)

			!Update S_emig (if this individual is a new species for the box emigrated to)
			new_immigrant = .true.
			do jj = 1, size(Z)
				if (Z(jj)%iBox .eq. LL .and. Z(jj)%Spname .eq. Z(j)%Spname) then
					new_immigrant = .false.
					exit
				endif
			enddo

			if (new_immigrant) S_emig(Z(j)%iBox) = S_emig(Z(j)%iBox) + 1

			!change the ibox for Z(j)
	        Z(j)%iBox = LL

			  !Update beta and delta
			  Z(j)%b_time  = rexp( beta(LL))
			  Z(j)%d_time  = rexp(delta(LL) + delta1(LL)*Abun_Box(LL))

		   Endif

	     IF (Z(j)%b_time .le. Z(j)%d_time) THEN
	        !if the clock exceeds birth time but not the death time, then it divides;
	        if (Z(j)%clock .ge. Z(j)%b_time) then
			
	           N_birth(Z(j)%iBox) = N_birth(Z(j)%iBox) + 1

	           !Reset the clock
	           Z(j)%clock = 0.d0
			
				ALLOCATE(Zc(1+size(Z)), stat=status)
				IF (Status /= 0) STOP "*** Problem in allocating Zc***"
				Zc(1:size(Z))   = Z

	           !Divide (copy all the attributes of this individual to a new individual)
	           Zc(size(Z) + 1) = Z(j)

	           call move_alloc(Zc, Z)   !Zc is deallocated and length of Z increases by 1

	           !Reset the birth and death time for both the previous and the new individual 
	           Z(j)%b_time = rexp( beta(Z(j)%ibox))
			   Z(j)%d_time = rexp(delta(Z(j)%ibox) + delta1(Z(j)%ibox)*Abun_Box(Z(j)%ibox) )
			
	           Z(size(Z))%b_time = rexp( beta(Z(j)%iBox))
	           Z(size(Z))%d_time = rexp(delta(Z(j)%iBox) + delta1(Z(j)%ibox)*Abun_Box(Z(j)%ibox))

	           !A small probability to generate a new species during birth
	           call random_number(rnd)
			
	           if (rnd < nu(Z(j)%ibox)) then
	             N_spec(Z(j)%iBox) = N_spec(Z(j)%iBox)  + 1
	             NewName           = maxval(AllSPnames) + 1
	             Z(size(Z))%spname = NewName

							 !get a new trait for the mutated individual
	             rnd2 = srand_mtGaus(1, [Z(size(Z))%Trait], TraitVAR)   
							 Z(size(Z))%Trait  = rnd2(1)

	             !Update species name
				 ALLOCATE(add_Spnames(size(AllSPnames) + 1), stat=status)
				 IF (Status /= 0) STOP "*** Problem in allocating add_Spnames***"
				 add_Spnames(1:size(AllSPnames))   = AllSPnames
				 add_Spnames(size(AllSPnames) + 1) = NewName
				 call move_alloc(add_Spnames, AllSPnames)  !AllSPnames has been updated and add_Spnames has been deallocated!
				
	           endif
	          endif
	        ELSE
	          !Otherwise, it dies
	          if (Z(j)%clock >= Z(j)%d_time) then
			
	            N_death(Z(j)%iBox)  = N_death(Z(j)%iBox) + 1

	            !Do not remove this individual now, but record this index (needs to do this because two deaths can occur at the same time step)
	            !Update ideath

				if (allocated(ideath)) then
			        allocate(add_death(size(ideath) + 1), stat=status)
				    IF (Status /= 0) STOP "*** Problem in allocating add_death***"
			        add_death(1:size(ideath))   = ideath
	                add_death(size(ideath) + 1) = j
	                call move_alloc(add_death, ideath)   !Update ideath and add_death has been deallocated
				else
	               allocate(ideath(1))
				   ideath(1)  = j 
				endif
		      endif
	        ENDIF
	  ENDDO

		!Remove the dead individuals
		IF (allocated(ideath) .and. size(ideath) > 0) THEN

	    !Construct a new Zd with length of size(Z)-length(ideath) 
			allocate(Zd(size(Z) - size(ideath)), stat=Status)
			IF (Status /= 0) STOP "*** Problem in allocating Zd ***"

			jj = 0   !Counting the index in Zd
			do j = 1, size(Z)
				match = .false.
				do m = 1, size(ideath)
					if (j .eq. ideath(m)) then
						match = .TRUE.

			      !Update the local extinction rate 
						!Check whether this individual is the last individual of the local community
						conspecific = .false.
						do LL = 1, size(Z)

							 !Find if there is any individual that is in the same box and has the same species name
							 !but is not itself
						   if (j .ne. LL .and. Z(LL)%ibox   .eq. Z(j)%ibox         &
						                 .and. Z(LL)%Spname .eq. Z(j)%Spname) then
								  
									!Double check whether this individual of the same species is not within ideath 
								  within_ideath = .false.
								  do LiD = 1, size(ideath)
										if (LL .eq. LiD)  then
								      within_ideath = .true.
											exit
										endif
									enddo
								  
									if (within_ideath) then
										cycle
									else
										conspecific = .true.
										exit
									endif
						   endif
						enddo

						if (.not. conspecific) then
							!Local extinction occurs
							N_ext(Z(j)%ibox) = N_ext(Z(j)%ibox) + 1
						endif

						exit
					endif
				enddo  !End of finding whether the individual has died or not

				if (.not. match) then   !Include this individual in the new community Zc
					jj     = jj + 1
					Zd(jj) = Z(j)
				endif
			enddo	!End of constructing Zd

			!Transfer Zd to Z
			call move_alloc(Zd, Z)  !Update Z and Zd is deallocated
			deallocate(ideath)
	  ENDIF  !End of removing dead individuals

		!Update N
	  N = size(Z)

		!Update Abun_box
		Abun_Box(:) = 0
		do i = 1, NBox
			do m = 1, N
				if (Z(m)%ibox .eq. i) Abun_Box(i) = Abun_Box(i) + 1
			enddo
		enddo

		!Save results at time2save
		If (mod(k, nsave) .eq. 0) Then
		  OPEN(indunit, file=indfile, status='old', action='write', position='append')
		  	do i = 1, size(Z)
		  		write(indunit, 1001) k, curr_day, Z(i)%ibox, Z(i)%spname, Z(i)%Trait
		  	enddo
		  CLOSE(indunit)

		  !Save rate data and reset N_birth etc.
	    OPEN(ratunit, file=ratfile, status='old', action='write', position='append')
	    DO i = 1, NBox
	    	WRITE(ratunit, 1003) k, curr_day, i, N_birth(i), N_death(i), N_spec(i), N_ext(i), N_emig(i), S_emig(i)
	    ENDDO
	    CLOSE(ratunit)

	    do i = 1, NBox
	    	N_birth(i) = 0
	    	N_death(i) = 0
	    	N_spec(i)  = 0
				N_ext(i)   = 0
	      N_emig(i)  = 0
	      S_emig(i)  = 0
	    enddo

		Endif

	ENDDO  !End of each individual simulation

	if (allocated(Z))          deallocate(Z)
	if (allocated(AllSpnames)) deallocate(AllSpnames)
ENDDO

1000 format(5(A8, 2x))
1001 format(I8, 2x, F10.1, 2x, 2(I0, 2x), F12.6)
1002 format(9(A8, 2x))
1003 format(I8, 2x, F10.1, 2x, I1, 2x, 6(I0, 2x))

END subroutine neutral_2box

real function rexp(u) result(y)
implicit none

real, intent(in) :: u   !Rate parameter
real             :: rnd !a random number between 0 and 1

y = 0.d0

if (u < 0.d0) then
        write(6, *) 'The rate parameter has to be positive!'
        return
endif

call random_number(rnd)

y = -log(rnd)/u
return
END function rexp

Integer function sample(nn, x) result(y)
!A function to randomly sample one element from a vector of integers
implicit none
integer, intent(in) :: nn
integer, intent(in) :: x(nn)

real    :: r = 0.d0
integer :: i

call random_number(r)

if (r .le. 1./dble(nn)) then
	y = x(1)
else if (r .gt. (1. - 1./dble(nn))) then
	y = x(nn)
else
	do i = 2, (nn-1)
		if (r .gt. dble(i-1)/dble(nn) .and. r .le. dble(i)/dble(nn)) then
			y = x(i)
			return
		endif
	enddo
endif

return
END function sample

integer function neighbor(N_box, curr_box)
implicit none
integer, intent(in)  :: N_box      !Total number of boxes
integer, intent(in)  :: curr_box   !Current box index

!Find a neiboring box (for dispersal)
if (curr_box .eq. N_box) then
	neighbor = 1
else
	neighbor = curr_box + 1
endif
return
end function neighbor