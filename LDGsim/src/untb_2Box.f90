Program UNTB_2Box
implicit none

integer, parameter :: NBox = 2    !Number of regions
integer   :: Nens     = 2
integer   :: Ninv     = 10
integer   :: Duration = 100  !Number of simulation days
real      :: dtdays   = 0.1  !Time step
real      :: theta0   = 2d0  !N*nu(mutation rate per capita)
real      :: beta0    = 1d0  !normalized birth rate (per day)
real      :: delta0   = 0.9d0  !normalized abundance-independent death rate (d-1)
real      :: delta01  = 0.001d0  !carrying capacity K = (beta - delta)/delta1
real      :: Ea       = 0.65d0
real      :: pmig     = 0.01d0
real      :: theta(NBox) = [0., 0.]
real      ::  beta(NBox) = [0., 0.] 
real      :: delta(NBox) = [0., 0.]
real      ::delta1(NBox) = [0., 0.]
real      ::  Temp(NBox) = [0., 0.]

integer :: i
integer, parameter  :: namlst = 20
real, external      :: TEMPBOL

!Recording time
real                 :: start, finish

namelist /parameters/ Nens, Ninv, Duration, dtdays, theta0, beta0, delta0, delta01, Ea, pmig, Temp
!End of declaration-----------------------------------

CALL cpu_time(start)

!  open the namelist file and read the parameters.
open(namlst,file='params.nml',status='old',action='read')
read(namlst,nml = parameters)
close(namlst)
  
!Calculate theta, beta, and delta based on temperature in two boxes
do  i = 1, NBox
    theta(i) = theta0 * TEMPBOL(Ea, Temp(i))
    beta(i)  =  beta0 * TEMPBOL(Ea, Temp(i))
    delta(i) = delta0 * TEMPBOL(Ea, Temp(i))
   delta1(i) =delta01 * TEMPBOL(Ea, Temp(i))
enddo

! Repeat the simulation for Nsim times

CALL neutral_2box(Nens, Ninv, NBox, Duration, dtdays, theta, beta, delta, delta1, pmig)

CALL CPU_TIME(finish)
PRINT '("Time = ",1pe12.2," hours.")', (finish-start)/3600.0 
END PROGRAM

PURE REAL FUNCTION TEMPBOL(Ea,tC)
implicit none
!DESCRIPTION:
!The temperature dependence of plankton rates are fomulated according to the Arrhenuis equation. 
! tC: in situ temperature
! Tr: reference temperature
!
!INPUT PARAMETERS:
real, intent (in) :: Ea, tC
! boltzman constant constant [ eV /K ]
real, parameter   :: kb = 8.62d-5, Tr = 15D0

TEMPBOL = exp(-(Ea/kb)*(1D0/(273.15 + tC)-1D0/(273.15 + Tr)))
return 
END FUNCTION TEMPBOL
