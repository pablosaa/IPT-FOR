! ========================================
! Variable initialization
! Program IPT
! ========================================

module variables
  implicit none
! *********************************************
!     1.) SETTINGS
! *********************************************
  ! Number of Days to analyze
  integer, parameter :: nt_days = 17

  !++++ Differentiations in cloud thickness
      !do not change unless other cloud model data used
  integer, parameter :: wd = 5
  real, parameter :: time_crit = -99.

  !++++ Maximum number of cloudy levels below 5 Km
    !do not change unless other cloud model data used
  integer, parameter :: pr = 6
          !number of frequencies
  integer, parameter :: n_freq_r = 14

  !++++ HATPRO frequencies used [GHz]
      ! exclude some if desired, change if other profiler is used
  real, dimension(n_freq_r), parameter :: freq = &
       (/ 22.240,23.040,23.840,25.440,26.240,27.840,&
          29.000,50.800,52.580,53.860,54.940,56.660,57.300,58.000 /)


  !++++ internal HATPRO numbering
      ! exclude some if desired, change if other profiler is used
  integer, parameter, dimension(n_freq_r) :: freq_num = &
       (/ 1,2,3,4,5,6,7,8,9,10,11,12,13,14 /)
  

  integer, parameter :: n_freq = size(freq_num)


  !++++ radar frequency [GHz]
  real, parameter :: rad_freq = 35.
      !change if desired and corresponding input variable

  !++++ number of T EOFs
  integer, parameter :: nT_eof = 9   !do not change

  !++++ number of q OEFs
  integer,parameter :: nq_eof = 5    !do not change
      
  !++++
  integer, parameter :: n_oef = nT_eof + nq_eof

  !++++ zenith angle
  real, parameter :: theta = 0.   !change if desired

  !++++ measurement error ground temp. (K)
  real, parameter :: e_T_gr = 0.5   !change if desired

  !++++ measurement error ground hum. (Kgm-3)
  real, parameter :: e_q_gr = 1e-4  !change if desired

  !++++ height resolution
  real,parameter :: deltaz_low = 250.  !do not change unless other cloud model data used

  real, dimension (44) :: F0O2_l92
  real, dimension (6,44) :: A_l92
  real, dimension (35) :: F0water_l92
  real, dimension (6,35) :: B_l92
  real(kind=8), dimension (44) :: F0O2_l93
  real(kind=8), dimension (6,44) :: A_l93
  real(kind=8), dimension (35) :: F0water_l93
  real(kind=8), dimension (6,35) :: B_l93


  !++++ Covariance fudging... for testing purposes if not EQ 1
  real :: fudge_sap
  real :: fudge_sap_lwc
  real :: fudge_sap_T
  real :: fudge_sap_q
  real :: fudge_se
  real :: fudge_se_tb

  !++++ path to IPT parameter data
  character(len=*), parameter :: path = '~pablo/retrieval/'
!  character(len=*), parameter :: path = '/usr/users/pablo/retrieval/' 
  character(len=*), parameter :: path_in = path//'databank/'
  character(len=*), parameter :: path_out = path//'output/'
  !++++ calculate absorption according to Liebe '92, '93 or Rosenkranz '98
  !     ('l92','l93','r98')
  character(len=*), parameter :: abs = 'r98'   !change if desired

  !++++ use dBZ to calculate LWC profile?
   !change between 'YES' or 'NO'
  character(len=*), parameter :: use_dbz = 'YES'

  !++++ initialization
  real :: d_i2_old          !convergence criterion (1e10)
  integer :: ab             !abort profile when ab = 1
  integer :: mat_inv
  integer :: cloud_prop
  integer, parameter :: n_Tq = 35
  integer, parameter :: n_Tq_2 = 2 * n_Tq
  integer, parameter :: n_Tq_new = 31
  integer, parameter :: n_Tq_new_2 = 2 * n_Tq_new
  integer, parameter :: z_level = 20

  ! ++++++++++++++++++
  ! Physical constants
  ! ++++++++++++++++++

  ! speed of light in vacuum (def) m/s
  real, parameter ::c_li = 0.299792458e9

  !Planck constant (40)J s
  real, parameter :: h = 6.6260755e-34
! Boltzmann constant
  real, parameter :: kB = 1.38066e-23 !Boltzmann constant [J/K]

  ! +++++++ Mathematical constants
  real, parameter :: pi = 3.141592      
  real, parameter :: n_e = 2.718281828459 
  
  ! --- other constants
  real, parameter :: Rw = 462., L = 2.5E6, T0 = 273.15, e0 = 611.

end module variables
