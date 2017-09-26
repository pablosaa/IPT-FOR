!     ---------------------------------
!     PROGRAM IPT_KNR_BBC2
!     ---------------------------------
!     
!     Translated from the original IDL version
!     ----------------------------------------

!     BBC2: use_dbz='YES'
!     KNMI radar
!     BBC2
!     FAB 'R98'
!     changes to IPT_KNI concerning
!     1.) ground values
!     2.) missing dbz
!     3.) fudge factors
!     4.) no skydip-cals for 50.8 GHz used!
!     5.) number of maximum iterations exceeded -> fudge_sap = 0.1
!     change to KN4: also calculate T, q profiles in case of cloud-free / pure-ice scenes
!     6.) ssat check included (6.11.03) + output to RH
!     R   RPG-HATPRO radiometric profiler
!
!     *************************************************
!     Integrated Profiling Technique
!     *************************************************
!

!     --> Retrieve profiles of T, q, lwc from real data <--
!     Date of start translation 10.02.05
!     Date of last change 08.04.05

!     ************************************
!     0.) OUTLINE
!     ************************************
!     -->
!     Löhnert, U., S. Crewell, C. Simmer, 2003
!     "An integrated approach towards retrieving physically consistent profiles of temperature, humidity and cloud liquid water"
!     -->
!     GENERAL
!     This routine calculate T, q and lwc profiles from microwave TB, radar dbz, radiosonde T & q, ground T & q measurements and a priori profiles of lwc.
!     The present form was developed for the Cabauw site (NL) using data from the BBC2 measurement campaign (May 2003). TBs are taken from HATPRO (RPG), radar data from KNMI (35 GHz) cloud radar, radiosondes from De Bilt (~30 Km distance from Cabauw), and ground T, q measurements from the Cabauw tower (~300m from site). The cloud model data, which was used for the lwc a priori profile and its covariance were determined using the Issig microphysical cloud model (Univ. Bonn) with radiosondes from Essen, Germany run on a vertical resolution of 250m, which determines the final resolution of the retrieved lwc profile. A priori data for T and q profiles are taken from interpolated De Bilt radiosondes (12-hourly). The a priori covariances (Sap) are calculated from the variance between Cabauw radiosonde measurements and the interpolate DeBilt sondes.
!     The uncertainty of the TB forward modelling is contained in the Se covariance matrix, which takes into account the variance between the Liebe (1993) and Rosenkranz (1998) water vapor absorption models. Se also contain an assumed calibration error of 0.5K and the uncertainty of the FAP approximation. Note all errors used to calculate Se are random errors: systematic errors cannot be described within Se. The user may choose between different absorption models and the differences in results will show the sensitivity of IPT towards systematic absorption model differences.
!     Systematic calibrations offsets must be considered before running IPT.

!     NOTIFI!ATION
!     The user principally has the possibility of adapting this routine to different sites or to different cloud model data. But beaware, the SETTING must be carefully changed and input files (infile) and program parameters (GET_MATRICES, GET_AP, GET_COEFF) must be adjusted accordingly.

!     1.) SETTING: algorithm setting -> please read carefully
!     2.) BEGIN LOOP OVER DAYS
!     3.) BEGIN LOOP OVER PROFILES
!     4.) BEGIN ITERATION: iteration loop to find best solution
!     5.) CALCULATE K_J BY COMBINING K1 AND K3: calculation of Jacobian matrix via pertubation
!     6.) OPTIMAL ESTIMATION: apply optimal estimation equations
!     7.) CHECK FOR CONVERGENCE CRITERION: do I need another iteratin or not?
!     8.) SOLUTION CHECKS: once covergence has been reached: is the solution reasonable?
!          (several error checks and calculations)
!     9.) WRITE RETRIEVAL RESULTS TO FILE

PROGRAM IPT_KNR_BBC2
  use variables
  implicit none

  !     *********************************************
  !     1.) SETTINGS
  !     *********************************************

  integer :: istat, status
  real :: t_start, t_end
  real, dimension(10, n_freq_r) :: coeff_fap
  real, dimension(n_freq_r, n_freq_r) :: Se_fap_org
  real, dimension(2, n_freq_r) :: bias_coeff_fap
  real(kind=8), dimension (7,44) :: array1
  real(kind=8), dimension (7,35) :: array2
  real(kind=8), dimension(wd, n_freq_r+pr, n_freq_r+pr) :: Se_a
  real(kind=8), dimension(wd, n_Tq_2+pr, n_Tq_2+pr) :: Sap_a
  real(kind=8), dimension(wd, n_Tq) :: T_ap_a,p_ap_a,q_ap_a
  real(kind=8), dimension(wd, pr) :: lwc_ap_a, c_a, d_a
  real(kind=8), dimension(n_Tq_2, n_Tq_2) :: cov_rs
  character(len=6), dimension(nt_days) :: date   !date 6 characters and NT_DAYS number of days
  character(len=66), dimension(2) :: aa

  call cpu_time(t_start)

  call open_file(path_in//'ipt_bbc2.list',1)
  read(unit=1,fmt='(A)') date
  close(unit=1)

  !     +++++++++++++++++++++++++++++++++++++
  !     ++++ READ ABSORTION PARAMETERS   ++++
  !     +++++++++++++++++++++++++++++++++++++

  call open_binary(path_in//'O2.dat',4,44*6*44)
  read(4,rec=1) F0O2_l92, A_l92
  close(4)

  call open_binary(path_in//'H2O.dat',4,35*6*35)
  read(4,rec=1) F0water_l92,B_l92
  close(4)
  !_________________________________________
  call open_file(path_in//'oxygen_l93.dat',16)
  read(unit=16,fmt='(A)') aa
  read(unit=16,fmt=10) array1
10 format(F12.6,3F12.3,F5.2,2F8.3)
  close(unit=16)
  F0O2_l93 = array1(1,:)
  A_l93 = array1(2:7,:)
  !__________________________________________
  call open_file(path_in//'water_l93.dat',17)
  read(unit=17,fmt='(A)') aa
  read(unit=17,fmt=20) array2
20 format(F11.6,F11.5,F7.3,F8.3,1X,3F6.2)
  close(unit=17)
  F0water_l93 = array2(1,:)
  B_l93 = array2(2:7,:)

  !     ++++
  !     Read covariance matrices
  call get_matrices(Se_a, Sap_a, path_in//'cov_35_r.dat')
  
  !     ----
  !     ++++
  !     Read fap coefficients
  call get_coeff_fap(path_in//'coeff_fap_'//abs//'_r.dat',coeff_fap)

  !     ----
  !     ++++
  !     Read fap forward model error covariance and bias_coeff
  call get_bias_coeff(path_in//'error_fap_'//abs//'_r.dat',&
       Se_fap_org,bias_coeff_fap)

  !     ----

  !     ++++
  !     Read a priori data in parameter space
  !     (Here: statistical mean Issig model output, only lwc data used)
  lwc_ap_a = -99.
  call get_ap(path_in//'ap.dat',T_ap_a,p_ap_a,q_ap_a,lwc_ap_a)

  !     ----

  !     ++++
  !     Read coefficients for Z-LWC relation
  c_a = -99.
  d_a = -99.
  call get_coeff(path_in//'coeff.dat',c_a,d_a)

  !     ----

  !     ++++
  !     Read RS covariance (variability from RS to RS)
  cov_rs = 0
  call open_file(path_in//'cov_rsrs_bbc2.dat',12)
  read(12,*) cov_rs
  close(12)
  !     ----

  !     ++++++++++++++++++++++++++++
  !     BEGIN LOOP OVER DAYS
  !     ++++++++++++++++++++++++++++

  call over_days(date,cov_rs,Sap_a,Se_a,Se_fap_org,c_a,d_a,lwc_ap_a,&
       coeff_fap, bias_coeff_fap)


  call cpu_time(t_end)


  if ((t_end-t_start).LT.60.) then
     print *, 'Total elapsed time: ', t_end-t_start, 'seconds.' 
  else if ((t_end-t_start).LT.3600.) then
     print *, 'Total elapsed time: ', (t_end-t_start)/60., 'minutes.'
  else
     print *, 'Total elapsed time: ', (t_end-t_start)/3600., 'hours.'
  end if
  !     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  STOP
END PROGRAM IPT_KNR_BBC2







!     =======================================
!     SUBROUTINES
!     ======================================

!     ++++
!     Subroutine to check the state of a file
subroutine allocate_lun(lun,fname)
  implicit none
  integer, intent(in) :: lun      !, MINLUN, MAXLUN
  character(len=*), intent(in) :: fname
  !      parameter   (MINLUN = 50, MAXLUN = 100)
  logical :: log
  !      do lun = MINLUN, MAXLUN 
  inquire(unit=lun, opened=log)
  if (.not. log) return 
  !     end do
  print *, 'Can''t allocate file: ',fname
  stop
end subroutine allocate_lun
!     ----

!     ++++
!     Subroutine to open an unformatted file. parameters: filename, unit=
subroutine open_binary(fname,lun,reco)
  implicit none
  character(len=*), intent(in) :: fname
  integer, intent(in) :: lun,reco
  integer :: istat
  call allocate_lun(lun,fname)
  open(unit=lun, file=fname, form='UNFORMATTED', status='old',&
        access='direct',recl=4*reco,iostat=istat)
  if (istat /= 0) then
     write (*,'(A,A,A,I3)') 'Input binary file: ', fname,&
           ', open failed = ', istat 
  endif
  return
end subroutine open_binary
!     ----

!     ++++
!     Subroutine to open a file. Parameters: filename, unit=
subroutine open_file(fname,lun)
  implicit none
  character(len=*), intent(in) ::  fname
  integer, intent(in) :: lun
  integer :: istat
  call allocate_lun(lun, fname)
  open(unit=lun,file=fname,form='formatted',status='old', action='read',&
       iostat=istat)
  if (istat /= 0) then
     write (*,'(A,A,A,I3)') 'Input file:', fname, ', open failed =', istat
  endif
  return
end subroutine open_file
!     ----

!     ++++
!     Subroutine GET MATRICES
subroutine get_matrices(Se,Sap,dat)
  use variables , only : pr, wd, n_freq_r, n_Tq_2
  implicit none
  character(len=*),intent(in) :: dat
             ! Se: Fehlermatrix TB, dBZ
  real(kind=8),dimension(wd,n_freq_r+pr,n_freq_r+pr), intent(out):: Se
           ! Sap: T,q,LWC Kovarianzen
  real(kind=8),dimension(wd,n_Tq_2+pr,n_Tq_2+pr), intent(out) :: Sap
  !--- local variables ---
  character(len=5) :: aa
  integer :: hnn, k

  hnn = n_Tq_2
  call open_file(dat,11)
  do k=4,0,-1
     read(11,'(A)') aa
     read(11,*) Sap(5-k,1:hnn+pr-k,1:hnn+pr-k)   !-Sap_0
     read(11,'(A)') aa
     read(11,*) Se(5-k,1:n_freq_r+pr-k,1:n_freq_r+pr-k)  !-Se_0
  end do
  close(11)
  return
end subroutine get_matrices
!     ----

!     ++++
!     GET_COEFF_FAP
subroutine get_coeff_fap(dat,coeff_fap)
  use variables, only : n_freq_r
  implicit none
  character(len=*),intent(in) :: dat
  real, dimension(10,n_freq_r), intent(out) :: coeff_fap
  !-- internal variables --
  character(len=100) :: aa
  integer :: k, j

  call open_file(dat,12)
  do k=1,n_freq_r
     read(12,'(A)') (aa, j=0,2)
     read(12,*) coeff_fap(:,k)  !-xx
  end do
  close(12)
  return
end subroutine get_coeff_fap
!     ----

!     ++++
!     GET_BIAS_COEFF
subroutine get_bias_coeff(dat,Se_fap_org,bias_coeff_fap)
  use variables , only : n_freq_r
  implicit none
  character(len=*), intent(in) :: dat
  real, dimension(n_freq_r,n_freq_r), intent(out) :: Se_fap_org
  real, dimension(2,n_freq_r), intent(out) :: bias_coeff_fap
  !-- internal variables --
  character(len=80) :: aa
  integer :: k

  call open_file(dat,12)
  read(12,'(A)') (aa ,k=0,24)
  read(12,*) Se_fap_org
  read(12,'(A)') aa
  read(12,*) bias_coeff_fap
  close(12)
  return
end subroutine get_bias_coeff
!     ----

!     ++++
!     Read a prior data (normal space) from files
subroutine get_ap(dat,T_ap_a,p_ap_a,q_ap_a,lwc_ap_a)
  use variables , only : pr, wd, n_Tq
  implicit none
  character(len=*), intent(in) :: dat
  real(kind=8),dimension(wd,n_Tq),intent(out) :: T_ap_a, p_ap_a , q_ap_a
  real(kind=8),dimension(wd,pr),intent(inout) :: lwc_ap_a

  integer :: k
  character(len=8) :: aa

  call open_file(dat,12)

  do k=1,5
     read(12,'(A)') aa
     read(12,*) T_ap_a(k,:)
     read(12,'(A)') aa
     read(12,*) q_ap_a(k,:)
     read(12,'(A)') aa
     read(12,*) p_ap_a(k,:)
     read(12,'(A)') aa
     read(12,*) lwc_ap_a(k,1:pr+k-5)
  end do
  close(12)
  return
end subroutine get_ap
!     ----

!     ++++
!     Read Z-LWC coefficients from files
subroutine get_coeff(dat,c_a,d_a)
  use variables , only : pr, wd
  implicit none
  character(len=*), intent(in) :: dat
  real(kind=8),dimension(wd,pr),intent(inout) :: c_a,d_a

  integer :: k
  character(len=10) :: aa

  call open_file(dat,12)

  do k=1,5
     read(12,'(A)') aa
     read(12,*) c_a(k,1:pr+k-5)
     read(12,'(A)') aa
     read(12,*) d_a(k,1:pr+k-5)
  end do
  close(12)
  return
end subroutine get_coeff
!     ----
