! ++++ DE-ATTENUATE dbz
! * q from radiosonde
! * lwc from radiometer

subroutine de_attenuate(c,d,dbz,dbz_da,h_true,&
     lwp_mic,T_rs,q_rs,p_rs,z_final,base_m,top_m,&
     dbz_hr,dbz_hr_da,z_hr,h_true_m,nn,z_lr_lwc)

  use variables, only: pi, theta, rad_freq
  use absorption, only: abwvr98, abliq
  IMPLICIT NONE
      real(kind=8), dimension(:),intent(in) :: c, d
      real, dimension(:), intent(out) :: dbz_da, dbz_hr_da
      integer, intent(in) :: h_true, nn
      real, intent(in) :: lwp_mic
      real, dimension(:), intent(in) :: dbz, T_rs, q_rs, p_rs, z_final
      real, dimension(:), intent(in) :: base_m, top_m, dbz_hr
      real, dimension(:), intent(in) :: z_hr, z_lr_lwc
      integer, dimension(:), intent(in) :: h_true_m
!     --
      ! Local variables
      integer :: i, j, k, nt
      integer :: n_q
      real :: mu, number_e, tau_gas, deltaz0, deltaz
      real :: lwp, corr, tau_under
      real(kind=8) :: lwp_at, absg, absc
      real(kind=8), dimension(h_true) :: lwc_at, lwc_corr
      real, allocatable, dimension(:) :: tau_gas_count

      dbz_da = 1e-10      ! de-attenuated dbz
      dbz_hr_da = dbz_hr
      lwp = lwp_mic
      mu = cos(theta*pi/180.)
      number_e = exp(1.)
      tau_gas = 0. ; absg = 0.
      n_q = size(q_rs);
      deltaz0 = z_final(2) - z_final(1)
      if (lwp .LT. 0.) lwp = 0.

      ! Water vapor: Rosenkranz '98
      j = 1 ; allocate(tau_gas_count(n_q)) ; tau_gas_count = 0.
      do i = 1, n_q
         if(i .LT. n_q) deltaz = z_final(i+1) - z_final(i)
         absg = abwvr98(real(1000.*q_rs(i),KIND=8), real(T_rs(i),KIND=8), &
              real(p_rs(i)/100.,KIND=8), rad_freq)
         absg = absg/1e3
         if(i .GT. 1) tau_gas = tau_gas + (absg/mu)*deltaz
         if(z_final(i) .GE. base_m(1)) then
            tau_gas_count(j) = tau_gas
            j = j + 1
         end if
      end do
      lwc_at = 0.
      do i = 1, h_true
         lwc_at(i) = (dbz(i) - c(i))/(10.*d(i))  ! log
         lwc_at(i) = (10.**(lwc_at(i)))        ! lwc_at linear
         lwc_at(i) = lwc_at(i)/1e3            ! [Kg/m3]
      end do
      lwp_at = sum(1e3*lwc_at)*deltaz0
      corr = lwp/lwp_at
      lwc_corr = corr*lwc_at     ! meets lwp
      tau_under = tau_gas_count(1)    ! initialize
      do i = 1, h_true
         dbz_da(i) = dbz(i) + 20.*tau_under*(log10(number_e))
         where((z_hr .GT. z_lr_lwc(i)-deltaz0/2.)&
              .AND.(z_hr .LT. z_lr_lwc(i)+deltaz0/2.)&
              .AND.(dbz_hr .GT. -100.)) dbz_hr_da = dbz_hr + 20.*tau_under*log10(number_e)
         ! ++ Liebe et al., 1993 ++

         absc = abliq(real(lwc_corr(i)*1e3,KIND=8), real(T_rs(i),KIND=8), rad_freq)
         absc = absc/1000.
         tau_under = absc*deltaz0 + tau_gas_count(i)
      end do
      deallocate(tau_gas_count)
return
end subroutine de_attenuate

