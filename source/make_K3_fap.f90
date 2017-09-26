! +++++++++++++++++++++++++++++++++++++++++++++++++
! SUBROUTINE MAKE_K3_FAP
! -------------------------------------------------
! Part of the program IPT_KNR_BBC2.F90
! This subroutine is called by OVER_DAYS.F90 program.
! ---
! Important: Use the character parameter 'R98' or 'L93'.
subroutine make_K3_fap(c_stp, z_lr, x_n, z_final, p_rs, h_true, inc, coeff, bias, cloud, diff)
  use variables, only: freq, n_freq, n_Tq, n_Tq_2, deltaz_low, theta, pi
  use absorption
  implicit none
  character(len=*), intent(in) :: c_stp
  integer, intent(in) :: h_true, cloud
  real(kind=8), intent(in) :: x_n(:)
  real, intent(in) :: z_lr(:),z_final(:),p_rs(:),inc(:),coeff(:,:),bias(:,:)
  real(kind=8), intent(out) :: diff(:,:)

  ! Local variables
  integer :: i, i_lwc, j, j_inc
  real(kind=8) :: inc_q_new
  real(kind=8) :: mu, tau_final(n_Tq-1,n_freq)
  real(kind=8), dimension(h_true) :: lwc
  real(kind=8), dimension(n_Tq) :: lwc_final, T_final, T_prepp, T_prepm, &
       q_final, q_prepp, q_prepm, lwc_prepp, lwc_prepm
  real(kind=8), dimension(n_freq) ::  TB_final, TB_bf, TB_p, TB_m
  real, dimension(size(p_rs)) :: p_final

  mu = COS(theta*pi/180.)
  T_final = x_n(1:n_Tq)
  q_final = x_n(n_Tq+1:n_Tq_2)
  p_final = p_rs            !real(p_rs, KIND=8)
  lwc_final = 0.
  if (cloud .GT. 0) then
     lwc = x_n(n_Tq_2+1:n_Tq_2+h_true)
     ! position lwc of final grid
     i_lwc = 1
     do i=1,n_Tq
        if(z_lr(i_lwc) .EQ. (z_final(i)+deltaz_low/2.)) then
           lwc_final(i) = lwc(i_lwc)
           i_lwc = i_lwc + 1
        end if
        if (i_lwc .GT. h_true) EXIT
     end do

     where(lwc_final .NE. 0.) lwc_final = (10.**(lwc_final/10.))*1e-3
  end if
  ! Determine components of K3 : diff matrix (2*n_Tq+h_true, n-freq)
  diff = 0.

  ! * T sensitivity
  T_sens:  do i = 1, n_Tq
     TB_final = 0.
     tau_final = 0.
     T_prepp = T_final
     T_prepm = T_final
     T_prepp(i) = T_final(i) + inc(i)
     T_prepm(i) = T_final(i) - inc(i)
     call tau_calc_fap(z_final,T_prepp,p_final,q_final,lwc_final,mu,&
          coeff,c_stp,tau_final)  ! c_stp: 'R98' or 'L93'

     !TB_final = tb_calc_pl(T_prepp,tau_final,mu)
     call tb_calc_pl(T_prepp,tau_final,mu,TB_final)
     TB_bf = (TB_final - bias(1,:))/bias(2,:)
     TB_p = TB_bf

     TB_final = 0.
     tau_final = 0.
     call tau_calc_fap(z_final,T_prepm,p_final,q_final,lwc_final,mu,&
          coeff,c_stp,tau_final)  ! c_stp: 'R98' or 'L93'
     call tb_calc_pl(T_prepm,tau_final,mu,TB_final)
     TB_bf = (TB_final - bias(1,:))/bias(2,:)
     TB_m = TB_bf
     diff(i,:) = (TB_p - TB_m)/(2.*inc(i))
  end do T_sens

  ! * q sensitivity
  Q_sens:  do i = 1, n_Tq
     TB_final = 0.
     tau_final = 0.
     q_prepp = q_final
     q_prepm = q_final
     inc_q_new = inc(n_Tq+i)
     ENDLESS: do
        if(((q_final(i)-inc_q_new).LT.0.) .AND. (q_final(i).GT.0.)) then
           inc_q_new = inc_q_new*0.1
        else
           exit
        end if
     end do ENDLESS
     q_prepp(i) = q_final(i) + inc_q_new
     q_prepm(i) = q_final(i) - inc_q_new
     if (q_final(i) .EQ. 0.) q_prepm(i) = 0.
     call tau_calc_fap(z_final,T_final,p_final,q_prepp,lwc_final,mu,&
          coeff,c_stp,tau_final)   ! c_stp: 'R98' or 'L93'

     call tb_calc_pl(T_final,tau_final,mu,TB_final)
     TB_bf = (TB_final - bias(1,:))/bias(2,:)
     TB_p = TB_bf

     TB_final = 0.
     tau_final = 0.
     call tau_calc_fap(z_final,T_final,p_final,q_prepm,lwc_final,mu,&
          coeff,c_stp,tau_final)   ! c_stp: 'R98' or 'L93'
     call tb_calc_pl(T_final,tau_final,mu,TB_final)
     TB_bf = (TB_final - bias(1,:))/bias(2,:)
     TB_m = TB_bf

     diff(n_Tq+i,:) = (TB_p - TB_m)/(2.*inc_q_new)
     if (q_final(i) .EQ. 0) diff(n_Tq+i,:) = (TB_p - TB_m)/(inc_q_new)
  end do Q_sens

  if (cloud .GT. 0) then
     where(lwc_final.NE.0)
        lwc_final = (10.*LOG10(lwc_final*1e3))
     end where
     ! * LWC sensitivity
     TB_final = 0.
     j_inc = n_Tq_2+1

     LWC_sens: do i=1,n_Tq
        TB_final = 0.
        tau_final = 0.

        if (lwc_final(i) .EQ. 0) cycle LWC_sens
        lwc_prepp = lwc_final
        lwc_prepm = lwc_final

        lwc_prepp(i) = lwc_final(i) + inc(j_inc)
        lwc_prepm(i) = lwc_final(i) - inc(j_inc)

        where(lwc_prepp.NE.0) lwc_prepp = (10.**(lwc_prepp/10.))*1e-3
        call tau_calc_fap(z_final,T_final,p_final,q_final,lwc_prepp,mu,&
             coeff,c_stp,tau_final)    ! c_stp: 'R98' or 'L93'
        call tb_calc_pl(T_final,tau_final,mu,TB_final)
        TB_bf = (TB_final - bias(1,:))/bias(2,:)
        TB_p = TB_bf

        TB_final = 0.
        tau_final = 0.
        where(lwc_prepm.NE.0) lwc_prepm = (10.**(lwc_prepm/10.))*1e-3
        call tau_calc_fap(z_final,T_final,p_final,q_final,lwc_prepm,mu,&
             coeff,c_stp,tau_final)
        call tb_calc_pl(T_final,tau_final,mu,TB_final)
        TB_bf = (TB_final - bias(1,:))/bias(2,:)
        TB_m = TB_bf

        diff(j_inc,:) = (TB_p - TB_m)/(2.*inc(j_inc))
        j_inc = j_inc + 1
     end do LWC_sens
     lwc = lwc*1000.
  end if
end subroutine make_K3_fap
