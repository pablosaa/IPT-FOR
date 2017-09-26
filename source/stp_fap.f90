subroutine stp_fap(c_stp,z_lr_lwc,x_fg,h_true,TB_n,z_final,p_rs,coeff_fap,bias_coeff_fap)

  use variables, only : theta, freq, n_freq, deltaz_low, n_Tq, n_Tq_2, pi
  use absorption
  implicit none
  character(len=*), intent(in) :: c_stp
  integer, intent(in) :: h_true
  real(kind=8), intent(in) :: x_fg(:)
  real, intent(in) :: z_lr_lwc(:),z_final(:),p_rs(:),coeff_fap(:,:),bias_coeff_fap(:,:)
  real(kind=8), intent(out) :: TB_n(:)

  ! Local variables
  integer :: i, i_lwc, j, k
  real(kind=8) :: mu, TB_bf(size(TB_n)), tau(n_Tq-1,n_freq)
  real(kind=8), dimension(h_true) :: lwc
  real(kind=8), dimension(n_Tq) :: lwc_final, T_final, q_final
  real, dimension(size(p_rs)) :: p_final

  mu = COS(theta*pi/180)
  T_final = x_fg(1:n_Tq)
  q_final = x_fg(n_Tq+1:n_Tq_2)
  p_final = p_rs                 ! REAL(p_rs, KIND=8)
  lwc_final = 0.
  TB_n = 0.
  tau = 0.
  if (h_true .NE. 0) then    ! If h_true equal to zero, then cloud-free scene.
     i_lwc = 1
     lwc = 10.**(x_fg(n_Tq_2+1:n_Tq_2+h_true)/10.d0)
     lwc = lwc/1000.d0          ! lwc in Kg/m3
     LWCF: do i=1,n_Tq
        if(z_lr_lwc(i_lwc) .EQ. (z_final(i)+deltaz_low/2.) ) then
           lwc_final(i) = lwc(i_lwc)
           i_lwc = i_lwc + 1
           if (i_lwc .GT. h_true) EXIT LWCF
        end if
     end do LWCF
  end if

  !++++ Radiative transfer. ++++
  select case (c_stp)
  case ('L93')
     ! processes L93
     call tau_calc_fap(z_final,T_final,p_final,q_final,lwc_final,mu,coeff_fap,'L93',tau)
  case ('R98')
     ! processes R98  (without the last 4 arguments!!!)
     call tau_calc_fap(z_final,T_final,p_final,q_final,lwc_final,mu,coeff_fap,'R98',tau)
  case default
     STOP 'Error selecting case R98 or L93 in STP_FAP subroutine.'
  end select

  !TB_n = tb_calc_pl(T_final,tau,mu)
  call tb_calc_pl(T_final,tau,mu,TB_n)

  ! calculate Bias-free fap predictor
  if (c_stp .EQ. 'L93') then
     TB_bf = (TB_n - bias_coeff_fap(1,:))/bias_coeff_fap(2,:)
     TB_n = TB_bf
  end if
  ! ++++ lwc -> back to g/m3
!  lwc = lwc*1000.
  return
end subroutine stp_fap

