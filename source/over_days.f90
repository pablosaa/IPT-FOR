!     ++++++++++++++++++
!     Begin Loop over days
!     ++++++++++++++++++

subroutine over_days(date, cov_rs, Sap_a,Se_a,Se_fap_org,&
     c_a,d_a,lwc_ap_a,coeff_fap,bias_coeff_fap)
  use variables
  use maths
  IMPLICIT NONE
  character(len=*), intent(in), dimension(17) :: date
  real(kind=8), intent(in), dimension(n_Tq_2,n_Tq_2) :: cov_rs
  real(kind=8), intent(in) :: Sap_a(wd,n_Tq_2+pr,n_Tq_2+pr),&
       Se_a(wd,n_freq_r+pr,n_freq_r+pr)
  real(kind=8), intent(in) :: c_a(wd,pr), d_a(wd,pr),lwc_ap_a(wd,pr)
  real, intent(in) :: Se_fap_org(n_freq_r,n_freq_r), coeff_fap(10,n_freq_r), &
       bias_coeff_fap(2,n_freq_r)


  !     Local variables
  ! Integers
  integer :: i,j,k,jh,nn,ss,i_days, k_pf
  integer :: profs, nd_rad, ip, it, allost
  integer :: h_true, i_neg, i_ssat, cloud
  integer :: kf, i_qn, chi_check, st_1, st_2, TB_faulty
  integer :: error_status
  ! Logicals
  logical :: REP_PROF
  logical :: SKIP_PROF
  logical :: SKIP_PROF_ICY
  ! Reals
  real :: ptime_s, ptime_e
  real :: T_gr, tau, faktor, q_gr
  real, parameter, dimension(n_Tq) :: z_final=&
       (/ (i,i=0,5000,250), (i,i=5500,10000,500), (i,i=15000,30000,5000) /)
  real, dimension(z_level) :: z_rad_liquid
  real, dimension(n_Tq) :: T_rs, q_rs, p_rs
  real(kind=8) :: d_i2
  real(kind=8) :: chi_1, chi_2, chi_3, chi_4, chi_test_x, chi_test_y
  real(kind=8), dimension(n_Tq) :: q_check, T_check, e_check, es_check, rh_check
  real(kind=8), dimension(n_freq) :: TB_n
  real(kind=8) :: deltaz_hr

  !     Allocatable Vectors
  integer, allocatable, dimension(:) :: ind, h_true_m, h_m
  real, allocatable, dimension(:) :: garbage, z_hr
  real, allocatable, dimension(:) :: nd_rad_x, time_a, tem_a
  real, allocatable, dimension(:) :: dpt_a, lwp_mic_a, base_hr_a
  real, allocatable, dimension(:) :: iwv_mic_a, code_c_a
  real, allocatable, dimension(:) :: top_m, base_m
  real, allocatable, dimension(:) :: dbz, z_lr_lwc, qc, e_qc, inc
  real, allocatable, dimension(:) :: dbz_hr,vel_hr,class_hr, y_gen
  real(kind=8), allocatable, dimension(:) :: garbage2, c, d, lwc_ap, x_fg, y_n
  real(kind=8), allocatable, dimension(:) :: x_n, x_ap, y_ap, x_n_s, x_ap_s, xd, bb
  !     --
  !     Allocatable Matrices
  real, allocatable, dimension(:,:) :: tri
  real, allocatable, dimension(:,:) :: TB_a, dbz_a, vel_a, sigma_a
  real, allocatable, dimension(:,:) :: code_a, T_rs_a 
  real, allocatable, dimension(:,:) :: p_rs_a, q_rs_a
  real, allocatable, dimension(:,:) :: dbz_hr_a, lin_hr_a, vel_hr_a
  real, allocatable, dimension(:,:) :: class_hr_a, base_m_a, top_m_a
  real, allocatable, dimension(:,:) :: dd_a, vv_a
  real(kind=8), allocatable, dimension(:,:) :: garbage3, Sap, Se, Se_fap, K_J, K3
  real(kind=8), allocatable, dimension(:,:) :: tK_J, Sd_m, Sd_m_inv, cc, Sdy, Sdy_inv
  real(kind=8), allocatable, dimension(:) :: dy
  ! Allocate variables for SOLUTION CHECKS part
  real(kind=8), allocatable, dimension(:) :: dF_s, dx
  real(kind=8), allocatable, dimension(:,:) :: Sd_s, Sd_s_inv, S_deyd, S_deyd_inv
  ! Allocate variables for THEORETICAL OPTIMAL ESTIMATION
  real(kind=8), allocatable, dimension(:) :: err_oe
  real(kind=8), allocatable, dimension(:,:) :: K3_e, K_e, tK_e, Sx, Sx_inv, Sd

  ! Strings
  character(len=70) :: infile, outfile
  character(len=8) :: d_today
  character(len=10) :: d_now
  character(len=9) :: fudge_cor
  character(len=7) :: mode, status


  interface
     subroutine de_attenuate(c,d,dbz,dbz_da,h_true,lwp_mic,&
          T_rs,q_rs,p_rs,z_final,base_m,top_m,dbz_hr,dbz_hr_da,z_hr,&
          h_true_m,nn,z_lr_lwc)
       real(kind=8), dimension(:),intent(in) :: c, d
       real, dimension(:), intent(out) :: dbz_da, dbz_hr_da
       integer, intent(in) :: h_true, nn
       real, intent(in) :: lwp_mic
       real, dimension(:), intent(in) :: dbz, T_rs, q_rs, p_rs, z_final
       real, dimension(:), intent(in) :: base_m, top_m, dbz_hr
       real, dimension(:), intent(in) :: z_hr, z_lr_lwc
       integer, dimension(:), intent(in) :: h_true_m
     end subroutine de_attenuate

     subroutine stp_fap(c_stp,z_lr_lwc,x_fg,h_true,TB_n,z_final,p_rs,coeff_fap,bias_coeff_fap)
       character(len=*), intent(in) :: c_stp
       integer, intent(in) :: h_true
       real(kind=8), intent(in) :: x_fg(:)
       real, intent(in) :: z_lr_lwc(:),z_final(:),p_rs(:),coeff_fap(:,:),bias_coeff_fap(:,:)
       real(kind=8), intent(out) :: TB_n(:)
     end subroutine stp_fap

     subroutine make_K3_fap(c_stp, z_lr, x_n, z_final, p_rs, h_true, inc, coeff, bias, cloud, diff)
       character(len=*), intent(in) :: c_stp
       integer, intent(in) :: h_true, cloud
       real(kind=8), intent(in) :: x_n(:)
       real, intent(in) :: z_lr(:),z_final(:),p_rs(:),inc(:),coeff(:,:),bias(:,:)
       real(kind=8), intent(out) :: diff(:,:)
     end subroutine make_K3_fap

     subroutine lwc_high_res(z_lr_lwc, h_true_m, lwc, xn, c, d, dbz_hr, z_hr, n_cl)
       integer, intent(in) :: h_true_m(:), n_cl
       real, intent(in) :: z_lr_lwc(:), dbz_hr(:), z_hr(:)
       real(kind=8), intent(in) :: lwc(:), xn(:), c(:), d(:)

     end subroutine lwc_high_res

  end interface

  fudge_sap = 1.
  fudge_sap_lwc = 1.
  fudge_sap_T = 1.
  fudge_sap_q = 1.
  fudge_se = 1.
  fudge_se_tb = 1.

  i = 0

  OVERDAYS: do i_days = 1,nt_days

     infile = path_in//'input_bbc2/'//date(i_days)//'r_bbc2_com.dat'
     outfile = path_out//'out_ipt_kn4_'//date(i_days)//'r.dat'

     call open_binary(infile, 12, 2)
     read(12,rec=1) profs,nd_rad
     close(12)

     ip = 61
     allocate(z_hr(nd_rad), tri(profs,261+ip+4*(nd_rad)+12))
     z_rad_liquid = -99.
     z_hr = -99. ; tri = -99.
     call open_binary(infile, 12, 50000000)
     read(12,rec=1) profs,nd_rad,z_rad_liquid,z_hr,tri
     close(12)

     write(*,*) profs, nd_rad      

     !     Start to allocate all variables
     allocate(nd_rad_x(profs), time_a(profs), TB_a(profs,25), &
          dbz_a(profs,z_level), vel_a(profs,z_level), &
          sigma_a(profs,z_level), code_a(profs,8), &
          tem_a(profs), dpt_a(profs), lwp_mic_a(profs), &
          T_rs_a(profs,35), p_rs_a(profs,35), q_rs_a(profs,35), &
          dbz_hr_a(profs,nd_rad), lin_hr_a(profs,nd_rad), &
          vel_hr_a(profs,nd_rad), class_hr_a(profs,nd_rad), &
          base_m_a(profs,5),top_m_a(profs,5), &
          base_hr_a(profs), iwv_mic_a(profs), code_c_a(profs), &
          dd_a(profs,35), vv_a(profs,35), stat=allost)

     if(allost /= 0) stop 'Error while allocating the firsts variables...!!!'

     nd_rad_x = tri(:,1)
     time_a = tri(:,2)
     TB_a = tri(:,5:29)
     it = 30
     dbz_a = tri(:,it:it+z_level-1)
     it = it + z_level
     vel_a = tri(:,it:it+z_level-1)
     it = it + z_level
     sigma_a = tri(:,it:it+z_level-1)
     code_a = tri(:,74+ip:81+ip) ! OJO CON IP....
     tem_a = tri(:,82+ip)
     dpt_a = tri(:,83+ip)
     lwp_mic_a = tri(:,84+ip)
     T_rs_a = tri(:,85+ip:119+ip)
     p_rs_a = tri(:,120+ip:154+ip)
     q_rs_a = tri(:,155+ip:189+ip)
     dbz_hr_a = tri(:,190+ip:190+ip+nd_rad-1)
     lin_hr_a = tri(:,190+ip+nd_rad:190+ip+2*nd_rad-1)
     vel_hr_a = tri(:,190+ip+2*nd_rad:190+ip+3*nd_rad-1)
     class_hr_a = tri(:,190+ip+3*nd_rad:190+ip+4*nd_rad-1)
     base_m_a = tri(:,190+ip+4*nd_rad:190+ip+4*nd_rad+4)
     top_m_a = tri(:,190+ip+4*nd_rad+5:190+ip+4*nd_rad+9)
     base_hr_a = tri(:,190+ip+4*nd_rad+10)
     iwv_mic_a = tri(:,190+ip+4*nd_rad+11)
     i = 190+ip+4*nd_rad+12
     dd_a = tri(:,i:i+34)
     vv_a = tri(:,i+35:i+69)
     code_c_a = tri(:,i+70)

     !     Number of applicable cases (cloudy & cloud-free)
     nn=count((sum(code_a(:,1:7),DIM=2).EQ.7) .OR. (code_c_a.EQ.1))
     write(*,*) 'Number of applicable cases on ', date(i_days),nn

     z_rad_liquid = 1000.*z_rad_liquid
     base_m_a = 1000.*base_m_a
     base_hr_a = 1000.*base_hr_a
     top_m_a = 1000.*top_m_a
     z_hr = 1000.*z_hr

     OPEN(UNIT=13,FILE=outfile,STATUS='unknown',ACTION='write',IOSTAT=allost)
     if(allost/=0) STOP 'Error opening output file.'

     write(13,'(A,A)') 'INFILE: ', infile
     call DATE_AND_TIME(d_today,d_now)
     write(13,'(A,A,''.'',A,''.'',A,4X,A,'':'',A,'':'',A)') 'Start: ', d_today(1:4), d_today(5:6),d_today(7:8),&
          d_now(1:2), d_now(3:4), d_now(5:10)
     write(13,'(A,A,I8)') 'Number of applicable cases on ', date(i_days), nn

     !     +++
     !     BEGIN LOOP OVER PROFILES
     !     ---
     PROFILE: do k_pf = 1,profs,2 ! Profile loop

        fudge_cor  = 'NO_FUDGE'
        REPEAT_PROF: do
           REP_PROF = .FALSE.
           SKIP_PROF = .FALSE.
           mode = 'XX'
           ab = 0 ; mat_inv = -99 ; cloud_prop = -99
           h_true = -99 ; i_neg = -99 ; i_ssat = -99

           T_rs = T_rs_a(k_pf,:)
           p_rs = p_rs_a(k_pf,:)*100.
           q_rs = q_rs_a(k_pf,:)

           !     ++++ Applicable cases ?

           if ((sum(code_a(k_pf,1:7)).LT.7) .AND. (code_c_a(k_pf).EQ.0)) then
              status = 'NO_APP'
              print *,'DATE/TIME: ', date(i_days),' ',time_a(k_pf)
              print *,'Not applicable.'
              ! start to write in output file
              write(13,'(A8,F10.6)') date(i_days),time_a(k_pf)
              write(13,*) mode
              write(13,*) status
              write(13,*) fudge_cor
              write(13,'(10G15.7)') TRANSPOSE(code_a(k_pf:k_pf,:))
              write(13,'(20G15.7)') T_rs_a(k_pf,:)
              write(13,'(20G15.7)') q_rs_a(k_pf,:)
              write(13,'(20G15.7)') p_rs_a(k_pf,:)
              write(13,'(20G15.7)') z_rad_liquid
              write(13,'(I6)') nd_rad
              write(13,'(20G15.7)') z_hr
              write(13,'(20G15.7)') dbz_hr_a(k_pf,:)
              write(13,'(20G15.7)') vel_hr_a(k_pf,:)
              write(13,'(20G15.7)') class_hr_a(k_pf,:)
              write(13,'(3G15.7)') lwp_mic_a(k_pf), iwv_mic_a(k_pf), base_hr_a(k_pf)
              ! end to write in output file
              cycle PROFILE
           endif
           status = ''
           ! Applicable cloudy ?
           APP_CLOUDY: if ((sum(code_a(k_pf,1:7)) .EQ. 7) .AND. &
                (time_a(k_pf) .GT. time_crit)) then
              cloud = 1              ! STP parameter
              ! abort parameters
              kf = -99 ; i_qn = -99 ; chi_1 = -99. ; chi_2 = -99. ; chi_3 = -99.
              chi_4 = -99. ; chi_test_x = -99. ; chi_test_y = -99. ; chi_check = -99
              st_1 = -99 ; st_2 = -99 ; it = -99 ; TB_faulty = -99
              ! Output info
              print *,'DATE/TIME: ', date(i_days),' ',time_a(k_pf)
              print *,'Cloudy'
              mode = 'CLOUDY'
              call cpu_time(ptime_s)
              ! base / top determination
              i = 0 ; j = 0 ; nn = 0
              nn = count((top_m_a(k_pf,:).GT.0 .AND. base_m_a(k_pf,:).GE.0) &
                   .AND. (top_m_a(k_pf,:).GT.base_m_a(k_pf,:) &
                   .AND. top_m_a(k_pf,:).LE.5000)) 

              BASE_TOP: if (nn.EQ.0) then    !(nn.LT.2) then  !     
                 write(*,*) 'Base/Top not applicable...!!!'
                 ab = 1 ; cloud_prop = 1
                 SKIP_PROF = .TRUE.
              else
                 allocate(top_m(nn), base_m(nn), h_true_m(nn), h_m(nn), stat=allost)
                 if(allost/=0) stop 'Allocate top_m, base_m, h_true_m failed...!!'
                 write(*,*) 'Number of cloud layers: ', nn

                 do i=1,5
                    if((top_m_a(k_pf,i).GT.0 .AND. base_m_a(k_pf,i).GE.0) &
                         .AND. (top_m_a(k_pf,i).GT.base_m_a(k_pf,i) .AND. &
                         top_m_a(k_pf,i).LE.5000)) then
                       j = j + 1
                       top_m(j) = top_m_a(k_pf,i) ; base_m(j) = base_m_a(k_pf,i)
                    end if
                 end do

                 h_true_m = int((top_m - base_m)/deltaz_low)
                 h_m = h_true_m
                 do k = 1,nn
                    i = minloc(z_rad_liquid,DIM=1, &
                         MASK=((z_rad_liquid - deltaz_low/2) .EQ. base_m(k)))
                    j = minloc(z_rad_liquid,DIM=1, &
                         MASK=((z_rad_liquid + deltaz_low/2) .EQ. top_m(k)))
                    ss = size(dbz_a(k_pf,i:j))
                    if (k == 1) then
                       allocate(dbz(ss))
                       dbz = dbz_a(k_pf,i:j)
                    else
                       allocate(garbage(size(dbz)))
                       garbage = dbz
                       deallocate(dbz)
                       allocate(dbz(size(garbage) + ss))
                       dbz = [garbage, dbz_a(k_pf,i:j)]
                       deallocate(garbage)
                    end if
                 end do
                 where(h_true_m .EQ. 1) h_m = 2
                 h_m = h_m + 1
                 h_true = size(dbz)

                 !     ++++ Calculate low-resolution height on which lwc is to be determined
                 allocate(z_lr_lwc(h_true), stat=allost)
                 if(allost/=0) stop 'Allocate z_lr_wc failed...!!'
                 z_lr_lwc = -99.
                 k = 1
                 do i = 1,nn
                    do j = 1,h_true_m(i)
                       z_lr_lwc(k) = base_m(i) + deltaz_low*(REAL(j)-.5)
                       k = k + 1
                    end do
                 end do

                 !     ++++ Check vertical cloud thickness
                 V_CLOUD: if (count(h_true_m .GT. pr) .GT. 0) then
                    ab = 1 ; cloud_prop = 1
                    print *, 'Vertical cloud thickness failed.'
                    SKIP_PROF = .TRUE.
                 else
                    !     ++++ specify T, q profiles to use
                    !     ++++ Define ground measurements
                    T_gr = tem_a(k_pf)
                    tau = dpt_a(k_pf)
                    faktor = L/(T0*Rw)
                    q_gr = (e0*EXP(faktor*(tau - T0)/tau))/(T_gr*Rw)
                    !     ++++ Calculate q saturation within cloud boundaries
                    allocate(qc(h_true), e_qc(h_true), ind(h_true), stat=allost)
                    if(allost/=0) stop 'Allocate qc, e_qc, ind failed...!!'
                    qc = -99. ; e_qc = -99. ; ind = -99
                    ! reference here: i = i_lwc ; j = multi ; k = level_count
                    i = 1 ; j = 1 ; k = 0
                    JHDO: do jh = 1,n_Tq
                       if (base_m(j).LE.z_final(jh)) then
                          qc(i) = (e0/(Rw*T_rs(jh)))* &
                               (EXP(faktor*(T_rs(jh) - T0)/T_rs(jh)))

                          !    Gaussian error propagation
                          !e_Tc = SQRT(cov_rs(jh,jh))
                          e_qc(i) = (e0/(Rw*(T_rs(jh))**2))* &
                               (EXP(faktor*(T_rs(jh) - T0)/T_rs(jh)))* &
                               (L/(Rw*T_rs(jh))-1.)*SQRT(cov_rs(jh,jh))
                          ind(i) = jh ; i =i + 1 ; k = k + 1
                          if (k .EQ. h_true_m(j)) then
                             j = j + 1 ; k = 0
                             if (j.GT.nn) exit JHDO     !(j.EQ.nn+1) 
                          end if
                       end if
                    end do JHDO

                    ! ++++ define
                    do i = 1,nn
                       j = h_m(i)-2 ; k = h_true_m(i)
                       ss = size(c_a(j,1:k))                  
                       if (i==1) then                     
                          allocate(c(ss),d(ss),lwc_ap(ss))
                          c = c_a(j,1:k)
                          d = d_a(j,1:k)
                          lwc_ap = lwc_ap_a(j,1:k)
                       else
                          ! Matrix c
                          allocate(garbage2(size(c)))
                          garbage2 = c
                          deallocate(c)
                          allocate(c(size(garbage2)+ss))
                          c = [garbage2,c_a(j,1:k)]
                          !  Matrix d
                          garbage2 = d
                          deallocate(d)
                          allocate(d(size(garbage2)+ss))
                          d = [garbage2,d_a(j,1:k)]
                          ! Matrix lwc_ap
                          garbage2 = lwc_ap
                          deallocate(lwc_ap)
                          allocate(lwc_ap(size(garbage2)+ss))
                          lwc_ap = [garbage2,lwc_ap_a(j,1:k)]
                          deallocate(garbage2)
                       end if
                    end do     ! end define

                    ! ++++  Define Sap
                    allocate(Sap(n_Tq_new_2+h_true, n_Tq_new_2+h_true), stat=allost)
                    if(allost/=0) stop 'Allocate Sap matrix failed...!!'
                    Sap = 0.d0
                    ! Define Sap_eof due to RS "a priori"
                    Sap(1:n_Tq_new, 1:n_Tq_new) = cov_rs(1:n_Tq_new, 1:n_Tq_new)
                    Sap(n_Tq_new+1:n_Tq_new_2,1:n_Tq_new) = cov_rs(n_Tq+1:n_Tq+n_Tq_new,1:n_Tq_new)
                    Sap(n_Tq_new+1:n_Tq_new_2,n_Tq_new+1:n_Tq_new_2) = cov_rs(n_Tq+1:n_Tq+n_Tq_new,n_Tq+1:n_Tq+n_Tq_new)
                    Sap(1:n_Tq_new, n_Tq_new+1:n_Tq_new_2) = cov_rs(1:n_Tq_new, n_Tq+1:n_Tq+n_Tq_new)


                    ! Consider cloud covariances
                    jh = n_Tq_new_2
                    do i = 1,nn
                       Sap(jh+1:jh+h_true_m(i),jh+1:jh+h_true_m(i)) = &
                            Sap_a(h_m(i)-2, n_Tq_2+1:n_Tq_2+h_true_m(i), n_Tq_2+1:n_Tq_2+h_true_m(i))
                       jh = jh + h_true_m(i)
                    end do
                    Sap = fudge_sap*Sap
                    Sap(n_Tq_new_2+1:n_Tq_new_2+h_true, &
                         n_Tq_new_2+1:n_Tq_new_2+h_true) = fudge_sap_lwc* &
                         Sap(n_Tq_new_2+1:n_Tq_new_2+h_true, &
                         n_Tq_new_2+1:n_Tq_new_2+h_true)

                    Sap(1:n_Tq_new,1:n_Tq_new) = fudge_sap_T* &
                         Sap(1:n_Tq_new,1:n_Tq_new)

                    Sap(n_Tq_new+1:n_Tq_new_2,n_Tq_new+1:n_Tq_new_2) = &
                         fudge_sap_q*Sap(n_Tq_new+1:n_Tq_new_2, &
                         n_Tq_new+1:n_Tq_new_2)

                    ! Define increments for Jacobian
                    allocate(inc(n_Tq_2+h_true), stat=allost)
                    if(allost/=0) stop 'Allocate inc failed...!!'
                    j = 1 ; k = n_Tq_new_2 + 1 
                    do i =1,n_Tq_2
                       if(cov_rs(i,i).NE.0) j=i
                       inc(i) = SQRT(cov_rs(j,j))/10.
                    end do
                    do i = n_Tq_2+1,n_Tq_2+h_true
                       inc(i) = SQRT(Sap(k,k))/10.
                       k = k + 1
                    end do

                    ! ++++ Determine Se
                    jh = MAXVAL(h_true_m)
                    if (jh.EQ.1) jh = 2
                    jh = jh + 1
                    select case (use_dbz)
                    case ('NO')
                       allocate(garbage3(n_freq_r+2+h_true,n_freq_r+2+h_true))
                       garbage3 = 0.d0
                       garbage3(1:n_freq_r,1:n_freq_r) = &
                            Se_a(jh-2,1:n_freq_r,1:n_freq_r)
                       garbage3(n_freq_r+1,n_freq_r+1) = e_T_gr**2
                       garbage3(n_freq_r+2,n_freq_r+2) = e_q_gr**2
                       do i=1,h_true
                          garbage3(n_freq_r+2+i,n_freq_r+2+i) = e_qc(i)**2
                       end do
                       allocate(Se(n_freq+2+h_true,n_freq+2+h_true))
                       Se = 0.d0
                       forall(j=1:n_freq,i=1:n_freq,.TRUE.)
                          Se(i,j) = garbage3(freq_num(i),freq_num(j))
                       end forall
                       Se(n_freq+1:n_freq+2+h_true, n_freq+1:n_freq+2+h_true) = &
                            garbage3(n_freq_r+1:n_freq_r+2+h_true, &
                            n_freq_r+1:n_freq_r+2+h_true)
                       deallocate(garbage3)
                       allocate(Se_fap(n_freq+2+h_true,n_freq+2+h_true), stat=allost)
                       if(allost/=0) stop 'Allocate Se_fap matrix failed...!!'
                       Se_fap = 0.d0
                    case ('YES')
                       allocate(garbage3(n_freq_r+h_true+2+h_true,n_freq_r+h_true+2+h_true))
                       garbage3 = 0.d0
                       k = n_freq_r         ! + 1
                       j = k
                       garbage3(1:n_freq_r,1:n_freq_r) = &
                            Se_a(jh-2,1:n_freq_r,1:n_freq_r)
                       do i = 1,nn
                          garbage3(k+1:k+h_true_m(i), k+1:k+h_true_m(i)) = &
                               Se_a(h_m(i)-2,j+1:j+h_true_m(i),j+1:j+h_true_m(i))
                          k = k + h_true_m(i)
                       end do

                       garbage3(n_freq_r+h_true+1,n_freq_r+h_true+1) = e_T_gr**2
                       garbage3(n_freq_r+h_true+2,n_freq_r+h_true+2) = e_q_gr**2
                       do i=1,h_true
                          garbage3(n_freq_r+h_true+2+i,n_freq_r+h_true+2+i)= e_qc(i)**2
                       end do
                       allocate(Se(n_freq+h_true+2+h_true,n_freq+h_true+2+h_true), stat=allost)
                       if(allost/=0) stop 'Allocate Se matrix failed...!!'
                       Se = 0.d0
                       forall(j=1:n_freq,i=1:n_freq,.TRUE.)
                          Se(i,j) = garbage3(freq_num(i),freq_num(j))
                       end forall
                       Se(n_freq+1:n_freq+h_true+2+h_true,n_freq+1:n_freq+h_true+2+h_true)=&
                            garbage3(n_freq_r+1:n_freq_r+h_true+2+h_true, &
                            n_freq_r+1:n_freq_r+h_true+2+h_true)
                       deallocate(garbage3)
                       allocate(Se_fap(n_freq+h_true+2+h_true,n_freq+h_true+2+h_true), stat=allost)
                       if(allost/=0) stop 'Allocate Se_fap matrix failed...!!'
                       Se_fap = 0.d0
                    case default
                       stop 'Error selecting case USE_DBZ.'
                    end select

                    Se_fap = Se
                    Se_fap(1:n_freq,1:n_freq)=Se(1:n_freq,1:n_freq) + Se_fap_org
                    Se_fap = fudge_se * Se_fap
                    Se_fap(1:n_freq,1:n_freq)=fudge_se_tb * Se_fap(1:n_freq,1:n_freq)

                    !  ++++ Invert matrix Se_fap
                    ss = size(Se_fap,DIM=1)
                    allocate(garbage3(ss,ss))
                    call invert_LU(Se_fap, garbage3, ss, error_status)
                    deallocate(garbage3)
                    if (error_status .EQ. 1) then   ! Invert matrix Se_fap
                       ab = 1 ; mat_inv = 1
                       print *, 'Out from inversion of Se_fap matrix.'
                       SKIP_PROF = .TRUE.
                       GOTO 666
                    end if
                    !  ++++ Invert matrix Sap
                    ss = size(Sap,DIM=1)
                    allocate(garbage3(ss,ss))
                    call invert_LU(Sap, garbage3, ss, error_status)
                    deallocate(garbage3)
                    if (error_status .EQ. 1) then   !! .EQ. 1 ! Invert matrix Sap
                       ab = 1 ; mat_inv = 2
                       print *, 'Out from inversion of Sap matrix.'
                       SKIP_PROF = .TRUE.
                       GOTO 666
                    end if
                    !  ++++ de-attenuate dbz considering q_rs and approx. of lwc
                    allocate(garbage(h_true), dbz_hr(nd_rad))
                    call de_attenuate(c,d,dbz,garbage,h_true, &
                         lwp_mic_a(k_pf),T_rs,q_rs,p_rs,z_final, &
                         base_m,top_m, dbz_hr_a(k_pf,:),dbz_hr,z_hr, &
                         h_true_m,nn,z_lr_lwc)
                    dbz = garbage
                    deallocate(garbage)

                    ! ++++ Initialization
                    select case (use_dbz)
                    case ('YES')
                       allocate(y_gen(n_freq + size(dbz) + 2 + h_true))
                       y_gen = [TB_a(k_pf, freq_num), dbz, T_gr, q_gr, qc]
                    case ('NO')
                       allocate(y_gen(n_freq + 2 + h_true))
                       y_gen = [TB_a(k_pf, freq_num), T_gr, q_gr, qc]
                    case default
                       STOP 'Error selecting case for y_gen.'
                    end select

                    TB_faulty = 0
                    if (count(TB_a(k_pf, freq_num).LT.0) .GT. 0 ) then
                       ab = 1 ; TB_faulty = 1
                       print *, 'NOT first guess.'
                       SKIP_PROF = .TRUE.
                       GOTO 666
                    end if
                    ! ++++ first guess
                    allocate(x_fg(n_Tq_2 + size(lwc_ap)))
                    x_fg = [real(T_rs,KIND=8), real(q_rs,KIND=8), lwc_ap]   ! First guess profile.

                    ! ++++ calculate y_n = F(x_n)

                    select case (abs)
                    case ('r98')
                       call stp_fap('R98',z_lr_lwc,x_fg,h_true,TB_n,z_final,p_rs,coeff_fap,&
                            bias_coeff_fap)
                    case ('l93')
                       call stp_fap('L93',z_lr_lwc,x_fg,h_true,TB_n,z_final,p_rs,coeff_fap,&
                            bias_coeff_fap)
                    case default
                       stop 'Error selecting ABS for y_n = F(x_n).!!!'
                    end select

                    ! +++++++ Assigning y_n and Dimension to  K_J and tK_J matrices
                    select case (use_dbz)
                    case ('YES')
                       !*Reference:  y_n = [TB_n, dbz, T_gr_n, q_gr_n, qc_n]  ;  T_gr_n = x_fg(1)
                       !             q_gr_n = x_fg(n_Tq+1)       ;     qc_n = x_fg(n_Tq+ind)
                       ss = size(c,DIM=1)
                       allocate(garbage2(ss))
                       garbage2 = c + d*x_fg(n_Tq_2+1:n_Tq_2+h_true)
                       allocate(y_n(n_freq+ss+2+h_true), K_J(n_Tq_new_2+h_true, n_freq+h_true+2+h_true), &
                            tK_J(n_freq+h_true+2+h_true, n_Tq_new_2+h_true),stat=allost)
                       if(allost/=0) stop 'Allocate K_J, tK_J failed...!!'
                       y_n = [TB_n, garbage2, x_fg(1), x_fg(n_Tq+1), x_fg(n_Tq+ind)]
                       K_J = 0.d0
                       deallocate(garbage2)
                    case ('NO')
                       allocate(y_n(n_freq+2+h_true), K_J(n_Tq_new_2+h_true, n_freq+2+h_true), &
                            tK_J(n_freq+2+h_true, n_Tq_new_2+h_true))
                       y_n = [TB_n, x_fg(1), x_fg(n_Tq+1), x_fg(n_Tq+ind)]
                       K_J = 0.d0
                    case default
                       stop 'Error selecting vector: y_n and K_J.'
                    end select

                    allocate(x_n(size(x_fg)), stat=allost)
                    if(allost/=0) stop 'Allocate x_n failed...!!'
                    x_n = x_fg

                    ! +++ connocate a priori.
                    allocate(x_ap(n_Tq_2 + size(lwc_ap)), y_ap(size(y_n)), stat=allost)
                    if(allost/=0) stop 'Allocate x_ap, y_ap failed...!!'
                    x_ap = [real(T_rs,KIND=8), real(q_rs,KIND=8), lwc_ap]
                    y_ap = y_n

                    ! +++++++++++++++++++++
                    ! 4.) BEGIN ITERATION
                    ! +++++++++++++++++++++
                    ITER_IT: do it = 0, 9
                       ! +++ Cost function check.
                       kf = -99
                       ! +++ q negative check.
                       i_neg = -99
                       ! inversion checks.
                       st_1 = -99   ;  st_2 = -99

                       ! +++++++++++++++++++++++++++++++++++++++++
                       ! 5.) CALCULATE K_J BY COMBINING K1 AND K3
                       ! +++++++++++++++++++++++++++++++++++++++++
                       ! (K3: de/dTB, df/dTB, dlwc/dTB)
                       ! (K1: dLWC/ddbz)
                       allocate(K3(n_Tq_2 + h_true, n_freq),stat=allost)
                       if(allost/=0) stop 'Allocate K3 failed...!!'
                       K3 = 0.d0
                       select case (abs)
                       case ('r98')
                          call make_K3_fap('R98', z_lr_lwc, x_n, z_final, p_rs, h_true, inc, &
                               coeff_fap, bias_coeff_fap, cloud, K3)
                       case ('l93')
                          call make_K3_fap('L93', z_lr_lwc, x_n, z_final, p_rs, h_true, inc, &
                               coeff_fap, bias_coeff_fap, cloud, K3)
                       case default
                          stop 'Error selecting make_K3_fap subroutine.'
                       end select

                       select case (use_dbz)
                       case ('YES')
                          ! ++++ Make K1 subsection
                          ! ++ sensitivity d(dBZ)/d(LWC)
                          ! Here, temporal variable "garbage3" is "K1". Instant of MAKE_K1 subroutine.
                          allocate(garbage3(h_true, h_true)) ; garbage3 = 0.
                          forall(j=1:h_true,i=1:h_true,j.EQ.i)
                             !  1. Ableitung der Z/LWC Beziehung
                             !   K1(i,j) = d(i)*10.*1./(lwc(i)*log(10.))
                             garbage3(i,j) = d(i)              ! K1(i,j) = d(i)
                          end forall
                          K_J(1:n_Tq_new, 1:n_freq) = K3(1:n_Tq_new, 1:n_freq)
                          K_J(n_Tq_new+1:n_Tq_new_2, 1:n_freq) = K3(n_Tq+1:n_Tq+n_Tq_new, 1:n_freq)
                          K_J(n_Tq_new_2+1:n_Tq_new_2+h_true, 1:n_freq) = K3(n_Tq_2+1:n_Tq_2+h_true, 1:n_freq)
                          K_J(n_Tq_new_2+1:n_Tq_new_2+h_true, n_freq+1:n_freq+h_true) = garbage3
                          K_J(1,n_freq+h_true+1) = 1.d0          ! T_gr
                          K_J(n_Tq_new+1,n_freq+h_true+2) = 1.d0 ! q_gr
                          do j=1,h_true
                             K_J(n_Tq_new+ind(j), n_freq+h_true+2+j) = 1.d0
                          end do
                          deallocate(garbage3)
                       case ('NO')
                          K_J(1:n_Tq_new, 1:n_freq) = K3(1:n_Tq_new, 1:n_freq)
                          K_J(n_Tq_new+1:n_Tq_new_2, 1:n_freq) = K3(n_Tq+1:n_Tq+n_Tq_new, 1:n_freq)
                          K_J(n_Tq_new_2+1:n_Tq_new_2+h_true, 1:n_freq) = K3(n_Tq_2+1:n_Tq_2+h_true, 1:n_freq)
                          K_J(1, n_freq+1) = 1.d0           ! T_gr
                          K_J(n_Tq_new+1, n_freq+2) = 1.d0  ! q_gr
                          do j=1,h_true
                             K_J(n_Tq_new+ind(j), n_freq+2+j) = 1.d0
                          end do
                       case default
                          stop 'Error selecting use_dbz: YES / NO for K_J.'
                       end select
                       deallocate(K3)

                       ! +++++++++++++++++++++++++++
                       ! 6.) OPTIMAL ESTIMATION
                       ! +++++++++++++++++++++++++++
                       ! Rodgers (5.8 / 5.9)
                       allocate(x_n_s(n_Tq_new_2+h_true), x_ap_s(n_Tq_new_2+h_true), &
                            xd(n_Tq_new_2+h_true), stat=allost)
                       if(allost/=0) stop 'Allocate x_n_s, x_ap_s, xd failed...!!'
                       x_n_s = [x_n(1:n_Tq_new),x_n(n_Tq+1:n_Tq+n_Tq_new),x_n(n_Tq_2+1:n_Tq_2+h_true)]
                       x_ap_s = [x_ap(1:n_Tq_new),x_ap(n_Tq+1:n_Tq+n_Tq_new),x_ap(n_Tq_2+1:n_Tq_2+h_true)]
                       tK_J = transpose(K_J)
                       ss = size(Se_fap,DIM=1)
                       allocate(Sd_m(ss,ss), Sd_m_inv(ss,ss), stat=allost)
                       if (allost/=0) stop 'Allocate Sd_m failed...!!!'

                       Sd_m = Se_fap + matmul(tK_J, matmul(Sap,K_J))

                       call invert_LU(Sd_m, Sd_m_inv, ss, error_status)
                       if (error_status .EQ. 1) then   !! .EQ. 1 ! Invert matrix Sap
                          ab = 1 ; mat_inv = 3
                          deallocate(x_n_s, x_ap_s, xd, Sd_m, Sd_m_inv)
                          print *, 'Out from inversion of Sd_m matrix.'
                          SKIP_PROF = .TRUE.
                          GOTO 666    ! SKIP_PROF
                       end if

                       allocate(bb(size(y_gen)), cc(size(tK_J,DIM=1), size(Sap,DIM=2)), garbage2(size(x_n_s)),stat=allost)
                       if(allost/=0) stop 'Allocate bb, cc failed...!!'
                       garbage2 = x_n_s - x_ap_s
                       bb = y_gen - y_n + MATMUL(garbage2, K_J)
                       cc = MATMUL(tK_J, Sap)

                       xd = x_ap_s + MATMUL(bb, MATMUL(Sd_m_inv, cc))
                       deallocate(bb,cc, garbage2)
                       x_n = [xd(1:n_Tq_new), x_fg(n_Tq_new+1:n_Tq), xd(n_Tq_new+1:n_Tq_new_2), &
                            x_fg(n_Tq+n_Tq_new+1:n_Tq_2), xd(n_Tq_new_2+1:n_Tq_new_2+h_true)]

                       ! Check for negatives values.
                       i_neg = 0
                       do i=1,n_Tq_2 
                          if (x_n(i) .LT. 0.d0) then
                             x_n(i) = x_ap(i)
                             i_neg = 1
                             if(x_n(i).LT.-1e-4) then
                                ab = 1
                                i_neg = 2
                                print*, 'skip_prof x_n < -1e-4.'
                                SKIP_PROF = .TRUE.
                                GOTO 666   ! GOTO SKIP_PROF
                             end if
                          end if
                       end do

                       ! Check for super-saturation above water.
                       i_ssat = 0
                       q_check = x_n(n_Tq+1:n_Tq_2)
                       T_check = x_n(1:n_Tq)
                       ! -- water vapor pressure
                       e_check = q_check*Rw*T_check
                       ! -- water vapor saturation pressure
                       es_check = e0*EXP(faktor*(T_check-T0)/T_check)
                       ! -- rel. humidity
                       rh_check = (e_check/es_check)*100.
                       if (COUNT(rh_check .GT. 110.) .GT. 0) i_ssat = 1

                       ! Calculate y_n = F(xd)
                       select case (abs)
                       case ('r98')
                          call stp_fap('R98',z_lr_lwc,x_n,h_true,TB_n,z_final,p_rs,coeff_fap,&
                               bias_coeff_fap)
                       case ('l93')
                          call stp_fap('L93',z_lr_lwc,x_n,h_true,TB_n,z_final,p_rs,coeff_fap,&
                               bias_coeff_fap)
                       case default
                          STOP 'Error selecting ABS for y_n=F(xd).'
                       end select
                       allocate(garbage2(size(y_n)), dy(size(y_n)))
                       garbage2 = y_n           ! y_n_old = garbage2
                       select case (use_dbz)
                       case ('NO')
                          y_n = [TB_n, x_n(1), x_n(n_Tq+1), x_n(n_Tq+ind)]
                       case ('YES')
                          y_n = [TB_n,c+d*x_n(n_Tq_2+1:n_Tq_2+h_true), x_n(1), x_n(n_Tq+1), x_n(n_Tq+ind)]
                       case default
                          STOP 'Error selecting USE_DBZ for y_n=F(xd)'
                       end select

                       ! +++++++++++++++++++++++++++++++++++++++
                       ! 7.) CHECK FOR CONVERGENCES CRITERION
                       ! +++++++++++++++++++++++++++++++++++++++
                       dy = y_n - garbage2
                       ss = size(Se_fap,DIM=1)
                       allocate(Sdy(ss,ss), Sdy_inv(ss,ss), stat=allost)
                       if(allost/=0) stop 'Allocate Sdy, Sdy_inv failed...!!'

                       Sdy = MATMUL(Se_fap, MATMUL(Sd_m_inv, Se_fap))
                       call invert_LU(Sdy, Sdy_inv, ss, error_status)
                       deallocate(Sdy, garbage2)
                       if (error_status .EQ. 1) then !Invert matrix Sap
                          ab = 1 ; mat_inv = 4
                          deallocate(x_n_s, x_ap_s, xd, Sd_m, Sd_m_inv, dy, Sdy_inv)
                          print*, 'Out from inversion of Sdy matrix.'
                          SKIP_PROF = .TRUE.
                          GOTO 666    ! SKIP_PROF
                       end if

                       d_i2 = minval(MATMUL(spread(dy,1,1), MATMUL(Sdy_inv,dy)))
                       kf = 0
                       if (d_i2_old.LT.d_i2) kf = 1
                       if (d_i2.LT.0.) then
                          kf = 2 ; ab = 1
                       end if
                       write(*,'(A,F13.5,F6.2)') 'd_i2:',d_i2,(n_freq+2*h_true+2)/10.
                       if (ab.EQ.1) then
                          deallocate(x_n_s, x_ap_s, xd, Sd_m, Sd_m_inv, dy, Sdy_inv)
                          SKIP_PROF = .TRUE.
                          GOTO 666    !SKIP_PROF
                       end if
                       k = 0
                       if (d_i2.LT.(n_freq+2.*h_true+2.)/10.) then
                          k = 1
                          deallocate(Sd_m, Sd_m_inv, dy, Sdy_inv)
                          EXIT ITER_IT
                       end if
                       d_i2_old = d_i2

                       deallocate(x_n_s, x_ap_s, xd, Sd_m, Sd_m_inv, dy, Sdy_inv)
                    end do ITER_IT

                    if (k.EQ.0) then
                       if (fudge_sap.EQ.1.) then
                          fudge_sap = 0.1
                          fudge_cor = 'FUDGE'
                          REP_PROF = .TRUE.
                          print*, 'Repeat Profile: ' 
                          GOTO 666 ! REPEAT_PROF
                       end if
                       ab = 1
                       fudge_sap = 1.
                       SKIP_PROF = .TRUE.
                       GOTO 666         ! GOTO, SKIP_PROF
                    end if
                    fudge_sap = 1.
                    ! +++++++++++++++++++++++++
                    ! 8.) SOLUTION CHECKS
                    ! +++++++++++++++++++++++++
                    ! ++  test  y_gen/y_ret
                    !     spurious convergence ?
                    !     (Rodgers 12.9)

                    ss = size(Se_fap, DIM=1)
                    allocate(dF_s(size(y_n)), Sd_s(ss,ss), Sd_s_inv(ss,ss),stat=allost)
                    if (allost/=0) stop 'Allocate dF_s, Sd_s, Sd_s_inv failed.'
                    dF_s = y_n - y_gen
                    Sd_s = MATMUL(tK_J, MATMUL(Sap, K_J)) + Se_fap

                    ! +++ Invert Sd_s
                    call invert_LU(Sd_s, Sd_s_inv, ss, error_status)
                    deallocate(Sd_s)
                    MAT_SD_S: if (error_status .EQ. 1) then
                       ab = 1
                       mat_inv = 5
                       print *, 'Out from inversion of Sd_s matrix.'
                       SKIP_PROF = .TRUE.
                    else
                       ss = size(Se_fap, DIM=1)
                       allocate(S_deyd(ss,ss), S_deyd_inv(ss,ss), stat=allost)
                       if (allost/=0) stop 'Allocate S_deyd, S_deyd_inv failed.'
                       S_deyd = MATMUL(Se_fap,MATMUL(Sd_s_inv,Se_fap))
                       ! +++ Invert S_deyd
                       call invert_LU(S_deyd, S_deyd_inv, ss, error_status)
                       deallocate(S_deyd)
                       MAT_S_DEYD: if (error_status .EQ. 1) then
                          ab = 1
                          mat_inv = 6
                          print *, 'Out from inversion of S_deyd matrix.'
                          SKIP_PROF = .TRUE.
                       else
                          chi_1 = minval(MATMUL(spread(dF_s,1,1),MATMUL(S_deyd_inv,dF_s)))

                          ! +++++ Test y_gen/y_ap
                          ! (Rodgers 12.3.3.1)
                          dF_s = y_gen - y_ap
                          chi_2 = minval(MATMUL(spread(dF_s,1,1),MATMUL(Sd_s_inv,dF_s)))

                          ! +++++ Test y_ret/y_ap
                          ! (Rodgers 12.16)
                          ss = size(Se_fap,DIM=1)
                          allocate(dy(size(y_n)), garbage3(ss,ss), stat=allost)
                          if (allost/=0) stop 'Allocate dy, garbage3 failed.'
                          dy = y_n - y_ap
                          garbage3 = MATMUL(tK_J, MATMUL(Sap, MATMUL(K_J, MATMUL(Sd_s_inv, &
                               MATMUL(tK_J, MATMUL(Sap, K_J))))))
                          chi_4 = minval(MATMUL(spread(dy,1,1), MATMUL(garbage3,dy)))
                          deallocate(dy, garbage3)

                          ! +++++ Test x_ret/x_ap
                          ! (Rodgers 12.12)
                          ss = size(Sap,DIM=1)
                          allocate(dx(size(xd)), garbage3(ss,ss), stat=allost)
                          if (allost/=0) stop 'Allocate dx, garbage3 failed.'
                          dx = x_ap_s - xd
                          garbage3 = MATMUL(Sap,MATMUL(K_J,MATMUL(Sd_s_inv,MATMUL(tK_J,Sap))))
                          chi_3 = minval(MATMUL(spread(dx,1,1),MATMUL(garbage3,dx)))
                          deallocate(dx, garbage3)

                          ! +++++  chi test value: y

                          chi_test_y = chisqr_cvf(0.05, size(y_gen))

                          ! +++++ chi test value: x

                          chi_test_x = chisqr_cvf(0.05, size(xd))
                          deallocate(xd, x_ap_s)

                          ! +++ Calculate theoretical optimal estimation errors
                          allocate(K3_e(n_Tq_2 + h_true, n_freq),stat=allost)
                          if(allost/=0) stop 'Allocate K3_e failed...!!'

                          select case (abs)
                          case ('r98')
                             call make_K3_fap('R98', z_lr_lwc, x_n, z_final, p_rs, h_true, inc, &
                                  coeff_fap, bias_coeff_fap, cloud, K3_e)
                          case ('l93')
                             call make_K3_fap('L93', z_lr_lwc, x_n, z_final, p_rs, h_true, inc, &
                                  coeff_fap, bias_coeff_fap, cloud, K3_e)
                          case default
                             stop 'Error selecting make_K3_fap subroutine for K3_e.'
                          end select

                          select case (use_dbz)
                          case ('YES')
                             ! ++++ Make K1_e subsection
                             ! 
                             ! Here, temporal variable "garbage3" is "K1_e".
                             ! Instant of MAKE_K1 subroutine.

                             allocate(garbage3(h_true, h_true)) ; garbage3 = 0.d0
                             forall(i=1:h_true,j=1:h_true,i.EQ.j)
                                !  1. Ableitung der Z/LWC Beziehung
                                !   K1(i,j) = d(i)*10.*1./(lwc(i)*log(10.))
                                garbage3(i,j) = d(i)              ! K1_e(i,j) = d(i)
                             end forall
                             allocate(K_e(n_Tq_new_2+h_true, n_freq+h_true+2+h_true), &
                                  tK_e(n_freq+h_true+2+h_true, n_Tq_new_2+h_true), stat=allost)
                             if (allost/=0) stop 'Allocate of K_e with ABS=YES failed.'
                             K_e = 0.d0 ; tK_e = 0.d0
                             K_e(1:n_Tq_new, 1:n_freq) = K3_e(1:n_Tq_new, 1:n_freq)
                             K_e(n_Tq_new+1:n_Tq_new_2, 1:n_freq) = K3_e(n_Tq+1:n_Tq+n_Tq_new, 1:n_freq)
                             K_e(n_Tq_new_2+1:n_Tq_new_2+h_true, 1:n_freq) = K3_e(n_Tq_2+1:n_Tq_2+h_true, 1:n_freq)
                             K_e(n_Tq_new_2+1:n_Tq_new_2+h_true, n_freq+1:n_freq+h_true) = garbage3
                             K_e(1,n_freq+h_true+1) = 1.d0       ! T_gr
                             K_e(n_Tq_new+1,n_freq+h_true+2) = 1.d0 !q_gr
                             do j=1,h_true
                                K_e(n_Tq_new+ind(j), n_freq+h_true+2+j) = 1.d0
                             end do
                             deallocate(garbage3)
                          case ('NO')
                             allocate(K_e(n_Tq_new_2+h_true, n_freq+2+h_true), &
                                  tK_e(n_freq+2+h_true, n_Tq_new_2+h_true), stat=allost)
                             if (allost/=0) stop 'Allocate of K_e with ABS=NO failed.'
                             K_e = 0.d0 ; tK_e = 0.d0
                             K_e(1:n_Tq_new, 1:n_freq) = K3_e(1:n_Tq_new, 1:n_freq)
                             K_e(n_Tq_new+1:n_Tq_new_2, 1:n_freq) = K3_e(n_Tq+1:n_Tq+n_Tq_new, 1:n_freq)
                             K_e(n_Tq_new_2+1:n_Tq_new_2+h_true, 1:n_freq) = K3_e(n_Tq_2+1:n_Tq_2+h_true, 1:n_freq)
                             K_e(1, n_freq+1) = 1.d0             ! T_gr
                             K_e(n_Tq_new+1, n_freq+2) = 1.d0    ! q_gr
                             do j=1,h_true
                                K_e(n_Tq_new+ind(j), n_freq+2+j) = 1.d0
                             end do
                          case default
                             STOP 'Error selecting use_dbz: YES / NO for K_e. '
                          end select

                          deallocate(K3_e)

                          tK_e = transpose(K_e)

                          ss = size(Se_fap, DIM=1)
                          allocate(Sx(ss,ss), Sx_inv(ss,ss), stat=allost)
                          if (allost/=0) stop 'Allocate Sx, Sx_inv matrices failed.'
                          Sx = Se_fap + MATMUL(tK_e,MATMUL(Sap,K_e))
                          ! ++++ Invert Sx matrix
                          call invert_LU(Sx, Sx_inv, ss, error_status)
                          deallocate(Sx)
                          if (error_status .EQ. 1) then
                             ab = 1
                             mat_inv = 7
                             SKIP_PROF = .TRUE.
                          else
                             ss = size(Sap, DIM=1)
                             allocate(Sd(ss,ss), err_oe(n_Tq_new_2+h_true), stat=allost)
                             if (allost/=0) stop 'Allocate Sd, err_oe failed.'

                             Sd = Sap - MATMUL(Sap,MATMUL(K_e,MATMUL(Sx_inv,MATMUL(tK_e,Sap))))
                             err_oe = 0.d0
                             do j = 1,n_Tq_new_2+h_true        ! j = i_err
                                err_oe(j) = SQRT(Sd(j,j))
                             end do
                             ! +++++++++++++++++++++++++++++++++++++
                             ! 9.) WRITE RETRIEVAL RESULTS TO FILE
                             ! +++++++++++++++++++++++++++++++++++++

                             ! ++++ reset convergence criterion
                             d_i2_old = 1e10

                          end if
                       end if MAT_S_DEYD       ! End matrix inversion S_deyd
                    end if MAT_SD_S          ! End matrix inversion Sd_s

666                 CONTINUE !SKIP_PROF or REPEAT_PROF:
                    !end if       ! End else first guess.
                    !end if           ! end if Invert matrix Sap
                    !end if           ! end if Invert matrix Se_fap
                 end if V_CLOUD              ! end check vertical cloud thickness
              end if BASE_TOP
              RESU: if (.NOT. REP_PROF) then
                 ! +++++++++++++++++++++++++++++++++++++
                 ! 9.) WRITE RETRIEVAL RESULTS TO FILE
                 ! +++++++++++++++++++++++++++++++++++++
                 if(ab.EQ.0) status='OK'
                 if(ab.EQ.1) status='ABORT'
                 write(13,'(A8,F10.6)') date(i_days),time_a(k_pf)
                 write(13,*) mode
                 write(13,*) status
                 write(13,*) fudge_cor
                 write(13,*) TRANSPOSE(code_a(k_pf:k_pf,:))
                 write(13,'(20G15.7)') T_rs_a(k_pf,:)
                 write(13,'(20G15.7)') q_rs_a(k_pf,:)
                 write(13,'(20G15.7)') p_rs_a(k_pf,:)
                 write(13,'(20G15.7)') z_rad_liquid
                 write(13,'(I6)') nd_rad
                 write(13,'(20G15.7)') z_hr
                 if (ab.EQ.0) write(13,'(20G15.7)') dbz_hr
		 if (ab.EQ.1) write(13,'(20G15.7)') dbz_hr_a(k_pf,:)
                 write(13,'(20G15.7)') vel_hr_a(k_pf,:)
                 write(13,'(20G15.7)') class_hr_a(k_pf,:)
                 write(13,'(3G15.7)') lwp_mic_a(k_pf), iwv_mic_a(k_pf), base_hr_a(k_pf)

                 write(13,'(A14,I4)') 'cloud layers:',nn
                 write(13,'(A14,I4)') 'cloud levels:',h_true
                 write(13,'(A6,5G15.7)') 'base:',base_m
                 write(13,'(A6,5G15.7)') 'top:',top_m
                 write(13,'(A10,I3)') 'TB check:',TB_faulty
                 write(13,'(A18,I3)') 'cloud properties:',cloud_prop
                 write(13,'(A18,I4)') 'matrix inversion:',mat_inv
                 write(13,'(A7,I3)') 'neg.q:',i_neg
                 write(13,'(A14,I3)') 'cost function:',kf
                 write(13,'(A10,I3)') 'num.its:',it
                 write(13,'(A10,I3)') 'ssat check:',i_ssat
                 write(13,'(A17,G15.3)') 'y-chi-test(95%):',chi_test_y
                 write(13,'(A13,G15.3)') 'y_ret/y_gen:',chi_1
                 write(13,'(A13,G15.3)') 'y_ap/y_gen:',chi_2
                 write(13,'(A17,G15.3)') 'x-chi-test(95%):',chi_test_x
                 write(13,'(A13,G15.3)') 'x_ret/x_ap:',chi_3

                 if (ab.EQ.0) then
                    write(13,*) 'results:'
                    write(13,'(20G15.7)') y_gen
                    write(13,'(20G15.7)') x_n
                    write(13,'(20G15.7)') y_n
                    write(13,'(20G15.7)') err_oe
                    write(13,'(20G15.7)') rh_check
                    write(13,'(20G15.7)') z_lr_lwc
                    write(13,'(20G15.7)') (10.**(x_n(n_Tq_2+1:n_Tq_2+h_true)/10.))
                    call lwc_high_res(z_lr_lwc, h_true_m, &
                         x_n_s(n_Tq_new_2+1:n_Tq_new_2+h_true), x_n(n_Tq_2+1:n_Tq_2+h_true), c, d, dbz_hr, z_hr, nn )

                    call cpu_time(ptime_e)
                    print '(A,F6.1,A)', 'Finish succesfully after',(ptime_e-ptime_s),' seconds.'
                 end if
              end if RESU
              ! DEALLOCATING ALL VARIABLES IN PROFILE LOOP
              if (allocated(x_n_s)) deallocate(x_n_s)
              if (allocated(Sd)) deallocate(Sd)
              if (allocated(err_oe)) deallocate(err_oe)
              if (allocated(Sx_inv) .AND. allocated(K_e) .AND. allocated(tK_e)) deallocate(Sx_inv, K_e, tK_e)
              if (allocated(S_deyd_inv)) deallocate(S_deyd_inv)
              if (allocated(dF_s) .AND. allocated(Sd_s_inv)) deallocate(dF_s, Sd_s_inv)
              if (allocated(x_fg) .AND. allocated(x_n)) deallocate(x_fg, x_n)
              if (allocated(y_n)) deallocate(y_n)
              if (allocated(K_J) .AND. allocated(tK_J)) deallocate(K_J, tK_J)
              if (allocated(x_ap) .AND. allocated(y_ap)) deallocate(x_ap, y_ap)
              if (allocated(dbz_hr) .AND. allocated(y_gen)) deallocate(dbz_hr, y_gen)
              if (allocated(c) .AND. allocated(d)) deallocate(c,d)
              if (allocated(lwc_ap) .AND. allocated(Sap)) deallocate(lwc_ap,Sap)
              if (allocated(Se) .AND. allocated(Se_fap) .AND. allocated(inc)) deallocate(Se,Se_fap,inc)
              if (allocated(qc) .AND. allocated(ind) .AND. allocated(e_qc)) deallocate(qc, ind, e_qc)
              if (allocated(base_m)) deallocate(base_m)
              if (allocated(top_m)) deallocate(top_m)
              if (allocated(h_true_m) .AND. allocated(h_m)) deallocate(h_true_m, h_m)
              if (allocated(dbz) .AND. allocated(z_lr_lwc)) deallocate(dbz, z_lr_lwc)
           endif APP_CLOUDY                      !Applicable cloudy..
           if (REP_PROF) CYCLE REPEAT_PROF

! *****************************************************************
! -----------------------------------------------------------------
!
!           CLOUD FREE PROCEDURE
! -----------------------------------------------------------------
! *****************************************************************

           status = ''
           SKIP_PROF_ICY = .FALSE.
           ! Free applicable cloudy ?
           APP_FREE_CLOUDY: if ((code_c_a(k_pf) .EQ. 1) .AND. (time_a(k_pf) .GT. time_crit)) then
              cloud = 0              ! STP parameter
              ! abort parameters
              kf = -99 ; i_qn = -99 ; chi_1 = -99. ; chi_2 = -99. ; chi_3 = -99.
              chi_4 = -99. ; chi_test_x = -99. ; chi_test_y = -99. ; chi_check = -99
              st_1 = -99 ; st_2 = -99 ; it = -99 ; TB_faulty = -99
              ! Output info
              print *,'DATE/TIME: ', date(i_days),' ',time_a(k_pf)
              print *,'Cloud-free/icy'
              mode = 'CF/ICY'
              call cpu_time(ptime_s)
              !     ++++ specify T, q profiles to use
              jh = 4
              !     ++++ Define ground measurements
              T_gr = tem_a(k_pf)
              tau = dpt_a(k_pf)
              faktor = L/(T0*Rw)
              q_gr = (e0*EXP(faktor*(tau - T0)/tau))/(T_gr*Rw)

              ! ++++  Define Sap
              allocate(Sap(n_Tq_new_2, n_Tq_new_2), stat=allost)
              if(allost/=0) stop 'Allocate Sap (icy) matrix failed...!!'
              Sap = 0.d0
              ! Define Sap_eof due to RS "a priori"
              Sap(1:n_Tq_new, 1:n_Tq_new) = cov_rs(1:n_Tq_new, 1:n_Tq_new)
              Sap(1:n_Tq_new, n_Tq_new+1:n_Tq_new_2) = cov_rs(1:n_Tq_new, n_Tq+1:n_Tq+n_Tq_new)
              Sap(n_Tq_new+1:n_Tq_new_2,1:n_Tq_new) = cov_rs(n_Tq+1:n_Tq+n_Tq_new,1:n_Tq_new)
              Sap(n_Tq_new+1:n_Tq_new_2,n_Tq_new+1:n_Tq_new_2) = cov_rs(n_Tq+1:n_Tq+n_Tq_new,n_Tq+1:n_Tq+n_Tq_new)
              ! Consider cloud covariances
              Sap = fudge_sap*Sap
              Sap(1:n_Tq_new,1:n_Tq_new) = fudge_sap_T*Sap(1:n_Tq_new,1:n_Tq_new)
              Sap(n_Tq_new+1:n_Tq_new_2,n_Tq_new+1:n_Tq_new_2) = &
                   fudge_sap_q*Sap(n_Tq_new+1:n_Tq_new_2, n_Tq_new+1:n_Tq_new_2)

              ! Define increments for Jacobian
              allocate(inc(n_Tq_2), stat=allost)
              if(allost/=0) stop 'Allocate inc (icy) failed...!!'
              j = 1 ; k = n_Tq_new_2 + 1 
              do i =1,n_Tq_2
                 if(cov_rs(i,i).NE.0) j=i
                 inc(i) = SQRT(cov_rs(j,j))/10.
              end do

              ! ++++ Determine Se
              allocate(garbage3(n_freq_r+2, n_freq_r+2))
              garbage3 = 0.d0
              garbage3(1:n_freq_r,1:n_freq_r) = Se_a(jh-2,1:n_freq_r,1:n_freq_r)
              garbage3(n_freq_r+1,n_freq_r+1) = e_T_gr**2
              garbage3(n_freq_r+2,n_freq_r+2) = e_q_gr**2
              
              allocate(Se(n_freq+2,n_freq+2))
              Se = 0.d0
              forall(i=1:n_freq,j=1:n_freq,.TRUE.)
                 Se(i,j) = garbage3(freq_num(i),freq_num(j))
              end forall
              Se(n_freq+1:n_freq+2, n_freq+1:n_freq+2) = &
                   garbage3(n_freq_r+1:n_freq_r+2, n_freq_r+1:n_freq_r+2)
              deallocate(garbage3)
              allocate(Se_fap(n_freq+2,n_freq+2), stat=allost)
              if(allost/=0) stop 'Allocate Se_fap (icy) matrix failed...!!'
              Se_fap = 0.d0
              
              Se_fap = Se
              Se_fap(1:n_freq,1:n_freq)=Se(1:n_freq,1:n_freq) + Se_fap_org
              Se_fap = fudge_se * Se_fap
              Se_fap(1:n_freq,1:n_freq)=fudge_se_tb * Se_fap(1:n_freq,1:n_freq)

              !  ++++ Invert matrix Se_fap
              ss = size(Se_fap,DIM=1)
              allocate(garbage3(ss,ss))
              call invert_LU(Se_fap, garbage3, ss, error_status)
              deallocate(garbage3)
              if (error_status .EQ. 1) then   ! Invert matrix Se_fap
                 ab = 1 ; mat_inv = 1
                 print *, 'Out from inversion of Se_fap (icy) matrix.'
                 SKIP_PROF_ICY = .TRUE.
                 GOTO 777
              end if

              !  ++++ Invert matrix Sap
              ss = size(Sap,DIM=1)
              allocate(garbage3(ss,ss))
              call invert_LU(Sap, garbage3, ss, error_status)
              deallocate(garbage3)
              if (error_status .EQ. 1) then   !! .EQ. 1 ! Invert matrix Sap
                 ab = 1 ; mat_inv = 2
                 print *, 'Out from inversion of Sap (icy) matrix.'
                 SKIP_PROF_ICY = .TRUE.
                 GOTO 777
              end if

              ! ++++ Initialization
              allocate(y_gen(n_freq + 2))
              y_gen = [TB_a(k_pf, freq_num), T_gr, q_gr]

              TB_faulty = 0
              if (count(TB_a(k_pf, freq_num).LT.0) .GT. 0 ) then
                 ab = 1 ; TB_faulty = 1
                 print *, 'NOT first guess.'
                 SKIP_PROF_ICY = .TRUE.
                 GOTO 777
              end if
              ! ++++ first guess
              allocate(x_fg(n_Tq_2))
              x_fg = [real(T_rs,KIND=8), real(q_rs,KIND=8)]   ! First guess profile.

              ! ++++ calculate y_n = F(x_n)

              select case (abs)
              case ('r98')
                 call stp_fap('R98',[0.],x_fg,0,TB_n,z_final,p_rs,coeff_fap,bias_coeff_fap)
              case ('l93')
                 call stp_fap('L93',[0.],x_fg,0,TB_n,z_final,p_rs,coeff_fap,bias_coeff_fap)
              case default
                 stop 'Error selecting ABS for y_n = F(x_n) (icy).!!!'
              end select
              
              ! ++++ Asigning y_n and Dimension to K_J and tK_J matrices.
              allocate(y_n(n_freq+2), K_J(n_Tq_new_2, n_freq+2), tK_J(n_freq+2, n_Tq_new_2))
              y_n = [TB_n, x_fg(1), x_fg(n_Tq+1)]
              K_J = 0.d0

              allocate(x_n(size(x_fg)), stat=allost)
              if(allost/=0) stop 'Allocate x_n (icy) failed...!!'
              x_n = x_fg

              ! +++ connocate a priori.
              allocate(x_ap(n_Tq_2), y_ap(size(y_n)), stat=allost)
              if(allost/=0) stop 'Allocate x_ap, y_ap (icy) failed...!!'
              x_ap = [real(T_rs,KIND=8), real(q_rs,KIND=8)]
              y_ap = y_n

              ! +++++++++++++++++++++
              ! 4.) BEGIN ITERATION
              ! +++++++++++++++++++++
              ITER_IT_ICY: do it = 0, 9
                 ! +++ Cost function check.
                 kf = -99
                 ! +++ q negative check.
                 i_neg = -99
                 ! inversion checks.
                 st_1 = -99   ;  st_2 = -99

                 ! +++++++++++++++++++++++++++++++++++++++++
                 ! 5.) CALCULATE K_J BY COMBINING K1 AND K3
                 ! +++++++++++++++++++++++++++++++++++++++++
                 ! (K3: de/dTB, df/dTB, dlwc/dTB)
                 ! (K1: dLWC/ddbz)
                 allocate(K3(n_Tq_2, n_freq),stat=allost)
                 if(allost/=0) stop 'Allocate K3 (icy) failed...!!'
                 K3 = 0.d0
                 select case (abs)
                 case ('r98')
                    call make_K3_fap('R98', [0.], x_n, z_final, p_rs, 0, inc, &
                         coeff_fap, bias_coeff_fap, cloud, K3)
                 case ('l93')
                    call make_K3_fap('L93', [0.], x_n, z_final, p_rs, 0, inc, &
                         coeff_fap, bias_coeff_fap, cloud, K3)
                 case default
                    stop 'Error selecting make_K3_fap (icy) subroutine.'
                 end select

                 K_J(1:n_Tq_new, 1:n_freq) = K3(1:n_Tq_new, 1:n_freq)
                 K_J(n_Tq_new+1:n_Tq_new_2, 1:n_freq) = K3(n_Tq+1:n_Tq+n_Tq_new, 1:n_freq)
                 K_J(n_Tq_new_2+1:n_Tq_new_2+h_true, 1:n_freq) = K3(n_Tq_2+1:n_Tq_2+h_true, 1:n_freq)
                 K_J(1, n_freq+1) = 1.d0           ! T_gr
                 K_J(n_Tq_new+1, n_freq+2) = 1.d0  ! q_gr
                 deallocate(K3)
                 
                 ! +++++++++++++++++++++++++++
                 ! 6.) OPTIMAL ESTIMATION
                 ! +++++++++++++++++++++++++++
                 ! Rodgers (5.8 / 5.9)
                 allocate(x_n_s(n_Tq_new_2), x_ap_s(n_Tq_new_2), &
                      xd(n_Tq_new_2), stat=allost)
                 if(allost/=0) stop 'Allocate x_n_s, x_ap_s, xd (icy) failed...!!'
                 x_n_s = [x_n(1:n_Tq_new),x_n(n_Tq+1:n_Tq+n_Tq_new)]
                 x_ap_s = [x_ap(1:n_Tq_new),x_ap(n_Tq+1:n_Tq+n_Tq_new)]
                 tK_J = transpose(K_J)
                 ss = size(Se_fap,DIM=1)
                 allocate(Sd_m(ss,ss), Sd_m_inv(ss,ss), stat=allost)
                 if (allost/=0) stop 'Allocate Sd_m (icy) failed...!!!'

                 Sd_m = Se_fap + matmul(tK_J, matmul(Sap,K_J))

                 call invert_LU(Sd_m, Sd_m_inv, ss, error_status)
                 if (error_status .EQ. 1) then   !! .EQ. 1 ! Invert matrix Sap
                    ab = 1 ; mat_inv = 3
                    deallocate(x_n_s, x_ap_s, xd, Sd_m, Sd_m_inv)
                    print *, 'Out from inversion of Sd_m (icy) matrix.'
                    SKIP_PROF_ICY = .TRUE.
                    GOTO 777    ! SKIP_PROF
                 end if

                 allocate(bb(size(y_gen)), cc(size(tK_J,DIM=1), size(Sap,DIM=2)), garbage2(size(x_n_s)),stat=allost)
                 if(allost/=0) stop 'Allocate bb, cc (icy) failed...!!'
                 garbage2 = x_n_s - x_ap_s
                 bb = y_gen - y_n + MATMUL(garbage2, K_J)
                 cc = MATMUL(tK_J, Sap)

                 xd = x_ap_s + MATMUL(bb, MATMUL(Sd_m_inv, cc))
                 deallocate(bb,cc, garbage2)
                 x_n = [xd(1:n_Tq_new), x_fg(n_Tq_new+1:n_Tq), xd(n_Tq_new+1:n_Tq_new_2), &
                      x_fg(n_Tq+n_Tq_new+1:n_Tq_2)]

                 ! Check for negatives values.
                 i_neg = 0
                 do i=1,n_Tq_2 
                    if (x_n(i) .LT. 0.d0) then
                       x_n(i) = x_ap(i)
                       i_neg = 1
                       if(x_n(i).LT.-1e-4) then
                          ab = 1
                          i_neg = 2
                          print*, 'skip_prof x_n < -1e-4 (icy).'
                          SKIP_PROF_ICY = .TRUE.
                          GOTO 777   ! GOTO SKIP_PROF
                       end if
                    end if
                 end do
                 ! Check for super-saturation above water.
                 i_ssat = 0
                 q_check = x_n(n_Tq+1:n_Tq_2)
                 T_check = x_n(1:n_Tq)
                 ! -- water vapor pressure
                 e_check = q_check*Rw*T_check
                 ! -- water vapor saturation pressure
                 es_check = e0*EXP(faktor*(T_check-T0)/T_check)
                 ! -- rel. humidity
                 rh_check = (e_check/es_check)*100.
                 if (COUNT(rh_check .GT. 110.) .GT. 0) then
                    if (MAXVAL(z_final,MASK=rh_check.GT.110.) .LT. 5000.) i_ssat = 1
                 end if
                 ! Calculate y_n = F(xd)
                 select case (abs)
                 case ('r98')
                    call stp_fap('R98',[0.],x_n,0,TB_n,z_final,p_rs,coeff_fap,&
                         bias_coeff_fap)
                 case ('l93')
                    call stp_fap('L93',[0.],x_n,0,TB_n,z_final,p_rs,coeff_fap,&
                         bias_coeff_fap)
                 case default
                    STOP 'Error selecting ABS (icy) for y_n=F(xd).'
                 end select
                 allocate(garbage2(size(y_n)), dy(size(y_n)))
                 garbage2 = y_n           ! y_n_old = garbage2

                 y_n = [TB_n, x_n(1), x_n(n_Tq+1)]

                 ! +++++++++++++++++++++++++++++++++++++++
                 ! 7.) CHECK FOR CONVERGENCES CRITERION
                 ! +++++++++++++++++++++++++++++++++++++++
                 dy = y_n - garbage2
                 ss = size(Se_fap,DIM=1)
                 allocate(Sdy(ss,ss), Sdy_inv(ss,ss), stat=allost)
                 if(allost/=0) stop 'Allocate Sdy, Sdy_inv (icy) failed...!!'

                 Sdy = MATMUL(Se_fap, MATMUL(Sd_m_inv, Se_fap))
                 call invert_LU(Sdy, Sdy_inv, ss, error_status)
                 deallocate(Sdy, garbage2)
                 if (error_status .EQ. 1) then !Invert matrix Sap
                    ab = 1 ; mat_inv = 4
                    deallocate(x_n_s, x_ap_s, xd, Sd_m, Sd_m_inv, dy, Sdy_inv)
                    print*, 'Out from inversion of Sdy (icy) matrix.'
                    SKIP_PROF_ICY = .TRUE.
                    GOTO 777    ! SKIP_PROF
                 end if

                 d_i2 = minval(MATMUL(spread(dy,1,1), MATMUL(Sdy_inv,dy)))

                 kf = 0
                 if (d_i2_old.LT.d_i2) kf = 1
                 if (d_i2.LT.0.) then
                    kf = 2 ; ab = 1
                 end if
                 write(*,'(A,F13.5,F6.2)') 'd_i2:',d_i2,(n_freq+2)/10.
                 if (ab.EQ.1) then
                    deallocate(x_n_s, x_ap_s, xd, Sd_m, Sd_m_inv, dy, Sdy_inv)
                    SKIP_PROF_ICY = .TRUE.
                    GOTO 777    !SKIP_PROF
                 end if
                 k = 0
                 if (d_i2.LT.(n_freq+2.)/10.) then
                    k = 1
                    deallocate(Sd_m, Sd_m_inv, dy, Sdy_inv)
                    EXIT ITER_IT_ICY
                 end if
                 d_i2_old = d_i2

                 deallocate(x_n_s, x_ap_s, xd, Sd_m, Sd_m_inv, dy, Sdy_inv)
              end do ITER_IT_ICY

              if (k.EQ.0) then
                 if (fudge_sap.EQ.1.) then
                    fudge_sap = 0.1
                    fudge_cor = 'FUDGE'
                    REP_PROF = .TRUE.
                    print*, 'Repeat Profile: ' 
                    GOTO 777 ! REPEAT_PROF
                 end if
                 ab = 1
                 fudge_sap = 1.
                 SKIP_PROF_ICY = .TRUE.
                 GOTO 777         ! GOTO, SKIP_PROF
              end if
              fudge_sap = 1.

              ! +++++++++++++++++++++++++
              ! 8.) SOLUTION CHECKS
              ! +++++++++++++++++++++++++
              ! ++  test  y_gen/y_ret
              !     spurious convergence ?
              !     (Rodgers 12.9)

              ss = size(Se_fap, DIM=1)
              allocate(dF_s(size(y_n)), Sd_s(ss,ss), Sd_s_inv(ss,ss),stat=allost)
              if (allost/=0) stop 'Allocate dF_s, Sd_s, Sd_s_inv (icy) failed.'
              dF_s = y_n - y_gen
              Sd_s = MATMUL(tK_J, MATMUL(Sap, K_J)) + Se_fap

              ! +++ Invert Sd_s
              call invert_LU(Sd_s, Sd_s_inv, ss, error_status)
              deallocate(Sd_s)
              MAT_SD_S2: if (error_status .EQ. 1) then
                 ab = 1
                 mat_inv = 5
                 print *, 'Out from inversion of Sd_s (icy) matrix.'
                 SKIP_PROF_ICY = .TRUE.
              else
                 ss = size(Se_fap, DIM=1)
                 allocate(S_deyd(ss,ss), S_deyd_inv(ss,ss), stat=allost)
                 if (allost/=0) stop 'Allocate S_deyd, S_deyd_inv (icy) failed.'
                 S_deyd = MATMUL(Se_fap,MATMUL(Sd_s_inv,Se_fap))
                 ! +++ Invert S_deyd
                 call invert_LU(S_deyd, S_deyd_inv, ss, error_status)
                 deallocate(S_deyd)
                 MAT_S_DEYD2: if (error_status .EQ. 1) then
                    ab = 1
                    mat_inv = 6
                    print *, 'Out from inversion of S_deyd (icy) matrix.'
                    SKIP_PROF_ICY = .TRUE.
                 else
                    chi_1 = minval(MATMUL(spread(dF_s,1,1),MATMUL(S_deyd_inv,dF_s)))
                    
                    ! +++++ Test y_gen/y_ap
                    ! (Rodgers 12.3.3.1)
                    dF_s = y_gen - y_ap
                    chi_2 = minval(MATMUL(spread(dF_s,1,1),MATMUL(Sd_s_inv,dF_s)))

                    ! +++++ Test y_ret/y_ap
                    ! (Rodgers 12.16)
                    ss = size(Se_fap,DIM=1)
                    allocate(dy(size(y_n)), garbage3(ss,ss), stat=allost)
                    if (allost/=0) stop 'Allocate dy, garbage3 (icy) failed.'
                    dy = y_n - y_ap
                    garbage3 = MATMUL(tK_J, MATMUL(Sap, MATMUL(K_J, MATMUL(Sd_s_inv, &
                         MATMUL(tK_J, MATMUL(Sap, K_J))))))
                    chi_4 = minval(MATMUL(spread(dy,1,1), MATMUL(garbage3,dy)))
                    deallocate(dy, garbage3)

                    ! +++++ Test x_ret/x_ap
                    ! (Rodgers 12.12)
                    ss = size(Sap,DIM=1)
                    allocate(dx(size(xd)), garbage3(ss,ss), stat=allost)
                    if (allost/=0) stop 'Allocate dx, garbage3 (icy) failed.'
                    dx = x_ap_s - xd
                    garbage3 = MATMUL(Sap,MATMUL(K_J,MATMUL(Sd_s_inv,MATMUL(tK_J,Sap))))
                    chi_3 = minval(MATMUL(spread(dx,1,1),MATMUL(garbage3,dx)))
                    deallocate(dx, garbage3)

                    ! +++++  chi test value: y

                    chi_test_y = chisqr_cvf(0.05, size(y_gen))

                    ! +++++ chi test value: x

                    chi_test_x = chisqr_cvf(0.05, size(xd))
                    deallocate(xd, x_ap_s)

                    ! +++ Calculate theoretical optimal estimation errors
                    allocate(K3_e(n_Tq_2, n_freq),stat=allost)
                    if(allost/=0) stop 'Allocate K3_e (icy) failed...!!'

                    select case (abs)
                    case ('r98')
                       call make_K3_fap('R98', [0.], x_n, z_final, p_rs, 0, inc, &
                            coeff_fap, bias_coeff_fap, cloud, K3_e)
                    case ('l93')
                       call make_K3_fap('L93', [0.], x_n, z_final, p_rs, 0, inc, &
                            coeff_fap, bias_coeff_fap, cloud, K3_e)
                    case default
                       stop 'Error selecting make_K3_fap subroutine for K3_e (icy).'
                    end select

                    allocate(K_e(n_Tq_new_2, n_freq+2), &
                         tK_e(n_freq+2, n_Tq_new_2), stat=allost)
                    if (allost/=0) stop 'Allocate of K_e with ABS=NO (icy) failed.'
                    K_e = 0.d0 ; tK_e = 0.d0
                    K_e(1:n_Tq_new, 1:n_freq) = K3_e(1:n_Tq_new, 1:n_freq)
                    K_e(n_Tq_new+1:n_Tq_new_2, 1:n_freq) = K3_e(n_Tq+1:n_Tq+n_Tq_new, 1:n_freq)
                    K_e(1, n_freq+1) = 1.d0             ! T_gr
                    K_e(n_Tq_new+1, n_freq+2) = 1.d0    ! q_gr

                    deallocate(K3_e)

                    tK_e = transpose(K_e)
                    
                    ss = size(Se_fap, DIM=1)
                    allocate(Sx(ss,ss), Sx_inv(ss,ss), stat=allost)
                    if (allost/=0) stop 'Allocate Sx, Sx_inv (icy) matrices failed.'
                    Sx = Se_fap + MATMUL(tK_e,MATMUL(Sap,K_e))
                    ! ++++ Invert Sx matrix
                    call invert_LU(Sx, Sx_inv, ss, error_status)
                    deallocate(Sx)
                    if (error_status .EQ. 1) then
                       ab = 1
                       mat_inv = 7
                       SKIP_PROF_ICY = .TRUE.
                    else
                       ss = size(Sap, DIM=1)
                       allocate(Sd(ss,ss), err_oe(n_Tq_new_2), stat=allost)
                       if (allost/=0) stop 'Allocate Sd, err_oe (icy) failed.'

                       Sd = Sap - MATMUL(Sap,MATMUL(K_e,MATMUL(Sx_inv,MATMUL(tK_e,Sap))))
                       err_oe = 0.d0
                       do j = 1,n_Tq_new_2        ! j = i_err
                          err_oe(j) = SQRT(Sd(j,j))
                       end do
                       ! ++++ reset convergence criterion
                       d_i2_old = 1e10

                    end if
                 end if MAT_S_DEYD2       ! End matrix inversion S_deyd
              end if MAT_SD_S2          ! End matrix inversion Sd_s

777           CONTINUE !SKIP_PROF_ICY or REPEAT_PROF:
              RESU2: if (.NOT. REP_PROF) then
                 ! +++++++++++++++++++++++++++++++++++++
                 ! 9.) WRITE RETRIEVAL RESULTS TO FILE
                 ! +++++++++++++++++++++++++++++++++++++
                 if(ab.EQ.0) status='OK'
                 if(ab.EQ.1) status='ABORT'
                 write(13,'(A8,F10.6)') date(i_days),time_a(k_pf)
                 write(13,*) mode
                 write(13,*) status
                 write(13,*) fudge_cor
                 write(13,*) TRANSPOSE(code_a(k_pf:k_pf,:))
                 write(13,'(20G15.7)') T_rs_a(k_pf,:)
                 write(13,'(20G15.7)') q_rs_a(k_pf,:)
                 write(13,'(20G15.7)') p_rs_a(k_pf,:)
                 write(13,'(20G15.7)') z_rad_liquid
                 write(13,'(I6)') nd_rad
                 write(13,'(20G15.7)') z_hr
                 write(13,'(20G15.7)') dbz_hr_a(k_pf,:)
                 write(13,'(20G15.7)') vel_hr_a(k_pf,:)
                 write(13,'(20G15.7)') class_hr_a(k_pf,:)
                 write(13,'(3G15.7)') lwp_mic_a(k_pf), iwv_mic_a(k_pf), base_hr_a(k_pf)

                 write(13,'(A10,I3)') 'TB check:',TB_faulty
                 write(13,'(A18,I4)') 'matrix inversion:',mat_inv
                 write(13,'(A7,I3)') 'neg.q:',i_neg
                 write(13,'(A14,I3)') 'cost function:',kf
                 write(13,'(A10,I3)') 'num.its:',it
                 write(13,'(A10,I3)') 'ssat check:',i_ssat
                 write(13,'(A17,G15.7)') 'y-chi-test(95%):',chi_test_y
                 write(13,'(A13,G15.7)') 'y_ret/y_gen:',chi_1
                 write(13,'(A13,G15.7)') 'y_ap/y_gen:',chi_2
                 write(13,'(A17,G15.7)') 'x-chi-test(95%):',chi_test_x
                 write(13,'(A13,G15.7)') 'x_ret/x_ap:',chi_3

                 if (ab.EQ.0) then
                    write(13,*) 'results:'
                    write(13,'(20G15.7)') y_gen
                    write(13,'(20G15.7)') x_n
                    write(13,'(20G15.7)') y_n
                    write(13,'(20G15.7)') err_oe
                    write(13,'(20G15.7)') rh_check
                    call cpu_time(ptime_e)
                    print '(A,F6.1,A)', 'Finish succesfully after',(ptime_e-ptime_s),' seconds.'
                 end if
              end if RESU2
              ! DEALLOCATING ALL VARIABLES IN PROFILE LOOP
              if (allocated(x_n_s)) deallocate(x_n_s)
              if (allocated(Sd)) deallocate(Sd)
              if (allocated(err_oe)) deallocate(err_oe)
              if (allocated(Sx_inv) .AND. allocated(K_e) .AND. allocated(tK_e)) deallocate(Sx_inv, K_e, tK_e)
              if (allocated(S_deyd_inv)) deallocate(S_deyd_inv)
              if (allocated(dF_s) .AND. allocated(Sd_s_inv)) deallocate(dF_s, Sd_s_inv)
              if (allocated(x_fg) .AND. allocated(x_n)) deallocate(x_fg, x_n)
              if (allocated(y_n)) deallocate(y_n)
              if (allocated(K_J) .AND. allocated(tK_J)) deallocate(K_J, tK_J)
              if (allocated(x_ap) .AND. allocated(y_ap)) deallocate(x_ap, y_ap)
              if (allocated(y_gen)) deallocate(y_gen)
              if (allocated(Sap)) deallocate(Sap)
              if (allocated(Se) .AND. allocated(Se_fap) .AND. allocated(inc)) deallocate(Se,Se_fap,inc)
           end if APP_FREE_CLOUDY             ! Applicable CF/ICY
           if (REP_PROF) CYCLE REPEAT_PROF
           if (.NOT. REP_PROF) EXIT REPEAT_PROF
        end do REPEAT_PROF        ! Never should be reached
     end do PROFILE            ! End Profile loop
     !     ++++ deallocating all variables +++++
     deallocate(z_hr, tri, nd_rad_x, time_a)
     deallocate(TB_a, dbz_a, vel_a, sigma_a, code_a, tem_a, dpt_a, lwp_mic_a, &
          T_rs_a, p_rs_a, q_rs_a, dbz_hr_a, lin_hr_a, vel_hr_a, class_hr_a, &
          base_m_a,top_m_a, base_hr_a, iwv_mic_a, code_c_a, dd_a, vv_a)

     CLOSE(13)
  end do OVERDAYS

  return
end subroutine over_days
