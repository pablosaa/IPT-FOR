! -----------------------------------------
! Subroutine: LWC_HIGH_RES.F90
! -----
! Part of the main program IPT_KNR_BBC2.F90
! -----------------------------------------
! OUTPUT Variables:
!       deltaz_hr
!       lwc_hr
!       z_hr_lwc
!       dbz_hr_lwc
!       re_hr
! These variables are allocateds in this subroutine and must be deallocated
! in the caller main program.
!

subroutine lwc_high_res(z_lr_lwc, h_true_m, lwc, xn, c, d, dbz_hr, z_hr, n_cl) !, &
   !  deltaz_hr, lwc_hr, z_hr_lwc, dbz_hr_lwc, re_hr)
  use variables, only: deltaz_low, pi
  implicit none

  integer, intent(in) :: h_true_m(:), n_cl
  real, intent(in) :: z_lr_lwc(:), dbz_hr(:), z_hr(:)
  real(kind=8), intent(in) :: lwc(:), xn(:), c(:), d(:)


  ! Local Variables.
  real(kind=8)  :: deltaz_hr
  real(kind=8), allocatable :: lwc_hr(:), z_hr_lwc(:), &
       dbz_hr_lwc(:), re_hr(:)

  integer :: ss, p, n_valid, i_tot, i, j, k, lwc_check, allost
  integer, allocatable, dimension(:) :: i_valid
  real(kind=8) :: f1, f2, lwp
  real(kind=8), allocatable, dimension(:) :: z_hr_v, dbz_hr_v, lwc_hr_v
  real(kind=8), allocatable, dimension(:) :: zz_hr_v, re_hr_v, garbage
  real(kind=8) :: sum_lwc, sum_zz, lwp_zz

  i_tot = 1
  lwc_check = 0
  f1 = (pi*1000./6.)**(1./3.)
  f2 = EXP(-2*(.32**2))

  do i = 1, n_cl
     do j = 1, h_true_m(i)
        lwp = 10.**(lwc(i_tot)/10.)*deltaz_low
        sum_lwc = 0.
        sum_zz = 0.
        lwp_zz = lwp/1e3
        p = 1
        n_valid = COUNT((z_hr.GT.(z_lr_lwc(i_tot)-deltaz_low/2.)).AND.&
             (z_hr.LT.(z_lr_lwc(i_tot)+deltaz_low/2.)).AND.(dbz_hr.GT.-70.))

        VALID: if (n_valid.GT.0) then
           allocate(z_hr_v(n_valid), dbz_hr_v(n_valid), lwc_hr_v(n_valid), &
                zz_hr_v(n_valid), re_hr_v(n_valid), i_valid(n_valid), stat=allost)
           if (allost/=0) STOP 'Allocate failed in LWC_HIGH_RES subroutine.'
           
           INDEX: do k=1,size(z_hr)
              if ((z_hr(k).GT.(z_lr_lwc(i_tot)-deltaz_low/2.)).AND.&
                  (z_hr(k).LT.(z_lr_lwc(i_tot)+deltaz_low/2.)).AND.&
                  (dbz_hr(k).GT.-70.)) then
                           
                 i_valid(p) = k
                 p = p + 1
              end if
           end do INDEX
                     
           z_hr_v = z_hr(i_valid)
           dbz_hr_v = dbz_hr(i_valid)
           lwc_hr_v = 0.
           zz_hr_v = (10.**(dbz_hr_v/10.))*1e-18
           
           do k=1,n_valid
              if (i_valid(k).EQ.1) deltaz_hr = z_hr(i_valid(k))
              if (i_valid(k).GT.1) deltaz_hr = z_hr(i_valid(k))-z_hr(i_valid(k)-1)
              sum_lwc = sum_lwc + 10.**((dbz_hr_v(k)-c(i_tot))/(10.*d(i_tot)))*deltaz_hr
              sum_zz = sum_zz + SQRT(zz_hr_v(k))*deltaz_hr
           end do
           sum_zz = sum_zz**(1./3.)

           lwc_hr_v = 10.**((dbz_hr_v-c(i_tot))/(10.*d(i_tot)))*lwp/sum_lwc
           re_hr_v = (zz_hr_v**(1./6.))/(2.*(lwp_zz**(1./3.)))*sum_zz*f1*f2

           select case (lwc_check)
              case (0)
                 allocate(lwc_hr(n_valid), z_hr_lwc(n_valid), &
                      dbz_hr_lwc(n_valid), re_hr(n_valid))
                 lwc_hr = lwc_hr_v
                 z_hr_lwc = z_hr_v
                 dbz_hr_lwc = dbz_hr_v
                 re_hr = re_hr_v
              case (1)
                 allocate(garbage(size(lwc_hr)))
                 garbage = lwc_hr
                 deallocate(lwc_hr)
                 allocate(lwc_hr(n_valid+ss))
                 lwc_hr = [garbage, lwc_hr_v]
                 garbage = z_hr_lwc
                 deallocate(z_hr_lwc)
                 allocate(z_hr_lwc(n_valid+ss))
                 z_hr_lwc = [garbage, z_hr_v]
                 garbage = dbz_hr_lwc
                 deallocate(dbz_hr_lwc)
                 allocate(dbz_hr_lwc(n_valid+ss))
                 dbz_hr_lwc = [garbage, dbz_hr_v]
                 garbage = re_hr
                 deallocate(re_hr)
                 allocate(re_hr(n_valid+ss))
                 re_hr = [garbage, re_hr_v]
                 deallocate(garbage)
              case default
                 stop 'Error selecting LWC_CHECK in LWC_HIGH_RES subroutine.'
              end select
              lwc_check = 1
        end if VALID
        i_tot = i_tot + 1
        ss = n_valid
        deallocate(z_hr_v, dbz_hr_v, lwc_hr_v, zz_hr_v, re_hr_v, i_valid, &
             stat=allost)
        if (allost/=0) STOP 'Deallocate failed in LWC_HIGH_RES.'
           
     end do
  end do
  write(13,'(I5)') size(z_hr_lwc)
  write(13,'(20F12.7)') z_hr_lwc
  write(13,'(20F12.7)') dbz_hr_lwc
  write(13,'(20F12.7)') lwc_hr
  write(13,'(20F12.7)') re_hr*1e6
  write(13,'(2F12.7)') SUM(10.**((xn)/10.))*deltaz_low, SUM(lwc_hr)*deltaz_hr
  return
end subroutine lwc_high_res


