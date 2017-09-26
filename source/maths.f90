MODULE maths
  implicit none
  public :: chisqr_cvf
  private :: LN_FactX
  private :: Chi_Square
  private :: Chi_Square_Cumul
  private :: Bisec_ChiSqr_pdf
  !!=============
  public :: invert_LU
  private :: LUDCMP
  private :: LUBKSB
CONTAINS

  !=========================================
  !! PURPOSE:
  !! This function computes the cutoff value (v) such that:
  !!  Probability(X > v) = p
  !! where X is a random vaiable from Chi-square distribtion
  !! with df degrees of freedom.
  !! CATEGORY:
  !! Statistics.
  !! INPUT:
  !! p : A non-negative scalar, in the interval [0.0,1.0], of type
  !! float that specifies the probability of occurence or success.
  !! df : A positive scalar of type integer, that specifies the degrees of
  !! freedom of the Chi-square distribution.
  !! ---------------------------------------------

  function chisqr_cvf(p,df) result(v)
    implicit none
    real :: p
    integer :: df
    real(8) :: up, below, pdf, v

    if(p .LT. 0. .OR. p .GT. 1.) STOP 'p must be in the interval [0.0,1.0].'
    if(p .EQ. 0.) then
       v = 1.0d12
       return
    end if
    if(p .EQ. 1) then
       v = 0.0
       return
    end if
    select case (df)
    case (-HUGE(1):0)
       STOP 'Degrees of freedom must be positive.'
    case (1)
       up = 300.0
    case (2)
       up = 100.0
    case (3:5)
       up = 30.0
    case (6:14)
       up = 20.0
    case default
       up = 12.0
    end select
    below = 0.
    UPBEL: do
       call Chi_Square_Cumul(df,up,pdf)
       if (pdf .LT. (1.-p)) then
          below = up
          up = 2*up
       else
          exit UPBEL
       end if
    end do UPBEL
    call Bisec_ChiSqr_pdf((1.-p), df, up, below, v)
    return
  end function chisqr_cvf

  !***************************************************
  !* Series approximation subroutine LN(X!)          *
  !* Accuracy better then 6 places for x>=3          *
  !* Accuracy better than 12 places for x>10         *
  !* Advantage is that very large values of the      *
  !* argument can be used without fear of over flow. *
  !* x is the input (real,kind=8),                   *
  !* y is the output (real, kind=8).                 *
  !* ----------------------------------------------- *
  !* Reference: CRC Math Tables.                     *
  !* *************************************************
  Subroutine LN_FactX(x,y)
    implicit none
    real(8) :: x, x1, y
    x1 = 1.d0 / (x * x)
    y = (x + 0.5d0) * dlog(x) - x * (1.d0 - x1 / 12 + x1 * x1 / 360.d0 - x1 * x1 * x1 / 1260.d0 + x1 * x1 * x1 * x1 / 1680.d0)
    y = y + 0.918938533205d0
    return
  end Subroutine LN_FactX

  !******************************************************
  !* Chi-square function subroutine. This program takes *
  !* a given degree of freedom, m and value, x, and     *
  !* calculates the chi-square density distribution     *
  !* function value, y. Subroutine used: LN(X!).        *
  !* All variables are Real, Kind=8.                    *
  !* -------------------------------------------------- *
  !* Reference: Texas Instruments SR-51 owners Manual,  *
  !* 1974.                                              *
  !******************************************************
  Subroutine Chi_Square(m,x,y)
    implicit none
    integer :: m
    real(8) :: c, m1, x, y
    !Save X
    m1 = x
    !Perform calculation
    x = (dfloat(m) / 2.d0 - 1.d0)
    !Call LN(X!) subroutine
    call LN_FactX(x,y)
    x = m1
    c = -x / 2 + (dfloat(m) / 2.d0 - 1.d0) * dlog(x) - (dfloat(m) / 2.d0) * dlog(2.d0) - y
    y = dexp(c)
    return
  end Subroutine Chi_Square

  !*****************************************************************
  !*       Chi-square cumulative distribution subroutine           *
  !* ------------------------------------------------------------- *
  !* The program is fairly accurate and calls upon the chi-square  *
  !* probability density function subroutine. The input parameter  *
  !* is m, the number of degrees of freedom. Also required is the  *
  !* ordinate value, x. The subroutine returns y, the cummulative  *
  !* distribution integral from 0 to x. This program also requires *
  !* an accuracy parameter, e, to determine the level of summation.*
  !* Calls Chi_Square() subroutine.                                *
  !* All variables are Real, Kind=8, except m and m2 (integer).    *
  !* ------------------------------------------------------------- *
  !* Reference: Hewlett-Packard statistics programs, 1974.         *
  !*****************************************************************
  Subroutine Chi_Square_Cumul(m,x,y)
    implicit none
    ! Labels: 100, 200
    integer :: m, m2
    real(8) x, x2, y, y1
    real(8), parameter :: e = 1e-6

    y1 = 1; x2 = x; m2 = m + 2
    x2 = x2 / m2
    OVER: DO
       y1 = y1 + x2
       if (x2 < e) EXIT OVER
       m2 = m2 + 2;
       ! This form is used to avoid overflow
       x2 = x2 * (x / m2);
       ! Loop to continue sum  
    END DO OVER
    ! Obtain y, the probability density function
    call Chi_Square(m,x,y)
    y = (y1 * y * 2) * (x / m)
    return
  end Subroutine Chi_Square_Cumul


  !*******************************************************
  !* PURPOSE:
  !* This function computes the cutoff value x such that the probability
  !* of an abservation from the given distribution, less than x, is p.
  !* up and low are the upper and lower limits for x, respectively.
  !* df is the degrees of freedom.
  !* p : real.
  !* df: integer.
  !* up, low : real kind 8.
  !* mid : output vriable, real kind 8.
  !*******************************************************
  subroutine Bisec_ChiSqr_pdf(p, df, up, low, mid)
    implicit none
    real(8), parameter :: del = 1.0e-6
    real :: p
    real(8) :: mid, up, low, z
    integer :: df, count

    mid = low + (up-low)*p
    CNT: do count=1,100
       if (ABS(up-low) .GT. del*mid) then
          call Chi_Square_Cumul(df,mid,z)

          if (z .GT. p) then
             up = mid
          else
             low = mid
          end if
          mid = (up + low)/2.
       else
          exit CNT
       endif
    end do CNT
    return
  end subroutine Bisec_ChiSqr_pdf
  !! ====================================================
  !!
  !! ====================================================
  ! ============================================================
  ! Subroutine for inversion of a real square matrix by the
  ! method LU decomposition.
  ! Parameters:
  ! mat     : NxN Matrix to be inverted. (kind=8)
  ! mat_inv : Inverted Matrix. (kind=8)
  ! nn      : rank of the matrix.
  ! singu   : Flag to verify if the inverted matrix exist or
  !           a singularity is found. 1=singular.
  ! 

  subroutine invert_LU(mat, mat_inv,nn, singu)
    implicit none
    integer, intent(in) :: nn
    real(kind=8), intent(in), dimension(nn,nn) :: mat
    real(kind=8), intent(out), dimension(nn,nn) :: mat_inv
    integer, intent(out) :: singu

    ! Local variables.
    integer :: i,j,k, D, rc
    real(kind=8), allocatable, dimension(:,:) :: A,Y
    real :: deter
    integer, allocatable, dimension(:) :: indx

    singu = 0
    allocate(A(nn,nn), Y(nn,nn), indx(nn)) ; Y = 0
    forall(i=1:nn,.TRUE.)
       Y(i,i) = 1
    end forall
    A = mat
    !call LU decomposition routine (only once)
    call ludcmp(A,nn,indx,D,rc)
    deter = real(D)
    do i=1,nn
       deter = deter*A(i,i)
    end do
    !if(abs(deter) .LT. sqrt(tiny(0.))) write(*,*) 'Matrix has not an invers. (det=',deter,')'

    if (rc .EQ. 1) then
       singu = 1
    else
       !call solver if previous return code is ok
       !to obtain inverse of A one column at a time
       do j=1,nn
          call LUBKSB(A,nn,indx,Y(:,j))
       end do
       !the inverse matrix is now in matrix Y
       !the original matrix A is destroyed
       mat_inv = Y
    end if
    deallocate(A,Y,indx)
    return
  end subroutine invert_LU


  !  ***************************************************************
  !  * Given an N x N matrix A, this routine replaces it by the LU *
  !  * decomposition of a rowwise permutation of itself. A and N   *
  !  * are input. INDX is an output vector which records the row   *
  !  * permutation effected by the partial pivoting; D is output   *
  !  * as -1 or 1, depending on whether the number of row inter-   *
  !  * changes was even or odd, respectively. This routine is used *
  !  * in combination with LUBKSB to solve linear equations or to  *
  !  * invert a matrix. Return code is 1, if matrix is singular.   *
  !  ***************************************************************
  Subroutine LUDCMP(A,N,INDX,D,CODE)
    implicit none
    integer, intent(in) :: N
    real(kind=8), intent(inout), dimension(N,N) :: A
    integer, intent(out), dimension(N) :: INDX
    integer, intent(out) :: D, CODE

    ! LOCAL VARIABLES
    integer, parameter :: NMAX=100
    real(kind=8), parameter :: TINY=1.5D-16  !PARAMETER(NMAX=100,TINY=1.5D-16)
    real(kind=8) ::  AMAX,DUM, SUM, VV(NMAX)
    integer :: I, J, K, IMAX

    D=1; CODE=0

    DO I=1,N
       AMAX=0.d0
       DO J=1,N
          IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
       END DO ! j loop
       IF(AMAX.LT.TINY) THEN
          CODE = 1
          RETURN
       END IF
       VV(I) = 1.d0 / AMAX
    END DO ! i loop

    DO J=1,N
       DO I=1,J-1
          SUM = A(I,J)
          DO K=1,I-1
             SUM = SUM - A(I,K)*A(K,J) 
          END DO ! k loop
          A(I,J) = SUM
       END DO ! i loop
       AMAX = 0.d0
       DO I=J,N
          SUM = A(I,J)
          DO K=1,J-1
             SUM = SUM - A(I,K)*A(K,J) 
          END DO ! k loop
          A(I,J) = SUM
          DUM = VV(I)*DABS(SUM)
          IF(DUM.GE.AMAX) THEN
             IMAX = I
             AMAX = DUM
          END IF
       END DO ! i loop  

       IF(J.NE.IMAX) THEN
          DO K=1,N
             DUM = A(IMAX,K)
             A(IMAX,K) = A(J,K)
             A(J,K) = DUM
          END DO ! k loop
          D = -D
          VV(IMAX) = VV(J)
       END IF

       INDX(J) = IMAX
       IF(DABS(A(J,J)) < TINY) A(J,J) = TINY

       IF(J.NE.N) THEN
          DUM = 1.d0 / A(J,J)
          DO I=J+1,N
             A(I,J) = A(I,J)*DUM
          END DO ! i loop
       END IF
    END DO ! j loop

    RETURN
  END subroutine LUDCMP


  !  ******************************************************************
  !  * Solves the set of N linear equations A . X = B.  Here A is     *
  !  * input, not as the matrix A but rather as its LU decomposition, *
  !  * determined by the routine LUDCMP. INDX is input as the permuta-*
  !  * tion vector returned by LUDCMP. B is input as the right-hand   *
  !  * side vector B, and returns with the solution vector X. A, N and*
  !  * INDX are not modified by this routine and can be used for suc- *
  !  * cessive calls with different right-hand sides. This routine is *
  !  * also efficient for plain matrix inversion.                     *
  !  ******************************************************************
  Subroutine LUBKSB(A,N,INDX,B)
    implicit none
    integer, intent(in) :: N
    REAL(kind=8), intent(in) ::  A(N,N)
    INTEGER, intent(in) :: INDX(N)
    real(kind=8), intent(inout) :: B(N)

    real(kind=8) :: SUM
    integer :: II, I, LL, J
    II = 0 

    DO I=1,N 
       LL = INDX(I)
       SUM = B(LL)
       B(LL) = B(I)
       IF(II.NE.0) THEN
          DO J=II,I-1
             SUM = SUM - A(I,J)*B(J)
          END DO ! j loop
       ELSE IF(SUM.NE.0.d0) THEN
          II = I
       END IF
       B(I) = SUM
    END DO ! i loop

    DO I=N,1,-1
       SUM = B(I)
       IF(I < N) THEN
          DO J=I+1,N
             SUM = SUM - A(I,J)*B(J)
          END DO ! j loop
       END IF
       B(I) = SUM / A(I,I)
    END DO ! i loop

    RETURN
  END subroutine LUBKSB


END MODULE MATHS
