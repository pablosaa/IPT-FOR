! Module containning ABS functions
! Part of the IPT main program
module absorption
  implicit none

  contains
    
    ! Function to calculate FAP optical thickness tau at height k (index counting from bottom of zgrid).
    ! Rayleight calculation
    ! ------------
    ! For use L93 must be included the flag 'L93' in the argument before-last.
    ! For use R98 must be included the flag 'R98' in the argument before-last.
    
    subroutine tau_calc_fap(z_final,T_final,p_final,q_final,lwc_final,mu,coeff_fap,char,tau)
      use variables, only : n_Tq, n_freq, freq, pi, n_freq_r, c_li
!      use absorption
      implicit none
      real, intent(in) :: z_final(:), coeff_fap(:,:), p_final(:)
      real(kind=8), intent(in) :: mu, T_final(:), q_final(:), lwc_final(:)
      character(len=*), intent(in) :: char        ! must be 'L93' or 'R98'.
      real(kind=8), intent(out) :: tau(:,:)

      ! local variables
      integer :: i, j
      real(kind=8), dimension(n_freq_r) :: absg, absc, AWV, AO2, AN2
      real(kind=8) :: T_mean, P_mean, q_mean, xp
      real :: deltaz
      real(kind=8) :: T_mean_reg, P_mean_reg, q_mean_reg
      real(kind=8), dimension(n_freq) :: tau_x
      real(kind=8), dimension(n_Tq-1,n_freq) :: abs_all

      tau = 0.
      abs_all = 0.
      do i = 1, n_Tq-1
         ! alles SI!!
         deltaz = z_final(n_Tq-i+1) - z_final(n_Tq-i)
         T_mean = (T_final(n_Tq-i+1) + T_final(n_Tq-i))/2.d0
         xp = -LOG(p_final(n_Tq-i+1)/p_final(n_Tq-i))/deltaz
         P_mean = -p_final(n_Tq-i)/xp*(EXP(-xp*deltaz)-1.d0)/deltaz
         q_mean = (q_final(n_Tq-i+1) + q_final(n_Tq-i))/2.d0
         ! gas absortion (FAP)
         absg = 0. ; absc = 0.
         if (z_final(n_Tq-i+1) .LT. 1e4) then
            T_mean_reg = T_mean-273.15
            q_mean_reg = q_mean*1000.
            P_mean_reg = P_mean/100.
            absg = coeff_fap(1,:) + &
                 T_mean_reg*coeff_fap(2,:) + &
                 (T_mean_reg**2)*coeff_fap(3,:) + &
                 (T_mean_reg**3)*coeff_fap(4,:) + &
                 (P_mean_reg)*coeff_fap(5,:) + &
                 (P_mean_reg**2)*coeff_fap(6,:) + &
                 (P_mean_reg**3)*coeff_fap(7,:) + &
                 (q_mean_reg)*coeff_fap(8,:) + &
                 (q_mean_reg**2)*coeff_fap(9,:) + &
                 (q_mean_reg**3)*coeff_fap(10,:)
         else
            select case (char)
            case ('L93')
               ! subroutines for L93.
               AWV = 0. ; AO2 = 0.
               AWV = abwvl93(q_mean, T_mean, P_mean, freq)
               AO2 = abo2l93(q_mean, T_mean, P_mean, freq)
               absg = AWV + AO2
            case ('R98')
               ! subroutines for R98.
               AWV = 0. ; AO2 = 0. ; AN2 = 0.
               AWV = abwvr98(q_mean*1000., T_mean, P_mean/100., freq)
               AO2 = abo2r98(T_mean, P_mean/100., q_mean*1000., freq)
               AN2 = absn2(T_mean, P_mean/100., freq)
               absg = (AWV + AO2 + AN2)/1000.
            case default
               STOP 'Error selecting l93 or r98 in TAU_CALC_FAP subroutine.' 
            end select
         end if

         ! Cloud absortion (Ulaby et al. 1982)
         ! (Liebe et al. 1993)
         absc = abliq(1e3*lwc_final(n_Tq-i), T_mean, freq)
         absc = absc/1000.
         abs_all(n_Tq-i,:) = absg + absc

         tau_x = 0.

         do j=1,i
            deltaz = z_final(n_Tq-j+1)-z_final(n_Tq-j)
            tau_x = tau_x+abs_all(n_Tq-j,:)*deltaz
         end do

         tau(n_Tq-i,:) = tau_x

      end do

      return
    end subroutine tau_calc_fap

    ! ------------------------
    ! Function: ABWVR98
    ! Original Subroutine IDL: ABWVR98.PRO
    ! Compute Absorption Coef in Atmosphere due to Water Vapor.
    ! Parameters:
    !          *   rho [g/m^3] ; water vapor density
    !          *   T [Kelvin]  ; Temperature
    !          *   P [mbar]    ; Presure       (0.1 to 1000)
    !          +   alpha [Nepers/Km]  ; Absorption coefficient.
    ! This function use the parameter FREQ from the module VARIABLES
    ! FREQ is the frequency in GHz with a valid range of 0 to 800
    ! REFERENCES: P.W. ROSENKRANZ, Radio Science V.33, pp. 919-928 (1998);
    !             V.34, P.1025 (1999).
    !
    ! LINE INTENSITIES SELECTION THRESHOLD=
    !     HALF OF CONTINUUM ABSORPTION AT 1000 MB.
    !     WIDTHS MEASURED AT 22, 183, 380 GHZ, OTHERS CALCULATED.
    !     A.BAUER ET AL.ASA WORKSHOP (SEPT. 1989) (380GHz).

    !   REVISION HISTORY-
    !    DATE- OCT.6, 1988  P.W.ROSENKRANZ - EQS AS PUBL. IN 1993.
    !          OCT.4, 1995  PWR- USE CLOUGH'S DEFINITION OF LOCAL LINE
    !                   CONTRIBUTION,  HITRAN INTENSITIES, ADD 7 LINES.
    !          OCT. 24, 95  PWR -ADD 1 LINE.
    !          JULY 7, 97   PWR -SEPARATE COEFF. FOR SELF-BROADENING,
    !                       REVISED CONTINUUM.
    !          DEC. 11, 98  PWR - ADDED COMMENTS

    elemental function abwvr98(rho, T , P, freq) result(alpha)
      implicit none
      real(kind=8), intent(in) :: rho, T, P
      real, intent(in) :: freq
      real(kind=8) :: alpha

      !local variables
      integer :: i, j, k
      real(kind=8) :: pvap, pda, den, TI, TI2, width, wsq, S, base
      real(kind=8) :: suma, con, res
      real(kind=8), dimension(2) :: DF
              ! Line Frequencies:
      real, parameter, dimension(15) :: FL = &
           (/22.2351, 183.3101, 321.2256, 325.1529, 380.1974, 439.1508,&
            443.0183, 448.0011, 470.8890, 474.6891, 488.4911, 556.9360,&
            620.7008, 752.0332, 916.1712/)
              ! Line Intensities at 300K:
      real, parameter, dimension(15) :: S1 = &
           (/.1310E-13, .2273E-11, .8036E-13, .2694E-11, .2438E-10,&
             .2179E-11, .4624E-12, .2562E-10, .8369E-12, .3263E-11,&
             .6659E-12, .1531E-08, .1707E-10, .1011E-08, .4227E-10/)
              ! T coeff. of Intensities:
      real, parameter, dimension(15) :: B2 = &
           (/2.144, .668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405,&
             3.597, 2.379, 2.852, .159, 2.391, .396, 1.441/)
              ! Air-Broadened width parameters at 300K:
      real, parameter, dimension(15) :: W3 = &
           (/.00281, .00281, .0023, .00278, .00287, .0021, .00186,&
             .00263, .00215, .00236, .0026, .00321, .00244, .00306, .00267/)
              ! T-Exponent of Air-Broadening:
      real, parameter, dimension(15) :: X = &
           (/.69, .64, .67, .68, .54, .63, .60, .66, .66, .65, .69, .69,&
             .71, .68, .70/)
              ! Self-Broadened Width Parameters at 300K
      real, parameter, dimension(15) :: WS = &
           (/.01349, .01491, .0108, .0135, .01541, .0090, .00788,&
             .01275, .00983, .01095, .01313, .01320, .01140, .01253, .01275/)
              ! T-Exponent of self-broadening:
      real, parameter, dimension(15) :: XS = &
           (/.61, .85, .54, .74, .89, .52, .50, .67, .65, .64, .72,&
             1.0, .68, .84, .78/)

      if(rho .LE. 0.) then
         alpha = 0.
      else
         pvap = rho*T/217.d0
         pda = P - pvap
         den = (3.335e16)*rho
         TI = 300.d0/T
         TI2 = TI**2.5
         ! Continuum terms
         con = ((5.43e-10)*pda*TI**3 + (1.8e-8)*pvap*TI**7.5)*pvap*freq*freq

         ! Add Resonances
         suma = 0.
         do i = 1, 15
            width = W3(i)*pda*TI**X(i) + WS(i)*pvap*TI**XS(i)
            wsq = width*width
            S = S1(i)*TI2*EXP(B2(i)*(1.-TI))
            DF(1) = freq - FL(i) ; DF(2) = freq + FL(i)
            ! Use clough's definition of local line contribution
            base = width/(562500.d0 + wsq)
            ! Do for positive and negative resonances
            res = 0.
            !forall(k=1:n_freq,j=1:2, abs(DF(j,k)) .LT. 750.)
            !forall(j=1:2, abs(DF(j)) .LT. 750.)
            do j=1,2
               if (ABS(DF(j)) .LT. 750.) res = res + width/(DF(j)*DF(j) + wsq) - base
            end do
            !end forall
            suma = suma + S*res*(freq/FL(i))*(freq/FL(i))
         end do
         alpha = (0.3183d-4)*den*suma + con
      end if
      return
    end function abwvr98


    ! ------------------------
    ! Function: ABO2R98
    ! Original Subroutine IDL: ABO2R98.PRO
    ! Return Absorption Coefficient due to Oxygen in air.
    ! Parameters:
    !          *   vapden [g/m^3] ; water vapor density (enters linewidth
    !                               calculation due to greater broadening
    !                               efficiency of H2O)
    !          *   Temp [Kelvin]  ; Temperature  (Uncertain, but belived to be
    !                               valid for atmosphere)
    !          *   Pres [mbar]    ; Presure       (3 to 1000)
    !          +   o2abs [Nepers/Km]  ; Absorption coefficient.
    ! This function use the parameter FREQ from the module VARIABLES
    ! FREQ is the frequency in GHz with a valid range of 0 to 900
    ! REFERENCES: P.W. ROSENKRANZ, Atmosphere Remote Sensing by Microwave
    !             Radiometry, Chap.2 and appendix. (M.A. Janssen,ed.,1993).
    !  H.J. Liebe et al, JQSRT V.48, PP.629-643 (1992).
    !  M.J. Schwartz, Ph.D. thesis, M.I.T. (1997).
    !  SUBMILLIMETER LINE INTENSITIES FROM HITRAN96.
    !  This version differs from Liebe's MPM92 in two significant respects:
    !  1. It uses the modification of the 1- line width temperature dependence
    !     recommended by Schwartz: (1/T).
    !  2. It uses the same temperature dependence (X) for submillimeter
    !  line widths as in the 60 GHz band: (1/T)^0.8
    
    elemental function abo2r98(temp, pres, vapden, freq) result(o2abs)
      use variables, only: pi
      implicit none
      real(kind=8), intent(in) :: temp, pres, vapden
      real, intent(in) :: freq
      real(kind=8) :: o2abs

      ! local variables
            ! LINES ARE ARRANGED 1-,1+,3-,3+,ETC. IN SPIN-ROTATION SPECTRUM
      real, parameter, dimension(40) :: F = &
          (/ 118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.5910,&
          59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002,&
          56.9682, 62.4112, 56.3634, 62.9980, 55.7838, 63.5685,&
          55.2214, 64.1278, 54.6712, 64.6789, 54.1300, 65.2241,&
          53.5957, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368,&
          52.0214, 67.3696, 51.5034, 67.9009, 368.4984, 424.7632,&
          487.2494, 715.3931, 773.8397, 834.1458/)

      real, parameter, dimension(40) :: S300 = &
           (/.2936E-14, .8079E-15, .2480E-14, .2228E-14,&
             .3351E-14, .3292E-14, .3721E-14, .3891E-14,&
             .3640E-14, .4005E-14, .3227E-14, .3715E-14,&
             .2627E-14, .3156E-14, .1982E-14, .2477E-14,&
             .1391E-14, .1808E-14, .9124E-15, .1230E-14,&
             .5603E-15, .7842E-15, .3228E-15, .4689E-15,&
             .1748E-15, .2632E-15, .8898E-16, .1389E-15,&
             .4264E-16, .6899E-16, .1924E-16, .3229E-16,&
             .8191E-17, .1423E-16, .6494E-15, .7083E-14,&
             .3025E-14, .1835E-14, .1158E-13, .3993E-14/)
      real, parameter, dimension(40) :: BE = &
           (/.009, .015, .083, .084, .212, .212, .391, .391, .626,&
           .626, .915, .915, 1.260, 1.260, 1.660, 1.665, 2.119,&
           2.115, 2.624, 2.625, 3.194, 3.194, 3.814, 3.814,&
           4.484, 4.484, 5.224, 5.224, 6.004, 6.004, 6.844,&
           6.844, 7.744, 7.744, .048, .044, .049, .145, .141, .145/)
      ! Width in MHz/MB
      real, parameter, dimension(40) :: W300 = &
           (/1.63, 1.646, 1.468, 1.449, 1.382, 1.360,&
             1.319, 1.297, 1.266, 1.248, 1.221, 1.207, 1.181, 1.171,&
             1.144, 1.139, 1.110, 1.108, 1.079, 1.078, 1.05, 1.05,&
             1.02, 1.02, 1.00, 1.00, .97, .97, .94, .94, .92, .92, .89,&
             .89, 1.92, 1.92, 1.92, 1.81, 1.81, 1.81/)

      real, parameter, dimension(40) :: Y300 = &
           (/-0.0233,  0.2408, -0.3486,  0.5227,&
           -0.5430,  0.5877, -0.3970,  0.3237, -0.1348,  0.0311,&
           0.0725, -0.1663,  0.2832, -0.3629,  0.3970, -0.4599,&
           0.4695, -0.5199,  0.5187, -0.5597,  0.5903, -0.6246,&
           0.6656, -0.6942,  0.7086, -0.7325,  0.7348, -0.7546,&
           0.7702, -0.7864,  0.8083, -0.8210,  0.8439, -0.8529,&
           0., 0., 0., 0., 0., 0./)

      real, parameter, dimension(40) :: V = &
           (/0.0079, -0.0978,  0.0844, -0.1273,&
           0.0699, -0.0776,  0.2309, -0.2825,  0.0436, -0.0584,&
           0.6056, -0.6619,  0.6451, -0.6759,  0.6547, -0.6675,&
           0.6135, -0.6139,  0.2952, -0.2895,  0.2654, -0.2590,&
           0.3750, -0.3680,  0.5085, -0.5002,  0.6206, -0.6091,&
           0.6526, -0.6393,  0.6640, -0.6475,  0.6729, -0.6545,&
           0., 0., 0., 0., 0., 0./)

      real, parameter :: WB300=0.56, X=0.8
      real(kind=8) :: TH, TH1, B, presWV, presDA, den
      real(kind=8) :: dens, dfnr, DF, Y, STR, SF1, SF2
      real(kind=8) :: suma
      integer :: k
      
      TH = 300.d0/temp
      TH1 = TH - 1.d0
      B = TH**X
      presWV = vapden*temp/217.d0
      presDA = pres - presWV

      den = 0.001*(presDA*B + 1.1*presWV*TH)
      dens = 0.001*(presDA + 1.1*presWV)*TH
      dfnr = WB300*den

      suma = 1.6e-17*freq*freq*dfnr/(TH*(freq*freq + dfnr*dfnr))

      do k = 1, 40
         if (k .EQ. 1) then 
            DF = W300(k)*dens
         else
            DF = W300(k)*den
         end if
         Y = 0.001*pres*B*(Y300(k) + V(k)*TH1)
         STR = S300(k)*EXP(-BE(k)*TH1)
         SF1 = (DF + (freq - F(k))*Y)/((freq - F(k))**2 + DF*DF)
         SF2 = (DF - (freq + F(k))*Y)/((freq + F(k))**2 + DF*DF)
         suma = suma + STR*(SF1 + SF2)*(freq/F(k))**2
      end do

      o2abs = 0.5034d12*suma*presDA*TH**3/pi

      return
    end function abo2r98

    ! ------------------------
    ! Function: ABSN2
    ! Original Subroutine IDL: ABSN2.PRO
    ! Return ABSORPTION COEFFICIENT DUE TO NITROGEN IN AIR (NEPER/KM)
    ! Parameters:
    !T = TEMPERATURE (K)
    !P = PRESSURE (MB)
    ! alpha = absortion coefficient (NEPER/KM)
    !freq = FREQUENCY (GHZ), from de module VARIABLES.
    elemental function absn2(T,P,freq) result(alpha)
      implicit none
      real(kind=8), intent(in) :: T, P
      real, intent(in) :: freq
      real(kind=8) :: alpha

      !Local Variables
      real(kind=8) :: TH

      TH = 300.d0/T
      alpha = 6.4d-14*P*P*freq*freq*TH**3.55

      return
    end function absn2


    ! ------------------------
    ! Function: ABLIQ
    ! Original Subroutine IDL: ABLIQ.PRO
    ! Return ABSORPTION IN NEPERS/KM BY SUSPENDED WATER DROPLETS
    ! Parameters:
    !           WATER [G/M**3]
    !           FREQ [GHZ]     (VALID FROM 0 TO 1000 GHZ)
    !           TEMP [KELVIN]
    !REFERENCES:
    !;LIEBE, HUFFORD AND MANABE, INT. J. IR & MM WAVES V.12, pp.659-675
    !;(1991);  Liebe et al, AGARD Conf. Proc. 542, May 1993.
    !;REVISION HISTORY:
    !;PWR 8/3/92   original version
    !;PWR 12/14/98 temp. dependence of EPS2 eliminated to agree
    !with MPM93
    !;pwr 2/27/02  use exponential dep. on T, eq. 2b instead of eq. 4a
    elemental function abliq(water, temp, freq) result(alpha)
      implicit none

      real(kind=8), intent(in) :: water, temp
      real, intent(in) :: freq
      real(kind=8) :: alpha

      ! Local Variables
      real(kind=8) :: THETA1, EPS0, EPS1, EPS2, FP, FS
      complex(kind=8) :: EPS, RE

!      if (water .LE. 0) then
       if (water .LE. 0.)  alpha = 0.
      !else

         THETA1 = 1.d0-300.d0/temp
         EPS0 = 77.66 - 103.3*THETA1
         EPS1 = .0671*EPS0
         EPS2 = 3.52                 ! from MPM93
         FP = 20.1*EXP(7.88*THETA1)  ! from eq. 2b
         FS = 39.8*FP
         EPS = (EPS0 - EPS1)/CMPLX(1.d0,freq/FP) + &
              (EPS1 - EPS2)/CMPLX(1.d0,freq/FS) + EPS2
         RE = (EPS - 1.d0)/(EPS + 2.d0)
         alpha = -.06286*AIMAG(RE)*freq*water
      !end if
      return
    end function abliq


    ! ------------------------
    ! Function: ABWVl93
    ! Original Subroutine IDL: ABWVl93.PRO
    ! Return WATER VAPOR ABSORPTION MODEL
    ! Parameters:
    !           RHO : water vapor [kgm-3]
    !             T : temperature [K]
    !            PP : pressure [Pa]
    !         ALPHA : absorption coefficient [1/m]
    ! Description:
    !SOURCE: MPM93, June 1993
    !--------
    !       ATMOSPHERIC ATTENUATION AND DELAY RATES UP TO 1000 GHz
    !       Hans J. Liebe     (303-497-3310)
    !       George A. Hufford (       -3457)
    !       Michael G. Cotton (       -7346)
    !       Institute for Telecommunication Sciences
    !       NTIA/ITS.S3
    !       325 BROADWAY
    !       Boulder, CO  80303,  USA
    !
    !       FAX   :  (303) 497-5993 (ITS), 497-3680 (ITS.S2)
    !       E-Mail:  HLIEBE@NTIA.ITS.BLDRDOC.GOV
    !--------
    elemental function abwvl93(rho, T, PP, freq) result(alpha)
      use variables, only: F0water_l93, B_l93
      implicit none
      real(kind=8), intent(in) :: rho, T , PP
      real, intent(in) :: freq
      real(kind=8) :: alpha

      ! Local Variables
      integer :: i
      real(kind=8) :: V, E, Pb, P, GAMH, GAMD2, S, DELH
      complex(kind=8) :: ZN, ZF

      V  = 300.d0 / T                            ! norm. temperature
      E  = rho * T * 461.52 * 1.0d-2           ! water vapor pressure hPa
      Pb = PP * 0.01                           ! pressure in hPa
      P  = Pb - E                              ! dry air pressure

      ZN = (0., 0.)
      do i = 1, 35
         GAMH = 0.
         S = B_l93(1,i)*E*(V**3.5)*EXP(B_l93(2,i)*(1. - V))
         !Doppler approximation
         GAMH = B_l93(3,i)*(P*(V**B_l93(5,i)) + &
                B_l93(4,i)*E*(V**B_l93(6,i)))*1.e-3
         GAMD2 = 1e-12/V*(1.46*F0water_l93(i))**2
         GAMH = 0.535*GAMH + (0.217*GAMH**2 + GAMD2)**0.5
         DELH = 0.
         ZF = freq/F0water_l93(i)*(&
              CMPLX(1.,-DELH)/CMPLX(F0water_l93(i)-freq,-1.*GAMH) - &
              CMPLX(1.,DELH)/CMPLX(F0water_l93(i)+freq,GAMH))
         if (i .EQ. 1) ZN = S*ZF
         ZN = ZN + S*ZF
      end do
      alpha = .000182*freq*AIMAG(ZN)/4.343
      return
    end function abwvl93


    ! ------------------------
    ! Function: ABO2l93
    ! Original Subroutine IDL: ABO2l93.PRO
    ! Return Absorption Coefficient
    ! Parameters:
    !           RHO : water vapor [kgm-3]
    !             T : temperature [K]
    !            PP : pressure [Pa]
    !         ALPHA : absorption coefficient [1/m]
    ! Description:
    !OXYGEN ABSORPTION MODEL
    !SOURCE: MPM93, June 1993
    !-------------
    !       ATMOSPHERIC ATTENUATION AND DELAY RATES UP TO 1000 GHz
    !       Hans J. Liebe     (303-497-3310)
    !       George A. Hufford (       -3457)
    !       Michael G. Cotton (       -7346)
    !       Institute for Telecommunication Sciences
    !       NTIA/ITS.S3
    !       325 BROADWAY
    !       Boulder, CO  80303,  USA
    !
    !       FAX   :  (303) 497-5993 (ITS), 497-3680 (ITS.S2)
    !       E-Mail:  HLIEBE@NTIA.ITS.BLDRDOC.GOV
    !-------------
    elemental function abo2l93(rho, T, PP, freq) result(alpha)
      use variables, only: F0O2_l93, A_l93
      implicit none
      real(kind=8), intent(in) :: rho, T, PP
      real, intent(in) :: freq
      real(kind=8) :: alpha

      ! Local Variable
      integer :: i, i_num
      real(kind=8) :: V, E, Pb, P, GAMMA, DELTA, S, So, GAMMAo, Sn, ALPHA1, ALPHA2
      complex(kind=8) :: ZF, ZN, ZFo, ZFn

      V = 300./T                            ! norm. temperature
      E = RHO*T*461.52*1.0d-2               ! water vapor pressure hPa
      Pb = PP*0.01                          ! pressure in hPa
      P = Pb - E                            ! dry air pressure

      do i = 1,44
         GAMMA = 0.
         S = A_l93(1,i)*P*V**3*EXP(A_l93(2,i)*(1. - V))*1.e-6
         GAMMA = A_l93(3,i)*(P*V**(0.8 - A_l93(4,i)) + 1.1*E*V)*1.e-3
         GAMMA = (GAMMA**2 + (25*0.6e-4)**2)**0.5
         DELTA = (A_l93(5,i) + A_l93(6,i)*V)*(P + E)*(V**0.8)*1.e-3
         ZF = freq/F0O2_l93(i)*(&
              CMPLX(1.,-DELTA)/CMPLX(F0O2_l93(i)-freq,-GAMMA) - &
              CMPLX(1.,DELTA)/CMPLX(F0O2_l93(i)+freq,GAMMA))
         if (i .EQ. 0) ZN = S*ZF
         ZN = ZN + S*ZF
      end do

      ! OXYGEN LINE ABSORPTION
      ALPHA1 = .000182*freq*AIMAG(ZN)/4.343
      ! Cannot be less than 0.  M
      if (ALPHA1 .LT. 0.) ALPHA1 = 0.
      ! ojo i_num = WHERE(ALPHA LT 0.)
      ! ojo  IF i_num(0) GT -1 THEN ALPHA(i_num) = 0.

      ! DRY AIR CONTINUUMM
      So = 6.14d-5*P*(V**2)
      GAMMAo = 0.56d-3*(P + E)*(V**0.8)
      ZFo = -1.d0*freq/CMPLX(freq,GAMMAo)
      Sn = 1.40e-12*(P**2)*(V**3.5)
      ZFn = CMPLX(0.,freq/(1.93e-5*(freq**1.5)+1.))
      ZN = So*ZFo + Sn*ZFn

      ! NONRESONAT DRY AIR ABSORPTION
      ALPHA2 = .000182*freq*AIMAG(ZN)/4.343

      alpha = ALPHA1 + ALPHA2

      return
    end function abo2l93




    ! ------------------------
    ! Function: TB_CALC_PL
    ! Original Subroutine IDL: TB_CALC_PL.PRO
    ! Description: Calculate Brightness temperatures without scattering
    !              according to Simmer (94) pp. 87-91 (alpha=1, no scattering)
    !              Planck/thermodynamic conform (28.05.03); UL
    ! Parameters:
    !             T [Kelvin] : Temperature
    !             tau        :
    !             mu         :
    ! This routine use the values FREQ and N_FREQ from the module VARIABLES.
    !-------------------
    subroutine tb_calc_pl(T,tau,mu,TB)
      use variables, only: n_Tq, n_freq, freq, pi, h, c_li, kB
      implicit none
      real(kind=8), intent(in) :: T(:), tau(:,:), mu
      real(kind=8), intent(out), dimension(:) :: TB
      ! Local Variables
      real(kind=8), dimension(n_freq) :: freq_si, lambda_si, IN, tau_bot, tau_top
      real(kind=8), dimension(n_freq) :: delta_tau, A, B, T_pl2, T_pl1, diffe
      integer :: i

      freq_si = 1d9*freq
      lambda_si = c_li/freq_si
      IN = 2.73         ! cosmic background
      IN = (2.*h*freq_si)/(lambda_si**2)/(EXP(h*freq_si/kB/IN) - 1.d0)

      TB = 0.
      tau_top = 0.
      tau_bot = tau(n_Tq-1,:)
      do i=1,n_Tq-1
         if(i .GT. 1) then
            tau_top = tau(n_Tq-i+1,:)
            tau_bot = tau(n_Tq-i,:)
         endif
         delta_tau = tau_bot - tau_top
         A = 1.d0 - EXP(-delta_tau/mu)
         B = delta_tau - mu*(1.d0 - EXP(-delta_tau/mu))
         T_pl2 = (2.d0*h*freq_si/(lambda_si**2))/(EXP(h*freq_si/kB/T(n_Tq-i)) - 1.d0)
         T_pl1 = (2.d0*h*freq_si/(lambda_si**2))/(EXP(h*freq_si/kB/T(n_Tq-i+1)) - 1.d0)
         diffe = (T_pl2 - T_pl1)/delta_tau
         IN = IN*EXP(-delta_tau/mu) + T_pl1*A + diffe*B
      end do
       
      TB = (h*freq_si/kB)/LOG((2.d0*h*freq_si/IN/(lambda_si**2)) + 1.d0)
      return
    end subroutine tb_calc_pl
  end module absorption



!!$
!!$!===============================================================
!!$    function tb_calc_pl(T,tau,mu) result(TB)
!!$      use variables, only: n_Tq, n_freq, freq, pi, h, c_li, kB
!!$      implicit none
!!$      real(kind=8), intent(in) :: T(:), tau(:,:), mu
!!$      real(kind=8), dimension(n_freq) :: TB
!!$      ! Local Variables
!!$      real(kind=8), dimension(n_freq) :: freq_si, lambda_si, IN, tau_bot, tau_top
!!$      real(kind=8), dimension(n_freq) :: delta_tau, A, B, T_pl2, T_pl1, diffe
!!$      integer :: i
!!$
!!$      freq_si = 1e9*freq
!!$      lambda_si = c_li/freq_si
!!$      IN = 2.73         ! cosmic background
!!$      IN = (2.*h*freq_si)/(lambda_si**2)/(EXP(h*freq_si/kB/IN) - 1.)
!!$
!!$      TB = 0.
!!$      tau_top = 0.
!!$      tau_bot = tau(n_Tq-1,:)
!!$      do i=1,n_Tq-1
!!$         if(i .GT. 1) then
!!$            tau_top = tau(n_Tq-i+1,:)
!!$            tau_bot = tau(n_Tq-i,:)
!!$         endif
!!$         delta_tau = tau_bot - tau_top
!!$         A = SPREAD(1.,1,n_freq) - EXP(-delta_tau/mu)
!!$         B = delta_tau - mu*(1. - EXP(-delta_tau/mu))
!!$         T_pl2 = (2.*h*freq_si)/(lambda_si**2)/(EXP(h*freq_si/kB/T(n_Tq-i)) - 1.)
!!$         T_pl1 = (2.*h*freq_si)/(lambda_si**2)/(EXP(h*freq_si/kB/T(n_Tq-i+1)) - 1.)
!!$         diffe = (T_pl2 - T_pl1)/delta_tau
!!$         IN = IN*EXP(-delta_tau/mu) + T_pl1*A + diffe*B
!!$      end do
!!$       
!!$      TB = (h*freq_si/kB)/LOG((2.*h*freq_si/IN/(lambda_si**2)) + 1.)
!!$      return
!!$    end function tb_calc_pl
!!$  end module absorption


