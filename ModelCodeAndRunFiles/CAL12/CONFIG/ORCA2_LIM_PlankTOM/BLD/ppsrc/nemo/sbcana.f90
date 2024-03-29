/* Copyright (C) 1991-2012 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http://www.gnu.org/licenses/>.  */


/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */

/* We do support the IEC 559 math functionality, real and complex.  */

/* wchar_t uses ISO/IEC 10646 (2nd ed., published 2011-03-15) /
   Unicode 6.0.  */

/* We do not support C11 <threads.h>.  */

MODULE sbcana
   !!======================================================================
   !!                       ***  MODULE  sbcana  ***
   !! Ocean forcing:  analytical momentum, heat and freshwater forcings
   !!=====================================================================
   !! History :  3.0   ! 2006-06  (G. Madec)  Original code
   !!            3.2   ! 2009-07  (G. Madec)  Style only
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_ana  : set an analytical ocean forcing
   !!   sbc_gyre : set the GYRE configuration analytical forcing
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_fortran

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_ana    ! routine called in sbcmod module
   PUBLIC   sbc_gyre   ! routine called in sbcmod module

   !                       !!* Namelist namsbc_ana *
   INTEGER  ::   nn_tau000  ! nb of time-step during which the surface stress
   !                        ! increase from 0 to its nominal value 
   REAL(wp) ::   rn_utau0   ! constant wind stress value in i-direction
   REAL(wp) ::   rn_vtau0   ! constant wind stress value in j-direction
   REAL(wp) ::   rn_qns0    ! non solar heat flux
   REAL(wp) ::   rn_qsr0    !     solar heat flux
   REAL(wp) ::   rn_emp0    ! net freshwater flux
   
   !! * Substitutions
   !!----------------------------------------------------------------------
   !!                    ***  domzgr_substitute.h90   ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
   !!      factors depending on the vertical coord. used, using CPP macro.
   !!----------------------------------------------------------------------
   !! History :  1.0  !  2005-10  (A. Beckmann, G. Madec) generalisation to all coord.
   !!            3.1  !  2009-02  (G. Madec, M. Leclair)  pure z* coordinate
   !!----------------------------------------------------------------------

! z- or s-coordinate (1D or 3D + no time dependency) use reference in all cases







! This part should be removed one day ...
! ... In that case all occurence of the above statement functions
!     have to be replaced in the code by xxx_n

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: domzgr_substitute.h90 4488 2014-02-06 10:43:09Z rfurner $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!                   ***  vectopt_loop_substitute  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute the inner loop start/end indices with CPP macro
   !!                allow unrolling of do-loop (useful with vector processors)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.7 , NEMO Consortium (2014)
   !! $Id: vectopt_loop_substitute.h90 4990 2014-12-15 16:42:49Z timgraham $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: sbcana.F90 4624 2014-04-28 12:09:03Z acc $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_ana( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_ana ***
      !!              
      !! ** Purpose :   provide at each time-step the ocean surface boundary
      !!              condition, i.e. the momentum, heat and freshwater fluxes.
      !!
      !! ** Method  :   Constant and uniform surface forcing specified from
      !!              namsbc_ana namelist parameters. All the fluxes are time
      !!              independant except the stresses which increase from zero
      !!              during the first nn_tau000 time-step 
      !!
      !! ** Action  : - set the ocean surface boundary condition, i.e.  
      !!                   utau, vtau, taum, wndm, qns, qsr, emp, sfx
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !
      INTEGER  ::   ios                   ! local integer
      REAL(wp) ::   zrhoa  = 1.22_wp      ! air density kg/m3
      REAL(wp) ::   zcdrag = 1.5e-3_wp    ! drag coefficient
      REAL(wp) ::   zfact, ztx            ! local scalars
      REAL(wp) ::   zcoef, zty, zmod      !   -      -
      !!
      NAMELIST/namsbc_ana/ nn_tau000, rn_utau0, rn_vtau0, rn_qns0, rn_qsr0, rn_emp0
      !!---------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         !
         REWIND( numnam_ref )              ! Namelist namsbc_ana in reference namelist : Analytical surface fluxes
         READ  ( numnam_ref, namsbc_ana, IOSTAT = ios, ERR = 901)
901      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_ana in reference namelist', lwp )

         REWIND( numnam_cfg )              ! Namelist namsbc_ana in configuration namelist : Analytical surface fluxes
         READ  ( numnam_cfg, namsbc_ana, IOSTAT = ios, ERR = 902 )
902      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_ana in configuration namelist', lwp )
         IF(lwm) WRITE ( numond, namsbc_ana )
         !
         IF(lwp) WRITE(numout,*)' '
         IF(lwp) WRITE(numout,*)' sbc_ana : Constant surface fluxes read in namsbc_ana namelist'
         IF(lwp) WRITE(numout,*)' ~~~~~~~ '
         IF(lwp) WRITE(numout,*)'              spin up of the stress  nn_tau000 = ', nn_tau000, ' time-steps'
         IF(lwp) WRITE(numout,*)'              constant i-stress      rn_utau0  = ', rn_utau0 , ' N/m2'
         IF(lwp) WRITE(numout,*)'              constant j-stress      rn_vtau0  = ', rn_vtau0 , ' N/m2'
         IF(lwp) WRITE(numout,*)'              non solar heat flux    rn_qns0   = ', rn_qns0  , ' W/m2'
         IF(lwp) WRITE(numout,*)'              solar heat flux        rn_qsr0   = ', rn_qsr0  , ' W/m2'
         IF(lwp) WRITE(numout,*)'              net heat flux          rn_emp0   = ', rn_emp0  , ' Kg/m2/s'
         !
         nn_tau000 = MAX( nn_tau000, 1 )     ! must be >= 1
         !
         utau(:,:) = rn_utau0
         vtau(:,:) = rn_vtau0
         taum(:,:) = SQRT ( rn_utau0 * rn_utau0 + rn_vtau0 * rn_vtau0 )
         wndm(:,:) = SQRT ( taum(1,1) /  ( zrhoa * zcdrag ) )
         !
         emp (:,:) = rn_emp0
         sfx (:,:) = 0.0_wp
         qns (:,:) = rn_qns0 - emp(:,:) * sst_m(:,:) * rcp      ! including heat content associated with mass flux at SST
         qsr (:,:) = rn_qsr0
         !         
      ENDIF

      IF( MOD( kt - 1, nn_fsbc ) == 0 ) THEN 
         !
         IF( kt <= nn_tau000 ) THEN       ! Increase the stress to its nominal value 
            !                             ! during the first nn_tau000 time-steps
            zfact = 0.5 * (  1. - COS( rpi * REAL( kt, wp ) / REAL( nn_tau000, wp ) )  )
            zcoef = 1. / ( zrhoa * zcdrag ) 
            ztx   = zfact * rn_utau0
            zty   = zfact * rn_vtau0
            zmod  = SQRT( ztx * ztx + zty * zty )
            utau(:,:) = ztx
            vtau(:,:) = zty
            taum(:,:) = zmod
            zmod = SQRT( zmod * zcoef )
            wndm(:,:) = zmod
         ENDIF
         !                                ! update heat and fresh water fluxes
         !                                ! as they may have been changed by sbcssr module
         emp (:,:) = rn_emp0              ! NB: qns changes with SST if emp /= 0
         sfx (:,:) = 0._wp
         qns (:,:) = rn_qns0 - emp(:,:) * sst_m(:,:) * rcp
         qsr (:,:) = rn_qsr0
         !
      ENDIF
      !
   END SUBROUTINE sbc_ana


   SUBROUTINE sbc_gyre( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_ana ***
      !!              
      !! ** Purpose :   provide at each time-step the GYRE surface boundary
      !!              condition, i.e. the momentum, heat and freshwater fluxes.
      !!
      !! ** Method  :   analytical seasonal cycle for GYRE configuration.
      !!                CAUTION : never mask the surface stress field !
      !!
      !! ** Action  : - set the ocean surface boundary condition, i.e.   
      !!                   utau, vtau, taum, wndm, qns, qsr, emp, sfx
      !!
      !! Reference : Hazeleger, W., and S. Drijfhout, JPO, 30, 677-695, 2000.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt          ! ocean time step
      !!
      INTEGER  ::   ji, jj                 ! dummy loop indices
      INTEGER  ::   zyear0                 ! initial year 
      INTEGER  ::   zmonth0                ! initial month
      INTEGER  ::   zday0                  ! initial day
      INTEGER  ::   zday_year0             ! initial day since january 1st
      REAL(wp) ::   ztau     , ztau_sais   ! wind intensity and of the seasonal cycle
      REAL(wp) ::   ztime                  ! time in hour
      REAL(wp) ::   ztimemax , ztimemin    ! 21th June, and 21th decem. if date0 = 1st january
      REAL(wp) ::   ztimemax1, ztimemin1   ! 21th June, and 21th decem. if date0 = 1st january
      REAL(wp) ::   ztimemax2, ztimemin2   ! 21th June, and 21th decem. if date0 = 1st january
      REAL(wp) ::   ztaun                  ! intensity
      REAL(wp) ::   zemp_s, zemp_n, zemp_sais, ztstar
      REAL(wp) ::   zcos_sais1, zcos_sais2, ztrp, zconv, t_star
      REAL(wp) ::   zsumemp, zsurf
      REAL(wp) ::   zrhoa  = 1.22         ! Air density kg/m3
      REAL(wp) ::   zcdrag = 1.5e-3       ! drag coefficient
      REAL(wp) ::   ztx, zty, zmod, zcoef ! temporary variables
      REAL(wp) ::   zyydd                 ! number of days in one year
      !!---------------------------------------------------------------------
      zyydd = REAL(nyear_len(1),wp)

      ! ---------------------------- !
      !  heat and freshwater fluxes  !
      ! ---------------------------- !
      !same temperature, E-P as in HAZELEGER 2000

      zyear0     =   ndate0 / 10000                             ! initial year
      zmonth0    = ( ndate0 - zyear0 * 10000 ) / 100            ! initial month
      zday0      =   ndate0 - zyear0 * 10000 - zmonth0 * 100    ! initial day betwen 1 and 30
      zday_year0 = ( zmonth0 - 1 ) * 30.+zday0                  ! initial day betwen 1 and 360

      ! current day (in hours) since january the 1st of the current year
      ztime = REAL( kt ) * rdt / (rmmss * rhhmm)   &       !  total incrementation (in hours)
         &      - (nyear  - 1) * rjjhh * zyydd             !  minus years since beginning of experiment (in hours)

      ztimemax1 = ((5.*30.)+21.)* 24.                      ! 21th june     at 24h in hours
      ztimemin1 = ztimemax1 + rjjhh * zyydd / 2            ! 21th december        in hours
      ztimemax2 = ((6.*30.)+21.)* 24.                      ! 21th july     at 24h in hours
      ztimemin2 = ztimemax2 - rjjhh * zyydd / 2            ! 21th january         in hours
      !                                                    ! NB: rjjhh * zyydd / 4 = one seasonal cycle in hours

      ! amplitudes
      zemp_S    = 0.7       ! intensity of COS in the South
      zemp_N    = 0.8       ! intensity of COS in the North
      zemp_sais = 0.1
      zTstar    = 28.3      ! intemsity from 28.3 a -5 deg

      ! 1/2 period between 21th June and 21th December and between 21th July and 21th January
      zcos_sais1 = COS( (ztime - ztimemax1) / (ztimemin1 - ztimemax1) * rpi ) 
      zcos_sais2 = COS( (ztime - ztimemax2) / (ztimemax2 - ztimemin2) * rpi )

      ztrp= - 40.e0        ! retroaction term on heat fluxes (W/m2/K)
      zconv = 3.16e-5      ! convertion factor: 1 m/yr => 3.16e-5 mm/s
      DO jj = 1, jpj
         DO ji = 1, jpi
            ! domain from 15 deg to 50 deg between 27 and 28  degC at 15N, -3
            ! and 13 degC at 50N 53.5 + or - 11 = 1/4 period :
            ! 64.5 in summer, 42.5 in winter
            t_star = zTstar * ( 1 + 1. / 50. * zcos_sais2 )                &
               &                    * COS( rpi * (gphit(ji,jj) - 5.)               &
               &                    / ( 53.5 * ( 1 + 11 / 53.5 * zcos_sais2 ) * 2.) )
            ! 23.5 deg : tropics
            qsr (ji,jj) =  230 * COS( 3.1415 * ( gphit(ji,jj) - 23.5 * zcos_sais1 ) / ( 0.9 * 180 ) )
            qns (ji,jj) = ztrp * ( tsb(ji,jj,1,jp_tem) - t_star ) - qsr(ji,jj)
            IF( gphit(ji,jj) >= 14.845 .AND. 37.2 >= gphit(ji,jj) ) THEN    ! zero at 37.8 deg, max at 24.6 deg
               emp  (ji,jj) =   zemp_S * zconv   &
                  &         * SIN( rpi / 2 * (gphit(ji,jj) - 37.2) / (24.6 - 37.2) )  &
                  &         * ( 1 - zemp_sais / zemp_S * zcos_sais1)
            ELSE
               emp (ji,jj) =  - zemp_N * zconv   &
                  &         * SIN( rpi / 2 * (gphit(ji,jj) - 37.2) / (46.8 - 37.2) )  &
                  &         * ( 1 - zemp_sais / zemp_N * zcos_sais1 )
            ENDIF
         END DO
      END DO

      ! Compute the emp flux such as its integration on the whole domain at each time is zero
      IF( nbench /= 1 ) THEN
         zsumemp = GLOB_SUM( emp(:,:) ) 
         zsurf   = GLOB_SUM( tmask(:,:,1) ) 
         ! Default GYRE configuration
         zsumemp = zsumemp / zsurf
      ELSE
         ! Benchmark GYRE configuration (to allow the bit to bit comparison between Mpp/Mono case)
         zsumemp = 0.e0   ;    zsurf = 0.e0
      ENDIF

      ! freshwater (mass flux) and update of qns with heat content of emp
      emp (:,:) = emp(:,:) - zsumemp * tmask(:,:,1)        ! freshwater flux (=0 in domain average)
      sfx (:,:) = 0.0_wp                                   ! no salt flux
      qns (:,:) = qns(:,:) - emp(:,:) * sst_m(:,:) * rcp   ! evap and precip are at SST


      ! ---------------------------- !
      !       momentum fluxes        !
      ! ---------------------------- !
      ! same wind as in Wico
      !test date0 : ndate0 = 010203
      zyear0  =   ndate0 / 10000
      zmonth0 = ( ndate0 - zyear0 * 10000 ) / 100
      zday0   =   ndate0 - zyear0 * 10000 - zmonth0 * 100
      !Calculates nday_year, day since january 1st
      zday_year0 = (zmonth0-1)*30.+zday0

      !accumulates days of previous months of this year
      ! day (in hours) since january the 1st
      ztime = FLOAT( kt ) * rdt / (rmmss * rhhmm)  &  ! incrementation in hour
         &     - (nyear - 1) * rjjhh * zyydd          !  - nber of hours the precedent years
      ztimemax = ((5.*30.)+21.)* 24.               ! 21th june     in hours
      ztimemin = ztimemax + rjjhh * zyydd / 2      ! 21th december in hours
      !                                            ! NB: rjjhh * zyydd / 4 = 1 seasonal cycle in hours

      ! mean intensity at 0.105 ; srqt(2) because projected with 45deg angle
      ztau = 0.105 / SQRT( 2. )
      ! seasonal oscillation intensity
      ztau_sais = 0.015
      ztaun = ztau - ztau_sais * COS( (ztime - ztimemax) / (ztimemin - ztimemax) * rpi )
      DO jj = 1, jpj
         DO ji = 1, jpi
           ! domain from 15deg to 50deg and 1/2 period along 14deg
           ! so 5/4 of half period with seasonal cycle
           utau(ji,jj) = - ztaun * SIN( rpi * (gphiu(ji,jj) - 15.) / (29.-15.) )
           vtau(ji,jj) =   ztaun * SIN( rpi * (gphiv(ji,jj) - 15.) / (29.-15.) )
         END DO
      END DO

      ! module of wind stress and wind speed at T-point
      zcoef = 1. / ( zrhoa * zcdrag ) 
!CDIR NOVERRCHK
      DO jj = 2, jpjm1
!CDIR NOVERRCHK
         DO ji = 2, jpim1   ! vect. opt.
            ztx = utau(ji-1,jj  ) + utau(ji,jj) 
            zty = vtau(ji  ,jj-1) + vtau(ji,jj) 
            zmod = 0.5 * SQRT( ztx * ztx + zty * zty )
            taum(ji,jj) = zmod
            wndm(ji,jj) = SQRT( zmod * zcoef )
         END DO
      END DO
      CALL lbc_lnk( taum(:,:), 'T', 1. )   ;   CALL lbc_lnk( wndm(:,:), 'T', 1. )

      ! ---------------------------------- !
      !  control print at first time-step  !
      ! ---------------------------------- !
      IF( kt == nit000 .AND. lwp ) THEN 
         WRITE(numout,*)
         WRITE(numout,*)'sbc_gyre : analytical surface fluxes for GYRE configuration'               
         WRITE(numout,*)'~~~~~~~~ ' 
         WRITE(numout,*)'           nyear      = ', nyear
         WRITE(numout,*)'           nmonth     = ', nmonth
         WRITE(numout,*)'           nday       = ', nday
         WRITE(numout,*)'           nday_year  = ', nday_year
         WRITE(numout,*)'           ztime      = ', ztime
         WRITE(numout,*)'           ztimemax   = ', ztimemax
         WRITE(numout,*)'           ztimemin   = ', ztimemin
         WRITE(numout,*)'           ztimemax1  = ', ztimemax1
         WRITE(numout,*)'           ztimemin1  = ', ztimemin1
         WRITE(numout,*)'           ztimemax2  = ', ztimemax2
         WRITE(numout,*)'           ztimemin2  = ', ztimemin2
         WRITE(numout,*)'           zyear0     = ', zyear0
         WRITE(numout,*)'           zmonth0    = ', zmonth0
         WRITE(numout,*)'           zday0      = ', zday0
         WRITE(numout,*)'           zday_year0 = ', zday_year0
         WRITE(numout,*)'           zyydd      = ', zyydd
         WRITE(numout,*)'           zemp_S     = ', zemp_S
         WRITE(numout,*)'           zemp_N     = ', zemp_N
         WRITE(numout,*)'           zemp_sais  = ', zemp_sais
         WRITE(numout,*)'           zTstar     = ', zTstar
         WRITE(numout,*)'           zsumemp    = ', zsumemp
         WRITE(numout,*)'           zsurf      = ', zsurf
         WRITE(numout,*)'           ztrp       = ', ztrp
         WRITE(numout,*)'           zconv      = ', zconv
         WRITE(numout,*)'           ndastp     = ', ndastp
         WRITE(numout,*)'           adatrj     = ', adatrj
      ENDIF
      !
   END SUBROUTINE sbc_gyre

   !!======================================================================
END MODULE sbcana
