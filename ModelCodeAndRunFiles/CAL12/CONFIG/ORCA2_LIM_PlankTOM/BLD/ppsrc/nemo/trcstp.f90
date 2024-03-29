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

MODULE trcstp
   !!======================================================================
   !!                       ***  MODULE trcstp  ***
   !! Time-stepping    : time loop of opa for passive tracer
   !!======================================================================
   !! History :  1.0  !  2004-03  (C. Ethe)  Original
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   trc_stp      : passive tracer system time-stepping
   !!----------------------------------------------------------------------
   USE oce_trc          ! ocean dynamics and active tracers variables
   USE sbc_oce
   USE trc
   USE trctrp           ! passive tracers transport
   USE trcsms           ! passive tracers sources and sinks
   USE prtctl_trc       ! Print control for debbuging
   USE trcdia
   USE trcwri
   USE trcrst
   USE trdtrc_oce
   USE trdmxl_trc
   USE iom
   USE in_out_manager
   USE trcsub

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_stp    ! called by step

   REAL(wp), DIMENSION(:,:,:), SAVE, ALLOCATABLE ::   qsr_arr ! save qsr during TOP time-step
   REAL(wp) :: rdt_sampl
   INTEGER  :: nb_rec_per_days
   INTEGER  :: isecfst, iseclast
   LOGICAL  :: llnew

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
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcstp.F90 5407 2015-06-11 19:13:22Z smasson $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_stp( kt )
      !!-------------------------------------------------------------------
      !!                     ***  ROUTINE trc_stp  ***
      !!                      
      !! ** Purpose : Time loop of opa for passive tracer
      !! 
      !! ** Method  : 
      !!              Compute the passive tracers trends 
      !!              Update the passive tracers
      !!-------------------------------------------------------------------
      INTEGER, INTENT( in ) ::  kt      ! ocean time-step index
      INTEGER               ::  jk, jn  ! dummy loop indices
      REAL(wp)              ::  ztrai
      CHARACTER (len=25)    ::  charout 

      !!-------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('trc_stp')
      !
      IF( kt == nittrc000 .AND. lk_trdmxl_trc )  CALL trd_mxl_trc_init    ! trends: Mixed-layer
      !
      IF( lk_vvl ) THEN                                                   ! update ocean volume due to ssh temporal evolution
         DO jk = 1, jpk
            cvol(:,:,jk) = e1e2t(:,:) * e3t_0(:,:,jk) * tmask(:,:,jk)
         END DO
         IF( lk_degrad )  cvol(:,:,:) = cvol(:,:,:) * facvol(:,:,:)       ! degrad option: reduction by facvol
         areatot         = glob_sum( cvol(:,:,:) )
      ENDIF
      !
      IF( l_trcdm2dc )   CALL trc_mean_qsr( kt )
      !    
      IF( nn_dttrc /= 1 )   CALL trc_sub_stp( kt )  ! averaging physical variables for sub-stepping
      !    
      IF( MOD( kt , nn_dttrc ) == 0 ) THEN      ! only every nn_dttrc time step
         !
         IF(ln_ctl) THEN
            WRITE(charout,FMT="('kt =', I4,'  d/m/y =',I2,I2,I4)") kt, nday, nmonth, nyear
            CALL prt_ctl_trc_info(charout)
         ENDIF
         !
         tra(:,:,:,:) = 0.e0
         !
                                   CALL trc_rst_opn  ( kt )       ! Open tracer restart file 
         IF( lrst_trc )            CALL trc_rst_cal  ( kt, 'WRITE' )   ! calendar
         IF( lk_iomput ) THEN  ;   CALL trc_wri      ( kt )       ! output of passive tracers with iom I/O manager
         ELSE                  ;   CALL trc_dia      ( kt )       ! output of passive tracers with old I/O manager
         ENDIF
                                   CALL trc_sms      ( kt )       ! tracers: sinks and sources
                                   CALL trc_trp      ( kt )       ! transport of passive tracers
         IF( kt == nittrc000 ) THEN
            CALL iom_close( numrtr )       ! close input tracer restart file
            IF(lwm) CALL FLUSH( numont )   ! flush namelist output
         ENDIF
         IF( lrst_trc )            CALL trc_rst_wri  ( kt )       ! write tracer restart file
         IF( lk_trdmxl_trc  )      CALL trd_mxl_trc  ( kt )       ! trends: Mixed-layer
         !
         IF( nn_dttrc /= 1   )     CALL trc_sub_reset( kt )       ! resetting physical variables when sub-stepping
         !
      ENDIF
      !
      ztrai = 0._wp                                                   !  content of all tracers
      DO jn = 1, jptra
         ztrai = ztrai + glob_sum( trn(:,:,:,jn) * cvol(:,:,:)   )
      END DO
      IF( lwp ) WRITE(numstr,9300) kt,  ztrai / areatot
9300  FORMAT(i10,e18.10)
      !
      IF( nn_timing == 1 )   CALL timing_stop('trc_stp')
      !
   END SUBROUTINE trc_stp

   SUBROUTINE trc_mean_qsr( kt )
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE trc_mean_qsr  ***
      !!
      !! ** Purpose :  Compute daily mean qsr for biogeochemical model in case
      !!               of diurnal cycle
      !!
      !! ** Method  : store in TOP the qsr every hour ( or every time-step the latter 
      !!              is greater than 1 hour ) and then, compute the  mean with 
      !!              a moving average over 24 hours. 
      !!              In coupled mode, the sampling is done at every coupling frequency 
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt
      INTEGER  :: jn

      IF( kt == nittrc000 ) THEN
         IF( ln_cpl )  THEN  
            rdt_sampl = 86400. / ncpl_qsr_freq
            nb_rec_per_days = ncpl_qsr_freq
         ELSE  
            rdt_sampl = MAX( 3600., rdt * nn_dttrc )
            nb_rec_per_days = INT( 86400 / rdt_sampl )
         ENDIF
         !
         IF( lwp ) THEN
            WRITE(numout,*) 
            WRITE(numout,*) ' Sampling frequency dt = ', rdt_sampl, 's','   Number of sampling per day  nrec = ', nb_rec_per_days
            WRITE(numout,*) 
         ENDIF
         !
         ALLOCATE( qsr_arr(jpi,jpj,nb_rec_per_days ) )
         DO jn = 1, nb_rec_per_days
            qsr_arr(:,:,jn) = qsr(:,:)
         ENDDO
         qsr_mean(:,:) = qsr(:,:)
         !
         isecfst  = nsec_year + nsec1jan000   !   number of seconds between Jan. 1st 00h of nit000 year and the middle of time step
         iseclast = isecfst
         !
      ENDIF
      !
      iseclast = nsec_year + nsec1jan000
      llnew   = ( iseclast - isecfst )  > INT( rdt_sampl )   !   new shortwave to store
      IF( kt /= nittrc000 .AND. llnew ) THEN
          IF( lwp ) WRITE(numout,*) ' New shortwave to sample for TOP at time kt = ', kt, &
             &                      ' time = ', (iseclast+rdt*nn_dttrc/2.)/3600.,'hours '
          isecfst = iseclast
          DO jn = 1, nb_rec_per_days - 1
             qsr_arr(:,:,jn) = qsr_arr(:,:,jn+1)
          ENDDO
          qsr_arr (:,:,nb_rec_per_days) = qsr(:,:)
          qsr_mean(:,:                ) = SUM( qsr_arr(:,:,:), 3 ) / nb_rec_per_days
      ENDIF
      !
   END SUBROUTINE trc_mean_qsr


   !!======================================================================
END MODULE trcstp
