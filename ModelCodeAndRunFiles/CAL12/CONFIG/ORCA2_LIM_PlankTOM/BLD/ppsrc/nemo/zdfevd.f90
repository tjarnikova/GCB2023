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

MODULE zdfevd
   !!======================================================================
   !!                       ***  MODULE  zdfevd  ***
   !! Ocean physics: parameterization of convection through an enhancement
   !!                of vertical eddy mixing coefficient
   !!======================================================================
   !! History :  OPA  !  1997-06  (G. Madec, A. Lazar)  Original code
   !!   NEMO     1.0  !  2002-06  (G. Madec)  F90: Free form and module
   !!             -   !  2005-06  (C. Ethe) KPP parameterization
   !!            3.2  !  2009-03  (M. Leclair, G. Madec, R. Benshila) test on both before & after
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   zdf_evd      : increase the momentum and tracer Kz at the location of
   !!                  statically unstable portion of the water column (ln_zdfevd=T)
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE zdf_oce         ! ocean vertical physics variables
   USE zdfkpp          ! KPP vertical mixing
   USE in_out_manager  ! I/O manager
   USE iom             ! for iom_put
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   zdf_evd    ! called by step.F90

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
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: zdfevd.F90 4990 2014-12-15 16:42:49Z timgraham $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE zdf_evd( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_evd  ***
      !!                   
      !! ** Purpose :   Local increased the vertical eddy viscosity and diffu-
      !!      sivity coefficients when a static instability is encountered.
      !!
      !! ** Method  :   avt, avm, and the 4 neighbouring avmu, avmv coefficients
      !!      are set to avevd (namelist parameter) if the water column is 
      !!      statically unstable (i.e. if rn2 < -1.e-12 )
      !!
      !! ** Action  :   avt, avm, avmu, avmv updted in static instability cases
      !!
      !! References :   Lazar, A., these de l'universite Paris VI, France, 1997
      !!----------------------------------------------------------------------
      USE oce,   zavt_evd => ua , zavm_evd => va  ! (ua,va) used ua workspace
      !
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step indexocean time step
      !
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('zdf_evd')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'zdf_evd : Enhanced Vertical Diffusion (evd)'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
         IF(lwp) WRITE(numout,*)
      ENDIF

      zavt_evd(:,:,:) = avt(:,:,:)           ! set avt prior to evd application

      SELECT CASE ( nn_evdm )
      !
      CASE ( 1 )           ! enhance vertical eddy viscosity and diffusivity (if rn2<-1.e-12)
         !
         zavm_evd(:,:,:) = avm(:,:,:)           ! set avm prior to evd application
         !
         DO jk = 1, jpkm1 
            DO jj = 2, jpj             ! no vector opt.
               DO ji = 2, jpi
                  IF(  MIN( rn2(ji,jj,jk), rn2b(ji,jj,jk) ) <= -1.e-12 ) THEN
                     avt (ji  ,jj  ,jk) = rn_avevd * tmask(ji  ,jj  ,jk)
                     avm (ji  ,jj  ,jk) = rn_avevd * tmask(ji  ,jj  ,jk)
                     avmu(ji  ,jj  ,jk) = rn_avevd * umask(ji  ,jj  ,jk)
                     avmu(ji-1,jj  ,jk) = rn_avevd * umask(ji-1,jj  ,jk)
                     avmv(ji  ,jj  ,jk) = rn_avevd * vmask(ji  ,jj  ,jk)
                     avmv(ji  ,jj-1,jk) = rn_avevd * vmask(ji  ,jj-1,jk)
                  ENDIF
               END DO
            END DO
         END DO 
         CALL lbc_lnk( avt , 'W', 1. )   ;   CALL lbc_lnk( avm , 'W', 1. )   ! Lateral boundary conditions
         CALL lbc_lnk( avmu, 'U', 1. )   ;   CALL lbc_lnk( avmv, 'V', 1. )
         !
         zavm_evd(:,:,:) = avm(:,:,:) - zavm_evd(:,:,:)   ! change in avm due to evd
         CALL iom_put( "avm_evd", zavm_evd )              ! output this change
         !
      CASE DEFAULT         ! enhance vertical eddy diffusivity only (if rn2<-1.e-12) 
         DO jk = 1, jpkm1
!!!         WHERE( rn2(:,:,jk) <= -1.e-12 ) avt(:,:,jk) = tmask(:,:,jk) * avevd   ! agissant sur T SEUL! 
            DO jj = 1, jpj             ! loop over the whole domain (no lbc_lnk call)
               DO ji = 1, jpi
                  IF(  MIN( rn2(ji,jj,jk), rn2b(ji,jj,jk) ) <= -1.e-12 )   &
                     avt(ji,jj,jk) = rn_avevd * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
         !
      END SELECT 

      zavt_evd(:,:,:) = avt(:,:,:) - zavt_evd(:,:,:)   ! change in avt due to evd
      CALL iom_put( "avt_evd", zavt_evd )              ! output this change
      !
      IF( nn_timing == 1 )  CALL timing_stop('zdf_evd')
      !
   END SUBROUTINE zdf_evd

   !!======================================================================
END MODULE zdfevd
