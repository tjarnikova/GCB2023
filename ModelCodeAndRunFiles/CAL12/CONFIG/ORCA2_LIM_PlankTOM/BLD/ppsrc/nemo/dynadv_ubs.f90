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

MODULE dynadv_ubs
   !!======================================================================
   !!                       ***  MODULE  dynadv_ubs  ***
   !! Ocean dynamics: Update the momentum trend with the flux form advection
   !!                 trend using a 3rd order upstream biased scheme
   !!======================================================================
   !! History :  2.0  ! 2006-08  (R. Benshila, L. Debreu)  Original code
   !!            3.2  ! 2009-07  (R. Benshila)  Suppression of rigid-lid option
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_adv_ubs   : flux form momentum advection using    (ln_dynadv=T)
   !!                   an 3rd order Upstream Biased Scheme or Quick scheme
   !!                   combined with 2nd or 4th order finite differences 
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics
   !
   USE in_out_manager ! I/O manager
   USE prtctl         ! Print control
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! MPP library
   USE wrk_nemo       ! Memory Allocation
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   REAL(wp), PARAMETER :: gamma1 = 1._wp/3._wp  ! =1/4 quick      ; =1/3  3rd order UBS
   REAL(wp), PARAMETER :: gamma2 = 1._wp/32._wp ! =0   2nd order  ; =1/32 4th order centred

   PUBLIC   dyn_adv_ubs   ! routine called by step.F90

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
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: dynadv_ubs.F90 5069 2015-02-09 10:08:53Z smasson $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_adv_ubs( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_adv_ubs  ***
      !!
      !! ** Purpose :   Compute the now momentum advection trend in flux form
      !!              and the general trend of the momentum equation.
      !!
      !! ** Method  :   The scheme is the one implemeted in ROMS. It depends 
      !!      on two parameter gamma1 and gamma2. The former control the 
      !!      upstream baised part of the scheme and the later the centred 
      !!      part:     gamma1 = 0    pure centered  (no diffusive part)
      !!                       = 1/4  Quick scheme
      !!                       = 1/3  3rd order Upstream biased scheme
      !!                gamma2 = 0    2nd order finite differencing 
      !!                       = 1/32 4th order finite differencing
      !!      For stability reasons, the first term of the fluxes which cor-
      !!      responds to a second order centered scheme is evaluated using  
      !!      the now velocity (centered in time) while the second term which  
      !!      is the diffusive part of the scheme, is evaluated using the 
      !!      before velocity (forward in time). 
      !!      Default value (hard coded in the begining of the module) are 
      !!      gamma1=1/3 and gamma2=1/32.
      !!
      !! ** Action : - (ua,va) updated with the 3D advective momentum trends
      !!
      !! Reference : Shchepetkin & McWilliams, 2005, Ocean Modelling. 
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk            ! dummy loop indices
      REAL(wp) ::   zbu, zbv    ! temporary scalars
      REAL(wp) ::   zui, zvj, zfuj, zfvi, zl_u, zl_v   ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:,:  ) ::  zfu, zfv
      REAL(wp), POINTER, DIMENSION(:,:,:  ) ::  zfu_t, zfv_t, zfu_f, zfv_f, zfu_uw, zfv_vw, zfw
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::  zlu_uu, zlv_vv, zlu_uv, zlv_vu
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_adv_ubs')
      !
      CALL wrk_alloc( jpi, jpj, jpk,       zfu_t , zfv_t , zfu_f , zfv_f, zfu_uw, zfv_vw, zfu, zfv, zfw )
      CALL wrk_alloc( jpi, jpj, jpk, jpts, zlu_uu, zlv_vv, zlu_uv, zlv_vu                               )
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_adv_ubs : UBS flux form momentum advection'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      zfu_t(:,:,:) = 0._wp
      zfv_t(:,:,:) = 0._wp
      zfu_f(:,:,:) = 0._wp
      zfv_f(:,:,:) = 0._wp
      !
      zlu_uu(:,:,:,:) = 0._wp
      zlv_vv(:,:,:,:) = 0._wp 
      zlu_uv(:,:,:,:) = 0._wp 
      zlv_vu(:,:,:,:) = 0._wp 

      IF( l_trddyn ) THEN           ! Save ua and va trends
         zfu_uw(:,:,:) = ua(:,:,:)
         zfv_vw(:,:,:) = va(:,:,:)
      ENDIF

      !                                      ! =========================== !
      DO jk = 1, jpkm1                       !  Laplacian of the velocity  !
         !                                   ! =========================== !
         !                                         ! horizontal volume fluxes
         zfu(:,:,jk) = e2u(:,:) * e3u_0(:,:,jk) * un(:,:,jk)
         zfv(:,:,jk) = e1v(:,:) * e3v_0(:,:,jk) * vn(:,:,jk)
         !            
         DO jj = 2, jpjm1                          ! laplacian
            DO ji = 2, jpim1   ! vector opt.
               !
               zlu_uu(ji,jj,jk,1) = ( ub (ji+1,jj  ,jk) - 2.*ub (ji,jj,jk) + ub (ji-1,jj  ,jk) ) * umask(ji,jj,jk)
               zlv_vv(ji,jj,jk,1) = ( vb (ji  ,jj+1,jk) - 2.*vb (ji,jj,jk) + vb (ji  ,jj-1,jk) ) * vmask(ji,jj,jk)
               zlu_uv(ji,jj,jk,1) = ( ub (ji  ,jj+1,jk) - ub (ji  ,jj  ,jk) ) * fmask(ji  ,jj  ,jk)   &
                  &               - ( ub (ji  ,jj  ,jk) - ub (ji  ,jj-1,jk) ) * fmask(ji  ,jj-1,jk)
               zlv_vu(ji,jj,jk,1) = ( vb (ji+1,jj  ,jk) - vb (ji  ,jj  ,jk) ) * fmask(ji  ,jj  ,jk)   &
                  &               - ( vb (ji  ,jj  ,jk) - vb (ji-1,jj  ,jk) ) * fmask(ji-1,jj  ,jk)
               !
               zlu_uu(ji,jj,jk,2) = ( zfu(ji+1,jj  ,jk) - 2.*zfu(ji,jj,jk) + zfu(ji-1,jj  ,jk) ) * umask(ji,jj,jk)
               zlv_vv(ji,jj,jk,2) = ( zfv(ji  ,jj+1,jk) - 2.*zfv(ji,jj,jk) + zfv(ji  ,jj-1,jk) ) * vmask(ji,jj,jk)
               zlu_uv(ji,jj,jk,2) = ( zfu(ji  ,jj+1,jk) - zfu(ji  ,jj  ,jk) ) * fmask(ji  ,jj  ,jk)   &
                  &               - ( zfu(ji  ,jj  ,jk) - zfu(ji  ,jj-1,jk) ) * fmask(ji  ,jj-1,jk)
               zlv_vu(ji,jj,jk,2) = ( zfv(ji+1,jj  ,jk) - zfv(ji  ,jj  ,jk) ) * fmask(ji  ,jj  ,jk)   &
                  &               - ( zfv(ji  ,jj  ,jk) - zfv(ji-1,jj  ,jk) ) * fmask(ji-1,jj  ,jk)
            END DO
         END DO
      END DO
      CALL lbc_lnk( zlu_uu(:,:,:,1), 'U', 1. )   ;   CALL lbc_lnk( zlu_uv(:,:,:,1), 'U', 1. )
      CALL lbc_lnk( zlu_uu(:,:,:,2), 'U', 1. )   ;   CALL lbc_lnk( zlu_uv(:,:,:,2), 'U', 1. )
      CALL lbc_lnk( zlv_vv(:,:,:,1), 'V', 1. )   ;   CALL lbc_lnk( zlv_vu(:,:,:,1), 'V', 1. )
      CALL lbc_lnk( zlv_vv(:,:,:,2), 'V', 1. )   ;   CALL lbc_lnk( zlv_vu(:,:,:,2), 'V', 1. ) 
      
      !                                      ! ====================== !
      !                                      !  Horizontal advection  !
      DO jk = 1, jpkm1                       ! ====================== !
         !                                         ! horizontal volume fluxes
         zfu(:,:,jk) = 0.25 * e2u(:,:) * e3u_0(:,:,jk) * un(:,:,jk)
         zfv(:,:,jk) = 0.25 * e1v(:,:) * e3v_0(:,:,jk) * vn(:,:,jk)
         !
         DO jj = 1, jpjm1                          ! horizontal momentum fluxes at T- and F-point
            DO ji = 1, jpim1   ! vector opt.
               zui = ( un(ji,jj,jk) + un(ji+1,jj  ,jk) )
               zvj = ( vn(ji,jj,jk) + vn(ji  ,jj+1,jk) )
               !
               IF (zui > 0) THEN   ;   zl_u = zlu_uu(ji  ,jj,jk,1)
               ELSE                ;   zl_u = zlu_uu(ji+1,jj,jk,1)
               ENDIF
               IF (zvj > 0) THEN   ;   zl_v = zlv_vv(ji,jj  ,jk,1)
               ELSE                ;   zl_v = zlv_vv(ji,jj+1,jk,1)
               ENDIF
               !
               zfu_t(ji+1,jj  ,jk) = ( zfu(ji,jj,jk) + zfu(ji+1,jj  ,jk)                               &
                  &                    - gamma2 * ( zlu_uu(ji,jj,jk,2) + zlu_uu(ji+1,jj  ,jk,2) )  )   &
                  &                * ( zui - gamma1 * zl_u)
               zfv_t(ji  ,jj+1,jk) = ( zfv(ji,jj,jk) + zfv(ji  ,jj+1,jk)                               &
                  &                    - gamma2 * ( zlv_vv(ji,jj,jk,2) + zlv_vv(ji  ,jj+1,jk,2) )  )   &
                  &                * ( zvj - gamma1 * zl_v)
               !
               zfuj = ( zfu(ji,jj,jk) + zfu(ji  ,jj+1,jk) )
               zfvi = ( zfv(ji,jj,jk) + zfv(ji+1,jj  ,jk) )
               IF (zfuj > 0) THEN   ;    zl_v = zlv_vu( ji  ,jj  ,jk,1)
               ELSE                 ;    zl_v = zlv_vu( ji+1,jj,jk,1)
               ENDIF
               IF (zfvi > 0) THEN   ;    zl_u = zlu_uv( ji,jj  ,jk,1)
               ELSE                 ;    zl_u = zlu_uv( ji,jj+1,jk,1)
               ENDIF
               !
               zfv_f(ji  ,jj  ,jk) = ( zfvi - gamma2 * ( zlv_vu(ji,jj,jk,2) + zlv_vu(ji+1,jj  ,jk,2) )  )   &
                  &                * ( un(ji,jj,jk) + un(ji  ,jj+1,jk) - gamma1 * zl_u )
               zfu_f(ji  ,jj  ,jk) = ( zfuj - gamma2 * ( zlu_uv(ji,jj,jk,2) + zlu_uv(ji  ,jj+1,jk,2) )  )   &
                  &                * ( vn(ji,jj,jk) + vn(ji+1,jj  ,jk) - gamma1 * zl_v )
            END DO
         END DO
         DO jj = 2, jpjm1                          ! divergence of horizontal momentum fluxes
            DO ji = 2, jpim1   ! vector opt.
               zbu = e1u(ji,jj) * e2u(ji,jj) * e3u_0(ji,jj,jk)
               zbv = e1v(ji,jj) * e2v(ji,jj) * e3v_0(ji,jj,jk)
               !
               ua(ji,jj,jk) = ua(ji,jj,jk) - (  zfu_t(ji+1,jj  ,jk) - zfu_t(ji  ,jj  ,jk)    &
                  &                           + zfv_f(ji  ,jj  ,jk) - zfv_f(ji  ,jj-1,jk)  ) / zbu
               va(ji,jj,jk) = va(ji,jj,jk) - (  zfu_f(ji  ,jj  ,jk) - zfu_f(ji-1,jj  ,jk)    &
                  &                           + zfv_t(ji  ,jj+1,jk) - zfv_t(ji  ,jj  ,jk)  ) / zbv
            END DO
         END DO
      END DO
      IF( l_trddyn ) THEN                          ! save the horizontal advection trend for diagnostic
         zfu_uw(:,:,:) = ua(:,:,:) - zfu_uw(:,:,:)
         zfv_vw(:,:,:) = va(:,:,:) - zfv_vw(:,:,:)
         CALL trd_dyn( zfu_uw, zfv_vw, jpdyn_keg, kt )
         zfu_t(:,:,:) = ua(:,:,:)
         zfv_t(:,:,:) = va(:,:,:)
      ENDIF

      !                                      ! ==================== !
      !                                      !  Vertical advection  !
      DO jk = 1, jpkm1                       ! ==================== !
         !                                         ! Vertical volume fluxes�
         zfw(:,:,jk) = 0.25 * e1t(:,:) * e2t(:,:) * wn(:,:,jk)
         !
         IF( jk == 1 ) THEN                        ! surface/bottom advective fluxes                   
            zfu_uw(:,:,jpk) = 0.e0                      ! Bottom  value : flux set to zero
            zfv_vw(:,:,jpk) = 0.e0
            !                                           ! Surface value :
            IF( lk_vvl ) THEN                                ! variable volume : flux set to zero
               zfu_uw(:,:, 1 ) = 0.e0    
               zfv_vw(:,:, 1 ) = 0.e0
            ELSE                                             ! constant volume : advection through the surface
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1
                     zfu_uw(ji,jj, 1 ) = 2.e0 * ( zfw(ji,jj,1) + zfw(ji+1,jj  ,1) ) * un(ji,jj,1)
                     zfv_vw(ji,jj, 1 ) = 2.e0 * ( zfw(ji,jj,1) + zfw(ji  ,jj+1,1) ) * vn(ji,jj,1)
                  END DO
               END DO
            ENDIF
         ELSE                                      ! interior fluxes
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zfu_uw(ji,jj,jk) = ( zfw(ji,jj,jk)+ zfw(ji+1,jj  ,jk) ) * ( un(ji,jj,jk) + un(ji,jj,jk-1) )
                  zfv_vw(ji,jj,jk) = ( zfw(ji,jj,jk)+ zfw(ji  ,jj+1,jk) ) * ( vn(ji,jj,jk) + vn(ji,jj,jk-1) )
               END DO
            END DO
         ENDIF
      END DO
      DO jk = 1, jpkm1                             ! divergence of vertical momentum flux divergence
         DO jj = 2, jpjm1 
            DO ji = 2, jpim1   ! vector opt.
               ua(ji,jj,jk) =  ua(ji,jj,jk) - ( zfu_uw(ji,jj,jk) - zfu_uw(ji,jj,jk+1) )    &
                  &  / ( e1u(ji,jj) * e2u(ji,jj) * e3u_0(ji,jj,jk) )
               va(ji,jj,jk) =  va(ji,jj,jk) - ( zfv_vw(ji,jj,jk) - zfv_vw(ji,jj,jk+1) )    &
                  &  / ( e1v(ji,jj) * e2v(ji,jj) * e3v_0(ji,jj,jk) )
            END DO
         END DO
      END DO
      !
      IF( l_trddyn ) THEN                          ! save the vertical advection trend for diagnostic
         zfu_t(:,:,:) = ua(:,:,:) - zfu_t(:,:,:)
         zfv_t(:,:,:) = va(:,:,:) - zfv_t(:,:,:)
         CALL trd_dyn( zfu_t, zfv_t, jpdyn_zad, kt )
      ENDIF
      !                                            ! Control print
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' ubs2 adv - Ua: ', mask1=umask,   &
         &                       tab3d_2=va, clinfo2=           ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      CALL wrk_dealloc( jpi, jpj, jpk,       zfu_t , zfv_t , zfu_f , zfv_f, zfu_uw, zfv_vw, zfu, zfv, zfw )
      CALL wrk_dealloc( jpi, jpj, jpk, jpts, zlu_uu, zlv_vv, zlu_uv, zlv_vu                               )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_adv_ubs')
      !
   END SUBROUTINE dyn_adv_ubs

   !!==============================================================================
END MODULE dynadv_ubs
