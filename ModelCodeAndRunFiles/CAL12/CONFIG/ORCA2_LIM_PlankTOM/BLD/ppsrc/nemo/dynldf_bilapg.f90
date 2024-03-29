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

MODULE dynldf_bilapg
   !!======================================================================
   !!                       ***  MODULE  dynldf_bilapg  ***
   !! Ocean dynamics:  lateral viscosity trend
   !!======================================================================
   !! History :  OPA  !  1997-07  (G. Madec)  Original code
   !!  NEMO      1.0  !  2002-08  (G. Madec)  F90: Free form and module
   !!            2.0  !  2004-08  (C. Talandier) New trends organization
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_ldfslp'                              Rotation of mixing tensor
   !!----------------------------------------------------------------------
   !!   dyn_ldf_bilapg : update the momentum trend with the horizontal part
   !!                    of the horizontal s-coord. bilaplacian diffusion
   !!   ldfguv         : 
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE ldfdyn_oce      ! ocean dynamics lateral physics
   USE zdf_oce         ! ocean vertical physics
   USE ldfslp          ! iso-neutral slopes available
   USE ldftra_oce, ONLY: ln_traldf_iso
   !
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_ldf_bilapg       ! called by step.F90

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  zfuw, zfvw , zdiu, zdiv   ! 2D workspace (ldfguv)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::  zdju, zdj1u, zdjv, zdj1v  ! 2D workspace (ldfguv)

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
   !!                    ***  ldfdyn_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsahm., the lateral eddy viscosity coeff. 
   !!      with a constant, or 1D, or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldfdyn_substitute.h90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   !!
   !! fsahmt, fsahmf - used for laplaian operator only
   !! fsahmu, fsahmv - used for bilaplacian operator only
   !!
!   ' key_dynldf_c3d' :                  3D coefficient
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: dynldf_bilapg.F90 4990 2014-12-15 16:42:49Z timgraham $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION dyn_ldf_bilapg_alloc()
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE dyn_ldf_bilapg_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( zfuw(jpi,jpk) , zfvw (jpi,jpk) , zdiu(jpi,jpk) , zdiv (jpi,jpk) ,     &
         &      zdju(jpi,jpk) , zdj1u(jpi,jpk) , zdjv(jpi,jpk) , zdj1v(jpi,jpk) , STAT=dyn_ldf_bilapg_alloc )
         !
      IF( dyn_ldf_bilapg_alloc /= 0 )   CALL ctl_warn('dyn_ldf_bilapg_alloc: failed to allocate arrays')
   END FUNCTION dyn_ldf_bilapg_alloc


   SUBROUTINE dyn_ldf_bilapg( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dyn_ldf_bilapg  ***
      !!                      
      !! ** Purpose :   Compute the before trend of the horizontal momentum
      !!      diffusion and add it to the general trend of momentum equation.
      !!
      !! ** Method  :   The lateral momentum diffusive trends is provided by a 
      !!      a 4th order operator rotated along geopotential surfaces. It is 
      !!      computed using before fields (forward in time) and geopotential
      !!      slopes computed in routine inildf.
      !!         -1- compute the geopotential harmonic operator applied to
      !!      (ub,vb) and multiply it by the eddy diffusivity coefficient
      !!      (done by a call to ldfgpu and ldfgpv routines) The result is in
      !!      (zwk1,zwk2) arrays. Applied the domain lateral boundary conditions
      !!      by call to lbc_lnk.
      !!         -2- applied to (zwk1,zwk2) the geopotential harmonic operator
      !!      by a second call to ldfgpu and ldfgpv routines respectively. The
      !!      result is in (zwk3,zwk4) arrays.
      !!         -3- Add this trend to the general trend (ta,sa):
      !!            (ua,va) = (ua,va) + (zwk3,zwk4)
      !!
      !! ** Action  : - Update (ua,va) arrays with the before geopotential
      !!                biharmonic mixing trend.
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
      !
      INTEGER ::   ji, jj, jk                 ! dummy loop indices
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zwk1, zwk2, zwk3, zwk4
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_ldf_bilapg')
      !
      CALL wrk_alloc( jpi, jpj, jpk, zwk1, zwk2, zwk3, zwk4 ) 
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_ldf_bilapg : horizontal biharmonic operator in s-coordinate'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'
         !                                      ! allocate dyn_ldf_bilapg arrays
         IF( dyn_ldf_bilapg_alloc() /= 0 )   CALL ctl_stop('STOP', 'dyn_ldf_bilapg: failed to allocate arrays')
      ENDIF

      ! s-coordinate: Iso-level diffusion on tracer, but geopotential level diffusion on momentum
      IF( ln_dynldf_hor .AND. ln_traldf_iso ) THEN
         !
         DO jk = 1, jpk         ! set the slopes of iso-level
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  uslp (ji,jj,jk) = -1./e1u(ji,jj) * ( gdept_0(ji+1,jj,jk) - gdept_0(ji ,jj ,jk) ) * umask(ji,jj,jk)
                  vslp (ji,jj,jk) = -1./e2v(ji,jj) * ( gdept_0(ji,jj+1,jk) - gdept_0(ji ,jj ,jk) ) * vmask(ji,jj,jk)
                  wslpi(ji,jj,jk) = -1./e1t(ji,jj) * ( gdepw_0(ji+1,jj,jk) - gdepw_0(ji-1,jj,jk) ) * tmask(ji,jj,jk) * 0.5
                  wslpj(ji,jj,jk) = -1./e2t(ji,jj) * ( gdepw_0(ji,jj+1,jk) - gdepw_0(ji,jj-1,jk) ) * tmask(ji,jj,jk) * 0.5
               END DO
            END DO
         END DO
         ! Lateral boundary conditions on the slopes
         CALL lbc_lnk( uslp , 'U', -1. )      ;      CALL lbc_lnk( vslp , 'V', -1. )
         CALL lbc_lnk( wslpi, 'W', -1. )      ;      CALL lbc_lnk( wslpj, 'W', -1. )
 
!!bug
         IF( kt == nit000 ) then
            IF(lwp) WRITE(numout,*) ' max slop: u', SQRT( MAXVAL(uslp*uslp)), ' v ', SQRT(MAXVAL(vslp)),  &
               &                             ' wi', sqrt(MAXVAL(wslpi))     , ' wj', sqrt(MAXVAL(wslpj))
         endif
!!end
      ENDIF

      zwk1(:,:,:) = 0.e0   ;   zwk3(:,:,:) = 0.e0
      zwk2(:,:,:) = 0.e0   ;   zwk4(:,:,:) = 0.e0

      ! Laplacian of (ub,vb) multiplied by ahm
      ! --------------------------------------  
      CALL ldfguv( ub, vb, zwk1, zwk2, 1 )      ! rotated harmonic operator applied to (ub,vb)
      !                                         ! and multiply by ahmu, ahmv (output in (zwk1,zwk2) )
      CALL lbc_lnk( zwk1, 'U', -1. )   ;   CALL lbc_lnk( zwk2, 'V', -1. )     ! Lateral boundary conditions

      ! Bilaplacian of (ub,vb)
      ! ---------------------- 
      CALL ldfguv( zwk1, zwk2, zwk3, zwk4, 2 )  ! rotated harmonic operator applied to (zwk1,zwk2) 
      !                                         ! (output in (zwk3,zwk4) )

      ! Update the momentum trends
      ! --------------------------
      DO jj = 2, jpjm1               ! add the diffusive trend to the general momentum trends
         DO jk = 1, jpkm1
            DO ji = 2, jpim1
               ua(ji,jj,jk) = ua(ji,jj,jk) + zwk3(ji,jj,jk)
               va(ji,jj,jk) = va(ji,jj,jk) + zwk4(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      CALL wrk_dealloc( jpi, jpj, jpk, zwk1, zwk2, zwk3, zwk4 ) 
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_ldf_bilapg')
      !
   END SUBROUTINE dyn_ldf_bilapg


   SUBROUTINE ldfguv( pu, pv, plu, plv, kahm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldfguv  ***
      !!      
      !! ** Purpose :   Apply a geopotential harmonic operator to (pu,pv)
      !!      (defined at u- and v-points) and multiply it by the eddy
      !!      viscosity coefficient (if kahm=1).
      !!
      !! ** Method  :   The harmonic operator rotated along geopotential 
      !!      surfaces is applied to (pu,pv) using the slopes of geopotential
      !!      surfaces computed in inildf routine. The result is provided in
      !!      (plu,plv) arrays. It is computed in 2 stepv:
      !!
      !!      First step: horizontal part of the operator. It is computed on
      !!      ==========  pu as follows (idem on pv)
      !!      horizontal fluxes :
      !!         zftu = e2u*e3u/e1u di[ pu ] - e2u*uslp dk[ mi(mk(pu)) ]
      !!         zftv = e1v*e3v/e2v dj[ pu ] - e1v*vslp dk[ mj(mk(pu)) ]
      !!      take the horizontal divergence of the fluxes (no divided by
      !!      the volume element :
      !!         plu  = di-1[ zftu ] +  dj-1[ zftv ]
      !!
      !!      Second step: vertical part of the operator. It is computed on
      !!      ===========  pu as follows (idem on pv)
      !!      vertical fluxes :
      !!         zftw = e1t*e2t/e3w * (wslpi^2+wslpj^2)  dk-1[ pu ]
      !!              -     e2t     *       wslpi        di[ mi(mk(pu)) ]
      !!              -     e1t     *       wslpj        dj[ mj(mk(pu)) ]
      !!      take the vertical divergence of the fluxes add it to the hori-
      !!      zontal component, divide the result by the volume element and
      !!      if kahm=1, multiply by the eddy diffusivity coefficient:
      !!         plu  = aht / (e1t*e2t*e3t) { plu + dk[ zftw ] }
      !!      else:
      !!         plu  =  1  / (e1t*e2t*e3t) { plu + dk[ zftw ] }
      !!
      !! ** Action :
      !!        plu, plv        : partial harmonic operator applied to
      !!                          pu and pv (all the components except
      !!                          second order vertical derivative term)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pu , pv    ! 1st call: before horizontal velocity 
      !                                                               ! 2nd call: ahm x these fields
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(  out) ::   plu, plv   ! partial harmonic operator applied to
      !                                                               ! pu and pv (all the components except
      !                                                               ! second order vertical derivative term)
      INTEGER                         , INTENT(in   ) ::   kahm       ! =1 1st call ; =2 2nd call
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zabe1 , zabe2 , zcof1 , zcof2        ! local scalar
      REAL(wp) ::   zcoef0, zcoef3, zcoef4               !   -      -
      REAL(wp) ::   zbur, zbvr, zmkt, zmkf, zuav, zvav   !   -      -
      REAL(wp) ::   zuwslpi, zuwslpj, zvwslpi, zvwslpj   !   -      -
      !
      REAL(wp), POINTER, DIMENSION(:,:) :: ziut, zjuf, zjvt, zivf, zdku, zdk1u, zdkv, zdk1v
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('ldfguv')
      !
      CALL wrk_alloc( jpi, jpj, ziut, zjuf, zjvt, zivf, zdku, zdk1u, zdkv, zdk1v ) 
      !
      !                               ! ********** !   ! ===============
      DO jk = 1, jpkm1                ! First step !   ! Horizontal slab
         !                            ! ********** !   ! ===============

         ! I.1 Vertical gradient of pu and pv at level jk and jk+1
         ! -------------------------------------------------------
         ! surface boundary condition: zdku(jk=1)=zdku(jk=2)
         !                             zdkv(jk=1)=zdkv(jk=2)

         zdk1u(:,:) = ( pu(:,:,jk) - pu(:,:,jk+1) ) * umask(:,:,jk+1)
         zdk1v(:,:) = ( pv(:,:,jk) - pv(:,:,jk+1) ) * vmask(:,:,jk+1)

         IF( jk == 1 ) THEN
            zdku(:,:) = zdk1u(:,:)
            zdkv(:,:) = zdk1v(:,:)
         ELSE
            zdku(:,:) = ( pu(:,:,jk-1) - pu(:,:,jk) ) * umask(:,:,jk)
            zdkv(:,:) = ( pv(:,:,jk-1) - pv(:,:,jk) ) * vmask(:,:,jk)
         ENDIF

         !                                -----f-----
         ! I.2 Horizontal fluxes on U          |
         ! ------------------------===     t   u   t
         !                                     |
         ! i-flux at t-point              -----f-----
         DO jj = 1, jpjm1
            DO ji = 2, jpi
               zabe1 = e2t(ji,jj) * e3t_0(ji,jj,jk) / e1t(ji,jj)

               zmkt  = 1./MAX(  umask(ji-1,jj,jk  )+umask(ji,jj,jk+1)   &
                              + umask(ji-1,jj,jk+1)+umask(ji,jj,jk  ), 1. )

               zcof1 = -e2t(ji,jj) * zmkt   &
                     * 0.5  * ( uslp(ji-1,jj,jk) + uslp(ji,jj,jk) )

               ziut(ji,jj) = tmask(ji,jj,jk) *   &
                           (  zabe1 * ( pu(ji,jj,jk) - pu(ji-1,jj,jk) )   &
                            + zcof1 * ( zdku (ji,jj) + zdk1u(ji-1,jj)     &
                                       +zdk1u(ji,jj) + zdku (ji-1,jj) )  )
            END DO
         END DO

         ! j-flux at f-point
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               zabe2 = e1f(ji,jj) * e3f_0(ji,jj,jk) / e2f(ji,jj)

               zmkf  = 1./MAX(  umask(ji,jj+1,jk  )+umask(ji,jj,jk+1)   &
                              + umask(ji,jj+1,jk+1)+umask(ji,jj,jk  ), 1. )

               zcof2 = -e1f(ji,jj) * zmkf   &
                     * 0.5  * ( vslp(ji+1,jj,jk) + vslp(ji,jj,jk) )

               zjuf(ji,jj) = fmask(ji,jj,jk) *   &
                           (  zabe2 * ( pu(ji,jj+1,jk) - pu(ji,jj,jk) )   &
                            + zcof2 * ( zdku (ji,jj+1) + zdk1u(ji,jj)     &
                                       +zdk1u(ji,jj+1) + zdku (ji,jj) )  )
            END DO
         END DO

         !                                 |   t   |
         ! I.3 Horizontal fluxes on V      |       |
         ! ------------------------===     f---v---f
         !                                 |       |
         ! i-flux at f-point               |   t   |
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               zabe1 = e2f(ji,jj) * e3f_0(ji,jj,jk) / e1f(ji,jj)

               zmkf  = 1./MAX(  vmask(ji+1,jj,jk  )+vmask(ji,jj,jk+1)   &
                              + vmask(ji+1,jj,jk+1)+vmask(ji,jj,jk  ), 1. )

               zcof1 = -e2f(ji,jj) * zmkf   &
                     * 0.5 * ( uslp(ji,jj+1,jk) + uslp(ji,jj,jk) )

               zivf(ji,jj) = fmask(ji,jj,jk) *   &
                           (  zabe1 * ( pu(ji+1,jj,jk) - pu(ji,jj,jk) )   &
                            + zcof1 * ( zdku (ji,jj) + zdk1u(ji+1,jj)     &
                                       +zdk1u(ji,jj) + zdku (ji+1,jj) )  )
            END DO
         END DO

         ! j-flux at t-point
         DO jj = 2, jpj
            DO ji = 1, jpim1
               zabe2 = e1t(ji,jj) * e3t_0(ji,jj,jk) / e2t(ji,jj)

               zmkt  = 1./MAX(  vmask(ji,jj-1,jk  )+vmask(ji,jj,jk+1)   &
                              + vmask(ji,jj-1,jk+1)+vmask(ji,jj,jk  ), 1. )

               zcof2 = -e1t(ji,jj) * zmkt   &
                     * 0.5 * ( vslp(ji,jj-1,jk) + vslp(ji,jj,jk) )

               zjvt(ji,jj) = tmask(ji,jj,jk) *   &
                           (  zabe2 * ( pu(ji,jj,jk) - pu(ji,jj-1,jk) )   &
                            + zcof2 * ( zdku (ji,jj-1) + zdk1u(ji,jj)     &
                                       +zdk1u(ji,jj-1) + zdku (ji,jj) )  )
            END DO
         END DO


         ! I.4 Second derivative (divergence) (not divided by the volume)
         ! ---------------------

         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               plu(ji,jj,jk) = ziut (ji+1,jj) - ziut (ji,jj  )   &
                             + zjuf (ji  ,jj) - zjuf (ji,jj-1)
               plv(ji,jj,jk) = zivf (ji,jj  ) - zivf (ji-1,jj)   &
                             + zjvt (ji,jj+1) - zjvt (ji,jj  ) 
            END DO
         END DO

         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,synchro,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

      !                             ! ************ !   ! ===============
      DO jj = 2, jpjm1              !  Second step !   ! Horizontal slab
         !                          ! ************ !   ! ===============

         ! II.1 horizontal (pu,pv) gradients
         ! ---------------------------------

         DO jk = 1, jpk
            DO ji = 2, jpi
               ! i-gradient of u at jj
               zdiu (ji,jk) = tmask(ji,jj  ,jk) * ( pu(ji,jj  ,jk) - pu(ji-1,jj  ,jk) )
               ! j-gradient of u and v at jj
               zdju (ji,jk) = fmask(ji,jj  ,jk) * ( pu(ji,jj+1,jk) - pu(ji  ,jj  ,jk) )
               zdjv (ji,jk) = tmask(ji,jj  ,jk) * ( pv(ji,jj  ,jk) - pv(ji  ,jj-1,jk) )
               ! j-gradient of u and v at jj+1
               zdj1u(ji,jk) = fmask(ji,jj-1,jk) * ( pu(ji,jj  ,jk) - pu(ji  ,jj-1,jk) )
               zdj1v(ji,jk) = tmask(ji,jj+1,jk) * ( pv(ji,jj+1,jk) - pv(ji  ,jj  ,jk) )
            END DO
         END DO
         DO jk = 1, jpk
            DO ji = 1, jpim1
               ! i-gradient of v at jj
               zdiv (ji,jk) = fmask(ji,jj  ,jk) * ( pv(ji+1,jj,jk) - pv(ji  ,jj  ,jk) )
            END DO
         END DO


         ! II.2 Vertical fluxes
         ! --------------------

         ! Surface and bottom vertical fluxes set to zero

         zfuw(:, 1 ) = 0.e0
         zfvw(:, 1 ) = 0.e0
         zfuw(:,jpk) = 0.e0
         zfvw(:,jpk) = 0.e0

         ! interior (2=<jk=<jpk-1) on pu field

         DO jk = 2, jpkm1
            DO ji = 2, jpim1
               ! i- and j-slopes at uw-point
               zuwslpi = 0.5 * ( wslpi(ji+1,jj,jk) + wslpi(ji,jj,jk) )
               zuwslpj = 0.5 * ( wslpj(ji+1,jj,jk) + wslpj(ji,jj,jk) )
               ! coef. for the vertical dirative
               zcoef0 = e1u(ji,jj) * e2u(ji,jj) / e3u_0(ji,jj,jk)   &
                      * ( zuwslpi * zuwslpi + zuwslpj * zuwslpj )
               ! weights for the i-k, j-k averaging at t- and f-points, resp.
               zmkt = 1./MAX(  tmask(ji,jj,jk-1)+tmask(ji+1,jj,jk-1)   &
                             + tmask(ji,jj,jk  )+tmask(ji+1,jj,jk  ), 1. )
               zmkf = 1./MAX(  fmask(ji,jj-1,jk-1)+fmask(ji,jj,jk-1)   &
                             + fmask(ji,jj-1,jk  )+fmask(ji,jj,jk  ), 1. )
               ! coef. for the horitontal derivative
               zcoef3 = - e2u(ji,jj) * zmkt * zuwslpi
               zcoef4 = - e1u(ji,jj) * zmkf * zuwslpj
               ! vertical flux on u field
               zfuw(ji,jk) = umask(ji,jj,jk) *   &
                           (  zcoef0 * ( pu  (ji,jj,jk-1) - pu  (ji,jj,jk) )   &
                            + zcoef3 * ( zdiu (ji,jk-1) + zdiu (ji+1,jk-1)     &
                                        +zdiu (ji,jk  ) + zdiu (ji+1,jk  ) )   &
                            + zcoef4 * ( zdj1u(ji,jk-1) + zdju (ji  ,jk-1)     &
                                        +zdj1u(ji,jk  ) + zdju (ji  ,jk  ) ) )
            END DO
         END DO

         ! interior (2=<jk=<jpk-1) on pv field

         DO jk = 2, jpkm1
            DO ji = 2, jpim1
               ! i- and j-slopes at vw-point
               zvwslpi = 0.5 * ( wslpi(ji,jj+1,jk) + wslpi(ji,jj,jk) )
               zvwslpj = 0.5 * ( wslpj(ji,jj+1,jk) + wslpj(ji,jj,jk) )
               ! coef. for the vertical derivative
               zcoef0 = e1v(ji,jj) * e2v(ji,jj) / e3v_0(ji,jj,jk)   &
                      * ( zvwslpi * zvwslpi + zvwslpj * zvwslpj )
               ! weights for the i-k, j-k averaging at f- and t-points, resp.
               zmkf = 1./MAX(  fmask(ji-1,jj,jk-1)+fmask(ji,jj,jk-1)   &
                             + fmask(ji-1,jj,jk  )+fmask(ji,jj,jk  ), 1. )
               zmkt = 1./MAX(  tmask(ji,jj,jk-1)+tmask(ji,jj+1,jk-1)   &
                             + tmask(ji,jj,jk  )+tmask(ji,jj+1,jk  ), 1. )
               ! coef. for the horizontal derivatives
               zcoef3 = - e2v(ji,jj) * zmkf * zvwslpi
               zcoef4 = - e1v(ji,jj) * zmkt * zvwslpj
               ! vertical flux on pv field
               zfvw(ji,jk) = vmask(ji,jj,jk) *   &
                           (  zcoef0 * ( pv  (ji,jj,jk-1) - pv  (ji,jj,jk) )   &
                            + zcoef3 * ( zdiv (ji,jk-1) + zdiv (ji-1,jk-1)     &
                                        +zdiv (ji,jk  ) + zdiv (ji-1,jk  ) )   &
                            + zcoef4 * ( zdjv (ji,jk-1) + zdj1v(ji  ,jk-1)     &
                                        +zdjv (ji,jk  ) + zdj1v(ji  ,jk  ) )  )
            END DO
         END DO


         ! II.3 Divergence of vertical fluxes added to the horizontal divergence
         ! ---------------------------------------------------------------------
         IF( (kahm -nkahm_smag) ==1 ) THEN
            ! multiply the laplacian by the eddy viscosity coefficient
            DO jk = 1, jpkm1
               DO ji = 2, jpim1
                  ! eddy coef. divided by the volume element
                  zbur = ahm3(ji,jj,jk) / ( e1u(ji,jj)*e2u(ji,jj)*e3u_0(ji,jj,jk) )
                  zbvr = ahm4(ji,jj,jk) / ( e1v(ji,jj)*e2v(ji,jj)*e3v_0(ji,jj,jk) )
                  ! vertical divergence
                  zuav = zfuw(ji,jk) - zfuw(ji,jk+1)
                  zvav = zfvw(ji,jk) - zfvw(ji,jk+1)
                  ! harmonic operator applied to (pu,pv) and multiply by ahm
                  plu(ji,jj,jk) = ( plu(ji,jj,jk) + zuav ) * zbur
                  plv(ji,jj,jk) = ( plv(ji,jj,jk) + zvav ) * zbvr
               END DO
            END DO
         ELSEIF( (kahm +nkahm_smag ) == 2 ) THEN
            ! second call, no multiplication
            DO jk = 1, jpkm1
               DO ji = 2, jpim1
                  ! inverse of the volume element
                  zbur = 1. / ( e1u(ji,jj)*e2u(ji,jj)*e3u_0(ji,jj,jk) )
                  zbvr = 1. / ( e1v(ji,jj)*e2v(ji,jj)*e3v_0(ji,jj,jk) )
                  ! vertical divergence
                  zuav = zfuw(ji,jk) - zfuw(ji,jk+1)
                  zvav = zfvw(ji,jk) - zfvw(ji,jk+1)
                  ! harmonic operator applied to (pu,pv) 
                  plu(ji,jj,jk) = ( plu(ji,jj,jk) + zuav ) * zbur
                  plv(ji,jj,jk) = ( plv(ji,jj,jk) + zvav ) * zbvr
               END DO
            END DO
         ELSE
            IF(lwp)WRITE(numout,*) ' ldfguv: kahm= 1 or 2, here =', kahm
            IF(lwp)WRITE(numout,*) '         We stop'
            STOP 'ldfguv'
         ENDIF
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      CALL wrk_dealloc( jpi, jpj, ziut, zjuf, zjvt, zivf, zdku, zdk1u, zdkv, zdk1v ) 
      !
      IF( nn_timing == 1 )  CALL timing_stop('ldfguv')
      !
   END SUBROUTINE ldfguv


   !!======================================================================
END MODULE dynldf_bilapg
