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

MODULE dynvor
   !!======================================================================
   !!                       ***  MODULE  dynvor  ***
   !! Ocean dynamics: Update the momentum trend with the relative and
   !!                 planetary vorticity trends
   !!======================================================================
   !! History :  OPA  ! 1989-12  (P. Andrich)  vor_ens: Original code
   !!            5.0  ! 1991-11  (G. Madec) vor_ene, vor_mix: Original code
   !!            6.0  ! 1996-01  (G. Madec)  s-coord, suppress work arrays
   !!   NEMO     0.5  ! 2002-08  (G. Madec)  F90: Free form and module
   !!            1.0  ! 2004-02  (G. Madec)  vor_een: Original code
   !!             -   ! 2003-08  (G. Madec)  add vor_ctl
   !!             -   ! 2005-11  (G. Madec)  add dyn_vor (new step architecture)
   !!            2.0  ! 2006-11  (G. Madec)  flux form advection: add metric term
   !!            3.2  ! 2009-04  (R. Benshila)  vvl: correction of een scheme
   !!            3.3  ! 2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!            3.7  ! 2014-04  (G. Madec) trend simplification: suppress jpdyn_trd_dat vorticity 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_vor      : Update the momentum trend with the vorticity trend
   !!       vor_ens  : enstrophy conserving scheme       (ln_dynvor_ens=T)
   !!       vor_ene  : energy conserving scheme          (ln_dynvor_ene=T)
   !!       vor_mix  : mixed enstrophy/energy conserving (ln_dynvor_mix=T)
   !!       vor_een  : energy and enstrophy conserving   (ln_dynvor_een=T)
   !!   dyn_vor_init : set and control of the different vorticity option
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE dommsk         ! ocean mask
   USE dynadv         ! momentum advection (use ln_dynadv_vec value)
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE prtctl         ! Print control
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE wrk_nemo       ! Memory Allocation
   USE timing         ! Timing


   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_vor        ! routine called by step.F90
   PUBLIC   dyn_vor_init   ! routine called by opa.F90

   !                                   !!* Namelist namdyn_vor: vorticity term
   LOGICAL, PUBLIC ::   ln_dynvor_ene   !: energy conserving scheme
   LOGICAL, PUBLIC ::   ln_dynvor_ens   !: enstrophy conserving scheme
   LOGICAL, PUBLIC ::   ln_dynvor_mix   !: mixed scheme
   LOGICAL, PUBLIC ::   ln_dynvor_een   !: energy and enstrophy conserving scheme
   LOGICAL, PUBLIC ::   ln_dynvor_een_old !: energy and enstrophy conserving scheme (original formulation)

   INTEGER ::   nvor = 0   ! type of vorticity trend used
   INTEGER ::   ncor = 1   ! coriolis
   INTEGER ::   nrvm = 2   ! =2 relative vorticity ; =3 metric term
   INTEGER ::   ntot = 4   ! =4 total vorticity (relative + planetary) ; =5 coriolis + metric term

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
   !! $Id: dynvor.F90 5029 2015-01-14 16:32:37Z davestorkey $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_vor( kt )
      !!----------------------------------------------------------------------
      !!
      !! ** Purpose :   compute the lateral ocean tracer physics.
      !!
      !! ** Action : - Update (ua,va) with the now vorticity term trend
      !!             - save the trends in (ztrdu,ztrdv) in 2 parts (relative
      !!               and planetary vorticity trends) and send them to trd_dyn 
      !!               for futher diagnostics (l_trddyn=T)
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ztrdu, ztrdv
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_vor')
      !
      IF( l_trddyn )   CALL wrk_alloc( jpi,jpj,jpk, ztrdu, ztrdv )
      !
      !                                          ! vorticity term 
      SELECT CASE ( nvor )                       ! compute the vorticity trend and add it to the general trend
      !
      CASE ( -1 )                                      ! esopa: test all possibility with control print
         CALL vor_ene( kt, ntot, ua, va )
         CALL prt_ctl( tab3d_1=ua, clinfo1=' vor0 - Ua: ', mask1=umask, &
            &          tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
         CALL vor_ens( kt, ntot, ua, va )
         CALL prt_ctl( tab3d_1=ua, clinfo1=' vor1 - Ua: ', mask1=umask, &
            &          tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
         CALL vor_mix( kt )
         CALL prt_ctl( tab3d_1=ua, clinfo1=' vor2 - Ua: ', mask1=umask, &
            &          tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
         CALL vor_een( kt, ntot, ua, va )
         CALL prt_ctl( tab3d_1=ua, clinfo1=' vor3 - Ua: ', mask1=umask, &
            &          tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
         !
      CASE ( 0 )                                       ! energy conserving scheme
         IF( l_trddyn )   THEN
            ztrdu(:,:,:) = ua(:,:,:)
            ztrdv(:,:,:) = va(:,:,:)
            CALL vor_ene( kt, nrvm, ua, va )                ! relative vorticity or metric trend
            ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
            ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
            CALL trd_dyn( ztrdu, ztrdv, jpdyn_rvo, kt )
            ztrdu(:,:,:) = ua(:,:,:)
            ztrdv(:,:,:) = va(:,:,:)
            CALL vor_ene( kt, ncor, ua, va )                ! planetary vorticity trend
            ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
            ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
            CALL trd_dyn( ztrdu, ztrdv, jpdyn_pvo, kt )
         ELSE
            CALL vor_ene( kt, ntot, ua, va )                ! total vorticity
         ENDIF
         !
      CASE ( 1 )                                       ! enstrophy conserving scheme
         IF( l_trddyn )   THEN    
            ztrdu(:,:,:) = ua(:,:,:)
            ztrdv(:,:,:) = va(:,:,:)
            CALL vor_ens( kt, nrvm, ua, va )                ! relative vorticity or metric trend
            ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
            ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
            CALL trd_dyn( ztrdu, ztrdv, jpdyn_rvo, kt )
            ztrdu(:,:,:) = ua(:,:,:)
            ztrdv(:,:,:) = va(:,:,:)
            CALL vor_ens( kt, ncor, ua, va )                ! planetary vorticity trend
            ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
            ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
            CALL trd_dyn( ztrdu, ztrdv, jpdyn_pvo, kt )
         ELSE
            CALL vor_ens( kt, ntot, ua, va )                ! total vorticity
         ENDIF
         !
      CASE ( 2 )                                       ! mixed ene-ens scheme
         IF( l_trddyn )   THEN
            ztrdu(:,:,:) = ua(:,:,:)
            ztrdv(:,:,:) = va(:,:,:)
            CALL vor_ens( kt, nrvm, ua, va )                ! relative vorticity or metric trend (ens)
            ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
            ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
            CALL trd_dyn( ztrdu, ztrdv, jpdyn_rvo, kt )
            ztrdu(:,:,:) = ua(:,:,:)
            ztrdv(:,:,:) = va(:,:,:)
            CALL vor_ene( kt, ncor, ua, va )                ! planetary vorticity trend (ene)
            ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
            ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
            CALL trd_dyn( ztrdu, ztrdv, jpdyn_pvo, kt )
         ELSE
            CALL vor_mix( kt )                               ! total vorticity (mix=ens-ene)
         ENDIF
         !
      CASE ( 3 )                                       ! energy and enstrophy conserving scheme
         IF( l_trddyn )   THEN
            ztrdu(:,:,:) = ua(:,:,:)
            ztrdv(:,:,:) = va(:,:,:)
            CALL vor_een( kt, nrvm, ua, va )                ! relative vorticity or metric trend
            ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
            ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
            CALL trd_dyn( ztrdu, ztrdv, jpdyn_rvo, kt )
            ztrdu(:,:,:) = ua(:,:,:)
            ztrdv(:,:,:) = va(:,:,:)
            CALL vor_een( kt, ncor, ua, va )                ! planetary vorticity trend
            ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
            ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
            CALL trd_dyn( ztrdu, ztrdv, jpdyn_pvo, kt )
         ELSE
            CALL vor_een( kt, ntot, ua, va )                ! total vorticity
         ENDIF
         !
      END SELECT
      !
      !                       ! print sum trends (used for debugging)
      IF(ln_ctl) CALL prt_ctl( tab3d_1=ua, clinfo1=' vor  - Ua: ', mask1=umask,               &
         &                     tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( l_trddyn )   CALL wrk_dealloc( jpi,jpj,jpk, ztrdu, ztrdv )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_vor')
      !
   END SUBROUTINE dyn_vor


   SUBROUTINE vor_ene( kt, kvor, pua, pva )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE vor_ene  ***
      !!
      !! ** Purpose :   Compute the now total vorticity trend and add it to 
      !!      the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time) 
      !!      and the Sadourny (1975) flux form formulation : conserves the
      !!      horizontal kinetic energy.
      !!      The trend of the vorticity term is given by:
      !!       * s-coordinate (ln_sco=T), the e3. are inside the derivatives:
      !!          voru = 1/e1u  mj-1[ (rotn+f)/e3f  mi(e1v*e3v vn) ]
      !!          vorv = 1/e2v  mi-1[ (rotn+f)/e3f  mj(e2u*e3u un) ]
      !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
      !!          voru = 1/e1u  mj-1[ (rotn+f)  mi(e1v vn) ]
      !!          vorv = 1/e2v  mi-1[ (rotn+f)  mj(e2u un) ]
      !!      Add this trend to the general momentum trend (ua,va):
      !!          (ua,va) = (ua,va) + ( voru , vorv )
      !!
      !! ** Action : - Update (ua,va) with the now vorticity term trend
      !!
      !! References : Sadourny, r., 1975, j. atmos. sciences, 32, 680-689.
      !!----------------------------------------------------------------------
      !
      INTEGER , INTENT(in   )                         ::   kt     ! ocean time-step index
      INTEGER , INTENT(in   )                         ::   kvor   ! =ncor (planetary) ; =ntot (total) ;
      !                                                           ! =nrvm (relative vorticity or metric)
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   pua    ! total u-trend
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   pva    ! total v-trend
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zx1, zy1, zfact2, zx2, zy2   ! local scalars
      REAL(wp), POINTER, DIMENSION(:,:) :: zwx, zwy, zwz
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('vor_ene')
      !
      CALL wrk_alloc( jpi, jpj, zwx, zwy, zwz ) 
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:vor_ene : vorticity term: energy conserving scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF

      zfact2 = 0.5 * 0.5      ! Local constant initialization

!CDIR PARALLEL DO PRIVATE( zwx, zwy, zwz )
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         !
         ! Potential vorticity and horizontal fluxes
         ! -----------------------------------------
         SELECT CASE( kvor )      ! vorticity considered
         CASE ( 1 )   ;   zwz(:,:) =                  ff(:,:)      ! planetary vorticity (Coriolis)
         CASE ( 2 )   ;   zwz(:,:) =   rotn(:,:,jk)                ! relative  vorticity
         CASE ( 3 )                                                ! metric term
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwz(ji,jj) = (   ( vn(ji+1,jj  ,jk) + vn (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                       &         - ( un(ji  ,jj+1,jk) + un (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                       &     * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) )
               END DO
            END DO
         CASE ( 4 )   ;   zwz(:,:) = ( rotn(:,:,jk) + ff(:,:) )    ! total (relative + planetary vorticity)
         CASE ( 5 )                                                ! total (coriolis + metric)
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwz(ji,jj) = ( ff (ji,jj)                                                                       &
                       &       + (   ( vn(ji+1,jj  ,jk) + vn (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                       &           - ( un(ji  ,jj+1,jk) + un (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                       &       * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) )                                               &
                       &       )
               END DO
            END DO
         END SELECT

         IF( ln_sco ) THEN
            zwz(:,:) = zwz(:,:) / e3f_0(:,:,jk)
            zwx(:,:) = e2u(:,:) * e3u_0(:,:,jk) * un(:,:,jk)
            zwy(:,:) = e1v(:,:) * e3v_0(:,:,jk) * vn(:,:,jk)
         ELSE
            zwx(:,:) = e2u(:,:) * un(:,:,jk)
            zwy(:,:) = e1v(:,:) * vn(:,:,jk)
         ENDIF

         ! Compute and add the vorticity term trend
         ! ----------------------------------------
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zy1 = zwy(ji,jj-1) + zwy(ji+1,jj-1)
               zy2 = zwy(ji,jj  ) + zwy(ji+1,jj  )
               zx1 = zwx(ji-1,jj) + zwx(ji-1,jj+1)
               zx2 = zwx(ji  ,jj) + zwx(ji  ,jj+1)
               pua(ji,jj,jk) = pua(ji,jj,jk) + zfact2 / e1u(ji,jj) * ( zwz(ji  ,jj-1) * zy1 + zwz(ji,jj) * zy2 )
               pva(ji,jj,jk) = pva(ji,jj,jk) - zfact2 / e2v(ji,jj) * ( zwz(ji-1,jj  ) * zx1 + zwz(ji,jj) * zx2 ) 
            END DO  
         END DO  
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      CALL wrk_dealloc( jpi, jpj, zwx, zwy, zwz ) 
      !
      IF( nn_timing == 1 )  CALL timing_stop('vor_ene')
      !
   END SUBROUTINE vor_ene


   SUBROUTINE vor_mix( kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE vor_mix  ***
      !!
      !! ** Purpose :   Compute the now total vorticity trend and add it to
      !!      the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time)
      !!      Mixte formulation : conserves the potential enstrophy of a hori-
      !!      zontally non-divergent flow for (rotzu x uh), the relative vor-
      !!      ticity term and the horizontal kinetic energy for (f x uh), the
      !!      coriolis term. the now trend of the vorticity term is given by:
      !!       * s-coordinate (ln_sco=T), the e3. are inside the derivatives:
      !!          voru = 1/e1u  mj-1(rotn/e3f) mj-1[ mi(e1v*e3v vn) ]
      !!              +1/e1u  mj-1[ f/e3f          mi(e1v*e3v vn) ]
      !!          vorv = 1/e2v  mi-1(rotn/e3f) mi-1[ mj(e2u*e3u un) ]
      !!              +1/e2v  mi-1[ f/e3f          mj(e2u*e3u un) ]
      !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
      !!          voru = 1/e1u  mj-1(rotn) mj-1[ mi(e1v vn) ]
      !!              +1/e1u  mj-1[ f          mi(e1v vn) ]
      !!          vorv = 1/e2v  mi-1(rotn) mi-1[ mj(e2u un) ]
      !!              +1/e2v  mi-1[ f          mj(e2u un) ]
      !!      Add this now trend to the general momentum trend (ua,va):
      !!          (ua,va) = (ua,va) + ( voru , vorv )
      !!
      !! ** Action : - Update (ua,va) arrays with the now vorticity term trend
      !!
      !! References : Sadourny, r., 1975, j. atmos. sciences, 32, 680-689.
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean timestep index
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zfact1, zua, zcua, zx1, zy1   ! local scalars
      REAL(wp) ::   zfact2, zva, zcva, zx2, zy2   !   -      -
      REAL(wp), POINTER, DIMENSION(:,:) :: zwx, zwy, zwz, zww
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('vor_mix')
      !
      CALL wrk_alloc( jpi, jpj, zwx, zwy, zwz, zww ) 
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:vor_mix : vorticity term: mixed energy/enstrophy conserving scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF

      zfact1 = 0.5 * 0.25      ! Local constant initialization
      zfact2 = 0.5 * 0.5

!CDIR PARALLEL DO PRIVATE( zwx, zwy, zwz, zww )
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         !
         ! Relative and planetary potential vorticity and horizontal fluxes
         ! ----------------------------------------------------------------
         IF( ln_sco ) THEN        
            IF( ln_dynadv_vec ) THEN
               zww(:,:) = rotn(:,:,jk) / e3f_0(:,:,jk)
            ELSE                       
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1   ! vector opt.
                     zww(ji,jj) = (   ( vn(ji+1,jj  ,jk) + vn (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                        &           - ( un(ji  ,jj+1,jk) + un (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                        &       * 0.5 / ( e1f(ji,jj) * e2f (ji,jj) * e3f_0(ji,jj,jk) )
                  END DO
               END DO
            ENDIF
            zwz(:,:) = ff  (:,:)    / e3f_0(:,:,jk)
            zwx(:,:) = e2u(:,:) * e3u_0(:,:,jk) * un(:,:,jk)
            zwy(:,:) = e1v(:,:) * e3v_0(:,:,jk) * vn(:,:,jk)
         ELSE
            IF( ln_dynadv_vec ) THEN
               zww(:,:) = rotn(:,:,jk)
            ELSE                       
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1   ! vector opt.
                     zww(ji,jj) = (   ( vn(ji+1,jj  ,jk) + vn (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                        &           - ( un(ji  ,jj+1,jk) + un (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                        &       * 0.5 / ( e1f(ji,jj) * e2f (ji,jj) )
                  END DO
               END DO
            ENDIF
            zwz(:,:) = ff (:,:)
            zwx(:,:) = e2u(:,:) * un(:,:,jk)
            zwy(:,:) = e1v(:,:) * vn(:,:,jk)
         ENDIF

         ! Compute and add the vorticity term trend
         ! ----------------------------------------
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zy1 = ( zwy(ji,jj-1) + zwy(ji+1,jj-1) ) / e1u(ji,jj)
               zy2 = ( zwy(ji,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
               zx1 = ( zwx(ji-1,jj) + zwx(ji-1,jj+1) ) / e2v(ji,jj)
               zx2 = ( zwx(ji  ,jj) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
               ! enstrophy conserving formulation for relative vorticity term
               zua = zfact1 * ( zww(ji  ,jj-1) + zww(ji,jj) ) * ( zy1 + zy2 )
               zva =-zfact1 * ( zww(ji-1,jj  ) + zww(ji,jj) ) * ( zx1 + zx2 )
               ! energy conserving formulation for planetary vorticity term
               zcua = zfact2 * ( zwz(ji  ,jj-1) * zy1 + zwz(ji,jj) * zy2 )
               zcva =-zfact2 * ( zwz(ji-1,jj  ) * zx1 + zwz(ji,jj) * zx2 )
               ! mixed vorticity trend added to the momentum trends
               ua(ji,jj,jk) = ua(ji,jj,jk) + zcua + zua
               va(ji,jj,jk) = va(ji,jj,jk) + zcva + zva
            END DO  
         END DO  
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      CALL wrk_dealloc( jpi, jpj, zwx, zwy, zwz, zww ) 
      !
      IF( nn_timing == 1 )  CALL timing_stop('vor_mix')
      !
   END SUBROUTINE vor_mix


   SUBROUTINE vor_ens( kt, kvor, pua, pva )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE vor_ens  ***
      !!
      !! ** Purpose :   Compute the now total vorticity trend and add it to
      !!      the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time)
      !!      and the Sadourny (1975) flux FORM formulation : conserves the
      !!      potential enstrophy of a horizontally non-divergent flow. the
      !!      trend of the vorticity term is given by:
      !!       * s-coordinate (ln_sco=T), the e3. are inside the derivative:
      !!          voru = 1/e1u  mj-1[ (rotn+f)/e3f ]  mj-1[ mi(e1v*e3v vn) ]
      !!          vorv = 1/e2v  mi-1[ (rotn+f)/e3f ]  mi-1[ mj(e2u*e3u un) ]
      !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
      !!          voru = 1/e1u  mj-1[ rotn+f ]  mj-1[ mi(e1v vn) ]
      !!          vorv = 1/e2v  mi-1[ rotn+f ]  mi-1[ mj(e2u un) ]
      !!      Add this trend to the general momentum trend (ua,va):
      !!          (ua,va) = (ua,va) + ( voru , vorv )
      !!
      !! ** Action : - Update (ua,va) arrays with the now vorticity term trend
      !!
      !! References : Sadourny, r., 1975, j. atmos. sciences, 32, 680-689.
      !!----------------------------------------------------------------------
      !
      INTEGER , INTENT(in   )                         ::   kt     ! ocean time-step index
      INTEGER , INTENT(in   )                         ::   kvor   ! =ncor (planetary) ; =ntot (total) ;
         !                                                        ! =nrvm (relative vorticity or metric)
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   pua    ! total u-trend
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   pva    ! total v-trend
      !
      INTEGER  ::   ji, jj, jk           ! dummy loop indices
      REAL(wp) ::   zfact1, zuav, zvau   ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:) :: zwx, zwy, zwz, zww
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('vor_ens')
      !
      CALL wrk_alloc( jpi, jpj, zwx, zwy, zwz ) 
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:vor_ens : vorticity term: enstrophy conserving scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF

      zfact1 = 0.5 * 0.25      ! Local constant initialization

!CDIR PARALLEL DO PRIVATE( zwx, zwy, zwz )
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         !
         ! Potential vorticity and horizontal fluxes
         ! -----------------------------------------
         SELECT CASE( kvor )      ! vorticity considered
         CASE ( 1 )   ;   zwz(:,:) =                  ff(:,:)      ! planetary vorticity (Coriolis)
         CASE ( 2 )   ;   zwz(:,:) =   rotn(:,:,jk)                ! relative  vorticity
         CASE ( 3 )                                                ! metric term
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwz(ji,jj) = (   ( vn(ji+1,jj  ,jk) + vn (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                       &         - ( un(ji  ,jj+1,jk) + un (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                       &     * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) )
               END DO
            END DO
         CASE ( 4 )   ;   zwz(:,:) = ( rotn(:,:,jk) + ff(:,:) )    ! total (relative + planetary vorticity)
         CASE ( 5 )                                                ! total (coriolis + metric)
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwz(ji,jj) = ( ff (ji,jj)                                                                       &
                       &       + (   ( vn(ji+1,jj  ,jk) + vn (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                       &           - ( un(ji  ,jj+1,jk) + un (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                       &       * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) )                                                &
                       &       )
               END DO
            END DO
         END SELECT
         !
         IF( ln_sco ) THEN
            DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop 
               DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
                  zwz(ji,jj) = zwz(ji,jj) / e3f_0(ji,jj,jk)
                  zwx(ji,jj) = e2u(ji,jj) * e3u_0(ji,jj,jk) * un(ji,jj,jk)
                  zwy(ji,jj) = e1v(ji,jj) * e3v_0(ji,jj,jk) * vn(ji,jj,jk)
               END DO
            END DO
         ELSE
            DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop 
               DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
                  zwx(ji,jj) = e2u(ji,jj) * un(ji,jj,jk)
                  zwy(ji,jj) = e1v(ji,jj) * vn(ji,jj,jk)
               END DO
            END DO
         ENDIF
         !
         ! Compute and add the vorticity term trend
         ! ----------------------------------------
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zuav = zfact1 / e1u(ji,jj) * ( zwy(ji  ,jj-1) + zwy(ji+1,jj-1)   &
                  &                         + zwy(ji  ,jj  ) + zwy(ji+1,jj  ) )
               zvau =-zfact1 / e2v(ji,jj) * ( zwx(ji-1,jj  ) + zwx(ji-1,jj+1)   &
                  &                         + zwx(ji  ,jj  ) + zwx(ji  ,jj+1) )
               pua(ji,jj,jk) = pua(ji,jj,jk) + zuav * ( zwz(ji  ,jj-1) + zwz(ji,jj) )
               pva(ji,jj,jk) = pva(ji,jj,jk) + zvau * ( zwz(ji-1,jj  ) + zwz(ji,jj) )
            END DO  
         END DO  
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      CALL wrk_dealloc( jpi, jpj, zwx, zwy, zwz ) 
      !
      IF( nn_timing == 1 )  CALL timing_stop('vor_ens')
      !
   END SUBROUTINE vor_ens


   SUBROUTINE vor_een( kt, kvor, pua, pva )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE vor_een  ***
      !!
      !! ** Purpose :   Compute the now total vorticity trend and add it to 
      !!      the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time) 
      !!      and the Arakawa and Lamb (1980) flux form formulation : conserves 
      !!      both the horizontal kinetic energy and the potential enstrophy
      !!      when horizontal divergence is zero (see the NEMO documentation)
      !!      Add this trend to the general momentum trend (ua,va).
      !!
      !! ** Action : - Update (ua,va) with the now vorticity term trend
      !!
      !! References : Arakawa and Lamb 1980, Mon. Wea. Rev., 109, 18-36
      !!----------------------------------------------------------------------
      !
      INTEGER , INTENT(in   )                         ::   kt     ! ocean time-step index
      INTEGER , INTENT(in   )                         ::   kvor   ! =ncor (planetary) ; =ntot (total) ;
      !                                                           ! =nrvm (relative vorticity or metric)
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   pua    ! total u-trend
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   pva    ! total v-trend
      !!
      INTEGER  ::   ji, jj, jk                                    ! dummy loop indices
      INTEGER  ::   ierr                                          ! local integer
      REAL(wp) ::   zfac12, zua, zva                              ! local scalars
      REAL(wp) ::   zmsk, ze3                                     ! local scalars
      !                                                           !  3D workspace 
      REAL(wp), POINTER    , DIMENSION(:,:  )         :: zwx, zwy, zwz
      REAL(wp), POINTER    , DIMENSION(:,:  )         :: ztnw, ztne, ztsw, ztse
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:), SAVE   :: ze3f     ! lk_vvl=F, ze3f=1/e3f saved one for all
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('vor_een')
      !
      CALL wrk_alloc( jpi, jpj,      zwx , zwy , zwz        ) 
      CALL wrk_alloc( jpi, jpj,      ztnw, ztne, ztsw, ztse ) 
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:vor_een : vorticity term: energy and enstrophy conserving scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
         IF( .NOT.ALLOCATED(ze3f) ) THEN
            ALLOCATE( ze3f(jpi,jpj,jpk) , STAT=ierr )
            IF( lk_mpp    )   CALL mpp_sum ( ierr )
            IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'dyn:vor_een : unable to allocate arrays' )
         ENDIF
         ze3f(:,:,:) = 0._wp
      ENDIF

      IF( kt == nit000 .OR. lk_vvl ) THEN      ! reciprocal of e3 at F-point (masked averaging of e3t over ocean points)

         IF( ln_dynvor_een_old ) THEN ! original formulation
            DO jk = 1, jpk
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     ze3  = ( e3t_0(ji,jj+1,jk)*tmask(ji,jj+1,jk) + e3t_0(ji+1,jj+1,jk)*tmask(ji+1,jj+1,jk)   &
                        &   + e3t_0(ji,jj  ,jk)*tmask(ji,jj  ,jk) + e3t_0(ji+1,jj  ,jk)*tmask(ji+1,jj  ,jk) )
                     IF( ze3 /= 0._wp )   ze3f(ji,jj,jk) = 4.0_wp / ze3
                  END DO
               END DO
            END DO
         ELSE ! new formulation from NEMO 3.6
            DO jk = 1, jpk
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     ze3  = ( e3t_0(ji,jj+1,jk)*tmask(ji,jj+1,jk) + e3t_0(ji+1,jj+1,jk)*tmask(ji+1,jj+1,jk)   &
                        &   + e3t_0(ji,jj  ,jk)*tmask(ji,jj  ,jk) + e3t_0(ji+1,jj  ,jk)*tmask(ji+1,jj  ,jk) )
                     zmsk = (                   tmask(ji,jj+1,jk) +                     tmask(ji+1,jj+1,jk)   &
                        &                     + tmask(ji,jj  ,jk) +                     tmask(ji+1,jj  ,jk) )
                     IF( ze3 /= 0._wp )   ze3f(ji,jj,jk) = zmsk / ze3
                  END DO
               END DO
            END DO
         ENDIF

         CALL lbc_lnk( ze3f, 'F', 1. )
      ENDIF

      zfac12 = 1._wp / 12._wp    ! Local constant initialization

      
!CDIR PARALLEL DO PRIVATE( zwx, zwy, zwz, ztnw, ztne, ztsw, ztse )
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         
         ! Potential vorticity and horizontal fluxes
         ! -----------------------------------------
         SELECT CASE( kvor )      ! vorticity considered
         CASE ( 1 )                                                ! planetary vorticity (Coriolis)
            zwz(:,:) = ff(:,:)      * ze3f(:,:,jk)
         CASE ( 2 )                                                ! relative  vorticity
            zwz(:,:) = rotn(:,:,jk) * ze3f(:,:,jk)
         CASE ( 3 )                                                ! metric term
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwz(ji,jj) = (   ( vn(ji+1,jj  ,jk) + vn (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                       &         - ( un(ji  ,jj+1,jk) + un (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                       &     * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) ) * ze3f(ji,jj,jk)
               END DO
            END DO
            CALL lbc_lnk( zwz, 'F', 1. )
        CASE ( 4 )                                                ! total (relative + planetary vorticity)
            zwz(:,:) = ( rotn(:,:,jk) + ff(:,:) ) * ze3f(:,:,jk)
         CASE ( 5 )                                                ! total (coriolis + metric)
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwz(ji,jj) = ( ff (ji,jj)                                                                       &
                       &       + (   ( vn(ji+1,jj  ,jk) + vn (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                       &           - ( un(ji  ,jj+1,jk) + un (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                       &       * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) )                                                &
                       &       ) * ze3f(ji,jj,jk)
               END DO
            END DO
            CALL lbc_lnk( zwz, 'F', 1. )
         END SELECT

         zwx(:,:) = e2u(:,:) * e3u_0(:,:,jk) * un(:,:,jk)
         zwy(:,:) = e1v(:,:) * e3v_0(:,:,jk) * vn(:,:,jk)

         ! Compute and add the vorticity term trend
         ! ----------------------------------------
         jj = 2
         ztne(1,:) = 0   ;   ztnw(1,:) = 0   ;   ztse(1,:) = 0   ;   ztsw(1,:) = 0
         DO ji = 2, jpi   
               ztne(ji,jj) = zwz(ji-1,jj  ) + zwz(ji  ,jj  ) + zwz(ji  ,jj-1)
               ztnw(ji,jj) = zwz(ji-1,jj-1) + zwz(ji-1,jj  ) + zwz(ji  ,jj  )
               ztse(ji,jj) = zwz(ji  ,jj  ) + zwz(ji  ,jj-1) + zwz(ji-1,jj-1)
               ztsw(ji,jj) = zwz(ji  ,jj-1) + zwz(ji-1,jj-1) + zwz(ji-1,jj  )
         END DO
         DO jj = 3, jpj
            DO ji = 2, jpi   ! vector opt. ok because we start at jj = 3
               ztne(ji,jj) = zwz(ji-1,jj  ) + zwz(ji  ,jj  ) + zwz(ji  ,jj-1)
               ztnw(ji,jj) = zwz(ji-1,jj-1) + zwz(ji-1,jj  ) + zwz(ji  ,jj  )
               ztse(ji,jj) = zwz(ji  ,jj  ) + zwz(ji  ,jj-1) + zwz(ji-1,jj-1)
               ztsw(ji,jj) = zwz(ji  ,jj-1) + zwz(ji-1,jj-1) + zwz(ji-1,jj  )
            END DO
         END DO
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zua = + zfac12 / e1u(ji,jj) * (  ztne(ji,jj  ) * zwy(ji  ,jj  ) + ztnw(ji+1,jj) * zwy(ji+1,jj  )   &
                  &                           + ztse(ji,jj  ) * zwy(ji  ,jj-1) + ztsw(ji+1,jj) * zwy(ji+1,jj-1) )
               zva = - zfac12 / e2v(ji,jj) * (  ztsw(ji,jj+1) * zwx(ji-1,jj+1) + ztse(ji,jj+1) * zwx(ji  ,jj+1)   &
                  &                           + ztnw(ji,jj  ) * zwx(ji-1,jj  ) + ztne(ji,jj  ) * zwx(ji  ,jj  ) )
               pua(ji,jj,jk) = pua(ji,jj,jk) + zua
               pva(ji,jj,jk) = pva(ji,jj,jk) + zva
            END DO  
         END DO  
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      CALL wrk_dealloc( jpi, jpj,      zwx , zwy , zwz        ) 
      CALL wrk_dealloc( jpi, jpj,      ztnw, ztne, ztsw, ztse ) 
      !
      IF( nn_timing == 1 )  CALL timing_stop('vor_een')
      !
   END SUBROUTINE vor_een


   SUBROUTINE dyn_vor_init
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_vor_init  ***
      !!
      !! ** Purpose :   Control the consistency between cpp options for
      !!              tracer advection schemes
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio          ! local integer
      INTEGER ::   ji, jj, jk      ! dummy loop indices
      INTEGER ::   ios             ! Local integer output status for namelist read
      !!
      NAMELIST/namdyn_vor/ ln_dynvor_ens, ln_dynvor_ene, ln_dynvor_mix, ln_dynvor_een, ln_dynvor_een_old
      !!----------------------------------------------------------------------

      REWIND( numnam_ref )              ! Namelist namdyn_vor in reference namelist : Vorticity scheme options
      READ  ( numnam_ref, namdyn_vor, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namdyn_vor in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namdyn_vor in configuration namelist : Vorticity scheme options
      READ  ( numnam_cfg, namdyn_vor, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namdyn_vor in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namdyn_vor )

      IF(lwp) THEN                    ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_vor_init : vorticity term : read namelist and control the consistency'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '        Namelist namdyn_vor : choice of the vorticity term scheme'
         WRITE(numout,*) '           energy    conserving scheme                ln_dynvor_ene = ', ln_dynvor_ene
         WRITE(numout,*) '           enstrophy conserving scheme                ln_dynvor_ens = ', ln_dynvor_ens
         WRITE(numout,*) '           mixed enstrophy/energy conserving scheme   ln_dynvor_mix = ', ln_dynvor_mix
         WRITE(numout,*) '           enstrophy and energy conserving scheme     ln_dynvor_een = ', ln_dynvor_een
         WRITE(numout,*) '           enstrophy and energy conserving scheme (old) ln_dynvor_een_old= ', ln_dynvor_een_old
      ENDIF

      ! If energy, enstrophy or mixed advection of momentum in vector form change the value for masks
      ! at angles with three ocean points and one land point
      IF( ln_vorlat .AND. ( ln_dynvor_ene .OR. ln_dynvor_ens .OR. ln_dynvor_mix ) ) THEN
         DO jk = 1, jpk
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  IF( tmask(ji,jj,jk)+tmask(ji+1,jj,jk)+tmask(ji,jj+1,jk)+tmask(ji+1,jj+1,jk) == 3._wp ) &
                      fmask(ji,jj,jk) = 1._wp
               END DO
            END DO
         END DO
          !
          CALL lbc_lnk( fmask, 'F', 1._wp )      ! Lateral boundary conditions on fmask
          !
      ENDIF

      ioptio = 0                     ! Control of vorticity scheme options
      IF( ln_dynvor_ene )   ioptio = ioptio + 1
      IF( ln_dynvor_ens )   ioptio = ioptio + 1
      IF( ln_dynvor_mix )   ioptio = ioptio + 1
      IF( ln_dynvor_een )   ioptio = ioptio + 1
      IF( ln_dynvor_een_old )   ioptio = ioptio + 1
      IF( lk_esopa      )   ioptio =          1

      IF( ioptio /= 1 ) CALL ctl_stop( ' use ONE and ONLY one vorticity scheme' )

      !                              ! Set nvor (type of scheme for vorticity)
      IF( ln_dynvor_ene )   nvor =  0
      IF( ln_dynvor_ens )   nvor =  1
      IF( ln_dynvor_mix )   nvor =  2
      IF( ln_dynvor_een .or. ln_dynvor_een_old )   nvor =  3
      IF( lk_esopa      )   nvor = -1
      
      !                              ! Set ncor, nrvm, ntot (type of vorticity)
      IF(lwp) WRITE(numout,*)
      ncor = 1
      IF( ln_dynadv_vec ) THEN     
         IF(lwp) WRITE(numout,*) '         Vector form advection : vorticity = Coriolis + relative vorticity'
         nrvm = 2
         ntot = 4
      ELSE                        
         IF(lwp) WRITE(numout,*) '         Flux form advection   : vorticity = Coriolis + metric term'
         nrvm = 3
         ntot = 5
      ENDIF
      
      IF(lwp) THEN                   ! Print the choice
         WRITE(numout,*)
         IF( nvor ==  0 )   WRITE(numout,*) '         vorticity scheme : energy conserving scheme'
         IF( nvor ==  1 )   WRITE(numout,*) '         vorticity scheme : enstrophy conserving scheme'
         IF( nvor ==  2 )   WRITE(numout,*) '         vorticity scheme : mixed enstrophy/energy conserving scheme'
         IF( nvor ==  3 )   WRITE(numout,*) '         vorticity scheme : energy and enstrophy conserving scheme'
         IF( nvor == -1 )   WRITE(numout,*) '         esopa test: use all lateral physics options'
      ENDIF
      !
   END SUBROUTINE dyn_vor_init

   !!==============================================================================
END MODULE dynvor
