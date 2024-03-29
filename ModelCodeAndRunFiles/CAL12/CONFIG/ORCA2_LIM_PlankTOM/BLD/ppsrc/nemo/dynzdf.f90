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

MODULE dynzdf
   !!==============================================================================
   !!                 ***  MODULE  dynzdf  ***
   !! Ocean dynamics :  vertical component of the momentum mixing trend
   !!==============================================================================
   !! History :  1.0  !  2005-11  (G. Madec)  Original code
   !!            3.3  !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_zdf      : Update the momentum trend with the vertical diffusion
   !!   dyn_zdf_init : initializations of the vertical diffusion scheme
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE zdf_oce         ! ocean vertical physics variables

   USE dynzdf_exp      ! vertical diffusion: explicit (dyn_zdf_exp     routine)
   USE dynzdf_imp      ! vertical diffusion: implicit (dyn_zdf_imp     routine)

   USE ldfdyn_oce      ! ocean dynamics: lateral physics
   USE trd_oce         ! trends: ocean variables
   USE trddyn          ! trend manager: dynamics
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   USE prtctl          ! Print control
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_zdf       !  routine called by step.F90
   PUBLIC   dyn_zdf_init  !  routine called by opa.F90

   INTEGER  ::   nzdf = 0   ! type vertical diffusion algorithm used, defined from ln_zdf... namlist logicals
   REAL(wp) ::   r2dt       ! time-step, = 2 rdttra except at nit000 (=rdttra) if neuler=0

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
   !!                    *** zdfddm_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
   !!      with a constant or 1D or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
!   'key_zdfddm' :                      avs: 3D array defined in zdfddm module
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: zdfddm_substitute.h90 4152 2013-11-05 11:59:53Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
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
   !! $Id: dynzdf.F90 4990 2014-12-15 16:42:49Z timgraham $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS
   
   SUBROUTINE dyn_zdf( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_zdf  ***
      !!
      !! ** Purpose :   compute the vertical ocean dynamics physics.
      !!---------------------------------------------------------------------
      !!
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ztrdu, ztrdv
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_zdf')
      !
      !                                          ! set time step
      IF( neuler == 0 .AND. kt == nit000     ) THEN   ;   r2dt =      rdt   ! = rdtra (restart with Euler time stepping)
      ELSEIF(               kt <= nit000 + 1 ) THEN   ;   r2dt = 2. * rdt   ! = 2 rdttra (leapfrog)
      ENDIF

      IF( l_trddyn )   THEN                      ! temporary save of ta and sa trends
         CALL wrk_alloc( jpi, jpj, jpk, ztrdu, ztrdv ) 
         ztrdu(:,:,:) = ua(:,:,:)
         ztrdv(:,:,:) = va(:,:,:)
      ENDIF

      SELECT CASE ( nzdf )                       ! compute lateral mixing trend and add it to the general trend
      !
      CASE ( 0 )   ;   CALL dyn_zdf_exp( kt, r2dt )      ! explicit scheme
      CASE ( 1 )   ;   CALL dyn_zdf_imp( kt, r2dt )      ! implicit scheme
      !
      CASE ( -1 )                                        ! esopa: test all possibility with control print
                       CALL dyn_zdf_exp( kt, r2dt )
                       CALL prt_ctl( tab3d_1=ua, clinfo1=' zdf0 - Ua: ', mask1=umask,               &
                          &          tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
                       CALL dyn_zdf_imp( kt, r2dt )
                       CALL prt_ctl( tab3d_1=ua, clinfo1=' zdf1 - Ua: ', mask1=umask,               &
                          &          tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      END SELECT

      IF( l_trddyn )   THEN                      ! save the vertical diffusive trends for further diagnostics
         ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
         ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_zdf, kt )
         CALL wrk_dealloc( jpi, jpj, jpk, ztrdu, ztrdv ) 
      ENDIF
      !                                          ! print mean trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' zdf  - Ua: ', mask1=umask,               &
            &                    tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_zdf')
      !
   END SUBROUTINE dyn_zdf


   SUBROUTINE dyn_zdf_init
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dyn_zdf_init  ***
      !!
      !! ** Purpose :   initializations of the vertical diffusion scheme
      !!
      !! ** Method  :   implicit (euler backward) scheme (default)
      !!                explicit (time-splitting) scheme if ln_zdfexp=T
      !!----------------------------------------------------------------------
      USE zdftke
      USE zdfgls
      USE zdfkpp
      !!----------------------------------------------------------------------
      !
      ! Choice from ln_zdfexp read in namelist in zdfini
      IF( ln_zdfexp ) THEN   ;   nzdf = 0           ! use explicit scheme
      ELSE                   ;   nzdf = 1           ! use implicit scheme
      ENDIF
      !
      ! Force implicit schemes
      IF( lk_zdftke .OR. lk_zdfgls .OR. lk_zdfkpp )   nzdf = 1   ! TKE, GLS or KPP physics
      IF( ln_dynldf_iso                           )   nzdf = 1   ! iso-neutral lateral physics
      IF( ln_dynldf_hor .AND. ln_sco              )   nzdf = 1   ! horizontal lateral physics in s-coordinate
      !
      IF( lk_esopa )    nzdf = -1                   ! Esopa key: All schemes used
      !
      IF(lwp) THEN                                  ! Print the choice
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_zdf_init : vertical dynamics physics scheme'
         WRITE(numout,*) '~~~~~~~~~~~'
         IF( nzdf == -1 )   WRITE(numout,*) '              ESOPA test All scheme used'
         IF( nzdf ==  0 )   WRITE(numout,*) '              Explicit time-splitting scheme'
         IF( nzdf ==  1 )   WRITE(numout,*) '              Implicit (euler backward) scheme'
      ENDIF
      !
   END SUBROUTINE dyn_zdf_init

   !!==============================================================================
END MODULE dynzdf
