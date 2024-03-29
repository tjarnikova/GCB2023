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

MODULE trazdf
   !!==============================================================================
   !!                 ***  MODULE  trazdf  ***
   !! Ocean active tracers:  vertical component of the tracer mixing trend
   !!==============================================================================
   !! History :  1.0  ! 2005-11  (G. Madec)  Original code
   !!            3.0  ! 2008-01  (C. Ethe, G. Madec)  merge TRC-TRA
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_zdf      : Update the tracer trend with the vertical diffusion
   !!   tra_zdf_init : initialisation of the computation
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE domvvl          ! variable volume
   USE phycst          ! physical constant
   USE zdf_oce         ! ocean vertical physics variables
   USE sbc_oce         ! surface boundary condition: ocean
   USE dynspg_oce
   USE trazdf_exp      ! vertical diffusion: explicit (tra_zdf_exp     routine)
   USE trazdf_imp      ! vertical diffusion: implicit (tra_zdf_imp     routine)
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE trd_oce         ! trends: ocean variables
   USE trdtra          ! trends manager: tracers 
   !
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! Memory allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_zdf        ! routine called by step.F90
   PUBLIC   tra_zdf_init   ! routine called by nemogcm.F90

   INTEGER ::   nzdf = 0   ! type vertical diffusion algorithm used (defined from ln_zdf...  namlist logicals)

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
   !! NEMO/OPA 3.7 , NEMO Consortium (2014)
   !! $Id: trazdf.F90 5385 2015-06-09 13:50:42Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_zdf( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_zdf  ***
      !!
      !! ** Purpose :   compute the vertical ocean tracer physics.
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !!
      INTEGER  ::   jk                   ! Dummy loop indices
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   ztrdt, ztrds   ! 3D workspace
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_zdf')
      !
      IF( neuler == 0 .AND. kt == nit000 ) THEN     ! at nit000
         r2dtra(:) =  rdttra(:)                          ! = rdtra (restarting with Euler time stepping)
      ELSEIF( kt <= nit000 + 1) THEN                ! at nit000 or nit000+1
         r2dtra(:) = 2. * rdttra(:)                      ! = 2 rdttra (leapfrog)
      ENDIF

      IF( l_trdtra )   THEN                    !* Save ta and sa trends
         CALL wrk_alloc( jpi, jpj, jpk, ztrdt, ztrds )
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem)
         ztrds(:,:,:) = tsa(:,:,:,jp_sal)
      ENDIF

      SELECT CASE ( nzdf )                       ! compute lateral mixing trend and add it to the general trend
      CASE ( 0 )    ;    CALL tra_zdf_exp( kt, nit000, 'TRA', r2dtra, nn_zdfexp, tsb, tsa, jpts )  !   explicit scheme 
      CASE ( 1 )    ;    CALL tra_zdf_imp( kt, nit000, 'TRA', r2dtra,            tsb, tsa, jpts )  !   implicit scheme 
      CASE ( -1 )                                       ! esopa: test all possibility with control print
         CALL tra_zdf_exp( kt, nit000, 'TRA', r2dtra, nn_zdfexp, tsb, tsa, jpts )
         CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' zdf0 - Ta: ', mask1=tmask,               &
         &             tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
         CALL tra_zdf_imp( kt, nit000, 'TRA', r2dtra,            tsb, tsa, jpts ) 
         CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' zdf1 - Ta: ', mask1=tmask,               &
         &             tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
      END SELECT
      ! DRAKKAR SSS control {
      ! JMM avoid negative salinities near river outlet ! Ugly fix
      ! JMM : restore negative salinities to small salinities:
      WHERE ( tsa(:,:,:,jp_sal) < 0._wp )   tsa(:,:,:,jp_sal) = 0.1_wp

      IF( l_trdtra )   THEN                      ! save the vertical diffusive trends for further diagnostics
         DO jk = 1, jpkm1
            ztrdt(:,:,jk) = ( ( tsa(:,:,jk,jp_tem) - tsb(:,:,jk,jp_tem) ) / r2dtra(jk) ) - ztrdt(:,:,jk)
            ztrds(:,:,jk) = ( ( tsa(:,:,jk,jp_sal) - tsb(:,:,jk,jp_sal) ) / r2dtra(jk) ) - ztrds(:,:,jk)
         END DO
         CALL lbc_lnk( ztrdt, 'T', 1. )
         CALL lbc_lnk( ztrds, 'T', 1. )
         CALL trd_tra( kt, 'TRA', jp_tem, jptra_zdf, ztrdt )
         CALL trd_tra( kt, 'TRA', jp_sal, jptra_zdf, ztrds )
         CALL wrk_dealloc( jpi, jpj, jpk, ztrdt, ztrds )
      ENDIF

      !                                          ! print mean trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' zdf  - Ta: ', mask1=tmask,               &
         &                       tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_zdf')
      !
   END SUBROUTINE tra_zdf


   SUBROUTINE tra_zdf_init
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE tra_zdf_init  ***
      !!
      !! ** Purpose :   Choose the vertical mixing scheme
      !!
      !! ** Method  :   Set nzdf from ln_zdfexp
      !!      nzdf = 0   explicit (time-splitting) scheme (ln_zdfexp=T)
      !!           = 1   implicit (euler backward) scheme (ln_zdfexp=F)
      !!      NB: rotation of lateral mixing operator or TKE or KPP scheme,
      !!      the implicit scheme is required.
      !!----------------------------------------------------------------------
      USE zdftke
      USE zdfgls
      USE zdfkpp
      !!----------------------------------------------------------------------

      ! Choice from ln_zdfexp already read in namelist in zdfini module
      IF( ln_zdfexp ) THEN   ;   nzdf = 0           ! use explicit scheme
      ELSE                   ;   nzdf = 1           ! use implicit scheme
      ENDIF

      ! Force implicit schemes
      IF( lk_zdftke .OR. lk_zdfgls .OR. lk_zdfkpp )   nzdf = 1      ! TKE, GLS or KPP physics
      IF( ln_traldf_iso                           )   nzdf = 1      ! iso-neutral lateral physics
      IF( ln_traldf_hor .AND. ln_sco              )   nzdf = 1      ! horizontal lateral physics in s-coordinate
      IF( ln_zdfexp .AND. nzdf == 1 )   CALL ctl_stop( 'tra_zdf : If using the rotation of lateral mixing operator',   &
            &                         ' TKE or KPP scheme, the implicit scheme is required, set ln_zdfexp = .false.' )

      ! Test: esopa
      IF( lk_esopa )    nzdf = -1                      ! All schemes used

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tra_zdf_init : vertical tracer physics scheme'
         WRITE(numout,*) '~~~~~~~~~~~'
         IF( nzdf == -1 )   WRITE(numout,*) '              ESOPA test All scheme used'
         IF( nzdf ==  0 )   WRITE(numout,*) '              Explicit time-splitting scheme'
         IF( nzdf ==  1 )   WRITE(numout,*) '              Implicit (euler backward) scheme'
      ENDIF
      !
   END SUBROUTINE tra_zdf_init

   !!==============================================================================
END MODULE trazdf
