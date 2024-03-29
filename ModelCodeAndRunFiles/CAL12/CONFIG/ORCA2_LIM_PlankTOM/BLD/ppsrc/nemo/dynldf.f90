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

MODULE dynldf
   !!======================================================================
   !!                       ***  MODULE  dynldf  ***
   !! Ocean physics:  lateral diffusivity trends 
   !!=====================================================================
   !! History :  9.0  !  05-11  (G. Madec)  Original code (new step architecture)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_ldf      : update the dynamics trend with the lateral diffusion
   !!   dyn_ldf_init : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE ldfdyn_oce     ! ocean dynamics lateral physics
   USE ldftra_oce     ! ocean tracers  lateral physics
   USE ldfslp         ! lateral mixing: slopes of mixing orientation
   USE dynldf_bilapg  ! lateral mixing            (dyn_ldf_bilapg routine)
   USE dynldf_bilap   ! lateral mixing            (dyn_ldf_bilap  routine)
   USE dynldf_iso     ! lateral mixing            (dyn_ldf_iso    routine)
   USE dynldf_lap     ! lateral mixing            (dyn_ldf_lap    routine)
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics   (trd_dyn        routine)
   !
   USE prtctl         ! Print control
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_ldf       ! called by step module 
   PUBLIC   dyn_ldf_init  ! called by opa  module 

   INTEGER ::   nldf = -2   ! type of lateral diffusion used defined from ln_dynldf_... namlist logicals)

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
   !! $Id: dynldf.F90 4990 2014-12-15 16:42:49Z timgraham $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_ldf( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_ldf  ***
      !! 
      !! ** Purpose :   compute the lateral ocean dynamics physics.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ztrdu, ztrdv
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_ldf')
      !
      IF( l_trddyn )   THEN                      ! temporary save of ta and sa trends
         CALL wrk_alloc( jpi, jpj, jpk, ztrdu, ztrdv )
         ztrdu(:,:,:) = ua(:,:,:) 
         ztrdv(:,:,:) = va(:,:,:) 
      ENDIF

      SELECT CASE ( nldf )                       ! compute lateral mixing trend and add it to the general trend
      !
      CASE ( 0 )    ;   CALL dyn_ldf_lap    ( kt )      ! iso-level laplacian
      CASE ( 1 )    ;   CALL dyn_ldf_iso    ( kt )      ! rotated laplacian (except dk[ dk[.] ] part)
      CASE ( 2 )    ;   CALL dyn_ldf_bilap  ( kt )      ! iso-level bilaplacian
      CASE ( 3 )    ;   CALL dyn_ldf_bilapg ( kt )      ! s-coord. horizontal bilaplacian
      CASE ( 4 )                                        ! iso-level laplacian + bilaplacian
         CALL dyn_ldf_lap    ( kt )
         CALL dyn_ldf_bilap  ( kt )
      CASE ( 5 )                                        ! rotated laplacian + bilaplacian (s-coord)
         CALL dyn_ldf_iso    ( kt )
         CALL dyn_ldf_bilapg ( kt )
      !
      CASE ( -1 )                                       ! esopa: test all possibility with control print
                        CALL dyn_ldf_lap    ( kt )
                        CALL prt_ctl( tab3d_1=ua, clinfo1=' ldf0 - Ua: ', mask1=umask,   &
            &                         tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
                        CALL dyn_ldf_iso    ( kt )
                        CALL prt_ctl( tab3d_1=ua, clinfo1=' ldf1 - Ua: ', mask1=umask,   &
            &                         tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
                        CALL dyn_ldf_bilap  ( kt )
                        CALL prt_ctl( tab3d_1=ua, clinfo1=' ldf2 - Ua: ', mask1=umask,   &
            &                         tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
                        CALL dyn_ldf_bilapg ( kt )
                        CALL prt_ctl( tab3d_1=ua, clinfo1=' ldf3 - Ua: ', mask1=umask,   &
            &                         tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      CASE ( -2 )                                       ! neither laplacian nor bilaplacian schemes used
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn_ldf : no lateral diffusion on momentum setup'
            IF(lwp) WRITE(numout,*) '~~~~~~~ '
         ENDIF
      END SELECT

      IF( l_trddyn ) THEN                        ! save the horizontal diffusive trends for further diagnostics
         ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
         ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_ldf, kt )
         CALL wrk_dealloc( jpi, jpj, jpk, ztrdu, ztrdv )
      ENDIF
      !                                          ! print sum trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' ldf  - Ua: ', mask1=umask,   &
         &                       tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_ldf')
      !
   END SUBROUTINE dyn_ldf


   SUBROUTINE dyn_ldf_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_ldf_init  ***
      !! 
      !! ** Purpose :   initializations of the horizontal ocean dynamics physics
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio, ierr         ! temporary integers 
      !!----------------------------------------------------------------------
    
      !                                   ! Namelist nam_dynldf: already read in ldfdyn module

      IF(lwp) THEN                        ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_ldf_init : Choice of the lateral diffusive operator on dynamics'
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*) '       Namelist nam_dynldf : set lateral mixing parameters (type, direction, coefficients)'
         WRITE(numout,*) '          laplacian operator          ln_dynldf_lap   = ', ln_dynldf_lap
         WRITE(numout,*) '          bilaplacian operator        ln_dynldf_bilap = ', ln_dynldf_bilap
         WRITE(numout,*) '          iso-level                   ln_dynldf_level = ', ln_dynldf_level
         WRITE(numout,*) '          horizontal (geopotential)   ln_dynldf_hor   = ', ln_dynldf_hor
         WRITE(numout,*) '          iso-neutral                 ln_dynldf_iso   = ', ln_dynldf_iso
      ENDIF

      !                                   ! control the consistency
      ioptio = 0
      IF( ln_dynldf_lap   )   ioptio = ioptio + 1
      IF( ln_dynldf_bilap )   ioptio = ioptio + 1
      IF( ioptio <  1 ) CALL ctl_warn( '          neither laplacian nor bilaplacian operator set for dynamics' )
      ioptio = 0
      IF( ln_dynldf_level )   ioptio = ioptio + 1
      IF( ln_dynldf_hor   )   ioptio = ioptio + 1
      IF( ln_dynldf_iso   )   ioptio = ioptio + 1
      IF( ioptio >  1 ) CALL ctl_stop( '          use only ONE direction (level/hor/iso)' )

      IF( ln_dynldf_iso .AND. ln_traldf_hor ) CALL ctl_stop &
      &   ( 'Not sensible to use geopotential diffusion for tracers with isoneutral diffusion for dynamics' )

      !                                   ! Set nldf, the type of lateral diffusion, from ln_dynldf_... logicals
      ierr = 0
      IF ( ln_dynldf_lap ) THEN      ! laplacian operator
         IF ( ln_zco ) THEN                ! z-coordinate
            IF ( ln_dynldf_level )   nldf = 0      ! iso-level  (no rotation)
            IF ( ln_dynldf_hor   )   nldf = 0      ! horizontal (no rotation)
            IF ( ln_dynldf_iso   )   nldf = 1      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_zps ) THEN             ! z-coordinate
            IF ( ln_dynldf_level )   ierr = 1      ! iso-level not allowed
            IF ( ln_dynldf_hor   )   nldf = 0      ! horizontal (no rotation)
            IF ( ln_dynldf_iso   )   nldf = 1      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_sco ) THEN             ! s-coordinate
            IF ( ln_dynldf_level )   nldf = 0      ! iso-level  (no rotation)
            IF ( ln_dynldf_hor   )   nldf = 1      ! horizontal (   rotation)
            IF ( ln_dynldf_iso   )   nldf = 1      ! isoneutral (   rotation)
         ENDIF
      ENDIF

      IF( ln_dynldf_bilap ) THEN      ! bilaplacian operator
         IF ( ln_zco ) THEN                ! z-coordinate
            IF ( ln_dynldf_level )   nldf = 2      ! iso-level  (no rotation)
            IF ( ln_dynldf_hor   )   nldf = 2      ! horizontal (no rotation)
            IF ( ln_dynldf_iso   )   ierr = 2      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_zps ) THEN             ! z-coordinate
            IF ( ln_dynldf_level )   ierr = 1      ! iso-level not allowed 
            IF ( ln_dynldf_hor   )   nldf = 2      ! horizontal (no rotation)
            IF ( ln_dynldf_iso   )   ierr = 2      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_sco ) THEN             ! s-coordinate
            IF ( ln_dynldf_level )   nldf = 2      ! iso-level  (no rotation)
            IF ( ln_dynldf_hor   )   nldf = 3      ! horizontal (   rotation)
            IF ( ln_dynldf_iso   )   ierr = 2      ! isoneutral (   rotation)
         ENDIF
      ENDIF
      
      IF( ln_dynldf_lap .AND. ln_dynldf_bilap ) THEN  ! mixed laplacian and bilaplacian operators
         IF ( ln_zco ) THEN                ! z-coordinate
            IF ( ln_dynldf_level )   nldf = 4      ! iso-level  (no rotation)
            IF ( ln_dynldf_hor   )   nldf = 4      ! horizontal (no rotation)
            IF ( ln_dynldf_iso   )   ierr = 2      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_zps ) THEN             ! z-coordinate
            IF ( ln_dynldf_level )   ierr = 1      ! iso-level not allowed 
            IF ( ln_dynldf_hor   )   nldf = 4      ! horizontal (no rotation)
            IF ( ln_dynldf_iso   )   ierr = 2      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_sco ) THEN             ! s-coordinate
            IF ( ln_dynldf_level )   nldf = 4      ! iso-level  (no rotation)
            IF ( ln_dynldf_hor   )   nldf = 5      ! horizontal (   rotation)
            IF ( ln_dynldf_iso   )   ierr = 2      ! isoneutral (   rotation)
         ENDIF
      ENDIF

      IF( lk_esopa )                 nldf = -1     ! esopa test

      IF( ierr == 1 )   CALL ctl_stop( 'iso-level in z-coordinate - partial step, not allowed' )
      IF( ierr == 2 )   CALL ctl_stop( 'isoneutral bilaplacian operator does not exist' )
      IF( nldf == 1 .OR. nldf == 3 ) THEN      ! rotation
         IF( .NOT.lk_ldfslp )   CALL ctl_stop( 'the rotation of the diffusive tensor require key_ldfslp' )
      ENDIF

      IF(lwp) THEN
         WRITE(numout,*)
         IF( nldf == -2 )   WRITE(numout,*) '              neither laplacian nor bilaplacian schemes used'
         IF( nldf == -1 )   WRITE(numout,*) '              ESOPA test All scheme used'
         IF( nldf ==  0 )   WRITE(numout,*) '              laplacian operator'
         IF( nldf ==  1 )   WRITE(numout,*) '              rotated laplacian operator'
         IF( nldf ==  2 )   WRITE(numout,*) '              bilaplacian operator'
         IF( nldf ==  3 )   WRITE(numout,*) '              rotated bilaplacian'
         IF( nldf ==  4 )   WRITE(numout,*) '              laplacian and bilaplacian operators'
         IF( nldf ==  5 )   WRITE(numout,*) '              rotated laplacian and bilaplacian operators'
      ENDIF
      !
   END SUBROUTINE dyn_ldf_init

   !!======================================================================
END MODULE dynldf
