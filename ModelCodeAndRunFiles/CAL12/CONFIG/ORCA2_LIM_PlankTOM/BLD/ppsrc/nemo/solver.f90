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

MODULE solver
   !!======================================================================
   !!                     ***  MODULE  solver  ***
   !! Ocean solver :  initialization of ocean solver
   !!=====================================================================
   !! History :  OPA  ! 1990-10  (G. Madec)  Original code           
   !!                 ! 1993-02  (O. Marti)                         
   !!                 ! 1997-02  (G. Madec)  local depth inverse computation
   !!                 ! 1998-10  (G. Roullet, G. Madec)  free surface 
   !!   NEMO     1.0  ! 2003-07  (G. Madec)  free form, F90
   !!            3.2  ! 2009-07  (R. Benshila) suppression of rigid-lid & FETI solver
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_dynspg_flt'                              filtered free surface
   !!----------------------------------------------------------------------
   !!   solver_init: solver initialization
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE zdf_oce         ! ocean vertical physics variables
   USE sol_oce         ! solver variables
   USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
   USE solmat          ! matrix of the solver
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! MPP library
   USE timing          ! timing

   IMPLICIT NONE

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: solver.F90 5533 2015-07-02 13:50:50Z jchanut $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE solver_init( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE solver_init  ***
      !!                   
      !! ** Purpose :   Initialization of the elliptic solver
      !!      
      !! ** Method  :   a solver is required when using the filtered free
      !!              surface. 
      !!
      !! ** Action  : - c_solver_pt : nature of the gridpoint at which the solver is applied
      !!
      !! References : Jensen, 1986: Adv. Phys. Oceanogr. Num. Mod.,Ed. O Brien,87-110.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt
      INTEGER             ::   ios       ! Local integer output status for namelist read
      !!
      NAMELIST/namsol/ nn_solv, nn_sol_arp, nn_nmin, nn_nmax, nn_nmod, rn_eps, rn_resmax, rn_sor
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('solver_init')
      !

      IF(lwp) THEN                  !* open elliptic solver statistics file (only on the printing processors)
         CALL ctl_opn( numsol, 'solver.stat', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
      ENDIF

      REWIND( numnam_ref )              ! Namelist namsol in reference namelist : Elliptic solver / free surface
      READ  ( numnam_ref, namsol, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsol in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namsol in configuration namelist : Elliptic solver / free surface
      READ  ( numnam_cfg, namsol, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsol in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namsol )

      IF(lwp) THEN                  !* Control print
         WRITE(numout,*)
         WRITE(numout,*) 'solver_init : solver to compute the surface pressure gradient'
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namsol : set solver parameters'
         WRITE(numout,*) '      type of elliptic solver            nn_solv    = ', nn_solv
         WRITE(numout,*) '      absolute/relative (0/1) precision  nn_sol_arp = ', nn_sol_arp
         WRITE(numout,*) '      minimum iterations for solver      nn_nmin    = ', nn_nmin
         WRITE(numout,*) '      maximum iterations for solver      nn_nmax    = ', nn_nmax
         WRITE(numout,*) '      frequency for test                 nn_nmod    = ', nn_nmod
         WRITE(numout,*) '      absolute precision of solver       rn_eps     = ', rn_eps
         WRITE(numout,*) '      absolute precision for SOR solver  rn_resmax  = ', rn_resmax
         WRITE(numout,*) '      optimal coefficient of sor         rn_sor     = ', rn_sor
         WRITE(numout,*)
      ENDIF
      eps = rn_eps

      !                              ! allocate solver arrays
      IF( .NOT. lk_agrif .OR. .NOT. ln_rstart) THEN
         IF( sol_oce_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'solver_init : unable to allocate sol_oce arrays' )
         gcx (:,:) = 0.e0
         gcxb(:,:) = 0.e0
      ENDIF

      SELECT CASE( nn_solv )          !* parameter check
      !
      CASE ( 1 )                          ! preconditioned conjugate gradient solver
         IF(lwp) WRITE(numout,*) '   a preconditioned conjugate gradient solver is used'
         IF( jpr2di /= 0 .AND. jpr2dj /= 0 )   CALL ctl_stop( ' jpr2di and jpr2dj should be equal to zero' )
         !
      CASE ( 2 )                          ! successive-over-relaxation solver
         IF(lwp) WRITE(numout,*) '   a successive-over-relaxation solver with extra outer halo is used'
         IF(lwp) WRITE(numout,*) '   with jpr2di =', jpr2di, ' and  jpr2dj =', jpr2dj
         IF( .NOT. lk_mpp .AND. jpr2di /= 0 .AND. jpr2dj /= 0 ) THEN
             CALL ctl_stop( 'jpr2di and jpr2dj are not equal to zero',   &
             &              'In this case the algorithm should be used only with the key_mpp_... option' )
         ELSE
            IF( ( ( jperio == 1 .OR. jperio == 4 .OR. jperio == 6 ) .OR. ( jpni /= 1 ) ) &
              &  .AND. ( jpr2di /= jpr2dj ) )   CALL ctl_stop( 'jpr2di should be equal to jpr2dj' )
         ENDIF
         !
      CASE DEFAULT                        ! error in parameter
         WRITE(ctmp1,*) '          bad flag value for nn_solv = ', nn_solv
         CALL ctl_stop( ctmp1 )
      END SELECT
      !

      !                             !* Grid-point at which the solver is applied
!!gm  c_solver_pt should be removed: nomore bsf, only T-point is used
      c_solver_pt = 'T'                   ! always T-point (ssh solver only, not anymore bsf)

      CALL sol_mat( kt )            !* Construction of the elliptic system matrix
      !
      IF( nn_timing == 1 )  CALL timing_stop('solver_init')
      !
   END SUBROUTINE solver_init

   !!======================================================================
END MODULE solver
