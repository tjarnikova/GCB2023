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

MODULE limdyn_2
   !!======================================================================
   !!                     ***  MODULE  limdyn_2  ***
   !!   Sea-Ice dynamics :  
   !!======================================================================
   !! History :  1.0  ! 2001-04  (LIM)  Original code
   !!            2.0  ! 2002-08  (C. Ethe, G. Madec)  F90, mpp
   !!            2.0  ! 2003-08  (C. Ethe) add lim_dyn_init
   !!            2.0  ! 2006-07  (G. Madec)  Surface module
   !!            3.3  ! 2009-05 (G. Garric, C. Bricaud) addition of the lim2_evp case
   !!---------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_lim2' :                                  LIM 2.0 sea-ice model
   !!----------------------------------------------------------------------
   !!    lim_dyn_2      : computes ice velocities
   !!    lim_dyn_init_2 : initialization and namelist read
   !!----------------------------------------------------------------------
   USE dom_oce          ! ocean space and time domain
   USE sbc_oce          ! ocean surface boundary condition
   USE phycst           ! physical constant
   USE ice_2            ! LIM-2: ice variables
   USE sbc_ice          ! Surface boundary condition: sea-ice fields
   USE dom_ice_2        ! LIM-2: ice domain
   USE limistate_2      ! LIM-2: initial state
   USE limrhg_2         ! LIM-2: VP  ice rheology
   USE limrhg           ! LIM  : EVP ice rheology
   USE lbclnk           ! lateral boundary condition - MPP link
   USE lib_mpp          ! MPP library
   USE wrk_nemo         ! work arrays
   USE in_out_manager   ! I/O manager
   USE prtctl           ! Print control
   USE lib_fortran      ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   lim_dyn_2   ! routine called by sbc_ice_lim

   !! * Substitutions
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
   !! NEMO/LIM2 3.3 , UCL - NEMO Consortium (2010)
   !! $Id: limdyn_2.F90 5123 2015-03-04 16:06:03Z clem $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE lim_dyn_2( kt )
      !!-------------------------------------------------------------------
      !!               ***  ROUTINE lim_dyn_2  ***
      !!               
      !! ** Purpose :   compute ice velocity and ocean-ice friction velocity
      !!                
      !! ** Method  : 
      !!
      !! ** Action  : - Initialisation
      !!              - Call of the dynamic routine for each hemisphere
      !!              - computation of the friction velocity at the sea-ice base
      !!              - treatment of the case if no ice dynamic
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! number of iteration
      !!
      INTEGER  ::   ji, jj             ! dummy loop indices
      INTEGER  ::   i_j1, i_jpj        ! Starting/ending j-indices for rheology
      REAL(wp) ::   zcoef              ! temporary scalar
      REAL(wp), POINTER, DIMENSION(:  ) ::   zind           ! i-averaged indicator of sea-ice
      REAL(wp), POINTER, DIMENSION(:  ) ::   zmsk           ! i-averaged of tmask
      REAL(wp), POINTER, DIMENSION(:,:) ::   zu_io, zv_io   ! ice-ocean velocity
      !!---------------------------------------------------------------------

      CALL wrk_alloc( jpi, jpj, zu_io, zv_io )
      CALL wrk_alloc(      jpj, zind , zmsk  )

      IF( kt == nit000 )   CALL lim_dyn_init_2   ! Initialization (first time-step only)
      
      IF( ln_limdyn ) THEN
         !
         ! Mean ice and snow thicknesses.          
         hsnm(:,:)  = ( 1.0 - frld(:,:) ) * hsnif(:,:)
         hicm(:,:)  = ( 1.0 - frld(:,:) ) * hicif(:,:)
         !
         !                                     ! Rheology (ice dynamics)
         !                                     ! ========
         
         !  Define the j-limits where ice rheology is computed
         ! ---------------------------------------------------
         
         IF( lk_mpp .OR. lk_mpp_rep ) THEN                    ! mpp: compute over the whole domain
            i_j1 = 1   
            i_jpj = jpj
            IF(ln_ctl)   CALL prt_ctl_info( 'lim_dyn  :    i_j1 = ', ivar1=i_j1, clinfo2=' ij_jpj = ', ivar2=i_jpj )
            IF( lk_lim2_vp )   THEN   ;   CALL lim_rhg_2( i_j1, i_jpj )             !  VP rheology
            ELSE                      ;   CALL lim_rhg  ( i_j1, i_jpj )             ! EVP rheology
            ENDIF
            !
         ELSE                                 ! optimization of the computational area
            !
            DO jj = 1, jpj
               zind(jj) = SUM( frld (:,jj  ) )   ! = REAL(jpj) if ocean everywhere on a j-line
               zmsk(jj) = SUM( tmask(:,jj,1) )   ! = 0         if land  everywhere on a j-line
            END DO
            !
            IF( l_jeq ) THEN                     ! local domain include both hemisphere
               !                                 ! Rheology is computed in each hemisphere
               !                                 ! only over the ice cover latitude strip
               ! Northern hemisphere
               i_j1  = njeq
               i_jpj = jpj
               DO WHILE ( i_j1 <= jpj .AND. zind(i_j1) == FLOAT(jpi) .AND. zmsk(i_j1) /=0 )
                  i_j1 = i_j1 + 1
               END DO
               IF( lk_lim2_vp )   THEN             ! VP  rheology
                  i_j1 = MAX( 1, i_j1-1 )
                  CALL lim_rhg_2( i_j1, i_jpj )
               ELSE                                ! EVP rheology
                  i_j1 = MAX( 1, i_j1-2 )
                  CALL lim_rhg( i_j1, i_jpj )
               ENDIF
               IF(ln_ctl)   WRITE(numout,*) 'lim_dyn : NH i_j1 = ', i_j1, 'ij_jpj = ', i_jpj
               !
               ! Southern hemisphere
               i_j1  =  1 
               i_jpj = njeq
               DO WHILE ( i_jpj >= 1 .AND. zind(i_jpj) == FLOAT(jpi) .AND. zmsk(i_jpj) /=0 )
                  i_jpj = i_jpj - 1
               END DO
               IF( lk_lim2_vp )   THEN             ! VP  rheology
                  i_jpj = MIN( jpj, i_jpj+2 )
                  CALL lim_rhg_2( i_j1, i_jpj )
               ELSE                                ! EVP rheology
                  i_jpj = MIN( jpj, i_jpj+1 )
                  CALL lim_rhg( i_j1, i_jpj )
               ENDIF
               IF(ln_ctl)   WRITE(numout,*) 'lim_dyn : SH i_j1 = ', i_j1, 'ij_jpj = ', i_jpj
               !
            ELSE                                 ! local domain extends over one hemisphere only
               !                                 ! Rheology is computed only over the ice cover
               !                                 ! latitude strip
               i_j1  = 1
               DO WHILE ( i_j1 <= jpj .AND. zind(i_j1) == FLOAT(jpi) .AND. zmsk(i_j1) /=0 )
                  i_j1 = i_j1 + 1
               END DO
               i_j1 = MAX( 1, i_j1-1 )
    
               i_jpj  = jpj
               DO WHILE ( i_jpj >= 1  .AND. zind(i_jpj) == FLOAT(jpi) .AND. zmsk(i_jpj) /=0 )
                  i_jpj = i_jpj - 1
               END DO
               i_jpj = MIN( jpj, i_jpj+2 )
               ! 
               IF( lk_lim2_vp )   THEN             ! VP  rheology
                  i_jpj = MIN( jpj, i_jpj+2 )
                  CALL lim_rhg_2( i_j1, i_jpj )                !  VP rheology
               ELSE                                ! EVP rheology
                  i_j1  = MAX( 1  , i_j1-2  )
                  i_jpj = MIN( jpj, i_jpj+1 )
                  CALL lim_rhg  ( i_j1, i_jpj )                ! EVP rheology
               ENDIF
               IF(ln_ctl)   WRITE(numout,*) 'lim_dyn : one hemisphere: i_j1 = ', i_j1, ' ij_jpj = ', i_jpj
               !
            ENDIF
            !
         ENDIF

         IF(ln_ctl)   CALL prt_ctl(tab2d_1=u_ice , clinfo1=' lim_dyn  : u_ice :', tab2d_2=v_ice , clinfo2=' v_ice :')
         
         ! computation of friction velocity
         ! --------------------------------
         SELECT CASE( cp_ice_msh )           ! ice-ocean relative velocity at u- & v-pts
         CASE( 'C' )                               ! EVP : C-grid ice dynamics
            zu_io(:,:) = u_ice(:,:) - ssu_m(:,:)           ! ice-ocean & ice velocity at ocean velocity points
            zv_io(:,:) = v_ice(:,:) - ssv_m(:,:)
         CASE( 'I' )                               ! VP  : B-grid ice dynamics (I-point) 
            DO jj = 1, jpjm1                               ! u_ice v_ice at I-point ; ssu_m, ssv_m at U- & V-points
               DO ji = 1, jpim1   ! NO vector opt.         !
                  zu_io(ji,jj) = 0.5_wp * ( u_ice(ji+1,jj+1) + u_ice(ji+1,jj  ) ) - ssu_m(ji,jj)
                  zv_io(ji,jj) = 0.5_wp * ( v_ice(ji+1,jj+1) + v_ice(ji  ,jj+1) ) - ssv_m(ji,jj)
               END DO
            END DO
         END SELECT

         ! frictional velocity at T-point
         zcoef = 0.5_wp * cw
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! NO vector opt. because of zu_io
               ust2s(ji,jj) = zcoef * (  zu_io(ji,jj) * zu_io(ji,jj) + zu_io(ji-1,jj) * zu_io(ji-1,jj)   &
                  &                    + zv_io(ji,jj) * zv_io(ji,jj) + zv_io(ji,jj-1) * zv_io(ji,jj-1)   ) * tms(ji,jj)
            END DO
         END DO
         !
      ELSE      ! no ice dynamics : transmit directly the atmospheric stress to the ocean
         !
         zcoef = SQRT( 0.5 ) / rau0
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               ust2s(ji,jj) = zcoef * SQRT(  utau(ji,jj) * utau(ji,jj) + utau(ji-1,jj) * utau(ji-1,jj)   &
                  &                        + vtau(ji,jj) * vtau(ji,jj) + vtau(ji,jj-1) * vtau(ji,jj-1)   ) * tms(ji,jj)
            END DO
         END DO
         !
      ENDIF
      !
      CALL lbc_lnk( ust2s, 'T',  1. )   ! T-point
      !
      IF(ln_ctl)   CALL prt_ctl(tab2d_1=ust2s , clinfo1=' lim_dyn  : ust2s :')
      !
      CALL wrk_dealloc( jpi, jpj, zu_io, zv_io )
      CALL wrk_dealloc(      jpj, zind , zmsk  )
      !
   END SUBROUTINE lim_dyn_2


   SUBROUTINE lim_dyn_init_2
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE lim_dyn_init_2  ***
      !!
      !! ** Purpose :   Physical constants and parameters linked to the ice
      !!              dynamics
      !!
      !! ** Method  :   Read the namicedyn namelist and check the ice-dynamic
      !!              parameter values
      !!
      !! ** input   :   Namelist namicedyn
      !!-------------------------------------------------------------------
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      NAMELIST/namicedyn/ epsd, alpha,     &
         &                dm, nbiter, nbitdr, om, resl, cw, angvg, pstar,   &
         &                c_rhg, etamn, rn_creepl, rn_ecc, ahi0,                  &
         &                nn_nevp, telast, alphaevp
      !!-------------------------------------------------------------------
                    
      REWIND( numnam_ice_ref )              ! Namelist namicedyn in reference namelist : Ice dynamics
      READ  ( numnam_ice_ref, namicedyn, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namicedyn in reference namelist', lwp )

      REWIND( numnam_ice_cfg )              ! Namelist namicedyn in configuration namelist : Ice dynamics
      READ  ( numnam_ice_cfg, namicedyn, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namicedyn in configuration namelist', lwp )
      IF(lwm) WRITE ( numoni, namicedyn )

      IF(lwp) THEN                                ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'lim_dyn_init_2: ice parameters for ice dynamics '
         WRITE(numout,*) '~~~~~~~~~~~~~~'
         WRITE(numout,*) '       tolerance parameter                              epsd   = ', epsd
         WRITE(numout,*) '       coefficient for semi-implicit coriolis           alpha  = ', alpha
         WRITE(numout,*) '       diffusion constant for dynamics                  dm     = ', dm
         WRITE(numout,*) '       number of sub-time steps for relaxation          nbiter = ', nbiter
         WRITE(numout,*) '       maximum number of iterations for relaxation      nbitdr = ', nbitdr
         WRITE(numout,*) '       relaxation constant                              om     = ', om
         WRITE(numout,*) '       maximum value for the residual of relaxation     resl   = ', resl
         WRITE(numout,*) '       drag coefficient for oceanic stress              cw     = ', cw
         WRITE(numout,*) '       turning angle for oceanic stress                 angvg  = ', angvg, ' degrees'
         WRITE(numout,*) '       first bulk-rheology parameter                    pstar  = ', pstar
         WRITE(numout,*) '       second bulk-rhelogy parameter                    c_rhg  = ', c_rhg
         WRITE(numout,*) '       minimun value for viscosity                      etamn  = ', etamn
         WRITE(numout,*) '       creep limit                                      rn_creepl = ', rn_creepl
         WRITE(numout,*) '       eccentricity of the elliptical yield curve       rn_ecc = ', rn_ecc
         WRITE(numout,*) '       horizontal diffusivity coeff. for sea-ice        ahi0   = ', ahi0
         WRITE(numout,*) '       number of iterations for subcycling              nn_nevp= ', nn_nevp
         WRITE(numout,*) '       timescale for elastic waves telast = ', telast
         WRITE(numout,*) '       coefficient for the solution of int. stresses alphaevp = ', alphaevp
      ENDIF
      !
      IF( angvg /= 0._wp .AND. .NOT.lk_lim2_vp ) THEN
         CALL ctl_warn( 'lim_dyn_init_2: turning angle for oceanic stress not properly coded for EVP ',   &
            &           '(see limsbc_2 module). We force  angvg = 0._wp'  )
         angvg = 0._wp
      ENDIF

      !  Initialization
      usecc2 = 1.0 / ( rn_ecc * rn_ecc )
      rhoco  = rau0 * cw
      angvg  = angvg * rad      ! convert angvg from degree to radian
      sangvg = SIN( angvg )
      cangvg = COS( angvg )
      pstarh = pstar / 2.0
      !
      ahiu(:,:) = ahi0 * umask(:,:,1)            ! Ice eddy Diffusivity coefficients.
      ahiv(:,:) = ahi0 * vmask(:,:,1)
      !
   END SUBROUTINE lim_dyn_init_2


   !!======================================================================
END MODULE limdyn_2
