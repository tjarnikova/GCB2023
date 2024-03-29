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

MODULE sbcfwb
   !!======================================================================
   !!                       ***  MODULE  sbcfwb  ***
   !! Ocean fluxes   : domain averaged freshwater budget
   !!======================================================================
   !! History :  OPA  ! 2001-02  (E. Durand)  Original code
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module
   !!            3.0  ! 2006-08  (G. Madec)  Surface module
   !!            3.2  ! 2009-07  (C. Talandier) emp mean s spread over erp area 
   !!            3.6  ! 2014-11  (P. Mathiot  ) add ice shelf melting
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_fwb      : freshwater budget for global ocean configurations
   !!                  in free surface and forced mode
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! surface ocean boundary condition
   USE phycst          ! physical constants
   USE sbcrnf          ! ocean runoffs
   USE sbcisf          ! ice shelf melting contribution
   USE sbcssr          ! SS damping terms
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
   USE wrk_nemo        ! work arrays
   USE timing          ! Timing
   USE lbclnk          ! ocean lateral boundary conditions
   USE lib_fortran

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_fwb    ! routine called by step

   REAL(wp) ::   a_fwb_b   ! annual domain averaged freshwater budget
   REAL(wp) ::   a_fwb     ! for 2 year before (_b) and before year.
   REAL(wp) ::   fwfold    ! fwfold to be suppressed
   REAL(wp) ::   area      ! global mean ocean surface (interior domain)

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
   !! $Id: sbcfwb.F90 5120 2015-03-03 16:11:55Z acc $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_fwb( kt, kn_fwb, kn_fsbc )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_fwb  ***
      !!
      !! ** Purpose :   Control the mean sea surface drift
      !!
      !! ** Method  :   several ways  depending on kn_fwb
      !!                =0 no control 
      !!                =1 global mean of emp set to zero at each nn_fsbc time step
      !!                =2 annual global mean corrected from previous year
      !!                =3 global mean of emp set to zero at each nn_fsbc time step
      !!                   & spread out over erp area depending its sign
      !! Note: if sea ice is embedded it is taken into account when computing the budget 
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
      INTEGER, INTENT( in ) ::   kn_fsbc  ! 
      INTEGER, INTENT( in ) ::   kn_fwb   ! ocean time-step index
      !
      INTEGER  ::   inum, ikty, iyear     ! local integers
      REAL(wp) ::   z_fwf, z_fwf_nsrf, zsum_fwf, zsum_erp                ! local scalars
      REAL(wp) ::   zsurf_neg, zsurf_pos, zsurf_tospread, zcoef          !   -      -
      REAL(wp), POINTER, DIMENSION(:,:) ::   ztmsk_neg, ztmsk_pos, z_wgt ! 2D workspaces
      REAL(wp), POINTER, DIMENSION(:,:) ::   ztmsk_tospread, zerp_cor    !   -      -
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('sbc_fwb')
      !
      CALL wrk_alloc( jpi,jpj, ztmsk_neg, ztmsk_pos, ztmsk_tospread, z_wgt, zerp_cor )
      !
      IF( kt == nit000 ) THEN
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'sbc_fwb : FreshWater Budget correction'
            WRITE(numout,*) '~~~~~~~'
            IF( kn_fwb == 1 )   WRITE(numout,*) '          instantaneously set to zero'
            IF( kn_fwb == 2 )   WRITE(numout,*) '          adjusted from previous year budget'
            IF( kn_fwb == 3 )   WRITE(numout,*) '          fwf set to zero and spread out over erp area'
         ENDIF
         !
         IF( kn_fwb == 3 .AND. nn_sssr /= 2 )   CALL ctl_stop( 'sbc_fwb: nn_fwb = 3 requires nn_sssr = 2, we stop ' )
         IF( kn_fwb == 3 .AND. ln_isfcav    )   CALL ctl_stop( 'sbc_fwb: nn_fwb = 3 with ln_isfcav = .TRUE. not working, we stop ' )
         !
         area = glob_sum( e1e2t(:,:) * tmask(:,:,1))           ! interior global domain surface
         ! isf cavities are excluded because it can feedback to the melting with generation of inhibition of plumes
         ! and in case of no melt, it can generate HSSW.
         !
         !
      ENDIF
      

      SELECT CASE ( kn_fwb )
      !
      CASE ( 1 )                             !==  global mean fwf set to zero  ==!
         !
         IF( MOD( kt-1, kn_fsbc ) == 0 ) THEN
            z_fwf = glob_sum( e1e2t(:,:) * ( emp(:,:) - rnf(:,:) + rdivisf * fwfisf(:,:) -  snwice_fmass(:,:) ) ) / area   ! sum over the global domain
            zcoef = z_fwf * rcp
            emp(:,:) = emp(:,:) - z_fwf              * tmask(:,:,1)
            qns(:,:) = qns(:,:) + zcoef * sst_m(:,:) * tmask(:,:,1) ! account for change to the heat budget due to fw correction
         ENDIF
         !
      CASE ( 2 )                             !==  fwf budget adjusted from the previous year  ==!
         !
         IF( kt == nit000 ) THEN                      ! initialisation
            !                                         ! Read the corrective factor on precipitations (fwfold)
            CALL ctl_opn( inum, 'EMPave_old.dat', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
            READ ( inum, "(24X,I8,2ES24.16)" ) iyear, a_fwb_b, a_fwb
            CLOSE( inum )
            fwfold = a_fwb                            ! current year freshwater budget correction
            !                                         ! estimate from the previous year budget
            IF(lwp)WRITE(numout,*)
            IF(lwp)WRITE(numout,*)'sbc_fwb : year = ',iyear  , ' freshwater budget correction = ', fwfold
            IF(lwp)WRITE(numout,*)'          year = ',iyear-1, ' freshwater budget read       = ', a_fwb
            IF(lwp)WRITE(numout,*)'          year = ',iyear-2, ' freshwater budget read       = ', a_fwb_b
         ENDIF   
         !                                         ! Update fwfold if new year start
         ikty = 365 * 86400 / rdttra(1)    !!bug  use of 365 days leap year or 360d year !!!!!!!
         IF( MOD( kt, ikty ) == 0 ) THEN
            a_fwb_b = a_fwb                           ! mean sea level taking into account the ice+snow
                                                      ! sum over the global domain
            a_fwb   = glob_sum( e1e2t(:,:) * ( sshn(:,:) + snwice_mass(:,:) * r1_rau0 ) )
            a_fwb   = a_fwb * 1.e+3 / ( area * rday * 365. )     ! convert in Kg/m3/s = mm/s
!!gm        !                                                      !!bug 365d year 
            fwfold =  a_fwb                           ! current year freshwater budget correction
            !                                         ! estimate from the previous year budget
         ENDIF
         ! 
         IF( MOD( kt-1, kn_fsbc ) == 0 ) THEN         ! correct the freshwater fluxes
            zcoef = fwfold * rcp
            emp(:,:) = emp(:,:) + fwfold             * tmask(:,:,1)
            qns(:,:) = qns(:,:) - zcoef * sst_m(:,:) * tmask(:,:,1) ! account for change to the heat budget due to fw correction
         ENDIF
         !
         IF( kt == nitend .AND. lwp ) THEN            ! save fwfold value in a file
            CALL ctl_opn( inum, 'EMPave.dat', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE., narea )
            WRITE( inum, "(24X,I8,2ES24.16)" ) nyear, a_fwb_b, a_fwb
            CLOSE( inum )
         ENDIF
         !
      CASE ( 3 )                             !==  global fwf set to zero and spread out over erp area  ==!
         !
         IF( MOD( kt-1, kn_fsbc ) == 0 ) THEN
            ztmsk_pos(:,:) = tmask_i(:,:)                      ! Select <0 and >0 area of erp
            WHERE( erp < 0._wp )   ztmsk_pos = 0._wp
            ztmsk_neg(:,:) = tmask_i(:,:) - ztmsk_pos(:,:)
            !
            zsurf_neg = glob_sum( e1e2t(:,:)*ztmsk_neg(:,:) )  ! Area filled by <0 and >0 erp 
            zsurf_pos = glob_sum( e1e2t(:,:)*ztmsk_pos(:,:) )
            !                                                  ! fwf global mean (excluding ocean to ice/snow exchanges) 
            z_fwf     = glob_sum( e1e2t(:,:) * ( emp(:,:) - rnf(:,:) + rdivisf * fwfisf(:,:) - snwice_fmass(:,:) ) ) / area
            !            
            IF( z_fwf < 0._wp ) THEN         ! spread out over >0 erp area to increase evaporation
                zsurf_tospread      = zsurf_pos
                ztmsk_tospread(:,:) = ztmsk_pos(:,:)
            ELSE                             ! spread out over <0 erp area to increase precipitation
                zsurf_tospread      = zsurf_neg
                ztmsk_tospread(:,:) = ztmsk_neg(:,:)
            ENDIF
            !
            zsum_fwf   = glob_sum( e1e2t(:,:) * z_fwf )         ! fwf global mean over <0 or >0 erp area
!!gm :  zsum_fwf   = z_fwf * area   ???  it is right?  I think so....
            z_fwf_nsrf =  zsum_fwf / ( zsurf_tospread + rsmall )
            !                                                  ! weight to respect erp field 2D structure 
            zsum_erp   = glob_sum( ztmsk_tospread(:,:) * erp(:,:) * e1e2t(:,:) )
            z_wgt(:,:) = ztmsk_tospread(:,:) * erp(:,:) / ( zsum_erp + rsmall )
            !                                                  ! final correction term to apply
            zerp_cor(:,:) = -1. * z_fwf_nsrf * zsurf_tospread * z_wgt(:,:)
            !
!!gm   ===>>>>  lbc_lnk should be useless as all the computation is done over the whole domain !
            CALL lbc_lnk( zerp_cor, 'T', 1. )
            !
            emp(:,:) = emp(:,:) + zerp_cor(:,:)
            qns(:,:) = qns(:,:) - zerp_cor(:,:) * rcp * sst_m(:,:)  ! account for change to the heat budget due to fw correction
            erp(:,:) = erp(:,:) + zerp_cor(:,:)
            !
            IF( nprint == 1 .AND. lwp ) THEN                   ! control print
               IF( z_fwf < 0._wp ) THEN
                  WRITE(numout,*)'   z_fwf < 0'
                  WRITE(numout,*)'   SUM(erp+)     = ', SUM( ztmsk_tospread(:,:)*erp(:,:)*e1e2t(:,:) )*1.e-9,' Sv'
               ELSE
                  WRITE(numout,*)'   z_fwf >= 0'
                  WRITE(numout,*)'   SUM(erp-)     = ', SUM( ztmsk_tospread(:,:)*erp(:,:)*e1e2t(:,:) )*1.e-9,' Sv'
               ENDIF
               WRITE(numout,*)'   SUM(empG)     = ', SUM( z_fwf*e1e2t(:,:) )*1.e-9,' Sv'
               WRITE(numout,*)'   z_fwf         = ', z_fwf      ,' Kg/m2/s'
               WRITE(numout,*)'   z_fwf_nsrf    = ', z_fwf_nsrf ,' Kg/m2/s'
               WRITE(numout,*)'   MIN(zerp_cor) = ', MINVAL(zerp_cor) 
               WRITE(numout,*)'   MAX(zerp_cor) = ', MAXVAL(zerp_cor) 
            ENDIF
         ENDIF
         !
      CASE DEFAULT                           !==  you should never be there  ==!
         CALL ctl_stop( 'sbc_fwb : wrong nn_fwb value for the FreshWater Budget correction, choose either 1, 2 or 3' )
         !
      END SELECT
      !
      CALL wrk_dealloc( jpi,jpj, ztmsk_neg, ztmsk_pos, ztmsk_tospread, z_wgt, zerp_cor )
      !
      IF( nn_timing == 1 )  CALL timing_stop('sbc_fwb')
      !
   END SUBROUTINE sbc_fwb

   !!======================================================================
END MODULE sbcfwb
