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

MODULE trabbc
   !!==============================================================================
   !!                       ***  MODULE  trabbc  ***
   !! Ocean active tracers:  bottom boundary condition (geothermal heat flux)
   !!==============================================================================
   !! History :  OPA  ! 1999-10 (G. Madec)  original code
   !!   NEMO     1.0  ! 2002-08 (G. Madec)  free form + modules
   !!             -   ! 2002-11 (A. Bozec)  tra_bbc_init: original code
   !!            3.3  ! 2010-10 (G. Madec)  dynamical allocation + suppression of key_trabbc
   !!             -   ! 2010-11 (G. Madec)  use mbkt array (deepest ocean t-level)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_bbc      : update the tracer trend at ocean bottom 
   !!   tra_bbc_init : initialization of geothermal heat flux trend
   !!----------------------------------------------------------------------
   USE oce             ! ocean variables
   USE dom_oce         ! domain: ocean
   USE phycst          ! physical constants
   USE trd_oce         ! trends: ocean variables
   USE trdtra          ! trends manager: tracers 
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O manager
   USE fldread         ! read input fields
   USE lbclnk            ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp           ! distributed memory computing library
   USE prtctl          ! Print control
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC tra_bbc          ! routine called by step.F90
   PUBLIC tra_bbc_init     ! routine called by opa.F90

   !                                 !!* Namelist nambbc: bottom boundary condition *
   LOGICAL, PUBLIC ::   ln_trabbc     !: Geothermal heat flux flag
   INTEGER         ::   nn_geoflx     !  Geothermal flux (=1:constant flux, =2:read in file )
   REAL(wp)        ::   rn_geoflx_cst !  Constant value of geothermal heat flux

   REAL(wp), PUBLIC, DIMENSION(:,:), ALLOCATABLE ::   qgh_trd0   ! geothermal heating trend
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_qgh              ! structure of input qgh (file informations, fields read)
 
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
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: trabbc.F90 5397 2015-06-10 14:00:44Z cbricaud $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_bbc( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_bbc  ***
      !!
      !! ** Purpose :   Compute the bottom boundary contition on temperature 
      !!              associated with geothermal heating and add it to the 
      !!              general trend of temperature equations.
      !!
      !! ** Method  :   The geothermal heat flux set to its constant value of 
      !!              86.4 mW/m2 (Stein and Stein 1992, Huang 1999).
      !!       The temperature trend associated to this heat flux through the
      !!       ocean bottom can be computed once and is added to the temperature
      !!       trend juste above the bottom at each time step:
      !!            ta = ta + Qsf / (rau0 rcp e3T) for k= mbkt
      !!       Where Qsf is the geothermal heat flux.
      !!
      !! ** Action  : - update the temperature trends (ta) with the trend of
      !!                the ocean bottom boundary condition
      !!
      !! References : Stein, C. A., and S. Stein, 1992, Nature, 359, 123-129.
      !!              Emile-Geay and Madec, 2009, Ocean Science.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, ik    ! dummy loop indices
      REAL(wp) ::   zqgh_trd      ! geothermal heat flux trend
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   ztrdt
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_bbc')
      !
      IF( l_trdtra )   THEN         ! Save ta and sa trends
         CALL wrk_alloc( jpi, jpj, jpk, ztrdt )
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem)
      ENDIF
      !
      !                             !  Add the geothermal heat flux trend on temperature
      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            ik = mbkt(ji,jj)
            zqgh_trd = qgh_trd0(ji,jj) / e3t_0(ji,jj,ik)
            tsa(ji,jj,ik,jp_tem) = tsa(ji,jj,ik,jp_tem) + zqgh_trd
         END DO
      END DO
      !
      CALL lbc_lnk( tsa(:,:,:,jp_tem) , 'T', 1. )
      !
      IF( l_trdtra ) THEN        ! Save the geothermal heat flux trend for diagnostics
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem) - ztrdt(:,:,:)
         CALL trd_tra( kt, 'TRA', jp_tem, jptra_bbc, ztrdt )
         CALL wrk_dealloc( jpi, jpj, jpk, ztrdt )
      ENDIF
      !
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' bbc  - Ta: ', mask1=tmask, clinfo3='tra-ta' )
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_bbc')
      !
   END SUBROUTINE tra_bbc


   SUBROUTINE tra_bbc_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_bbc_init  ***
      !!
      !! ** Purpose :   Compute once for all the trend associated with geothermal
      !!              heating that will be applied at each time step at the
      !!              last ocean level
      !!
      !! ** Method  :   Read the nambbc namelist and check the parameters.
      !!
      !! ** Input   : - Namlist nambbc
      !!              - NetCDF file  : geothermal_heating.nc ( if necessary )
      !!
      !! ** Action  : - read/fix the geothermal heat qgh_trd0
      !!----------------------------------------------------------------------
      USE iom
      !!
      INTEGER  ::   ji, jj              ! dummy loop indices
      INTEGER  ::   inum                ! temporary logical unit
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      INTEGER  ::   ierror              ! local integer
      !
      TYPE(FLD_N)        ::   sn_qgh    ! informations about the geotherm. field to be read
      CHARACTER(len=256) ::   cn_dir    ! Root directory for location of ssr files
      !
      NAMELIST/nambbc/ln_trabbc, nn_geoflx, rn_geoflx_cst, sn_qgh, cn_dir 
      !!----------------------------------------------------------------------

      REWIND( numnam_ref )              ! Namelist nambbc in reference namelist : Bottom momentum boundary condition
      READ  ( numnam_ref, nambbc, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nambbc in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist nambbc in configuration namelist : Bottom momentum boundary condition
      READ  ( numnam_cfg, nambbc, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nambbc in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, nambbc )

      IF(lwp) THEN                     ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'tra_bbc : Bottom Boundary Condition (bbc), apply a Geothermal heating'
         WRITE(numout,*) '~~~~~~~   '
         WRITE(numout,*) '   Namelist nambbc : set bbc parameters'
         WRITE(numout,*) '      Apply a geothermal heating at ocean bottom   ln_trabbc     = ', ln_trabbc
         WRITE(numout,*) '      type of geothermal flux                      nn_geoflx     = ', nn_geoflx
         WRITE(numout,*) '      Constant geothermal flux value               rn_geoflx_cst = ', rn_geoflx_cst
         WRITE(numout,*)
      ENDIF

      IF( ln_trabbc ) THEN             !==  geothermal heating  ==!
         !
         ALLOCATE( qgh_trd0(jpi,jpj) )    ! allocation
         !
         SELECT CASE ( nn_geoflx )        ! geothermal heat flux / (rauO * Cp)
         !
         CASE ( 1 )                          !* constant flux
            IF(lwp) WRITE(numout,*) '      *** constant heat flux  =   ', rn_geoflx_cst
            qgh_trd0(:,:) = r1_rau0_rcp * rn_geoflx_cst
            !
         CASE ( 2 )                          !* variable geothermal heat flux : read the geothermal fluxes in mW/m2
            IF(lwp) WRITE(numout,*) '      *** variable geothermal heat flux'
            !
            ALLOCATE( sf_qgh(1), STAT=ierror )
            IF( ierror > 0 ) THEN
               CALL ctl_stop( 'tra_bbc_init: unable to allocate sf_qgh structure' )   ;
               RETURN
            ENDIF
            ALLOCATE( sf_qgh(1)%fnow(jpi,jpj,1)   )
            IF( sn_qgh%ln_tint )ALLOCATE( sf_qgh(1)%fdta(jpi,jpj,1,2) )
            ! fill sf_chl with sn_chl and control print
            CALL fld_fill( sf_qgh, (/ sn_qgh /), cn_dir, 'tra_bbc_init',   &
               &          'bottom temperature boundary condition', 'nambbc' )

            CALL fld_read( nit000, 1, sf_qgh )                         ! Read qgh data
            qgh_trd0(:,:) = r1_rau0_rcp * sf_qgh(1)%fnow(:,:,1) * 1.e-3 ! conversion in W/m2
            !
         CASE DEFAULT
            WRITE(ctmp1,*) '     bad flag value for nn_geoflx = ', nn_geoflx
            CALL ctl_stop( ctmp1 )
            !
         END SELECT
         !
      ELSE
         IF(lwp) WRITE(numout,*) '      *** no geothermal heat flux'
      ENDIF
      !
   END SUBROUTINE tra_bbc_init

   !!======================================================================
END MODULE trabbc
