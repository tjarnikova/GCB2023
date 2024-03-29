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

MODULE iceini_2
   !!======================================================================
   !!                       ***  MODULE iceini   ***
   !!   Sea-ice model : LIM 2.0 Sea ice model Initialization
   !!======================================================================
   !! History :  1.0  ! 2002-08  (G. Madec)  F90: Free form and modules
   !!            2.0  ! 2003-08  (C. Ethe)  add ice_run
   !!            3.3  ! 2009-05  (G. Garric, C. Bricaud) addition of the lim2_evp case
   !!            4.0  ! 2011-02  (G. Madec) dynamical allocation
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_lim2' :                                  LIM 2.0 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_init_2    : sea-ice model initialization
   !!   ice_run_2     : Definition some run parameter for ice model
   !!----------------------------------------------------------------------
   USE phycst         ! physical constants
   USE dom_oce        ! ocean domain
   USE sbc_oce        ! surface boundary condition: ocean
   USE sbc_ice        ! LIM2 surface boundary condition
   USE dom_ice_2      ! LIM2 ice domain
   USE par_ice_2      ! LIM2 parameters
   USE thd_ice_2      ! LIM2 thermodynamical variables
   USE ice_2          ! LIM2 ice variable
   USE limmsh_2       ! LIM2 mesh
   USE limistate_2    ! LIM2 initial state
   USE limrst_2       ! LIM2 restart
   USE limsbc_2       ! LIM2 surface boundary condition
   USE limdia_2
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_init_2               ! called by sbcice_lim_2.F90

   !!----------------------------------------------------------------------
   !! NEMO/LIM2 4.0 , UCL - NEMO Consortium (2011)
   !! $Id: iceini_2.F90 5385 2015-06-09 13:50:42Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_init_2
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ice_init_2  ***
      !!
      !! ** purpose :   initialisation of LIM-2 domain and variables  
      !!----------------------------------------------------------------------
      INTEGER :: ierr
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ice_init_2 : LIM-2 sea-ice - initialization'
         WRITE(numout,*) '~~~~~~~~~~~   '
      ENDIF
      !                                
      !                                ! Allocate the ice arrays
      ierr =        ice_alloc_2    ()       ! ice variables
      ierr = ierr + dom_ice_alloc_2()       ! domain
      ierr = ierr + sbc_ice_alloc  ()       ! surface forcing
      ierr = ierr + thd_ice_alloc_2()       ! thermodynamics

      IF( lk_mpp    )   CALL mpp_sum( ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'ice_init_2 : unable to allocate ice arrays' )

      !                                ! adequation jpk versus ice/snow layers
      IF( jpl > jpk  .OR.  jplayersp1 > jpk  )   CALL ctl_stop( 'STOP',           &
         &     'ice_init: the 3rd dimension of workspace arrays is too small.',   &
         &     'use more ocean levels or less ice layers/categories.' )

      !                                ! Open the reference and configuration namelist files and namelist output file 
      CALL ctl_opn( numnam_ice_ref, 'namelist_ice_ref',    'OLD',     'FORMATTED', 'SEQUENTIAL', -1, numout, lwp ) 
      CALL ctl_opn( numnam_ice_cfg, 'namelist_ice_cfg',    'OLD',     'FORMATTED', 'SEQUENTIAL', -1, numout, lwp )
      IF(lwm) CALL ctl_opn( numoni, 'output.namelist.ice', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, 1 )

      !                                ! Open the namelist file 
      !    
      CALL ice_run_2                   ! read in namelist some run parameters
      !          
      rdt_ice = nn_fsbc * rdttra(1)    ! sea-ice time step
      numit   = nit000 - 1
      !
      CALL lim_msh_2                   ! ice mesh initialization
      !
      !                                ! Initial sea-ice state
      IF( .NOT.ln_rstart ) THEN   ;   CALL lim_istate_2     ! start from rest: sea-ice deduced from sst
      ELSE                        ;   CALL lim_rst_read_2   ! start from a restart file
      ENDIF
      !
      IF( .NOT.lk_mpp )               CALL lim_dia_init_2  ! online diagnostics in mono proc only
      !
      tn_ice(:,:,1) = sist(:,:)        ! ice temperature  known by the ocean
      fr_i  (:,:)   = 1.0 - frld(:,:)  ! sea-ice fraction known by the ocean
      !
      CALL lim_sbc_init_2              ! ice surface boundary condition   
      !
      IF( lk_lim2_vp )   THEN   ;   IF(lwp) WRITE(numout,*) '                VP  rheology - B-grid case'
      ELSE                      ;   IF(lwp) WRITE(numout,*) '                EVP rheology - C-grid case'
      ENDIF
      !
   END SUBROUTINE ice_init_2


   SUBROUTINE ice_run_2
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_run_2 ***
      !!                 
      !! ** Purpose :   Definition some run parameter for ice model
      !!
      !! ** Method  :   Read the namicerun namelist and check the parameter 
      !!       values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namicerun
      !!-------------------------------------------------------------------
      NAMELIST/namicerun/ cn_icerst_in, cn_icerst_indir, cn_icerst_out, cn_icerst_outdir, &
                          ln_limdyn, ln_limdmp, acrit, hsndif, hicdif
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      !!-------------------------------------------------------------------
      !                    
      REWIND( numnam_ice_ref )              ! Namelist namicerun in reference namelist : Parameters for ice
      READ  ( numnam_ice_ref, namicerun, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namicerun in reference namelist', lwp )

      REWIND( numnam_ice_cfg )              ! Namelist namicerun in configuration namelist : Parameters for ice
      READ  ( numnam_ice_cfg, namicerun, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namicerun in configuration namelist', lwp )
      IF(lwm) WRITE ( numoni, namicerun )
      !
      IF(lwp) THEN                              ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_run : ice share parameters for dynamics/advection/thermo of sea-ice'
         WRITE(numout,*) ' ~~~~~~'
         WRITE(numout,*) '   switch for ice dynamics (1) or not (0)      ln_limdyn   = ', ln_limdyn
         WRITE(numout,*) '   Ice damping                                 ln_limdmp   = ', ln_limdmp
         WRITE(numout,*) '   minimum fraction for leads in the NH (SH)  acrit(1/2)   = ', acrit(:)
         WRITE(numout,*) '   computation of temp. in snow (=0) or not (=9999) hsndif = ', hsndif
         WRITE(numout,*) '   computation of temp. in ice  (=0) or not (=9999) hicdif = ', hicdif
      ENDIF
      !
   END SUBROUTINE ice_run_2


   !!======================================================================
END MODULE iceini_2
