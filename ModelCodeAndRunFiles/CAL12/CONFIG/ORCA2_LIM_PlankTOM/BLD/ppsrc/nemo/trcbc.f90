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

MODULE trcbc
   !!======================================================================
   !!                     ***  MODULE  trcdta  ***
   !! TOP :  module for passive tracer boundary conditions
   !!=====================================================================
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP model 
   !!----------------------------------------------------------------------
   !!   trc_dta    : read and time interpolated passive tracer data
   !!----------------------------------------------------------------------
   USE par_trc       !  passive tracers parameters
   USE oce_trc       !  shared variables between ocean and passive tracers
   USE trc           !  passive tracers common variables
   USE iom           !  I/O manager
   USE lib_mpp       !  MPP library
   USE fldread       !  read input fields

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_bc_init    ! called in trcini.F90 
   PUBLIC   trc_bc_read    ! called in trcstp.F90 or within

   INTEGER  , SAVE, PUBLIC                             :: nb_trcobc   ! number of tracers with open BC
   INTEGER  , SAVE, PUBLIC                             :: nb_trcsbc   ! number of tracers with surface BC
   INTEGER  , SAVE, PUBLIC                             :: nb_trccbc   ! number of tracers with coastal BC
   INTEGER  , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: n_trc_indobc ! index of tracer with OBC data
   INTEGER  , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: n_trc_indsbc ! index of tracer with SBC data
   INTEGER  , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: n_trc_indcbc ! index of tracer with CBC data
   INTEGER  , SAVE, PUBLIC                             :: ntra_obc     ! MAX( 1, nb_trcxxx ) to avoid compilation error with bounds checking
   INTEGER  , SAVE, PUBLIC                             :: ntra_sbc     ! MAX( 1, nb_trcxxx ) to avoid compilation error with bounds checking
   INTEGER  , SAVE, PUBLIC                             :: ntra_cbc     ! MAX( 1, nb_trcxxx ) to avoid compilation error with bounds checking
   REAL(wp) , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: rf_trofac   ! multiplicative factor for OBCtracer values
   TYPE(FLD), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: sf_trcobc   ! structure of data input OBC (file informations, fields read)
   REAL(wp) , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: rf_trsfac   ! multiplicative factor for SBC tracer values
   TYPE(FLD), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: sf_trcsbc   ! structure of data input SBC (file informations, fields read)
   REAL(wp) , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: rf_trcfac   ! multiplicative factor for CBC tracer values
   TYPE(FLD), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: sf_trccbc   ! structure of data input CBC (file informations, fields read)

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
   !! $Id: trcbc.F90 5215 2015-04-15 16:11:56Z nicolasmartin $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_bc_init(ntrc)
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_bc_init  ***
      !!                    
      !! ** Purpose :   initialisation of passive tracer BC data 
      !! 
      !! ** Method  : - Read namtsd namelist
      !!              - allocates passive tracer BC data structure 
      !!----------------------------------------------------------------------
      !
      INTEGER,INTENT(IN) :: ntrc                           ! number of tracers
      INTEGER            :: jl, jn                         ! dummy loop indices
      INTEGER            :: ierr0, ierr1, ierr2, ierr3     ! temporary integers
      INTEGER            ::  ios                           ! Local integer output status for namelist read
      CHARACTER(len=100) :: clndta, clntrc
      !
      CHARACTER(len=100) :: cn_dir
      TYPE(FLD_N), ALLOCATABLE, DIMENSION(:) :: slf_i  ! local array of namelist informations on the fields to read
      TYPE(FLD_N), DIMENSION(jpmaxtrc) :: sn_trcobc    ! open
      TYPE(FLD_N), DIMENSION(jpmaxtrc) :: sn_trcsbc    ! surface
      TYPE(FLD_N), DIMENSION(jpmaxtrc) :: sn_trccbc    ! coastal
      REAL(wp)   , DIMENSION(jpmaxtrc) :: rn_trofac    ! multiplicative factor for tracer values
      REAL(wp)   , DIMENSION(jpmaxtrc) :: rn_trsfac    ! multiplicative factor for tracer values
      REAL(wp)   , DIMENSION(jpmaxtrc) :: rn_trcfac    ! multiplicative factor for tracer values
      !!
      NAMELIST/namtrc_bc/ cn_dir, sn_trcobc, rn_trofac, sn_trcsbc, rn_trsfac, sn_trccbc, rn_trcfac 
      !!----------------------------------------------------------------------
      IF( nn_timing == 1 )  CALL timing_start('trc_bc_init')
      !
      !  Initialisation and local array allocation
      ierr0 = 0  ;  ierr1 = 0  ;  ierr2 = 0  ;  ierr3 = 0  
      ALLOCATE( slf_i(ntrc), STAT=ierr0 )
      IF( ierr0 > 0 ) THEN
         CALL ctl_stop( 'trc_bc_init: unable to allocate local slf_i' )   ;   RETURN
      ENDIF

      ! Compute the number of tracers to be initialised with open, surface and boundary data
      ALLOCATE( n_trc_indobc(ntrc), STAT=ierr0 )
      IF( ierr0 > 0 ) THEN
         CALL ctl_stop( 'trc_bc_init: unable to allocate n_trc_indobc' )   ;   RETURN
      ENDIF
      nb_trcobc      = 0
      n_trc_indobc(:) = 0
      !
      ALLOCATE( n_trc_indsbc(ntrc), STAT=ierr0 )
      IF( ierr0 > 0 ) THEN
         CALL ctl_stop( 'trc_bc_init: unable to allocate n_trc_indsbc' )   ;   RETURN
      ENDIF
      nb_trcsbc      = 0
      n_trc_indsbc(:) = 0
      !
      ALLOCATE( n_trc_indcbc(ntrc), STAT=ierr0 )
      IF( ierr0 > 0 ) THEN
         CALL ctl_stop( 'trc_bc_init: unable to allocate n_trc_indcbc' )   ;   RETURN
      ENDIF
      nb_trccbc      = 0
      n_trc_indcbc(:) = 0
      !
      DO jn = 1, ntrc
         IF( ln_trc_obc(jn) ) THEN
             nb_trcobc       = nb_trcobc + 1 
             n_trc_indobc(jn) = nb_trcobc 
         ENDIF
         IF( ln_trc_sbc(jn) ) THEN
             nb_trcsbc       = nb_trcsbc + 1
             n_trc_indsbc(jn) = nb_trcsbc
         ENDIF
         IF( ln_trc_cbc(jn) ) THEN
             nb_trccbc       = nb_trccbc + 1
             n_trc_indcbc(jn) = nb_trccbc
         ENDIF
      ENDDO
      ntra_obc = MAX( 1, nb_trcobc )   ! To avoid compilation error with bounds checking
      IF( lwp ) WRITE(numout,*) ' '
      IF( lwp ) WRITE(numout,*) ' Number of passive tracers to be initialized with open boundary data :', nb_trcobc
      IF( lwp ) WRITE(numout,*) ' '
      ntra_sbc = MAX( 1, nb_trcsbc )   ! To avoid compilation error with bounds checking
      IF( lwp ) WRITE(numout,*) ' '
      IF( lwp ) WRITE(numout,*) ' Number of passive tracers to be initialized with surface boundary data :', nb_trcsbc
      IF( lwp ) WRITE(numout,*) ' '
      ntra_cbc = MAX( 1, nb_trccbc )   ! To avoid compilation error with bounds checking
      IF( lwp ) WRITE(numout,*) ' '
      IF( lwp ) WRITE(numout,*) ' Number of passive tracers to be initialized with coastal boundary data :', nb_trccbc
      IF( lwp ) WRITE(numout,*) ' '

      REWIND( numnat_ref )              ! Namelist namtrc_bc in reference namelist : Passive tracer data structure
      READ  ( numnat_ref, namtrc_bc, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc_bc in reference namelist', lwp )

      REWIND( numnat_cfg )              ! Namelist namtrc_bc in configuration namelist : Passive tracer data structure
      READ  ( numnat_cfg, namtrc_bc, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc_bc in configuration namelist', lwp )
      IF(lwm) WRITE ( numont, namtrc_bc )

      ! print some information for each 
      IF( lwp ) THEN
         DO jn = 1, ntrc
            IF( ln_trc_obc(jn) )  THEN    
               clndta = TRIM( sn_trcobc(jn)%clvar ) 
               IF(lwp) WRITE(numout,*) 'Preparing to read OBC data file for passive tracer number :', jn, ' name : ', clndta, & 
               &               ' multiplicative factor : ', rn_trofac(jn)
            ENDIF
            IF( ln_trc_sbc(jn) )  THEN    
               clndta = TRIM( sn_trcsbc(jn)%clvar ) 
               IF(lwp) WRITE(numout,*) 'Preparing to read SBC data file for passive tracer number :', jn, ' name : ', clndta, & 
               &               ' multiplicative factor : ', rn_trsfac(jn)
            ENDIF
            IF( ln_trc_cbc(jn) )  THEN    
               clndta = TRIM( sn_trccbc(jn)%clvar ) 
               IF(lwp) WRITE(numout,*) 'Preparing to read CBC data file for passive tracer number :', jn, ' name : ', clndta, & 
               &               ' multiplicative factor : ', rn_trcfac(jn)
            ENDIF
         END DO
      ENDIF
      !
      ! The following code is written this way to reduce memory usage and repeated for each boundary data
      ! MAV: note that this is just a placeholder and the dimensions must be changed according to 
      !      what will be done with BDY. A new structure will probably need to be included
      !
      ! OPEN Lateral boundary conditions
      IF( nb_trcobc > 0 ) THEN       !  allocate only if the number of tracer to initialise is greater than zero
         ALLOCATE( sf_trcobc(nb_trcobc), rf_trofac(nb_trcobc), STAT=ierr1 )
         IF( ierr1 > 0 ) THEN
            CALL ctl_stop( 'trc_bc_init: unable to allocate  sf_trcobc structure' )   ;   RETURN
         ENDIF
         !
         DO jn = 1, ntrc
            IF( ln_trc_obc(jn) ) THEN      ! update passive tracers arrays with input data read from file
               jl = n_trc_indobc(jn)
               slf_i(jl)    = sn_trcobc(jn)
               rf_trofac(jl) = rn_trofac(jn)
                                            ALLOCATE( sf_trcobc(jl)%fnow(jpi,jpj,jpk)   , STAT=ierr2 )
               IF( sn_trcobc(jn)%ln_tint )  ALLOCATE( sf_trcobc(jl)%fdta(jpi,jpj,jpk,2) , STAT=ierr3 )
               IF( ierr2 + ierr3 > 0 ) THEN
                 CALL ctl_stop( 'trc_bc_init : unable to allocate passive tracer OBC data arrays' )   ;   RETURN
               ENDIF
            ENDIF
            !   
         ENDDO
         !                         ! fill sf_trcdta with slf_i and control print
         CALL fld_fill( sf_trcobc, slf_i, cn_dir, 'trc_bc_init', 'Passive tracer OBC data', 'namtrc_bc' )
         !
      ENDIF
      !
      ! SURFACE Boundary conditions
      IF( nb_trcsbc > 0 ) THEN       !  allocate only if the number of tracer to initialise is greater than zero
         ALLOCATE( sf_trcsbc(nb_trcsbc), rf_trsfac(nb_trcsbc), STAT=ierr1 )
         IF( ierr1 > 0 ) THEN
            CALL ctl_stop( 'trc_bc_init: unable to allocate  sf_trcsbc structure' )   ;   RETURN
         ENDIF
         !
         DO jn = 1, ntrc
            IF( ln_trc_sbc(jn) ) THEN      ! update passive tracers arrays with input data read from file
               jl = n_trc_indsbc(jn)
               slf_i(jl)    = sn_trcsbc(jn)
               rf_trsfac(jl) = rn_trsfac(jn)
                                            ALLOCATE( sf_trcsbc(jl)%fnow(jpi,jpj,1)   , STAT=ierr2 )
               IF( sn_trcsbc(jn)%ln_tint )  ALLOCATE( sf_trcsbc(jl)%fdta(jpi,jpj,1,2) , STAT=ierr3 )
               IF( ierr2 + ierr3 > 0 ) THEN
                 CALL ctl_stop( 'trc_bc_init : unable to allocate passive tracer SBC data arrays' )   ;   RETURN
               ENDIF
            ENDIF
            !   
         ENDDO
         !                         ! fill sf_trcsbc with slf_i and control print
         CALL fld_fill( sf_trcsbc, slf_i, cn_dir, 'trc_bc_init', 'Passive tracer SBC data', 'namtrc_bc' )
         !
      ENDIF
      !
      ! COSTAL Boundary conditions
      IF( nb_trccbc > 0 ) THEN       !  allocate only if the number of tracer to initialise is greater than zero
         ALLOCATE( sf_trccbc(nb_trccbc), rf_trcfac(nb_trccbc), STAT=ierr1 )
         IF( ierr1 > 0 ) THEN
            CALL ctl_stop( 'trc_bc_ini: unable to allocate  sf_trccbc structure' )   ;   RETURN
         ENDIF
         !
         DO jn = 1, ntrc
            IF( ln_trc_cbc(jn) ) THEN      ! update passive tracers arrays with input data read from file
               jl = n_trc_indcbc(jn)
               slf_i(jl)    = sn_trccbc(jn)
               rf_trcfac(jl) = rn_trcfac(jn)
                                            ALLOCATE( sf_trccbc(jl)%fnow(jpi,jpj,1)   , STAT=ierr2 )
               IF( sn_trccbc(jn)%ln_tint )  ALLOCATE( sf_trccbc(jl)%fdta(jpi,jpj,1,2) , STAT=ierr3 )
               IF( ierr2 + ierr3 > 0 ) THEN
                 CALL ctl_stop( 'trc_bc_ini : unable to allocate passive tracer CBC data arrays' )   ;   RETURN
               ENDIF
            ENDIF
            !   
         ENDDO
         !                         ! fill sf_trccbc with slf_i and control print
         CALL fld_fill( sf_trccbc, slf_i, cn_dir, 'trc_bc_init', 'Passive tracer CBC data', 'namtrc_bc' )
         !
      ENDIF
 
      DEALLOCATE( slf_i )          ! deallocate local field structure
      IF( nn_timing == 1 )  CALL timing_stop('trc_bc_init')

   END SUBROUTINE trc_bc_init


   SUBROUTINE trc_bc_read(kt)
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_bc_init  ***
      !!
      !! ** Purpose :  Read passive tracer Boundary Conditions data
      !!
      !! ** Method  :  Read BC inputs and update data structures using fldread
      !!              
      !!----------------------------------------------------------------------
   
      ! NEMO
      USE fldread
      
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index

      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_bc_read')

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_bc_read : Surface boundary conditions for passive tracers.'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
      ENDIF

      ! OPEN boundary conditions: DOES NOT WORK. Waiting for stable BDY
      IF( nb_trcobc > 0 ) THEN
        if (lwp) write(numout,'(a,i5,a,i5)') '   reading OBC data for ', nb_trcobc ,' variables at step ', kt
        CALL fld_read(kt,1,sf_trcobc)
        ! vertical interpolation on s-grid and partial step to be added
      ENDIF

      ! SURFACE boundary conditions       
      IF( nb_trcsbc > 0 ) THEN
        if (lwp) write(numout,'(a,i5,a,i5)') '   reading SBC data for ', nb_trcsbc ,' variables at step ', kt
        CALL fld_read(kt,1,sf_trcsbc)
      ENDIF

      ! COASTAL boundary conditions       
      IF( nb_trccbc > 0 ) THEN
        if (lwp) write(numout,'(a,i5,a,i5)') '   reading CBC data for ', nb_trccbc ,' variables at step ', kt
        CALL fld_read(kt,1,sf_trccbc)
      ENDIF   
      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_bc_read')
      !       

   END SUBROUTINE trc_bc_read

   !!======================================================================
END MODULE trcbc
