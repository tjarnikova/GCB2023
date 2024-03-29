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

MODULE dtauvd
   !!======================================================================
   !!                     ***  MODULE  dtauvd  ***
   !! Ocean data  :  read ocean U & V current data from gridded data
   !!======================================================================
   !! History :  3.5   ! 2013-08  (D. Calvert)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dta_uvd_init   : read namelist and allocate data structures
   !!   dta_uvd        : read and time-interpolate ocean U & V current data
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE fldread         ! read input fields
   USE in_out_manager  ! I/O manager
   USE phycst          ! physical constants
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! Memory allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dta_uvd_init   ! called by nemogcm.F90
   PUBLIC   dta_uvd        ! called by istate.F90 and dyndmp.90

   LOGICAL , PUBLIC ::   ln_uvd_init         ! Flag to initialise with U & V current data
   LOGICAL , PUBLIC ::   ln_uvd_dyndmp       ! Flag for Newtonian damping toward U & V current data

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_uvd   ! structure for input U & V current (file information and data)

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
   !! $Id: dtauvd.F90 5215 2015-04-15 16:11:56Z nicolasmartin $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dta_uvd_init( ld_dyndmp )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dta_uvd_init  ***
      !!                    
      !! ** Purpose :   initialization of U & V current input data 
      !! 
      !! ** Method  : - read namc1d_uvd namelist
      !!              - allocate U & V current data structure
      !!              - fld_fill data structure with namelist information
      !!----------------------------------------------------------------------
      LOGICAL, INTENT(in), OPTIONAL ::   ld_dyndmp         ! force the initialization when dyndmp is used
      !
      INTEGER ::   ierr0, ierr1, ierr2, ierr3              ! temporary integers
      !
      CHARACTER(len=100)            ::   cn_dir            ! Root directory for location of files to be used
      TYPE(FLD_N), DIMENSION(2)     ::   suv_i             ! Combined U & V namelist information
      TYPE(FLD_N)                   ::   sn_ucur, sn_vcur  ! U & V data namelist information
      !!
      NAMELIST/namc1d_uvd/ ln_uvd_init, ln_uvd_dyndmp, cn_dir, sn_ucur, sn_vcur
      INTEGER  ::   ios
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dta_uvd_init')
      !
      ierr0 = 0  ;  ierr1 = 0  ;  ierr2 = 0  ;  ierr3 = 0

      REWIND( numnam_ref )              ! Namelist namc1d_uvd in reference namelist : 
      READ  ( numnam_ref, namc1d_uvd, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namc1d_uvd in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namc1d_uvd in configuration namelist : Parameters of the run
      READ  ( numnam_cfg, namc1d_uvd, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namc1d_uvd in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namc1d_uvd )

      !                             ! force the initialization when dyndmp is used
      IF( PRESENT( ld_dyndmp ) )   ln_uvd_dyndmp = .TRUE.
      
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'dta_uvd_init : U & V current data '
         WRITE(numout,*) '~~~~~~~~~~~~ '
         WRITE(numout,*) '   Namelist namc1d_uvd : Set flags'
         WRITE(numout,*) '      Initialization of ocean U & V current with input data   ln_uvd_init   = ', ln_uvd_init
         WRITE(numout,*) '      Damping of ocean U & V current toward input data        ln_uvd_dyndmp = ', ln_uvd_dyndmp
         WRITE(numout,*)
         IF( .NOT. ln_uvd_init .AND. .NOT. ln_uvd_dyndmp ) THEN
            WRITE(numout,*)
            WRITE(numout,*) '   U & V current data not used'
         ENDIF
      ENDIF
      !                             ! no initialization when restarting
      IF( ln_rstart .AND. ln_uvd_init ) THEN
         CALL ctl_warn( 'dta_uvd_init: ocean restart and U & V current data initialization, ',   &
            &           'we keep the restart U & V current values and set ln_uvd_init to FALSE' )
         ln_uvd_init = .FALSE.
      ENDIF

      !
      IF(  ln_uvd_init .OR. ln_uvd_dyndmp  ) THEN
         !                          !==   allocate the data arrays   ==!
         ALLOCATE( sf_uvd(2), STAT=ierr0 )
         IF( ierr0 > 0 ) THEN
            CALL ctl_stop( 'dta_uvd_init: unable to allocate sf_uvd structure' )             ;   RETURN
         ENDIF
         !
                                 ALLOCATE( sf_uvd(1)%fnow(jpi,jpj,jpk)   , STAT=ierr0 )
         IF( sn_ucur%ln_tint )   ALLOCATE( sf_uvd(1)%fdta(jpi,jpj,jpk,2) , STAT=ierr1 )
                                 ALLOCATE( sf_uvd(2)%fnow(jpi,jpj,jpk)   , STAT=ierr2 )
         IF( sn_vcur%ln_tint )   ALLOCATE( sf_uvd(2)%fdta(jpi,jpj,jpk,2) , STAT=ierr3 )
         !
         IF( ierr0 + ierr1 + ierr2 + ierr3 > 0 ) THEN
            CALL ctl_stop( 'dta_uvd_init : unable to allocate U & V current data arrays' )   ;   RETURN
         ENDIF
         !                          !==   fill sf_uvd with sn_ucur, sn_vcur and control print   ==!
         suv_i(1) = sn_ucur   ;   suv_i(2) = sn_vcur
         CALL fld_fill( sf_uvd, suv_i, cn_dir, 'dta_uvd', 'U & V current data', 'namc1d_uvd' )
         !
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('dta_uvd_init')
      !
   END SUBROUTINE dta_uvd_init


   SUBROUTINE dta_uvd( kt, puvd )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dta_uvd  ***
      !!                    
      !! ** Purpose :   provides U & V current data at time step kt
      !! 
      !! ** Method  : - call fldread routine
      !!              - ORCA_R2: make some hand made alterations to the data (EMPTY)
      !!              - s- or mixed s-zps coordinate: vertical interpolation onto model mesh
      !!              - zps coordinate: vertical interpolation onto last partial level
      !!              - ln_uvd_dyndmp=False: deallocate the U & V current data structure,
      !!                                     as the data is no longer used
      !!
      !! ** Action  :   puvd,  U & V current data interpolated onto model mesh at time-step kt
      !!----------------------------------------------------------------------
      INTEGER                           , INTENT(in   ) ::   kt     ! ocean time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk,2), INTENT(  out) ::   puvd   ! U & V current data
      !
      INTEGER ::   ji, jj, jk, jl, jkk               ! dummy loop indicies
      INTEGER ::   ik, il0, il1, ii0, ii1, ij0, ij1  ! local integers
      REAL(wp)::   zl, zi                            ! local floats
      REAL(wp), POINTER, DIMENSION(:) ::  zup, zvp   ! 1D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dta_uvd')
      !
      CALL fld_read( kt, 1, sf_uvd )      !==   read U & V current data at time step kt   ==!
      !
      !
      !                                   !==   ORCA_R2 configuration and U & V current damping   ==! 
      IF( cp_cfg == "orca" .AND. jp_cfg == 2 .AND. ln_uvd_dyndmp ) THEN    ! some hand made alterations
         !!! EMPTY- to be added for running in 3D context !!!
      ENDIF
      !
      puvd(:,:,:,1) = sf_uvd(1)%fnow(:,:,:)                 ! NO mask
      puvd(:,:,:,2) = sf_uvd(2)%fnow(:,:,:) 
      !
      IF( ln_sco ) THEN                   !==   s- or mixed s-zps-coordinate   ==!
         !
         CALL wrk_alloc( jpk, zup, zvp )
         !
         IF( kt == nit000 .AND. lwp )THEN
            WRITE(numout,*)
            WRITE(numout,*) 'dta_uvd: interpolate U & V current data onto the s- or mixed s-z-coordinate mesh'
         ENDIF
         !
         DO jj = 1, jpj                   ! vertical interpolation of U & V current:
            DO ji = 1, jpi                ! determines the interpolated U & V current profiles at each (i,j) point
               DO jk = 1, jpk
                  zl = gdept_0(ji,jj,jk)
                  IF    ( zl < gdept_1d(1  ) ) THEN          ! extrapolate above the first level of data
                     zup(jk) =  puvd(ji,jj,1    ,1)
                     zvp(jk) =  puvd(ji,jj,1    ,2)
                  ELSEIF( zl > gdept_1d(jpk) ) THEN          ! extrapolate below the last level of data
                     zup(jk) =  puvd(ji,jj,jpkm1,1)
                     zvp(jk) =  puvd(ji,jj,jpkm1,2)
                  ELSE                                      ! inbetween : vertical interpolation between jkk & jkk+1
                     DO jkk = 1, jpkm1                      ! when  gdept(jkk) < zl < gdept(jkk+1)
                        IF( (zl-gdept_1d(jkk)) * (zl-gdept_1d(jkk+1)) <= 0._wp ) THEN
                           zi = ( zl - gdept_1d(jkk) ) / (gdept_1d(jkk+1)-gdept_1d(jkk))
                           zup(jk) = puvd(ji,jj,jkk,1) + ( puvd(ji,jj,jkk+1,1 ) - puvd(ji,jj,jkk,1) ) * zi 
                           zvp(jk) = puvd(ji,jj,jkk,2) + ( puvd(ji,jj,jkk+1,2 ) - puvd(ji,jj,jkk,2) ) * zi
                        ENDIF
                     END DO
                  ENDIF
               END DO
               DO jk = 1, jpkm1           ! apply mask
                  puvd(ji,jj,jk,1) = zup(jk) * umask(ji,jj,jk)
                  puvd(ji,jj,jk,2) = zvp(jk) * vmask(ji,jj,jk)
               END DO
               puvd(ji,jj,jpk,1) = 0._wp
               puvd(ji,jj,jpk,2) = 0._wp
            END DO
         END DO
         ! 
         CALL wrk_dealloc( jpk, zup, zvp )
         ! 
      ELSE                                !==   z- or zps- coordinate   ==!
         !                             
         puvd(:,:,:,1) = puvd(:,:,:,1) * umask(:,:,:)       ! apply mask
         puvd(:,:,:,2) = puvd(:,:,:,2) * vmask(:,:,:)
         !
         IF( ln_zps ) THEN                ! zps-coordinate (partial steps) interpolation at the last ocean level
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ik = mbkt(ji,jj) 
                  IF( ik > 1 ) THEN
                     zl = ( gdept_1d(ik) - gdept_0(ji,jj,ik) ) / ( gdept_1d(ik) - gdept_1d(ik-1) )
                     puvd(ji,jj,ik,1) = (1.-zl) * puvd(ji,jj,ik,1) + zl * puvd(ji,jj,ik-1,1)
                     puvd(ji,jj,ik,2) = (1.-zl) * puvd(ji,jj,ik,2) + zl * puvd(ji,jj,ik-1,2)
                  ENDIF
               END DO
            END DO
         ENDIF
         !
      ENDIF
      !
      IF( lwp .AND. kt == nit000 ) THEN   ! control print
         WRITE(numout,*) ' U current '
         WRITE(numout,*)
         WRITE(numout,*)'  level = 1'
         CALL prihre( puvd(:,:,1    ,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
         WRITE(numout,*)'  level = ', jpk/2
         CALL prihre( puvd(:,:,jpk/2,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
         WRITE(numout,*)'  level = ', jpkm1
         CALL prihre( puvd(:,:,jpkm1,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
         WRITE(numout,*)
         WRITE(numout,*) ' V current '
         WRITE(numout,*)
         WRITE(numout,*)'  level = 1'
         CALL prihre( puvd(:,:,1    ,2), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
         WRITE(numout,*)'  level = ', jpk/2
         CALL prihre( puvd(:,:,jpk/2,2), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
         WRITE(numout,*)'  level = ', jpkm1
         CALL prihre( puvd(:,:,jpkm1,2), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
         WRITE(numout,*)
      ENDIF
      !
      IF( .NOT. ln_uvd_dyndmp    ) THEN   !==   deallocate U & V current structure   ==! 
         !                                !==   (data used only for initialization)  ==!
         IF(lwp) WRITE(numout,*) 'dta_uvd: deallocate U & V current arrays as they are only used to initialize the run'
                                   DEALLOCATE( sf_uvd(1)%fnow )     ! U current arrays in the structure
         IF( sf_uvd(1)%ln_tint )   DEALLOCATE( sf_uvd(1)%fdta )
                                   DEALLOCATE( sf_uvd(2)%fnow )     ! V current arrays in the structure
         IF( sf_uvd(2)%ln_tint )   DEALLOCATE( sf_uvd(2)%fdta )
                                   DEALLOCATE( sf_uvd         )     ! the structure itself
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('dta_uvd')
      !
   END SUBROUTINE dta_uvd

   !!======================================================================
END MODULE dtauvd
