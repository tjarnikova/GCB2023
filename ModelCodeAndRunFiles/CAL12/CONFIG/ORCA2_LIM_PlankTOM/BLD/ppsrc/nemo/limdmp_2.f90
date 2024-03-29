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

MODULE limdmp_2
   !!======================================================================
   !!                       ***  MODULE limdmp_2   ***
   !!  LIM-2 ice model : restoring Ice thickness and Fraction leads
   !!======================================================================
   !! History :   2.0  !  2004-04 (S. Theetten) Original code
   !!             3.3  !  2010-06 (J.-M. Molines) use of fldread
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_lim2'                                    LIM 2.0 sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_dmp_2     : ice model damping
   !!----------------------------------------------------------------------
   USE ice_2          ! ice variables 
   USE sbc_oce, ONLY : nn_fsbc ! for fldread
   USE dom_oce        ! for mi0; mi1 etc ...
   USE fldread        ! read input fields
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   lim_dmp_2     ! called by sbc_ice_lim2

   INTEGER  , PARAMETER :: jp_hicif = 1 , jp_frld = 2
   REAL(wp) , ALLOCATABLE, DIMENSION(:,:,:) ::   resto_ice   ! restoring coeff. on ICE   [s-1]
   TYPE(FLD), ALLOCATABLE, DIMENSION(:)     ::   sf_icedmp   ! structure of ice damping input
   
   !! * Substitution
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
   !! NEMO/LIM 3.3 , UCL-NEMO-consortium (2010) 
   !! $Id: limdmp_2.F90 4624 2014-04-28 12:09:03Z acc $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE lim_dmp_2( kt )
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE lim_dmp_2  ***
      !!
      !! ** purpose :   restore ice thickness and lead fraction
      !!
      !! ** method  :   restore ice thickness and lead fraction using a restoring
      !!              coefficient defined by the user in lim_dmp_init
      !!
      !! ** Action  : - update hicif and frld  
      !!
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step
      !
      INTEGER  ::   ji, jj         ! dummy loop indices
      REAL(wp) ::   zfrld, zhice   ! local scalars
      !!---------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN 
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'lim_dmp_2 : Ice thickness and ice concentration restoring'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
         !
         ! ice_resto_init create resto_ice (in 1/s) for restoring ice parameters near open boundaries.
         ! Double check this routine to verify if it corresponds to your config
         CALL lim_dmp_init
      ENDIF
      !
      IF( ln_limdmp ) THEN   ! ice restoring in this case
         !
         CALL fld_read( kt, nn_fsbc, sf_icedmp )
         !
!CDIR COLLAPSE
         hicif(:,:) = MAX( 0._wp,                     &        ! h >= 0         avoid spurious out of physical range
            &         hicif(:,:) - rdt_ice * resto_ice(:,:,1) * ( hicif(:,:) - sf_icedmp(jp_hicif)%fnow(:,:,1) )  ) 
!CDIR COLLAPSE
         frld (:,:) = MAX( 0._wp, MIN( 1._wp,         &        ! 0<= frld<=1    values which blow the run up
            &         frld (:,:) - rdt_ice * resto_ice(:,:,1) * ( frld (:,:) - sf_icedmp(jp_frld )%fnow(:,:,1) )  )  )
         !
      ENDIF
      !
   END SUBROUTINE lim_dmp_2


   SUBROUTINE lim_dmp_init
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE lim_dmp_init  ***
      !!
      !! ** Purpose :   set the coefficient for the ice thickness and lead fraction restoring
      !!
      !! ** Method  :   restoring is used to mimic ice open boundaries.
      !!              the restoring coef. (a 2D array) has to be defined by the user.
      !!              here is given as an example a restoring along north and south boundaries
      !!      
      !! ** Action  :   define resto_ice(:,:,1)
      !!----------------------------------------------------------------------
      INTEGER  :: ji, jj, jk       ! dummy loop indices
      INTEGER  :: irelax, ierror   ! error flag for allocation
      INTEGER  ::   ios            ! Local integer output status for namelist read
      !
      REAL(wp) :: zdmpmax, zdmpmin, zfactor, zreltim ! temporary scalar
      !
      CHARACTER(len=100)           ::   cn_dir       ! Root directory for location of ssr files
      TYPE(FLD_N), DIMENSION (2)   ::   sl_icedmp    ! informations about the icedmp  field to be read
      TYPE(FLD_N)                  ::   sn_hicif     ! 
      TYPE(FLD_N)                  ::   sn_frld      ! 
      NAMELIST/namice_dmp/ cn_dir, ln_limdmp, sn_hicif, sn_frld
      !!----------------------------------------------------------------------
      !
      ! 1)  initialize fld read structure for input data 
      !     --------------------------------------------
                  
      REWIND( numnam_ice_ref )              ! Namelist namice_dmp in reference namelist : Ice restoring
      READ  ( numnam_ice_ref, namice_dmp, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namice_dmp in reference namelist', lwp )

      REWIND( numnam_ice_cfg )              ! Namelist  namice_dmp in configuration namelist : Ice restoring
      READ  ( numnam_ice_cfg, namice_dmp, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namice_dmp in configuration namelist', lwp )
      IF(lwm) WRITE ( numoni, namice_dmp )
      !
      IF ( lwp ) THEN                     !* control print
         WRITE (numout,*)'     lim_dmp_init : lim_dmp initialization ' 
         WRITE (numout,*)'       Namelist namicedmp read '
         WRITE (numout,*)'         Ice restoring (T) or not (F) ln_limdmp =', ln_limdmp 
         WRITE (numout,*)
         WRITE (numout,*)'     CAUTION : here hard coded ice restoring along northern and southern boundaries'
         WRITE (numout,*)'               adapt the lim_dmp_init routine to your needs'
      ENDIF

      ! 2)  initialise resto_ice    ==>  config dependant !
      !     --------------------         ++++++++++++++++
      !
      IF( ln_limdmp ) THEN                !* ice restoring is used, follow initialization
         ! 
         sl_icedmp ( jp_hicif ) = sn_hicif
         sl_icedmp ( jp_frld  ) = sn_frld
         ALLOCATE ( sf_icedmp (2) , resto_ice(jpi,jpj,1), STAT=ierror )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'lim_dmp_init: unable to allocate sf_icedmp structure or resto_ice array' )   ;   RETURN
         ENDIF
         ALLOCATE( sf_icedmp(jp_hicif)%fnow(jpi,jpj,1) , sf_icedmp(jp_hicif)%fdta(jpi,jpj,1,2) )
         ALLOCATE( sf_icedmp(jp_frld )%fnow(jpi,jpj,1) , sf_icedmp(jp_frld )%fdta(jpi,jpj,1,2) )
         !                         ! fill sf_icedmp with sn_icedmp and control print
         CALL fld_fill( sf_icedmp, sl_icedmp, cn_dir, 'lim_dmp_init', 'Ice  restoring input data', 'namicedmp' )
      
         resto_ice(:,:,:) = 0._wp
         !      Re-calculate the North and South boundary restoring term
         !      because those boundaries may change with the prescribed zoom area.
         !
         irelax  = 16                     ! width of buffer zone with respect to close boundary
         zdmpmax = 10._wp                 ! max restoring time scale  (days) (low restoring)
         zdmpmin = rdt_ice / 86400._wp    ! min restoring time scale  (days) (high restoring)
         !                                ! days / grid-point
         zfactor = ( zdmpmax - zdmpmin ) / REAL( irelax, wp )

         !    South boundary restoring term
         ! REM: if there is no ice in the model and in the data, 
         !      no restoring even with non zero resto_ice
         DO jj = mj0(jpjzoom - 1 + 1), mj1(jpjzoom -1 + irelax)
            zreltim = zdmpmin + zfactor * ( mjg(jj) - jpjzoom + 1 )
            resto_ice(:,jj,:) = 1._wp / ( zreltim * 86400._wp )
         END DO

         ! North boundary restoring term
         DO jj =  mj0(jpjzoom -1 + jpjglo - irelax), mj1(jpjzoom - 1 + jpjglo)
            zreltim = zdmpmin + zfactor * (jpjglo - ( mjg(jj) - jpjzoom + 1 ))
            resto_ice(:,jj,:) = 1.e0 / ( zreltim * 86400 )
         END DO
      ENDIF
      !
   END SUBROUTINE lim_dmp_init
   

   !!======================================================================
END MODULE limdmp_2
