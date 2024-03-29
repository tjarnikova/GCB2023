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

MODULE icbutl
   !!======================================================================
   !!                       ***  MODULE  icbutl  ***
   !! Icebergs:  various iceberg utility routines
   !!======================================================================
   !! History : 3.3.1 !  2010-01  (Martin&Adcroft) Original code
   !!            -    !  2011-03  (Madec)          Part conversion to NEMO form
   !!            -    !                            Removal of mapping from another grid
   !!            -    !  2011-04  (Alderson)       Split into separate modules
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   icb_utl_interp   :
   !!   icb_utl_bilin    :
   !!   icb_utl_bilin_e  :
   !!----------------------------------------------------------------------
   USE par_oce                             ! ocean parameters
   USE dom_oce                             ! ocean domain
   USE in_out_manager                      ! IO parameters
   USE lbclnk                              ! lateral boundary condition
   USE lib_mpp                             ! MPI code and lk_mpp in particular
   USE icb_oce                             ! define iceberg arrays
   USE sbc_oce                             ! ocean surface boundary conditions
   USE ice_2,         ONLY: u_ice, v_ice   ! LIM-2 ice velocities  (CAUTION in C-grid do not use key_vp option)
   USE ice_2,         ONLY: hicif          ! LIM-2 ice thickness

   IMPLICIT NONE
   PRIVATE

   PUBLIC   icb_utl_copy          ! routine called in icbstp module
   PUBLIC   icb_utl_interp        ! routine called in icbdyn, icbthm modules
   PUBLIC   icb_utl_bilin         ! routine called in icbini, icbdyn modules
   PUBLIC   icb_utl_bilin_x       ! routine called in icbdyn module
   PUBLIC   icb_utl_add           ! routine called in icbini.F90, icbclv, icblbc and icbrst modules
   PUBLIC   icb_utl_delete        ! routine called in icblbc, icbthm modules
   PUBLIC   icb_utl_destroy       ! routine called in icbstp module
   PUBLIC   icb_utl_track         ! routine not currently used, retain just in case
   PUBLIC   icb_utl_print_berg    ! routine called in icbthm module
   PUBLIC   icb_utl_print         ! routine called in icbini, icbstp module
   PUBLIC   icb_utl_count         ! routine called in icbdia, icbini, icblbc, icbrst modules
   PUBLIC   icb_utl_incr          ! routine called in icbini, icbclv modules
   PUBLIC   icb_utl_yearday       ! routine called in icbclv, icbstp module
   PUBLIC   icb_utl_mass          ! routine called in icbdia module
   PUBLIC   icb_utl_heat          ! routine called in icbdia module

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2011)
   !! $Id: icbutl.F90 5215 2015-04-15 16:11:56Z nicolasmartin $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!-------------------------------------------------------------------------

CONTAINS

   SUBROUTINE icb_utl_copy()
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE icb_utl_copy  ***
      !!
      !! ** Purpose :   iceberg initialization.
      !!
      !! ** Method  : - blah blah
      !!----------------------------------------------------------------------

      ! copy nemo forcing arrays into iceberg versions with extra halo
      ! only necessary for variables not on T points
      ! and ssh which is used to calculate gradients

      uo_e(:,:) = 0._wp ;   uo_e(1:jpi, 1:jpj) = ssu_m(:,:) * umask(:,:,1)
      vo_e(:,:) = 0._wp ;   vo_e(1:jpi, 1:jpj) = ssv_m(:,:) * vmask(:,:,1)
      ff_e(:,:) = 0._wp ;   ff_e(1:jpi, 1:jpj) = ff   (:,:) 
      tt_e(:,:) = 0._wp ;   tt_e(1:jpi, 1:jpj) = sst_m(:,:)
      fr_e(:,:) = 0._wp ;   fr_e(1:jpi, 1:jpj) = fr_i (:,:)
      ua_e(:,:) = 0._wp ;   ua_e(1:jpi, 1:jpj) = utau (:,:) * umask(:,:,1) ! maybe mask useless because mask applied in sbcblk
      va_e(:,:) = 0._wp ;   va_e(1:jpi, 1:jpj) = vtau (:,:) * vmask(:,:,1) ! maybe mask useless because mask applied in sbcblk

      CALL lbc_lnk_icb( uo_e, 'U', -1._wp, 1, 1 )
      CALL lbc_lnk_icb( vo_e, 'V', -1._wp, 1, 1 )
      CALL lbc_lnk_icb( ff_e, 'F', +1._wp, 1, 1 )
      CALL lbc_lnk_icb( ua_e, 'U', -1._wp, 1, 1 )
      CALL lbc_lnk_icb( va_e, 'V', -1._wp, 1, 1 )
      CALL lbc_lnk_icb( fr_e, 'T', +1._wp, 1, 1 )
      CALL lbc_lnk_icb( tt_e, 'T', +1._wp, 1, 1 )
      hicth(:,:) = 0._wp ;  hicth(1:jpi,1:jpj) = hicif(:,:)  
      CALL lbc_lnk_icb(hicth, 'T', +1._wp, 1, 1 )  

      ui_e(:,:) = 0._wp ;   ui_e(1:jpi, 1:jpj) = u_ice(:,:)
      vi_e(:,:) = 0._wp ;   vi_e(1:jpi, 1:jpj) = v_ice(:,:)

      CALL lbc_lnk_icb( ui_e, 'U', -1._wp, 1, 1 )
      CALL lbc_lnk_icb( vi_e, 'V', -1._wp, 1, 1 )

      !! special for ssh which is used to calculate slope
      !! so fudge some numbers all the way around the boundary

      ssh_e(:,:) = 0._wp ;   ssh_e(1:jpi, 1:jpj) = ssh_m(:,:) * tmask(:,:,1)
      ssh_e(0    ,    :) = ssh_e(1  ,  :)
      ssh_e(jpi+1,    :) = ssh_e(jpi,  :)
      ssh_e(:    ,    0) = ssh_e(:  ,  1)
      ssh_e(:    ,jpj+1) = ssh_e(:  ,jpj)
      ssh_e(0,0)         = ssh_e(1,1)
      ssh_e(jpi+1,0)     = ssh_e(jpi,1)
      ssh_e(0,jpj+1)     = ssh_e(1,jpj)
      ssh_e(jpi+1,jpj+1) = ssh_e(jpi,jpj)
      CALL lbc_lnk_icb( ssh_e, 'T', +1._wp, 1, 1 )
      !
   END SUBROUTINE icb_utl_copy


   SUBROUTINE icb_utl_interp( pi, pe1, puo, pui, pua, pssh_i,   &
      &                       pj, pe2, pvo, pvi, pva, pssh_j,   &
      &                       psst, pcn, phi, pff            )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE icb_utl_interp  ***
      !!
      !! ** Purpose :   interpolation
      !!
      !! ** Method  : - interpolate from various ocean arrays onto iceberg position
      !!
      !!       !!gm  CAUTION here I do not care of the slip/no-slip conditions
      !!             this can be done later (not that easy to do...)
      !!             right now, U is 0 in land so that the coastal value of velocity parallel to the coast
      !!             is half the off shore value, wile the normal-to-the-coast value is zero.
      !!             This is OK as a starting point.
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   ) ::   pi , pj                        ! position in (i,j) referential
      REAL(wp), INTENT(  out) ::   pe1, pe2                       ! i- and j scale factors
      REAL(wp), INTENT(  out) ::   puo, pvo, pui, pvi, pua, pva   ! ocean, ice and wind speeds
      REAL(wp), INTENT(  out) ::   pssh_i, pssh_j                 ! ssh i- & j-gradients
      REAL(wp), INTENT(  out) ::   psst, pcn, phi, pff            ! SST, ice concentration, ice thickness, Coriolis
      !
      REAL(wp) ::   zcd, zmod       ! local scalars
      !!----------------------------------------------------------------------

      pe1 = icb_utl_bilin_e( e1t, e1u, e1v, e1f, pi, pj )     ! scale factors
      pe2 = icb_utl_bilin_e( e2t, e2u, e2v, e2f, pi, pj )
      !
      puo  = icb_utl_bilin_h( uo_e, pi, pj, 'U' )             ! ocean velocities
      pvo  = icb_utl_bilin_h( vo_e, pi, pj, 'V' )
      psst = icb_utl_bilin_h( tt_e, pi, pj, 'T' )             ! SST
      pcn  = icb_utl_bilin_h( fr_e , pi, pj, 'T' )            ! ice concentration
      pff  = icb_utl_bilin_h( ff_e , pi, pj, 'F' )            ! Coriolis parameter
      !
      pua  = icb_utl_bilin_h( ua_e , pi, pj, 'U' )            ! 10m wind
      pva  = icb_utl_bilin_h( va_e , pi, pj, 'V' )            ! here (ua,va) are stress => rough conversion from stress to speed
      zcd  = 1.22_wp * 1.5e-3_wp                              ! air density * drag coefficient
      zmod = 1._wp / MAX(  1.e-20, SQRT(  zcd * SQRT( pua*pua + pva*pva)  )  )
      pua  = pua * zmod                                       ! note: stress module=0 necessarly implies ua=va=0
      pva  = pva * zmod

      pui = icb_utl_bilin_h( ui_e, pi, pj, 'U' )              ! sea-ice velocities
      pvi = icb_utl_bilin_h( vi_e, pi, pj, 'V' )
      phi = icb_utl_bilin_h(hicth, pi, pj, 'T' )              ! ice thickness

      ! Estimate SSH gradient in i- and j-direction (centred evaluation)
      pssh_i = ( icb_utl_bilin_h( ssh_e, pi+0.1_wp, pj, 'T' ) -   &
          &      icb_utl_bilin_h( ssh_e, pi-0.1_wp, pj, 'T' )  ) / ( 0.2_wp * pe1 )
      pssh_j = ( icb_utl_bilin_h( ssh_e, pi, pj+0.1_wp, 'T' ) -   &
          &      icb_utl_bilin_h( ssh_e, pi, pj-0.1_wp, 'T' )  ) / ( 0.2_wp * pe2 )
      !
   END SUBROUTINE icb_utl_interp


   REAL(wp) FUNCTION icb_utl_bilin_h( pfld, pi, pj, cd_type )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION icb_utl_bilin  ***
      !!
      !! ** Purpose :   bilinear interpolation at berg location depending on the grid-point type
      !!                this version deals with extra halo points
      !!
      !!       !!gm  CAUTION an optional argument should be added to handle
      !!             the slip/no-slip conditions  ==>>> to be done later
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(0:jpi+1,0:jpj+1), INTENT(in) ::   pfld      ! field to be interpolated
      REAL(wp)                            , INTENT(in) ::   pi, pj    ! targeted coordinates in (i,j) referential
      CHARACTER(len=1)                    , INTENT(in) ::   cd_type   ! type of pfld array grid-points: = T , U , V or F points
      !
      INTEGER  ::   ii, ij   ! local integer
      REAL(wp) ::   zi, zj   ! local real
      !!----------------------------------------------------------------------
      !
      SELECT CASE ( cd_type )
         CASE ( 'T' )
            ! note that here there is no +0.5 added
            ! since we're looking for four T points containing quadrant we're in of 
            ! current T cell
            ii = MAX(1, INT( pi     ))
            ij = MAX(1, INT( pj     ))    ! T-point
            zi = pi - REAL(ii,wp)
            zj = pj - REAL(ij,wp)
         CASE ( 'U' )
            ii = MAX(1, INT( pi-0.5 ))
            ij = MAX(1, INT( pj     ))    ! U-point
            zi = pi - 0.5 - REAL(ii,wp)
            zj = pj - REAL(ij,wp)
         CASE ( 'V' )
            ii = MAX(1, INT( pi     ))
            ij = MAX(1, INT( pj-0.5 ))    ! V-point
            zi = pi - REAL(ii,wp)
            zj = pj - 0.5 - REAL(ij,wp)
         CASE ( 'F' )
            ii = MAX(1, INT( pi-0.5 ))
            ij = MAX(1, INT( pj-0.5 ))    ! F-point
            zi = pi - 0.5 - REAL(ii,wp)
            zj = pj - 0.5 - REAL(ij,wp)
      END SELECT
      !
      ! find position in this processor. Prevent near edge problems (see #1389)

      if (ii.lt.mig(1)) then
        ii = 1
      else if (ii.gt.mig(jpi)) then
        ii = jpi
      else
        ii  = mi1( ii  )
      end if

      if (ij.lt.mjg(1)) then
        ij = 1
      else if (ij.gt.mjg(jpj)) then
        ij = jpj
      else
        ij  = mj1( ij  )
      end if

      if (ij.eq.jpj) ij=ij-1
      if (ii.eq.jpi) ii=ii-1      

      !
      icb_utl_bilin_h = ( pfld(ii,ij  ) * (1.-zi) + pfld(ii+1,ij  ) * zi ) * (1.-zj)   &
         &            + ( pfld(ii,ij+1) * (1.-zi) + pfld(ii+1,ij+1) * zi ) *     zj
      !
   END FUNCTION icb_utl_bilin_h


   REAL(wp) FUNCTION icb_utl_bilin( pfld, pi, pj, cd_type )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION icb_utl_bilin  ***
      !!
      !! ** Purpose :   bilinear interpolation at berg location depending on the grid-point type
      !!
      !!       !!gm  CAUTION an optional argument should be added to handle
      !!             the slip/no-slip conditions  ==>>> to be done later
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pfld      ! field to be interpolated
      REAL(wp)                    , INTENT(in) ::   pi, pj    ! targeted coordinates in (i,j) referential
      CHARACTER(len=1)            , INTENT(in) ::   cd_type   ! type of pfld array grid-points: = T , U , V or F points
      !
      INTEGER  ::   ii, ij   ! local integer
      REAL(wp) ::   zi, zj   ! local real
      !!----------------------------------------------------------------------
      !
      SELECT CASE ( cd_type )
         CASE ( 'T' )
            ! note that here there is no +0.5 added
            ! since we're looking for four T points containing quadrant we're in of 
            ! current T cell
            ii = MAX(1, INT( pi     ))
            ij = MAX(1, INT( pj     ))    ! T-point
            zi = pi - REAL(ii,wp)
            zj = pj - REAL(ij,wp)
         CASE ( 'U' )
            ii = MAX(1, INT( pi-0.5 ))
            ij = MAX(1, INT( pj     ))    ! U-point
            zi = pi - 0.5 - REAL(ii,wp)
            zj = pj - REAL(ij,wp)
         CASE ( 'V' )
            ii = MAX(1, INT( pi     ))
            ij = MAX(1, INT( pj-0.5 ))    ! V-point
            zi = pi - REAL(ii,wp)
            zj = pj - 0.5 - REAL(ij,wp)
         CASE ( 'F' )
            ii = MAX(1, INT( pi-0.5 ))
            ij = MAX(1, INT( pj-0.5 ))    ! F-point
            zi = pi - 0.5 - REAL(ii,wp)
            zj = pj - 0.5 - REAL(ij,wp)
      END SELECT
      !
      ! find position in this processor. Prevent near edge problems (see #1389)

      if (ii.lt.mig(1)) then
        ii = 1
      else if (ii.gt.mig(jpi)) then
        ii = jpi
      else
        ii  = mi1( ii  )
      end if

      if (ij.lt.mjg(1)) then
        ij = 1
      else if (ij.gt.mjg(jpj)) then
        ij = jpj
      else
        ij  = mj1( ij  )
      end if

      if (ij.eq.jpj) ij=ij-1
      if (ii.eq.jpi) ii=ii-1

      icb_utl_bilin = ( pfld(ii,ij  ) * (1.-zi) + pfld(ii+1,ij  ) * zi ) * (1.-zj)   &
         &          + ( pfld(ii,ij+1) * (1.-zi) + pfld(ii+1,ij+1) * zi ) *     zj
      !
   END FUNCTION icb_utl_bilin


   REAL(wp) FUNCTION icb_utl_bilin_x( pfld, pi, pj )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION icb_utl_bilin_x  ***
      !!
      !! ** Purpose :   bilinear interpolation at berg location depending on the grid-point type
      !!                Special case for interpolating longitude
      !!
      !!       !!gm  CAUTION an optional argument should be added to handle
      !!             the slip/no-slip conditions  ==>>> to be done later
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pfld      ! field to be interpolated
      REAL(wp)                    , INTENT(in) ::   pi, pj    ! targeted coordinates in (i,j) referential
      !
      INTEGER                                  ::   ii, ij   ! local integer
      REAL(wp)                                 ::   zi, zj   ! local real
      REAL(wp)                                 ::   zret     ! local real
      REAL(wp), DIMENSION(4)                   ::   z4
      !!----------------------------------------------------------------------
      !
      ! note that here there is no +0.5 added
      ! since we're looking for four T points containing quadrant we're in of 
      ! current T cell
      ii = MAX(1, INT( pi     ))
      ij = MAX(1, INT( pj     ))    ! T-point
      zi = pi - REAL(ii,wp)
      zj = pj - REAL(ij,wp)
      !
      ! find position in this processor. Prevent near edge problems (see #1389)

      if (ii.lt.mig(1)) then
        ii = 1
      else if (ii.gt.mig(jpi)) then
        ii = jpi
      else
        ii  = mi1( ii  )
      end if

      if (ij.lt.mjg(1)) then
        ij = 1
      else if (ij.gt.mjg(jpj)) then
        ij = jpj
      else
        ij  = mj1( ij  )
      end if

      if (ij.eq.jpj) ij=ij-1
      if (ii.eq.jpi) ii=ii-1

      z4(1) = pfld(ii  ,ij  )
      z4(2) = pfld(ii+1,ij  )
      z4(3) = pfld(ii  ,ij+1)
      z4(4) = pfld(ii+1,ij+1)
      IF( MAXVAL(z4) - MINVAL(z4) > 90._wp ) THEN
         WHERE( z4 < 0._wp ) z4 = z4 + 360._wp
      ENDIF
      !
      zret = (z4(1) * (1.-zi) + z4(2) * zi) * (1.-zj) + (z4(3) * (1.-zi) + z4(4) * zi) * zj
      IF( zret > 180._wp ) zret = zret - 360._wp
      icb_utl_bilin_x = zret
      !
   END FUNCTION icb_utl_bilin_x


   REAL(wp) FUNCTION icb_utl_bilin_e( pet, peu, pev, pef, pi, pj )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION dom_init  ***
      !!
      !! ** Purpose :   bilinear interpolation at berg location of horizontal scale factor
      !! ** Method  :   interpolation done using the 4 nearest grid points among
      !!                t-, u-, v-, and f-points.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in) ::   pet, peu, pev, pef   ! horizontal scale factor to be interpolated at t-,u-,v- & f-pts
      REAL(wp)                , INTENT(in) ::   pi, pj               ! targeted coordinates in (i,j) referential
      !
      INTEGER  ::   ii, ij, icase   ! local integer
      !
      ! weights corresponding to corner points of a T cell quadrant
      REAL(wp) ::   zi, zj          ! local real
      !
      ! values at corner points of a T cell quadrant
      ! 00 = bottom left, 10 = bottom right, 01 = top left, 11 = top right
      REAL(wp) ::   ze00, ze10, ze01, ze11
      !!----------------------------------------------------------------------
      !
      ii = MAX(1, INT( pi ))   ;   ij = MAX(1, INT( pj ))            ! left bottom T-point (i,j) indices

      ! fractional box spacing
      ! 0   <= zi < 0.5, 0   <= zj < 0.5   -->  NW quadrant of current T cell
      ! 0.5 <= zi < 1  , 0   <= zj < 0.5   -->  NE quadrant
      ! 0   <= zi < 0.5, 0.5 <= zj < 1     -->  SE quadrant
      ! 0.5 <= zi < 1  , 0.5 <= zj < 1     -->  SW quadrant

      zi = pi - REAL(ii,wp)          !!gm use here mig, mjg arrays
      zj = pj - REAL(ij,wp)

      ! find position in this processor. Prevent near edge problems (see #1389)

      if (ii.lt.mig(1)) then
        ii = 1
      else if (ii.gt.mig(jpi)) then
        ii = jpi
      else
        ii  = mi1( ii  )
      end if

      if (ij.lt.mjg(1)) then
        ij = 1
      else if (ij.gt.mjg(jpj)) then
        ij = jpj
      else
        ij  = mj1( ij  )
      end if

      if (ij.eq.jpj) ij=ij-1
      if (ii.eq.jpi) ii=ii-1

      IF(    0.0_wp <= zi .AND. zi < 0.5_wp   ) THEN
         IF( 0.0_wp <= zj .AND. zj < 0.5_wp        )   THEN        !  NE quadrant
            !                                                      !             i=I       i=I+1/2
            ze01 = pev(ii  ,ij  )   ;   ze11 = pef(ii  ,ij  )      !   j=J+1/2    V ------- F
            ze00 = pet(ii  ,ij  )   ;   ze10 = peu(ii  ,ij  )      !   j=J        T ------- U
            zi = 2._wp * zi
            zj = 2._wp * zj
         ELSE                                                      !  SE quadrant
            !                                                                    !             i=I       i=I+1/2
            ze01 = pet(ii  ,ij+1)   ;   ze11 = peu(ii  ,ij+1)      !   j=J+1      T ------- U
            ze00 = pev(ii  ,ij  )   ;   ze10 = pef(ii  ,ij  )      !   j=J+1/2    V ------- F
            zi = 2._wp *  zi
            zj = 2._wp * (zj-0.5_wp)
         ENDIF
      ELSE
         IF( 0.0_wp <= zj .AND. zj < 0.5_wp        )   THEN        !  NW quadrant
            !                                                                    !             i=I       i=I+1/2
            ze01 = pef(ii  ,ij  )   ;   ze11 = pev(ii+1,ij)        !   j=J+1/2    F ------- V
            ze00 = peu(ii  ,ij  )   ;   ze10 = pet(ii+1,ij)        !   j=J        U ------- T
            zi = 2._wp * (zi-0.5_wp)
            zj = 2._wp *  zj
         ELSE                                                      !  SW quadrant
            !                                                                    !             i=I+1/2   i=I+1
            ze01 = peu(ii  ,ij+1)   ;   ze11 = pet(ii+1,ij+1)      !   j=J+1      U ------- T
            ze00 = pef(ii  ,ij  )   ;   ze10 = pev(ii+1,ij  )      !   j=J+1/2    F ------- V
            zi = 2._wp * (zi-0.5_wp)
            zj = 2._wp * (zj-0.5_wp)
         ENDIF
      ENDIF
      !
      icb_utl_bilin_e = ( ze01 * (1.-zi) + ze11 * zi ) *     zj    &
         &            + ( ze00 * (1.-zi) + ze10 * zi ) * (1.-zj)
      !
   END FUNCTION icb_utl_bilin_e


   SUBROUTINE icb_utl_add( bergvals, ptvals )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE icb_utl_add           ***
      !!
      !! ** Purpose :   add a new berg to the iceberg list
      !!
      !!----------------------------------------------------------------------
      TYPE(iceberg), INTENT(in)           ::   bergvals
      TYPE(point)  , INTENT(in)           ::   ptvals
      !
      TYPE(iceberg), POINTER ::   new => NULL()
      !!----------------------------------------------------------------------
      !
      new => NULL()
      CALL icb_utl_create( new, bergvals, ptvals )
      CALL icb_utl_insert( new )
      new => NULL()     ! Clear new
      !
   END SUBROUTINE icb_utl_add         


   SUBROUTINE icb_utl_create( berg, bergvals, ptvals )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE icb_utl_create  ***
      !!
      !! ** Purpose :   add a new berg to the iceberg list
      !!
      !!----------------------------------------------------------------------
      TYPE(iceberg), INTENT(in) ::   bergvals
      TYPE(point)  , INTENT(in) ::   ptvals
      TYPE(iceberg), POINTER    ::   berg
      !
      TYPE(point)  , POINTER    ::   pt
      INTEGER                   ::   istat
      !!----------------------------------------------------------------------
      !
      IF( ASSOCIATED(berg) )   CALL ctl_stop( 'icebergs, icb_utl_create: berg already associated' )
      ALLOCATE(berg, STAT=istat)
      IF( istat /= 0 ) CALL ctl_stop( 'failed to allocate iceberg' )
      berg%number(:) = bergvals%number(:)
      berg%mass_scaling = bergvals%mass_scaling
      berg%prev => NULL()
      berg%next => NULL()
      !
      ALLOCATE(pt, STAT=istat)
      IF( istat /= 0 ) CALL ctl_stop( 'failed to allocate first iceberg point' )
      pt = ptvals
      berg%current_point => pt
      !
   END SUBROUTINE icb_utl_create


   SUBROUTINE icb_utl_insert( newberg )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE icb_utl_insert  ***
      !!
      !! ** Purpose :   add a new berg to the iceberg list
      !!
      !!----------------------------------------------------------------------
      TYPE(iceberg), POINTER  ::   newberg
      !
      TYPE(iceberg), POINTER  ::   this, prev, last
      !!----------------------------------------------------------------------
      !
      IF( ASSOCIATED( first_berg ) ) THEN
         last => first_berg
         DO WHILE (ASSOCIATED(last%next))
            last => last%next
         ENDDO
         newberg%prev => last
         last%next    => newberg
         last         => newberg
      ELSE                       ! list is empty so create it
         first_berg => newberg
      ENDIF
      !
   END SUBROUTINE icb_utl_insert


   REAL(wp) FUNCTION icb_utl_yearday(kmon, kday, khr, kmin, ksec)
      !!----------------------------------------------------------------------
      !!                 ***  FUNCTION icb_utl_yearday  ***
      !!
      !! ** Purpose :   
      !!
      ! sga - improved but still only applies to 365 day year, need to do this properly
      !
      !!gm  all these info are already known in daymod, no???
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in)     :: kmon, kday, khr, kmin, ksec
      ! 
      INTEGER, DIMENSION(12)  :: imonths = (/ 0,31,28,31,30,31,30,31,31,30,31,30 /)
      !!----------------------------------------------------------------------
      !
      icb_utl_yearday = REAL( SUM( imonths(1:kmon) ), wp )
      icb_utl_yearday = icb_utl_yearday + REAL(kday-1,wp) + (REAL(khr,wp) + (REAL(kmin,wp) + REAL(ksec,wp)/60.)/60.)/24.
      !
   END FUNCTION icb_utl_yearday

   !!-------------------------------------------------------------------------

   SUBROUTINE icb_utl_delete( first, berg )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE icb_utl_delete  ***
      !!
      !! ** Purpose :   
      !!
      !!----------------------------------------------------------------------
      TYPE(iceberg), POINTER :: first, berg
      !!----------------------------------------------------------------------
      ! Connect neighbors to each other
      IF ( ASSOCIATED(berg%prev) ) THEN
        berg%prev%next => berg%next
      ELSE
        first => berg%next
      ENDIF
      IF (ASSOCIATED(berg%next)) berg%next%prev => berg%prev
      !
      CALL icb_utl_destroy(berg)
      !
   END SUBROUTINE icb_utl_delete


   SUBROUTINE icb_utl_destroy( berg )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE icb_utl_destroy  ***
      !!
      !! ** Purpose :   remove a single iceberg instance
      !!
      !!----------------------------------------------------------------------
      TYPE(iceberg), POINTER :: berg
      !!----------------------------------------------------------------------
      !
      ! Remove any points
      IF( ASSOCIATED( berg%current_point ) )   DEALLOCATE( berg%current_point )
      !
      DEALLOCATE(berg)
      !
   END SUBROUTINE icb_utl_destroy


   SUBROUTINE icb_utl_track( knum, cd_label, kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE icb_utl_track  ***
      !!
      !! ** Purpose :   
      !!
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(nkounts)    :: knum       ! iceberg number
      CHARACTER(len=*)               :: cd_label   ! 
      INTEGER                        :: kt         ! timestep number
      ! 
      TYPE(iceberg), POINTER         :: this
      LOGICAL                        :: match
      INTEGER                        :: k
      !!----------------------------------------------------------------------
      !
      this => first_berg
      DO WHILE( ASSOCIATED(this) )
         match = .TRUE.
         DO k = 1, nkounts
            IF( this%number(k) /= knum(k) ) match = .FALSE.
         END DO
         IF( match )   CALL icb_utl_print_berg(this, kt)
         this => this%next
      END DO
      !
   END SUBROUTINE icb_utl_track


   SUBROUTINE icb_utl_print_berg( berg, kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE icb_utl_print_berg  ***
      !!
      !! ** Purpose :   print one
      !!
      !!----------------------------------------------------------------------
      TYPE(iceberg), POINTER :: berg
      TYPE(point)  , POINTER :: pt
      INTEGER                :: kt      ! timestep number
      !!----------------------------------------------------------------------
      !
      pt => berg%current_point
      WRITE(numicb, 9200) kt, berg%number(1), &
                   pt%xi, pt%yj, pt%lon, pt%lat, pt%uvel, pt%vvel,  &
                   pt%uo, pt%vo, pt%ua, pt%va, pt%ui, pt%vi
      CALL flush( numicb )
 9200 FORMAT(5x,i5,2x,i10,6(2x,2f10.4))
      !
   END SUBROUTINE icb_utl_print_berg


   SUBROUTINE icb_utl_print( cd_label, kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE icb_utl_print  ***
      !!
      !! ** Purpose :   print many
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=*)       :: cd_label
      INTEGER                :: kt             ! timestep number
      ! 
      INTEGER                :: ibergs, inbergs
      TYPE(iceberg), POINTER :: this
      !!----------------------------------------------------------------------
      !
      this => first_berg
      IF( ASSOCIATED(this) ) THEN
         WRITE(numicb,'(a," pe=(",i3,")")' ) cd_label, narea
         WRITE(numicb,'(a8,4x,a6,12x,a5,15x,a7,19x,a3,17x,a5,17x,a5,17x,a5)' )   &
            &         'timestep', 'number', 'xi,yj','lon,lat','u,v','uo,vo','ua,va','ui,vi'
      ENDIF
      DO WHILE( ASSOCIATED(this) )
        CALL icb_utl_print_berg(this, kt)
        this => this%next
      END DO
      ibergs = icb_utl_count()
      inbergs = ibergs
      IF( lk_mpp )   CALL mpp_sum(inbergs)
      IF( ibergs > 0 )   WRITE(numicb,'(a," there are",i5," bergs out of",i6," on PE ",i4)')   &
         &                                  cd_label, ibergs, inbergs, narea
      !
   END SUBROUTINE icb_utl_print


   SUBROUTINE icb_utl_incr()
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE icb_utl_incr  ***
      !!
      !! ** Purpose :   
      !!
      ! Small routine for coping with very large integer values labelling icebergs
      ! num_bergs is a array of integers
      ! the first member is incremented in steps of jpnij starting from narea
      ! this means each iceberg is labelled with a unique number
      ! when this gets to the maximum allowed integer the second and subsequent members are 
      ! used to count how many times the member before cycles
      !!----------------------------------------------------------------------
      INTEGER ::   ii, ibig
      !!----------------------------------------------------------------------

      ibig = HUGE(num_bergs(1))
      IF( ibig-jpnij < num_bergs(1) ) THEN
         num_bergs(1) = narea
         DO ii = 2,nkounts
            IF( num_bergs(ii) == ibig ) THEN
               num_bergs(ii) = 0
               IF( ii == nkounts ) CALL ctl_stop('Sorry, run out of iceberg number space')
            ELSE
               num_bergs(ii) = num_bergs(ii) + 1
               EXIT
            ENDIF
         END DO
      ELSE
         num_bergs(1) = num_bergs(1) + jpnij
      ENDIF
      !
   END SUBROUTINE icb_utl_incr


   INTEGER FUNCTION icb_utl_count()
      !!----------------------------------------------------------------------
      !!                 ***  FUNCTION icb_utl_count  ***
      !!
      !! ** Purpose :   
      !!----------------------------------------------------------------------
      TYPE(iceberg), POINTER :: this
      !!----------------------------------------------------------------------
      !
      icb_utl_count = 0
      this => first_berg
      DO WHILE( ASSOCIATED(this) )
         icb_utl_count = icb_utl_count+1
         this => this%next
      END DO
      !
   END FUNCTION icb_utl_count


   REAL(wp) FUNCTION icb_utl_mass( first, justbits, justbergs )
      !!----------------------------------------------------------------------
      !!                 ***  FUNCTION icb_utl_mass  ***
      !!
      !! ** Purpose :   compute the mass all iceberg, all berg bits or all bergs.
      !!----------------------------------------------------------------------
      TYPE(iceberg)      , POINTER  ::   first
      TYPE(point)        , POINTER  ::   pt
      LOGICAL, INTENT(in), OPTIONAL ::   justbits, justbergs
      !
      TYPE(iceberg), POINTER ::   this
      !!----------------------------------------------------------------------
      icb_utl_mass = 0._wp
      this => first
      !
      IF( PRESENT( justbergs  ) ) THEN
         DO WHILE( ASSOCIATED( this ) )
            pt => this%current_point
            icb_utl_mass = icb_utl_mass + pt%mass         * this%mass_scaling
            this => this%next
         END DO
      ELSEIF( PRESENT(justbits) ) THEN
         DO WHILE( ASSOCIATED( this ) )
            pt => this%current_point
            icb_utl_mass = icb_utl_mass + pt%mass_of_bits * this%mass_scaling
            this => this%next
         END DO
      ELSE
         DO WHILE( ASSOCIATED( this ) )
            pt => this%current_point
            icb_utl_mass = icb_utl_mass + ( pt%mass + pt%mass_of_bits ) * this%mass_scaling
            this => this%next
         END DO
      ENDIF
      !
   END FUNCTION icb_utl_mass


   REAL(wp) FUNCTION icb_utl_heat( first, justbits, justbergs )
      !!----------------------------------------------------------------------
      !!                 ***  FUNCTION icb_utl_heat  ***
      !!
      !! ** Purpose :   compute the heat in all iceberg, all bergies or all bergs.
      !!----------------------------------------------------------------------
      TYPE(iceberg)      , POINTER  ::   first
      LOGICAL, INTENT(in), OPTIONAL ::   justbits, justbergs
      !
      TYPE(iceberg)      , POINTER  ::   this
      TYPE(point)        , POINTER  ::   pt
      !!----------------------------------------------------------------------
      icb_utl_heat = 0._wp
      this => first
      !
      IF( PRESENT( justbergs  ) ) THEN
         DO WHILE( ASSOCIATED( this ) )
            pt => this%current_point
            icb_utl_heat = icb_utl_heat + pt%mass         * this%mass_scaling * pt%heat_density
            this => this%next
         END DO
      ELSEIF( PRESENT(justbits) ) THEN
         DO WHILE( ASSOCIATED( this ) )
            pt => this%current_point
            icb_utl_heat = icb_utl_heat + pt%mass_of_bits * this%mass_scaling * pt%heat_density
            this => this%next
         END DO
      ELSE
         DO WHILE( ASSOCIATED( this ) )
            pt => this%current_point
            icb_utl_heat = icb_utl_heat + ( pt%mass + pt%mass_of_bits ) * this%mass_scaling * pt%heat_density
            this => this%next
         END DO
      ENDIF
      !
   END FUNCTION icb_utl_heat

   !!======================================================================
END MODULE icbutl
