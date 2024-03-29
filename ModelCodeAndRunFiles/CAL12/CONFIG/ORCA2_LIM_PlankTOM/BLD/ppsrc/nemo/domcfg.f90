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

MODULE domcfg
   !!==============================================================================
   !!                       ***  MODULE domcfg   ***
   !! Ocean initialization : domain configuration initialization
   !!==============================================================================
   !! History :  1.0  ! 2003-09  (G. Madec)  Original code
   !!            3.2  ! 2009-07  (R. Benshila) Suppression of rigid-lid option
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dom_cfg        : initialize the domain configuration
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distributed memory computing library
   USE timing          ! Timing
   USE c1d             ! 1D configuration
   USE domc1d          ! 1D configuration: column location

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dom_cfg    ! called by opa.F90

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.2 , LODYC-IPSL  (2009)
   !! $Id: domcfg.F90 4667 2014-06-12 08:47:52Z cbricaud $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dom_cfg
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dom_cfg  ***
      !!                    
      !! ** Purpose :   set the domain configuration
      !!
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dom_cfg')
      !
      IF(lwp) THEN                   ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'dom_cfg : set the ocean configuration'
         WRITE(numout,*) '~~~~~~~ '
         WRITE(numout,*) '   ocean model configuration used :   cp_cfg = ', cp_cfg, ' jp_cfg = ', jp_cfg
         !
         WRITE(numout,*) '   global domain lateral boundaries'
         !
         IF( jperio == 0 )   WRITE(numout,*) '      jperio= 0, closed'
         IF( jperio == 1 )   WRITE(numout,*) '      jperio= 1, cyclic east-west'
         IF( jperio == 2 )   WRITE(numout,*) '      jperio= 2, equatorial symmetric'
         IF( jperio == 3 )   WRITE(numout,*) '      jperio= 3, north fold with T-point pivot'
         IF( jperio == 4 )   WRITE(numout,*) '      jperio= 4, cyclic east-west and north fold with T-point pivot'
         IF( jperio == 5 )   WRITE(numout,*) '      jperio= 5, north fold with F-point pivot'
         IF( jperio == 6 )   WRITE(numout,*) '      jperio= 6, cyclic east-west and north fold with F-point pivot'
      ENDIF
      !
      IF( jperio <  0 .OR. jperio > 6 )   CALL ctl_stop( 'jperio is out of range' )
      !
      CALL dom_glo                   ! global domain versus zoom and/or local domain
      !
      IF( nn_timing == 1 )  CALL timing_stop('dom_cfg')
      !
   END SUBROUTINE dom_cfg


   SUBROUTINE dom_glo
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_glo  ***
      !!
      !! ** Purpose :   initialization for global domain, zoom and local domain
      !!
      !! ** Method  :   
      !!
      !! ** Action  : - mig  , mjg : 
      !!              - mi0  , mi1   :
      !!              - mj0, , mj1   :
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj   ! dummy loop argument
      !!----------------------------------------------------------------------
      !                              ! recalculate jpizoom/jpjzoom given lat/lon
      IF( lk_c1d .AND. ln_c1d_locpt )  CALL dom_c1d( rn_lat1d, rn_lon1d )
      !
      !                        ! ============== !
      !                        !  Local domain  ! 
      !                        ! ============== !
      DO ji = 1, jpi                 ! local domain indices ==> data domain indices
        mig(ji) = ji + jpizoom - 1 + nimpp - 1
      END DO
      DO jj = 1, jpj
        mjg(jj) = jj + jpjzoom - 1 + njmpp - 1
      END DO
      !
      !                              ! data domain indices ==> local domain indices
      !                                   ! (return (m.0,m.1)=(1,0) if data domain gridpoint is to the west/south of the 
      !                                   !local domain, or (m.0,m.1)=(jp.+1,jp.) to the east/north of local domain. 
      DO ji = 1, jpidta
        mi0(ji) = MAX( 1, MIN( ji - jpizoom + 1 - nimpp + 1, jpi+1 ) )
        mi1(ji) = MAX( 0, MIN( ji - jpizoom + 1 - nimpp + 1, jpi   ) )
      END DO
      DO jj = 1, jpjdta
        mj0(jj) = MAX( 1, MIN( jj - jpjzoom + 1 - njmpp + 1, jpj+1 ) )
        mj1(jj) = MAX( 0, MIN( jj - jpjzoom + 1 - njmpp + 1, jpj   ) )
      END DO
      IF(lwp) THEN                   ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'dom_glo : domain: data / local '
         WRITE(numout,*) '~~~~~~~ '
         WRITE(numout,*) '          data input domain    : jpidta = ', jpidta,   &
            &                                            ' jpjdta = ', jpjdta, ' jpkdta = ', jpkdta
         WRITE(numout,*) '          global or zoom domain: jpiglo = ', jpiglo,   &
            &                                            ' jpjglo = ', jpjglo, ' jpk    = ', jpk
         WRITE(numout,*) '          local domain         : jpi    = ', jpi   ,   &
            &                                            ' jpj    = ', jpj   , ' jpk    = ', jpk
         WRITE(numout,*)
         WRITE(numout,*) '          south-west indices    jpizoom = ', jpizoom,   &
            &                                           ' jpjzoom = ', jpjzoom
         WRITE(numout,*)
         WRITE(numout,*) '          conversion local  ==> data i-index domain'
         WRITE(numout,25)              (mig(ji),ji = 1,jpi)
         WRITE(numout,*)
         WRITE(numout,*) '          conversion data   ==> local  i-index domain'
         WRITE(numout,*) '             starting index'
         WRITE(numout,25)              (mi0(ji),ji = 1,jpidta)
         WRITE(numout,*) '             ending index'
         WRITE(numout,25)              (mi1(ji),ji = 1,jpidta)
         WRITE(numout,*)
         WRITE(numout,*) '          conversion local  ==> data j-index domain'
         WRITE(numout,25)              (mjg(jj),jj = 1,jpj)
         WRITE(numout,*)
         WRITE(numout,*) '          conversion data  ==> local j-index domain'
         WRITE(numout,*) '             starting index'
         WRITE(numout,25)              (mj0(jj),jj = 1,jpjdta)
         WRITE(numout,*) '             ending index'
         WRITE(numout,25)              (mj1(jj),jj = 1,jpjdta)
      ENDIF
 25   FORMAT( 100(10x,19i4,/) )

      !                        ! ============== !
      !                        !  Zoom domain   !
      !                        ! ============== !
      !                              ! zoom control
      IF( jpiglo + jpizoom - 1  >  jpidta .OR.   &
          jpjglo + jpjzoom - 1  >  jpjdta      ) &
          &   CALL ctl_stop( ' global or zoom domain exceed the data domain ! ' )

      !                              ! set zoom flag
      IF( jpiglo < jpidta .OR. jpjglo < jpjdta )   lzoom = .TRUE.

      !                              ! set zoom type flags
      IF( lzoom .AND. jpizoom /= 1 )   lzoom_w = .TRUE.                     ! 
      IF( lzoom .AND. jpjzoom /= 1 )   lzoom_s = .TRUE.
      IF( lzoom .AND. jpiglo + jpizoom -1 /= jpidta )   lzoom_e = .TRUE.
      IF( lzoom .AND. jpjglo + jpjzoom -1 /= jpjdta )   lzoom_n = .TRUE.
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '          zoom flags : '
         WRITE(numout,*) '             lzoom   = ', lzoom  , ' (T = zoom, F = global )'
         WRITE(numout,*) '             lzoom_e = ', lzoom_e, ' (T = forced closed east  boundary)'
         WRITE(numout,*) '             lzoom_w = ', lzoom_w, ' (T = forced closed west  boundary)'
         WRITE(numout,*) '             lzoom_s = ', lzoom_s, ' (T = forced closed South boundary)'
         WRITE(numout,*) '             lzoom_n = ', lzoom_n, ' (T = forced closed North boundary)'
      ENDIF
      IF(  ( lzoom_e .OR. lzoom_w )  .AND.  ( jperio == 1 .OR. jperio == 4 .OR. jperio == 6 )  )   &
           &   CALL ctl_stop( ' Your zoom choice is inconsistent with east-west cyclic boundary condition' )
      IF(  lzoom_n  .AND.  ( 3 <= jperio .AND. jperio <= 6 )  )   &
           &   CALL ctl_stop( ' Your zoom choice is inconsistent with North fold boundary condition' )

      !                              ! Pre-defined arctic/antarctic zoom of ORCA configuration flag
      IF( cp_cfg == "orca" ) THEN
         SELECT CASE ( jp_cfg )
         CASE ( 2 )                               !  ORCA_R2 configuration
            IF(  cp_cfz == "arctic"    .AND. jpiglo  == 142    .AND. jpjglo  ==  53 .AND.   &
               & jpizoom ==  21    .AND. jpjzoom ==  97         )   THEN
              IF(lwp) WRITE(numout,*) '          ORCA configuration: arctic zoom '
            ENDIF
            IF(  cp_cfz == "antarctic" .AND. jpiglo  == jpidta .AND. jpjglo  ==  50 .AND.   &
               & jpizoom ==   1    .AND. jpjzoom ==   1         )   THEN
              IF(lwp) WRITE(numout,*) '          ORCA configuration: antarctic zoom '
            ENDIF
            !                             
         CASE ( 05 )                              !  ORCA_R05 configuration
            IF(    cp_cfz == "arctic"    .AND. jpiglo  == 562    .AND. jpjglo  == 202 .AND.   &
               & jpizoom ==  81    .AND. jpjzoom == 301         )   THEN
              IF(lwp) WRITE(numout,*) '          ORCA configuration: arctic zoom '
            ENDIF
            IF(    cp_cfz == "antarctic" .AND. jpiglo  == jpidta .AND. jpjglo  == 187 .AND.   &
               & jpizoom ==   1    .AND. jpjzoom ==   1         )   THEN
              IF(lwp) WRITE(numout,*) '          ORCA configuration: antarctic zoom '
            ENDIF
         END SELECT
         !
      ENDIF
      !
   END SUBROUTINE dom_glo

   !!======================================================================
END MODULE domcfg
