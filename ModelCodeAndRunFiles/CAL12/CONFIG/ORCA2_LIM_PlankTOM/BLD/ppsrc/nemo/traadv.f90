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

MODULE traadv
   !!==============================================================================
   !!                       ***  MODULE  traadv  ***
   !! Ocean active tracers:  advection trend 
   !!==============================================================================
   !! History :  2.0  !  2005-11  (G. Madec)  Original code
   !!            3.3  !  2010-09  (C. Ethe, G. Madec)  merge TRC-TRA + switch from velocity to transport
   !!            4.0  !  2011-06  (G. Madec)  Addition of Mixed Layer Eddy parameterisation
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_adv      : compute ocean tracer advection trend
   !!   tra_adv_ctl  : control the different options of advection scheme
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space and time domain
   USE domvvl          ! variable vertical scale factors
   USE traadv_cen2     ! 2nd order centered scheme (tra_adv_cen2   routine)
   USE traadv_tvd      ! TVD      scheme           (tra_adv_tvd    routine)
   USE traadv_muscl    ! MUSCL    scheme           (tra_adv_muscl  routine)
   USE traadv_muscl2   ! MUSCL2   scheme           (tra_adv_muscl2 routine)
   USE traadv_ubs      ! UBS      scheme           (tra_adv_ubs    routine)
   USE traadv_qck      ! QUICKEST scheme           (tra_adv_qck    routine)
   USE traadv_eiv      ! eddy induced velocity     (tra_adv_eiv    routine)
   USE traadv_mle      ! ML eddy induced velocity  (tra_adv_mle    routine)
   USE cla             ! cross land advection      (cla_traadv     routine)
   USE ldftra_oce      ! lateral diffusion coefficient on tracers
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O module
   USE prtctl          ! Print control
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing
   USE sbc_oce
   USE diaptr          ! Poleward heat transport 


   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_adv        ! routine called by step module
   PUBLIC   tra_adv_init   ! routine called by opa module

   !                              !!* Namelist namtra_adv *
   LOGICAL ::   ln_traadv_cen2     ! 2nd order centered scheme flag
   LOGICAL ::   ln_traadv_tvd      ! TVD scheme flag
   LOGICAL ::   ln_traadv_tvd_zts  ! TVD scheme flag with vertical sub time-stepping
   LOGICAL ::   ln_traadv_muscl    ! MUSCL scheme flag
   LOGICAL ::   ln_traadv_muscl2   ! MUSCL2 scheme flag
   LOGICAL ::   ln_traadv_ubs      ! UBS scheme flag
   LOGICAL ::   ln_traadv_qck      ! QUICKEST scheme flag
   LOGICAL ::   ln_traadv_msc_ups  ! use upstream scheme within muscl


   INTEGER ::   nadv   ! choice of the type of advection scheme

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
   !! $Id: traadv.F90 5147 2015-03-13 10:01:32Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_adv( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv  ***
      !!
      !! ** Purpose :   compute the ocean tracer advection trend.
      !!
      !! ** Method  : - Update (ua,va) with the advection term following nadv
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !
      INTEGER ::   jk   ! dummy loop index
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zun, zvn, zwn
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_adv')
      !
      CALL wrk_alloc( jpi, jpj, jpk, zun, zvn, zwn )
      !                                          ! set time step
      IF( neuler == 0 .AND. kt == nit000 ) THEN     ! at nit000
         r2dtra(:) =  rdttra(:)                          ! = rdtra (restarting with Euler time stepping)
      ELSEIF( kt <= nit000 + 1) THEN                ! at nit000 or nit000+1
         r2dtra(:) = 2._wp * rdttra(:)                   ! = 2 rdttra (leapfrog)
      ENDIF
      !
      IF( nn_cla == 1 .AND. cp_cfg == 'orca' .AND. jp_cfg == 2 )   CALL cla_traadv( kt )       !==  Cross Land Advection  ==! (hor. advection)
      !
      !                                               !==  effective transport  ==!
      DO jk = 1, jpkm1
         zun(:,:,jk) = e2u(:,:) * e3u_0(:,:,jk) * un(:,:,jk)                  ! eulerian transport only
         zvn(:,:,jk) = e1v(:,:) * e3v_0(:,:,jk) * vn(:,:,jk)
         zwn(:,:,jk) = e1t(:,:) * e2t(:,:)      * wn(:,:,jk)
      END DO
      !
      IF( ln_vvl_ztilde .OR. ln_vvl_layer ) THEN
         zun(:,:,:) = zun(:,:,:) + un_td(:,:,:)
         zvn(:,:,:) = zvn(:,:,:) + vn_td(:,:,:)
      ENDIF
      !
      zun(:,:,jpk) = 0._wp                                                     ! no transport trough the bottom
      zvn(:,:,jpk) = 0._wp                                                     ! no transport trough the bottom
      zwn(:,:,jpk) = 0._wp                                                     ! no transport trough the bottom
      !
      IF( lk_traldf_eiv .AND. .NOT. ln_traldf_grif )   &
         &              CALL tra_adv_eiv( kt, nit000, zun, zvn, zwn, 'TRA' )    ! add the eiv transport (if necessary)
      !
      IF( ln_mle    )   CALL tra_adv_mle( kt, nit000, zun, zvn, zwn, 'TRA' )    ! add the mle transport (if necessary)
      !
      CALL iom_put( "uocetr_eff", zun )                                         ! output effective transport      
      CALL iom_put( "vocetr_eff", zvn )
      CALL iom_put( "wocetr_eff", zwn )
      !
      IF( ln_diaptr )   CALL dia_ptr( zvn )                                     ! diagnose the effective MSF 
      !
   
      SELECT CASE ( nadv )                            !==  compute advection trend and add it to general trend  ==!
      CASE ( 1 )   ;    CALL tra_adv_cen2   ( kt, nit000, 'TRA',         zun, zvn, zwn, tsb, tsn, tsa, jpts )   !  2nd order centered
      CASE ( 2 )   ;    CALL tra_adv_tvd    ( kt, nit000, 'TRA', r2dtra, zun, zvn, zwn, tsb, tsn, tsa, jpts )   !  TVD 
      CASE ( 3 )   ;    CALL tra_adv_muscl  ( kt, nit000, 'TRA', r2dtra, zun, zvn, zwn, tsb,      tsa, jpts, ln_traadv_msc_ups )   !  MUSCL 
      CASE ( 4 )   ;    CALL tra_adv_muscl2 ( kt, nit000, 'TRA', r2dtra, zun, zvn, zwn, tsb, tsn, tsa, jpts )   !  MUSCL2 
      CASE ( 5 )   ;    CALL tra_adv_ubs    ( kt, nit000, 'TRA', r2dtra, zun, zvn, zwn, tsb, tsn, tsa, jpts )   !  UBS 
      CASE ( 6 )   ;    CALL tra_adv_qck    ( kt, nit000, 'TRA', r2dtra, zun, zvn, zwn, tsb, tsn, tsa, jpts )   !  QUICKEST 
      CASE ( 7 )   ;    CALL tra_adv_tvd_zts( kt, nit000, 'TRA', r2dtra, zun, zvn, zwn, tsb, tsn, tsa, jpts )   !  TVD ZTS
      !
      CASE (-1 )                                      !==  esopa: test all possibility with control print  ==!
         CALL tra_adv_cen2  ( kt, nit000, 'TRA',         zun, zvn, zwn, tsb, tsn, tsa, jpts )          
         CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' adv0 - Ta: ', mask1=tmask,               &
            &          tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
         CALL tra_adv_tvd   ( kt, nit000, 'TRA', r2dtra, zun, zvn, zwn, tsb, tsn, tsa, jpts )          
         CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' adv1 - Ta: ', mask1=tmask,               &
            &          tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
         CALL tra_adv_muscl ( kt, nit000, 'TRA', r2dtra, zun, zvn, zwn, tsb,      tsa, jpts, ln_traadv_msc_ups )          
         CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' adv3 - Ta: ', mask1=tmask,               &
            &          tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
         CALL tra_adv_muscl2( kt, nit000, 'TRA', r2dtra, zun, zvn, zwn, tsb, tsn, tsa, jpts )          
         CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' adv4 - Ta: ', mask1=tmask,               &
            &          tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
         CALL tra_adv_ubs   ( kt, nit000, 'TRA', r2dtra, zun, zvn, zwn, tsb, tsn, tsa, jpts )          
         CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' adv5 - Ta: ', mask1=tmask,               &
            &          tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
         CALL tra_adv_qck   ( kt, nit000, 'TRA', r2dtra, zun, zvn, zwn, tsb, tsn, tsa, jpts )          
         CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' adv6 - Ta: ', mask1=tmask,               &
            &          tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
      END SELECT
      !
      !                                              ! print mean trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' adv  - Ta: ', mask1=tmask,               &
         &                       tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'tra_adv' )
      !
      CALL wrk_dealloc( jpi, jpj, jpk, zun, zvn, zwn )
      !                                          
   END SUBROUTINE tra_adv


   SUBROUTINE tra_adv_init
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_init  ***
      !!                
      !! ** Purpose :   Control the consistency between namelist options for 
      !!              tracer advection schemes and set nadv
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio
      INTEGER ::   ios                 ! Local integer output status for namelist read
      !!
      NAMELIST/namtra_adv/ ln_traadv_cen2 , ln_traadv_tvd,     &
         &                 ln_traadv_muscl, ln_traadv_muscl2,  &
         &                 ln_traadv_ubs  , ln_traadv_qck,     &
         &                 ln_traadv_msc_ups, ln_traadv_tvd_zts
      !!----------------------------------------------------------------------

      REWIND( numnam_ref )              ! Namelist namtra_adv in reference namelist : Tracer advection scheme
      READ  ( numnam_ref, namtra_adv, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtra_adv in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namtra_adv in configuration namelist : Tracer advection scheme
      READ  ( numnam_cfg, namtra_adv, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtra_adv in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namtra_adv )

      IF(lwp) THEN                    ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'tra_adv_init : choice/control of the tracer advection scheme'
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtra_adv : chose a advection scheme for tracers'
         WRITE(numout,*) '      2nd order advection scheme     ln_traadv_cen2    = ', ln_traadv_cen2
         WRITE(numout,*) '      TVD advection scheme           ln_traadv_tvd     = ', ln_traadv_tvd
         WRITE(numout,*) '      MUSCL  advection scheme        ln_traadv_muscl   = ', ln_traadv_muscl
         WRITE(numout,*) '      MUSCL2 advection scheme        ln_traadv_muscl2  = ', ln_traadv_muscl2
         WRITE(numout,*) '      UBS    advection scheme        ln_traadv_ubs     = ', ln_traadv_ubs
         WRITE(numout,*) '      QUICKEST advection scheme      ln_traadv_qck     = ', ln_traadv_qck
         WRITE(numout,*) '      upstream scheme within muscl   ln_traadv_msc_ups = ', ln_traadv_msc_ups
         WRITE(numout,*) '      TVD advection scheme with zts  ln_traadv_tvd_zts = ', ln_traadv_tvd_zts
      ENDIF

      ioptio = 0                      ! Parameter control
      IF( ln_traadv_cen2   )   ioptio = ioptio + 1
      IF( ln_traadv_tvd    )   ioptio = ioptio + 1
      IF( ln_traadv_muscl  )   ioptio = ioptio + 1
      IF( ln_traadv_muscl2 )   ioptio = ioptio + 1
      IF( ln_traadv_ubs    )   ioptio = ioptio + 1
      IF( ln_traadv_qck    )   ioptio = ioptio + 1
      IF( ln_traadv_tvd_zts)   ioptio = ioptio + 1
      IF( lk_esopa         )   ioptio =          1

      IF( ( ln_traadv_muscl .OR. ln_traadv_muscl2 .OR. ln_traadv_ubs .OR. ln_traadv_qck .OR. ln_traadv_tvd_zts )   &
         .AND. ln_isfcav )   CALL ctl_stop( 'Only traadv_cen2 and traadv_tvd is compatible with ice shelf cavity')

      IF( ioptio /= 1 )   CALL ctl_stop( 'Choose ONE advection scheme in namelist namtra_adv' )

      !                              ! Set nadv
      IF( ln_traadv_cen2   )   nadv =  1
      IF( ln_traadv_tvd    )   nadv =  2
      IF( ln_traadv_muscl  )   nadv =  3
      IF( ln_traadv_muscl2 )   nadv =  4
      IF( ln_traadv_ubs    )   nadv =  5
      IF( ln_traadv_qck    )   nadv =  6
      IF( ln_traadv_tvd_zts)   nadv =  7
      IF( lk_esopa         )   nadv = -1

      IF(lwp) THEN                   ! Print the choice
         WRITE(numout,*)
         IF( nadv ==  1 )   WRITE(numout,*) '         2nd order scheme is used'
         IF( nadv ==  2 )   WRITE(numout,*) '         TVD       scheme is used'
         IF( nadv ==  3 )   WRITE(numout,*) '         MUSCL     scheme is used'
         IF( nadv ==  4 )   WRITE(numout,*) '         MUSCL2    scheme is used'
         IF( nadv ==  5 )   WRITE(numout,*) '         UBS       scheme is used'
         IF( nadv ==  6 )   WRITE(numout,*) '         QUICKEST  scheme is used'
         IF( nadv ==  7 )   WRITE(numout,*) '         TVD ZTS   scheme is used'
         IF( nadv == -1 )   WRITE(numout,*) '         esopa test: use all advection scheme'
      ENDIF
      !
      CALL tra_adv_mle_init          ! initialisation of the Mixed Layer Eddy parametrisation (MLE)
      !
   END SUBROUTINE tra_adv_init

  !!======================================================================
END MODULE traadv
