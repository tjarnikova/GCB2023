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

MODULE sbcice_lim_2
   !!======================================================================
   !!                       ***  MODULE  sbcice_lim_2  ***
   !! Surface module :  update surface ocean boundary condition over ice covered area using LIM sea-ice model
   !! Sea-Ice model  :  LIM-2 Sea ice model time-stepping
   !!======================================================================
   !! History :  1.0   !  06-2006  (G. Madec)  from icestp_2.F90
   !!            3.0   !  08-2008  (S. Masson, E. .... ) coupled interface
   !!            3.3   !  05-2009  (G.Garric) addition of the lim2_evp case
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_lim2' :                                    LIM-2 sea-ice model
   !!----------------------------------------------------------------------
   !!   sbc_ice_lim_2   : sea-ice model time-stepping and update ocean sbc over ice-covered area
   !!----------------------------------------------------------------------
   USE oce              ! ocean dynamics and tracers
   USE dom_oce          ! ocean space and time domain
   USE ice_2
   USE par_ice_2
   USE iceini_2
   USE dom_ice_2

   USE sbc_oce          ! Surface boundary condition: ocean fields
   USE sbc_ice          ! Surface boundary condition: ice   fields
   USE sbcblk_core      ! Surface boundary condition: CORE bulk
   USE sbcblk_clio      ! Surface boundary condition: CLIO bulk
   USE sbccpl           ! Surface boundary condition: coupled interface
   USE albedo

   USE phycst           ! Define parameters for the routines
   USE eosbn2           ! equation of state
   USE limdyn_2
   USE limtrp_2
   USE limdmp_2
   USE limthd_2
   USE limsbc_2         ! sea surface boundary condition
   USE limdia_2
   USE limwri_2
   USE limrst_2

   USE c1d              ! 1D vertical configuration

   USE lbclnk           ! lateral boundary condition - MPP link
   USE lib_mpp          ! MPP library
   USE wrk_nemo         ! work arrays
   USE iom              ! I/O manager library
   USE in_out_manager   ! I/O manager
   USE prtctl           ! Print control



   IMPLICIT NONE
   PRIVATE

   PUBLIC sbc_ice_lim_2 ! routine called by sbcmod.F90

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
   !! $Id: sbcice_lim_2.F90 5540 2015-07-02 15:11:23Z jchanut $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_ice_lim_2( kt, ksbc )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_ice_lim_2  ***
      !!                   
      !! ** Purpose :   update the ocean surface boundary condition via the 
      !!                Louvain la Neuve Sea Ice Model time stepping 
      !!
      !! ** Method  :   ice model time stepping
      !!              - call the ice dynamics routine 
      !!              - call the ice advection/diffusion routine 
      !!              - call the ice thermodynamics routine 
      !!              - call the routine that computes mass and 
      !!                heat fluxes at the ice/ocean interface
      !!              - save the outputs 
      !!              - save the outputs for restart when necessary
      !!
      !! ** Action  : - time evolution of the LIM sea-ice model
      !!              - update all sbc variables below sea-ice:
      !!                utau, vtau, taum, wndm, qns , qsr, emp , sfx 
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt      ! ocean time step
      INTEGER, INTENT(in) ::   ksbc    ! type of sbc ( =3 CLIO bulk ; =4 CORE bulk ; =5 coupled )
      !!
      INTEGER  ::   ji, jj   ! dummy loop indices
      REAL(wp), DIMENSION(:,:,:), POINTER :: zalb_os   ! ice albedo under overcast sky
      REAL(wp), DIMENSION(:,:,:), POINTER :: zalb_cs   ! ice albedo under clear sky
      REAL(wp), DIMENSION(:,:,:), POINTER :: zalb_ice  ! mean ice albedo
      REAL(wp), DIMENSION(:,:,:), POINTER :: zsist     ! ice surface temperature (K)
      REAL(wp), DIMENSION(:,:  ), POINTER :: zutau_ice, zvtau_ice 
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'sbc_ice_lim_2 : update ocean surface boudary condition' 
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~   via Louvain la Neuve Ice Model (LIM) time stepping'
         !
         CALL ice_init_2
         !
      ENDIF

      !                                        !----------------------!
      IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN     !  Ice time-step only  !
         !                                     !----------------------!

         CALL wrk_alloc( jpi,jpj  , zutau_ice, zvtau_ice)
         CALL wrk_alloc( jpi,jpj,1, zalb_os, zalb_cs, zalb_ice, zsist )

         !  Bulk Formulea !
         !----------------!
         ! ... mean surface ocean current at ice dynamics point
         SELECT CASE( cp_ice_msh )
         CASE( 'I' )                  !== B-grid ice dynamics :   I-point (i.e. F-point with sea-ice indexation)
            DO jj = 2, jpj
               DO ji = 2, jpi   ! NO vector opt. possible
                  u_oce(ji,jj) = 0.5_wp * ( ssu_m(ji-1,jj  ) * umask(ji-1,jj  ,1) &
                     &                    + ssu_m(ji-1,jj-1) * umask(ji-1,jj-1,1) ) * tmu(ji,jj)
                  v_oce(ji,jj) = 0.5_wp * ( ssv_m(ji  ,jj-1) * vmask(ji  ,jj-1,1) &
                     &                    + ssv_m(ji-1,jj-1) * vmask(ji-1,jj-1,1) ) * tmu(ji,jj)
               END DO
            END DO
            CALL lbc_lnk( u_oce, 'I', -1. )   ! I-point (i.e. F-point with ice indices)
            CALL lbc_lnk( v_oce, 'I', -1. )   ! I-point (i.e. F-point with ice indices)
            !
         CASE( 'C' )                  !== C-grid ice dynamics :   U & V-points (same as ocean)
            u_oce(:,:) = ssu_m(:,:) * umask(:,:,1)                     ! mean surface ocean current at ice velocity point
            v_oce(:,:) = ssv_m(:,:) * vmask(:,:,1)
            !
         END SELECT

         ! ... masked sea surface freezing temperature [Kelvin] (set to rt0 over land)
         CALL eos_fzp( sss_m(:,:), tfu(:,:) )
         tfu(:,:) = tfu(:,:) + rt0

         zsist (:,:,1) = sist (:,:) + rt0 * ( 1. - tmask(:,:,1) )

         ! Ice albedo

         CALL albedo_ice( zsist, reshape( hicif, (/jpi,jpj,1/) ), &
                                 reshape( hsnif, (/jpi,jpj,1/) ), &
                          zalb_cs, zalb_os )

         SELECT CASE( ksbc )
         CASE( jp_core , jp_purecpl )   ! CORE and COUPLED bulk formulations

            ! albedo depends on cloud fraction because of non-linear spectral effects
            zalb_ice(:,:,:) = ( 1. - cldf_ice ) * zalb_cs(:,:,:) + cldf_ice * zalb_os(:,:,:)
            ! In CLIO the cloud fraction is read in the climatology and the all-sky albedo 
            ! (zalb_ice) is computed within the bulk routine

         END SELECT

         ! ... Sea-ice surface boundary conditions output from bulk formulae :
         !     - utau_ice   ! surface ice stress i-component (I-point)   [N/m2]
         !     - vtau_ice   ! surface ice stress j-component (I-point)   [N/m2]
         !     - qns_ice    ! non solar heat flux over ice   (T-point)   [W/m2]
         !     - qsr_ice    !     solar heat flux over ice   (T-point)   [W/m2]
         !     - qla_ice    ! latent    heat flux over ice   (T-point)   [W/m2]
         !     - dqns_ice   ! non solar heat sensistivity    (T-point)   [W/m2]
         !     - dqla_ice   ! latent    heat sensistivity    (T-point)   [W/m2]
         !     - tprecip    ! total precipitation            (T-point)   [Kg/m2/s]
         !     - sprecip    ! solid precipitation            (T-point)   [Kg/m2/s]
         !     - fr1_i0     ! 1sr fraction of qsr penetration in ice     [%]
         !     - fr2_i0     ! 2nd fraction of qsr penetration in ice     [%]
         !
         SELECT CASE( ksbc )
         CASE( jp_clio )           ! CLIO bulk formulation
!           CALL blk_ice_clio( zsist, zalb_cs    , zalb_os    , zalb_ice   ,            &
!              &                      utau_ice   , vtau_ice   , qns_ice    , qsr_ice,   &
!              &                      qla_ice    , dqns_ice   , dqla_ice   ,            &
!              &                      tprecip    , sprecip    ,                         &
!              &                      fr1_i0     , fr2_i0     , cp_ice_msh , jpl  )
            CALL blk_ice_clio_tau
            CALL blk_ice_clio_flx( zsist, zalb_cs, zalb_os, zalb_ice )

         CASE( jp_core )           ! CORE bulk formulation
            CALL blk_ice_core_tau
            CALL blk_ice_core_flx( zsist, zalb_ice )

         CASE( jp_purecpl )            ! Coupled formulation : atmosphere-ice stress only (fluxes provided after ice dynamics)
            CALL sbc_cpl_ice_tau( utau_ice , vtau_ice )
         END SELECT
         
         IF( ln_mixcpl) THEN
            CALL sbc_cpl_ice_tau( zutau_ice , zvtau_ice )
            utau_ice(:,:) = utau_ice(:,:) * xcplmask(:,:,0) + zutau_ice(:,:) * ( 1. - xcplmask(:,:,0) )
            vtau_ice(:,:) = vtau_ice(:,:) * xcplmask(:,:,0) + zvtau_ice(:,:) * ( 1. - xcplmask(:,:,0) )
         ENDIF

         CALL iom_put( 'utau_ice', utau_ice )     ! Wind stress over ice along i-axis at I-point
         CALL iom_put( 'vtau_ice', vtau_ice )     ! Wind stress over ice along j-axis at I-point

         IF(ln_ctl) THEN         ! print mean trends (used for debugging)
            CALL prt_ctl_info( 'Ice Forcings ' )
            CALL prt_ctl( tab2d_1=tprecip ,clinfo1=' sbc_ice_lim: precip  : ', tab2d_2=sprecip , clinfo2=' Snow    : ' )
            CALL prt_ctl( tab2d_1=utau_ice,clinfo1=' sbc_ice_lim: utau_ice: ', tab2d_2=vtau_ice, clinfo2=' vtau_ice: ' )
            CALL prt_ctl( tab2d_1=sst_m   ,clinfo1=' sbc_ice_lim: sst     : ', tab2d_2=sss_m   , clinfo2=' sss     : ' )
            CALL prt_ctl( tab2d_1=u_oce   ,clinfo1=' sbc_ice_lim: u_io    : ', tab2d_2=v_oce   , clinfo2=' v_io    : ' )
            CALL prt_ctl( tab2d_1=hsnif   ,clinfo1=' sbc_ice_lim: hsnif  1: ', tab2d_2=hicif   , clinfo2=' hicif   : ' )
            CALL prt_ctl( tab2d_1=frld    ,clinfo1=' sbc_ice_lim: frld   1: ', tab2d_2=sist    , clinfo2=' sist    : ' )
         ENDIF

         ! ---------------- !
         !  Ice model step  !
         ! ---------------- !
         numit = numit + nn_fsbc                           ! Ice model time step

                           CALL lim_rst_opn_2  ( kt )  ! Open Ice restart file
         IF( .NOT. lk_c1d ) THEN                       ! Ice dynamics & transport (except in 1D case)
                           CALL lim_dyn_2      ( kt )      ! Ice dynamics    ( rheology/dynamics )
                           CALL lim_trp_2      ( kt )      ! Ice transport   ( Advection/diffusion )
           IF( ln_limdmp ) CALL lim_dmp_2      ( kt )      ! Ice damping 
         END IF
         !                                             ! Ice surface fluxes in coupled mode 
         IF( ln_cpl ) THEN   ! pure coupled and mixed forced-coupled configurations
            a_i(:,:,1)=fr_i
            CALL sbc_cpl_ice_flx( frld,                                              &
            !                                optional arguments, used only in 'mixed oce-ice' case
            &                                             palbi=zalb_ice, psst=sst_m, pist=zsist )
            sprecip(:,:) = - emp_ice(:,:)   ! Ugly patch, WARNING, in coupled mode, sublimation included in snow (parsub = 0.)
         ENDIF
                           CALL lim_thd_2      ( kt )      ! Ice thermodynamics 
                           CALL lim_sbc_flx_2  ( kt )      ! update surface ocean mass, heat & salt fluxes 

         IF(  .NOT. lk_mpp )THEN
            IF( MOD( kt+nn_fsbc-1, ninfo ) == 0 .OR. ntmoy == 1 )   &
            &              CALL lim_dia_2      ( kt )      ! Ice Diagnostics
         ENDIF
         IF( lrst_ice  )   CALL lim_rst_write_2( kt )      ! Ice restart file
         !
         !
         CALL wrk_dealloc( jpi,jpj  , zutau_ice, zvtau_ice)
         CALL wrk_dealloc( jpi,jpj,1, zalb_os, zalb_cs, zalb_ice, zsist )
         !
      ENDIF                                    ! End sea-ice time step only
      !
      !                                        !--------------------------!
      !                                        !  at all ocean time step  !
      !                                        !--------------------------!
      !                                               
      !                                              ! Update surface ocean stresses (only in ice-dynamic case)
      !                                                   ! otherwise the atm.-ocean stresses are used everywhere
      IF( ln_limdyn    )   CALL lim_sbc_tau_2( kt, ub(:,:,1), vb(:,:,1) )  ! using before instantaneous surf. currents
      !
   END SUBROUTINE sbc_ice_lim_2


   !!======================================================================
END MODULE sbcice_lim_2
