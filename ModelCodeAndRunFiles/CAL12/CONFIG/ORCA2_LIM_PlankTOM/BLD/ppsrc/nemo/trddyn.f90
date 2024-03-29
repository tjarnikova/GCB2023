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

MODULE trddyn
   !!======================================================================
   !!                       ***  MODULE  trddyn  ***
   !! Ocean diagnostics:  ocean dynamic trends
   !!=====================================================================
   !! History :  3.5  !  2012-02  (G. Madec) creation from trdmod: split DYN and TRA trends
   !!                                        and manage  3D trends output for U, V, and KE
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   trd_dyn       : manage the type of momentum trend diagnostics (3D I/O, domain averaged, KE)
   !!   trd_dyn_iom   : output 3D momentum and/or tracer trends using IOM
   !!   trd_dyn_init  : initialization step
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain variables
   USE zdf_oce        ! ocean vertical physics variables
   USE trd_oce        ! trends: ocean variables
   USE zdfbfr         ! bottom friction
   USE sbc_oce        ! surface boundary condition: ocean
   USE phycst         ! physical constants
   USE trdken         ! trends: Kinetic ENergy 
   USE trdglo         ! trends: global domain averaged
   USE trdvor         ! trends: vertical averaged vorticity 
   USE trdmxl         ! trends: mixed layer averaged 
   USE in_out_manager ! I/O manager
   USE lbclnk         ! lateral boundary condition 
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE wrk_nemo       ! Memory allocation

   IMPLICIT NONE
   PRIVATE

   PUBLIC trd_dyn        ! called by all dynXX modules

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
   !! $Id: trddyn.F90 5215 2015-04-15 16:11:56Z nicolasmartin $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trd_dyn( putrd, pvtrd, ktrd, kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_mod  ***
      !! 
      !! ** Purpose :   Dispatch momentum trend computation, e.g. 3D output, 
      !!              integral constraints, barotropic vorticity, kinetic enrgy, 
      !!              and/or mixed layer budget.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   putrd, pvtrd   ! U and V trends 
      INTEGER                   , INTENT(in   ) ::   ktrd           ! trend index
      INTEGER                   , INTENT(in   ) ::   kt             ! time step
      !!----------------------------------------------------------------------
      !
      putrd(:,:,:) = putrd(:,:,:) * umask(:,:,:)                       ! mask the trends
      pvtrd(:,:,:) = pvtrd(:,:,:) * vmask(:,:,:)
      !

!!gm NB : here a lbc_lnk should probably be added

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !   3D output of momentum and/or tracers trends using IOM interface
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( ln_dyn_trd )   CALL trd_dyn_iom( putrd, pvtrd, ktrd, kt )
         
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !  Integral Constraints Properties for momentum and/or tracers trends
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( ln_glo_trd )   CALL trd_glo( putrd, pvtrd, ktrd, 'DYN', kt )

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !  Kinetic Energy trends
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF( ln_KE_trd  )   CALL trd_ken( putrd, pvtrd, ktrd, kt )

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !  Vorticity trends
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF( ln_vor_trd )   CALL trd_vor( putrd, pvtrd, ktrd, kt )

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !  Mixed layer trends for active tracers
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!gm      IF( ln_dyn_mxl )   CALL trd_mxl_dyn   
      !
   END SUBROUTINE trd_dyn


   SUBROUTINE trd_dyn_iom( putrd, pvtrd, ktrd, kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_dyn_iom  ***
      !! 
      !! ** Purpose :   output 3D trends using IOM
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   putrd, pvtrd   ! U and V trends
      INTEGER                   , INTENT(in   ) ::   ktrd           ! trend index
      INTEGER                   , INTENT(in   ) ::   kt             ! time step
      !
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      INTEGER ::   ikbu, ikbv   ! local integers
      REAL(wp), POINTER, DIMENSION(:,:)   ::   z2dx, z2dy   ! 2D workspace 
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   z3dx, z3dy   ! 3D workspace 
      !!----------------------------------------------------------------------
      !
      SELECT CASE( ktrd )
      CASE( jpdyn_hpg )   ;   CALL iom_put( "utrd_hpg", putrd )    ! hydrostatic pressure gradient
                              CALL iom_put( "vtrd_hpg", pvtrd )
      CASE( jpdyn_spg )   ;   CALL iom_put( "utrd_spg", putrd )    ! surface pressure gradient
                              CALL iom_put( "vtrd_spg", pvtrd )
      CASE( jpdyn_spgexp );   CALL iom_put( "utrd_spgexp", putrd ) ! surface pressure gradient (explicit)
                              CALL iom_put( "vtrd_spgexp", pvtrd )
      CASE( jpdyn_spgflt );   CALL iom_put( "utrd_spgflt", putrd ) ! surface pressure gradient (filtered)
                              CALL iom_put( "vtrd_spgflt", pvtrd )
      CASE( jpdyn_pvo )   ;   CALL iom_put( "utrd_pvo", putrd )    ! planetary vorticity
                              CALL iom_put( "vtrd_pvo", pvtrd )
      CASE( jpdyn_rvo )   ;   CALL iom_put( "utrd_rvo", putrd )    ! relative  vorticity     (or metric term)
                              CALL iom_put( "vtrd_rvo", pvtrd )
      CASE( jpdyn_keg )   ;   CALL iom_put( "utrd_keg", putrd )    ! Kinetic Energy gradient (or had)
                              CALL iom_put( "vtrd_keg", pvtrd )
                              CALL wrk_alloc( jpi, jpj, jpk, z3dx, z3dy )
                              z3dx(:,:,:) = 0._wp                  ! U.dxU & V.dyV (approximation)
                              z3dy(:,:,:) = 0._wp
                              DO jk = 1, jpkm1					   ! no mask as un,vn are masked
                                 DO jj = 2, jpjm1
                                    DO ji = 2, jpim1
                                       z3dx(ji,jj,jk) = un(ji,jj,jk) * ( un(ji+1,jj,jk) - un(ji-1,jj,jk) ) / ( 2._wp * e1u(ji,jj) )
                                       z3dy(ji,jj,jk) = vn(ji,jj,jk) * ( vn(ji,jj+1,jk) - vn(ji,jj-1,jk) ) / ( 2._wp * e2v(ji,jj) )
                                    END DO
                                 END DO
                              END DO
                              CALL lbc_lnk( z3dx, 'U', -1. )
                              CALL lbc_lnk( z3dy, 'V', -1. )
                              CALL iom_put( "utrd_udx", z3dx  )
                              CALL iom_put( "vtrd_vdy", z3dy  )
                              CALL wrk_dealloc( jpi, jpj, jpk, z3dx, z3dy )
      CASE( jpdyn_zad )   ;   CALL iom_put( "utrd_zad", putrd )    ! vertical   advection
                              CALL iom_put( "vtrd_zad", pvtrd )
      CASE( jpdyn_ldf )   ;   CALL iom_put( "utrd_ldf", putrd )    ! lateral diffusion
                              CALL iom_put( "vtrd_ldf", pvtrd )
      CASE( jpdyn_zdf )   ;   CALL iom_put( "utrd_zdf", putrd )    ! vertical diffusion 
                              CALL iom_put( "vtrd_zdf", pvtrd )
                              !                                    ! wind stress trends
                              CALL wrk_alloc( jpi, jpj, z2dx, z2dy )
                              z2dx(:,:) = ( utau_b(:,:) + utau(:,:) ) / ( e3u_0(:,:,1) * rau0 )
                              z2dy(:,:) = ( vtau_b(:,:) + vtau(:,:) ) / ( e3v_0(:,:,1) * rau0 )
                              CALL iom_put( "utrd_tau", z2dx )
                              CALL iom_put( "vtrd_tau", z2dy )
                              CALL wrk_dealloc( jpi, jpj, z2dx, z2dy )
      CASE( jpdyn_bfr )       ! called if ln_bfrimp=T
                              CALL iom_put( "utrd_bfr", putrd )    ! bottom friction (explicit case)
                              CALL iom_put( "vtrd_bfr", pvtrd )
      CASE( jpdyn_atf )   ;   CALL iom_put( "utrd_atf", putrd )        ! asselin filter trends 
                              CALL iom_put( "vtrd_atf", pvtrd )
      CASE( jpdyn_bfri )  ;   IF( ln_bfrimp ) THEN                     ! bottom friction (implicit case)
                              	CALL wrk_alloc( jpi, jpj, jpk, z3dx, z3dy )
                              	z3dx(:,:,:) = 0._wp   ;   z3dy(:,:,:) = 0._wp  ! after velocity known (now filed at this stage)
                              	DO jk = 1, jpkm1
                                    DO jj = 2, jpjm1
                                       DO ji = 2, jpim1
                              	         ikbu = mbku(ji,jj)          ! deepest ocean u- & v-levels
                                          ikbv = mbkv(ji,jj)
                                          z3dx(ji,jj,jk) = bfrua(ji,jj) * un(ji,jj,ikbu) / e3u_0(ji,jj,ikbu)
                                          z3dy(ji,jj,jk) = bfrva(ji,jj) * vn(ji,jj,ikbv) / e3v_0(ji,jj,ikbv)
                              	      END DO
                              	   END DO
                              	END DO
                              	CALL lbc_lnk( z3dx, 'U', -1. ) ; CALL lbc_lnk( z3dy, 'V', -1. )
                              	CALL iom_put( "utrd_bfri", z3dx )
                              	CALL iom_put( "vtrd_bfri", z3dy )
                              	CALL wrk_dealloc( jpi, jpj, jpk, z3dx, z3dy )
         					      ENDIF
      END SELECT
      !
   END SUBROUTINE trd_dyn_iom

   !!======================================================================
END MODULE trddyn
