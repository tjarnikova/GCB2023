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

MODULE dynzdf_exp
   !!==============================================================================
   !!                     ***  MODULE  dynzdf_exp  ***
   !! Ocean dynamics:  vertical component(s) of the momentum mixing trend
   !!==============================================================================
   !! History :  OPA  !  1990-10  (B. Blanke)  Original code
   !!            8.0  !  1997-05  (G. Madec)  vertical component of isopycnal
   !!   NEMO     0.5  !  2002-08  (G. Madec)  F90: Free form and module
   !!            3.3  !  2010-04  (M. Leclair, G. Madec)  Forcing averaged over 2 time steps
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_zdf_exp  : update the momentum trend with the vertical diffu-
   !!                  sion using an explicit time-stepping scheme.
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE zdf_oce         ! ocean vertical physics
   USE sbc_oce         ! surface boundary condition: ocean
   USE lib_mpp         ! MPP library
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing


   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_zdf_exp   ! called by step.F90
   
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
   !! $Id: dynzdf_exp.F90 3625 2012-11-21 13:19:18Z acc $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_zdf_exp( kt, p2dt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_zdf_exp  ***
      !!                   
      !! ** Purpose :   Compute the trend due to the vert. momentum diffusion
      !!
      !! ** Method  :   Explicit forward time stepping with a time splitting
      !!      technique. The vertical diffusion of momentum is given by:
      !!         diffu = dz( avmu dz(u) ) = 1/e3u dk+1( avmu/e3uw dk(ub) )
      !!      Surface boundary conditions: wind stress input (averaged over kt-1/2 & kt+1/2)
      !!      Bottom boundary conditions : bottom stress (cf zdfbfr.F90)
      !!      Add this trend to the general trend ua :
      !!         ua = ua + dz( avmu dz(u) )
      !!
      !! ** Action : - Update (ua,va) with the vertical diffusive trend
      !!---------------------------------------------------------------------
      INTEGER , INTENT(in) ::   kt     ! ocean time-step index
      REAL(wp), INTENT(in) ::   p2dt   ! time-step 
      !
      INTEGER  ::   ji, jj, jk, jl   ! dummy loop indices
      REAL(wp) ::   zlavmr, zua, zva   ! local scalars
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zwx, zwy, zwz, zww
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_zdf_exp')
      !
      CALL wrk_alloc( jpi,jpj,jpk, zwx, zwy, zwz, zww ) 
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_zdf_exp : vertical momentum diffusion - explicit operator'
         WRITE(numout,*) '~~~~~~~~~~~ '
      ENDIF

      zlavmr = 1. / REAL( nn_zdfexp )


      DO jj = 2, jpjm1                 ! Surface boundary condition
         DO ji = 2, jpim1
            zwy(ji,jj,1) = ( utau_b(ji,jj) + utau(ji,jj) ) * r1_rau0
            zww(ji,jj,1) = ( vtau_b(ji,jj) + vtau(ji,jj) ) * r1_rau0
         END DO  
      END DO  
      DO jk = 1, jpk                   ! Initialization of x, z and contingently trends array
         DO jj = 2, jpjm1 
            DO ji = 2, jpim1
               zwx(ji,jj,jk) = ub(ji,jj,jk)
               zwz(ji,jj,jk) = vb(ji,jj,jk)
            END DO  
         END DO  
      END DO  
      !
      DO jl = 1, nn_zdfexp             ! Time splitting loop
         !
         DO jk = 2, jpk                      ! First vertical derivative
            DO jj = 2, jpjm1 
               DO ji = 2, jpim1
                  zwy(ji,jj,jk) = avmu(ji,jj,jk) * ( zwx(ji,jj,jk-1) - zwx(ji,jj,jk) ) / e3uw_0(ji,jj,jk) 
                  zww(ji,jj,jk) = avmv(ji,jj,jk) * ( zwz(ji,jj,jk-1) - zwz(ji,jj,jk) ) / e3vw_0(ji,jj,jk)
               END DO  
            END DO  
         END DO  
         DO jk = 1, jpkm1                    ! Second vertical derivative and trend estimation at kt+l*rdt/nn_zdfexp
            DO jj = 2, jpjm1 
               DO ji = 2, jpim1
                  zua = zlavmr * ( zwy(ji,jj,jk) - zwy(ji,jj,jk+1) ) / e3u_0(ji,jj,jk)
                  zva = zlavmr * ( zww(ji,jj,jk) - zww(ji,jj,jk+1) ) / e3v_0(ji,jj,jk)
                  ua(ji,jj,jk) = ua(ji,jj,jk) + zua
                  va(ji,jj,jk) = va(ji,jj,jk) + zva
                  !
                  zwx(ji,jj,jk) = zwx(ji,jj,jk) + p2dt * zua * umask(ji,jj,jk)
                  zwz(ji,jj,jk) = zwz(ji,jj,jk) + p2dt * zva * vmask(ji,jj,jk)
               END DO  
            END DO  
         END DO  
         !
      END DO                           ! End of time splitting
      !
      CALL wrk_dealloc( jpi,jpj,jpk, zwx, zwy, zwz, zww ) 
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_zdf_exp')
      !
   END SUBROUTINE dyn_zdf_exp

   !!==============================================================================
END MODULE dynzdf_exp
