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

MODULE dynkeg
   !!======================================================================
   !!                       ***  MODULE  dynkeg  ***
   !! Ocean dynamics:  kinetic energy gradient trend
   !!======================================================================
   !! History :  1.0  !  1987-09  (P. Andrich, M.-A. Foujols)  Original code
   !!            7.0  !  1997-05  (G. Madec)  Split dynber into dynkeg and dynhpg
   !!  NEMO      1.0  !  2002-07  (G. Madec)  F90: Free form and module
   !!            3.6  !  2015-05  (N. Ducousso, G. Madec)  add Hollingsworth scheme as an option 
   !!----------------------------------------------------------------------
   
   !!----------------------------------------------------------------------
   !!   dyn_keg      : update the momentum trend with the horizontal tke
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE trd_oce         ! trends: ocean variables
   USE trddyn          ! trend manager: dynamics
   !
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! MPP library
   USE prtctl          ! Print control
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_keg    ! routine called by step module
   
   INTEGER, PARAMETER, PUBLIC  ::   nkeg_C2  = 0   !: 2nd order centered scheme (standard scheme)
   INTEGER, PARAMETER, PUBLIC  ::   nkeg_HW  = 1   !: Hollingsworth et al., QJRMS, 1983
   !
   REAL(wp) ::   r1_48 = 1._wp / 48._wp   !: =1/(4*2*6)
   
   !! * Substitutions
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
   !! NEMO/OPA 3.6 , NEMO Consortium (2015)
   !! $Id: dynkeg.F90 5328 2015-06-01 13:07:18Z pabouttier $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_keg( kt, kscheme )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_keg  ***
      !!
      !! ** Purpose :   Compute the now momentum trend due to the horizontal
      !!      gradient of the horizontal kinetic energy and add it to the 
      !!      general momentum trend.
      !!
      !! ** Method  : * kscheme = nkeg_C2 : 2nd order centered scheme that 
      !!      conserve kinetic energy. Compute the now horizontal kinetic energy 
      !!         zhke = 1/2 [ mi-1( un^2 ) + mj-1( vn^2 ) ]
      !!              * kscheme = nkeg_HW : Hollingsworth correction following
      !!      Arakawa (2001). The now horizontal kinetic energy is given by:
      !!         zhke = 1/6 [ mi-1(  2 * un^2 + ((un(j+1)+un(j-1))/2)^2  )
      !!                    + mj-1(  2 * vn^2 + ((vn(i+1)+vn(i-1))/2)^2  ) ]
      !!      
      !!      Take its horizontal gradient and add it to the general momentum
      !!      trend (ua,va).
      !!         ua = ua - 1/e1u di[ zhke ]
      !!         va = va - 1/e2v dj[ zhke ]
      !!
      !! ** Action : - Update the (ua, va) with the hor. ke gradient trend
      !!             - send this trends to trd_dyn (l_trddyn=T) for post-processing
      !!
      !! ** References : Arakawa, A., International Geophysics 2001.
      !!                 Hollingsworth et al., Quart. J. Roy. Meteor. Soc., 1983.
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt        ! ocean time-step index
      INTEGER, INTENT( in ) ::   kscheme   ! =0/1   type of KEG scheme 
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zu, zv       ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zhke
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ztrdu, ztrdv 
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('dyn_keg')
      !
      CALL wrk_alloc( jpi,jpj,jpk,   zhke )
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_keg : kinetic energy gradient trend, scheme number=', kscheme
         IF(lwp) WRITE(numout,*) '~~~~~~~'
      ENDIF

      IF( l_trddyn ) THEN           ! Save ua and va trends
         CALL wrk_alloc( jpi,jpj,jpk,   ztrdu, ztrdv )
         ztrdu(:,:,:) = ua(:,:,:) 
         ztrdv(:,:,:) = va(:,:,:) 
      ENDIF
      
      zhke(:,:,jpk) = 0._wp
      
      SELECT CASE ( kscheme )             !== Horizontal kinetic energy at T-point  ==!
      !
      CASE ( nkeg_C2 )                          !--  Standard scheme  --!
         DO jk = 1, jpkm1
            DO jj = 2, jpj
               DO ji = 2, jpi   ! vector opt.
                  zu =    un(ji-1,jj  ,jk) * un(ji-1,jj  ,jk)   &
                     &  + un(ji  ,jj  ,jk) * un(ji  ,jj  ,jk)
                  zv =    vn(ji  ,jj-1,jk) * vn(ji  ,jj-1,jk)   &
                     &  + vn(ji  ,jj  ,jk) * vn(ji  ,jj  ,jk)
                  zhke(ji,jj,jk) = 0.25_wp * ( zv + zu )
               END DO  
            END DO
         END DO
         !
      CASE ( nkeg_HW )                          !--  Hollingsworth scheme  --!
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1       
               DO ji = 2, jpim1   ! vector opt.
                  zu = 8._wp * ( un(ji-1,jj  ,jk) * un(ji-1,jj  ,jk)    &
                     &         + un(ji  ,jj  ,jk) * un(ji  ,jj  ,jk) )  &
                     &   +     ( un(ji-1,jj-1,jk) + un(ji-1,jj+1,jk) ) * ( un(ji-1,jj-1,jk) + un(ji-1,jj+1,jk) )   &
                     &   +     ( un(ji  ,jj-1,jk) + un(ji  ,jj+1,jk) ) * ( un(ji  ,jj-1,jk) + un(ji  ,jj+1,jk) )
                     !
                  zv = 8._wp * ( vn(ji  ,jj-1,jk) * vn(ji  ,jj-1,jk)    &
                     &         + vn(ji  ,jj  ,jk) * vn(ji  ,jj  ,jk) )  &
                     &  +      ( vn(ji-1,jj-1,jk) + vn(ji+1,jj-1,jk) ) * ( vn(ji-1,jj-1,jk) + vn(ji+1,jj-1,jk) )   &
                     &  +      ( vn(ji-1,jj  ,jk) + vn(ji+1,jj  ,jk) ) * ( vn(ji-1,jj  ,jk) + vn(ji+1,jj  ,jk) )
                  zhke(ji,jj,jk) = r1_48 * ( zv + zu )
               END DO  
            END DO
         END DO
         CALL lbc_lnk( zhke, 'T', 1. )
         !
      END SELECT
      !
      DO jk = 1, jpkm1                    !==  grad( KE ) added to the general momentum trends  ==!
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               ua(ji,jj,jk) = ua(ji,jj,jk) - ( zhke(ji+1,jj  ,jk) - zhke(ji,jj,jk) ) / e1u(ji,jj)
               va(ji,jj,jk) = va(ji,jj,jk) - ( zhke(ji  ,jj+1,jk) - zhke(ji,jj,jk) ) / e2v(ji,jj)
            END DO 
         END DO
      END DO
      !
      IF( l_trddyn ) THEN                 ! save the Kinetic Energy trends for diagnostic
         ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
         ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_keg, kt )
         CALL wrk_dealloc( jpi,jpj,jpk,   ztrdu, ztrdv )
      ENDIF
      !
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' keg  - Ua: ', mask1=umask,   &
         &                       tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      CALL wrk_dealloc( jpi,jpj,jpk,   zhke )
      !
      IF( nn_timing == 1 )   CALL timing_stop('dyn_keg')
      !
   END SUBROUTINE dyn_keg

   !!======================================================================
END MODULE dynkeg
