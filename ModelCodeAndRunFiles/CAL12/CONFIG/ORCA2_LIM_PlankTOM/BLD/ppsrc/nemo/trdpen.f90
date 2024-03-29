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

MODULE trdpen
   !!======================================================================
   !!                       ***  MODULE  trdpen  ***
   !! Ocean diagnostics:  Potential Energy trends
   !!=====================================================================
   !! History :  3.5  !  2012-02  (G. Madec) original code 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   trd_pen       : compute and output Potential Energy trends from T & S trends
   !!   trd_pen_init  : initialisation of PE trends
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean domain 
   USE sbc_oce        ! surface boundary condition: ocean
   USE zdf_oce        ! ocean vertical physics
   USE trd_oce        ! trends: ocean variables
   USE eosbn2         ! equation of state and related derivatives
   USE ldftra_oce     ! ocean active tracers lateral physics
   USE zdfddm         ! vertical physics: double diffusion
   USE phycst         ! physical constants
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE wrk_nemo       ! Memory allocation

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trd_pen        ! called by all trdtra module
   PUBLIC   trd_pen_init   ! called by all nemogcm module

   INTEGER ::   nkstp   ! current time step 

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   rab_pe   ! partial derivatives of PE anomaly with respect to T and S

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
   !!                    *** zdfddm_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
   !!      with a constant or 1D or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
!   'key_zdfddm' :                      avs: 3D array defined in zdfddm module
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: zdfddm_substitute.h90 4152 2013-11-05 11:59:53Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
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
   !! $Id: trdpen.F90 5215 2015-04-15 16:11:56Z nicolasmartin $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION trd_pen_alloc()
      !!---------------------------------------------------------------------
      !!                  ***  FUNCTION trd_tra_alloc  ***
      !!---------------------------------------------------------------------
      ALLOCATE( rab_pe(jpi,jpj,jpk,jpts) , STAT= trd_pen_alloc )
      !
      IF( lk_mpp             )   CALL mpp_sum ( trd_pen_alloc )
      IF( trd_pen_alloc /= 0 )   CALL ctl_warn( 'trd_pen_alloc: failed to allocate arrays' )
   END FUNCTION trd_pen_alloc


   SUBROUTINE trd_pen( ptrdx, ptrdy, ktrd, kt, pdt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_tra_mng  ***
      !! 
      !! ** Purpose :   Dispatch all trends computation, e.g. 3D output, integral
      !!                constraints, barotropic vorticity, kinetic enrgy, 
      !!                potential energy, and/or mixed layer budget.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in) ::   ptrdx, ptrdy   ! Temperature & Salinity trends
      INTEGER                   , INTENT(in) ::   ktrd           ! tracer trend index
      INTEGER                   , INTENT(in) ::   kt             ! time step index
      REAL(wp)                  , INTENT(in) ::   pdt            ! time step [s]
      !
      INTEGER ::   jk                                            ! dummy loop indices
      REAL(wp), POINTER, DIMENSION(:,:)      ::   z2d            ! 2D workspace 
      REAL(wp), POINTER, DIMENSION(:,:,:)    ::   zpe            ! 3D workspace 
      !!----------------------------------------------------------------------
      !
      CALL wrk_alloc( jpi, jpj, jpk, zpe )
      zpe(:,:,:) = 0._wp
      !
      IF ( kt /= nkstp ) THEN   ! full eos: set partial derivatives at the 1st call of kt time step
         nkstp = kt
         CALL eos_pen( tsn, rab_PE, zpe )
         CALL iom_put( "alphaPE", rab_pe(:,:,:,jp_tem) )
         CALL iom_put( "betaPE" , rab_pe(:,:,:,jp_sal) )
         CALL iom_put( "PEanom" , zpe )
      ENDIF
      !
      zpe(:,:,jpk) = 0._wp
      DO jk = 1, jpkm1
         zpe(:,:,jk) = ( - ( rab_n(:,:,jk,jp_tem) + rab_pe(:,:,jk,jp_tem) ) * ptrdx(:,:,jk)   &
            &            + ( rab_n(:,:,jk,jp_sal) + rab_pe(:,:,jk,jp_sal) ) * ptrdy(:,:,jk)  )
      END DO

      SELECT CASE ( ktrd )
      CASE ( jptra_xad  )   ;   CALL iom_put( "petrd_xad", zpe )   ! zonal    advection
      CASE ( jptra_yad  )   ;   CALL iom_put( "petrd_yad", zpe )   ! merid.   advection
      CASE ( jptra_zad  )   ;   CALL iom_put( "petrd_zad", zpe )   ! vertical advection
                                IF( .NOT.lk_vvl ) THEN                   ! cst volume : adv flux through z=0 surface
                                   CALL wrk_alloc( jpi, jpj, z2d )
                                   z2d(:,:) = wn(:,:,1) * ( &
                                   	  &   - ( rab_n(:,:,1,jp_tem) + rab_pe(:,:,1,jp_tem) ) * tsn(:,:,1,jp_tem)    &
                                   	  &   + ( rab_n(:,:,1,jp_sal) + rab_pe(:,:,1,jp_sal) ) * tsn(:,:,1,jp_sal)    &
                                      &   					) / e3t_0(:,:,1)
                                   CALL iom_put( "petrd_sad" , z2d )
                                   CALL wrk_dealloc( jpi, jpj, z2d )
                                ENDIF
      CASE ( jptra_ldf  )   ;   CALL iom_put( "petrd_ldf" , zpe )   ! lateral  diffusion
      CASE ( jptra_zdf  )   ;   CALL iom_put( "petrd_zdf" , zpe )   ! lateral  diffusion (K_z)
      CASE ( jptra_zdfp )   ;   CALL iom_put( "petrd_zdfp", zpe )   ! vertical diffusion (K_z)
      CASE ( jptra_dmp  )   ;   CALL iom_put( "petrd_dmp" , zpe )   ! internal 3D restoring (tradmp)
      CASE ( jptra_bbl  )   ;   CALL iom_put( "petrd_bbl" , zpe )   ! bottom boundary layer
      CASE ( jptra_npc  )   ;   CALL iom_put( "petrd_npc" , zpe )   ! non penetr convect adjustment
      CASE ( jptra_nsr  )   ;   CALL iom_put( "petrd_nsr" , zpe )   ! surface forcing + runoff (ln_rnf=T)
      CASE ( jptra_qsr  )   ;   CALL iom_put( "petrd_qsr" , zpe )   ! air-sea : penetrative sol radiat
      CASE ( jptra_bbc  )   ;   CALL iom_put( "petrd_bbc" , zpe )   ! bottom bound cond (geoth flux)
      CASE ( jptra_atf  )   ;   CALL iom_put( "petrd_atf" , zpe )   ! asselin time filter (last trend)
                                !IF( .NOT.lk_vvl ) THEN                   ! cst volume : ssh term (otherwise include in e3t variation)
                                !   CALL wrk_alloc( jpi, jpj, z2d )
                                !   z2d(:,:) = ( ssha(:,:) - sshb(:,:) )                 &
                                !      &     * (   dPE_dt(:,:,1) * tsn(:,:,1,jp_tem)    &
                                !      &         + dPE_ds(:,:,1) * tsn(:,:,1,jp_sal)  ) / ( e3t_0(:,:,1) * pdt )
                                !   CALL iom_put( "petrd_sad" , z2d )
                                !   CALL wrk_dealloc( jpi, jpj, z2d )
                                !ENDIF
         !
      END SELECT
      !
      CALL wrk_dealloc( jpi, jpj, jpk, zpe )
      !
   END SUBROUTINE trd_pen


   SUBROUTINE trd_pen_init
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_pen_init  ***
      !! 
      !! ** Purpose :   initialisation of 3D Kinetic Energy trend diagnostic
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trd_pen_init : 3D Potential ENergy trends'
         WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF
      !                           ! allocate box volume arrays
      IF ( trd_pen_alloc() /= 0 )   CALL ctl_stop('trd_pen_alloc: failed to allocate arrays')
      !
      rab_pe(:,:,:,:) = 0._wp
      !
      IF ( lk_vvl               )   CALL ctl_stop('trd_pen_init : PE trends not coded for variable volume')
      !
      nkstp     = nit000 - 1
      !
   END SUBROUTINE trd_pen_init

   !!======================================================================
END MODULE trdpen
