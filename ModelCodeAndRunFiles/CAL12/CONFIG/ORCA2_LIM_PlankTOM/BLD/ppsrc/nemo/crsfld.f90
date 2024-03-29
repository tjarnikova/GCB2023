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

MODULE crsfld
   !!======================================================================
   !!                     ***  MODULE  crsdfld  ***
   !!  Ocean coarsening :  coarse ocean fields
   !!=====================================================================
   !!   2012-07  (J. Simeon, C. Calone, G. Madec, C. Ethe)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   crs_fld       : create the standard output files for coarse grid and prep
   !!                       other variables needed to be passed to TOP
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE zdf_oce         ! vertical  physics: ocean fields
   USE zdfddm          ! vertical  physics: double diffusion
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  ! I/O manager
   USE timing          ! preformance summary
   USE wrk_nemo        ! working array
   USE crs
   USE crsdom
   USE crslbclnk
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   crs_fld                 ! routines called by step.F90


   !! * Substitutions
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
   !! $Id: crsfld.F90 5217 2015-04-16 10:24:12Z nicolasmartin $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE crs_fld( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE crs_fld  ***
      !!                   
      !! ** Purpose :   Basic output of coarsened dynamics and tracer fields 
      !!      NETCDF format is used by default 
      !!      1. Accumulate in time the dimensionally-weighted fields
      !!      2. At time of output, rescale [1] by dimension and time
      !!         to yield the spatial and temporal average.
      !!  See. diawri_dimg.h90, sbcmod.F90
      !!
      !! ** Method  :  
      !!----------------------------------------------------------------------
      !!
      
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !!
      INTEGER               ::   ji, jj, jk              ! dummy loop indices
      !!
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zfse3t, zfse3u, zfse3v, zfse3w ! 3D workspace for e3
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zt, zs 
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zt_crs, zs_crs !
      REAL(wp)       :: z2dcrsu, z2dcrsv
      !!
       !!----------------------------------------------------------------------
      ! 

      IF( nn_timing == 1 )   CALL timing_start('crs_fld')

      !  Initialize arrays
      CALL wrk_alloc( jpi, jpj, jpk, zfse3t, zfse3w )
      CALL wrk_alloc( jpi, jpj, jpk, zfse3u, zfse3v )
      CALL wrk_alloc( jpi, jpj, jpk, zt, zs       )
      !
      CALL wrk_alloc( jpi_crs, jpj_crs, jpk, zt_crs, zs_crs )

      ! Depth work arrrays
      zfse3t(:,:,:) = e3t_0(:,:,:)
      zfse3u(:,:,:) = e3u_0(:,:,:)
      zfse3v(:,:,:) = e3v_0(:,:,:)
      zfse3w(:,:,:) = e3w_0(:,:,:)

      IF( kt == nit000  ) THEN
         tsn_crs  (:,:,:,:) = 0._wp    ! temp/sal  array, now 
         un_crs   (:,:,:  ) = 0._wp    ! u-velocity
         vn_crs   (:,:,:  ) = 0._wp    ! v-velocity
         wn_crs   (:,:,:  ) = 0._wp    ! w
         avt_crs  (:,:,:  ) = 0._wp    ! avt
         hdivn_crs(:,:,:  ) = 0._wp    ! hdiv
         rke_crs  (:,:,:  ) = 0._wp    ! rke
         sshn_crs (:,:    ) = 0._wp    ! ssh
         utau_crs (:,:    ) = 0._wp    ! taux
         vtau_crs (:,:    ) = 0._wp    ! tauy
         wndm_crs (:,:    ) = 0._wp    ! wind speed
         qsr_crs  (:,:    ) = 0._wp    ! qsr
         emp_crs  (:,:    ) = 0._wp    ! emp
         emp_b_crs(:,:    ) = 0._wp    ! emp
         rnf_crs  (:,:    ) = 0._wp    ! runoff
         fr_i_crs (:,:    ) = 0._wp    ! ice cover
      ENDIF

      CALL iom_swap( "nemo_crs" )    ! swap on the coarse grid

      ! 2. Coarsen fields at each time step
      ! --------------------------------------------------------

      !  Temperature
      zt(:,:,:) = tsn(:,:,:,jp_tem)  ;      zt_crs(:,:,:) = 0._wp
      CALL crs_dom_ope( zt, 'VOL', 'T', tmask, zt_crs, p_e12=e1e2t, p_e3=zfse3t, psgn=1.0 )
      tsn_crs(:,:,:,jp_tem) = zt_crs(:,:,:)

      CALL iom_put( "toce", tsn_crs(:,:,:,jp_tem) )    ! temp
      CALL iom_put( "sst" , tsn_crs(:,:,1,jp_tem) )    ! sst

      
      !  Salinity
      zs(:,:,:) = tsn(:,:,:,jp_sal)  ;      zs_crs(:,:,:) = 0._wp
      CALL crs_dom_ope( zs, 'VOL', 'T', tmask, zs_crs, p_e12=e1e2t, p_e3=zfse3t, psgn=1.0 )
      tsn_crs(:,:,:,jp_sal) = zt_crs(:,:,:)

      CALL iom_put( "soce" , tsn_crs(:,:,:,jp_sal) )    ! sal
      CALL iom_put( "sss"  , tsn_crs(:,:,1,jp_sal) )    ! sss

      !  U-velocity
      CALL crs_dom_ope( un, 'SUM', 'U', umask, un_crs, p_e12=e2u, p_e3=zfse3u, p_surf_crs=e2e3u_msk, psgn=-1.0 )
      !
      zt(:,:,:) = 0._wp     ;    zs(:,:,:) = 0._wp  ;   zt_crs(:,:,:) = 0._wp   ;    zs_crs(:,:,:) = 0._wp
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   
               zt(ji,jj,jk)  = un(ji,jj,jk) * 0.5 * ( tsn(ji,jj,jk,jp_tem) + tsn(ji+1,jj,jk,jp_tem) ) 
               zs(ji,jj,jk)  = un(ji,jj,jk) * 0.5 * ( tsn(ji,jj,jk,jp_sal) + tsn(ji+1,jj,jk,jp_sal) ) 
            END DO
         END DO
      END DO
      CALL crs_dom_ope( zt, 'SUM', 'U', umask, zt_crs, p_e12=e2u, p_e3=zfse3u, p_surf_crs=e2e3u_msk, psgn=-1.0 )
      CALL crs_dom_ope( zs, 'SUM', 'U', umask, zs_crs, p_e12=e2u, p_e3=zfse3u, p_surf_crs=e2e3u_msk, psgn=-1.0 )

      CALL iom_put( "uoce"  , un_crs )   ! i-current 
      CALL iom_put( "uocet" , zt_crs )   ! uT
      CALL iom_put( "uoces" , zs_crs )   ! uS

      !  V-velocity
      CALL crs_dom_ope( vn, 'SUM', 'V', vmask, vn_crs, p_e12=e1v, p_e3=zfse3v, p_surf_crs=e1e3v_msk, psgn=-1.0 )
      !                                                                                 
      zt(:,:,:) = 0._wp     ;    zs(:,:,:) = 0._wp  ;   zt_crs(:,:,:) = 0._wp   ;    zs_crs(:,:,:) = 0._wp
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   
               zt(ji,jj,jk)  = vn(ji,jj,jk) * 0.5 * ( tsn(ji,jj,jk,jp_tem) + tsn(ji,jj+1,jk,jp_tem) ) 
               zs(ji,jj,jk)  = vn(ji,jj,jk) * 0.5 * ( tsn(ji,jj,jk,jp_sal) + tsn(ji,jj+1,jk,jp_sal) ) 
            END DO
         END DO
      END DO
      CALL crs_dom_ope( zt, 'SUM', 'V', vmask, zt_crs, p_e12=e1v, p_e3=zfse3v, p_surf_crs=e1e3v_msk, psgn=-1.0 )
      CALL crs_dom_ope( zs, 'SUM', 'V', vmask, zs_crs, p_e12=e1v, p_e3=zfse3v, p_surf_crs=e1e3v_msk, psgn=-1.0 )
 
      CALL iom_put( "voce"  , vn_crs )   ! i-current 
      CALL iom_put( "vocet" , zt_crs )   ! vT
      CALL iom_put( "voces" , zs_crs )   ! vS

     
      !  Kinetic energy
      CALL crs_dom_ope( rke, 'VOL', 'T', tmask, rke_crs, p_e12=e1e2t, p_e3=zfse3t, psgn=1.0 )
      CALL iom_put( "eken", rke_crs )

      !  Horizontal divergence ( following OPA_SRC/DYN/divcur.F90 ) 
      DO jk = 1, jpkm1
         DO ji = 2, jpi_crsm1
            DO jj = 2, jpj_crsm1
               IF( tmask_crs(ji,jj,jk ) > 0 ) THEN
                   z2dcrsu =  ( un_crs(ji  ,jj  ,jk) * crs_surfu_wgt(ji  ,jj  ,jk) ) &
                      &     - ( un_crs(ji-1,jj  ,jk) * crs_surfu_wgt(ji-1,jj  ,jk) )
                   z2dcrsv =  ( vn_crs(ji  ,jj  ,jk) * crs_surfv_wgt(ji  ,jj  ,jk) ) &
                      &     - ( vn_crs(ji  ,jj-1,jk) * crs_surfv_wgt(ji  ,jj-1,jk) )
                   !
                   hdivn_crs(ji,jj,jk) = ( z2dcrsu + z2dcrsv ) / crs_volt_wgt(ji,jj,jk) 
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      CALL crs_lbc_lnk( hdivn_crs, 'T', 1.0 )
      !
      CALL iom_put( "hdiv", hdivn_crs )  


      !  W-velocity
      IF( ln_crs_wn ) THEN
         CALL crs_dom_ope( wn, 'SUM', 'W', tmask, wn_crs, p_e12=e1e2t, p_surf_crs=e1e2w_msk, psgn=1.0 )
       !  CALL crs_dom_ope( wn, 'VOL', 'W', tmask, wn_crs, p_e12=e1e2t, p_e3=zfse3w )
      ELSE
        wn_crs(:,:,jpk) = 0._wp
        DO jk = jpkm1, 1, -1
           wn_crs(:,:,jk) = wn_crs(:,:,jk+1) - e3t_crs(:,:,jk) * hdivn_crs(:,:,jk)
        ENDDO
      ENDIF
      CALL iom_put( "woce", wn_crs  )   ! vertical velocity
      !  free memory

      !  avt, avs
      SELECT CASE ( nn_crs_kz )
         CASE ( 0 )
            CALL crs_dom_ope( avt, 'VOL', 'W', tmask, avt_crs, p_e12=e1e2t, p_e3=zfse3w, psgn=1.0 )
         CASE ( 1 )
            CALL crs_dom_ope( avt, 'MAX', 'W', tmask, avt_crs, p_e12=e1e2t, p_e3=zfse3w, psgn=1.0 )
         CASE ( 2 )
            CALL crs_dom_ope( avt, 'MIN', 'W', tmask, avt_crs, p_e12=e1e2t, p_e3=zfse3w, psgn=1.0 )
      END SELECT
      !
      CALL iom_put( "avt", avt_crs )   !  Kz
      
      !  sbc fields  
      CALL crs_dom_ope( sshn , 'VOL', 'T', tmask, sshn_crs , p_e12=e1e2t, p_e3=zfse3t         , psgn=1.0 )  
      CALL crs_dom_ope( utau , 'SUM', 'U', umask, utau_crs , p_e12=e2u  , p_surf_crs=e2u_crs  , psgn=1.0 )
      CALL crs_dom_ope( vtau , 'SUM', 'V', vmask, vtau_crs , p_e12=e1v  , p_surf_crs=e1v_crs  , psgn=1.0 )
      CALL crs_dom_ope( wndm , 'SUM', 'T', tmask, wndm_crs , p_e12=e1e2t, p_surf_crs=e1e2t_crs, psgn=1.0 )
      CALL crs_dom_ope( rnf  , 'MAX', 'T', tmask, rnf_crs                                     , psgn=1.0 )
      CALL crs_dom_ope( qsr  , 'SUM', 'T', tmask, qsr_crs  , p_e12=e1e2t, p_surf_crs=e1e2t_crs, psgn=1.0 )
      CALL crs_dom_ope( emp_b, 'SUM', 'T', tmask, emp_b_crs, p_e12=e1e2t, p_surf_crs=e1e2t_crs, psgn=1.0 )
      CALL crs_dom_ope( emp  , 'SUM', 'T', tmask, emp_crs  , p_e12=e1e2t, p_surf_crs=e1e2t_crs, psgn=1.0 )
      CALL crs_dom_ope( sfx  , 'SUM', 'T', tmask, sfx_crs  , p_e12=e1e2t, p_surf_crs=e1e2t_crs, psgn=1.0 )
      CALL crs_dom_ope( fr_i , 'SUM', 'T', tmask, fr_i_crs , p_e12=e1e2t, p_surf_crs=e1e2t_crs, psgn=1.0 )

      CALL iom_put( "ssh"      , sshn_crs )   ! ssh output 
      CALL iom_put( "utau"     , utau_crs )   ! i-tau output 
      CALL iom_put( "vtau"     , vtau_crs )   ! j-tau output 
      CALL iom_put( "wspd"     , wndm_crs )   ! wind speed output 
      CALL iom_put( "runoffs"  , rnf_crs  )   ! runoff output 
      CALL iom_put( "qsr"      , qsr_crs  )   ! qsr output 
      CALL iom_put( "empmr"    , emp_crs  )   ! water flux output 
      CALL iom_put( "saltflx"  , sfx_crs  )   ! salt flux output 
      CALL iom_put( "ice_cover", fr_i_crs )   ! ice cover output 

      !  free memory
      CALL wrk_dealloc( jpi, jpj, jpk, zfse3t, zfse3w )
      CALL wrk_dealloc( jpi, jpj, jpk, zfse3u, zfse3v )
      CALL wrk_dealloc( jpi, jpj, jpk, zt, zs       )
      CALL wrk_dealloc( jpi_crs, jpj_crs, jpk, zt_crs, zs_crs )
      !
      CALL iom_swap( "nemo" )     ! return back on high-resolution grid
      !
      IF( nn_timing == 1 )   CALL timing_stop('crs_fld')
      !
   END SUBROUTINE crs_fld

   !!======================================================================
END MODULE crsfld
