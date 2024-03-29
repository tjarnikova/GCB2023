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

MODULE zdfddm
   !!======================================================================
   !!                       ***  MODULE  zdfddm  ***
   !! Ocean physics : double diffusion mixing parameterization
   !!======================================================================
   !! History :  OPA  ! 2000-08  (G. Madec)  double diffusive mixing
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module
   !!            3.3  ! 2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!            3.6  ! 2013-04  (G. Madec, F. Roquet) zrau compute locally using interpolation of alpha & beta
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_zdfddm' :                                     double diffusion
   !!----------------------------------------------------------------------
   !!   zdf_ddm       : compute the Ks for salinity
   !!   zdf_ddm_init  : read namelist and control the parameters
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE zdf_oce         ! ocean vertical physics variables
   USE eosbn2         ! equation of state
   !
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! work arrays
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   zdf_ddm       ! called by step.F90
   PUBLIC   zdf_ddm_init  ! called by opa.F90
   PUBLIC   zdf_ddm_alloc ! called by nemogcm.F90

   LOGICAL , PUBLIC, PARAMETER ::   lk_zdfddm = .TRUE.  !: double diffusive mixing flag

   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::   avs   !: salinity vertical diffusivity coeff. at w-point

   !                       !!* Namelist namzdf_ddm : double diffusive mixing *
   REAL(wp) ::   rn_avts    ! maximum value of avs for salt fingering
   REAL(wp) ::   rn_hsbfr   ! heat/salt buoyancy flux ratio

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
   !! NEMO/OPA 3.7 , NEMO Consortium (2014)
   !! $Id: zdfddm.F90 5120 2015-03-03 16:11:55Z acc $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION zdf_ddm_alloc()
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE zdf_ddm_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( avs(jpi,jpj,jpk) , STAT= zdf_ddm_alloc )
      IF( lk_mpp             )   CALL mpp_sum ( zdf_ddm_alloc )
      IF( zdf_ddm_alloc /= 0 )   CALL ctl_warn('zdf_ddm_alloc: failed to allocate arrays')
   END FUNCTION zdf_ddm_alloc


   SUBROUTINE zdf_ddm( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_ddm  ***
      !!                    
      !! ** Purpose :   Add to the vertical eddy diffusivity coefficient the 
      !!              effect of salt fingering and diffusive convection. 
      !!
      !! ** Method  :   Diapycnal mixing is increased in case of double
      !!      diffusive mixing (i.e. salt fingering and diffusive layering)
      !!      following Merryfield et al. (1999). The rate of double diffusive 
      !!      mixing depend on the buoyancy ratio (R=alpha/beta dk[T]/dk[S]):
      !!         * salt fingering (Schmitt 1981):
      !!      for R > 1 and rn2 > 0 : zavfs = rn_avts / ( 1 + (R/rn_hsbfr)^6 )
      !!      for R > 1 and rn2 > 0 : zavfs = O
      !!      otherwise                : zavft = 0.7 zavs / R
      !!         * diffusive layering (Federov 1988):
      !!      for 0< R < 1 and N^2 > 0 : zavdt = 1.3635e-6 * exp( 4.6 exp(-0.54 (1/R-1) ) )
      !!      otherwise                   : zavdt = 0 
      !!      for .5 < R < 1 and N^2 > 0 : zavds = zavdt (1.885 R -0.85)
      !!      for  0 < R <.5 and N^2 > 0 : zavds = zavdt 0.15 R      
      !!      otherwise                     : zavds = 0 
      !!         * update the eddy diffusivity:
      !!      avt = avt + zavft + zavdt
      !!      avs = avs + zavfs + zavds
      !!      avmu, avmv are required to remain at least above avt and avs.
      !!      
      !! ** Action  :   avt, avs : updated vertical eddy diffusivity coef. for T & S
      !!
      !! References :   Merryfield et al., JPO, 29, 1124-1142, 1999.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step indexocean time step
      !
      INTEGER  ::   ji, jj , jk     ! dummy loop indices
      REAL(wp) ::   zaw, zbw, zrw   ! local scalars
      REAL(wp) ::   zdt, zds
      REAL(wp) ::   zinr, zrr       !   -      -
      REAL(wp) ::   zavft, zavfs    !   -      -
      REAL(wp) ::   zavdt, zavds    !   -      -
      REAL(wp), POINTER, DIMENSION(:,:) ::   zrau, zmsks, zmskf, zmskd1, zmskd2, zmskd3
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('zdf_ddm')
      !
      CALL wrk_alloc( jpi,jpj, zrau, zmsks, zmskf, zmskd1, zmskd2, zmskd3 )
      !
      !                                                ! ===============
      DO jk = 2, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         ! Define the mask 
         ! ---------------
         DO jj = 1, jpj                                ! R=zrau = (alpha / beta) (dk[t] / dk[s])
            DO ji = 1, jpi
               zrw =   ( gdepw_0(ji,jj,jk  ) - gdept_0(ji,jj,jk) )   &
                  &  / ( gdept_0(ji,jj,jk-1) - gdept_0(ji,jj,jk) ) 
               !
               zaw = (  rab_n(ji,jj,jk,jp_tem) * (1. - zrw) + rab_n(ji,jj,jk-1,jp_tem) * zrw  )  &
                   &    * tmask(ji,jj,jk) * tmask(ji,jj,jk-1)
               zbw = (  rab_n(ji,jj,jk,jp_sal) * (1. - zrw) + rab_n(ji,jj,jk-1,jp_sal) * zrw  )  &
                   &    * tmask(ji,jj,jk) * tmask(ji,jj,jk-1)
               !
               zdt = zaw * ( tsn(ji,jj,jk-1,jp_tem) - tsn(ji,jj,jk,jp_tem) )
               zds = zbw * ( tsn(ji,jj,jk-1,jp_sal) - tsn(ji,jj,jk,jp_sal) ) 
               IF( ABS( zds) <= 1.e-20_wp )   zds = 1.e-20_wp
               zrau(ji,jj) = MAX(  1.e-20, zdt / zds  )    ! only retains positive value of zrau
            END DO
         END DO

         DO jj = 1, jpj                                     ! indicators:
            DO ji = 1, jpi
               ! stability indicator: msks=1 if rn2>0; 0 elsewhere
               IF( rn2(ji,jj,jk) + 1.e-12  <= 0. ) THEN   ;   zmsks(ji,jj) = 0._wp
               ELSE                                       ;   zmsks(ji,jj) = 1._wp
               ENDIF
               ! salt fingering indicator: msksf=1 if R>1; 0 elsewhere            
               IF( zrau(ji,jj) <= 1.             ) THEN   ;   zmskf(ji,jj) = 0._wp
               ELSE                                       ;   zmskf(ji,jj) = 1._wp
               ENDIF
               ! diffusive layering indicators: 
               !     ! mskdl1=1 if 0< R <1; 0 elsewhere
               IF( zrau(ji,jj) >= 1.             ) THEN   ;   zmskd1(ji,jj) = 0._wp
               ELSE                                       ;   zmskd1(ji,jj) = 1._wp
               ENDIF
               !     ! mskdl2=1 if 0< R <0.5; 0 elsewhere
               IF( zrau(ji,jj) >= 0.5            ) THEN   ;   zmskd2(ji,jj) = 0._wp
               ELSE                                       ;   zmskd2(ji,jj) = 1._wp
               ENDIF
               !   mskdl3=1 if 0.5< R <1; 0 elsewhere
               IF( zrau(ji,jj) <= 0.5 .OR. zrau(ji,jj) >= 1. ) THEN   ;   zmskd3(ji,jj) = 0._wp
               ELSE                                                   ;   zmskd3(ji,jj) = 1._wp
               ENDIF
            END DO
         END DO
         ! mask zmsk in order to have avt and avs masked
         zmsks(:,:) = zmsks(:,:) * wmask(:,:,jk)


         ! Update avt and avs
         ! ------------------
         ! Constant eddy coefficient: reset to the background value
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               zinr = 1._wp / zrau(ji,jj)
               ! salt fingering
               zrr = zrau(ji,jj) / rn_hsbfr
               zrr = zrr * zrr
               zavfs = rn_avts / ( 1 + zrr*zrr*zrr ) * zmsks(ji,jj) * zmskf(ji,jj)
               zavft = 0.7 * zavfs * zinr
               ! diffusive layering
               zavdt = 1.3635e-6 * EXP(  4.6 * EXP( -0.54*(zinr-1.) )  ) * zmsks(ji,jj) * zmskd1(ji,jj)
               zavds = zavdt * zmsks(ji,jj) * (  ( 1.85 * zrau(ji,jj) - 0.85 ) * zmskd3(ji,jj)   &
                  &                             +  0.15 * zrau(ji,jj)          * zmskd2(ji,jj)  )
               ! add to the eddy viscosity coef. previously computed
               avs (ji,jj,jk) = avt(ji,jj,jk) + zavfs + zavds
               avt (ji,jj,jk) = avt(ji,jj,jk) + zavft + zavdt
               avm (ji,jj,jk) = avm(ji,jj,jk) + MAX( zavft + zavdt, zavfs + zavds )
            END DO
         END DO


         ! Increase avmu, avmv if necessary
         ! --------------------------------
!!gm to be changed following the definition of avm.
         DO jj = 1, jpjm1
            DO ji = 1, jpim1   ! vector opt.
               avmu(ji,jj,jk) = MAX( avmu(ji,jj,jk),    &
                  &                  avt(ji,jj,jk), avt(ji+1,jj,jk),   &
                  &                  avs(ji,jj,jk), avs(ji+1,jj,jk) )  * wumask(ji,jj,jk)
               avmv(ji,jj,jk) = MAX( avmv(ji,jj,jk),    &
                  &                  avt(ji,jj,jk), avt(ji,jj+1,jk),   &
                  &                  avs(ji,jj,jk), avs(ji,jj+1,jk) )  * wvmask(ji,jj,jk)
            END DO
         END DO
         !                                                ! ===============
      END DO                                              !   End of slab
      !                                                   ! ===============
      !
      CALL lbc_lnk( avt , 'W', 1._wp )     ! Lateral boundary conditions   (unchanged sign)
      CALL lbc_lnk( avs , 'W', 1._wp )
      CALL lbc_lnk( avm , 'W', 1._wp )
      CALL lbc_lnk( avmu, 'U', 1._wp ) 
      CALL lbc_lnk( avmv, 'V', 1._wp )

      IF(ln_ctl) THEN
         CALL prt_ctl(tab3d_1=avt , clinfo1=' ddm  - t: ', tab3d_2=avs , clinfo2=' s: ', ovlap=1, kdim=jpk)
         CALL prt_ctl(tab3d_1=avmu, clinfo1=' ddm  - u: ', mask1=umask, &
            &         tab3d_2=avmv, clinfo2=       ' v: ', mask2=vmask, ovlap=1, kdim=jpk)
      ENDIF
      !
      CALL wrk_dealloc( jpi,jpj, zrau, zmsks, zmskf, zmskd1, zmskd2, zmskd3 )
      !
      IF( nn_timing == 1 )  CALL timing_stop('zdf_ddm')
      !
   END SUBROUTINE zdf_ddm
   
   
   SUBROUTINE zdf_ddm_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_ddm_init  ***
      !!
      !! ** Purpose :   Initialization of double diffusion mixing scheme
      !!
      !! ** Method  :   Read the namzdf_ddm namelist and check the parameter values
      !!              called by zdf_ddm at the first timestep (nit000)
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! local integer
      !!
      NAMELIST/namzdf_ddm/ rn_avts, rn_hsbfr
      !!----------------------------------------------------------------------
      !
      REWIND( numnam_ref )              ! Namelist namzdf_ddm in reference namelist : Double diffusion mixing scheme
      READ  ( numnam_ref, namzdf_ddm, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namzdf_ddm in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namzdf_ddm in configuration namelist : Double diffusion mixing scheme
      READ  ( numnam_cfg, namzdf_ddm, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namzdf_ddm in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namzdf_ddm )
      !
      IF(lwp) THEN                    ! Parameter print
         WRITE(numout,*)
         WRITE(numout,*) 'zdf_ddm : double diffusive mixing'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '   Namelist namzdf_ddm : set dd mixing parameter'
         WRITE(numout,*) '      maximum avs for dd mixing      rn_avts   = ', rn_avts
         WRITE(numout,*) '      heat/salt buoyancy flux ratio  rn_hsbfr  = ', rn_hsbfr
      ENDIF
      !
      !                               ! allocate zdfddm arrays
      IF( zdf_ddm_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'zdf_ddm_init : unable to allocate arrays' )
      !                               ! initialization to masked Kz
      avs(:,:,:) = rn_avt0 * wmask(:,:,:) 
      !
   END SUBROUTINE zdf_ddm_init


   !!======================================================================
END MODULE zdfddm
