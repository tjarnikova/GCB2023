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

MODULE crsini   
   !!======================================================================
   !!                         ***  MODULE crsini  ***
   !!            Manage the grid coarsening module initialization
   !!======================================================================
   !!  History     2012-05   (J. Simeon, G. Madec, C. Ethe, C. Calone) Original code
   !!----------------------------------------------------------------------

   USE timing                   ! Timing
   USE par_oce                  ! For parameter jpi,jpj,jphgr_msh
   USE dom_oce                  ! For parameters in par_oce (jperio, lk_vvl)
   USE crs                  ! Coarse grid domain
   USE phycst, ONLY: omega, rad ! physical constants
   USE wrk_nemo 
   USE in_out_manager
   USE par_kind, ONLY: wp
   USE iom
   USE crsdom
   USE crsdomwri
   USE crslbclnk
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   PUBLIC  crs_init

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

   !! $Id: crsini.F90 5215 2015-04-15 16:11:56Z nicolasmartin $
CONTAINS
   
   SUBROUTINE crs_init 
      !!-------------------------------------------------------------------
      !!                     *** SUBROUTINE crs_init
      !!  ** Purpose : Initialization of the grid coarsening module  
      !!               1. Read namelist
      !!               X2. MOVE TO crs_dom.F90 Set the domain definitions for coarse grid
      !!                 a. Define the coarse grid starting/ending indices on parent grid
      !!                    Here is where the T-pivot or F-pivot grids are discerned
      !!                 b. TODO.  Include option for north-centric or equator-centric binning.
      !!                 (centered only for odd reduction factors; even reduction bins bias equator to the south)
      !!               3. Mask and mesh creation. => calls to crsfun
      !!                  a. Use crsfun_mask to generate tmask,umask, vmask, fmask.
      !!                  b. Use crsfun_coordinates to get coordinates
      !!                  c. Use crsfun_UV to get horizontal scale factors
      !!                  d. Use crsfun_TW to get initial vertical scale factors   
      !!               4. Volumes and weights jes.... TODO. Updates for vvl? Where to do this? crsstp.F90?
      !!                  a. Calculate initial coarse grid box volumes: pvol_T, pvol_W
      !!                  b. Calculate initial coarse grid surface-averaging weights
      !!                  c. Calculate initial coarse grid volume-averaging weights
      !!                  
      !!               X5. MOVE TO crs_dom_wri.F90 Using iom_rstput output the initial meshmask.
      !!               ?. Another set of "masks" to generate
      !!                  are the u- and v- surface areas for U- and V- area-weighted means.
      !!                  Need to put this somewhere in section 3?
      !! jes. What do to about the vvl?  GM.  could separate the weighting (denominator), so
      !! output C*dA or C*dV as summation not mran, then do mean (division) at moment of output.
      !! As is, crsfun takes into account vvl.   
      !!      Talked about pre-setting the surface array to avoid IF/ENDIFS and division.
      !!      But have then to make that preset array here and elsewhere.
      !!      that is called every timestep...
      !! 
      !!               - Read in pertinent data ?
      !!-------------------------------------------------------------------
      !! Local variables
      INTEGER  :: ji,jj,jk      ! dummy indices
      INTEGER  :: ierr                                ! allocation error status
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      REAL(wp), DIMENSION(:,:,:), POINTER :: zfse3t, zfse3u, zfse3v, zfse3w

      NAMELIST/namcrs/ nn_factx, nn_facty, nn_binref, nn_msh_crs, nn_crs_kz, ln_crs_wn
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('crs_init')
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'crs_init : Initializing the grid coarsening module '
      ENDIF

     !---------------------------------------------------------
     ! 1. Read Namelist file
     !---------------------------------------------------------
     !

      REWIND( numnam_ref )              ! Namelist namrun in reference namelist : Parameters of the run
      READ  ( numnam_ref, namcrs, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namcrs in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namrun in configuration namelist : Parameters of the run
      READ  ( numnam_cfg, namcrs, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namcrs in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namcrs )

     IF(lwp) THEN
        WRITE(numout,*)
        WRITE(numout,*) 'crs_init: Namelist namcrs '
        WRITE(numout,*) '   coarsening factor in i-direction      nn_factx   = ', nn_factx
        WRITE(numout,*) '   coarsening factor in j-direction      nn_facty   = ', nn_facty
        WRITE(numout,*) '   bin centering preference              nn_binref  = ', nn_binref
        WRITE(numout,*) '   create (=1) a mesh file or not (=0)   nn_msh_crs = ', nn_msh_crs
        WRITE(numout,*) '   type of Kz coarsening (0,1,2)         nn_crs_kz  = ', nn_crs_kz
        WRITE(numout,*) '   wn coarsened or computed using hdivn  ln_crs_wn  = ', ln_crs_wn
     ENDIF
              
     rfactx_r = 1. / nn_factx
     rfacty_r = 1. / nn_facty

     !---------------------------------------------------------
     ! 2. Define Global Dimensions of the coarsened grid
     !---------------------------------------------------------
     CALL crs_dom_def      

     !---------------------------------------------------------
     ! 3. Mask and Mesh
     !---------------------------------------------------------

     !     Set up the masks and meshes     

     !  3.a. Get the masks   
 
     CALL crs_dom_msk


     !  3.b. Get the coordinates
     !      Odd-numbered reduction factor, center coordinate on T-cell
     !      Even-numbered reduction factor, center coordinate on U-,V- faces or f-corner.
     !      
     IF ( nresty /= 0 .AND. nrestx /= 0 ) THEN
        CALL crs_dom_coordinates( gphit, glamt, 'T', gphit_crs, glamt_crs ) 
        CALL crs_dom_coordinates( gphiu, glamu, 'U', gphiu_crs, glamu_crs )       
        CALL crs_dom_coordinates( gphiv, glamv, 'V', gphiv_crs, glamv_crs ) 
        CALL crs_dom_coordinates( gphif, glamf, 'F', gphif_crs, glamf_crs ) 
     ELSEIF ( nresty /= 0 .AND. nrestx == 0 ) THEN
        CALL crs_dom_coordinates( gphiu, glamu, 'T', gphit_crs, glamt_crs )
        CALL crs_dom_coordinates( gphiu, glamu, 'U', gphiu_crs, glamu_crs )
        CALL crs_dom_coordinates( gphif, glamf, 'V', gphiv_crs, glamv_crs )
        CALL crs_dom_coordinates( gphif, glamf, 'F', gphif_crs, glamf_crs )
     ELSEIF ( nresty == 0 .AND. nrestx /= 0 ) THEN
        CALL crs_dom_coordinates( gphiv, glamv, 'T', gphit_crs, glamt_crs )
        CALL crs_dom_coordinates( gphif, glamf, 'U', gphiu_crs, glamu_crs )
        CALL crs_dom_coordinates( gphiv, glamv, 'V', gphiv_crs, glamv_crs )
        CALL crs_dom_coordinates( gphif, glamf, 'F', gphif_crs, glamf_crs )
     ELSE 
        CALL crs_dom_coordinates( gphif, glamf, 'T', gphit_crs, glamt_crs )
        CALL crs_dom_coordinates( gphif, glamf, 'U', gphiu_crs, glamu_crs )
        CALL crs_dom_coordinates( gphif, glamf, 'V', gphiv_crs, glamv_crs )
        CALL crs_dom_coordinates( gphif, glamf, 'F', gphif_crs, glamf_crs )
     ENDIF


     !  3.c. Get the horizontal mesh

     !      3.c.1 Horizontal scale factors

     CALL crs_dom_hgr( e1t, e2t, 'T', e1t_crs, e2t_crs )
     CALL crs_dom_hgr( e1u, e2u, 'U', e1u_crs, e2u_crs )
     CALL crs_dom_hgr( e1v, e2v, 'V', e1v_crs, e2v_crs )
     CALL crs_dom_hgr( e1f, e2f, 'F', e1f_crs, e2f_crs )

     e1e2t_crs(:,:) = e1t_crs(:,:) * e2t_crs(:,:)
     
     
     !      3.c.2 Coriolis factor  

      SELECT CASE( jphgr_msh )   ! type of horizontal mesh

      CASE ( 0, 1, 4 )           ! mesh on the sphere

         ff_crs(:,:) = 2. * omega * SIN( rad * gphif_crs(:,:) )

      CASE DEFAULT 

       IF(lwp)    WRITE(numout,*) 'crsini.F90. crs_init. Only jphgr_msh = 0, 1 or 4 supported' 
 
      END SELECT 

     !    3.d.1 mbathy ( vertical k-levels of bathymetry )     

     CALL crs_dom_bat
     
     !
     CALL wrk_alloc(jpi, jpj, jpk, zfse3t, zfse3u, zfse3v, zfse3w )
     !
     zfse3t(:,:,:) = e3t_0(:,:,:)
     zfse3u(:,:,:) = e3u_0(:,:,:)
     zfse3v(:,:,:) = e3v_0(:,:,:)
     zfse3w(:,:,:) = e3w_0(:,:,:)

     !    3.d.2   Surfaces 
     CALL crs_dom_sfc( tmask, 'W', e1e2w_crs, e1e2w_msk, p_e1=e1t, p_e2=e2t    )
     CALL crs_dom_sfc( umask, 'U', e2e3u_crs, e2e3u_msk, p_e2=e2u, p_e3=zfse3u )
     CALL crs_dom_sfc( vmask, 'V', e1e3v_crs, e1e3v_msk, p_e1=e1v, p_e3=zfse3v )
   
     facsurfu(:,:,:) = umask_crs(:,:,:) * e2e3u_msk(:,:,:) / e2e3u_crs(:,:,:)
     facsurfv(:,:,:) = vmask_crs(:,:,:) * e1e3v_msk(:,:,:) / e1e3v_crs(:,:,:)

     !    3.d.3   Vertical scale factors
     !
   
  
     CALL crs_dom_e3( e1t, e2t, zfse3t, e1e2w_crs, 'T', tmask, e3t_crs, e3t_max_crs)
     CALL crs_dom_e3( e1u, e2u, zfse3u, e2e3u_crs, 'U', umask, e3u_crs, e3u_max_crs)
     CALL crs_dom_e3( e1v, e2v, zfse3v, e1e3v_crs, 'V', vmask, e3v_crs, e3v_max_crs)
     CALL crs_dom_e3( e1t, e2t, zfse3w, e1e2w_crs, 'W', tmask, e3w_crs, e3w_max_crs)

     ! Reset 0 to e3t_0 or e3w_0
     DO jk = 1, jpk
        DO ji = 1, jpi_crs
           DO jj = 1, jpj_crs
              IF( e3t_crs(ji,jj,jk) == 0._wp ) e3t_crs(ji,jj,jk) = e3t_1d(jk)
              IF( e3w_crs(ji,jj,jk) == 0._wp ) e3w_crs(ji,jj,jk) = e3w_1d(jk)
              IF( e3u_crs(ji,jj,jk) == 0._wp ) e3u_crs(ji,jj,jk) = e3t_1d(jk)
              IF( e3v_crs(ji,jj,jk) == 0._wp ) e3v_crs(ji,jj,jk) = e3t_1d(jk)
           ENDDO
        ENDDO
     ENDDO

     !    3.d.3   Vertical depth (meters)
     CALL crs_dom_ope( gdept_0, 'MAX', 'T', tmask, gdept_crs, p_e3=zfse3t, psgn=1.0 ) 
     CALL crs_dom_ope( gdepw_0, 'MAX', 'W', tmask, gdepw_crs, p_e3=zfse3w, psgn=1.0 )


     !---------------------------------------------------------
     ! 4.  Coarse grid ocean volume and averaging weights
     !---------------------------------------------------------
     ! 4.a. Ocean volume or area unmasked and masked
     CALL crs_dom_facvol( tmask, 'T', e1t, e2t, zfse3t, ocean_volume_crs_t, facvol_t )
     !
     bt_crs(:,:,:) = ocean_volume_crs_t(:,:,:) * facvol_t(:,:,:)
     !
     r1_bt_crs(:,:,:) = 0._wp 
     WHERE( bt_crs /= 0._wp ) r1_bt_crs(:,:,:) = 1._wp / bt_crs(:,:,:)

     CALL crs_dom_facvol( tmask, 'W', e1t, e2t, zfse3w, ocean_volume_crs_w, facvol_w )
     !
     !---------------------------------------------------------
     ! 5.  Write out coarse meshmask  (see OPA_SRC/DOM/domwri.F90 for ideas later)
     !---------------------------------------------------------

     IF( nn_msh_crs > 0 ) THEN 
        CALL dom_grid_crs   ! Save the parent grid information  & Switch to coarse grid domain
        CALL crs_dom_wri     
        CALL dom_grid_glo   ! Return to parent grid domain
     ENDIF
     
     !---------------------------------------------------------
     ! 7. Finish and clean-up
     !---------------------------------------------------------
     CALL wrk_dealloc(jpi, jpj, jpk, zfse3t, zfse3u, zfse3v, zfse3w )


   END SUBROUTINE crs_init
    
   !!======================================================================

END MODULE crsini
