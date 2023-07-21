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

MODULE trcnam
   !!======================================================================
   !!                       ***  MODULE trcnam  ***
   !! TOP :   Read and print options for the passive tracer run (namelist)
   !!======================================================================
   !! History :    -   !  1996-11  (M.A. Foujols, M. Levy)  original code
   !!              -   !  1998-04  (M.A Foujols, L. Bopp) ahtrb0 for isopycnal mixing
   !!              -   !  1999-10  (M.A. Foujols, M. Levy) separation of sms
   !!              -   !  2000-07  (A. Estublier) add TVD and MUSCL : Tests on ndttrc
   !!              -   !  2000-11  (M.A Foujols, E Kestenare) trcrat, ahtrc0 and aeivtr0
   !!              -   !  2001-01 (E Kestenare) suppress ndttrc=1 for CEN2 and TVD schemes
   !!             1.0  !  2005-03 (O. Aumont, A. El Moussaoui) F90
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_nam    :  Read and print options for the passive tracer run (namelist)
   !!----------------------------------------------------------------------
   USE oce_trc           ! shared variables between ocean and passive tracers
   USE trc               ! passive tracers common variables
   USE trcnam_trp        ! Transport namelist
   USE trcnam_pisces     ! PISCES namelist
   USE trcnam_cfc        ! CFC SMS namelist
   USE trcnam_c14b       ! C14 SMS namelist
   USE trcnam_my_trc     ! MY_TRC SMS namelist
   USE trd_oce       
   USE trdtrc_oce
   USE iom               ! I/O manager
   USE trcnam_planktom   ! PlankTOM namelist
   USE par_planktom, only : lk_planktom 

   IMPLICIT NONE
   PRIVATE 

   PUBLIC trc_nam_run  ! called in trcini
   PUBLIC trc_nam      ! called in trcini

   !! * Substitutions
   !!----------------------------------------------------------------------
   !!                    ***  top_substitute.h90   ***
   !!----------------------------------------------------------------------
   !! ** purpose : Statement function file: to be include in all passive tracer modules
   !!----------------------------------------------------------------------
   !! History :   1.0  !  2004-03 (C. Ethe) Original code
   !!             2.0  !  2007-12 (C. Ethe, G. Madec) new architecture
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
   !!                   ***  ldfeiv_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsaei. the eddy induced velocity coeff.
   !!      with a constant or 1D or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldfeiv_substitute.h90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
!   'traldf_c2d' :                           eiv: 2D coefficient
   !!----------------------------------------------------------------------
   !!                    *** ldftra_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
   !!      with a constant or 1D or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldftra_substitute.h90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
!   'key_traldf_c2d' :                 aht: 2D coefficient
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
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: top_substitute.h90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnam.F90 5411 2015-06-12 17:39:14Z cetlod $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS


   SUBROUTINE trc_nam
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_nam  ***
      !!
      !! ** Purpose :   READ and PRINT options for the passive tracer run (namelist) 
      !!
      !! ** Method  : - read passive tracer namelist 
      !!              - read namelist of each defined SMS model
      !!                ( (PISCES, CFC, MY_TRC )
      !!---------------------------------------------------------------------
      INTEGER  ::   jn                  ! dummy loop indice
      !                                        !   Parameters of the run 
      IF( .NOT. lk_offline ) CALL trc_nam_run
      
      !                                        !  passive tracer informations
      CALL trc_nam_trc
      
      !                                        !   Parameters of additional diagnostics
      CALL trc_nam_dia

      !                                        !   namelist of transport
      CALL trc_nam_trp


      IF( ln_rsttr )                      ln_trcdta = .FALSE.   ! restart : no need of clim data
      !
      IF( ln_trcdmp .OR. ln_trcdmp_clo )  ln_trcdta = .TRUE.   ! damping : need to have clim data
      !
      IF( .NOT.ln_trcdta ) THEN
         ln_trc_ini(:) = .FALSE.
      ENDIF

     IF(lwp) THEN                   ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist : namtrc'
         WRITE(numout,*) '   Read inputs data from file (y/n)             ln_trcdta     = ', ln_trcdta
         WRITE(numout,*) '   Damping of passive tracer (y/n)              ln_trcdmp     = ', ln_trcdmp
         WRITE(numout,*) '   Restoring of tracer on closed seas           ln_trcdmp_clo = ', ln_trcdmp_clo
         WRITE(numout,*) ' '
         DO jn = 1, jptra
            WRITE(numout,*) '  tracer nb : ', jn, '    short name : ', ctrcnm(jn)
         END DO
         WRITE(numout,*) ' '
      ENDIF

      IF(lwp) THEN                   ! control print
         IF( ln_rsttr ) THEN
            WRITE(numout,*)
            WRITE(numout,*) '  Read a restart file for passive tracer : ', TRIM( cn_trcrst_in )
            WRITE(numout,*)
         ENDIF
         IF( ln_trcdta .AND. .NOT.ln_rsttr ) THEN
            WRITE(numout,*)
            WRITE(numout,*) '  Some of the passive tracers are initialised from climatologies '
            WRITE(numout,*)
         ENDIF
         IF( .NOT.ln_trcdta ) THEN
            WRITE(numout,*)
            WRITE(numout,*) '  All the passive tracers are initialised with constant values '
            WRITE(numout,*)
         ENDIF
      ENDIF

      
      rdttrc(:) = rdttra(:) * FLOAT( nn_dttrc )   ! vertical profile of passive tracer time-step
  
      IF(lwp) THEN                   ! control print
        WRITE(numout,*) 
        WRITE(numout,*) '    Passive Tracer  time step    rdttrc  = ', rdttrc(1)
        WRITE(numout,*) 
        CALL FLUSH(numout)
      ENDIF
      ! Call the ice module for tracers
      ! -------------------------------
      CALL trc_nam_ice

      ! namelist of SMS
      ! ---------------      
      IF( lk_pisces  ) THEN   ;   CALL trc_nam_pisces      ! PISCES  bio-model
      ELSE                    ;   IF(lwp) WRITE(numout,*) '          PISCES not used'
      ENDIF

      IF( lk_cfc     ) THEN   ;   CALL trc_nam_cfc         ! CFC     tracers
      ELSE                    ;   IF(lwp) WRITE(numout,*) '          CFC not used'
      ENDIF

      IF( lk_c14b     ) THEN   ;   CALL trc_nam_c14b         ! C14 bomb     tracers
      ELSE                    ;   IF(lwp) WRITE(numout,*) '          C14 not used'
      ENDIF

      IF( lk_my_trc  ) THEN   ;   CALL trc_nam_my_trc      ! MY_TRC  tracers
      ELSE                    ;   IF(lwp) WRITE(numout,*) '          MY_TRC not used'
      IF(lwp) CALL FLUSH(numout)
      ENDIF
      IF( lk_planktom ) CALL trc_nam_planktom
      !
   END SUBROUTINE trc_nam

   SUBROUTINE trc_nam_run
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_nam  ***
      !!
      !! ** Purpose :   read options for the passive tracer run (namelist) 
      !!
      !!---------------------------------------------------------------------
      NAMELIST/namtrc_run/ nn_dttrc, nn_writetrc, ln_rsttr, nn_rsttr, ln_top_euler, &
        &                  cn_trcrst_indir, cn_trcrst_outdir, cn_trcrst_in, cn_trcrst_out


      INTEGER  ::   ios                 ! Local integer output status for namelist read

      !!---------------------------------------------------------------------


      IF(lwp) WRITE(numout,*) 'trc_nam : read the passive tracer namelists'
      IF(lwp) WRITE(numout,*) '~~~~~~~'
      IF(lwp) CALL FLUSH(numout)

      CALL ctl_opn( numnat_ref, 'namelist_top_ref'   , 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      CALL ctl_opn( numnat_cfg, 'namelist_top_cfg'   , 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      IF(lwm) CALL ctl_opn( numont, 'output.namelist.top', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE., 1 )

      REWIND( numnat_ref )              ! Namelist namtrc in reference namelist : Passive tracer variables
      READ  ( numnat_ref, namtrc_run, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc in reference namelist', lwp )

      REWIND( numnat_cfg )              ! Namelist namtrc in configuration namelist : Passive tracer variables
      READ  ( numnat_cfg, namtrc_run, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc in configuration namelist', lwp )
      IF(lwm) WRITE ( numont, namtrc_run )

      !  computes the first time step of tracer model
      nittrc000 = nit000 + nn_dttrc - 1

      IF(lwp) THEN                   ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist : namtrc_run'
         WRITE(numout,*) '   time step freq. for passive tracer           nn_dttrc      = ', nn_dttrc
         WRITE(numout,*) '   restart  for passive tracer                  ln_rsttr      = ', ln_rsttr
         WRITE(numout,*) '   control of time step for passive tracer      nn_rsttr      = ', nn_rsttr
         WRITE(numout,*) '   first time step for pass. trac.              nittrc000     = ', nittrc000
         WRITE(numout,*) '   frequency of outputs for passive tracers     nn_writetrc   = ', nn_writetrc 
         WRITE(numout,*) '   Use euler integration for TRC (y/n)          ln_top_euler  = ', ln_top_euler
         WRITE(numout,*) ' '
         CALL FLUSH(numout)
      ENDIF
      !
    END SUBROUTINE trc_nam_run

   SUBROUTINE trc_nam_ice
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_nam_ice ***
      !!
      !! ** Purpose :   Read the namelist for the ice effect on tracers
      !!
      !! ** Method  : -
      !!
      !!---------------------------------------------------------------------
      ! --- Variable declarations --- !
      INTEGER :: jn      ! dummy loop indices
      INTEGER :: ios     ! Local integer output status for namelist read

      ! --- Namelist declarations --- !
      TYPE(TRC_I_NML), DIMENSION(jptra) :: sn_tri_tracer
      NAMELIST/namtrc_ice/ nn_ice_tr, sn_tri_tracer

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_nam_ice : Read the namelist for trc_ice'
         WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      ENDIF

      IF( nn_timing == 1 )  CALL timing_start('trc_nam_ice')

      !
      REWIND( numnat_ref )              ! Namelist namtrc_ice in reference namelist : Passive tracer input data
      READ  ( numnat_ref, namtrc_ice, IOSTAT = ios, ERR = 901)
 901  IF( ios /= 0 ) CALL ctl_nam ( ios , ' namtrc_ice in reference namelist ', lwp )

      REWIND( numnat_cfg )              ! Namelist namtrc_ice in configuration namelist : Pisces external sources of nutrients
      READ  ( numnat_cfg, namtrc_ice, IOSTAT = ios, ERR = 902 )
 902  IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc_ice in configuration namelist', lwp )

      IF( lwp ) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Sea ice tracers option (nn_ice_tr) : ', nn_ice_tr
         WRITE(numout,*) ' '
         CALL FLUSH(numout)
      ENDIF

      ! Assign namelist stuff
      DO jn = 1, jptra
         trc_ice_ratio(jn)  = sn_tri_tracer(jn)%trc_ratio
         trc_ice_prescr(jn) = sn_tri_tracer(jn)%trc_prescr
         cn_trc_o      (jn) = sn_tri_tracer(jn)%ctrc_o
      END DO

      IF( nn_timing == 1 )   CALL timing_stop('trc_nam_ice')
      !
   END SUBROUTINE trc_nam_ice

   SUBROUTINE trc_nam_trc
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_nam  ***
      !!
      !! ** Purpose :   read options for the passive tracer run (namelist) 
      !!
      !!---------------------------------------------------------------------
      TYPE(PTRACER), DIMENSION(jptra) :: sn_tracer  ! type of tracer for saving if not 1
      !!
      NAMELIST/namtrc/ sn_tracer, ln_trcdta,ln_trcdmp, ln_trcdmp_clo
  
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      INTEGER  ::   jn                  ! dummy loop indice
      !!---------------------------------------------------------------------
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'trc_nam : read the passive tracer namelists'
      IF(lwp) WRITE(numout,*) '~~~~~~~'
      IF(lwp) CALL FLUSH(numout)


      REWIND( numnat_ref )              ! Namelist namtrc in reference namelist : Passive tracer variables
      READ  ( numnat_ref, namtrc, IOSTAT = ios)!, ERR = 901)
!      IF(lwp) WRITE(numout,*) 'trc_nam ', ln_trcdta,ln_trcdmp, ln_trcdmp_clo,jptra,sn_tracer
!901   
      IF( ios /= 0 ) THEN
       IF(lwp) WRITE(numout,*) 'trc_nam ', ln_trcdta,ln_trcdmp, ln_trcdmp_clo,jptra,sn_tracer 
       CALL ctl_nam ( ios , 'namtrc in reference namelist', lwp )
      ENDIF

      REWIND( numnat_cfg )              ! Namelist namtrc in configuration namelist : Passive tracer variables
      READ  ( numnat_cfg, namtrc, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc in configuration namelist', lwp )
      IF(lwm) WRITE ( numont, namtrc )

      DO jn = 1, jptra
         ctrcnm    (jn) = TRIM( sn_tracer(jn)%clsname )
         ctrcln    (jn) = TRIM( sn_tracer(jn)%cllname )
         ctrcun    (jn) = TRIM( sn_tracer(jn)%clunit  )
         ln_trc_ini(jn) =       sn_tracer(jn)%llinit
         ln_trc_wri(jn) =       sn_tracer(jn)%llsave
      END DO
      
    END SUBROUTINE trc_nam_trc


   SUBROUTINE trc_nam_dia
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_nam_dia  ***
      !!
      !! ** Purpose :   read options for the passive tracer diagnostics
      !!
      !! ** Method  : - read passive tracer namelist 
      !!              - read namelist of each defined SMS model
      !!                ( (PISCES, CFC, MY_TRC )
      !!---------------------------------------------------------------------
      INTEGER ::  ierr,jn
      NAMELIST/namtrc_dia/ ln_diatrc, ln_diabio, nn_writedia, nn_writebio

      INTEGER  ::   ios                 ! Local integer output status for namelist read
      !!---------------------------------------------------------------------
      IF(lwp) WRITE(numout,*) 
      IF(lwp) WRITE(numout,*) 'trc_nam_dia : read the passive tracer diagnostics options'
      IF(lwp) WRITE(numout,*) '~~~~~~~'
      IF(lwp) CALL FLUSH(numout)

      REWIND( numnat_ref )              ! Namelist namtrc_dia in reference namelist : Passive tracer diagnostics
      READ  ( numnat_ref, namtrc_dia, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc_dia in reference namelist', lwp )

      REWIND( numnat_cfg )              ! Namelist namtrc_dia in configuration namelist : Passive tracer diagnostics
      READ  ( numnat_cfg, namtrc_dia, IOSTAT = ios, ERR = 904 )
904   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtrc_dia in configuration namelist', lwp )
      IF(lwm) WRITE ( numont, namtrc_dia )

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist : namtrc_dia'
         WRITE(numout,*) '    save additionnal diagnostics arrays         ln_diatrc   = ', ln_diatrc
         WRITE(numout,*) '    save additionnal biology diagnostics arrays ln_diabio   = ', ln_diabio
         WRITE(numout,*) '    frequency of outputs for additional arrays  nn_writedia = ', nn_writedia
         WRITE(numout,*) '    frequency of outputs for biological trends  nn_writebio = ', nn_writebio
         WRITE(numout,*) ' '
      ENDIF

      IF( ln_diatrc .AND. .NOT. lk_iomput ) THEN 
         ALLOCATE( trc2d(jpi,jpj,jpdia2d), trc3d(jpi,jpj,jpk,jpdia3d),  &
           &       ctrc2d(jpdia2d), ctrc2l(jpdia2d), ctrc2u(jpdia2d) ,  & 
           &       ctrc3d(jpdia3d), ctrc3l(jpdia3d), ctrc3u(jpdia3d) ,  STAT = ierr ) 
         IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'trcnam: unable to allocate add. diag. array' )
         !
         trc2d(:,:,:  ) = 0._wp  ;   ctrc2d(:) = ' '   ;   ctrc2l(:) = ' '    ;    ctrc2u(:) = ' ' 
         trc3d(:,:,:,:) = 0._wp  ;   ctrc3d(:) = ' '   ;   ctrc3l(:) = ' '    ;    ctrc3u(:) = ' ' 
         !
      ENDIF

      IF( ( ln_diabio .AND. .NOT. lk_iomput ) .OR. l_trdtrc ) THEN
         ALLOCATE( trbio (jpi,jpj,jpk,jpdiabio) , &
           &       ctrbio(jpdiabio), ctrbil(jpdiabio), ctrbiu(jpdiabio), STAT = ierr ) 
         IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'trcnam: unable to allocate bio. diag. array' )
         !
         trbio(:,:,:,:) = 0._wp  ;   ctrbio(:) = ' '   ;   ctrbil(:) = ' '    ;    ctrbiu(:) = ' ' 
         !
      ENDIF
      !
   END SUBROUTINE trc_nam_dia


   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnam.F90 5411 2015-06-12 17:39:14Z cetlod $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE  trcnam
