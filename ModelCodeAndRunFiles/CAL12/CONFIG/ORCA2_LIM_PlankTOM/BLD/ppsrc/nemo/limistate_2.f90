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

MODULE limistate_2
   !!======================================================================
   !!                     ***  MODULE  limistate_2  ***
   !!              Initialisation of diagnostics ice variables
   !!======================================================================
   !! History :   1.0  !  01-04  (C. Ethe, G. Madec)  Original code
   !!             2.0  !  03-08  (G. Madec)  add lim_istate_init
   !!                  !  04-04  (S. Theetten) initialization from a file
   !!                  !  06-07  (S. Masson)  IOM to read the restart
   !!                  !  07-10  (G. Madec)  surface module
   !!--------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_lim2' :                                  LIM 2.0 sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_istate_2      :  Initialisation of diagnostics ice variables
   !!   lim_istate_init_2 :  initialization of ice state and namelist read
   !!----------------------------------------------------------------------
   USE phycst
   USE par_ice_2       ! ice parameters
   USE dom_ice_2
   USE eosbn2          ! equation of state
   USE lbclnk
   USE oce
   USE ice_2
   USE iom
   USE in_out_manager
   USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC lim_istate_2      ! routine called by lim_init_2.F90

   !                        !! **  namelist (namiceini) **
   LOGICAL  ::   ln_limini   ! Ice initialization state
   REAL(wp) ::   ttest       ! threshold water temperature for initial sea ice
   REAL(wp) ::   hninn       ! initial snow thickness in the north
   REAL(wp) ::   hginn       ! initial ice thickness in the north
   REAL(wp) ::   alinn       ! initial leads area in the north
   REAL(wp) ::   hnins       ! initial snow thickness in the south
   REAL(wp) ::   hgins       ! initial ice thickness in the south
   REAL(wp) ::   alins       ! initial leads area in the south
   
   REAL(wp) ::   zero      = 0.e0     ! constant value = 0
   REAL(wp) ::   zone      = 1.e0     ! constant value = 1
   !!----------------------------------------------------------------------
   !! NEMO/LIM2 3.3 , UCL - NEMO Consortium (2010)
   !! $Id: limistate_2.F90 5540 2015-07-02 15:11:23Z jchanut $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE lim_istate_2
      !!-------------------------------------------------------------------
      !!                    ***  ROUTINE lim_istate_2  ***
      !!
      !! ** Purpose :   defined the sea-ice initial state
      !!
      !! ** Method  :   restart from a state defined in a binary file
      !!                or from arbitrary sea-ice conditions
      !!--------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk                ! dummy loop indices
      REAL(wp) ::   zidto                     ! temporary scalar
      !--------------------------------------------------------------------
 
      CALL lim_istate_init_2     !  reading the initials parameters of the ice

      IF( .NOT. ln_limini ) THEN  
         
         CALL eos_fzp( tsn(:,:,1,jp_sal), tfu(:,:) )       ! freezing/melting point of sea water [Celcius]
         tfu(:,:) = tfu(:,:) *  tmask(:,:,1)

         DO jj = 1, jpj
            DO ji = 1, jpi
               !                     ! ice if sst <= t-freez + ttest
               IF( tsn(ji,jj,1,jp_tem)  - tfu(ji,jj) >= ttest ) THEN   ;   zidto = 0.e0      ! no ice
               ELSE                                                    ;   zidto = 1.e0      !    ice
               ENDIF
               !
               IF( fcor(ji,jj) >= 0.e0 ) THEN     !--  Northern hemisphere.
                  hicif(ji,jj)   = zidto * hginn
                  frld(ji,jj)    = zidto * alinn + ( 1.0 - zidto ) * 1.0
                  hsnif(ji,jj)   = zidto * hninn
               ELSE                               !---  Southern hemisphere.
                  hicif(ji,jj)   = zidto * hgins
                  frld(ji,jj)    = zidto * alins + ( 1.0 - zidto ) * 1.0
                  hsnif(ji,jj)   = zidto * hnins
               ENDIF
            END DO
         END DO

         tfu(:,:) = tfu(:,:) + rt0       ! ftu converted from Celsius to Kelvin (rt0 over land)
         
         sist  (:,:)   = tfu(:,:)
         tbif  (:,:,1) = tfu(:,:)
         tbif  (:,:,2) = tfu(:,:)
         tbif  (:,:,3) = tfu(:,:)

      ENDIF
     
      fsbbq (:,:)   = 0.e0
      qstoif(:,:)   = 0.e0
      u_ice (:,:)   = 0.e0
      v_ice (:,:)   = 0.e0

      !---  Moments for advection.             

      sxice (:,:)  = 0.e0   ;   sxsn (:,:)  = 0.e0   ;   sxa  (:,:)  = 0.e0
      syice (:,:)  = 0.e0   ;   sysn (:,:)  = 0.e0   ;   sya  (:,:)  = 0.e0
      sxxice(:,:)  = 0.e0   ;   sxxsn(:,:)  = 0.e0   ;   sxxa (:,:)  = 0.e0
      syyice(:,:)  = 0.e0   ;   syysn(:,:)  = 0.e0   ;   syya (:,:)  = 0.e0
      sxyice(:,:)  = 0.e0   ;   sxysn(:,:)  = 0.e0   ;   sxya (:,:)  = 0.e0

      sxc0  (:,:)  = 0.e0   ;   sxc1 (:,:)  = 0.e0   ;   sxc2 (:,:)  = 0.e0
      syc0  (:,:)  = 0.e0   ;   syc1 (:,:)  = 0.e0   ;   syc2 (:,:)  = 0.e0
      sxxc0 (:,:)  = 0.e0   ;   sxxc1(:,:)  = 0.e0   ;   sxxc2(:,:)  = 0.e0
      syyc0 (:,:)  = 0.e0   ;   syyc1(:,:)  = 0.e0   ;   syyc2(:,:)  = 0.e0
      sxyc0 (:,:)  = 0.e0   ;   sxyc1(:,:)  = 0.e0   ;   sxyc2(:,:)  = 0.e0

      sxst  (:,:)  = 0.e0
      syst  (:,:)  = 0.e0
      sxxst (:,:)  = 0.e0
      syyst (:,:)  = 0.e0
      sxyst (:,:)  = 0.e0
      stress1_i (:,:) = 0._wp                          ! EVP rheology
      stress2_i (:,:) = 0._wp
      stress12_i(:,:) = 0._wp

      !-- lateral boundary conditions
      CALL lbc_lnk( hicif, 'T', 1. )
      CALL lbc_lnk( frld , 'T', 1. )

      ! C A U T I O N  frld = 1 over land and lbc_lnk put zero along 
      ! *************  closed boundaries herefore we force to one over land
      frld(:,:) = tms(:,:) * frld(:,:) + ( 1. - tms(:,:) )   

      CALL lbc_lnk( hsnif, 'T', 1. )
      CALL lbc_lnk( sist , 'T', 1. , pval = rt0 )      ! set rt0 on closed boundary (required by bulk formulation)
      DO jk = 1, jplayersp1
         CALL lbc_lnk(tbif(:,:,jk), 'T', 1. )
      END DO
      CALL lbc_lnk( fsbbq  , 'T', 1. )
      CALL lbc_lnk( qstoif , 'T', 1. )

   END SUBROUTINE lim_istate_2

   
   SUBROUTINE lim_istate_init_2
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE lim_istate_init_2  ***
      !!        
      !! ** Purpose :   Definition of initial state of the ice 
      !!
      !! ** Method  :   Read the namiceini namelist and check the parameter 
      !!       values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namiceini
      !!-------------------------------------------------------------------
      INTEGER :: inum_ice
      INTEGER :: ji,jj
      INTEGER :: ios                 ! Local integer output status for namelist read

      NAMELIST/namiceini/ ln_limini, ttest, hninn, hginn, alinn, &
         &                hnins, hgins, alins
      !!-------------------------------------------------------------------
                   
      REWIND( numnam_ice_ref )              ! Namelist namiceini in reference namelist : Ice initial state
      READ  ( numnam_ice_ref, namiceini, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namiceini in reference namelist', lwp )

      REWIND( numnam_ice_cfg )              ! Namelist namiceini in configuration namelist : Ice initial state
      READ  ( numnam_ice_cfg, namiceini, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namiceini in configuration namelist', lwp )
      IF(lwm) WRITE ( numoni, namiceini )
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'lim_istate_init_2 : ice parameters inititialisation '
         WRITE(numout,*) '~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '         threshold water temp. for initial sea-ice    ttest      = ', ttest
         WRITE(numout,*) '         initial snow thickness in the north          hninn      = ', hninn
         WRITE(numout,*) '         initial ice thickness in the north           hginn      = ', hginn 
         WRITE(numout,*) '         initial leads area in the north              alinn      = ', alinn            
         WRITE(numout,*) '         initial snow thickness in the south          hnins      = ', hnins 
         WRITE(numout,*) '         initial ice thickness in the south           hgins      = ', hgins
         WRITE(numout,*) '         initial leads area in the south              alins      = ', alins
         WRITE(numout,*) '         Ice state initialization using input file    ln_limini  = ', ln_limini
      ENDIF

      IF( ln_limini ) THEN                      ! Ice initialization using input file
         !
         CALL iom_open( 'Ice_initialization.nc', inum_ice )
         !
         IF( inum_ice > 0 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '                  ice state initialization with : Ice_initialization.nc'
            
            CALL iom_get( inum_ice, jpdom_data, 'hicif', hicif )      
            CALL iom_get( inum_ice, jpdom_data, 'hsnif', hsnif )      
            CALL iom_get( inum_ice, jpdom_data, 'frld' , frld  )     
            CALL iom_get( inum_ice, jpdom_data, 'ts'   , sist  )
            CALL iom_get( inum_ice, jpdom_unknown, 'tbif', tbif(1:nlci,1:nlcj,:),   &
                 &        kstart = (/ mig(1),mjg(1),1 /), kcount = (/ nlci,nlcj,jplayersp1 /) )
            ! put some values in the extra-halo...
            DO jj = nlcj+1, jpj   ;   tbif(1:nlci,jj,:) = tbif(1:nlci,nlej,:)   ;   END DO
            DO ji = nlci+1, jpi   ;   tbif(ji    ,: ,:) = tbif(nlei  ,:   ,:)   ;   END DO

            CALL iom_close( inum_ice)
            !
         ENDIF
      ENDIF
      !     
   END SUBROUTINE lim_istate_init_2


   !!======================================================================
END MODULE limistate_2
