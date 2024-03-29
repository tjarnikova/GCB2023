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

MODULE limdia_2
   !!======================================================================
   !!                       ***  MODULE limdia_2   ***
   !!                      diagnostics of ice model 
   !!======================================================================
   !! History :   8.0  !  97-06  (Louvain-La-Neuve)  Original code
   !!             8.5  !  02-09  (C. Ethe , G. Madec )  F90: Free form and module
   !!             9.0  !  06-08  (S. Masson)  change frequency output control
   !!-------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_lim2' :                                  LIM 2.0 sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_dia_2      : computation of the time evolution of keys var.
   !!   lim_dia_init_2 : initialization and namelist read
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! 
   USE par_ice_2       ! ice parameters
   USE sbc_oce         ! surface boundary condition variables
   USE dom_ice_2       !
   USE ice_2           !
   USE limistate_2     !
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC               lim_dia_2          ! called by sbc_ice_lim_2
   PUBLIC               lim_dia_init_2     ! called by sbc_ice_lim_2

   INTEGER, PUBLIC ::   ntmoy  ,   &  !: instantaneous values of ice evolution or averaging ntmoy
      &                 ninfo         !: frequency of ouputs on file ice_evolu in case of averaging

   INTEGER, PARAMETER ::   &  ! Parameters for outputs to files "evolu"
      jpinfmx = 100         ,    &  ! maximum number of key variables
      jpchinf = 5           ,    &  ! ???
      jpchsep = jpchinf + 2         ! ???

   INTEGER ::   &
      nfrinf ,          &  ! number of variables written in one line 
      nferme ,          &  ! last time step at which the var. are written on file
      nvinfo ,          &  ! number of total variables 
      nbvt   ,          &  ! number of time variables
      naveg                ! number of step for accumulation before averaging

   CHARACTER(len= 8) ::   fmtinf   ! format of the output values  
   CHARACTER(len=30) ::   fmtw  ,           &  ! formats
      &                   fmtr  ,           &  ! ???
      &                   fmtitr               ! ???
   CHARACTER(len=jpchsep), DIMENSION(jpinfmx) ::   titvar               ! title of key variables
 
   REAL(wp)                     ::   epsi06 = 1.e-06      ! ???
   REAL(wp), DIMENSION(jpinfmx) ::   vinfom               ! temporary working space
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   aire                 ! masked grid cell area

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
   !! NEMO/LIM2 3.3 , UCL - NEMO Consortium (2010)
   !! $Id: limdia_2.F90 4624 2014-04-28 12:09:03Z acc $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE lim_dia_2( kt )
      !!--------------------------------------------------------------------
      !!                  ***  ROUTINE lim_dia_2  ***
      !!   
      !! ** Purpose : Computation and outputs on file ice.evolu 
      !!      the temporal evolution of some key variables
      !!-------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! number of iteration
      !!
      INTEGER  ::   jv,ji, jj   ! dummy loop indices
      INTEGER  ::   nv          ! indice of variable 
      REAL(wp) ::   zarea    , zldarea  ,    &  ! sea-ice and leads area
         &          zextent15, zextent85,    &  ! sea-ice extent (15% and 85%)
         &          zicevol  , zsnwvol  ,    &  ! sea-ice and snow volume volume
         &          zicespd                     ! sea-ice velocity
      REAL(wp), DIMENSION(jpinfmx) ::   vinfor  ! temporary working space 
      !!-------------------------------------------------------------------

      ! computation of key variables at each time step   

      nv = 1 
      vinfor(nv) = REAL( kt + nn_fsbc - 1 )
      nv = nv + 1
      vinfor(nv) = nyear
 
      DO jv = nbvt + 1, nvinfo
         vinfor(jv) = 0.e0
      END DO

      zextent15 = 0.e0
      zextent85 = 0.e0
      ! variables in northern Hemis
      DO jj = njeq, jpjm1
         DO ji = 2, jpim1   ! vector opt.
            IF( tms(ji,jj) == 1 ) THEN
               zarea = ( 1.0 - frld(ji,jj) ) * aire(ji,jj)
               IF (frld(ji,jj) <= 0.15 ) zextent15 = aire(ji,jj)    
               IF (frld(ji,jj) <= 0.85 ) zextent85 = aire(ji,jj)   
               zldarea = zarea   / MAX( ( 1 - frld(ji,jj) ) , epsi06 )
               zicevol = zarea   * hicif(ji,jj)
               zsnwvol = zarea   * hsnif(ji,jj)
               zicespd = zicevol * ( u_ice(ji,jj) * u_ice(ji,jj) + v_ice(ji,jj) * v_ice(ji,jj) )
               vinfor(nv+ 1) = vinfor(nv+ 1) + zarea
               vinfor(nv+ 3) = vinfor(nv+ 3) + zextent15
               vinfor(nv+ 5) = vinfor(nv+ 5) + zextent85
               vinfor(nv+ 7) = vinfor(nv+ 7) + zldarea
               vinfor(nv+ 9) = vinfor(nv+ 9) + zicevol
               vinfor(nv+11) = vinfor(nv+11) + zsnwvol
               vinfor(nv+13) = vinfor(nv+13) + zicespd
            ENDIF
         END DO
      END DO
      vinfor(nv+13) = SQRT( vinfor(nv+13) / MAX( vinfor(nv+9) , epsi06 ) )


      ! variables in southern Hemis
       nv = nv + 1
       DO jj = 2, njeqm1
          DO ji = 2, jpim1   ! vector opt.
             IF( tms(ji,jj) == 1 ) THEN
                zarea = ( 1.0 - frld(ji,jj) ) * aire(ji,jj)
                IF (frld(ji,jj) <= 0.15 ) zextent15 = aire(ji,jj)    
                IF (frld(ji,jj) <= 0.85 ) zextent85 = aire(ji,jj)   
                zldarea = zarea   / MAX( ( 1 - frld(ji,jj) ) , epsi06 )
                zicevol = zarea   * hicif(ji,jj)
                zsnwvol = zarea   * hsnif(ji,jj)
                zicespd = zicevol * ( u_ice(ji,jj) * u_ice(ji,jj) + v_ice(ji,jj) * v_ice(ji,jj) )
                vinfor(nv+ 1) = vinfor(nv+ 1) + zarea
                vinfor(nv+ 3) = vinfor(nv+ 3) + zextent15
                vinfor(nv+ 5) = vinfor(nv+ 5) + zextent85
                vinfor(nv+ 7) = vinfor(nv+ 7) + zldarea
                vinfor(nv+ 9) = vinfor(nv+ 9) + zicevol
                vinfor(nv+11) = vinfor(nv+11) + zsnwvol
                vinfor(nv+13) = vinfor(nv+13) + zicespd
             ENDIF
          END DO
       END DO
       vinfor(nv+13) = SQRT( vinfor(nv+13) / MAX( vinfor(nv+9) , epsi06 ) )    

       !  Accumulation before averaging 
       DO jv = 1, nvinfo
          vinfom(jv) = vinfom(jv) + vinfor(jv)
       END DO
       naveg = naveg + 1  
    
       ! oututs on file ice_evolu    
       IF( MOD( kt + nn_fsbc - 1, ninfo ) == 0 ) THEN
          WRITE(numevo_ice,fmtw) ( titvar(jv), vinfom(jv)/naveg, jv = 1, nvinfo )
          naveg = 0
          DO jv = 1, nvinfo
             vinfom(jv) = 0.e0
          END DO
       ENDIF
       !
    END SUBROUTINE lim_dia_2
 

    SUBROUTINE lim_dia_init_2
       !!-------------------------------------------------------------------
       !!                  ***  ROUTINE lim_dia_init_2  ***
       !!             
       !! ** Purpose : Preparation of the file ice_evolu for the output of
       !!      the temporal evolution of key variables
       !!
       !! ** input   : Namelist namicedia
       !!-------------------------------------------------------------------
       CHARACTER(len=jpchinf) ::   titinf
       INTEGER  ::   jv   ! dummy loop indice
       INTEGER  ::   ntot , ndeb, nv, ierr   ! local integer
       INTEGER  ::   ios                     ! Local integer output status for namelist read
       REAL(wp) ::   zxx0, zxx1              ! local scalars

       NAMELIST/namicedia/fmtinf, nfrinf, ninfo, ntmoy
       !!-------------------------------------------------------------------
                    
       REWIND( numnam_ice_ref )              ! Namelist namicedia in reference namelist : Ice diagnostics in ice_evolu
       READ  ( numnam_ice_ref, namicedia, IOSTAT = ios, ERR = 901)
901    IF( ios /= 0 ) CALL ctl_nam ( ios , 'namicedia in reference namelist', lwp )

       REWIND( numnam_ice_cfg )              ! Namelist namicedia in configuration namelist : Ice diagnostics in ice_evolu 
       READ  ( numnam_ice_cfg, namicedia, IOSTAT = ios, ERR = 902 )
902    IF( ios /= 0 ) CALL ctl_nam ( ios , 'namicedia in configuration namelist', lwp )
       IF(lwm) WRITE ( numoni, namicedia )
       
       IF(lwp) THEN                             ! control print
          WRITE(numout,*)
          WRITE(numout,*) 'lim_dia_init_2 : ice parameters for ice diagnostics '
          WRITE(numout,*) '~~~~~~~~~~~~~~'
          WRITE(numout,*) '   format of the output values                                 fmtinf = ', fmtinf
          WRITE(numout,*) '   number of variables written in one line                     nfrinf = ', nfrinf 
          WRITE(numout,*) '   Instantaneous values of ice evolution or averaging          ntmoy  = ', ntmoy
          WRITE(numout,*) '   frequency of ouputs on file ice_evolu in case of averaging  ninfo  = ', ninfo
       ENDIF

       ALLOCATE( aire(jpi,jpj) , STAT=ierr )    ! masked grid cell area
       IF( lk_mpp    )   CALL mpp_sum( ierr )
       IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'lim_dia_init_2 : unable to allocate standard arrays' )
       aire(:,:) = area(:,:) * tms(:,:)

       nv = 1                                   ! Titles of ice key variables
       titvar(nv) = 'NoIt'  ! iteration number
       nv = nv + 1
       titvar(nv) = 'T yr'  ! time step in years
       nbvt = nv - 1
       nv = nv + 1   ;   titvar(nv) = 'AEFN' ! sea ice area in the northern Hemisp.(10^12 km2)
       nv = nv + 1   ;   titvar(nv) = 'AEFS' ! sea ice area in the southern Hemisp.(10^12 km2)
       nv = nv + 1   ;   titvar(nv) = 'A15N'  ! sea ice extent (15%) in the northern Hemisp.(10^12 km2)
       nv = nv + 1   ;   titvar(nv) = 'A15S'  ! sea ice extent (15%) in the southern Hemisp.(10^12 km2)
       nv = nv + 1   ;   titvar(nv) = 'A85N'  ! sea ice extent (85%) in the northern Hemisp.(10^12 km2)
       nv = nv + 1   ;   titvar(nv) = 'A85S'  ! sea ice extent (85%) in the southern Hemisp.(10^12 km2)
       nv = nv + 1   ;   titvar(nv) = 'ALEN'  ! leads area in the northern Hemisp.(10^12 km2)
       nv = nv + 1   ;   titvar(nv) = 'ALES'  ! leads area in the southern Hemisp.(10^12 km2)
       nv = nv + 1   ;   titvar(nv) = 'VOLN'  ! sea ice volume in the northern Hemisp.(10^3 km3)
       nv = nv + 1   ;   titvar(nv) = 'VOLS'  ! sea ice volume in the southern Hemisp.(10^3 km3)
       nv = nv + 1   ;   titvar(nv) = 'VONN'  ! snow volume over sea ice in the northern Hemisp.(10^3 km3)
       nv = nv + 1   ;   titvar(nv) = 'VONS'  ! snow volume over sea ice in the southern Hemisp.(10^3 km3)
       nv = nv + 1   ;   titvar(nv) = 'ECGN'  ! mean sea ice velocity in the northern Hemisp.(m/s)
       nv = nv + 1   ;   titvar(nv) = 'ECGS'  ! mean sea ice velocity in the southern Hemisp.(m/s)

       nvinfo = nv

       ! Definition et Ecriture de l'entete : nombre d'enregistrements 
       ndeb   = ( nit000 - 1 + nn_fsbc - 1 ) / ninfo
       IF( nit000 - 1 + nn_fsbc == 1 ) ndeb = -1

       nferme = ( nitend + nn_fsbc - 1 ) / ninfo ! nit000 - 1 + nn_fsbc - 1 + nitend - nit000 + 1
       ntot   = nferme - ndeb
       ndeb   = ninfo * ( 1 + ndeb )
       nferme = ninfo * nferme

       ! definition of formats 
       WRITE( fmtw  , '(A,I3,A2,I1,A)' )  '(', nfrinf, '(A', jpchsep, ','//fmtinf//'))'
       WRITE( fmtr  , '(A,I3,A,I1,A)'  )  '(', nfrinf, '(', jpchsep, 'X,'//fmtinf//'))'
       WRITE( fmtitr, '(A,I3,A,I1,A)'  )  '(', nvinfo, 'A', jpchinf, ')'

       ! opening  "ice_evolu" file
       CALL ctl_opn( numevo_ice, 'ice_evolu', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp )

       !- ecriture de 2 lignes d''entete :
       WRITE(numevo_ice,1000) fmtr, fmtw, fmtitr, nvinfo, ntot, 0, nfrinf
       zxx0 = 0.001 * REAL( ninfo )
       zxx1 = 0.001 * REAL( ndeb  )
       WRITE(numevo_ice,1111) REAL(jpchinf), 0., zxx1, zxx0, 0., 0., 0

       !- ecriture de 2 lignes de titre :
       WRITE(numevo_ice,'(A,I8,A,I8,A,I5)')                 &
          'Evolution chronologique - Experience '//cexper   &
          //'   de', ndeb, ' a', nferme, ' pas', ninfo
       WRITE(numevo_ice,fmtitr) ( titvar(jv), jv = 1, nvinfo )


       !--preparation de "titvar" pour l''ecriture parmi les valeurs numeriques :
       DO  jv = 2 , nvinfo
          titinf     = titvar(jv)(:jpchinf)
          titvar(jv) = '  '//titinf
       END DO

       !--Initialisation of the arrays for the accumulation
       DO  jv = 1, nvinfo
          vinfom(jv) = 0.
       END DO
       naveg = 0

1000   FORMAT( 3(A20),4(1x,I6) )
1111   FORMAT( 3(F7.1,1X,F7.3,1X),I3,A )  
      !
    END SUBROUTINE lim_dia_init_2


   !!======================================================================
END MODULE limdia_2
