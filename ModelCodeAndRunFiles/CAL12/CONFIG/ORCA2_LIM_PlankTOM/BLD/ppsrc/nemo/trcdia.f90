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

MODULE trcdia
   !!======================================================================
   !!                       *** MODULE trcdia ***
   !! TOP :   Output of passive tracers
   !!======================================================================
   !! History :   OPA  !  1995-01 (M. Levy)  Original code
   !!              -   !  1998-01 (C. Levy) NETCDF format using ioipsl interface
   !!              -   !  1999-01 (M.A. Foujols) adapted for passive tracer
   !!              -   !  1999-09 (M.A. Foujols) split into three parts
   !!   NEMO      1.0  !  2005-03 (O. Aumont, A. El Moussaoui) F90
   !!                  !  2008-05 (C. Ethe re-organization)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !! trc_dia     : main routine of output passive tracer
   !! trcdit_wr   : outputs of concentration fields
   !! trcdii_wr   : outputs of additional 2D/3D diagnostics
   !! trcdib_wr   : outputs of biological fields
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain variables 
   USE oce_trc
   USE trc
   USE par_trc
   USE dianam    ! build name of file (routine)
   USE ioipsl    ! I/O manager
   USE iom       ! I/O manager
   USE lib_mpp   ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_dia        ! called by XXX module 

   INTEGER  ::   nit5      !: id for tracer output file
   INTEGER  ::   ndepit5   !: id for depth mesh
   INTEGER  ::   nhorit5   !: id for horizontal mesh
   INTEGER  ::   ndimt50   !: number of ocean points in index array
   INTEGER  ::   ndimt51   !: number of ocean points in index array
   REAL(wp) ::   zjulian   !: ????   not DOCTOR !
   INTEGER , ALLOCATABLE, SAVE, DIMENSION (:) ::   ndext50   !: integer arrays for ocean 3D index
   INTEGER , ALLOCATABLE, SAVE, DIMENSION (:) ::   ndext51   !: integer arrays for ocean surface index

   INTEGER  ::   nitd      !: id for additional array output file
   INTEGER  ::   ndepitd   !: id for depth mesh
   INTEGER  ::   nhoritd   !: id for horizontal mesh

   INTEGER  ::   nitb        !:         id.         for additional array output file
   INTEGER  ::   ndepitb   !:  id for depth mesh
   INTEGER  ::   nhoritb   !:  id for horizontal mesh

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
   !! $Id: trcdia.F90 4292 2013-11-20 16:28:04Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_dia( kt )  
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_dia  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! ocean time-step
      !
      INTEGER             ::  ierr   ! local integer
      !!---------------------------------------------------------------------
      !
      IF( kt == nittrc000 )  THEN
         ALLOCATE( ndext50(jpij*jpk), ndext51(jpij), STAT=ierr )
         IF( ierr > 0 ) THEN
            CALL ctl_stop( 'STOP', 'trc_diat: unable to allocate arrays' )  ;   RETURN
         ENDIF
      ENDIF
      !
      IF( .NOT.lk_iomput ) THEN
                          CALL trcdit_wr( kt )      ! outputs for tracer concentration
         IF( ln_diatrc )  CALL trcdii_wr( kt )      ! outputs for additional arrays
         IF( ln_diabio )  CALL trcdib_wr( kt )      ! outputs for biological trends
      ENDIF
      !
   END SUBROUTINE trc_dia


   SUBROUTINE trcdit_wr( kt )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE trcdit_wr  ***
      !!
      !! ** Purpose :   Standard output of passive tracer : concentration fields
      !!
      !! ** Method  :   At the beginning of the first time step (nittrc000), define all
      !!             the NETCDF files and fields for concentration of passive tracer
      !!
      !!        At each time step call histdef to compute the mean if necessary
      !!        Each nwritetrc time step, output the instantaneous or mean fields
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt       ! ocean time-step
      !
      INTEGER ::   jn
      LOGICAL ::   ll_print = .FALSE.
      CHARACTER (len=40) :: clhstnam, clop
      INTEGER ::   inum = 11             ! temporary logical unit
      CHARACTER (len=20) :: cltra, cltrau
      CHARACTER (len=80) :: cltral
      REAL(wp) :: zsto, zout, zdt
      INTEGER  :: iimi, iima, ijmi, ijma, ipk, it, itmod, iiter
      !!----------------------------------------------------------------------

      ! Initialisation
      ! --------------

      ! local variable for debugging
      ll_print = .FALSE.                  ! change it to true for more control print
      ll_print = ll_print .AND. lwp

      ! Define frequency of output and means
      zdt = rdt
      IF( ln_mskland )   THEN   ;   clop = "only(x)"   ! put 1.e+20 on land (very expensive!!)
      ELSE                      ;   clop = "x"         ! no use of the mask value (require less cpu time)
      ENDIF
      zsto = zdt
      clop = "ave("//TRIM(clop)//")"
      zout = nn_writetrc * zdt

      ! Define indices of the horizontal output zoom and vertical limit storage
      iimi = 1      ;      iima = jpi
      ijmi = 1      ;      ijma = jpj
      ipk = jpk

      ! define time axis
      itmod = kt - nittrc000 + 1
      it    = kt
      iiter = ( nittrc000 - 1 ) / nn_dttrc

      ! Define NETCDF files and fields at beginning of first time step
      ! --------------------------------------------------------------

      IF(ll_print)WRITE(numout,*)'trcdit_wr kt=',kt
      
      IF( kt == nittrc000 ) THEN

         IF(lwp) THEN                   ! control print
            WRITE(numout,*)
            WRITE(numout,*) '    frequency of outputs for passive tracers nn_writetrc = ', nn_writetrc
            DO jn = 1, jptra
               IF( ln_trc_wri(jn) )  WRITE(numout,*) ' ouput tracer nb : ', jn, '    short name : ', ctrcnm(jn) 
            END DO
            WRITE(numout,*) ' '
         ENDIF

         ! Compute julian date from starting date of the run
         CALL ymds2ju( nyear, nmonth, nday, rdt, zjulian )
         zjulian = zjulian - adatrj   !   set calendar origin to the beginning of the experiment
         IF(lwp)WRITE(numout,*)' '  
         IF(lwp)WRITE(numout,*)' Date 0 used :', nittrc000                         &
            &                 ,' YEAR ', nyear, ' MONTH ', nmonth, ' DAY ', nday   &
            &                 ,'Julian day : ', zjulian  
  
         IF(lwp) WRITE(numout,*) ' indexes of zoom = ', iimi, iima, ijmi, ijma,  &
            &                    ' limit storage in depth = ', ipk

         IF( lk_offline .AND. lwp ) THEN
            CALL dia_nam( clhstnam, nn_writetrc,' ' )
            CALL ctl_opn( inum, 'date.file', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', 1, numout, lwp, narea )
            WRITE(inum,*) clhstnam
            CLOSE(inum)
         ENDIF

         ! Define the NETCDF files for passive tracer concentration
         CALL dia_nam( clhstnam, nn_writetrc, 'ptrc_T' )
         IF(lwp)WRITE(numout,*)" Name of NETCDF file ", clhstnam

         ! Horizontal grid : glamt and gphit
         CALL histbeg( clhstnam, jpi, glamt, jpj, gphit,     &
            &          iimi, iima-iimi+1, ijmi, ijma-ijmi+1,         & 
            &          iiter, zjulian, zdt, nhorit5, nit5 , domain_id=nidom, snc4chunks=snc4set)

         ! Vertical grid for tracer : gdept
         CALL histvert( nit5, 'deptht', 'Vertical T levels', 'm', ipk, gdept_1d, ndepit5)

         ! Index of ocean points in 3D and 2D (surface)
         CALL wheneq( jpi*jpj*ipk, tmask, 1, 1., ndext50, ndimt50 )
         CALL wheneq( jpi*jpj    , tmask, 1, 1., ndext51, ndimt51 )

         ! Declare all the output fields as NETCDF variables
         DO jn = 1, jptra
            IF( ln_trc_wri(jn) ) THEN
               cltra  = TRIM( ctrcnm(jn) )   ! short title for tracer
               cltral = TRIM( ctrcln(jn) )   ! long title for tracer
               cltrau = TRIM( ctrcun(jn) )   ! UNIT for tracer
               CALL histdef( nit5, cltra, cltral, cltrau, jpi, jpj, nhorit5,  &
                  &          ipk, 1, ipk,  ndepit5, 32, clop, zsto, zout ) 
            ENDIF
         END DO

         ! end netcdf files header
         CALL histend( nit5, snc4set )
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'End of NetCDF Initialization in trcdit_wr'
         IF( ll_print )   CALL FLUSH(numout )

      ENDIF

      ! Start writing the tracer concentrations
      ! ---------------------------------------

      IF( lwp .AND. MOD( itmod, nn_writetrc ) == 0 ) THEN
         WRITE(numout,*) 'trcdit_wr : write NetCDF passive tracer concentrations at ', kt, 'time-step'
         WRITE(numout,*) '~~~~~~~~~ '
      ENDIF

      DO jn = 1, jptra
         cltra  = TRIM( ctrcnm(jn) )   ! short title for tracer
         IF( ln_trc_wri(jn) ) CALL histwrite( nit5, cltra, it, trn(:,:,:,jn), ndimt50, ndext50 )
      END DO

      ! close the file 
      ! --------------
      IF( kt == nitend )   CALL histclo( nit5 )
      !
   END SUBROUTINE trcdit_wr

   SUBROUTINE trcdii_wr( kt )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE trcdii_wr  ***
      !!
      !! ** Purpose :   output of passive tracer : additional 2D and 3D arrays
      !!
      !! ** Method  :   At the beginning of the first time step (nittrc000), define all
      !!             the NETCDF files and fields for concentration of passive tracer
      !!
      !!        At each time step call histdef to compute the mean if necessary
      !!        Each nn_writedia time step, output the instantaneous or mean fields
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt       ! ocean time-step
      !!
      LOGICAL ::   ll_print = .FALSE.
      CHARACTER (len=40) ::   clhstnam, clop
      CHARACTER (len=20) ::   cltra, cltrau
      CHARACTER (len=80) ::   cltral
      INTEGER  ::   jl
      INTEGER  ::   iimi, iima, ijmi, ijma, ipk, it, itmod, iiter
      REAL(wp) ::   zsto, zout, zdt
      !!----------------------------------------------------------------------

      ! Initialisation
      ! --------------
      
      ! local variable for debugging
      ll_print = .FALSE.
      ll_print = ll_print .AND. lwp
      !
      ! Define frequency of output and means
      zdt = rdt
      IF( ln_mskland )   THEN   ;   clop = "only(x)"   ! put 1.e+20 on land (very expensive!!)
      ELSE                      ;   clop = "x"         ! no use of the mask value (require less cpu time)
      ENDIF
      zsto = zdt
      clop = "ave("//TRIM(clop)//")"
      zout = nn_writedia * zdt

      ! Define indices of the horizontal output zoom and vertical limit storage
      iimi = 1      ;      iima = jpi
      ijmi = 1      ;      ijma = jpj
      ipk = jpk

      ! define time axis
      itmod = kt - nittrc000 + 1
      it    = kt
      iiter = ( nittrc000 - 1 ) / nn_dttrc

      ! 1. Define NETCDF files and fields at beginning of first time step
      ! -----------------------------------------------------------------

      IF( ll_print ) WRITE(numout,*) 'trcdii_wr kt=', kt

      IF( kt == nittrc000 ) THEN

         ! Define the NETCDF files for additional arrays : 2D or 3D

         ! Define the T grid file for tracer auxiliary files

         CALL dia_nam( clhstnam, nn_writedia, 'diad_T' )
         IF(lwp) WRITE(numout,*) " Name of NETCDF file ", clhstnam

         ! Define a netcdf FILE for 2d and 3d arrays

         CALL histbeg( clhstnam, jpi, glamt, jpj, gphit,             &
            &          iimi, iima-iimi+1, ijmi, ijma-ijmi+1,         &
            &          iiter, zjulian, zdt, nhoritd, nitd , domain_id=nidom, snc4chunks=snc4set )

         ! Vertical grid for 2d and 3d arrays

         CALL histvert( nitd, 'deptht', 'Vertical T levels','m', ipk, gdept_1d, ndepitd)

         ! Declare all the output fields as NETCDF variables

         ! more 3D horizontal arrays
         DO jl = 1, jpdia3d
            cltra  = TRIM( ctrc3d(jl) )   ! short title for 3D diagnostic
            cltral = TRIM( ctrc3l(jl) )  ! long title for 3D diagnostic
            cltrau = TRIM( ctrc3u(jl) )  ! UNIT for 3D diagnostic
            CALL histdef( nitd, cltra, cltral, cltrau, jpi, jpj, nhoritd,   &
               &          ipk, 1, ipk,  ndepitd, 32, clop, zsto, zout )
         END DO

         ! more 2D horizontal arrays
         DO jl = 1, jpdia2d
            cltra  = TRIM( ctrc2d(jl) )   ! short title for 2D diagnostic
            cltral = TRIM( ctrc2l(jl) )  ! long title for 2D diagnostic
            cltrau = TRIM( ctrc2u(jl) )  ! UNIT for 2D diagnostic
            CALL histdef( nitd, cltra, cltral, cltrau, jpi, jpj, nhoritd,  &
               &          1, 1, 1,  -99, 32, clop, zsto, zout )
         END DO

         ! TODO: more 2D vertical sections arrays : I or J indice fixed

         ! CLOSE netcdf Files
         CALL histend( nitd, snc4set )

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'End of NetCDF Initialization in trcdii_wr'
         IF( ll_print )   CALL FLUSH(numout )
         !
      ENDIF

      ! 2. Start writing data
      ! ---------------------

      IF( lwp .AND. MOD( itmod, nn_writedia ) == 0 ) THEN
         WRITE(numout,*) 'trcdii_wr : write NetCDF additional arrays at ', kt, 'time-step'
         WRITE(numout,*) '~~~~~~ '
      ENDIF

      ! more 3D horizontal arrays
      DO jl = 1, jpdia3d
         cltra  = TRIM( ctrc3d(jl) )   ! short title for 3D diagnostic
         CALL histwrite( nitd, cltra, it, trc3d(:,:,:,jl), ndimt50 ,ndext50)
      END DO

      ! more 2D horizontal arrays
      DO jl = 1, jpdia2d
         cltra  = TRIM( ctrc2d(jl) )   ! short title for 2D diagnostic
         CALL histwrite(nitd, cltra, it, trc2d(:,:,jl), ndimt51  ,ndext51)
      END DO

      ! Closing all files
      ! -----------------
      IF( kt == nitend )   CALL histclo(nitd)
      !

   END SUBROUTINE trcdii_wr

   SUBROUTINE trcdib_wr( kt )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE trcdib_wr  ***
      !!
      !! ** Purpose :   output of passive tracer : biological fields
      !!
      !! ** Method  :   At the beginning of the first time step (nittrc000), define all
      !!             the NETCDF files and fields for concentration of passive tracer
      !!
      !!        At each time step call histdef to compute the mean if necessary
      !!        Each nn_writebio time step, output the instantaneous or mean fields
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt          ! ocean time-step
      !!
      LOGICAL ::   ll_print = .FALSE.
      CHARACTER (len=40) ::   clhstnam, clop
      CHARACTER (len=20) ::   cltra, cltrau
      CHARACTER (len=80) ::   cltral
      INTEGER  ::   ji, jj, jk, jl
      INTEGER  ::   iimi, iima, ijmi, ijma, ipk, it, itmod, iiter
      REAL(wp) ::   zsto, zout, zdt
      !!----------------------------------------------------------------------

      ! Initialisation
      ! --------------
      
      ! local variable for debugging
      ll_print = .FALSE.
      ll_print = ll_print .AND. lwp

      ! Define frequency of output and means
      zdt = rdt
      IF( ln_mskland )   THEN   ;   clop = "only(x)"   ! put 1.e+20 on land (very expensive!!)
      ELSE                      ;   clop = "x"         ! no use of the mask value (require less cpu time)
      ENDIF
      zsto = zdt
      clop = "ave("//TRIM(clop)//")"
      zout = nn_writebio * zdt

      ! Define indices of the horizontal output zoom and vertical limit storage
      iimi = 1      ;      iima = jpi
      ijmi = 1      ;      ijma = jpj
      ipk = jpk

      ! define time axis
      itmod = kt - nittrc000 + 1
      it    = kt
      iiter = ( nittrc000 - 1 ) / nn_dttrc

      ! Define NETCDF files and fields at beginning of first time step
      ! --------------------------------------------------------------

      IF(ll_print) WRITE(numout,*)'trcdib_wr kt=',kt

      IF( kt == nittrc000 ) THEN

         ! Define the NETCDF files for biological trends

         CALL dia_nam(clhstnam,nn_writebio,'biolog')
         IF(lwp)WRITE(numout,*) " Name of NETCDF file for biological trends ", clhstnam
         ! Horizontal grid : glamt and gphit
         CALL histbeg( clhstnam, jpi, glamt, jpj, gphit,      &
            &    iimi, iima-iimi+1, ijmi, ijma-ijmi+1,          &
            &    iiter, zjulian, zdt, nhoritb, nitb , domain_id=nidom, snc4chunks=snc4set )
         ! Vertical grid for biological trends
         CALL histvert(nitb, 'deptht', 'Vertical T levels', 'm', ipk, gdept_1d, ndepitb)

         ! Declare all the output fields as NETCDF variables
         ! biological trends
         DO jl = 1, jpdiabio
            cltra  = TRIM( ctrbio(jl) )   ! short title for biological diagnostic
            cltral = TRIM( ctrbil(jl) )  ! long title for biological diagnostic
            cltrau = TRIM( ctrbiu(jl) )  ! UNIT for biological diagnostic
            CALL histdef( nitb, cltra, cltral, cltrau, jpi, jpj, nhoritb,  &
               &         ipk, 1, ipk,  ndepitb, 32, clop, zsto, zout)
         END DO

         ! CLOSE netcdf Files
          CALL histend( nitb, snc4set )

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'End of NetCDF Initialization in trcdib_wr'
         IF(ll_print) CALL FLUSH(numout )
         !
      ENDIF

      ! Start writing data
      ! ------------------

      ! biological trends
      IF( lwp .AND. MOD( itmod, nn_writebio ) == 0 ) THEN
         WRITE(numout,*) 'trcdit_wr : write NetCDF biological trends at ', kt, 'time-step'
         WRITE(numout,*) '~~~~~~ '
      ENDIF

      DO jl = 1, jpdiabio
         cltra  = TRIM( ctrbio(jl) )   ! short title for biological diagnostic
         CALL histwrite(nitb, cltra, it, trbio(:,:,:,jl), ndimt50,ndext50)
      END DO

      ! Closing all files
      ! -----------------
      IF( kt == nitend )   CALL histclo( nitb )
      !
   END SUBROUTINE trcdib_wr


   !!======================================================================
END MODULE trcdia
