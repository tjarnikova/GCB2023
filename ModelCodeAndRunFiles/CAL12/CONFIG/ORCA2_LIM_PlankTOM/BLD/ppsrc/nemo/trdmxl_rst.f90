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

MODULE trdmxl_rst
   !!=================================================================================
   !!                       ***  MODULE  trdmxl_rst  ***
   !! Ocean dynamic :  Input/Output files for restart on mixed-layer diagnostics
   !!=================================================================================
   !! History :  1.0  ! 2005-05 (C. Deltel)  Original code
   !!---------------------------------------------------------------------------------

   !!---------------------------------------------------------------------------------
   !!  trd_mxl_rst_write : write mixed layer trend restart
   !!  trd_mxl_rst_read  : read  mixed layer trend restart
   !!---------------------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE trd_oce         ! trends: ocean variables
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O module
   USE restart         ! only for lrst_oce

   IMPLICIT NONE
   PRIVATE
  
   PUBLIC   trd_mxl_rst_read    ! routine called by trd_mxl_init
   PUBLIC   trd_mxl_rst_write   ! routine called by step.F90
  
   INTEGER ::   nummxlw         ! logical unit for mxl restart

   !!---------------------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: trdmxl_rst.F90 5341 2015-06-03 14:59:46Z davestorkey $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!---------------------------------------------------------------------------------
CONTAINS
  
   SUBROUTINE trd_mxl_rst_write( kt )
      !!--------------------------------------------------------------------------------
      !!                  ***  SUBROUTINE trd_mxl_rst_wri  ***
      !!                
      !! ** Purpose :   Write mixed-layer diagnostics restart fields.
      !!--------------------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !
      CHARACTER (len=35) :: charout
      INTEGER ::   jk                 ! loop indice
      CHARACTER(LEN=20)   ::   clkt     ! ocean time-step deine as a character
      CHARACTER(LEN=50)   ::   clname   ! output restart file name
      CHARACTER(LEN=256)  ::   clpath   ! full path to restart file
      !!--------------------------------------------------------------------------------

      ! to get better performances with NetCDF format:
      ! we open and define the ocean restart_mxl file one time step before writing the data (-> at nitrst - 1)
      ! except if we write ocean restart_mxl files every time step or if an ocean restart_mxl file was writen at nitend - 1
      IF( kt == nitrst - 1 .OR. nstock == 1 .OR. ( kt == nitend .AND. MOD( nitend - 1, nstock ) == 0 ) ) THEN
         ! beware of the format used to write kt (default is i8.8, that should be large enough...)
         IF( nitrst > 999999999 ) THEN   ;   WRITE(clkt, *       ) nitrst
         ELSE                            ;   WRITE(clkt, '(i8.8)') nitrst
         ENDIF
         ! create the file
         clname = TRIM(cexper)//"_"//TRIM(ADJUSTL(clkt))//"_"//TRIM(cn_trdrst_out)
         clpath = TRIM(cn_ocerst_outdir)
         IF( clpath(LEN_TRIM(clpath):) /= '/' ) clpath = TRIM(clpath) // '/'
         IF(lwp) THEN
            WRITE(numout,*)
            SELECT CASE ( jprstlib )
            CASE ( jprstdimg )   ;   WRITE(numout,*) '             open ocean restart_mxl binary file: '//clname
            CASE DEFAULT         ;   WRITE(numout,*) '             open ocean restart_mxl NetCDF file: '//clname
            END SELECT
            IF( kt == nitrst - 1 ) THEN   ;   WRITE(numout,*) '             kt = nitrst - 1 = ', kt,' date= ', ndastp
            ELSE                          ;   WRITE(numout,*) '             kt = '             , kt,' date= ', ndastp
            ENDIF
         ENDIF

         CALL iom_open( TRIM(clpath)//TRIM(clname), nummxlw, ldwrt = .TRUE., kiolib = jprstlib )
      ENDIF

      IF( kt == nitrst .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trdmxl_rst: output for ML diags. restart, with trd_mxl_rst_write routine kt =', kt
         WRITE(numout,*) '~~~~~~~~~~'
         WRITE(numout,*)
      ENDIF

      IF( ln_trdmxl_instant ) THEN 
         !-- Temperature
         CALL iom_rstput( kt, nitrst, nummxlw, 'tmlbb'           , tmlbb           )
         CALL iom_rstput( kt, nitrst, nummxlw, 'tmlbn'           , tmlbn           )
         CALL iom_rstput( kt, nitrst, nummxlw, 'tmlatfb'         , tmlatfb         )

         !-- Salinity
         CALL iom_rstput( kt, nitrst, nummxlw, 'smlbb'           , smlbb           )
         CALL iom_rstput( kt, nitrst, nummxlw, 'smlbn'           , smlbn           )
         CALL iom_rstput( kt, nitrst, nummxlw, 'smlatfb'         , smlatfb         )
      ELSE
         CALL iom_rstput( kt, nitrst, nummxlw, 'hmxlbn'          , hmxlbn          )

         !-- Temperature
         CALL iom_rstput( kt, nitrst, nummxlw, 'tmlbn'           , tmlbn           )
         CALL iom_rstput( kt, nitrst, nummxlw, 'tml_sumb'        , tml_sumb        )
         DO jk = 1, jpltrd
            IF( jk < 10 ) THEN   ;   WRITE(charout,FMT="('tmltrd_csum_ub_', I1)")   jk
            ELSE                 ;   WRITE(charout,FMT="('tmltrd_csum_ub_', I2)")   jk
            ENDIF
            CALL iom_rstput( kt, nitrst, nummxlw, charout,  tmltrd_csum_ub(:,:,jk) )
         ENDDO
         CALL iom_rstput( kt, nitrst, nummxlw, 'tmltrd_atf_sumb' , tmltrd_atf_sumb )

         !-- Salinity
         CALL iom_rstput( kt, nitrst, nummxlw, 'smlbn'           , smlbn           )
         CALL iom_rstput( kt, nitrst, nummxlw, 'sml_sumb'        , sml_sumb        )
         DO jk = 1, jpltrd
            IF( jk < 10 ) THEN   ;   WRITE(charout,FMT="('smltrd_csum_ub_', I1)")   jk
            ELSE                 ;   WRITE(charout,FMT="('smltrd_csum_ub_', I2)")   jk
            ENDIF
            CALL iom_rstput( kt, nitrst, nummxlw, charout , smltrd_csum_ub(:,:,jk) )
         ENDDO
         CALL iom_rstput( kt, nitrst, nummxlw, 'smltrd_atf_sumb' , smltrd_atf_sumb )
      ENDIF
      !
      IF( kt == nitrst ) THEN
         CALL iom_close( nummxlw )     ! close the restart file (only at last time step)
         lrst_oce = .FALSE.
      ENDIF
      ! 
   END SUBROUTINE trd_mxl_rst_write


   SUBROUTINE trd_mxl_rst_read
      !!----------------------------------------------------------------------------
      !!                   ***  SUBROUTINE trd_mxl_rst_lec  ***
      !!                   
      !! ** Purpose :   Read file for mixed-layer diagnostics restart.
      !!----------------------------------------------------------------------------
      INTEGER  ::  inum       ! temporary logical unit
      !
      CHARACTER (len=35) :: charout
      INTEGER ::   jk         ! loop indice
      INTEGER ::   jlibalt = jprstlib
      LOGICAL ::   llok
      CHARACTER(LEN=256)  ::   clpath   ! full path to restart file
      !!-----------------------------------------------------------------------------

      IF(lwp)  THEN
         WRITE(numout,*)
         WRITE(numout,*) ' trd_mxl_rst_read : read the NetCDF mixed layer trend restart file'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~'
      ENDIF

      clpath = TRIM(cn_ocerst_indir)
      IF( clpath(LEN_TRIM(clpath):) /= '/' ) clpath = TRIM(clpath) // '/'

      IF ( jprstlib == jprstdimg ) THEN
         ! eventually read netcdf file (monobloc)  for restarting on different number of processors
         ! if {cn_trdrst_in}.nc exists, then set jlibalt to jpnf90
         INQUIRE( FILE = TRIM(clpath)//TRIM(cn_trdrst_in)//'.nc', EXIST = llok )
         IF ( llok ) THEN   ;   jlibalt = jpnf90   
         ELSE               ;   jlibalt = jprstlib   
         ENDIF
      ENDIF

      CALL iom_open( TRIM(clpath)//TRIM(cn_trdrst_in), inum, kiolib = jlibalt ) 

      IF( ln_trdmxl_instant ) THEN 
         !-- Temperature
         CALL iom_get( inum, jpdom_autoglo, 'tmlbb'           , tmlbb          )
         CALL iom_get( inum, jpdom_autoglo, 'tmlbn'           , tmlbn          )
         CALL iom_get( inum, jpdom_autoglo, 'tmlatfb'         , tmlatfb        )
         !
         !-- Salinity
         CALL iom_get( inum, jpdom_autoglo, 'smlbb'           , smlbb          )
         CALL iom_get( inum, jpdom_autoglo, 'smlbn'           , smlbn          )
         CALL iom_get( inum, jpdom_autoglo, 'smlatfb'         , smlatfb        )
      ELSE
         CALL iom_get( inum, jpdom_autoglo, 'hmxlbn'          , hmxlbn         ) ! needed for hmxl_sum
         !
         !-- Temperature
         CALL iom_get( inum, jpdom_autoglo, 'tmlbn'           , tmlbn          ) ! needed for tml_sum
         CALL iom_get( inum, jpdom_autoglo, 'tml_sumb'        , tml_sumb       )
         DO jk = 1, jpltrd
            IF( jk < 10 ) THEN   ;   WRITE(charout,FMT="('tmltrd_csum_ub_', I1)")   jk
            ELSE                 ;   WRITE(charout,FMT="('tmltrd_csum_ub_', I2)")   jk
            ENDIF
            CALL iom_get( inum, jpdom_autoglo, charout, tmltrd_csum_ub(:,:,jk) )
         END DO
         CALL iom_get( inum, jpdom_autoglo, 'tmltrd_atf_sumb' , tmltrd_atf_sumb)
         !
         !-- Salinity
         CALL iom_get( inum, jpdom_autoglo, 'smlbn'           , smlbn          ) ! needed for sml_sum
         CALL iom_get( inum, jpdom_autoglo, 'sml_sumb'        , sml_sumb       )
         DO jk = 1, jpltrd
            IF( jk < 10 ) THEN   ;   WRITE(charout,FMT="('smltrd_csum_ub_', I1)")   jk
            ELSE                 ;   WRITE(charout,FMT="('smltrd_csum_ub_', I2)")   jk
            ENDIF
            CALL iom_get( inum, jpdom_autoglo, charout, smltrd_csum_ub(:,:,jk) )
         END DO
         CALL iom_get( inum, jpdom_autoglo, 'smltrd_atf_sumb' , smltrd_atf_sumb)
         !
         CALL iom_close( inum )
      ENDIF
      !
   END SUBROUTINE trd_mxl_rst_read
  
  !!=================================================================================
END MODULE trdmxl_rst
