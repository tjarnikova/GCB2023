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

MODULE domstp
   !!==============================================================================
   !!                       ***  MODULE domstp   ***
   !! Ocean initialization : time domain
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   dom_stp        : ocean time domain initialization
   !!----------------------------------------------------------------------
   !! History :  OPA  ! 1990-10  (O. Marti)  Original code
   !!                 ! 1996-01  (G. Madec)  terrain following coordinates
   !!   NEMO     1.0  ! 2002-08  (G. Madec)  F90: Free form and module
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dom_stp   ! routine called by inidom.F90

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
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: domstp.F90 4292 2013-11-20 16:28:04Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dom_stp
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dom_stp  ***
      !!          
      !! ** Purpose :   Intialize ocean time step for the run
      !!
      !! ** Method  : - Initialization of a coef. use in the Asselin time
      !!      filter:  atfp1 = 1 - 2 * atfp  where atfp is the Asselin time
      !!      filter parameter read in namelist
      !!              - Model time step:
      !!      nacc = 0 : synchronous time intergration. 
      !!      There is one time step only, defined by: rdt, rdttra(k)=rdt
      !!      nacc = 1 : accelerating the convergence. There is 2 different
      !!      time steps for dynamics and tracers:
      !!        rdt      : dynamical part
      !!        rdttra(k): temperature and salinity
      !!      The tracer time step is a function of vertical level. the model
      !!      reference time step ( i.e. for wind stress, surface heat and
      !!      salt fluxes) is the surface tracer time step is rdttra(1).
      !!         N.B. depth dependent acceleration of convergence is not im-
      !!      plemented for s-coordinate.
      !!
      !! ** Action  : - rdttra   : vertical profile of tracer time step
      !!              - atfp1    : = 1 - 2*atfp
      !!
      !! References :   Bryan, K., 1984, J. Phys. Oceanogr., 14, 666-673.
      !!----------------------------------------------------------------------
      INTEGER ::   jk              ! dummy loop indice
      !!----------------------------------------------------------------------

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dom_stp : time stepping setting'
         WRITE(numout,*) '~~~~~~~'
      ENDIF

      ! 0. Asselin Time filter
      ! ----------------------
      
      atfp1 = 1. - 2. * atfp

      SELECT CASE ( nacc )

         CASE ( 0 )                ! Synchronous time stepping
            IF(lwp) WRITE(numout,*)'               synchronous time stepping'
            IF(lwp) WRITE(numout,*)'               dynamics and tracer time step = ', rdt/3600., ' hours'

            rdttra(:) = rdt

         CASE ( 1 )                ! Accelerating the convergence
            IF(lwp) WRITE(numout,*) '              no tracer damping in the turbocline'
            IF(lwp) WRITE(numout,*)'               accelerating the convergence'
            IF(lwp) WRITE(numout,*)'               dynamics time step = ', rdt/3600., ' hours'
            IF( ln_sco .AND. rdtmin /= rdtmax .AND. lk_vvl )   &
                 & CALL ctl_stop ( ' depth dependent acceleration of convergence not implemented in s-coordinates &
                 &                   nor in variable volume' )
            IF(lwp) WRITE(numout,*)'         tracers   time step :  dt (hours)  level'

            DO jk = 1, jpk
               IF( gdept_1d(jk) <= rdth ) rdttra(jk) = rdtmin
               IF( gdept_1d(jk) >  rdth ) THEN
                  rdttra(jk) = rdtmin + ( rdtmax - rdtmin )   &
                                      * ( EXP( ( gdept_1d(jk ) - rdth ) / rdth ) - 1. )   &
                                      / ( EXP( ( gdept_1d(jpk) - rdth ) / rdth ) - 1. )
               ENDIF
               IF(lwp) WRITE(numout,"(36x,f5.2,5x,i3)") rdttra(jk)/3600., jk
            END DO  

         CASE DEFAULT              ! E R R O R 

            WRITE(ctmp1,*) ' nacc value e r r o r, nacc= ',nacc
            CALL ctl_stop( ctmp1 )

      END SELECT

   END SUBROUTINE dom_stp

   !!======================================================================
END MODULE domstp
