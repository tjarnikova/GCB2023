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

MODULE par_planktom
!!======================================================================
!!
!!                         PARAMETER SMS
!!                       *******************************
!!
!!  purpose :
!!  ---------
!!     PARAMETER FILE for biogeochemical model
!!
!!======================================================================
!!  TOP 1.0,  LOCEAN-IPSL (2005)
!! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
!!----------------------------------------------------------------------
   IMPLICIT NONE
   PUBLIC
!!
!!    ASSIGN AN INTEGER TO NAME INDIVIDUAL TRACERS.
!!    FOR THE DGOM, THE TRACERS SHOULD BE IN THE ORDER:
!!    chemical tracers, detritus, zooplankton,
!!    nutrients, phytoplankton (C,Fe,Chl,Si), isotopes
!!    For documentation of variables see:
!!    http://lgmacweb.env.uea.ac.uk/green_ocean/model/code_description/var_dgom.html
!!----------------------------------------------------------------------
      LOGICAL, PUBLIC, PARAMETER ::   lk_planktom     = .TRUE.  !: PlankTOM flag 
      LOGICAL, PUBLIC, PARAMETER ::   lk_kriest       = .FALSE. !: Kriest flag 
      INTEGER, PUBLIC, PARAMETER :: jpzft=5
      INTEGER, PUBLIC, PARAMETER :: jppft=6,jpfmx=jpzft*(3+jpzft+jppft)
      INTEGER, PUBLIC, PARAMETER :: jptal=1,jpoxy=2,jpdic=3
      INTEGER, PUBLIC, PARAMETER :: jppiic=4,jpdin=5
      INTEGER, PUBLIC, PARAMETER :: jpsil=jpdin+1,jppo4=jpdin+2
      INTEGER, PUBLIC, PARAMETER :: jpfer=jpdin+3,jpdoc=jpdin+4,jpdsi=jpdin+5,jpcal=jpdin+6
      INTEGER, PUBLIC, PARAMETER :: jpara=jpdin+7,jpgon=jpdin+8,jpsfe=jpdin+9
      INTEGER, PUBLIC, PARAMETER :: jpbfe=jpsfe+1,jppoc=jpsfe+2,jpgoc=jpsfe+3,jpbac=jpsfe+4
      INTEGER, PUBLIC, PARAMETER :: jpmic=jpsfe+5
      INTEGER, PUBLIC, PARAMETER :: jppte=jpsfe+6,jpmes=jpsfe+7,jpgel=jpsfe+8,jpmac=jpsfe+9
      INTEGER, PUBLIC, PARAMETER :: jpdia=jpmac+1 ,jpmix=jpmac+2 ,jpcoc=jpmac+3 ,jppic=jpmac+4 ,jppha=jpmac+5 ,jpfix=jpmac+6
      INTEGER, PUBLIC, PARAMETER :: jpdfe=jpmac+7 ,jpnfe=jpmac+8 ,jpcfe=jpmac+9 ,jppfe=jpmac+10,jphfe=jpmac+11,jpffe=jpmac+12
      INTEGER, PUBLIC, PARAMETER :: jpdch=jpmac+13,jpnch=jpmac+14,jpcch=jpmac+15,jppch=jpmac+16,jphch=jpmac+17,jpfch=jpmac+18
      INTEGER, PUBLIC, PARAMETER :: jpbsi=jpmac+19
      INTEGER, PUBLIC, PARAMETER :: jpc11=jpmac+19
      INTEGER, PUBLIC, PARAMETER :: jpdmd=jpc11
      INTEGER, PUBLIC, PARAMETER :: jpn2s=jpdmd
      INTEGER, PUBLIC, PARAMETER :: jp_planktom=jpn2s
END MODULE par_planktom
