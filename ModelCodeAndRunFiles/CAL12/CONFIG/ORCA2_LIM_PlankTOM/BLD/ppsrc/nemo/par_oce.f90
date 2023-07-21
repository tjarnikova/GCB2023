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

MODULE par_oce
   !!======================================================================
   !!                        ***  par_oce  ***
   !! Ocean :   set the ocean parameters
   !!======================================================================
   !! History :  OPA  !  1991     (Imbard, Levy, Madec)  Original code
   !!   NEMO     1.0  !  2004-01  (G. Madec, J.-M. Molines)  Free form and module
   !!            3.3  !  2010-09  (C. Ethe) TRA-TRC merge: add jpts, jp_tem & jp_sal
   !!----------------------------------------------------------------------
   USE par_kind          ! kind parameters

   IMPLICIT NONE
   PUBLIC

   !!----------------------------------------------------------------------
   !!   Domain decomposition
   !!----------------------------------------------------------------------
   !! if we dont use massively parallel computer (parameters jpni=jpnj=1) so jpiglo=jpi and jpjglo=jpj
   INTEGER, PUBLIC            ::   jpni         !: number of processors following i 
   INTEGER, PUBLIC            ::   jpnj         !: number of processors following j
   INTEGER, PUBLIC            ::   jpnij        !: nb of local domain = nb of processors ( <= jpni x jpnj )
   INTEGER, PUBLIC, PARAMETER ::   jpr2di = 0   !: number of columns for extra outer halo 
   INTEGER, PUBLIC, PARAMETER ::   jpr2dj = 0   !: number of rows    for extra outer halo 
   INTEGER, PUBLIC, PARAMETER ::   jpreci = 1   !: number of columns for overlap 
   INTEGER, PUBLIC, PARAMETER ::   jprecj = 1   !: number of rows    for overlap 

   !!----------------------------------------------------------------------
   !!                   namcfg namelist parameters
   !!----------------------------------------------------------------------
   CHARACTER(lc) ::   cp_cfg           !: name of the configuration
   CHARACTER(lc) ::   cp_cfz           !: name of the zoom of configuration
   INTEGER       ::   jp_cfg           !: resolution of the configuration

   ! data size                                       !!! * size of all input files *
   INTEGER       ::   jpidta           !: 1st lateral dimension ( >= jpi )
   INTEGER       ::   jpjdta           !: 2nd    "         "    ( >= jpj )
   INTEGER       ::   jpkdta           !: number of levels      ( >= jpk )

   ! global or zoom domain size                      !!! * computational domain *
   INTEGER       ::   jpiglo           !: 1st dimension of global domain --> i
   INTEGER       ::   jpjglo           !: 2nd    -                  -    --> j

   ! zoom starting position 
   INTEGER       ::   jpizoom          !: left bottom (i,j) indices of the zoom
   INTEGER       ::   jpjzoom          !: in data domain indices

   ! Domain characteristics
   INTEGER       ::   jperio           !: lateral cond. type (between 0 and 6)
   !                                       !  = 0 closed                 ;   = 1 cyclic East-West
   !                                       !  = 2 equatorial symmetric   ;   = 3 North fold T-point pivot
   !                                       !  = 4 cyclic East-West AND North fold T-point pivot
   !                                       !  = 5 North fold F-point pivot
   !                                       !  = 6 cyclic East-West AND North fold F-point pivot

   ! Input file read offset
   LOGICAL       ::   ln_use_jattr     !: Use file global attribute: open_ocean_jstart to determine start j-row 
                                           ! when reading input from those netcdf files that have the 
                                           ! attribute defined. This is designed to enable input files associated 
                                           ! with the extended grids used in the under ice shelf configurations to 
                                           ! be used without redundant rows when the ice shelves are not in use.

   !!  Values set to pp_not_used indicates that this parameter is not used in THIS config.
   !!  Values set to pp_to_be_computed  indicates that variables will be computed in domzgr
   REAL(wp)      ::   pp_not_used       = 999999._wp   !: vertical grid parameter
   REAL(wp)      ::   pp_to_be_computed = 999999._wp   !:    -      -       -




   !!---------------------------------------------------------------------
   !! Active tracer parameters
   !!---------------------------------------------------------------------
   INTEGER, PUBLIC, PARAMETER ::   jpts   = 2    !: Number of active tracers (=2, i.e. T & S )
   INTEGER, PUBLIC, PARAMETER ::   jp_tem = 1    !: indice for temperature
   INTEGER, PUBLIC, PARAMETER ::   jp_sal = 2    !: indice for salinity

   !!---------------------------------------------------------------------
   !! Domain Matrix size  (if AGRIF, they are not all parameters)
   !!---------------------------------------------------------------------
   INTEGER, PUBLIC  ::   jpi   ! = ( jpiglo-2*jpreci + (jpni-1) ) / jpni + 2*jpreci   !: first  dimension
   INTEGER, PUBLIC  ::   jpj   ! = ( jpjglo-2*jprecj + (jpnj-1) ) / jpnj + 2*jprecj   !: second dimension
   INTEGER, PUBLIC  ::   jpk   ! = jpkdta
   INTEGER, PUBLIC  ::   jpim1 ! = jpi-1                                            !: inner domain indices
   INTEGER, PUBLIC  ::   jpjm1 ! = jpj-1                                            !:   -     -      -
   INTEGER, PUBLIC  ::   jpkm1 ! = jpk-1                                            !:   -     -      -
   INTEGER, PUBLIC  ::   jpij  ! = jpi*jpj                                          !:  jpi x jpj

   !!---------------------------------------------------------------------
   !! Optimization/control flags
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_esopa     = .FALSE.  !: flag to activate the all options

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: par_oce.F90 5118 2015-03-03 13:36:04Z acc $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE par_oce
