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

MODULE trp_trc

   !!======================================================================
   !! Module trp_trc
   !!======================================================================
   !!  TOP 1.0,  LOCEAN-IPSL (2005)
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
   !! passive tracers number
   USE par_trc , ONLY : &
      jptra    =>   jptra       !!: number of passive tracers

   USE par_trc , ONLY : &
      jpdia2d  =>  jpdia2d , &  !!: number of 2d outputs
      jpdia3d  =>  jpdia3d

   !! passive tracers fields 
   USE trc , ONLY :  &
      trai   =>   trai , &  !!: initial total tracer
      trb    =>   trb  , &  !!: tracer field (before)
      tra    =>   tra  , &  !!: tracer field (now)
      trn    =>   trn       !!: tracer field (after)

   !! time step - not used in PlankTOM so not required
!!   USE trc , ONLY :  &
!!      nn_dttrc =>   nn_dttrc    !!: frequency of step on passive tracers (NAMELIST)

   !! non-centered advection scheme (smolarkiewicz)
   USE trc , ONLY : &
      rtrn   =>   rtrn      !!: value for truncation (NAMELIST)

   USE trc , ONLY : &
      ctrcnm   =>   ctrcnm      !!: value for truncation (NAMELIST)

END MODULE trp_trc
