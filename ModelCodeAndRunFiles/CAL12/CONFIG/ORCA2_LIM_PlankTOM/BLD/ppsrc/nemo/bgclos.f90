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

      SUBROUTINE bgclos
!CC   ------------------------------------------------------------------
!CC   
!CC   ROUTINE bgclos: DGOM MODEL
!CC   ****************************
!CC   
!C
!C     PURPOSE.
!C     --------
!C                   Calculates loss terms of all plankton
!C     EXTERNALS.
!C     ----------
!C          NONE.
!C
!C   MODIFICATIONS:
!C   --------------
!C      original p4zloss : 2002    E. Buitenhuis
!C ---------------------------------------------------------------------
      USE trp_trc
      
      USE sms_planktom
      USE oce_trc
      IMPLICIT NONE
!
! local variables
!
      INTEGER ji, jj, jk, jl, jm, jn
      REAL cmpfoo(jppoc:jpdia+jppft-1),compabio
      REAL stosil
      REAL graze(jpzft),zdenom(jpzft),zprozoo(jpzft)
      REAL zmokmes,zmokmac,zmokgel
!
! local dgom variables
!
      REAL zlosphy
      REAL zgrdmp,zstre2,zstre3
      grazoc = 0.
      grazof = 0.
! 
! -------------------------------------------------------------------
! Big loop to compute the various sink and source terms 
! for phytoplankton, zooplankton and organic matter reservoirs. 
! -------------------------------------------------------------------
!     
! no grazing if concentration goes below rn_gramin (Strom et al., MEPS 2000)
!
      DO jk = 1, jpk-1
        DO jj = 2, nlcj-1
          DO ji = 2, nlci-1
!     
! Protect plankton from extinction
!     
            DO jl = jppoc, jpdia+jppft-1
              cmpfoo(jl) = max((trn(ji,jj,jk,jl)-rn_gramin),0.)
            END DO
!
! sum all the phyto + zoo biomass 
!
            compabio=0. 
            DO jl = jpmic, jpdia+jppft-1
              compabio = compabio + trn(ji,jj,jk,jl)
            ENDDO
!
! k-half for top mortality
!
!           zmokmes = 20.e-6
            zmokmac = 20.e-6
            zmokgel = 20.e-6
!
! -------------------------------------------------------------------
! Respiration rates of phytoplankton as a fraction of growth
! -------------------------------------------------------------------
            stosil=trn(ji,jj,jk,jpbsi)/(trn(ji,jj,jk,jpdia)+rtrn)
            DO jl = jppoc, jpgoc
              stofoo(ji,jj,jk,jl,2) = trn(ji,jj,jk,jl-2)          &
     &          /max(trn(ji,jj,jk,jl),rtrn)
            END DO
            DO jl = jpdia, jpdia+jppft-1
              zlosphy = rn_resphy(jl)*rn_mumpft(jl)*tgfunc(ji,jj,jk,jl)
              resphy(ji,jj,jk,jl,1) = zlosphy*cmpfoo(jl) &
     &          /rjjss*rfact*tmask(ji,jj,jk)
              resphy(ji,jj,jk,jl,2) = resphy(ji,jj,jk,jl,1)      &
     &           *stofoo(ji,jj,jk,jl,2)
              resphy(ji,jj,jk,jl,3) = resphy(ji,jj,jk,jl,1)      &
     &           *stofoo(ji,jj,jk,jl,3)
            END DO
            DO jm = 1, jpzft
! -------------------------------------------------------------------
! Respiration rates of zooplankton
! -------------------------------------------------------------------
              reszoo(ji,jj,jk,jm) = rn_reszoo(jm)*rn_retzoo(jm)**tsn(ji,jj,jk,1) &
     &          /rjjss*rfact*cmpfoo(jpbac+jm)*tmask(ji,jj,jk)
! -------------------------------------------------------------------
! Grazing rates of zooplankton as a function of temperature
! -------------------------------------------------------------------
              graze(jm) = rn_grazoo(jm)/rjjss*rfact*tmask(ji,jj,jk)     &
     &          *tgfunc(ji,jj,jk,jpbac+jm)
! compute denominator
              zdenom(jm) = rn_grkzoo(jm)
            END DO
            DO jn = 1, jpfoo
              jm = grizoo(1,jn)
              jl = grizoo(2,jn)
              zdenom(jm) = zdenom(jm)+rn_prfzoo(jm,jl)*trn(ji,jj,jk,jl)
            END DO
!
          losbsi(ji,jj,jk) = resphy(ji,jj,jk,jpdia,1)*stosil
!     
! -------------------------------------------------------------------
! Zooplankton mortality
! -------------------------------------------------------------------
!
! add mortality by top predators on macrozoo, using global biomass as a proxy 
!
          tormac(ji,jj,jk) =                &
     &       rn_mormac/rjjss*rfact*cmpfoo(jpmac)/(zmokmac+cmpfoo(jpmac))*compabio &
     &       *tmask(ji,jj,jk)*(rn_motmac**tsn(ji,jj,jk,1))
!
! macrozooplankton predation mortality turns to near zero when ice is present. 
!     this is a protection mechanisms specific to krill. 
!
          tormac(ji,jj,jk) = tormac(ji,jj,jk) *                           &
     &       max( (1.- fr_i(ji,jj)/(fr_i(ji,jj)+rtrn)) ,0.01)
! add mortality by top predators on gelatinouszoo, using global biomass as a proxy
          torgel(ji,jj,jk) =                &
!    &       rn_morgel/rjjss*rfact*cmpfoo(jpgel)*morfsh(ji,jj)     &
! changing proxy biomass mortality to fish mortality ^^
     &       rn_morgel/rjjss*rfact*cmpfoo(jpgel)/(zmokgel+cmpfoo(jpgel))*compabio &
     &       *tmask(ji,jj,jk)*(rn_motgel**tsn(ji,jj,jk,1))
!     
! -------------------------------------------------------------------
! Preference of zooplankton for Phyto and POC (Fasham et al. 1990)
! -------------------------------------------------------------------
            DO jm = 1, jpzft
              zprozoo(jm)=graze(jm)*trn(ji,jj,jk,jpbac+jm)/zdenom(jm)
            END DO
            DO jn = 1, jpfoo
              jm = grizoo(1,jn)
              jl = grizoo(2,jn)
              grazoo(ji,jj,jk,jm,jl) = zprozoo(jm)                      &
     &          *rn_prfzoo(jm,jl)*cmpfoo(jl)
              grazoc(ji,jj,jk,jm)=grazoc(ji,jj,jk,jm)+grazoo(ji,jj,jk,jm,jl)
              grazof(ji,jj,jk,jm)=grazof(ji,jj,jk,jm)+grazoo(ji,jj,jk,jm,jl)  &
     &          *stofoo(ji,jj,jk,jl,2)
            END DO
!
! zooplankton model growth efficiency
!      if meso food respiration is lower than basal respiration increase GE
!      if there's not enough iron in food, decrease GE
!
            DO jm = 1, jpzft
          mgezoo(ji,jj,jk,jm)=min(1.-rn_unazoo(jm),                     &
     &      rn_ggezoo(jm)+reszoo(ji,jj,jk,jm)/max(rtrn,grazoc(ji,jj,jk,jm)), &
     &      grazof(ji,jj,jk,jm)*(1.-rn_unazoo(jm))/                     &
     &      max(grazoc(ji,jj,jk,jm)*ferat3,minfer))
          grapoc(ji,jj,jk,jm)=grazoc(ji,jj,jk,jm)*rn_unazoo(jm)
          grarem(ji,jj,jk,jm)=grazoc(ji,jj,jk,jm)*                      &
     &      (1.-mgezoo(ji,jj,jk,jm)-rn_unazoo(jm))
          grafer(ji,jj,jk,jm)=grazof(ji,jj,jk,jm)                       &
     &      *(1.-rn_unazoo(jm))-ferat3*mgezoo(ji,jj,jk,jm)*grazoc(ji,jj,jk,jm)
! grazing on diatom Si
              losbsi(ji,jj,jk) = losbsi(ji,jj,jk)                       &
     &          +grazoo(ji,jj,jk,jm,jpdia)*stosil
            END DO
          END DO
        END DO
      END DO
      RETURN
      END
