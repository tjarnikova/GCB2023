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

      SUBROUTINE bgcbio(kt)
!CC   ------------------------------------------------------------------
!CC   
!CC   ROUTINE bgcbio : DGOM MODEL
!CC   *****************************
!CC   
!C
!C     PURPOSE.
!C     --------
!C          ECOSYSTEM MODEL
!C
!C     EXTERNALS.
!C     ----------
!C          NONE.
!C
!C   MODIFICATIONS:
!C   --------------
!C      original      : 2001    O. Aumont
!C      modifications : 2002    E. Buitenhuis
!C ---------------------------------------------------------------------
      USE trc
      USE trp_trc
      USE sms_planktom
      USE oce_trc
      USE lbclnk
      USE lib_mpp
      USE iom
      USE traqsr
      USE trcsms_c14b , ONLY : rn_bfrd14
      IMPLICIT NONE
! local declarations
      INTEGER ji, jj, jk, ikt, jl, jm, jn, iostat, kt,level
      REAL consum,consump,gongoc,dgongoc,mfecal,dumch4,temch4
      REAL omeara,omecal,remara,remco3
! For output 
      REAL zgrdmp
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
!
! Initialisation of variables
! ------------------------------------------------------------------
!
      grazing=0.
      ppint  =0.
      ppt=0.
      pptdoc =0.
      tchl=0.
      trophic=0.
      out2d=0.
      out3d=0.
      DO jk = 1, jpk-1
       DO jj = 2, nlcj-1
        DO ji = 2, nlci-1
!
! Iron and Si deposition at the surface
! the dust variable is in kgdust/m2/s (from Jickells et al. 2005). 
! We use 0.035gFe/gdust, 0.308 gSi/gdust (or 8.8gSi/gFe). 
! The solubility of Fe in dust is rn_fersol (usually 2%, set in namelist). 
! The solubility of Si in dust is 7.5%. 
! Variables irondep and sidep are in mol/L/time_step for Fe and Si, respectivel
! kg/m3=g/L
! The molecular weight of Fe and Si is 55.85 and 28.01, respectively.
! 0.3314 is the correction so that the total iron added is the same as when
! iron was only added at the surface (0.80037 kmol/s), and is approximately
! the integral of
! (gdept(jk)/gdept(1))**(-0.858) over 30 layers
!
            irondep(ji,jj,jk) = dust(ji,jj)*rn_fersol*rfact             &
     &        *(gdept_0(ji,jj,jk)/gdept_0(ji,jj,1))**(-0.858)*0.3314    &
     &        /(e3t_0(ji,jj,jk)*55.85) 
!            if ( rn_fersol .gt. .999 ) then
!                    sidep(ji,jj,1) = irondep(ji,jj,1)*             &
!     &                   (0.308/0.035)*(0.075/0.0125)*(55.85/28.01)
!            else 
            sidep(ji,jj,jk)=irondep(ji,jj,jk)*                     &
     &                (0.308/0.035)*(rn_silsol/rn_fersol)* (55.85/28.01)
!            endif
        END DO
       END DO
      END DO
!
! Call optical /production routine to compute phytoplankton growth rate
! over the global ocean. Growth rates for each element are 
! computed (C, Si, Fe, Chl)
! -------------------------------------------------------------------
!
      IF ( kt .EQ. nit000 ) CALL tra_qsr(kt)
      CALL bgcpro
!
! Calculate loss rates
! -------------------------------------------------------------------
!     
      CALL bgclos
!
! Call subroutine for computation of the vertical flux 
! of particulate organic matter
! -------------------------------------------------------------------
!
      CALL bgcsnk
!
! pre-compute concentration changes of the (rapidly
! varying) tracers for preventing them to fall below 0
! -------------------------------------------------------------------
!
      CALL bgcnul
!
! Recompute the SMS related to zooplankton grazing
! -------------------------------------------------------------------
!
      grazoc = 0.
      grazof = 0.
      DO jk = 1, jpk-1
        DO jj = 2, nlcj-1
          DO ji = 2, nlci-1
!
! assign grazing to organic matter 
! 
            DO jn = 1, jpfoo
              jm = grizoo(1,jn)
              jl = grizoo(2,jn)
              grazoc(ji,jj,jk,jm) = grazoc(ji,jj,jk,jm)+grazoo(ji,jj,jk,jm,jl)
              grazof(ji,jj,jk,jm) = grazof(ji,jj,jk,jm)+grazoo(ji,jj,jk,jm,jl)*stofoo(ji,jj,jk,jl,2)
              level=min(max(int(float(jl-jpbac)/float(jpzft)+0.81),0),2)+1 !1=detriti,2=carni,3=herbivory
              trophic(ji,jj,jk,level) = trophic(ji,jj,jk,level)+grazoo(ji,jj,jk,jm,jl)
            END DO
!
! calculate GGE zooplankton
!
            DO jm = 1, jpzft
              mgezoo(ji,jj,jk,jm) = min(1.-rn_unazoo(jm),rn_ggezoo(jm)  &
     &          +reszoo(ji,jj,jk,jm)/max(rtrn,grazoc(ji,jj,jk,jm)),     &
     &          grazof(ji,jj,jk,jm)*(1.-rn_unazoo(jm))                  &
     &          /max(grazoc(ji,jj,jk,jm)*ferat3,minfer))
              grafer(ji,jj,jk,jm)=grazof(ji,jj,jk,jm)*(1.-rn_unazoo(jm))&
     &          -ferat3*mgezoo(ji,jj,jk,jm)*grazoc(ji,jj,jk,jm)
              grarem(ji,jj,jk,jm) = grazoc(ji,jj,jk,jm)                 &
     &          *(1.-mgezoo(ji,jj,jk,jm)-rn_unazoo(jm))
              grapoc(ji,jj,jk,jm) = grazoc(ji,jj,jk,jm)*rn_unazoo(jm)
            END DO
!
            bactge(ji,jj,jk) = min(bactge(ji,jj,jk),                   &
     &        (ubafer(ji,jj,jk)+remsfe(ji,jj,jk)+rembfe(ji,jj,jk))        &
     &        /max((remdoc(ji,jj,jk)+rempoc(ji,jj,jk)+remgoc(ji,jj,jk))    &
     &        *ferat3,minfer))
!     
! Determination of tracer concentrations as a function of 
! biological sources and sinks
! -------------------------------------------------------------------
!
! Evolution of O2
! ---------------
!
            consum = (1.-bactge(ji,jj,jk))*(remdoc(ji,jj,jk)            &
     &        +rempoc(ji,jj,jk))+resbac(ji,jj,jk)/3.0
            DO jl = jpdia, jpdia+jppft-1
              consum = consum +resphy(ji,jj,jk,jl,1)
            END DO
            DO jm = 1, jpzft
              consum = consum+grarem(ji,jj,jk,jm)*rn_sigzoo(jm)         &
     &          +reszoo(ji,jj,jk,jm)
            END DO
            consump= consum+remgon(ji,jj,jk)*ratc2n-bactge(ji,jj,jk)*remgoc(ji,jj,jk)
! =(1.-bactge(ji,jj,jk))*remgon(ji,jj,jk)*ratc2n&
!     &        +bactge(ji,jj,jk)*(remgon(ji,jj,jk)*ratc2n-remgoc(ji,jj,jk))
            consum = consum+(1.-bactge(ji,jj,jk))*remgoc(ji,jj,jk)
            delo2(ji,jj,jk) = min(rato2c*consum,delo2(ji,jj,jk))
!
            denitr(ji,jj,jk) = 0.8*(rato2c*consum-delo2(ji,jj,jk))
!            prbn2s = denitr(ji,jj,jk)*rn_denn2o
!            trn(ji,jj,jk,jpn2o)= trn(ji,jj,jk,jpn2o)+prbn2s
!            trc3d(ji,jj,jk,17) = prbn2s*rfactr*1e3
!            denno3=max(-0.5,(rn_doxn2o-trno2(ji,jj,jk))     &
!       &      /(2.*rn_doxn2o+trno2(ji,jj,jk)))
!            denno3=(sin(denno3*rpi)+1.)/2.
!            degn2o=min(trn(ji,jj,jk,jpn2o), &
!       &      0.8*rato2c*consum*denno3*rn_degn2o)
!            trn(ji,jj,jk,jpn2o)= trn(ji,jj,jk,jpn2o)-degn2o
!            trc3d(ji,jj,jk,18) = degn2o*rfactr*1e3
!
! 0.5 N2 + 1.25 O2 + 0.5 H2O -> HNO3
!
            trn(ji,jj,jk,jpoxy) = trn(ji,jj,jk,jpoxy)-delo2(ji,jj,jk)  &
     &        -(prophy(ji,jj,jk,jpfix,1)+docphy(ji,jj,jk,jpfix))*ratn2c* &
     &        (1.-dinpft(ji,jj,jk,jpfix))*1.25
!
            prodt(ji,jj,jk) = 0.
!
! Evolution of phytoplankton
! --------------------------
!
!     Loop over phytoplankton types
!
            DO jl = jpdia, jpdia+jppft-1
!
              trn(ji,jj,jk,jppo4) = trn(ji,jj,jk,jppo4)                &
     &          -prophy(ji,jj,jk,jl,1)-docphy(ji,jj,jk,jl)
!
              trn(ji,jj,jk,jpdin) = trn(ji,jj,jk,jpdin)                &
     &          -(prophy(ji,jj,jk,jl,1)+docphy(ji,jj,jk,jl))           &
     &          *ratn2c*dinpft(ji,jj,jk,jl)
!
              prodt(ji,jj,jk) = prodt(ji,jj,jk)+prophy(ji,jj,jk,jl,1)  &
     &          +docphy(ji,jj,jk,jl)
              ppint(ji,jj) = ppint(ji,jj)+prophy(ji,jj,jk,jl,1)*e3t_0(ji,jj,jk)*1e3*rfactr
              pptdoc(ji,jj,jk) = pptdoc(ji,jj,jk)+docphy(ji,jj,jk,jl)*1e3*rfactr
!
! Loop over elements (C, Fe, Chl)
!
              DO jm = 0, 2
!
                trn(ji,jj,jk,jl+jm*jppft) = trn(ji,jj,jk,jl+jm*jppft)  &
     &            +prophy(ji,jj,jk,jl,jm+1)                            &
     &            -resphy(ji,jj,jk,jl,jm+1)
                DO jn = 1, jpzft
                trn(ji,jj,jk,jl+jm*jppft) = trn(ji,jj,jk,jl+jm*jppft)  &
     &            -grazoo(ji,jj,jk,jn,jl)*stofoo(ji,jj,jk,jl,jm+1)
                END DO
              END DO
!
              trn(ji,jj,jk,jpdoc) = trn(ji,jj,jk,jpdoc)+docphy(ji,jj,jk,jl)
!
              trn(ji,jj,jk,jpoxy) = trn(ji,jj,jk,jpoxy)                &
     &          +rato2c*(prophy(ji,jj,jk,jl,1)+docphy(ji,jj,jk,jl))
!
              trn(ji,jj,jk,jpfer) = trn(ji,jj,jk,jpfer)                &
     &          -prophy(ji,jj,jk,jl,2)+resphy(ji,jj,jk,jl,2)
            END DO
!
! Evolution of Calcite
! --------------------
!
! Production of attached CaCO3 (no tracer, but subtracted from DIC and Alk)
!
            prcaca(ji,jj,jk) =                                          &
     &        rn_coccal*prophy(ji,jj,jk,jpcoc,1)
            proara(ji,jj,jk) = rn_pteara*mgezoo(ji,jj,jk,2)             &
     &         *grazoc(ji,jj,jk,2)
!
! Transfer of attached CaCO3 to sinking CaCO3 by loss processes.
!
            losara(ji,jj,jk) = rn_pteara*reszoo(ji,jj,jk,2)
            DO jm = 3, jpzft
              losara(ji,jj,jk) = losara(ji,jj,jk)+rn_pteara*grazoo(ji,jj,jk,jm,jppte)
            END DO
            trn(ji,jj,jk,jpara) = trn(ji,jj,jk,jpara)+(1.-rn_disara)*losara(ji,jj,jk)  &
     &        +(snkara(ji,jj,jk)-snkara(ji,jj,jk+1))/e3t_0(ji,jj,jk)
            trn(ji,jj,jk,jpdic) = trn(ji,jj,jk,jpdic)-proara(ji,jj,jk)  &
     &        +rn_disara*losara(ji,jj,jk)
            trn(ji,jj,jk,jppiic) = trn(ji,jj,jk,jppiic)-proara(ji,jj,jk)  &
     &        +rn_disara*losara(ji,jj,jk)
            trn(ji,jj,jk,jptal) = trn(ji,jj,jk,jptal)-2.*proara(ji,jj,jk) &
     &        +2.*rn_disara*losara(ji,jj,jk)
!
! Transfer of attached CaCO3 to sinking CaCO3 by loss processes.
!
            loscal(ji,jj,jk) = rn_coccal*resphy(ji,jj,jk,jpcoc,1)
            DO jm = 1, jpzft
              loscal(ji,jj,jk) = loscal(ji,jj,jk)+rn_coccal*grazoo(ji,jj,jk,jm,jpcoc)
            END DO
!
            trn(ji,jj,jk,jpcal) = trn(ji,jj,jk,jpcal)+(1.-rn_discal)*loscal(ji,jj,jk) &
     &        +(snkcal(ji,jj,jk)-snkcal(ji,jj,jk+1))/e3t_0(ji,jj,jk)
            d14pro(ji,jj,jk) = rn_bfrd14/10.*prodt(ji,jj,jk)/max(trn(ji,jj,jk,jpdic),rtrn)
            d14res(ji,jj,jk) = (trn(ji,jj,jk,jpd14)-(trn(ji,jj,1,jpd14)-rn_bfrd14/10.)) &
     &        *consum/max(trn(ji,jj,jk,jpdic),rtrn)
            trn(ji,jj,jk,jpd14) = trn(ji,jj,jk,jpd14)+(d14pro(ji,jj,jk)-d14res(ji,jj,jk))*tmask(ji,jj,jk)
!
! Evolution of PO4
! ----------------
!
            trn(ji,jj,jk,jppo4) = trn(ji,jj,jk,jppo4)+consump           &
     &        +deppo4(ji,jj,jk)
!
            trn(ji,jj,jk,jpdin) = trn(ji,jj,jk,jpdin)-denitr(ji,jj,jk) &
     &        +consump*ratn2c+depnit(ji,jj,jk)+atmdin(ji,jj,jk)
!
            prodt(ji,jj,jk) = prodt(ji,jj,jk)-consum
!    
! Evolution of zooplankton
! -----------------------------
!     
            DO jm = 1, jpzft
              trn(ji,jj,jk,jpbac+jm) = trn(ji,jj,jk,jpbac+jm)           &
     &      -reszoo(ji,jj,jk,jm)+mgezoo(ji,jj,jk,jm)*grazoc(ji,jj,jk,jm)
            END DO
            DO jm = 0, jpzft
              DO jn = 1, jpzft
                trn(ji,jj,jk,jpbac+jm) = trn(ji,jj,jk,jpbac+jm)         &
     &            -grazoo(ji,jj,jk,jn,jpbac+jm)
              END DO
            END DO
!
! Evolution of Macrozooplankton
! -----------------------------
!
            trn(ji,jj,jk,jpmac) = trn(ji,jj,jk,jpmac)-tormac(ji,jj,jk)
            trn(ji,jj,jk,jpgel) = trn(ji,jj,jk,jpgel)-torgel(ji,jj,jk)
!
! Evolution of Bacteria
! ---------------------
!
            trn(ji,jj,jk,jpbac) = trn(ji,jj,jk,jpbac)+bactge(ji,jj,jk)*&
     &        (remdoc(ji,jj,jk)+rempoc(ji,jj,jk)+remgoc(ji,jj,jk))         &
     &        -resbac(ji,jj,jk)
!    
! Evolution of DOC
! ----------------
!     
            trn(ji,jj,jk,jpdoc) = trn(ji,jj,jk,jpdoc)-remdoc(ji,jj,jk)  &
     &        +resbac(ji,jj,jk)/3.0  &
     &        -xaggdoc(ji,jj,jk)-xaggdoc2(ji,jj,jk)+depdoc(ji,jj,jk)    
            gongoc = trn(ji,jj,jk,jpgon)/max(trn(ji,jj,jk,jpgoc),rtrn)
            DO jm = 1, jpzft
              trn(ji,jj,jk,jpdoc) = trn(ji,jj,jk,jpdoc)                 &
     &          +grarem(ji,jj,jk,jm)*(1.-rn_sigzoo(jm))
              trn(ji,jj,jk,jppoc) = trn(ji,jj,jk,jppoc)                 &
     &          -grazoo(ji,jj,jk,jm,jppoc)
              trn(ji,jj,jk,jpgoc) = trn(ji,jj,jk,jpgoc)                 &
     &          -grazoo(ji,jj,jk,jm,jpgoc)
              trn(ji,jj,jk,jpgon) = trn(ji,jj,jk,jpgon)                 &
     &          -grazoo(ji,jj,jk,jm,jpgoc)*gongoc
              trn(ji,jj,jk,jppoc+nn_sizzoo(jm))=                        &
     &          trn(ji,jj,jk,jppoc+nn_sizzoo(jm))+grapoc(ji,jj,jk,jm)
              trn(ji,jj,jk,jpfer) = trn(ji,jj,jk,jpfer)                  &
     &          +grafer(ji,jj,jk,jm)+ferat3*reszoo(ji,jj,jk,jm)
              trn(ji,jj,jk,jpsfe) = trn(ji,jj,jk,jpsfe)                 &
     &          -grazoo(ji,jj,jk,jm,jppoc)*stofoo(ji,jj,jk,jppoc,2)
              trn(ji,jj,jk,jpbfe) = trn(ji,jj,jk,jpbfe)                 &
     &          -grazoo(ji,jj,jk,jm,jpgoc)*stofoo(ji,jj,jk,jpgoc,2)
              trn(ji,jj,jk,jpsfe+nn_sizzoo(jm))=                        &
     &          trn(ji,jj,jk,jpsfe+nn_sizzoo(jm))                       &
     &          +rn_unazoo(jm)*grazof(ji,jj,jk,jm)
            END DO
            dgongoc=gongoc-ratn2c
            DO jm = 2, jpzft
              trn(ji,jj,jk,jpgon)=trn(ji,jj,jk,jpgon)                   &
     &          +float(nn_sizzoo(jm))*(grapoc(ji,jj,jk,jm)*ratn2c       &
     &          +grazoo(ji,jj,jk,jm,jpgoc)*rn_unazoo(jm)*dgongoc)
            END DO
!
! Evolution of small POC
! ----------------------
!     
            trn(ji,jj,jk,jppoc) = trn(ji,jj,jk,jppoc)-rempoc(ji,jj,jk)   &
     &        +resbac(ji,jj,jk)/3.0  &
     &        +(snkpoc(ji,jj,jk)-snkpoc(ji,jj,jk+1))/e3t_0(ji,jj,jk) &
     &        -xagg(ji,jj,jk)+xaggdoc(ji,jj,jk)+deppoc(ji,jj,jk)
!    
! Evolution of big POC
! --------------------
!
            zgrdmp =                                                    &
     &        torgel(ji,jj,jk)+                                         &
     &        tormac(ji,jj,jk)+xagg(ji,jj,jk)+xaggdoc2(ji,jj,jk)
            trn(ji,jj,jk,jpgoc) = trn(ji,jj,jk,jpgoc)-remgoc(ji,jj,jk)   &
     &        +zgrdmp                                                   &
     &        +(snkgoc(ji,jj,jk)-snkgoc(ji,jj,jk+1))/e3t_0(ji,jj,jk)
            trn(ji,jj,jk,jpgon) = trn(ji,jj,jk,jpgon)-remgon(ji,jj,jk)   &
     &        +zgrdmp*ratn2c                                            &
     &        +(snkgon(ji,jj,jk)-snkgon(ji,jj,jk+1))/e3t_0(ji,jj,jk)
!
! Evolution of dissolved Iron
! ---------------------------
!
            trn(ji,jj,jk,jpfer) = trn(ji,jj,jk,jpfer)                  &
     &        +rbafer(ji,jj,jk)-ubafer(ji,jj,jk)                       &
     &        -xscave(ji,jj,jk)+irondep(ji,jj,jk)+depfer(ji,jj,jk)
!
! Evolution of small biogenic Iron
! --------------------------------
!
            trn(ji,jj,jk,jpsfe) = trn(ji,jj,jk,jpsfe)-remsfe(ji,jj,jk)    &
     &        -xaggfe(ji,jj,jk)+scasfe(ji,jj,jk)                       &
     &        +(snksfe(ji,jj,jk)-snksfe(ji,jj,jk+1))/e3t_0(ji,jj,jk) &
     &        +deppoc(ji,jj,jk)*ferat3
!
! Evolution of big biogenic Iron
! ------------------------------
!
            trn(ji,jj,jk,jpbfe) = trn(ji,jj,jk,jpbfe)                  &
     &        +ferat3*tormac(ji,jj,jk)              &
     &        +xaggfe(ji,jj,jk)+scabfe(ji,jj,jk)-rembfe(ji,jj,jk)       &
     &        +(snkbfe(ji,jj,jk)-snkbfe(ji,jj,jk+1))/e3t_0(ji,jj,jk)
!
! Evolution of biogenic Silica
! ----------------------------
!
            trn(ji,jj,jk,jpbsi) = trn(ji,jj,jk,jpbsi)                  &
     &        +prorca3(ji,jj,jk)-losbsi(ji,jj,jk)
!
! Evolution of sinking biogenic silica
! ------------------------------------
!
            trn(ji,jj,jk,jpdsi)=trn(ji,jj,jk,jpdsi)     &
     &        +losbsi(ji,jj,jk)-remdsi(ji,jj,jk)          &
     &        +(snkdsi(ji,jj,jk)-snkdsi(ji,jj,jk+1))/e3t_0(ji,jj,jk)
!
! Evolution of dissolved Silica
! -----------------------------
!
            trn(ji,jj,jk,jpsil) = trn(ji,jj,jk,jpsil)                  &
     &        -prorca3(ji,jj,jk)+remdsi(ji,jj,jk)                        &
     &        +depsil(ji,jj,jk)+sidep(ji,jj,jk)
!
! Consumption of Total (12C)O2
! ----------------------------
!     
            trn(ji,jj,jk,jpdic) = trn(ji,jj,jk,jpdic)                   &
     &        -prodt(ji,jj,jk)-prcaca(ji,jj,jk)+rn_discal*loscal(ji,jj,jk) &
     &        +depdic(ji,jj,jk)
            trn(ji,jj,jk,jppiic) = trn(ji,jj,jk,jppiic)                   &
     &        -prodt(ji,jj,jk)-prcaca(ji,jj,jk)+rn_discal*loscal(ji,jj,jk) &
     &        +depdic(ji,jj,jk)
!     
! Consumption of alkalinity due to ca++ uptake and increase 
! of alkalinity due to nitrate consumption during organic 
! soft tissue production
! ------------------------------------------------------------------
!     
            trn(ji,jj,jk,jptal) = trn(ji,jj,jk,jptal)                  &
     &        +alknut*prodt(ji,jj,jk)-2.*prcaca(ji,jj,jk)+             &
     &        2.*rn_discal*loscal(ji,jj,jk)+depdic(ji,jj,jk)+denitr(ji,jj,jk)
!            deln2o = nitrif*rn_aoun2o
!            trn(ji,jj,jk,jpn2o)= trn(ji,jj,jk,jpn2o)+deln2o
!            trc3d(ji,jj,jk,16) = deln2o*rfactr*1e3
          END DO
        END DO
      END DO
!
! Save additional diagnostics

        out3d(:,:,:)=(snkgoc(:,:,:)+snkpoc(:,:,:))*1e3*rfactr   
!        where (tmask(:,:,:) .eq. 0 )
!          out3d = ncf_fill
!        end where
        CALL iom_put("EXP", out3d )

!        where (tmask(:,:,:) .eq. 0 )
!          out3d = ncf_fill
!        else where
          out3d = grazoc(:,:,:,1)*1e3*rfactr
!        end where
        CALL iom_put("GRAMIC", out3d(:,:,:) )
!        
        DO jl = jpdia, jpdia+jppft-1   
            ppt(2:nlci-1,2:nlcj-1,1:jpk-1) = &
     &           ppt(2:nlci-1,2:nlcj-1,1:jpk-1)+ &
     &           prophy(2:nlci-1,2:nlcj-1,1:jpk-1,jl,1)*rfactr*1e3
            tchl(2:nlci-1,2:nlcj-1,1:jpk-1) = &
     &           tchl(2:nlci-1,2:nlcj-1,1:jpk-1)+ &
     &           trn(2:nlci-1,2:nlcj-1,1:jpk-1,jl+2*jppft)
            grazing(2:nlci-1,2:nlcj-1,1:jpk-1,1) = & 
     &           grazing(2:nlci-1,2:nlcj-1,1:jpk-1,1)+ &
     &           grazoo(2:nlci-1,2:nlcj-1,1:jpk-1,1,jl)*rfactr*1e3
            grazing(2:nlci-1,2:nlcj-1,1:jpk-1,2) = &
     &           grazing(2:nlci-1,2:nlcj-1,1:jpk-1,2) +  &
     &           grazoo(2:nlci-1,2:nlcj-1,1:jpk-1,3,jl)*rfactr*1e3
            grazing(2:nlci-1,2:nlcj-1,1:jpk-1,3) = &
     &           grazing(2:nlci-1,2:nlcj-1,1:jpk-1,3) + &
     &           grazoo(2:nlci-1,2:nlcj-1,1:jpk-1,5,jl)*rfactr*1e3
        END DO
!        where (tmask(:,:,:) .eq. 0 )
!          ppt(:,:,:) = ncf_fill
!          tchl(:,:,:) = ncf_fill
!          grazing(:,:,:,1) = ncf_fill
!          grazing(:,:,:,2) = ncf_fill
!          grazing(:,:,:,3) = ncf_fill
!        end where
        CALL iom_put("PPT", ppt(:,:,:) )
        CALL iom_put("PPTDOC",pptdoc)
        CALL iom_put("TChl",tchl)
        CALL iom_put("Detrit",trophic(:,:,:,1))
        CALL iom_put("Carniv",trophic(:,:,:,2))
        CALL iom_put("Herbiv",trophic(:,:,:,3))
        CALL iom_put("GRAMICPHY", grazing(:,:,:,1) )
        CALL iom_put("GRAMESPHY", grazing(:,:,:,2) )
        CALL iom_put("GRAMACPHY", grazing(:,:,:,3) )
      
!        where (tmask(:,:,:) .eq. 0 )
!          out3d = ncf_fill 
!        else where
          out3d = (prophy(:,:,:,jpfix,1)             &
     &        +docphy(:,:,:,jpfix))*ratn2c*(1.-dinpft(:,:,:,jpfix))*rfactr*1.e3
!        end where    
        CALL iom_put("nitrfix", out3d) 
  
!        where (tmask(:,:,:) .eq. 0 )
!          out3d = ncf_fill
!        else where
          out3d = denitr(:,:,:)*rfactr*1e3 
!        end where    
        CALL iom_put("denitr", out3d )

!        where (tmask(:,:,:) .eq. 0 )
!          out3d = ncf_fill
!        else where
          out3d = delo2(:,:,:)*1e3*rfactr
!        end where 
        CALL iom_put("DELO2", out3d)
          out3d = snkdsi(:,:,:)*1e3*rfactr
        CALL iom_put("sinksil", out3d)
          out3d = (snkcal(:,:,:)+snkara(:,:,:))*1e3*rfactr
        CALL iom_put("ExpCO3", out3d )
          out3d = snkara(:,:,:)*1e3*rfactr
        CALL iom_put("ExpARA", out3d )
          out3d = grazoc(:,:,:,5)*1e3*rfactr
        CALL iom_put("GRAMAC", out3d(:,:,:) )
          out3d = proara(:,:,:)*rfactr*1e3
        CALL iom_put("proara", out3d )
          out3d = grazoc(:,:,:,2)*1e3*rfactr
        CALL iom_put("GRAPTE", out3d )
          out3d = grazoc(:,:,:,3)*1e3*rfactr
        CALL iom_put("GRAMES", out3d(:,:,:) )
          out3d = grazoc(:,:,:,4)*1e3*rfactr
        CALL iom_put( "GRAGEL", out3d(:,:,:) )
          out3d = rn_coccal*prophy(:,:,:,jpcoc,1)*rfactr*1e3
        CALL iom_put("prococ", out3d )
        CALL iom_put("probsi",prorca3*rfactr*1e3)
        CALL iom_put("losbsi",losbsi*rfactr*1e3)
        CALL iom_put("remdsi",remdsi*rfactr*1e3)
        CALL iom_put("D14PRO",d14pro)
        CALL iom_put("D14RES",d14res)
        CALL iom_put("PPINT",ppint)
      RETURN
      END
