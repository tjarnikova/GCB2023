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

MODULE step_oce
   !!======================================================================
   !!                       ***  MODULE step_oce  ***
   !! Ocean time-stepping : module used in both initialisation phase and time stepping
   !!======================================================================
   !! History :   3.3  ! 2010-08  (C. Ethe)  Original code - reorganisation of the initial phase
   !!----------------------------------------------------------------------
   USE oce              ! ocean dynamics and tracers variables
   USE dom_oce          ! ocean space and time domain variables
   USE zdf_oce          ! ocean vertical physics variables
   USE ldftra_oce       ! ocean tracer   - trends
   USE ldfdyn_oce       ! ocean dynamics - trends
   USE divcur           ! hor. divergence and curl      (div & cur routines)
   USE in_out_manager   ! I/O manager
   USE iom              !
   USE lbclnk
   USE restart          ! restart
   USE xios

   USE daymod           ! calendar                         (day     routine)

   USE sbcmod           ! surface boundary condition       (sbc     routine)
   USE sbcrnf           ! surface boundary condition: runoff variables
   USE sbccpl           ! surface boundary condition: coupled formulation (call send at end of step)
   USE sbc_oce          ! surface boundary condition: ocean
   USE sbctide          ! Tide initialisation
   USE sbcapr           ! surface boundary condition: ssh_ib required by bdydta 

   USE traqsr           ! solar radiation penetration      (tra_qsr routine)
   USE trasbc           ! surface boundary condition       (tra_sbc routine)
   USE trabbc           ! bottom boundary condition        (tra_bbc routine)
   USE trabbl           ! bottom boundary layer            (tra_bbl routine)
   USE tradmp           ! internal damping                 (tra_dmp routine)
   USE traadv           ! advection scheme control     (tra_adv_ctl routine)
   USE traldf           ! lateral mixing                   (tra_ldf routine)
   !   zdfkpp           ! KPP non-local tracer fluxes      (tra_kpp routine)
   USE trazdf           ! vertical mixing                  (tra_zdf routine)
   USE tranxt           ! time-stepping                    (tra_nxt routine)
   USE tranpc           ! non-penetrative convection       (tra_npc routine)

   USE eosbn2           ! equation of state                (eos_bn2 routine)

   USE dynadv           ! advection                        (dyn_adv routine)
   USE dynbfr           ! Bottom friction terms            (dyn_bfr routine)
   USE dynvor           ! vorticity term                   (dyn_vor routine)
   USE dynhpg           ! hydrostatic pressure grad.       (dyn_hpg routine)
   USE dynldf           ! lateral momentum diffusion       (dyn_ldf routine)
   USE dynzdf           ! vertical diffusion               (dyn_zdf routine)
   USE dynspg_oce       ! surface pressure gradient        (dyn_spg routine)
   USE dynspg           ! surface pressure gradient        (dyn_spg routine)
   USE dynnept          ! simp. form of Neptune effect(dyn_nept_cor routine)

   USE dynnxt           ! time-stepping                    (dyn_nxt routine)

   USE stopar           ! Stochastic parametrization       (sto_par routine)
   USE stopts 

   USE bdy_par          ! for lk_bdy
   USE bdy_oce          ! for dmp logical
   USE bdydta           ! open boundary condition data     (bdy_dta routine)
   USE bdytra           ! bdy cond. for tracers            (bdy_tra routine)
   USE bdydyn3d         ! bdy cond. for baroclinic vel.  (bdy_dyn3d routine)

   USE sshwzv           ! vertical velocity and ssh        (ssh_nxt routine)
   !                                                       (ssh_swp routine)
   !                                                       (wzv     routine)
   USE domvvl           ! variable vertical scale factors  (dom_vvl_sf_nxt routine)
   !                                                       (dom_vvl_sf_swp routine)

   USE ldfslp           ! iso-neutral slopes               (ldf_slp routine)
   USE ldfeiv           ! eddy induced velocity coef.      (ldf_eiv routine)
   USE ldftra_smag      ! Smagirinsky diffusion            (ldftra_smag routine)
   USE ldfdyn_smag      ! Smagorinsky viscosity            (ldfdyn_smag routine) 

   USE zdftmx           ! tide-induced vertical mixing     (zdf_tmx routine)
   USE zdfbfr           ! bottom friction                  (zdf_bfr routine)
   USE zdftke           ! TKE vertical mixing              (zdf_tke routine)
   USE zdfgls           ! GLS vertical mixing              (zdf_gls routine)
   USE zdfkpp           ! KPP vertical mixing              (zdf_kpp routine)
   USE zdfddm           ! double diffusion mixing          (zdf_ddm routine)
   USE zdfevd           ! enhanced vertical diffusion      (zdf_evd routine)
   USE zdfric           ! Richardson vertical mixing       (zdf_ric routine)
   USE zdfmxl           ! Mixed-layer depth                (zdf_mxl routine)

   USE zpshde           ! partial step: hor. derivative     (zps_hde routine)

   USE diawri           ! Standard run outputs             (dia_wri routine)
   USE diaptr           ! poleward transports              (dia_ptr routine)
   USE diadct           ! sections transports              (dia_dct routine)
   USE diaar5           ! AR5 diagnosics                   (dia_ar5 routine)
   USE diahth           ! thermocline depth                (dia_hth routine)
   USE diafwb           ! freshwater budget                (dia_fwb routine)
   USE diahsb           ! heat, salt and volume budgets    (dia_hsb routine)
   USE diaharm
   USE flo_oce          ! floats variables
   USE floats           ! floats computation               (flo_stp routine)

   USE crsfld           ! Standard output on coarse grid   (crs_fld routine)

   USE asminc           ! assimilation increments      (tra_asm_inc routine)
   !                                                   (dyn_asm_inc routine)
   USE asmbkg
   USE stpctl           ! time stepping control            (stp_ctl routine)
   USE prtctl           ! Print control                    (prt_ctl routine)

   USE diaobs           ! Observation operator

   USE timing           ! Timing

   USE trcstp           ! passive tracer time-stepping      (trc_stp routine)
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: step_oce.F90 5501 2015-06-29 10:08:15Z deazer $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE step_oce
