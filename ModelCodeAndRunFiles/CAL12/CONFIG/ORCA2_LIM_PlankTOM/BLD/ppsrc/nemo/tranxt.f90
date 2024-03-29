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

MODULE tranxt
   !!======================================================================
   !!                       ***  MODULE  tranxt  ***
   !! Ocean active tracers:  time stepping on temperature and salinity
   !!======================================================================
   !! History :  OPA  !  1991-11  (G. Madec)  Original code
   !!            7.0  !  1993-03  (M. Guyon)  symetrical conditions
   !!            8.0  !  1996-02  (G. Madec & M. Imbard)  opa release 8.0
   !!             -   !  1996-04  (A. Weaver)  Euler forward step
   !!            8.2  !  1999-02  (G. Madec, N. Grima)  semi-implicit pressure grad.
   !!  NEMO      1.0  !  2002-08  (G. Madec)  F90: Free form and module
   !!             -   !  2002-11  (C. Talandier, A-M Treguier) Open boundaries
   !!             -   !  2005-04  (C. Deltel) Add Asselin trend in the ML budget
   !!            2.0  !  2006-02  (L. Debreu, C. Mazauric) Agrif implementation
   !!            3.0  !  2008-06  (G. Madec)  time stepping always done in trazdf
   !!            3.1  !  2009-02  (G. Madec, R. Benshila)  re-introduce the vvl option
   !!            3.3  !  2010-04  (M. Leclair, G. Madec)  semi-implicit hpg with asselin filter + modified LF-RA
   !!             -   !  2010-05  (C. Ethe, G. Madec)  merge TRC-TRA
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_nxt       : time stepping on tracers
   !!   tra_nxt_fix   : time stepping on tracers : fixed    volume case
   !!   tra_nxt_vvl   : time stepping on tracers : variable volume case
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE sbc_oce         ! surface boundary condition: ocean
   USE sbcrnf          ! river runoffs
   USE zdf_oce         ! ocean vertical mixing
   USE domvvl          ! variable volume
   USE dynspg_oce      ! surface     pressure gradient variables
   USE dynhpg          ! hydrostatic pressure gradient 
   USE trd_oce         ! trends: ocean variables
   USE trdtra          ! trends manager: tracers 
   USE traqsr          ! penetrative solar radiation (needed for nksr)
   USE phycst          ! physical constant
   USE ldftra_oce      ! lateral physics on tracers
   USE bdy_oce         ! BDY open boundary condition variables
   USE bdytra          ! open boundary condition (bdy_tra routine)
   !
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control
   USE wrk_nemo        ! Memory allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_nxt       ! routine called by step.F90
   PUBLIC   tra_nxt_fix   ! to be used in trcnxt
   PUBLIC   tra_nxt_vvl   ! to be used in trcnxt

   REAL(wp) ::   rbcp   ! Brown & Campana parameters for semi-implicit hpg

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
   !! NEMO/OPA 3.3 , NEMO-Consortium (2010) 
   !! $Id: tranxt.F90 5467 2015-06-24 10:07:54Z jchanut $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_nxt( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tranxt  ***
      !!
      !! ** Purpose :   Apply the boundary condition on the after temperature  
      !!             and salinity fields, achieved the time stepping by adding
      !!             the Asselin filter on now fields and swapping the fields.
      !! 
      !! ** Method  :   At this stage of the computation, ta and sa are the 
      !!             after temperature and salinity as the time stepping has
      !!             been performed in trazdf_imp or trazdf_exp module.
      !!
      !!              - Apply lateral boundary conditions on (ta,sa) 
      !!             at the local domain   boundaries through lbc_lnk call, 
      !!             at the one-way open boundaries (lk_bdy=T), 
      !!             at the AGRIF zoom   boundaries (lk_agrif=T)
      !!
      !!              - Update lateral boundary conditions on AGRIF children
      !!             domains (lk_agrif=T)
      !!
      !! ** Action  : - (tb,sb) and (tn,sn) ready for the next time step
      !!              - (ta,sa) time averaged (t,s)   (ln_dynhpg_imp = T)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !!
      INTEGER  ::   jk, jn    ! dummy loop indices
      REAL(wp) ::   zfact     ! local scalars
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ztrdt, ztrds
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'tra_nxt')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_nxt : achieve the time stepping by Asselin filter and array swap'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
         !
         rbcp = 0.25_wp * (1._wp + atfp) * (1._wp + atfp) * ( 1._wp - atfp)      ! Brown & Campana parameter for semi-implicit hpg
      ENDIF

      ! Update after tracer on domain lateral boundaries
      ! 
      CALL lbc_lnk( tsa(:,:,:,jp_tem), 'T', 1._wp )      ! local domain boundaries  (T-point, unchanged sign)
      CALL lbc_lnk( tsa(:,:,:,jp_sal), 'T', 1._wp )
      !
 
      ! set time step size (Euler/Leapfrog)
      IF( neuler == 0 .AND. kt == nit000 ) THEN   ;   r2dtra(:) =     rdttra(:)      ! at nit000             (Euler)
      ELSEIF( kt <= nit000 + 1 )           THEN   ;   r2dtra(:) = 2._wp* rdttra(:)      ! at nit000 or nit000+1 (Leapfrog)
      ENDIF

      ! trends computation initialisation
      IF( l_trdtra )   THEN                    ! store now fields before applying the Asselin filter
         CALL wrk_alloc( jpi, jpj, jpk, ztrdt, ztrds )
         ztrdt(:,:,:) = tsn(:,:,:,jp_tem) 
         ztrds(:,:,:) = tsn(:,:,:,jp_sal)
         IF( ln_traldf_iso ) THEN              ! diagnose the "pure" Kz diffusive trend 
            CALL trd_tra( kt, 'TRA', jp_tem, jptra_zdfp, ztrdt )
            CALL trd_tra( kt, 'TRA', jp_sal, jptra_zdfp, ztrds )
         ENDIF
      ENDIF

      IF( neuler == 0 .AND. kt == nit000 ) THEN       ! Euler time-stepping at first time-step (only swap)
         DO jn = 1, jpts
            DO jk = 1, jpkm1
               tsn(:,:,jk,jn) = tsa(:,:,jk,jn)    
            END DO
         END DO
      ELSE                                            ! Leap-Frog + Asselin filter time stepping
         !
         IF( lk_vvl )  THEN   ;   CALL tra_nxt_vvl( kt, nit000, rdttra, 'TRA', tsb, tsn, tsa,   &
           &                                                              sbc_tsc, sbc_tsc_b, jpts )  ! variable volume level (vvl) 
         ELSE                 ;   CALL tra_nxt_fix( kt, nit000,         'TRA', tsb, tsn, tsa, jpts )  ! fixed    volume level 
         ENDIF
      ENDIF 
      !
      !
      ! trends computation
      IF( l_trdtra ) THEN      ! trend of the Asselin filter (tb filtered - tb)/dt     
         DO jk = 1, jpkm1
            zfact = 1._wp / r2dtra(jk)             
            ztrdt(:,:,jk) = ( tsb(:,:,jk,jp_tem) - ztrdt(:,:,jk) ) * zfact
            ztrds(:,:,jk) = ( tsb(:,:,jk,jp_sal) - ztrds(:,:,jk) ) * zfact
         END DO
         CALL trd_tra( kt, 'TRA', jp_tem, jptra_atf, ztrdt )
         CALL trd_tra( kt, 'TRA', jp_sal, jptra_atf, ztrds )
         CALL wrk_dealloc( jpi, jpj, jpk, ztrdt, ztrds )
      END IF
      !
      !                        ! control print
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=tsn(:,:,:,jp_tem), clinfo1=' nxt  - Tn: ', mask1=tmask,   &
         &                       tab3d_2=tsn(:,:,:,jp_sal), clinfo2=       ' Sn: ', mask2=tmask )
      !
      IF( nn_timing == 1 )   CALL timing_stop('tra_nxt')
      !
   END SUBROUTINE tra_nxt


   SUBROUTINE tra_nxt_fix( kt, kit000, cdtype, ptb, ptn, pta, kjpt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tra_nxt_fix  ***
      !!
      !! ** Purpose :   fixed volume: apply the Asselin time filter and 
      !!                swap the tracer fields.
      !! 
      !! ** Method  : - Apply a Asselin time filter on now fields.
      !!              - save in (ta,sa) an average over the three time levels 
      !!             which will be used to compute rdn and thus the semi-implicit
      !!             hydrostatic pressure gradient (ln_dynhpg_imp = T)
      !!              - swap tracer fields to prepare the next time_step.
      !!                This can be summurized for tempearture as:
      !!             ztm = tn + rbcp * [ta -2 tn + tb ]       ln_dynhpg_imp = T
      !!             ztm = 0                                   otherwise
      !!                   with rbcp=1/4 * (1-atfp^4) / (1-atfp)
      !!             tb  = tn + atfp*[ tb - 2 tn + ta ]
      !!             tn  = ta  
      !!             ta  = ztm       (NB: reset to 0 after eos_bn2 call)
      !!
      !! ** Action  : - (tb,sb) and (tn,sn) ready for the next time step
      !!              - (ta,sa) time averaged (t,s)   (ln_dynhpg_imp = T)
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in   )                               ::   kt       ! ocean time-step index
      INTEGER         , INTENT(in   )                               ::   kit000   ! first time step index
      CHARACTER(len=3), INTENT(in   )                               ::   cdtype   ! =TRA or TRC (tracer indicator)
      INTEGER         , INTENT(in   )                               ::   kjpt     ! number of tracers
      REAL(wp)        , INTENT(inout), DIMENSION(jpi,jpj,jpk,kjpt)  ::   ptb      ! before tracer fields
      REAL(wp)        , INTENT(inout), DIMENSION(jpi,jpj,jpk,kjpt)  ::   ptn      ! now tracer fields
      REAL(wp)        , INTENT(inout), DIMENSION(jpi,jpj,jpk,kjpt)  ::   pta      ! tracer trend
      !
      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      LOGICAL  ::   ll_tra_hpg       ! local logical
      REAL(wp) ::   ztn, ztd         ! local scalars
      !!----------------------------------------------------------------------

      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_nxt_fix : time stepping', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      IF( cdtype == 'TRA' )  THEN   ;   ll_tra_hpg = ln_dynhpg_imp    ! active  tracers case  and  semi-implicit hpg    
      ELSE                          ;   ll_tra_hpg = .FALSE.          ! passive tracers case or NO semi-implicit hpg
      ENDIF
      !
      DO jn = 1, kjpt
         !
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ztn = ptn(ji,jj,jk,jn)                                    
                  ztd = pta(ji,jj,jk,jn) - 2. * ztn + ptb(ji,jj,jk,jn)      !  time laplacian on tracers
                  !
                  ptb(ji,jj,jk,jn) = ztn + atfp * ztd                       ! ptb <-- filtered ptn 
                  ptn(ji,jj,jk,jn) = pta(ji,jj,jk,jn)                       ! ptn <-- pta
                  !
                  IF( ll_tra_hpg )   pta(ji,jj,jk,jn) = ztn + rbcp * ztd    ! pta <-- Brown & Campana average
               END DO
           END DO
         END DO
         !
      END DO
      !
   END SUBROUTINE tra_nxt_fix


   SUBROUTINE tra_nxt_vvl( kt, kit000, p2dt, cdtype, ptb, ptn, pta, psbc_tc, psbc_tc_b, kjpt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tra_nxt_vvl  ***
      !!
      !! ** Purpose :   Time varying volume: apply the Asselin time filter  
      !!                and swap the tracer fields.
      !! 
      !! ** Method  : - Apply a thickness weighted Asselin time filter on now fields.
      !!              - save in (ta,sa) a thickness weighted average over the three 
      !!             time levels which will be used to compute rdn and thus the semi-
      !!             implicit hydrostatic pressure gradient (ln_dynhpg_imp = T)
      !!              - swap tracer fields to prepare the next time_step.
      !!                This can be summurized for tempearture as:
      !!             ztm = ( e3t_n*tn + rbcp*[ e3t_b*tb - 2 e3t_n*tn + e3t_a*ta ] )   ln_dynhpg_imp = T
      !!                  /( e3t_n    + rbcp*[ e3t_b    - 2 e3t_n    + e3t_a    ] )   
      !!             ztm = 0                                                       otherwise
      !!             tb  = ( e3t_n*tn + atfp*[ e3t_b*tb - 2 e3t_n*tn + e3t_a*ta ] )
      !!                  /( e3t_n    + atfp*[ e3t_b    - 2 e3t_n    + e3t_a    ] )
      !!             tn  = ta 
      !!             ta  = zt        (NB: reset to 0 after eos_bn2 call)
      !!
      !! ** Action  : - (tb,sb) and (tn,sn) ready for the next time step
      !!              - (ta,sa) time averaged (t,s)   (ln_dynhpg_imp = T)
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in   )                               ::  kt       ! ocean time-step index
      INTEGER         , INTENT(in   )                               ::  kit000   ! first time step index
      REAL(wp)        , INTENT(in   ), DIMENSION(jpk)               ::  p2dt     ! time-step
      CHARACTER(len=3), INTENT(in   )                               ::  cdtype   ! =TRA or TRC (tracer indicator)
      INTEGER         , INTENT(in   )                               ::  kjpt     ! number of tracers
      REAL(wp)        , INTENT(inout), DIMENSION(jpi,jpj,jpk,kjpt)  ::  ptb      ! before tracer fields
      REAL(wp)        , INTENT(inout), DIMENSION(jpi,jpj,jpk,kjpt)  ::  ptn      ! now tracer fields
      REAL(wp)        , INTENT(inout), DIMENSION(jpi,jpj,jpk,kjpt)  ::  pta      ! tracer trend
      REAL(wp)        , INTENT(in   ), DIMENSION(jpi,jpj,kjpt)      ::  psbc_tc   ! surface tracer content
      REAL(wp)        , INTENT(in   ), DIMENSION(jpi,jpj,kjpt)      ::  psbc_tc_b ! before surface tracer content

      !!     
      LOGICAL  ::   ll_tra_hpg, ll_traqsr, ll_rnf   ! local logical
      INTEGER  ::   ji, jj, jk, jn              ! dummy loop indices
      REAL(wp) ::   zfact1, ztc_a , ztc_n , ztc_b , ztc_f , ztc_d    ! local scalar
      REAL(wp) ::   zfact2, ze3t_b, ze3t_n, ze3t_a, ze3t_f, ze3t_d   !   -      -
      !!----------------------------------------------------------------------
      !
      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_nxt_vvl : time stepping', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      IF( cdtype == 'TRA' )  THEN   
         ll_tra_hpg = ln_dynhpg_imp    ! active  tracers case  and  semi-implicit hpg
         ll_traqsr  = ln_traqsr        ! active  tracers case  and  solar penetration
         ll_rnf     = ln_rnf           ! active  tracers case  and  river runoffs
      ELSE                          
         ll_tra_hpg = .FALSE.          ! passive tracers case or NO semi-implicit hpg
         ll_traqsr  = .FALSE.          ! active  tracers case and NO solar penetration
         ll_rnf     = .FALSE.          ! passive tracers or NO river runoffs
      ENDIF
      !
      DO jn = 1, kjpt      
         DO jk = 1, jpkm1
            zfact1 = atfp * p2dt(jk)
            zfact2 = zfact1 / rau0
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ze3t_b = e3t_0(ji,jj,jk)
                  ze3t_n = e3t_0(ji,jj,jk)
                  ze3t_a = e3t_0(ji,jj,jk)
                  !                                         ! tracer content at Before, now and after
                  ztc_b  = ptb(ji,jj,jk,jn) * ze3t_b
                  ztc_n  = ptn(ji,jj,jk,jn) * ze3t_n
                  ztc_a  = pta(ji,jj,jk,jn) * ze3t_a
                  !
                  ze3t_d = ze3t_a - 2. * ze3t_n + ze3t_b
                  ztc_d  = ztc_a  - 2. * ztc_n  + ztc_b
                  !
                  ze3t_f = ze3t_n + atfp * ze3t_d
                  ztc_f  = ztc_n  + atfp * ztc_d
                  !
                  IF( jk == 1 ) THEN           ! first level 
                     ze3t_f = ze3t_f - zfact2 * ( emp_b(ji,jj) - emp(ji,jj) + rnf(ji,jj) - rnf_b(ji,jj) )
                     ztc_f  = ztc_f  - zfact1 * ( psbc_tc(ji,jj,jn) - psbc_tc_b(ji,jj,jn) )
                  ENDIF

                  IF( ll_traqsr .AND. jn == jp_tem .AND. jk <= nksr )   &     ! solar penetration (temperature only)
                     &     ztc_f  = ztc_f  - zfact1 * ( qsr_hc(ji,jj,jk) - qsr_hc_b(ji,jj,jk) ) 

                  IF( ll_rnf .AND. jk <= nk_rnf(ji,jj) )   & 		      ! river runoffs
                     &     ztc_f  = ztc_f  - zfact1 * ( rnf_tsc(ji,jj,jn) - rnf_tsc_b(ji,jj,jn) ) & 
                     &                              * e3t_0(ji,jj,jk) / h_rnf(ji,jj)

                  ze3t_f = 1.e0 / ze3t_f
                  ptb(ji,jj,jk,jn) = ztc_f * ze3t_f       ! ptb <-- ptn filtered
                  ptn(ji,jj,jk,jn) = pta(ji,jj,jk,jn)     ! ptn <-- pta
                  !
                  IF( ll_tra_hpg ) THEN        ! semi-implicit hpg (T & S only)
                     ze3t_d           = 1.e0   / ( ze3t_n + rbcp * ze3t_d )
                     pta(ji,jj,jk,jn) = ze3t_d * ( ztc_n  + rbcp * ztc_d  )   ! ta <-- Brown & Campana average
                  ENDIF
               END DO
            END DO
         END DO
         ! 
      END DO
      !
   END SUBROUTINE tra_nxt_vvl

   !!======================================================================
END MODULE tranxt
