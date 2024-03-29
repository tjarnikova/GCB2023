MODULE p4zmort
   !!======================================================================
   !!                         ***  MODULE p4zmort  ***
   !! TOP :   PISCES Compute the mortality terms for phytoplankton
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont)  Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_mort       :   Compute the mortality terms for phytoplankton
   !!   p4z_mort_init  :   Initialize the mortality params for phytoplankton
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zsink         !  vertical flux of particulate matter due to sinking
   USE p4zprod         !  Primary productivity 
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_mort    
   PUBLIC   p4z_mort_init    

   !! * Shared module variables
   REAL(wp), PUBLIC :: wchl    !:
   REAL(wp), PUBLIC :: wchld   !:
   REAL(wp), PUBLIC :: wchldm  !:
   REAL(wp), PUBLIC :: mprat   !:
   REAL(wp), PUBLIC :: mprat2  !:


   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zmort.F90 5385 2015-06-09 13:50:42Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_mort( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_mort  ***
      !!
      !! ** Purpose :   Calls the different subroutine to initialize and compute
      !!                the different phytoplankton mortality terms
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt ! ocean time step
      !!---------------------------------------------------------------------

      CALL p4z_nano            ! nanophytoplankton

      CALL p4z_diat            ! diatoms

   END SUBROUTINE p4z_mort


   SUBROUTINE p4z_nano
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_nano  ***
      !!
      !! ** Purpose :   Compute the mortality terms for nanophytoplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zsizerat, zcompaph
      REAL(wp) :: zfactfe, zfactch, zprcaca, zfracal
      REAL(wp) :: ztortp , zrespp , zmortp , zstep
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_nano')
      !
      prodcal(:,:,:) = 0.  !: calcite production variable set to zero
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zcompaph = MAX( ( trb(ji,jj,jk,jpphy) - 1e-8 ), 0.e0 )
               zstep    = xstep
# if defined key_degrad
               zstep    = zstep * facvol(ji,jj,jk)
# endif
               !     When highly limited by macronutrients, very small cells 
               !     dominate the community. As a consequence, aggregation
               !     due to turbulence is negligible. Mortality is also set
               !     to 0
               zsizerat = MIN(1., MAX( 0., (quotan(ji,jj,jk) - 0.2) / 0.3) ) * trb(ji,jj,jk,jpphy)
               !     Squared mortality of Phyto similar to a sedimentation term during
               !     blooms (Doney et al. 1996)
               zrespp = wchl * 1.e6 * zstep * xdiss(ji,jj,jk) * zcompaph * zsizerat 

               !     Phytoplankton mortality. This mortality loss is slightly
               !     increased when nutrients are limiting phytoplankton growth
               !     as observed for instance in case of iron limitation.
               ztortp = mprat * xstep * zcompaph / ( xkmort + trb(ji,jj,jk,jpphy) ) * zsizerat

               zmortp = zrespp + ztortp

               !   Update the arrays TRA which contains the biological sources and sinks

               zfactfe = trb(ji,jj,jk,jpnfe)/(trb(ji,jj,jk,jpphy)+rtrn)
               zfactch = trb(ji,jj,jk,jpnch)/(trb(ji,jj,jk,jpphy)+rtrn)
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zmortp
               tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) - zmortp * zfactch
               tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) - zmortp * zfactfe
               zprcaca = xfracal(ji,jj,jk) * zmortp
               !
               prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
               !
               zfracal = 0.5 * xfracal(ji,jj,jk)
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - zprcaca
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2. * zprcaca
               tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal) + zprcaca
#if defined key_kriest
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zmortp
               tra(ji,jj,jk,jpnum) = tra(ji,jj,jk,jpnum) + ztortp * xkr_dnano + zrespp * xkr_ddiat
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zmortp * zfactfe
#else
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zfracal * zmortp
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + ( 1. - zfracal ) * zmortp
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + ( 1. - zfracal ) * zmortp * zfactfe
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zfracal * zmortp * zfactfe
#endif
            END DO
         END DO
      END DO
      !
       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('nano')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_nano')
      !
   END SUBROUTINE p4z_nano

   SUBROUTINE p4z_diat
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_diat  ***
      !!
      !! ** Purpose :   Compute the mortality terms for diatoms
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER  ::  ji, jj, jk
      REAL(wp) ::  zfactfe,zfactsi,zfactch, zcompadi
      REAL(wp) ::  zrespp2, ztortp2, zmortp2, zstep
      REAL(wp) ::  zlim2, zlim1
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_diat')
      !

      !    Aggregation term for diatoms is increased in case of nutrient
      !    stress as observed in reality. The stressed cells become more
      !    sticky and coagulate to sink quickly out of the euphotic zone
      !     ------------------------------------------------------------

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               zcompadi = MAX( ( trb(ji,jj,jk,jpdia) - 1e-9), 0. )

               !    Aggregation term for diatoms is increased in case of nutrient
               !    stress as observed in reality. The stressed cells become more
               !    sticky and coagulate to sink quickly out of the euphotic zone
               !     ------------------------------------------------------------
               zstep   = xstep
# if defined key_degrad
               zstep = zstep * facvol(ji,jj,jk)
# endif
               !  Phytoplankton respiration 
               !     ------------------------
               zlim2   = xlimdia(ji,jj,jk) * xlimdia(ji,jj,jk)
               zlim1   = 0.25 * ( 1. - zlim2 ) / ( 0.25 + zlim2 ) 
               zrespp2 = 1.e6 * zstep * (  wchld + wchldm * zlim1 ) * xdiss(ji,jj,jk) * zcompadi * trb(ji,jj,jk,jpdia)

               !     Phytoplankton mortality. 
               !     ------------------------
               ztortp2 = mprat2 * zstep * trb(ji,jj,jk,jpdia)  / ( xkmort + trb(ji,jj,jk,jpdia) ) * zcompadi 

               zmortp2 = zrespp2 + ztortp2

               !   Update the arrays tra which contains the biological sources and sinks
               !   ---------------------------------------------------------------------
               zfactch = trb(ji,jj,jk,jpdch) / ( trb(ji,jj,jk,jpdia) + rtrn )
               zfactfe = trb(ji,jj,jk,jpdfe) / ( trb(ji,jj,jk,jpdia) + rtrn )
               zfactsi = trb(ji,jj,jk,jpdsi) / ( trb(ji,jj,jk,jpdia) + rtrn )
               tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - zmortp2 
               tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) - zmortp2 * zfactch
               tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) - zmortp2 * zfactfe
               tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) - zmortp2 * zfactsi
               tra(ji,jj,jk,jpgsi) = tra(ji,jj,jk,jpgsi) + zmortp2 * zfactsi
#if defined key_kriest
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zmortp2  
               tra(ji,jj,jk,jpnum) = tra(ji,jj,jk,jpnum) + ztortp2 * xkr_ddiat + zrespp2 * xkr_daggr
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zmortp2 * zfactfe
#else
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zrespp2 + 0.5 * ztortp2
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + 0.5 * ztortp2
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + 0.5 * ztortp2 * zfactfe
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + ( zrespp2 + 0.5 * ztortp2 ) * zfactfe
#endif
            END DO
         END DO
      END DO
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('diat')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_diat')
      !
   END SUBROUTINE p4z_diat

   SUBROUTINE p4z_mort_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_mort_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton parameters
      !!
      !! ** Method  :   Read the nampismort namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampismort
      !!
      !!----------------------------------------------------------------------

      NAMELIST/nampismort/ wchl, wchld, wchldm, mprat, mprat2
      INTEGER :: ios                 ! Local integer output status for namelist read

      REWIND( numnatp_ref )              ! Namelist nampismort in reference namelist : Pisces phytoplankton
      READ  ( numnatp_ref, nampismort, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismort in reference namelist', lwp )

      REWIND( numnatp_cfg )              ! Namelist nampismort in configuration namelist : Pisces phytoplankton
      READ  ( numnatp_cfg, nampismort, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismort in configuration namelist', lwp )
      IF(lwm) WRITE ( numonp, nampismort )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for phytoplankton mortality, nampismort'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    quadratic mortality of phytoplankton      wchl      =', wchl
         WRITE(numout,*) '    maximum quadratic mortality of diatoms    wchld     =', wchld
         WRITE(numout,*) '    maximum quadratic mortality of diatoms    wchldm    =', wchldm
         WRITE(numout,*) '    phytoplankton mortality rate              mprat     =', mprat
         WRITE(numout,*) '    Diatoms mortality rate                    mprat2    =', mprat2
      ENDIF

   END SUBROUTINE p4z_mort_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_mort                    ! Empty routine
   END SUBROUTINE p4z_mort
#endif 

   !!======================================================================
END MODULE  p4zmort
