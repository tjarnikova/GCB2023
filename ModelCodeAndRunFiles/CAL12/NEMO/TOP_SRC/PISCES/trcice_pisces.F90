MODULE trcice_pisces
   !!======================================================================
   !!                         ***  MODULE trcice_pisces  ***
   !! TOP :   initialisation of the PISCES biochemical model
   !!======================================================================
   !! History :  3.5  ! 2013    (M. Vancoppenolle, O. Aumont, G. Madec), original code
   !! Comment ! probably not properly done when the second particle export
   !! scheme (kriest) is used
   !!----------------------------------------------------------------------
#if defined key_pisces || defined key_pisces_reduced
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !! trc_ice_pisces   : PISCES fake sea ice model setting
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE par_pisces      ! PISCES parameters
   USE oce_trc         ! Shared variables between ocean and passive tracers
   USE trc             ! Passive tracers common variables 
   USE phycst          ! Ocean physics parameters
   USE sms_pisces      ! PISCES Source Minus Sink variables
   USE in_out_manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ice_ini_pisces ! called by trcini.F90 module

CONTAINS

   SUBROUTINE trc_ice_ini_pisces
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_ice_ini_pisces ***
      !!
      !! ** Purpose :   PISCES fake sea ice model setting
      !!    Method  :   Assign prescribe values to tracer concentrations in sea ice
      !!
      !! For levitating sea ice, constant ocean tracer concentrations also have to be defined. 
      !! This is done specifically for Global, Arctic, Antarctic and Baltic regions
      !! 
      !! Sea ice concentrations are by default prescribed as follows
      !!  trc_i = zratio * trc_o
      !!
      !! This formulation is modulated by the namelist parameter trc_ice_ratio
      !!
      !! trc_ice_ratio  * betw 0 and 1: prescribed ice/ocean tracer concentration ratio
      !!                * -1 => the ice-ocean tracer concentration ratio follows the
      !!                         ice-ocean salinity ratio
      !!                * -2 => no ice-ocean tracer concentration is used
      !!                        instead, the tracer concentration in sea ice 
      !!                        is prescribed to trc_ice_prescr
      !! 
      !! cn_trc_o  specifies which disinctions are made for prescribed tracer concentration
      !!                * 'GL' use global ocean values making distinction for Baltic Sea only
      !!                * 'AA' use Arctic/Antarctic contrasted values, + Baltic
      !!
      !!----------------------------------------------------------------------

                                        !--- Dummy variables
      REAL(wp), DIMENSION(jp_pisces,2)  :: zratio  ! effective ice-ocean tracer cc ratio
      REAL(wp), DIMENSION(jp_pisces,4)  :: zpisc   ! prescribes concentration 
      !                                            !  1:global, 2:Arctic, 3:Antarctic, 4:Baltic

      REAL(wp), DIMENSION(2) :: zrs  ! ice-ocean salinity ratio, 1 - global, 2- Baltic
      REAL(wp) :: zsice_bal          ! prescribed ice salinity in the Baltic
      REAL(wp) :: zsoce_bal          ! prescribed ocean salinity in the Baltic
      REAL(wp) :: zfeoce_glo         ! prescribed iron concentration in the global ocean
      REAL(wp) :: zfeoce_bal         ! prescribed iron concentration in the global ocean
      INTEGER  :: jn                 ! dummy loop index

      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ice_ini_pisces: Prescribed sea ice biogeochemistry '
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~'

      !--------------------------------------------
      ! Initialize ocean prescribed concentrations
      !--------------------------------------------
      ! values taken from a 500 yr equilibrium run
      ! used only in the levitating sea ice case with virtual salt / tracer
      ! fluxes

      !--- Global values
      zpisc(jpdic,1) =  1.99e-3_wp 
      zpisc(jpdoc,1) =  2.04e-5_wp 
      zpisc(jptal,1) =  2.31e-3_wp 
      zpisc(jpoxy,1) =  2.47e-4_wp
      zpisc(jpcal,1) =  1.04e-8_wp
      zpisc(jppo4,1) =  5.77e-7_wp / po4r 
      zpisc(jppoc,1) =  1.27e-6_wp  
#  if ! defined key_kriest
      zpisc(jpgoc,1) =  5.23e-8_wp  
      zpisc(jpbfe,1) =  9.84e-13_wp 
#  else
      zpisc(jpnum,1) = 0. ! could not get this value since did not use it
#  endif
      zpisc(jpsil,1) =  7.36e-6_wp  
      zpisc(jpdsi,1) =  1.07e-7_wp 
      zpisc(jpgsi,1) =  1.53e-8_wp
      zpisc(jpphy,1) =  9.57e-8_wp
      zpisc(jpdia,1) =  4.24e-7_wp
      zpisc(jpzoo,1) =  6.07e-7_wp
      zpisc(jpmes,1) =  3.44e-7_wp
      zpisc(jpfer,1) =  4.06e-10_wp
      zpisc(jpsfe,1) =  2.51e-11_wp
      zpisc(jpdfe,1) =  6.57e-12_wp
      zpisc(jpnfe,1) =  1.76e-11_wp
      zpisc(jpnch,1) =  1.67e-7_wp
      zpisc(jpdch,1) =  1.02e-7_wp
      zpisc(jpno3,1) =  5.79e-6_wp / rno3 
      zpisc(jpnh4,1) =  3.22e-7_wp / rno3

      !--- Arctic specificities (dissolved inorganic & DOM)
      zpisc(jpdic,2) =  1.98e-3_wp 
      zpisc(jpdoc,2) =  6.00e-6_wp 
      zpisc(jptal,2) =  2.13e-3_wp 
      zpisc(jpoxy,2) =  3.65e-4_wp  
      zpisc(jpcal,2) =  1.50e-9_wp  
      zpisc(jppo4,2) =  4.09e-7_wp / po4r 
      zpisc(jppoc,2) =  4.05e-7_wp  
#  if ! defined key_kriest
      zpisc(jpgoc,2) =  2.84e-8_wp  
      zpisc(jpbfe,2) =  7.03e-13_wp 
#  else
      zpisc(jpnum,2) =  0.00e-00_wp 
#  endif
      zpisc(jpsil,2) =  6.87e-6_wp  
      zpisc(jpdsi,2) =  1.73e-7_wp 
      zpisc(jpgsi,2) =  7.93e-9_wp
      zpisc(jpphy,2) =  5.25e-7_wp  
      zpisc(jpdia,2) =  7.75e-7_wp 
      zpisc(jpzoo,2) =  3.34e-7_wp
      zpisc(jpmes,2) =  2.49e-7_wp  
      zpisc(jpfer,2) =  1.43e-9_wp 
      zpisc(jpsfe,2) =  2.21e-11_wp 
      zpisc(jpdfe,2) =  2.04e-11_wp 
      zpisc(jpnfe,2) =  1.75e-11_wp 
      zpisc(jpnch,2) =  1.46e-07_wp 
      zpisc(jpdch,2) =  2.36e-07_wp 
      zpisc(jpno3,2) =  3.51e-06_wp / rno3 
      zpisc(jpnh4,2) =  6.15e-08_wp / rno3 

      !--- Antarctic specificities (dissolved inorganic & DOM)
      zpisc(jpdic,3) =  2.20e-3_wp  
      zpisc(jpdoc,3) =  7.02e-6_wp  
      zpisc(jptal,3) =  2.37e-3_wp  
      zpisc(jpoxy,3) =  3.42e-4_wp  
      zpisc(jpcal,3) =  3.17e-9_wp  
      zpisc(jppo4,3) =  1.88e-6_wp / po4r  
      zpisc(jppoc,3) =  1.13e-6_wp  
#  if ! defined key_kriest
      zpisc(jpgoc,3) =  2.89e-8_wp  
      zpisc(jpbfe,3) =  5.63e-13_wp 
#  else
      zpisc(jpnum,3) =  0.00e-00_wp 
#  endif
      zpisc(jpsil,3) =  4.96e-5_wp  
      zpisc(jpdsi,3) =  5.63e-7_wp 
      zpisc(jpgsi,3) =  5.35e-8_wp
      zpisc(jpphy,3) =  8.10e-7_wp  
      zpisc(jpdia,3) =  5.77e-7_wp 
      zpisc(jpzoo,3) =  6.68e-7_wp
      zpisc(jpmes,3) =  3.55e-7_wp  
      zpisc(jpfer,3) =  1.62e-10_wp
      zpisc(jpsfe,3) =  2.29e-11_wp 
      zpisc(jpdfe,3) =  8.75e-12_wp
      zpisc(jpnfe,3) =  1.48e-11_wp 
      zpisc(jpnch,3) =  2.02e-7_wp  
      zpisc(jpdch,3) =  1.60e-7_wp  
      zpisc(jpno3,3) =  2.64e-5_wp / rno3  
      zpisc(jpnh4,3) =  3.39e-7_wp / rno3  

      !--- Baltic Sea particular case for ORCA configurations
      zpisc(jpdic,4) = 1.14e-3_wp
      zpisc(jpdoc,4) = 1.06e-5_wp
      zpisc(jptal,4) = 1.16e-3_wp
      zpisc(jpoxy,4) = 3.71e-4_wp
      zpisc(jpcal,4) = 1.51e-9_wp
      zpisc(jppo4,4) = 2.85e-9_wp / po4r
      zpisc(jppoc,4) = 4.84e-7_wp
#  if ! defined key_kriest
      zpisc(jpgoc,4) = 1.05e-8_wp
      zpisc(jpbfe,4) = 4.97e-13_wp
#  else
      zpisc(jpnum,4) = 0. ! could not get this value
#  endif
      zpisc(jpsil,4) = 4.91e-5_wp
      zpisc(jpdsi,4) = 3.25e-7_wp
      zpisc(jpgsi,4) = 1.93e-8_wp
      zpisc(jpphy,4) = 6.64e-7_wp
      zpisc(jpdia,4) = 3.41e-7_wp
      zpisc(jpzoo,4) = 3.83e-7_wp
      zpisc(jpmes,4) = 0.225e-6_wp
      zpisc(jpfer,4) = 2.45e-9_wp
      zpisc(jpsfe,4) = 3.89e-11_wp
      zpisc(jpdfe,4) = 1.33e-11_wp
      zpisc(jpnfe,4) = 2.62e-11_wp
      zpisc(jpnch,4) = 1.17e-7_wp
      zpisc(jpdch,4) = 9.69e-8_wp
      zpisc(jpno3,4) = 5.36e-5_wp / rno3
      zpisc(jpnh4,4) = 7.18e-7_wp / rno3
 
      DO jn = jp_pcs0, jp_pcs1
         IF( cn_trc_o(jn) == 'GL ' ) trc_o(:,:,jn) = zpisc(jn,1)  ! Global case
         IF( cn_trc_o(jn) == 'AA ' ) THEN 
            WHERE( gphit(:,:) >= 0._wp ) ; trc_o(:,:,jn) = zpisc(jn,2) ; END WHERE ! Arctic 
            WHERE( gphit(:,:) <  0._wp ) ; trc_o(:,:,jn) = zpisc(jn,3) ; END WHERE ! Antarctic 
         ENDIF
         IF( cp_cfg == "orca" ) THEN     !  Baltic Sea particular case for ORCA configurations
             WHERE( 14._wp <= glamt(:,:) .AND. glamt(:,:) <= 32._wp .AND.    &
                    54._wp <= gphit(:,:) .AND. gphit(:,:) <= 66._wp )
                    trc_o(:,:,jn) = zpisc(jn,4)
            END WHERE
         ENDIF 
      ENDDO



      !-----------------------------
      ! Assign ice-ocean cc ratios 
      !-----------------------------
      ! 0 means zero concentration in sea ice
      ! 1 means same concentration in the sea ice as in the ocean

      ! Ice ocean salinity ratio
      zsoce_bal   = 4. ; zsice_bal   = 2. !! Baltic ocean and sea ice salinities
      zrs(1) = sice / soce                !! ice-ocean salinity ratio, global case
      zrs(2) = zsice_bal / zsoce_bal      !! ice-ocean salinity ratio, Baltic case

      DO jn = jp_pcs0, jp_pcs1
         IF( trc_ice_ratio(jn) >= 0._wp )  zratio(jn,:) = trc_ice_ratio(jn)
         IF( trc_ice_ratio(jn) == -1._wp ) zratio(jn,:) = zrs(:)
         IF( trc_ice_ratio(jn) == -2._wp ) zratio(jn,:) = -9999.99_wp
      END DO

      !-------------------------------
      ! Sea ice tracer concentrations
      !-------------------------------
      DO jn = jp_pcs0, jp_pcs1
         !-- Everywhere but in the Baltic
         IF ( trc_ice_ratio(jn) >= -1._wp ) THEN ! no prescribed conc. ; typically everything but iron) 
            trc_i(:,:,jn) = zratio(jn,1) * trc_o(:,:,jn) 
         ELSE                                    ! prescribed concentration
            trc_i(:,:,jn) = trc_ice_prescr(jn)
         ENDIF
       
         !-- Baltic
         IF( cp_cfg == "orca" ) THEN  ! Baltic treated seperately for ORCA configs
            IF ( trc_ice_ratio(jn) >= - 1._wp ) THEN ! no prescribed conc. ; typically everything but iron) 
               WHERE( 14._wp <= glamt(:,:) .AND. glamt(:,:) <= 32._wp .AND.    &
                      54._wp <= gphit(:,:) .AND. gphit(:,:) <= 66._wp )
                     trc_i(:,:,jn) = zratio(jn,2) * trc_o(:,:,jn) 
               END WHERE
            ELSE                                 ! prescribed tracer concentration in ice
               WHERE( 14._wp <= glamt(:,:) .AND. glamt(:,:) <= 32._wp .AND.    &
                   54._wp <= gphit(:,:) .AND. gphit(:,:) <= 66._wp )
                     trc_i(:,:,jn) = trc_ice_prescr(jn)
               END WHERE
            ENDIF ! trc_ice_ratio
         ENDIF
      !
      END DO ! jn

   END SUBROUTINE trc_ice_ini_pisces

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                            No PISCES biochemical model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ice_ini_pisces         ! Empty routine
   END SUBROUTINE trc_ice_ini_pisces
#endif

   !!======================================================================
END MODULE trcice_pisces
