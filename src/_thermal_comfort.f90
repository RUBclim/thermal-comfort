    !~ UTCI, Version a 0.002, October 2009
    !~ Copyright (C) 2009  Peter Broede

    !~ Program for calculating UTCI Temperature (UTCI)
    !~ released for public use after termination of COST Action 730

    !~ replaces Version a 0.001, from September 2009

    !~ Version History:

    !~ pre-alpha version  October 2008:
    !~ based on calculations as presented at the MC/WG meeting in Eilat, 2008-09-10

    !~ Version a 0.001, September 2009:
    !~ Changes compared to the pre-alpha version were based on decisions after the
    !~ WG1 meeting in Stuttgart 2009-02-24:
     !~ - Changed clothing model with respect to the maximum insulation and
       !~ to reduction of clothing insulation by wind speed
     !~ - Changed physiological model with respect to hidromeiosis
     !~ - Equivalent temperature calculated by explicit comparisons
       !~ to values of a 'Response Index' computed
       !~ as first principal component from values after 30 and 120 min of
       !~ rectal temperature, mean skin temperature, face skin temperature,
       !~ sweat rate, wetted skin area, shivering, skin blood flow.
     !~ - as in the pre-alpha version, UTCI values are approximated by a 6th order polynomial regression function
     !~ - The limits of mean radiant temperatures were extended to 30 degC below
       !~ and 70 degC above air temperature
      !~ Please note: The given polynomial approximation limits the application
      !~ of this procedure to values of wind speed between 0.5 and 17 m/s!

    !~ Version a 0.002, October 2009:
    !~ Changed ReadMe text and program messages for public release

    !~ Copyright (C) 2009  Peter Broede

    !~ This program is distributed in the hope that it will be useful,
    !~ but WITHOUT ANY WARRANTY; without even the implied warranty of
    !~ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

!~ Disclaimer of Warranty.

!~ THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN OTHERWISE STATED IN WRITING
!~ THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM �AS IS� WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED
!~ OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
!~ THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE DEFECTIVE, YOU
!~ ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.


!~ Limitation of Liability.

!~ IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO
!~ MODIFIES AND/OR CONVEYS THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL,
!~ INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM
!~ (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES
!~ OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED
!~ OF THE POSSIBILITY OF SUCH DAMAGES.

module thermal_comfort_mod
    IMPLICIT NONE

    PRIVATE
!
!   Internal constants:
! ----------------
    REAL(kind=8), PARAMETER :: po   = 1013.25   !< air pressure at sea level (hPa)
    REAL(kind=8), PARAMETER :: rob  =    1.06   !<
    REAL(kind=8), PARAMETER :: cb   = 3.64 * 1000.  !<
    REAL(kind=8), PARAMETER :: food =    0.         !< heat gain by food        (W)
    REAL(kind=8), PARAMETER :: emsk =    0.99   !< longwave emission coef. of skin
    REAL(kind=8), PARAMETER :: emcl =    0.95   !< longwave emission coef. of cloth
    REAL(kind=8), PARAMETER :: evap = 2.42 * 10. ** 6.  !<
    REAL(kind=8), PARAMETER :: sigm = 5.67 * 10. **(-8.)  !< Stefan-Boltzmann-Const.
    REAL(kind=8), PARAMETER :: cair = 1.01 * 1000.      !<
!
!   Internal variables:
! ----------------
    REAL(kind=8) :: h           !< internal heat        (W)

!   MEMI configuration
    REAL(kind=8) :: age         !< persons age          (a)
    REAL(kind=8) :: mbody       !< persons body mass    (kg)
    REAL(kind=8) :: ht          !< persons height       (m)
    REAL(kind=8) :: work        !< work load            (W)
    REAL(kind=8) :: eta         !< work efficiency      (dimensionless)
    REAL(kind=8) :: icl         !< clothing insulation index (clo)
    REAL(kind=8) :: fcl         !< surface area modification by clothing (factor)
    INTEGER :: pos     !< posture: 1 = standing, 2 = sitting
    INTEGER :: sex     !< sex: 1 = male, 2 = female

    PUBLIC calculate_pet_static
    PUBLIC UTCI_approx
    PUBLIC UTCI_approx_vectorized
    !~ TODO: they should likely be private
    PUBLIC es
    PUBLIC es_vectorized

  CONTAINS

!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Physiologically Equivalent Temperature (PET),
!  stationary (calculated based on MEMI),
!  Subroutine based on PETBER vers. 1.5.1996 by P. Hoeppe
!------------------------------------------------------------------------------!
  SUBROUTINE calculate_pet_static( ta, vpa, v, tmrt, p, tx )
    !-- Configure sample person (optional)
    !  age, mbody, ht, work, eta, icl, fcl, pos, sex )

       IMPLICIT NONE
    !
    !  Input arguments:
    ! ----------------
       REAL(kind=8), INTENT( IN ) :: ta     !< Air temperature 		(°C)	REAL(kind=8)
       REAL(kind=8), INTENT( IN ) :: tmrt   !< Mean radiant temperature 	(°C)	REAL(kind=8)
       REAL(kind=8), INTENT( IN ) :: v      !< Wind speed 			(m/s)	REAL(kind=8)
       REAL(kind=8), INTENT( IN ) :: vpa    !< Vapor pressure 		(hPa)	REAL(kind=8)
       REAL(kind=8), INTENT( IN ) :: p      !< Air pressure                  (hPa)   REAL(kind=8)
    !
    !  Output arguments:
    ! ----------------
       REAL(kind=8), INTENT ( OUT ) :: tx   !< PET 				(°C)	REAL(kind=8)
    !  former intent (out), disabled:
       !  - tsk        : Skin temperature              (°C)    real
       !  - tcl        : Clothing temperature          (°C)    real
       !  - ws         :                                       real
       !  - wetsk      : Fraction of wet skin                  real

    !
    !  Internal variables:
    ! ----------------
       REAL(kind=8) :: acl, adu, aeff, ere, erel, esw, facl, feff, rdcl,      &
              rdsk, rtv, vpts, tsk, tcl, wetsk
    !
    !  Optional arguments not supported, removed
    !  REAL(kind=8), INTENT ( in ), optional :: age, mbody, ht, work, eta, icl, fcl
    !  INTEGER, INTENT ( in ), optional :: pos, sex

    !--      Person data
    !   IF ( .NOT. present( age ) )   age   = 35.
    !   IF ( .NOT. present( mbody ) ) mbody = 75.
    !   IF ( .NOT. present( ht ) )    ht    =  1.75
    !   IF ( .NOT. present( work ) )  work  = 80.
    !   IF ( .NOT. present( eta ) )   eta   =  0.
    !   IF ( .NOT. present( icl ) )   icl   =  0.9
    !   IF ( .NOT. present( fcl ) )   fcl   =  1.15
    !   IF ( .NOT. present( pos ) )   pos   =  1
    !   IF ( .NOT. present( sex ) )   sex   =  1

    !   MEMI configuration
        age   = 35.
        mbody = 75.
        ht    =  1.75
        work  = 80.
        eta   =  0.
        icl   =  0.9
        fcl   =  1.15
        pos   =  1
        sex   =  1


    ! !-- call subfunctions
       CALL in_body ( ere, erel, p, rtv, ta, vpa )

       CALL heat_exch ( acl, adu, aeff, ere, erel, esw, facl, feff,       &
                p, rdcl, rdsk, ta, tcl, tmrt, tsk, v, vpa, vpts, wetsk )

       CALL pet_iteration ( acl, adu, aeff, esw, facl, feff, p, rdcl,     &
                rdsk, rtv, ta, tcl, tsk, tx, vpts, wetsk )


    END SUBROUTINE calculate_pet_static

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate internal energy ballance
!------------------------------------------------------------------------------!
    SUBROUTINE in_body ( ere, erel, p, rtv, ta, vpa )
      !
      !  Input arguments:
      ! ----------------
         REAL(kind=8), INTENT( IN )  :: p         !< air pressure     (hPa)
         REAL(kind=8), INTENT( IN )  :: ta        !< air temperature  (°C)
         REAL(kind=8), INTENT( IN )  :: vpa       !< vapor pressure   (hPa)
      !
      !  Output arguments:
      ! ----------------
         REAL(kind=8), INTENT( OUT ) :: ere       !< energy ballance          (W)
         REAL(kind=8), INTENT( OUT ) :: erel      !< latent energy ballance   (W)
         REAL(kind=8), INTENT( OUT ) :: rtv       !<
      !
      !  Internal variables:
      ! ----------------
         REAL(kind=8) :: eres, met, tex, vpex

         met = 3.45 * mbody ** ( 3. / 4. ) * (1. + 0.004 *                           &
               ( 30. - age) + 0.010 * ( ( ht * 100. /                                &
               ( mbody ** ( 1. / 3. ) ) ) - 43.4 ) )
      !    IF ( sex .EQ. 2 ) THEN
      !      met = 3.19 * mbody ** ( 3. / 4. ) * ( 1. + 0.004 *                        &
      !            ( 30. - age ) + 0.018 * ( ( ht * 100. / ( mbody **                  &
      !            ( 1. / 3. ) ) ) - 42.1 ) )
      !   END IF
         met = work + met

         h = met * (1. - eta)

      !-- SENSIBLE RESPIRATION ENERGY

         tex  = 0.47 * ta + 21.0
         rtv  = 1.44 * 10. ** (-6.) * met
         eres = cair * (ta - tex) * rtv

      !-- LATENT RESPIRATION ENERGY

         vpex = 6.11 * 10. ** ( 7.45 * tex / ( 235. + tex ) )
         erel = 0.623 * evap / p * ( vpa - vpex ) * rtv

      !-- SUM OF RESULTS

         ere = eres + erel

         RETURN
       END SUBROUTINE in_body


      !------------------------------------------------------------------------------!
      ! Description:
      ! ------------
      !> Calculate heat gain or loss
      !------------------------------------------------------------------------------!
       SUBROUTINE heat_exch ( acl, adu, aeff, ere, erel, esw, facl, feff,   &
              p, rdcl, rdsk, ta, tcl, tmrt, tsk, v, vpa, vpts, wetsk )

      !
      !  Input arguments:
      ! ----------------
         REAL(kind=8), INTENT( IN )  :: ere       !< energy ballance          (W)
         REAL(kind=8), INTENT( IN )  :: erel      !< latent energy ballance   (W)
         REAL(kind=8), INTENT( IN )  :: p         !< air pressure             (hPa)
         REAL(kind=8), INTENT( IN )  :: ta        !< air temperature          (°C)
         REAL(kind=8), INTENT( IN )  :: tmrt      !< mean radiant temperature (°C)
         REAL(kind=8), INTENT( IN )  :: v         !< Wind speed               (m/s)
         REAL(kind=8), INTENT( IN )  :: vpa       !< vapor pressure           (hPa)
      !
      !  Output arguments:
      ! ----------------
         REAL(kind=8), INTENT( OUT ) :: acl       !< clothing surface area        (m²)
         REAL(kind=8), INTENT( OUT ) :: adu       !< Du-Bois area                 (m²)
         REAL(kind=8), INTENT( OUT ) :: aeff      !< effective surface area       (m²)
         REAL(kind=8), INTENT( OUT ) :: esw       !< energy-loss through sweat evap. (W)
         REAL(kind=8), INTENT( OUT ) :: facl      !< surface area extension through clothing (factor)
         REAL(kind=8), INTENT( OUT ) :: feff      !< surface modification by posture (factor)
         REAL(kind=8), INTENT( OUT ) :: rdcl      !<
         REAL(kind=8), INTENT( OUT ) :: rdsk      !<
         REAL(kind=8), INTENT( OUT ) :: tcl       !< clothing temperature         (°C)
         REAL(kind=8), INTENT( OUT ) :: tsk       !< skin temperature             (°C)
         REAL(kind=8), INTENT( OUT ) :: vpts      !<
         REAL(kind=8), INTENT( OUT ) :: wetsk     !< fraction of wet skin (dimensionless)


         REAL(kind=8) :: c(0:10), cbare, cclo, csum, di, ed, enbal, enbal2, eswdif,      &
              eswphy, eswpot, fec, hc, he, htcl, r1, r2, rbare, rcl, rclo, rclo2,    &
              rsum, sw, swf, swm, tbody, tcore(1:7), vb, vb1, vb2, wd, wr, ws, wsum, &
              xx, y

         INTEGER :: count1, count3, j
         logical :: skipIncreaseCount

         wetsk = 0.
         adu = 0.203 * mbody ** 0.425 * ht ** 0.725

         hc = 2.67 + ( 6.5 * v ** 0.67)
         hc   = hc * (p /po) ** 0.55
         feff = 0.725                     !< Posture: 0.725 for stading
         IF ( pos .EQ. 2 ) feff = 0.696   !<          0.696 for sitting
         facl = (- 2.36 + 173.51 * icl - 100.76 * icl * icl + 19.28                  &
               * (icl ** 3.)) / 100.

         IF ( facl .GT. 1. )   facl = 1.
         rcl = ( icl / 6.45) / facl
         IF ( icl .GE. 2. )  y  = 1.

         IF ( ( icl .GT. 0.6 ) .AND. ( icl .LT. 2. ) )  y = ( ht - 0.2 ) / ht
         IF ( ( icl .LE. 0.6 ) .AND. ( icl .GT. 0.3 ) ) y = 0.5
         IF ( ( icl .LE. 0.3 ) .AND. ( icl .GT. 0. ) )  y = 0.1

         r2   = adu * (fcl - 1. + facl) / (2. * 3.14 * ht * y)
         r1   = facl * adu / (2. * 3.14 * ht * y)

         di   = r2 - r1

      !-- SKIN TEMPERATURE

         DO j = 1,7

           tsk    = 34.
           count1 = 0
           tcl    = ( ta + tmrt + tsk ) / 3.
           count3 = 1
           enbal2 = 0.

           DO
             acl   = adu * facl + adu * ( fcl - 1. )
             rclo2 = emcl * sigm * ( ( tcl + 273.2 )**4. -                           &
               ( tmrt + 273.2 )** 4. ) * feff
             htcl  = 6.28 * ht * y * di / ( rcl * LOG( r2 / r1 ) * acl )
             tsk   = 1. / htcl * ( hc * ( tcl - ta ) + rclo2 ) + tcl

      !--    RADIATION SALDO

             aeff  = adu * feff
             rbare = aeff * ( 1. - facl ) * emsk * sigm *                            &
               ( ( tmrt + 273.2 )** 4. - ( tsk + 273.2 )** 4. )
             rclo  = feff * acl * emcl * sigm *                                      &
               ( ( tmrt + 273.2 )** 4. - ( tcl + 273.2 )** 4. )
             rsum  = rbare + rclo

      !--    CONVECTION

             cbare = hc * ( ta - tsk ) * adu * ( 1. - facl )
             cclo  = hc * ( ta - tcl ) * acl
             csum  = cbare + cclo

       !--   CORE TEMPERATUR

             c(0)  = h + ere
             c(1)  = adu * rob * cb
             c(2)  = 18. - 0.5 * tsk
             c(3)  = 5.28 * adu * c(2)
             c(4)  = 0.0208 * c(1)
             c(5)  = 0.76075 * c(1)
             c(6)  = c(3) - c(5) - tsk * c(4)
             c(7)  = - c(0) * c(2) - tsk * c(3) + tsk * c(5)
             c(8)  = c(6) * c(6) - 4. * c(4) * c(7)
             c(9)  = 5.28 * adu - c(5) - c(4) * tsk
             c(10) = c(9) * c(9) - 4. * c(4) *                                       &
               ( c(5) * tsk - c(0) - 5.28 * adu * tsk )

             IF ( tsk .EQ. 36. ) tsk = 36.01
             tcore(7) = c(0) / ( 5.28 * adu + c(1) * 6.3 / 3600. ) + tsk
             tcore(3) = c(0) / ( 5.28 * adu + ( c(1) * 6.3 / 3600. ) /               &
               ( 1. + 0.5 * ( 34. - tsk ) ) ) + tsk
             IF ( c(10) .GE. 0.) THEN
               tcore(6) = ( - c(9) - c(10)**0.5 ) / ( 2. * c(4) )
               tcore(1) = ( - c(9) + c(10)**0.5 ) / ( 2. * c(4) )
             END IF

             IF ( c(8) .GE. 0. ) THEN
               tcore(2) = ( - c(6) + ABS( c(8) ) ** 0.5 ) / ( 2. * c(4) )
               tcore(5) = ( - c(6) - ABS( c(8) ) ** 0.5 ) / ( 2. * c(4) )
               tcore(4) = c(0) / ( 5.28 * adu + c(1) * 1. / 40. ) + tsk
             END IF

      !--    TRANSPIRATION

             tbody = 0.1 * tsk + 0.9 * tcore(j)
             swm   = 304.94 * ( tbody - 36.6 ) * adu / 3600000.
             vpts  = 6.11 * 10.**( 7.45 * tsk / ( 235. + tsk ) )

             IF ( tbody .LE. 36.6 ) swm = 0.

             sw = swm                                 !< make sure sw is initialized
             IF ( sex .EQ. 2 ) sw = 0.7 * swm         !< reduce if female
             eswphy = - sw * evap
             he     = 0.633 * hc / ( p * cair )
             fec    = 1. / ( 1. + 0.92 * hc * rcl )
             eswpot = he * ( vpa - vpts ) * adu * evap * fec
             wetsk  = eswphy / eswpot

             IF ( wetsk .GT. 1. ) wetsk = 1.

             eswdif = eswphy - eswpot

             IF ( eswdif .LE. 0. ) esw = eswpot
             IF ( eswdif .GT. 0. ) esw = eswphy
             IF ( esw  .GT. 0. )   esw = 0.

      !--    DIFFUSION

             rdsk = 0.79 * 10. ** 7.
             rdcl = 0.
             ed   = evap / ( rdsk + rdcl ) * adu * ( 1. - wetsk ) * ( vpa - vpts )

      !--    MAX VB

             vb1 = 34. - tsk
             vb2 = tcore(j) - 36.6

             IF ( vb2 .LT. 0. ) vb2 = 0.
             IF ( vb1 .LT. 0. ) vb1 = 0.
             vb = ( 6.3 + 75. * vb2 ) / ( 1. + 0.5 * vb1 )

      !--    ENERGY BALLANCE

             enbal = h + ed + ere + esw + csum + rsum + food


      !--    CLOTHING TEMPERATURE

             xx = 0.001
             IF ( count1 .EQ. 0 ) xx = 1.
             IF ( count1 .EQ. 1 ) xx = 0.1
             IF ( count1 .EQ. 2 ) xx = 0.01
             IF ( count1 .EQ. 3 ) xx = 0.001

             IF ( enbal .GT. 0. ) tcl = tcl + xx
             IF ( enbal .LT. 0. ) tcl = tcl - xx

             skipIncreaseCount = .FALSE.
             IF ( ( (enbal .LE. 0.) .AND. (enbal2 .GT. 0.) ) .OR.                    &
                ( ( enbal .GE. 0. ) .AND. ( enbal2 .LT. 0. ) ) ) THEN
               skipIncreaseCount = .TRUE.
             ELSE
               enbal2 = enbal
               count3 = count3 + 1
             END IF

             IF ( ( count3 .GT. 200 ) .OR. skipIncreaseCount ) THEN
               IF ( count1 .LT. 3 ) THEN
                 count1 = count1 + 1
                 enbal2 = 0.
               ELSE
                 EXIT
               END IF
             END IF
           END DO

           IF ( count1 .EQ. 3 ) THEN
             SELECT CASE ( j )
               CASE ( 2, 5)
                 IF ( .NOT. ( ( tcore(j) .GE. 36.6 ) .AND.                           &
                   ( tsk .LE. 34.050 ) ) ) CYCLE
               CASE ( 6, 1 )
                 IF ( c(10) .LT. 0. ) CYCLE
                 IF ( .NOT. ( ( tcore(j) .GE. 36.6 ) .AND.                           &
                   ( tsk .GT. 33.850 ) ) ) CYCLE
               CASE ( 3 )
                 IF ( .NOT. ( ( tcore(j) .LT. 36.6 ) .AND.                           &
                   ( tsk .LE. 34.000 ) ) ) CYCLE
               CASE ( 7 )
                 IF ( .NOT. ( ( tcore(j) .LT. 36.6 ) .AND.                           &
                   ( tsk .GT. 34.000 ) ) ) CYCLE
               CASE default  !< same as CASE ( 4 ), does actually nothing
             END SELECT
           END IF

           IF ( ( j .NE. 4 ) .AND. ( vb .GE. 91. ) ) CYCLE
           IF ( ( j .EQ. 4 ) .AND. ( vb .LT. 89. ) ) CYCLE
           IF ( vb .GT. 90.) vb = 90.

      !--  LOSSES BY WATER

           ws = sw * 3600. * 1000.
           IF ( ws .GT. 2000. ) ws = 2000.
           wd = ed / evap * 3600. * ( -1000. )
           wr = erel / evap * 3600. * ( -1000. )

           wsum = ws + wr + wd

           RETURN
         END DO
         RETURN
       END SUBROUTINE heat_exch



      !------------------------------------------------------------------------------!
      ! Description:
      ! ------------
      !> Calculate PET
      !------------------------------------------------------------------------------!
       SUBROUTINE pet_iteration ( acl, adu, aeff, esw, facl, feff, p, rdcl,          &
              rdsk, rtv, ta, tcl, tsk, tx, vpts, wetsk )
      !
      !  Input arguments:
      ! ----------------
         REAL(kind=8), INTENT( IN )  :: acl       !< clothing surface area        (m²)
         REAL(kind=8), INTENT( IN )  :: adu       !< Du-Bois area                 (m²)
         REAL(kind=8), INTENT( IN )  :: esw       !< energy-loss through sweat evap. (W)
         REAL(kind=8), INTENT( IN )  :: facl      !< surface area extension through clothing (factor)
         REAL(kind=8), INTENT( IN )  :: feff      !< surface modification by posture (factor)
         REAL(kind=8), INTENT( IN )  :: p         !< air pressure                 (hPa)
         REAL(kind=8), INTENT( IN )  :: rdcl      !<
         REAL(kind=8), INTENT( IN )  :: rdsk      !<
         REAL(kind=8), INTENT( IN )  :: rtv       !<
         REAL(kind=8), INTENT( IN )  :: ta        !< air temperature              (°C)
         REAL(kind=8), INTENT( IN )  :: tcl       !< clothing temperature         (°C)
         REAL(kind=8), INTENT( IN )  :: tsk       !< skin temperature             (°C)
         REAL(kind=8), INTENT( IN )  :: vpts      !<
         REAL(kind=8), INTENT( IN )  :: wetsk     !< fraction of wet skin (dimensionless)
      !
      !  Output arguments:
      ! ----------------
         REAL(kind=8), INTENT( OUT ) :: aeff      !< effective surface area       (m²)
         REAL(kind=8), INTENT( OUT ) :: tx        !< PET                          (°C)



         REAL (kind=8)  :: cbare, cclo, csum, ed, enbal, enbal2, ere, erel, eres, hc,  &
              rbare, rclo, rsum, tex, vpex, xx

         INTEGER :: count1

         tx = ta
         enbal2 = 0.

          DO count1 = 0, 3
           DO
             hc = 2.67 + 6.5 * 0.1 ** 0.67
             hc = hc * ( p / po ) ** 0.55

      !--    Radiation

             aeff  = adu * feff
             rbare = aeff * ( 1. - facl ) * emsk * sigm *                            &
                 ( ( tx + 273.2 ) ** 4. - ( tsk + 273.2 ) ** 4. )
             rclo  = feff * acl * emcl * sigm *                                      &
                 ( ( tx + 273.2 ) ** 4. - ( tcl + 273.2 ) ** 4. )
             rsum  = rbare + rclo

      !--    Covection

             cbare = hc * ( tx - tsk ) * adu * ( 1. - facl )
             cclo  = hc * ( tx - tcl ) * acl
             csum  = cbare + cclo

      !--    Diffusion

             ed = evap / ( rdsk + rdcl ) * adu * ( 1. - wetsk ) * ( 12. - vpts )

      !--    Respiration

             tex  = 0.47 * tx + 21.
             eres = cair * ( tx - tex ) * rtv
             vpex = 6.11 * 10. ** ( 7.45 * tex / ( 235. + tex ) )
             erel = 0.623 * evap / p * ( 12. - vpex ) * rtv
             ere  = eres + erel

      !--    Energy ballance

             enbal = h + ed + ere + esw + csum + rsum

      !--    Iteration concerning ta

             IF ( count1 .EQ. 0 )  xx = 1.
             IF ( count1 .EQ. 1 )  xx = 0.1
             IF ( count1 .EQ. 2 )  xx = 0.01
             IF ( count1 .EQ. 3 )  xx = 0.001
             IF ( enbal .GT. 0. )  tx = tx - xx
             IF ( enbal .LT. 0. )  tx = tx + xx
             IF ( ( enbal .LE. 0. ) .AND. ( enbal2 .GT. 0. ) ) EXIT
             IF ( ( enbal .GE. 0. ) .AND. ( enbal2 .LT. 0. ) ) EXIT

             enbal2 = enbal
           END DO
         END DO
         RETURN
       END SUBROUTINE pet_iteration

  DOUBLE precision function UTCI_approx(Ta,Tmrt,va,rh)
  !~ **********************************************
  !~ DOUBLE PRECISION Function value is the UTCI in degree Celsius
  !~ computed by a 6th order approximating polynomial from the 4 Input paramters
  !~
  !~ Input parameters (all of type DOUBLE PRECISION)
  !~ - Ta   : air temperature, degree Celsius
  !~ - Tmrt : mean radiant temperature, degree Celsius
  !~ - va   : wind speed 10 m above ground level in m/s
  !~ - rh   : relative humidity in %

  implicit none
      !~ type of input of the argument list
      DOUBLE PRECISION Ta,va,Tmrt,Pa,D_Tmrt,rh;
      D_TMRT=Tmrt-Ta
      !~ calculate vapour pressure from relative humidity
      PA=(rh * es(Ta)) / 1000.0 !~ use vapour pressure in kPa
      !~ calculate 6th order polynomial as approximation
    UTCI_approx=Ta+&

  ( 6.07562052D-01 )   + &
  ( -2.27712343D-02 ) * Ta + &
  ( 8.06470249D-04 ) * Ta*Ta + &
  ( -1.54271372D-04 ) * Ta*Ta*Ta + &
  ( -3.24651735D-06 ) * Ta*Ta*Ta*Ta + &
  ( 7.32602852D-08 ) * Ta*Ta*Ta*Ta*Ta + &
  ( 1.35959073D-09 ) * Ta*Ta*Ta*Ta*Ta*Ta + &
  ( -2.25836520D+00 ) * va + &
  ( 8.80326035D-02 ) * Ta*va + &
  ( 2.16844454D-03 ) * Ta*Ta*va + &
  ( -1.53347087D-05 ) * Ta*Ta*Ta*va + &
  ( -5.72983704D-07 ) * Ta*Ta*Ta*Ta*va + &
  ( -2.55090145D-09 ) * Ta*Ta*Ta*Ta*Ta*va + &
  ( -7.51269505D-01 ) * va*va + &
  ( -4.08350271D-03 ) * Ta*va*va + &
  ( -5.21670675D-05 ) * Ta*Ta*va*va + &
  ( 1.94544667D-06 ) * Ta*Ta*Ta*va*va + &
  ( 1.14099531D-08 ) * Ta*Ta*Ta*Ta*va*va + &
  ( 1.58137256D-01 ) * va*va*va + &
  ( -6.57263143D-05 ) * Ta*va*va*va + &
  ( 2.22697524D-07 ) * Ta*Ta*va*va*va + &
  ( -4.16117031D-08 ) * Ta*Ta*Ta*va*va*va + &
  ( -1.27762753D-02 ) * va*va*va*va + &
  ( 9.66891875D-06 ) * Ta*va*va*va*va + &
  ( 2.52785852D-09 ) * Ta*Ta*va*va*va*va + &
  ( 4.56306672D-04 ) * va*va*va*va*va + &
  ( -1.74202546D-07 ) * Ta*va*va*va*va*va + &
  ( -5.91491269D-06 ) * va*va*va*va*va*va + &
  ( 3.98374029D-01 ) * D_Tmrt + &
  ( 1.83945314D-04 ) * Ta*D_Tmrt + &
  ( -1.73754510D-04 ) * Ta*Ta*D_Tmrt + &
  ( -7.60781159D-07 ) * Ta*Ta*Ta*D_Tmrt + &
  ( 3.77830287D-08 ) * Ta*Ta*Ta*Ta*D_Tmrt + &
  ( 5.43079673D-10 ) * Ta*Ta*Ta*Ta*Ta*D_Tmrt + &
  ( -2.00518269D-02 ) * va*D_Tmrt + &
  ( 8.92859837D-04 ) * Ta*va*D_Tmrt + &
  ( 3.45433048D-06 ) * Ta*Ta*va*D_Tmrt + &
  ( -3.77925774D-07 ) * Ta*Ta*Ta*va*D_Tmrt + &
  ( -1.69699377D-09 ) * Ta*Ta*Ta*Ta*va*D_Tmrt + &
  ( 1.69992415D-04 ) * va*va*D_Tmrt + &
  ( -4.99204314D-05 ) * Ta*va*va*D_Tmrt + &
  ( 2.47417178D-07 ) * Ta*Ta*va*va*D_Tmrt + &
  ( 1.07596466D-08 ) * Ta*Ta*Ta*va*va*D_Tmrt + &
  ( 8.49242932D-05 ) * va*va*va*D_Tmrt + &
  ( 1.35191328D-06 ) * Ta*va*va*va*D_Tmrt + &
  ( -6.21531254D-09 ) * Ta*Ta*va*va*va*D_Tmrt + &
  ( -4.99410301D-06 ) * va*va*va*va*D_Tmrt + &
  ( -1.89489258D-08 ) * Ta*va*va*va*va*D_Tmrt + &
  ( 8.15300114D-08 ) * va*va*va*va*va*D_Tmrt + &
  ( 7.55043090D-04 ) * D_Tmrt*D_Tmrt + &
  ( -5.65095215D-05 ) * Ta*D_Tmrt*D_Tmrt + &
  ( -4.52166564D-07 ) * Ta*Ta*D_Tmrt*D_Tmrt + &
  ( 2.46688878D-08 ) * Ta*Ta*Ta*D_Tmrt*D_Tmrt + &
  ( 2.42674348D-10 ) * Ta*Ta*Ta*Ta*D_Tmrt*D_Tmrt + &
  ( 1.54547250D-04 ) * va*D_Tmrt*D_Tmrt + &
  ( 5.24110970D-06 ) * Ta*va*D_Tmrt*D_Tmrt + &
  ( -8.75874982D-08 ) * Ta*Ta*va*D_Tmrt*D_Tmrt + &
  ( -1.50743064D-09 ) * Ta*Ta*Ta*va*D_Tmrt*D_Tmrt + &
  ( -1.56236307D-05 ) * va*va*D_Tmrt*D_Tmrt + &
  ( -1.33895614D-07 ) * Ta*va*va*D_Tmrt*D_Tmrt + &
  ( 2.49709824D-09 ) * Ta*Ta*va*va*D_Tmrt*D_Tmrt + &
  ( 6.51711721D-07 ) * va*va*va*D_Tmrt*D_Tmrt + &
  ( 1.94960053D-09 ) * Ta*va*va*va*D_Tmrt*D_Tmrt + &
  ( -1.00361113D-08 ) * va*va*va*va*D_Tmrt*D_Tmrt + &
  ( -1.21206673D-05 ) * D_Tmrt*D_Tmrt*D_Tmrt + &
  ( -2.18203660D-07 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 7.51269482D-09 ) * Ta*Ta*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 9.79063848D-11 ) * Ta*Ta*Ta*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 1.25006734D-06 ) * va*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( -1.81584736D-09 ) * Ta*va*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( -3.52197671D-10 ) * Ta*Ta*va*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( -3.36514630D-08 ) * va*va*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 1.35908359D-10 ) * Ta*va*va*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 4.17032620D-10 ) * va*va*va*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( -1.30369025D-09 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 4.13908461D-10 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 9.22652254D-12 ) * Ta*Ta*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( -5.08220384D-09 ) * va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( -2.24730961D-11 ) * Ta*va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 1.17139133D-10 ) * va*va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 6.62154879D-10 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 4.03863260D-13 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 1.95087203D-12 ) * va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( -4.73602469D-12 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 5.12733497D+00 ) * Pa + &
  ( -3.12788561D-01 ) * Ta*Pa + &
  ( -1.96701861D-02 ) * Ta*Ta*Pa + &
  ( 9.99690870D-04 ) * Ta*Ta*Ta*Pa + &
  ( 9.51738512D-06 ) * Ta*Ta*Ta*Ta*Pa + &
  ( -4.66426341D-07 ) * Ta*Ta*Ta*Ta*Ta*Pa + &
  ( 5.48050612D-01 ) * va*Pa + &
  ( -3.30552823D-03 ) * Ta*va*Pa + &
  ( -1.64119440D-03 ) * Ta*Ta*va*Pa + &
  ( -5.16670694D-06 ) * Ta*Ta*Ta*va*Pa + &
  ( 9.52692432D-07 ) * Ta*Ta*Ta*Ta*va*Pa + &
  ( -4.29223622D-02 ) * va*va*Pa + &
  ( 5.00845667D-03 ) * Ta*va*va*Pa + &
  ( 1.00601257D-06 ) * Ta*Ta*va*va*Pa + &
  ( -1.81748644D-06 ) * Ta*Ta*Ta*va*va*Pa + &
  ( -1.25813502D-03 ) * va*va*va*Pa + &
  ( -1.79330391D-04 ) * Ta*va*va*va*Pa + &
  ( 2.34994441D-06 ) * Ta*Ta*va*va*va*Pa + &
  ( 1.29735808D-04 ) * va*va*va*va*Pa + &
  ( 1.29064870D-06 ) * Ta*va*va*va*va*Pa + &
  ( -2.28558686D-06 ) * va*va*va*va*va*Pa + &
  ( -3.69476348D-02 ) * D_Tmrt*Pa + &
  ( 1.62325322D-03 ) * Ta*D_Tmrt*Pa + &
  ( -3.14279680D-05 ) * Ta*Ta*D_Tmrt*Pa + &
  ( 2.59835559D-06 ) * Ta*Ta*Ta*D_Tmrt*Pa + &
  ( -4.77136523D-08 ) * Ta*Ta*Ta*Ta*D_Tmrt*Pa + &
  ( 8.64203390D-03 ) * va*D_Tmrt*Pa + &
  ( -6.87405181D-04 ) * Ta*va*D_Tmrt*Pa + &
  ( -9.13863872D-06 ) * Ta*Ta*va*D_Tmrt*Pa + &
  ( 5.15916806D-07 ) * Ta*Ta*Ta*va*D_Tmrt*Pa + &
  ( -3.59217476D-05 ) * va*va*D_Tmrt*Pa + &
  ( 3.28696511D-05 ) * Ta*va*va*D_Tmrt*Pa + &
  ( -7.10542454D-07 ) * Ta*Ta*va*va*D_Tmrt*Pa + &
  ( -1.24382300D-05 ) * va*va*va*D_Tmrt*Pa + &
  ( -7.38584400D-09 ) * Ta*va*va*va*D_Tmrt*Pa + &
  ( 2.20609296D-07 ) * va*va*va*va*D_Tmrt*Pa + &
  ( -7.32469180D-04 ) * D_Tmrt*D_Tmrt*Pa + &
  ( -1.87381964D-05 ) * Ta*D_Tmrt*D_Tmrt*Pa + &
  ( 4.80925239D-06 ) * Ta*Ta*D_Tmrt*D_Tmrt*Pa + &
  ( -8.75492040D-08 ) * Ta*Ta*Ta*D_Tmrt*D_Tmrt*Pa + &
  ( 2.77862930D-05 ) * va*D_Tmrt*D_Tmrt*Pa + &
  ( -5.06004592D-06 ) * Ta*va*D_Tmrt*D_Tmrt*Pa + &
  ( 1.14325367D-07 ) * Ta*Ta*va*D_Tmrt*D_Tmrt*Pa + &
  ( 2.53016723D-06 ) * va*va*D_Tmrt*D_Tmrt*Pa + &
  ( -1.72857035D-08 ) * Ta*va*va*D_Tmrt*D_Tmrt*Pa + &
  ( -3.95079398D-08 ) * va*va*va*D_Tmrt*D_Tmrt*Pa + &
  ( -3.59413173D-07 ) * D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
  ( 7.04388046D-07 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
  ( -1.89309167D-08 ) * Ta*Ta*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
  ( -4.79768731D-07 ) * va*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
  ( 7.96079978D-09 ) * Ta*va*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
  ( 1.62897058D-09 ) * va*va*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
  ( 3.94367674D-08 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
  ( -1.18566247D-09 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
  ( 3.34678041D-10 ) * va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
  ( -1.15606447D-10 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
  ( -2.80626406D+00 ) * Pa*Pa + &
  ( 5.48712484D-01 ) * Ta*Pa*Pa + &
  ( -3.99428410D-03 ) * Ta*Ta*Pa*Pa + &
  ( -9.54009191D-04 ) * Ta*Ta*Ta*Pa*Pa + &
  ( 1.93090978D-05 ) * Ta*Ta*Ta*Ta*Pa*Pa + &
  ( -3.08806365D-01 ) * va*Pa*Pa + &
  ( 1.16952364D-02 ) * Ta*va*Pa*Pa + &
  ( 4.95271903D-04 ) * Ta*Ta*va*Pa*Pa + &
  ( -1.90710882D-05 ) * Ta*Ta*Ta*va*Pa*Pa + &
  ( 2.10787756D-03 ) * va*va*Pa*Pa + &
  ( -6.98445738D-04 ) * Ta*va*va*Pa*Pa + &
  ( 2.30109073D-05 ) * Ta*Ta*va*va*Pa*Pa + &
  ( 4.17856590D-04 ) * va*va*va*Pa*Pa + &
  ( -1.27043871D-05 ) * Ta*va*va*va*Pa*Pa + &
  ( -3.04620472D-06 ) * va*va*va*va*Pa*Pa + &
  ( 5.14507424D-02 ) * D_Tmrt*Pa*Pa + &
  ( -4.32510997D-03 ) * Ta*D_Tmrt*Pa*Pa + &
  ( 8.99281156D-05 ) * Ta*Ta*D_Tmrt*Pa*Pa + &
  ( -7.14663943D-07 ) * Ta*Ta*Ta*D_Tmrt*Pa*Pa + &
  ( -2.66016305D-04 ) * va*D_Tmrt*Pa*Pa + &
  ( 2.63789586D-04 ) * Ta*va*D_Tmrt*Pa*Pa + &
  ( -7.01199003D-06 ) * Ta*Ta*va*D_Tmrt*Pa*Pa + &
  ( -1.06823306D-04 ) * va*va*D_Tmrt*Pa*Pa + &
  ( 3.61341136D-06 ) * Ta*va*va*D_Tmrt*Pa*Pa + &
  ( 2.29748967D-07 ) * va*va*va*D_Tmrt*Pa*Pa + &
  ( 3.04788893D-04 ) * D_Tmrt*D_Tmrt*Pa*Pa + &
  ( -6.42070836D-05 ) * Ta*D_Tmrt*D_Tmrt*Pa*Pa + &
  ( 1.16257971D-06 ) * Ta*Ta*D_Tmrt*D_Tmrt*Pa*Pa + &
  ( 7.68023384D-06 ) * va*D_Tmrt*D_Tmrt*Pa*Pa + &
  ( -5.47446896D-07 ) * Ta*va*D_Tmrt*D_Tmrt*Pa*Pa + &
  ( -3.59937910D-08 ) * va*va*D_Tmrt*D_Tmrt*Pa*Pa + &
  ( -4.36497725D-06 ) * D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa + &
  ( 1.68737969D-07 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa + &
  ( 2.67489271D-08 ) * va*D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa + &
  ( 3.23926897D-09 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa + &
  ( -3.53874123D-02 ) * Pa*Pa*Pa + &
  ( -2.21201190D-01 ) * Ta*Pa*Pa*Pa + &
  ( 1.55126038D-02 ) * Ta*Ta*Pa*Pa*Pa + &
  ( -2.63917279D-04 ) * Ta*Ta*Ta*Pa*Pa*Pa + &
  ( 4.53433455D-02 ) * va*Pa*Pa*Pa + &
  ( -4.32943862D-03 ) * Ta*va*Pa*Pa*Pa + &
  ( 1.45389826D-04 ) * Ta*Ta*va*Pa*Pa*Pa + &
  ( 2.17508610D-04 ) * va*va*Pa*Pa*Pa + &
  ( -6.66724702D-05 ) * Ta*va*va*Pa*Pa*Pa + &
  ( 3.33217140D-05 ) * va*va*va*Pa*Pa*Pa + &
  ( -2.26921615D-03 ) * D_Tmrt*Pa*Pa*Pa + &
  ( 3.80261982D-04 ) * Ta*D_Tmrt*Pa*Pa*Pa + &
  ( -5.45314314D-09 ) * Ta*Ta*D_Tmrt*Pa*Pa*Pa + &
  ( -7.96355448D-04 ) * va*D_Tmrt*Pa*Pa*Pa + &
  ( 2.53458034D-05 ) * Ta*va*D_Tmrt*Pa*Pa*Pa + &
  ( -6.31223658D-06 ) * va*va*D_Tmrt*Pa*Pa*Pa + &
  ( 3.02122035D-04 ) * D_Tmrt*D_Tmrt*Pa*Pa*Pa + &
  ( -4.77403547D-06 ) * Ta*D_Tmrt*D_Tmrt*Pa*Pa*Pa + &
  ( 1.73825715D-06 ) * va*D_Tmrt*D_Tmrt*Pa*Pa*Pa + &
  ( -4.09087898D-07 ) * D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa*Pa + &
  ( 6.14155345D-01 ) * Pa*Pa*Pa*Pa + &
  ( -6.16755931D-02 ) * Ta*Pa*Pa*Pa*Pa + &
  ( 1.33374846D-03 ) * Ta*Ta*Pa*Pa*Pa*Pa + &
  ( 3.55375387D-03 ) * va*Pa*Pa*Pa*Pa + &
  ( -5.13027851D-04 ) * Ta*va*Pa*Pa*Pa*Pa + &
  ( 1.02449757D-04 ) * va*va*Pa*Pa*Pa*Pa + &
  ( -1.48526421D-03 ) * D_Tmrt*Pa*Pa*Pa*Pa + &
  ( -4.11469183D-05 ) * Ta*D_Tmrt*Pa*Pa*Pa*Pa + &
  ( -6.80434415D-06 ) * va*D_Tmrt*Pa*Pa*Pa*Pa + &
  ( -9.77675906D-06 ) * D_Tmrt*D_Tmrt*Pa*Pa*Pa*Pa + &
  ( 8.82773108D-02 ) * Pa*Pa*Pa*Pa*Pa + &
  ( -3.01859306D-03 ) * Ta*Pa*Pa*Pa*Pa*Pa + &
  ( 1.04452989D-03 ) * va*Pa*Pa*Pa*Pa*Pa + &
  ( 2.47090539D-04 ) * D_Tmrt*Pa*Pa*Pa*Pa*Pa + &
  ( 1.48348065D-03 ) * Pa*Pa*Pa*Pa*Pa*Pa
    return
    END

!~ **********************************************
  DOUBLE precision function es(ta)
!~ **********************************************
!~ calculates saturation vapour pressure over water in hPa for input air temperature (ta) in celsius according to:
!~ Hardy, R.; ITS-90 Formulations for Vapor Pressure, Frostpoint Temperature, Dewpoint Temperature and Enhancement Factors in the Range -100 to 100 �C;
!~ Proceedings of Third International Symposium on Humidity and Moisture; edited by National Physical Laboratory (NPL), London, 1998, pp. 214-221
!~ http://www.thunderscientific.com/tech_info/reflibrary/its90formulas.pdf (retrieved 2008-10-01)

    implicit none
!
    DOUBLE precision ta, tk
    INTEGER I
    REAL :: g(0:7)=(/&
      -2.8365744E3,&
      -6.028076559E3,&
      1.954263612E1,&
      -2.737830188E-2,&
      1.6261698E-5,&
      7.0229056E-10,&
      -1.8680009E-13,&
      2.7150305 /)
!
    tk=ta+273.15 ! air temp in K
    es=g(7)*log(tk)
    do i=0,6
      es=es+g(i)*tk**(i-2)
    end do
    es=exp(es)*0.01 ! *0.01: convert Pa to hPa
!
    return
  END

  function es_vectorized(array) result(new_array)
      ! Subroutine to add 10 to every element in the array
      implicit none
      real(kind=8), intent(in) :: array(:)  ! Assumes a 1D array of any size
      real(kind=8) :: new_array(size(array))
      integer :: i

      do i = 1, size(array)
          new_array(i) = 6.112 * exp((17.67 * array(i)) / (array(i) + 243.5))
      end do
  end function es_vectorized

  function UTCI_approx_vectorized(Ta,Tmrt,va,rh) result(result_array)
  !~ **********************************************
  !~ DOUBLE PRECISION Function value is the UTCI in degree Celsius
  !~ computed by a 6th order approximating polynomial from the 4 Input paramters
  !~
  !~ Input parameters (all of type DOUBLE PRECISION)
  !~ - Ta   : air temperature, degree Celsius
  !~ - Tmrt : mean radiant temperature, degree Celsius
  !~ - va   : wind speed 10 m above ground level in m/s
  !~ - rh   : relative humidity in %

  implicit none
      real(kind=8), intent(in) :: Ta(:), va(:), Tmrt(:), rh(:)
      real(kind=8) :: D_TMRT(size(Ta)), PA(size(Ta))
      real(kind=8) :: result_array(size(Ta))

      !~ type of input of the argument list
      ! DOUBLE PRECISION Ta,va,Tmrt,ehPa,Pa,D_Tmrt,rh;
      D_TMRT=Tmrt-Ta
  !     !~ calculate vapour pressure from relative humidity
      PA=(rh * es_vectorized(Ta)) / 1000.0  !~ use vapour pressure in kPa
  !     !~ calculate 6th order polynomial as approximation
    result_array=Ta+&
  ( 6.07562052D-01 )   + &
  ( -2.27712343D-02 ) * Ta + &
  ( 8.06470249D-04 ) * Ta*Ta + &
  ( -1.54271372D-04 ) * Ta*Ta*Ta + &
  ( -3.24651735D-06 ) * Ta*Ta*Ta*Ta + &
  ( 7.32602852D-08 ) * Ta*Ta*Ta*Ta*Ta + &
  ( 1.35959073D-09 ) * Ta*Ta*Ta*Ta*Ta*Ta + &
  ( -2.25836520D+00 ) * va + &
  ( 8.80326035D-02 ) * Ta*va + &
  ( 2.16844454D-03 ) * Ta*Ta*va + &
  ( -1.53347087D-05 ) * Ta*Ta*Ta*va + &
  ( -5.72983704D-07 ) * Ta*Ta*Ta*Ta*va + &
  ( -2.55090145D-09 ) * Ta*Ta*Ta*Ta*Ta*va + &
  ( -7.51269505D-01 ) * va*va + &
  ( -4.08350271D-03 ) * Ta*va*va + &
  ( -5.21670675D-05 ) * Ta*Ta*va*va + &
  ( 1.94544667D-06 ) * Ta*Ta*Ta*va*va + &
  ( 1.14099531D-08 ) * Ta*Ta*Ta*Ta*va*va + &
  ( 1.58137256D-01 ) * va*va*va + &
  ( -6.57263143D-05 ) * Ta*va*va*va + &
  ( 2.22697524D-07 ) * Ta*Ta*va*va*va + &
  ( -4.16117031D-08 ) * Ta*Ta*Ta*va*va*va + &
  ( -1.27762753D-02 ) * va*va*va*va + &
  ( 9.66891875D-06 ) * Ta*va*va*va*va + &
  ( 2.52785852D-09 ) * Ta*Ta*va*va*va*va + &
  ( 4.56306672D-04 ) * va*va*va*va*va + &
  ( -1.74202546D-07 ) * Ta*va*va*va*va*va + &
  ( -5.91491269D-06 ) * va*va*va*va*va*va + &
  ( 3.98374029D-01 ) * D_Tmrt + &
  ( 1.83945314D-04 ) * Ta*D_Tmrt + &
  ( -1.73754510D-04 ) * Ta*Ta*D_Tmrt + &
  ( -7.60781159D-07 ) * Ta*Ta*Ta*D_Tmrt + &
  ( 3.77830287D-08 ) * Ta*Ta*Ta*Ta*D_Tmrt + &
  ( 5.43079673D-10 ) * Ta*Ta*Ta*Ta*Ta*D_Tmrt + &
  ( -2.00518269D-02 ) * va*D_Tmrt + &
  ( 8.92859837D-04 ) * Ta*va*D_Tmrt + &
  ( 3.45433048D-06 ) * Ta*Ta*va*D_Tmrt + &
  ( -3.77925774D-07 ) * Ta*Ta*Ta*va*D_Tmrt + &
  ( -1.69699377D-09 ) * Ta*Ta*Ta*Ta*va*D_Tmrt + &
  ( 1.69992415D-04 ) * va*va*D_Tmrt + &
  ( -4.99204314D-05 ) * Ta*va*va*D_Tmrt + &
  ( 2.47417178D-07 ) * Ta*Ta*va*va*D_Tmrt + &
  ( 1.07596466D-08 ) * Ta*Ta*Ta*va*va*D_Tmrt + &
  ( 8.49242932D-05 ) * va*va*va*D_Tmrt + &
  ( 1.35191328D-06 ) * Ta*va*va*va*D_Tmrt + &
  ( -6.21531254D-09 ) * Ta*Ta*va*va*va*D_Tmrt + &
  ( -4.99410301D-06 ) * va*va*va*va*D_Tmrt + &
  ( -1.89489258D-08 ) * Ta*va*va*va*va*D_Tmrt + &
  ( 8.15300114D-08 ) * va*va*va*va*va*D_Tmrt + &
  ( 7.55043090D-04 ) * D_Tmrt*D_Tmrt + &
  ( -5.65095215D-05 ) * Ta*D_Tmrt*D_Tmrt + &
  ( -4.52166564D-07 ) * Ta*Ta*D_Tmrt*D_Tmrt + &
  ( 2.46688878D-08 ) * Ta*Ta*Ta*D_Tmrt*D_Tmrt + &
  ( 2.42674348D-10 ) * Ta*Ta*Ta*Ta*D_Tmrt*D_Tmrt + &
  ( 1.54547250D-04 ) * va*D_Tmrt*D_Tmrt + &
  ( 5.24110970D-06 ) * Ta*va*D_Tmrt*D_Tmrt + &
  ( -8.75874982D-08 ) * Ta*Ta*va*D_Tmrt*D_Tmrt + &
  ( -1.50743064D-09 ) * Ta*Ta*Ta*va*D_Tmrt*D_Tmrt + &
  ( -1.56236307D-05 ) * va*va*D_Tmrt*D_Tmrt + &
  ( -1.33895614D-07 ) * Ta*va*va*D_Tmrt*D_Tmrt + &
  ( 2.49709824D-09 ) * Ta*Ta*va*va*D_Tmrt*D_Tmrt + &
  ( 6.51711721D-07 ) * va*va*va*D_Tmrt*D_Tmrt + &
  ( 1.94960053D-09 ) * Ta*va*va*va*D_Tmrt*D_Tmrt + &
  ( -1.00361113D-08 ) * va*va*va*va*D_Tmrt*D_Tmrt + &
  ( -1.21206673D-05 ) * D_Tmrt*D_Tmrt*D_Tmrt + &
  ( -2.18203660D-07 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 7.51269482D-09 ) * Ta*Ta*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 9.79063848D-11 ) * Ta*Ta*Ta*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 1.25006734D-06 ) * va*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( -1.81584736D-09 ) * Ta*va*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( -3.52197671D-10 ) * Ta*Ta*va*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( -3.36514630D-08 ) * va*va*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 1.35908359D-10 ) * Ta*va*va*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 4.17032620D-10 ) * va*va*va*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( -1.30369025D-09 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 4.13908461D-10 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 9.22652254D-12 ) * Ta*Ta*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( -5.08220384D-09 ) * va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( -2.24730961D-11 ) * Ta*va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 1.17139133D-10 ) * va*va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 6.62154879D-10 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 4.03863260D-13 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 1.95087203D-12 ) * va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( -4.73602469D-12 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
  ( 5.12733497D+00 ) * Pa + &
  ( -3.12788561D-01 ) * Ta*Pa + &
  ( -1.96701861D-02 ) * Ta*Ta*Pa + &
  ( 9.99690870D-04 ) * Ta*Ta*Ta*Pa + &
  ( 9.51738512D-06 ) * Ta*Ta*Ta*Ta*Pa + &
  ( -4.66426341D-07 ) * Ta*Ta*Ta*Ta*Ta*Pa + &
  ( 5.48050612D-01 ) * va*Pa + &
  ( -3.30552823D-03 ) * Ta*va*Pa + &
  ( -1.64119440D-03 ) * Ta*Ta*va*Pa + &
  ( -5.16670694D-06 ) * Ta*Ta*Ta*va*Pa + &
  ( 9.52692432D-07 ) * Ta*Ta*Ta*Ta*va*Pa + &
  ( -4.29223622D-02 ) * va*va*Pa + &
  ( 5.00845667D-03 ) * Ta*va*va*Pa + &
  ( 1.00601257D-06 ) * Ta*Ta*va*va*Pa + &
  ( -1.81748644D-06 ) * Ta*Ta*Ta*va*va*Pa + &
  ( -1.25813502D-03 ) * va*va*va*Pa + &
  ( -1.79330391D-04 ) * Ta*va*va*va*Pa + &
  ( 2.34994441D-06 ) * Ta*Ta*va*va*va*Pa + &
  ( 1.29735808D-04 ) * va*va*va*va*Pa + &
  ( 1.29064870D-06 ) * Ta*va*va*va*va*Pa + &
  ( -2.28558686D-06 ) * va*va*va*va*va*Pa + &
  ( -3.69476348D-02 ) * D_Tmrt*Pa + &
  ( 1.62325322D-03 ) * Ta*D_Tmrt*Pa + &
  ( -3.14279680D-05 ) * Ta*Ta*D_Tmrt*Pa + &
  ( 2.59835559D-06 ) * Ta*Ta*Ta*D_Tmrt*Pa + &
  ( -4.77136523D-08 ) * Ta*Ta*Ta*Ta*D_Tmrt*Pa + &
  ( 8.64203390D-03 ) * va*D_Tmrt*Pa + &
  ( -6.87405181D-04 ) * Ta*va*D_Tmrt*Pa + &
  ( -9.13863872D-06 ) * Ta*Ta*va*D_Tmrt*Pa + &
  ( 5.15916806D-07 ) * Ta*Ta*Ta*va*D_Tmrt*Pa + &
  ( -3.59217476D-05 ) * va*va*D_Tmrt*Pa + &
  ( 3.28696511D-05 ) * Ta*va*va*D_Tmrt*Pa + &
  ( -7.10542454D-07 ) * Ta*Ta*va*va*D_Tmrt*Pa + &
  ( -1.24382300D-05 ) * va*va*va*D_Tmrt*Pa + &
  ( -7.38584400D-09 ) * Ta*va*va*va*D_Tmrt*Pa + &
  ( 2.20609296D-07 ) * va*va*va*va*D_Tmrt*Pa + &
  ( -7.32469180D-04 ) * D_Tmrt*D_Tmrt*Pa + &
  ( -1.87381964D-05 ) * Ta*D_Tmrt*D_Tmrt*Pa + &
  ( 4.80925239D-06 ) * Ta*Ta*D_Tmrt*D_Tmrt*Pa + &
  ( -8.75492040D-08 ) * Ta*Ta*Ta*D_Tmrt*D_Tmrt*Pa + &
  ( 2.77862930D-05 ) * va*D_Tmrt*D_Tmrt*Pa + &
  ( -5.06004592D-06 ) * Ta*va*D_Tmrt*D_Tmrt*Pa + &
  ( 1.14325367D-07 ) * Ta*Ta*va*D_Tmrt*D_Tmrt*Pa + &
  ( 2.53016723D-06 ) * va*va*D_Tmrt*D_Tmrt*Pa + &
  ( -1.72857035D-08 ) * Ta*va*va*D_Tmrt*D_Tmrt*Pa + &
  ( -3.95079398D-08 ) * va*va*va*D_Tmrt*D_Tmrt*Pa + &
  ( -3.59413173D-07 ) * D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
  ( 7.04388046D-07 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
  ( -1.89309167D-08 ) * Ta*Ta*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
  ( -4.79768731D-07 ) * va*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
  ( 7.96079978D-09 ) * Ta*va*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
  ( 1.62897058D-09 ) * va*va*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
  ( 3.94367674D-08 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
  ( -1.18566247D-09 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
  ( 3.34678041D-10 ) * va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
  ( -1.15606447D-10 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
  ( -2.80626406D+00 ) * Pa*Pa + &
  ( 5.48712484D-01 ) * Ta*Pa*Pa + &
  ( -3.99428410D-03 ) * Ta*Ta*Pa*Pa + &
  ( -9.54009191D-04 ) * Ta*Ta*Ta*Pa*Pa + &
  ( 1.93090978D-05 ) * Ta*Ta*Ta*Ta*Pa*Pa + &
  ( -3.08806365D-01 ) * va*Pa*Pa + &
  ( 1.16952364D-02 ) * Ta*va*Pa*Pa + &
  ( 4.95271903D-04 ) * Ta*Ta*va*Pa*Pa + &
  ( -1.90710882D-05 ) * Ta*Ta*Ta*va*Pa*Pa + &
  ( 2.10787756D-03 ) * va*va*Pa*Pa + &
  ( -6.98445738D-04 ) * Ta*va*va*Pa*Pa + &
  ( 2.30109073D-05 ) * Ta*Ta*va*va*Pa*Pa + &
  ( 4.17856590D-04 ) * va*va*va*Pa*Pa + &
  ( -1.27043871D-05 ) * Ta*va*va*va*Pa*Pa + &
  ( -3.04620472D-06 ) * va*va*va*va*Pa*Pa + &
  ( 5.14507424D-02 ) * D_Tmrt*Pa*Pa + &
  ( -4.32510997D-03 ) * Ta*D_Tmrt*Pa*Pa + &
  ( 8.99281156D-05 ) * Ta*Ta*D_Tmrt*Pa*Pa + &
  ( -7.14663943D-07 ) * Ta*Ta*Ta*D_Tmrt*Pa*Pa + &
  ( -2.66016305D-04 ) * va*D_Tmrt*Pa*Pa + &
  ( 2.63789586D-04 ) * Ta*va*D_Tmrt*Pa*Pa + &
  ( -7.01199003D-06 ) * Ta*Ta*va*D_Tmrt*Pa*Pa + &
  ( -1.06823306D-04 ) * va*va*D_Tmrt*Pa*Pa + &
  ( 3.61341136D-06 ) * Ta*va*va*D_Tmrt*Pa*Pa + &
  ( 2.29748967D-07 ) * va*va*va*D_Tmrt*Pa*Pa + &
  ( 3.04788893D-04 ) * D_Tmrt*D_Tmrt*Pa*Pa + &
  ( -6.42070836D-05 ) * Ta*D_Tmrt*D_Tmrt*Pa*Pa + &
  ( 1.16257971D-06 ) * Ta*Ta*D_Tmrt*D_Tmrt*Pa*Pa + &
  ( 7.68023384D-06 ) * va*D_Tmrt*D_Tmrt*Pa*Pa + &
  ( -5.47446896D-07 ) * Ta*va*D_Tmrt*D_Tmrt*Pa*Pa + &
  ( -3.59937910D-08 ) * va*va*D_Tmrt*D_Tmrt*Pa*Pa + &
  ( -4.36497725D-06 ) * D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa + &
  ( 1.68737969D-07 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa + &
  ( 2.67489271D-08 ) * va*D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa + &
  ( 3.23926897D-09 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa + &
  ( -3.53874123D-02 ) * Pa*Pa*Pa + &
  ( -2.21201190D-01 ) * Ta*Pa*Pa*Pa + &
  ( 1.55126038D-02 ) * Ta*Ta*Pa*Pa*Pa + &
  ( -2.63917279D-04 ) * Ta*Ta*Ta*Pa*Pa*Pa + &
  ( 4.53433455D-02 ) * va*Pa*Pa*Pa + &
  ( -4.32943862D-03 ) * Ta*va*Pa*Pa*Pa + &
  ( 1.45389826D-04 ) * Ta*Ta*va*Pa*Pa*Pa + &
  ( 2.17508610D-04 ) * va*va*Pa*Pa*Pa + &
  ( -6.66724702D-05 ) * Ta*va*va*Pa*Pa*Pa + &
  ( 3.33217140D-05 ) * va*va*va*Pa*Pa*Pa + &
  ( -2.26921615D-03 ) * D_Tmrt*Pa*Pa*Pa + &
  ( 3.80261982D-04 ) * Ta*D_Tmrt*Pa*Pa*Pa + &
  ( -5.45314314D-09 ) * Ta*Ta*D_Tmrt*Pa*Pa*Pa + &
  ( -7.96355448D-04 ) * va*D_Tmrt*Pa*Pa*Pa + &
  ( 2.53458034D-05 ) * Ta*va*D_Tmrt*Pa*Pa*Pa + &
  ( -6.31223658D-06 ) * va*va*D_Tmrt*Pa*Pa*Pa + &
  ( 3.02122035D-04 ) * D_Tmrt*D_Tmrt*Pa*Pa*Pa + &
  ( -4.77403547D-06 ) * Ta*D_Tmrt*D_Tmrt*Pa*Pa*Pa + &
  ( 1.73825715D-06 ) * va*D_Tmrt*D_Tmrt*Pa*Pa*Pa + &
  ( -4.09087898D-07 ) * D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa*Pa + &
  ( 6.14155345D-01 ) * Pa*Pa*Pa*Pa + &
  ( -6.16755931D-02 ) * Ta*Pa*Pa*Pa*Pa + &
  ( 1.33374846D-03 ) * Ta*Ta*Pa*Pa*Pa*Pa + &
  ( 3.55375387D-03 ) * va*Pa*Pa*Pa*Pa + &
  ( -5.13027851D-04 ) * Ta*va*Pa*Pa*Pa*Pa + &
  ( 1.02449757D-04 ) * va*va*Pa*Pa*Pa*Pa + &
  ( -1.48526421D-03 ) * D_Tmrt*Pa*Pa*Pa*Pa + &
  ( -4.11469183D-05 ) * Ta*D_Tmrt*Pa*Pa*Pa*Pa + &
  ( -6.80434415D-06 ) * va*D_Tmrt*Pa*Pa*Pa*Pa + &
  ( -9.77675906D-06 ) * D_Tmrt*D_Tmrt*Pa*Pa*Pa*Pa + &
  ( 8.82773108D-02 ) * Pa*Pa*Pa*Pa*Pa + &
  ( -3.01859306D-03 ) * Ta*Pa*Pa*Pa*Pa*Pa + &
  ( 1.04452989D-03 ) * va*Pa*Pa*Pa*Pa*Pa + &
  ( 2.47090539D-04 ) * D_Tmrt*Pa*Pa*Pa*Pa*Pa + &
  ( 1.48348065D-03 ) * Pa*Pa*Pa*Pa*Pa*Pa
    END



end module thermal_comfort_mod
