!===========================================================================
!===========================================================================
!This file is part of ModelLib.
!
!    ModelLib is a free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ModelLib is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ModelLib.  If not, see <http://www.gnu.org/licenses/>.
!
!    Copyright 2016  David LAUVERGNAT, Félix MOUHAT and Liang LIANG
!
!===========================================================================
!===========================================================================

!> @brief Module which makes the initialization, calculation of the PSB3 potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!> @author Emanuele Marsili 
!! @date 26/04/2019
!!
MODULE mod_PSB3Pot

  USE mod_NumParameters

  IMPLICIT NONE

!> @brief Derived type in which the PSB3 parameters are set-up.
!!
!! @param option                  integer: it enables to chose between the 2 models (default 1, not published model)

  TYPE Param_PSB3

     PRIVATE

      real (kind=Rkind) :: blamin   = 0.0912615_Rkind
      real (kind=Rkind) :: blaTSdir = 0.025079167_Rkind
      real (kind=Rkind) :: deepth   = 2000_Rkind

      real (kind=Rkind) :: kf1 = 3733.5_Rkind
      real (kind=Rkind) :: d2  = 54.633_Rkind
      real (kind=Rkind) :: d3  = 3.8008_Rkind
      real (kind=Rkind) :: kf4 = 1097.19_Rkind

      real (kind=Rkind) :: c1  = 437.068_Rkind
      real (kind=Rkind) :: c2  = 16.7317_Rkind
      real (kind=Rkind) :: c3  = 7.35468_Rkind
      real (kind=Rkind) :: c4  = 88.517_Rkind
      real (kind=Rkind) :: c5  = 5.95221_Rkind

      real (kind=Rkind) :: h1  = 155.749_Rkind

      real (kind=Rkind) :: k1  = 24.0358_Rkind

!---Additional parameter for the Published model-------!
      real (kind=Rkind) :: d1  = 2571.3_Rkind
      real (kind=Rkind) :: d4  = 594.97_Rkind

      real (kind=Rkind) :: hd1  = 126.22_Rkind
      real (kind=Rkind) :: hd2  = 75.003_Rkind
      real (kind=Rkind) :: hc1  = 86.464_Rkind
      real (kind=Rkind) :: hc2  = 43.461_Rkind

      real (kind=Rkind) :: k3  = 1.1347_Rkind
!-----------------------------------------------------!

     ! Warning the parameters are given as in the publication.
     !   Therefore, the BLA(=Q(1)) is in Angstrom and the energy is in kcal.mol^-1.

      integer :: option    = 1
      logical :: PubliUnit = .FALSE.
 
  END TYPE Param_PSB3
 
  PRIVATE Read_PSB3Pot,Init0_PSB3Pot,eval_PSB3Pot1,eval_PSB3Pot2
 
  CONTAINS
!> @brief Subroutine which makes the initialization of the PSB3 parameters.
!!
!! @param Para_PSB3          TYPE(Param_PSB3):   derived type in which the parameters are set-up.
!! @param option             integer:            to be able to chose between the 3 models (default 1, Simple avoided crossing).
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.

  SUBROUTINE Init_PSB3Pot(Para_PSB3,option,nio,read_param,PubliUnit,blamin,blaTSdir,deepth,&
                          kf1,d1,d2,d3,d4,kf4,c1,c2,c3,c4,c5,h1,hd1,hd2,hc1,hc2,k1,k3)
    IMPLICIT NONE

    TYPE (Param_PSB3),           intent(inout) :: Para_PSB3
    integer,                     intent(in)    :: option
    integer,           optional, intent(in)    :: nio
    logical,           optional, intent(in)    :: read_param,PubliUnit              
    real (kind=Rkind), optional, intent(in)    :: blamin,blaTSdir,deepth,kf1,d1,d2,d3,d4,kf4, &
                                                  c1,c2,c3,c4,c5,h1,hd1,hd2,hc1,hc2,k1,k3
    logical :: read_param_loc

    IF (present(PubliUnit)) Para_PSB3%PubliUnit = PubliUnit

    read_param_loc = .FALSE.
    IF (present(read_param)) read_param_loc = read_param
    IF (read_param_loc .AND. .NOT. present(nio)) THEN
       write(out_unitp,*) ' ERROR in Init_PSB3Pot '
       write(out_unitp,*) ' read_param = t and The file unit (nio) is not present '
       write(out_unitp,*) ' => impossible to read the input file '
       STOP 'STOP in Init_PSB3Pot: impossible to read the input file. The file unit (nio) is not present'
    END IF

    Para_PSB3%option = option

    IF (Para_PSB3%option < 1 .OR. Para_PSB3%option > 2) Para_PSB3%option = 1

    IF (read_param_loc) THEN
      CALL Read_PSB3Pot(Para_PSB3,nio)
      CALL Init0_PSB3Pot(Para_PSB3)

    ELSE

      Para_PSB3%blamin   = Huge(ONE)
      Para_PSB3%blaTSdir = Huge(ONE)
      Para_PSB3%deepth   = Huge(ONE)
      Para_PSB3%kf1      = Huge(ONE)
      Para_PSB3%d2       = Huge(ONE)
      Para_PSB3%d3       = Huge(ONE)
      Para_PSB3%kf4      = Huge(ONE)
      Para_PSB3%d1       = Huge(ONE)
      Para_PSB3%d4       = Huge(ONE)
      Para_PSB3%c1       = Huge(ONE)
      Para_PSB3%c2       = Huge(ONE)
      Para_PSB3%c3       = Huge(ONE)
      Para_PSB3%c4       = Huge(ONE)
      Para_PSB3%c5       = Huge(ONE)
      Para_PSB3%h1       = Huge(ONE)
      Para_PSB3%hc1      = Huge(ONE)
      Para_PSB3%hc2      = Huge(ONE)
      Para_PSB3%hd1      = Huge(ONE)
      Para_PSB3%hd2      = Huge(ONE)
      Para_PSB3%k1       = Huge(ONE)
      Para_PSB3%k3       = Huge(ONE)

      CALL Init0_PSB3Pot(Para_PSB3)

      SELECT CASE (Para_PSB3%option)
      CASE (1) ! Not published potential 

        IF (present(blamin))   Para_PSB3%blamin   = blamin
        IF (present(blaTSdir)) Para_PSB3%blaTSdir = blaTSdir
        IF (present(deepth))   Para_PSB3%deepth   = deepth
        IF (present(kf1)) Para_PSB3%kf1    = kf1 
        IF (present(d2))  Para_PSB3%d2     = d2
        IF (present(d3))  Para_PSB3%d3     = d3
        IF (present(kf4)) Para_PSB3%kf4    = kf4
        IF (present(c1))  Para_PSB3%c1     = c1
        IF (present(c2))  Para_PSB3%c2     = c2
        IF (present(c3))  Para_PSB3%c3     = c3
        IF (present(c4))  Para_PSB3%c4     = c4
        IF (present(c5))  Para_PSB3%c5     = c5
        IF (present(h1))  Para_PSB3%h1     = h1
        IF (present(k1))  Para_PSB3%k1     = k1

      CASE (2) ! Published potential 

        IF (present(blamin))   Para_PSB3%blamin   = blamin
        IF (present(kf1)) Para_PSB3%d1     = d1  
        IF (present(d2))  Para_PSB3%d2     = d2
        IF (present(d3))  Para_PSB3%d3     = d3
        IF (present(kf4)) Para_PSB3%d4     = d4
        IF (present(c1))  Para_PSB3%c1     = c1
        IF (present(c2))  Para_PSB3%c2     = c2
        IF (present(c3))  Para_PSB3%c3     = c3
        IF (present(c4))  Para_PSB3%c4     = c4
        IF (present(c5))  Para_PSB3%c5     = c5
        IF (present(h1))  Para_PSB3%hd1    = hd1
        IF (present(h1))  Para_PSB3%hd2    = hd2
        IF (present(h1))  Para_PSB3%hc1    = hc1
        IF (present(h1))  Para_PSB3%hc2    = hc2
        IF (present(k1))  Para_PSB3%k1     = k1
        IF (present(k1))  Para_PSB3%k3     = k3

      CASE Default

          write(out_unitp,*) ' ERROR in Init_PSB3Pot '
          write(out_unitp,*) ' This option is not possible. option: ',Para_PSB3%option
          write(out_unitp,*) ' Its value MUST be 1 or 2'

          STOP
      END SELECT
    END IF

    IF (Para_PSB3%PubliUnit) THEN
      write(out_unitp,*) 'PubliUnit=.TRUE.,  Q:[Angs,Rad,Rad], Energy: [kcal.mol^-1]'
    ELSE
      write(out_unitp,*) 'PubliUnit=.FALSE., Q:[Bohr,Rad,Rad], Energy: [Hartree]'
    END IF

  END SUBROUTINE Init_PSB3Pot

!> @brief Subroutine which makes the initialization of the PSB3 parameters from his paper.
!!
!! @param Para_PSB3         TYPE(Param_PSB3):  derived type in which the parameters are set-up.

  SUBROUTINE Init0_PSB3Pot(Para_PSB3)
    
    IMPLICIT NONE

    TYPE (Param_PSB3), intent(inout) :: Para_PSB3

    character (len=*), parameter     :: name_sub='Init0_PSB3Pot'

    IF (Para_PSB3%option < 1 .OR. Para_PSB3%option > 2) Para_PSB3%option = 1

    SELECT CASE (Para_PSB3%option)
    CASE (1) ! unpublished model
     
      IF (Para_PSB3%blamin   == huge(ONE)) Para_PSB3%blamin   = 0.0912615_Rkind
      IF (Para_PSB3%blaTSdir == huge(ONE)) Para_PSB3%blaTSdir = 0.025079167_Rkind
      IF (Para_PSB3%deepth   == huge(ONE)) Para_PSB3%deepth   = 2000_Rkind

      IF (Para_PSB3%kf1 == huge(ONE)) Para_PSB3%kf1 = 3733.5_Rkind
      IF (Para_PSB3%d2  == huge(ONE)) Para_PSB3%d2  = 54.633_Rkind
      IF (Para_PSB3%d3  == huge(ONE)) Para_PSB3%d3  = 3.8008_Rkind
      IF (Para_PSB3%kf4 == huge(ONE)) Para_PSB3%kf4 = 1097.19_Rkind

      IF (Para_PSB3%c1  == huge(ONE)) Para_PSB3%c1 = 437.068_Rkind
      IF (Para_PSB3%c2  == huge(ONE)) Para_PSB3%c2 = 16.7317_Rkind
      IF (Para_PSB3%c3  == huge(ONE)) Para_PSB3%c3 = 7.35468_Rkind
      IF (Para_PSB3%c4  == huge(ONE)) Para_PSB3%c4 = 88.517_Rkind
      IF (Para_PSB3%c5  == huge(ONE)) Para_PSB3%c5 = 5.95221_Rkind

      IF (Para_PSB3%h1  == huge(ONE)) Para_PSB3%h1 = 155.749_Rkind
      IF (Para_PSB3%k1  == huge(ONE)) Para_PSB3%k1 = 24.0358_Rkind
    
    CASE (2) ! Published model 
   
      IF (Para_PSB3%blamin   == huge(ONE)) Para_PSB3%blamin   = 0.0912615_Rkind
      IF (Para_PSB3%d1  == huge(ONE)) Para_PSB3%d1  = 2571.3_Rkind
      IF (Para_PSB3%d2  == huge(ONE)) Para_PSB3%d2  = 54.633_Rkind
      IF (Para_PSB3%d3  == huge(ONE)) Para_PSB3%d3  = 3.8008_Rkind
      IF (Para_PSB3%d4  == huge(ONE)) Para_PSB3%d4  = 594.97_Rkind

      IF (Para_PSB3%c1  == huge(ONE)) Para_PSB3%c1 = 437.068_Rkind
      IF (Para_PSB3%c2  == huge(ONE)) Para_PSB3%c2 = 16.7317_Rkind
      IF (Para_PSB3%c3  == huge(ONE)) Para_PSB3%c3 = 7.35468_Rkind
      IF (Para_PSB3%c4  == huge(ONE)) Para_PSB3%c4 = 88.517_Rkind
      IF (Para_PSB3%c5  == huge(ONE)) Para_PSB3%c5 = 5.95221_Rkind
      
      IF (Para_PSB3%hd1 == huge(ONE)) Para_PSB3%hd1 = 126.22_Rkind
      IF (Para_PSB3%hd2 == huge(ONE)) Para_PSB3%hd2 = 75.003_Rkind
      IF (Para_PSB3%hc1 == huge(ONE)) Para_PSB3%hc1 = 86.464_Rkind
      IF (Para_PSB3%hc2 == huge(ONE)) Para_PSB3%hc2 = 43.461_Rkind

      IF (Para_PSB3%k1  == huge(ONE)) Para_PSB3%k1  = 13.699_Rkind
      IF (Para_PSB3%k3  == huge(ONE)) Para_PSB3%k3  = 1.1347_Rkind

    CASE Default
        write(out_unitp,*) ' ERROR in Init0_PSB3Pot'
        write(out_unitp,*) ' This option is not possible. option:',Para_PSB3%option
        write(out_unitp,*) ' Its value MUST be 1 or 2'
        STOP
    END SELECT

  END SUBROUTINE Init0_PSB3Pot

!> @brief Subroutine wich reads the PSB3 parameters with a namelist.
!!   This can be called only from the "Init_PSB3Pot" subroutine.
!!
!! @param Para_PSB3         TYPE(Param_PSB3):   derived type in which the parameters are set-up.
!! @param nio                integer:           file unit to read the parameters.
 
  SUBROUTINE Read_PSB3Pot(Para_PSB3,nio)
  
    TYPE (Param_PSB3),     intent(inout) :: Para_PSB3
    integer          ,     intent(in)    :: nio
    real (kind=Rkind)      :: blamin,blaTSdir,deepth,kf1,d1,d2,d3,d4,kf4,      &
                              c1,c2,c3,c4,c5,h1,hd1,hd2,hc1,hc2,k1,k3
    integer                :: err_read

    namelist /PSB3/ blamin,blaTSdir,deepth,kf1,d1,d2,d3,d4,kf4,c1,c2,c3,c4,c5, &
                    h1,hd1,hd2,hc1,hc2,k1,k3

     blamin   = Huge(ONE)
     blaTSdir = Huge(ONE)
     deepth   = Huge(ONE)
     kf1      = Huge(ONE)
     d2       = Huge(ONE)
     d3       = Huge(ONE)
     kf4      = Huge(ONE)
     c1       = Huge(ONE)
     c2       = Huge(ONE)
     c3       = Huge(ONE)
     c4       = Huge(ONE)
     c5       = Huge(ONE)
     h1       = Huge(ONE)
     k1       = Huge(ONE)

     d1       = Huge(ONE)
     d4       = Huge(ONE)
     hd1      = Huge(ONE)
     hd2      = Huge(ONE)
     hc1      = Huge(ONE)
     hc2      = Huge(ONE)
     k3       = Huge(ONE)


    read(nio,PSB3,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unitp,*) ' ERROR in Read_PSB3Pot'
      write(out_unitp,*) ' End-of-file or End-of-record'
      write(out_unitp,*) ' The namelist "PSB3" is probably absent'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,*)
      STOP ' ERROR in Read_PSB3Pot'
    ELSE IF (err_read > 0) THEN
      write(out_unitp,*) ' ERROR in Read_PSB3Pot'
      write(out_unitp,*) ' Some parameter names of the namelist "PSB3" are probaly wrong'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,nml=PSB3)
      STOP ' ERROR in Read_PSB3Pot'
    END IF

      Para_PSB3%blamin   = blamin
      Para_PSB3%blaTSdir = blaTSdir
      Para_PSB3%deepth   = deepth
      Para_PSB3%kf1      = kf1
      Para_PSB3%d2       = d2
      Para_PSB3%d3       = d3
      Para_PSB3%kf4      = kf4
      Para_PSB3%c1       = c1
      Para_PSB3%c2       = c2
      Para_PSB3%c3       = c3
      Para_PSB3%c4       = c4
      Para_PSB3%c5       = c5
      Para_PSB3%h1       = h1
      Para_PSB3%k1       = k1

      Para_PSB3%d1       = d1
      Para_PSB3%d4       = d4
      Para_PSB3%hd1      = hd1
      Para_PSB3%hd2      = hd2
      Para_PSB3%hc1      = hc1
      Para_PSB3%hc2      = hc2
      Para_PSB3%k3       = k3

  END SUBROUTINE Read_PSB3Pot
!> @brief Subroutine wich prints the PSB3 parameters.
!!
!! @param nio                integer:             file unit to print the parameters.
  SUBROUTINE Write0_PSB3Pot(nio)

    integer         , intent(in) :: nio

    write(nio,*) 'PSB3 default parameters'
    write(nio,*)
    write(nio,*) ' Warning the parameters are given as in the publication.'
    write(nio,*) '  Therefore, the BLA(=Q(1)) is in Angstrom and the energy is in kcal.mol^-1.'
    write(nio,*)
    write(nio,*) 'end PSB3 default parameters'

  END SUBROUTINE Write0_PSB3Pot

!> @brief Subroutine wich prints the PSB3 parameters.
!!
!! @param Para_PSB3         TYPE(Param_PSB3):   derived type in which the parameters are set-up.
!! @param nio                integer:             file unit to print the parameters.
  SUBROUTINE Write_PSB3Pot(Para_PSB3,nio)
    
    TYPE(Param_PSB3), intent(in) :: Para_PSB3
    integer         , intent(in) :: nio

    write(nio,*) 'PSB3 current parameters'
    write(nio,*)
    write(nio,*) '  PubliUnit:      ',Para_PSB3%PubliUnit
    write(nio,*)
    write(nio,*) '  Option   :      ',Para_PSB3%option
    write(nio,*)

    SELECT CASE (Para_PSB3%option)

    CASE (1)

    write(nio,*) '  ... for first Model:'
    write(nio,*) '  blamin   :      ',Para_PSB3%blamin
    write(nio,*) '  blaTSdir :      ',Para_PSB3%blaTSdir
    write(nio,*) '  deepth   :      ',Para_PSB3%deepth
    write(nio,*) '  kf1      :      ',Para_PSB3%kf1
    write(nio,*) '  d2       :      ',Para_PSB3%d2
    write(nio,*) '  d3       :      ',Para_PSB3%d3
    write(nio,*) '  kf4      :      ',Para_PSB3%kf4
    write(nio,*) '  c1       :      ',Para_PSB3%c1
    write(nio,*) '  c2       :      ',Para_PSB3%c2
    write(nio,*) '  c3       :      ',Para_PSB3%c3
    write(nio,*) '  c4       :      ',Para_PSB3%c4
    write(nio,*) '  c5       :      ',Para_PSB3%c5
    write(nio,*) '  h1       :      ',Para_PSB3%h1
    write(nio,*) '  k1       :      ',Para_PSB3%k1

    CASE (2)

    write(nio,*) '  ... for second Model:'
    write(nio,*) '  blamin   :      ',Para_PSB3%blamin
    write(nio,*) '  d1       :      ',Para_PSB3%d1
    write(nio,*) '  d2       :      ',Para_PSB3%d2
    write(nio,*) '  d3       :      ',Para_PSB3%d3
    write(nio,*) '  d4       :      ',Para_PSB3%d4
    write(nio,*) '  c1       :      ',Para_PSB3%c1
    write(nio,*) '  c2       :      ',Para_PSB3%c2
    write(nio,*) '  c3       :      ',Para_PSB3%c3
    write(nio,*) '  c4       :      ',Para_PSB3%c4
    write(nio,*) '  c5       :      ',Para_PSB3%c5
    write(nio,*) '  hd1      :      ',Para_PSB3%hd1
    write(nio,*) '  hd2      :      ',Para_PSB3%hd2
    write(nio,*) '  hc1      :      ',Para_PSB3%hc1
    write(nio,*) '  hc2      :      ',Para_PSB3%hc2
    write(nio,*) '  k1       :      ',Para_PSB3%k1
    write(nio,*) '  k3       :      ',Para_PSB3%k3

    CASE Default
        write(out_unitp,*) ' ERROR in write_PSB3Pot '
        write(out_unitp,*) ' This option is not possible. option: ',Para_PSB3%option
        write(out_unitp,*) ' Its value MUST be 1 or 2 '

        STOP
    END SELECT

    write(nio,*)
    write(nio,*) 'end PSB3 current parameters'

  END SUBROUTINE Write_PSB3Pot

!> @brief Subroutine wich calculates the PSB3 potential (for the 3 models) with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE(dnMatPot):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param Para_PSB3          TYPE(Param_PSB3):    derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_PSB3Pot(Mat_OF_PotDia,dnQ,Para_PSB3,nderiv)
    USE mod_dnS

    TYPE(Param_PSB3) , intent(in)    :: Para_PSB3
    TYPE(dnS),       intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE(dnS),       intent(in)    :: dnQ(:) ! BLA, Tors, HOOP
    integer          , intent(in)    :: nderiv

    SELECT CASE (Para_PSB3%option)

    CASE (1)
      CALL eval_PSB3Pot1(Mat_OF_PotDia,dnQ,Para_PSB3,nderiv)

    CASE (2)
      CALL eval_PSB3Pot2(Mat_OF_PotDia,dnQ,Para_PSB3,nderiv)

    CASE Default
        write(out_unitp,*) ' ERROR in eval_PSB3Pot '
        write(out_unitp,*) ' This option is not possible. option: ',Para_PSB3%option
        write(out_unitp,*) ' Its value MUST be 1 or 2 '

        STOP
    END SELECT

  END SUBROUTINE eval_PSB3Pot

!> @brief Subroutine wich calculates the PSB3 potential (Not published model) with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE(dnMatPot):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param Para_PSB3          TYPE(Param_PSB3):    derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).

  SUBROUTINE eval_PSB3Pot1(Mat_OF_PotDia,dnQ,Para_PSB3,nderiv)
    !Not Published model potential 
    USE mod_dnS

    real (kind=Rkind),  parameter     :: EnergyConv = 627.509_Rkind,LenghtConv = 0.52917721067121_Rkind
    TYPE(dnS),        intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE(dnS),        intent(in)    :: dnQ(:) ! BLA, Tors, HOOP
    TYPE(Param_PSB3) ,  intent(in)    :: Para_PSB3
    integer,            intent(in)    :: nderiv


    !local variables (derived type). They have to be deallocated
    TYPE(dnS)        :: dnPot,BLA,Tors,HOOP
    TYPE(dnS)        :: MorseBLAP,Morsemin,Overlap
    TYPE(dnS)        :: Hdir2D,Hct2D

    real (kind=Rkind)  :: d1,d4



    BLA      = dnQ(1)
    Tors     = dnQ(2)
    HOOP     = dnQ(3)

    !If PubliUnit is False conversion must be done, the potential is expressed in Angstrom 
    !and it requires the proper conversion into Bhor  
    IF(.NOT. Para_PSB3%PubliUnit) THEN
       BLA = BLA * LenghtConv
    END IF

!-----------------------------------------------------------------------!
    d1 = SQRT(Para_PSB3%kf1/(TWO * Para_PSB3%deepth))
    d4 = SQRT(Para_PSB3%kf4/(TWO * Para_PSB3%deepth))

    MorseBLAP = Para_PSB3%deepth * (-ONE + EXP(-d1 * (BLA - Para_PSB3%blaTSdir)))**2
    Morsemin  = Para_PSB3%deepth * (-ONE + EXP(-d4 * (BLA - Para_PSB3%blamin  )))**2

    Overlap   = Tors - HOOP/TWO

    Hdir2D = sin(Overlap)**2 * (MorseBLAP + Para_PSB3%d2) + Para_PSB3%d3 * cos(Overlap / TWO)**2 + Morsemin * Cos(Overlap)**2
    Hct2D  = (ONE + Para_PSB3%c5 * sin(Overlap)**2) * (Para_PSB3%c1 * BLA**2 + &
             Para_PSB3%c2 * BLA + Para_PSB3%c3) + Para_PSB3%c4 * cos(Overlap)**2
!-----------------------------------------------------------------------!

!-----------------------------------------------------------------------!
! V11 matrix element 

   dnPot =  Hdir2D + Para_PSB3%h1 * sin(HOOP/FOUR)**2

   IF(.NOT. Para_PSB3%PubliUnit) dnPot = dnPot/EnergyConv

   Mat_OF_PotDia(1,1) = dnPot
!-----------------------------------------------------------------------!
! V22 matrix element 

   dnPot =  Hct2D + Para_PSB3%h1 * sin(HOOP/FOUR)**2

   IF(.NOT. Para_PSB3%PubliUnit) dnPot = dnPot/EnergyConv

   Mat_OF_PotDia(2,2) = dnPot

!-----------------------------------------------------------------------!
!V12 = V21

   dnPot = Para_PSB3%k1 * sin((Overlap) * TWO)

   IF(.NOT. Para_PSB3%PubliUnit) dnPot = dnPot/EnergyConv

   Mat_OF_PotDia(1,2) = dnPot
   Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)


!-----------------------------------------------------------------------!
   CALL dealloc_dnS(BLA)
   CALL dealloc_dnS(Tors)
   CALL dealloc_dnS(HOOP)

   CALL dealloc_dnS(Overlap)
   CALL dealloc_dnS(MorseBLAP)
   CALL dealloc_dnS(Morsemin)
   CALL dealloc_dnS(Hdir2D)
   CALL dealloc_dnS(Hct2D)

   CALL dealloc_dnS(dnPot)

  END SUBROUTINE eval_PSB3Pot1

!> @brief Subroutine wich calculates the PSB3 potential (Published model) with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE(dnMatPot):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param Para_PSB3          TYPE(Param_PSB3):    derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).

  SUBROUTINE eval_PSB3Pot2(Mat_OF_PotDia,dnQ,Para_PSB3,nderiv) !Second PSB3's potential
  ! Published potential 
    USE mod_dnS

    real (kind=Rkind),  parameter      :: EnergyConv = 627.509_Rkind,LenghtConv = 0.52917721067121_Rkind 
    TYPE (Param_PSB3),  intent(in)     :: Para_PSB3
    TYPE(dnS),        intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE(dnS),        intent(in)     :: dnQ(:) ! BLA, Tors, HOOP
    integer, intent(in)                :: nderiv
    real (kind=Rkind)                  :: d1,d4

    !local variables (derived type). They have to be deallocated
    TYPE(dnS)     :: dnPot,BLA,Tors,HOOP
    TYPE(dnS)     :: Overlap
    TYPE(dnS)     :: Hdir2D,Hct2D



    BLA      = dnQ(1)
    Tors     = dnQ(2)
    HOOP     = dnQ(3)


    !If PubliUnit is False conversion must be done, the potential is expressed in Angstrom 
    !and it requires the proper conversion into Bhor  
    IF(.NOT. Para_PSB3%PubliUnit) THEN
       !BLA = BLA * LenghtConv            ! wrong derivative. Here with respect ot BLA
       BLA = d0S_TIME_R(BLA,LenghtConv) ! to set up the correct derivatives with respect to (R*LenghtConv)
    END IF

!-----------------------------------------------------------------------!
    Overlap   = Tors - HOOP/TWO

    Hdir2D = SIN(Tors)**2 * (Para_PSB3%d1 * BLA**2 + Para_PSB3%d2) + &
             Para_PSB3%d3 * COS(Tors / TWO)**2 + Para_PSB3%d4 * (BLA - Para_PSB3%blamin)**2
    Hct2D  = (ONE + Para_PSB3%c5 * SIN(Tors)**2) *                   &
             (Para_PSB3%c1 * BLA**2 + Para_PSB3%c2 * BLA + Para_PSB3%c3) + Para_PSB3%c4 * COS(Tors)**2
!-----------------------------------------------------------------------!

!-----------------------------------------------------------------------!
! V11 matrix element 

   dnPot =  Hdir2D + Para_PSB3%hd1 * sin(HOOP/FOUR)**2 - Para_PSB3%hd2 * Sin(HOOP/FOUR) * sin(Tors * TWO)

   IF(.NOT. Para_PSB3%PubliUnit) dnPot = dnPot/EnergyConv

   Mat_OF_PotDia(1,1) = dnPot

!-----------------------------------------------------------------------!
! V22 matrix element 

   dnPot =  Hct2D  + Para_PSB3%hc1 * Sin(HOOP/FOUR)**2 + Para_PSB3%hc2 * Sin(HOOP/FOUR) * sin(Tors * TWO)

   IF(.NOT. Para_PSB3%PubliUnit) dnPot = dnPot/EnergyConv

   Mat_OF_PotDia(2,2) = dnPot

!-----------------------------------------------------------------------!
! V12 = V21

   dnPot = (ONE + Para_PSB3%k3 * SIN(Tors)**2) *  Para_PSB3%k1 * SIN((Overlap) * TWO)

   IF(.NOT. Para_PSB3%PubliUnit) dnPot = dnPot/EnergyConv

   Mat_OF_PotDia(1,2) = dnPot
   Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)


!-----------------------------------------------------------------------!
   CALL dealloc_dnS(BLA)
   CALL dealloc_dnS(Tors)
   CALL dealloc_dnS(HOOP)

   CALL dealloc_dnS(Overlap)
   CALL dealloc_dnS(Hdir2D)
   CALL dealloc_dnS(Hct2D)

   CALL dealloc_dnS(dnPot)

  END SUBROUTINE eval_PSB3Pot2

END MODULE mod_PSB3Pot
