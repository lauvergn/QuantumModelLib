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
!    Copyright 2016 David Lauvergnat [1]
!      with contributions of:
!        Félix MOUHAT [2]
!        Liang LIANG [3]
!        Emanuele MARSILI [1,4]
!
![1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
![2]: Laboratoire PASTEUR, ENS-PSL-Sorbonne Université-CNRS, France
![3]: Maison de la Simulation, CEA-CNRS-Université Paris-Saclay,France
![4]: Durham University, Durham, UK
!* Originally, it has been developed during the Quantum-Dynamics E-CAM project :
!     https://www.e-cam2020.eu/quantum-dynamics
!
!===========================================================================
!===========================================================================

!> @brief Module which makes the initialization, calculation of the LinearHBond potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE mod_LinearHBondModel
  USE mod_NumParameters
  USE mod_EmptyModel
  USE mod_MorseModel
  USE mod_BuckModel
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the LinearHBond parameters are set-up.
  TYPE, EXTENDS (EmptyModel_t) :: LinearHBondModel_t
   PRIVATE

     TYPE (MorseModel_t)  :: Morse1
     TYPE (MorseModel_t)  :: Morse2
     real (kind=Rkind)    :: Eref2 = ZERO  ! energy reference for the second morse (-D*epsi^2)

     TYPE (BuckModel_t)   :: Buck

     real (kind=Rkind)    :: muQQ= 29156.946380706224_Rkind/TWO  ! reduced mass associated to QQ (O---O)
     real (kind=Rkind)    :: muq = 1837.1526464003414_Rkind      ! reduced mass associated to q (H atom)

   CONTAINS
    PROCEDURE :: Eval_QModel_Pot => eval_LinearHBondPot
    PROCEDURE :: Write_QModel    => Write_LinearHBondModel
    PROCEDURE :: Write0_QModel   => Write0_LinearHBondModel
  END TYPE LinearHBondModel_t

  PUBLIC :: LinearHBondModel_t,Init_LinearHBondModel

  CONTAINS
!> @brief Function which makes the initialization of the LinearHBond parameters.
!!when PubliUnit=.TRUE., the units (Angstrom and eV) are used. Default (atomic unit).
!!
!! @param QModel             TYPE(LinearHBondModel_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(EmptyModel_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
!! @param D,a,req            real (optional):     parameters for the first Morse
!! @param epsi               real (optional):     scaling parameters for the 2d Morse (using parameters of the first Morse)
!! @param Abuck,Bbuck,Cbuck  real (optional):     parameters for the Buckingham potential
  FUNCTION Init_LinearHBondModel(QModel_in,read_param,nio_param_file,   &
                                 D,a,req,epsi,Abuck,Bbuck,Cbuck) RESULT(QModel)
  IMPLICIT NONE

    TYPE (LinearHBondModel_t)                    :: QModel ! RESULT

    TYPE(EmptyModel_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param
    real (kind=Rkind), optional, intent(in)      :: D,a,req,Abuck,Bbuck,Cbuck,epsi


    real (kind=Rkind) :: D_loc,a_loc,req_loc,Abuck_loc,Bbuck_loc,Cbuck_loc,epsi_loc
    real (kind=Rkind), parameter  :: a0               = 0.52917720835354106_Rkind
    real (kind=Rkind), parameter  :: auTOkcalmol_inv  = 627.51_Rkind
    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_LinearHBondModel'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL Init0_EmptyModel(QModel%EmptyModel_t,QModel_in)

    QModel%nsurf    = 1
    QModel%ndim     = 2
    QModel%pot_name = 'hbond'

    ! initalization of the default values
    D_loc     = 60._Rkind
    a_loc     = 2.52_Rkind
    req_loc   = 0.95_Rkind

    epsi_loc  = ONE

    Abuck_loc = 2.32e5_Rkind
    Bbuck_loc = 3.15_Rkind
    Cbuck_loc = 2.31e4_Rkind

    IF (read_param) THEN
      CALL Read_LinearHBondModel(nio_param_file,                                &
                                 D_loc,a_loc,req_loc,epsi_loc,                  &
                                 Abuck_loc,Bbuck_loc,Cbuck_loc)
    ELSE

      IF (present(D))       D_loc     = D
      IF (present(a))       a_loc     = a
      IF (present(req))     req_loc   = req

      IF (present(epsi))    epsi_loc  = epsi

      IF (present(Abuck))   Abuck_loc = Abuck
      IF (present(Bbuck))   Bbuck_loc = Bbuck
      IF (present(Cbuck))   Cbuck_loc = Cbuck
    END IF

    CALL Init0_MorseModel(QModel%Morse1,D=D_loc,            a=a_loc,         req=req_loc,model_name='Morse1')
    CALL Init0_MorseModel(QModel%Morse2,D=D_loc*epsi_loc**2,a=a_loc/epsi_loc,req=req_loc,model_name='Morse2')
    CALL Init0_BuckModel(QModel%Buck,A=Abuck_loc,B=Bbuck_loc,C=Cbuck_loc,model_name='Buck')
    QModel%Eref2 = -D_loc*epsi_loc**2

    IF (QModel%PubliUnit) THEN
      write(out_unitp,*) 'PubliUnit=.TRUE.,  Q:[Angs,Angs], Energy: [kcal.mol^-1]'
    ELSE
      write(out_unitp,*) 'PubliUnit=.FALSE., Q:[Bohr,Bohr], Energy: [Hartree]'
    END IF


    IF (debug) write(out_unitp,*) 'init Q0 of LinearHBond'
    QModel%Q0 = [2.75_Rkind,ZERO] ! in Angstrom
    IF (.NOT. QModel%PubliUnit) QModel%Q0 = QModel%Q0/a0    ! in Bohr

    IF (debug) write(out_unitp,*) 'init d0GGdef of LinearHBond'
    CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)
    QModel%d0GGdef(1,1) = ONE/QModel%muQQ
    QModel%d0GGdef(2,2) = ONE/QModel%muq

    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_LinearHBondModel
!> @brief Subroutine wich reads the model parameters with a namelist.
!!   This can be called only from the "Init_LinearHBondModel" subroutine.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param nio                         integer:                    file unit to read the parameters.
!! @param Dsub,asub,reqsub            real (optional):            parameters for the first Morse
!! @param epsisub                     real (optional):            scaling parameters for the 2d Morse (using parameters of the first Morse)
!! @param Abucksub,Bbucksub,Cbucksub  real (optional):            parameters for the Buckingham potential
  SUBROUTINE Read_LinearHBondModel(nio,Dsub,asub,reqsub,epsisub,        &
                                   Abucksub,Bbucksub,Cbucksub)

    real (kind=Rkind),        intent(inout) :: Dsub,asub,reqsub,epsisub
    real (kind=Rkind),        intent(inout) :: Abucksub,Bbucksub,Cbucksub
    integer,                  intent(in)    :: nio

    real (kind=Rkind) :: D,a,req,epsi,Abuck,Bbuck,Cbuck ! for the namelist
    integer           :: err_read

    namelist /LinearHBond/ D,a,req,epsi,Abuck,Bbuck,Cbuck

    ! to recover the default value
    D     = Dsub
    a     = asub
    req   = reqsub
    epsi  = epsisub
    Abuck = Abucksub
    Bbuck = Bbucksub
    Cbuck = Cbucksub

    read(nio,nml=LinearHBond,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unitp,*) ' ERROR in Read_LinearHBondModel'
      write(out_unitp,*) ' End-of-file or End-of-record'
      write(out_unitp,*) ' The namelist "LinearHBond" is probably absent'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,*)
      STOP ' ERROR in Read_LinearHBondPot'
    ELSE IF (err_read > 0) THEN
      write(out_unitp,*) ' ERROR in Read_LinearHBondModel'
      write(out_unitp,*) ' Some parameter names of the namelist "LinearHBond" are probaly wrong'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,nml=LinearHBond)
      STOP ' ERROR in Read_LinearHBondModel'
    END IF
    !write(out_unitp,nml=LinearHBond)

    Dsub     = D
    asub     = a
    reqsub   = req
    epsisub  = epsi
    Abucksub = Abuck
    Bbucksub = Bbuck
    Cbucksub = Cbuck

  END SUBROUTINE Read_LinearHBondModel
!> @brief Subroutine wich prints the current LinearHBondModel parameters.
!!
!! @param QModel            CLASS(LinearHBondModel_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_LinearHBondModel(QModel,nio)
  IMPLICIT NONE

    CLASS(LinearHBondModel_t),   intent(in) :: QModel
    integer,                     intent(in) :: nio

    write(nio,*) 'LinearHBond current parameters:'
    write(nio,*)
    write(nio,*) 'PubliUnit: ',QModel%PubliUnit
    write(nio,*)
    write(nio,*) '   Morse parameters:   '
    CALL Write_MorseModel(QModel%Morse1,nio)
    CALL Write_MorseModel(QModel%Morse2,nio)
    write(nio,*) '   Buckingham parameters:   '
    CALL Write_BuckModel(QModel%Buck,nio)
    write(nio,*) '  Eref2 = ',QModel%Eref2

    write(nio,*) 'END LinearHBond current parameters'

  END SUBROUTINE Write_LinearHBondModel
!> @brief Subroutine wich prints the default LinearHBondModel parameters.
!!
!! @param QModel            CLASS(LinearHBondModel_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write0_LinearHBondModel(QModel,nio)
  IMPLICIT NONE

    CLASS(LinearHBondModel_t),   intent(in) :: QModel
    integer,                     intent(in) :: nio

    write(nio,*) 'LinearHBond default parameters:'
    write(nio,*)
    write(nio,*) ' Potential and parameters from:'
    write(nio,*) '  - Dana Codruta Marinica, Marie-Pierre Gaigeot, Daniel Borgis, '
    write(nio,*) '    Chemical Physics Letters 423 (2006) 390–394'
    write(nio,*) '    DOI: 10.1016/j.cplett.2006.04.007'
    write(nio,*) ' -  J. Beutier PhD thesis (Eq 3.79)'
    write(nio,*) '    Hydrogen motion between two atoms, A and B.'
    write(nio,*)
    write(nio,*) '    A--------------H-----X--------------------B'
    write(nio,*) '     <--------------QQ----------------------->'
    write(nio,*) '                    <-q->'
    write(nio,*)
    write(nio,*) '  V(Q,q) = Morse1(Q/2+q,param1)+Morse2(Q/2-q,param2)+Eref2+Buckingham(Q)'
    write(nio,*) '    Units: QQ and q in Angstrom, energy in kcal.mol^-1'
    write(nio,*)
    write(nio,*) '  Morse1 (D1*(1-exp(-a1*(r-Req1))**2) parameters:'
    write(nio,*) '    D1   = 60. kcal.mol^-1'
    write(nio,*) '    a1   = 2.52 Angs^-2'
    write(nio,*) '    Req1 = 0.95 Angs'
    write(nio,*)
    write(nio,*) '  Morse2 (D2*(1-exp(-a2*(r-Req2))**2) parameters:'
    write(nio,*) '    D2   = D1*epsi^2'
    write(nio,*) '    a2   = a1/epsi'
    write(nio,*) '    Req2 = Req1'
    write(nio,*) '       with espi=1'
    write(nio,*)
    write(nio,*) '  Buckingham (A.Exp(-B*r)-C/r^6) parameters:'
    write(nio,*) '    A   = 2.32 10^5 kcal.mol^-1'
    write(nio,*) '    B   = 3.15      Angs^-1'
    write(nio,*) '    C   = 2.31 10^4 kcal.mol^-1.Angs^-6'
    write(nio,*)
    write(nio,*) '  Eref2 = -D1*epsi^2'
    write(nio,*)
    write(nio,*)
    write(nio,*) 'Potential Value at: QQ=2.75 Angs and q=0.0 Angs'
    write(nio,*) 'V = -21.433838 kcal.mol^-1'
    write(nio,*) 'gradient (kcal.mol^-1 Angs^-1) = [ 58.250585, 0.000000]'
    write(nio,*)
    write(nio,*) 'end LinearHBond default parameters'


  END SUBROUTINE Write0_LinearHBondModel

!> @brief Subroutine wich calculates the LinearHBond potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(LinearHBondModel_t):   derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_LinearHBondPot(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE mod_dnS
  IMPLICIT NONE

    CLASS(LinearHBondModel_t),  intent(in)    :: QModel
    TYPE (dnS_t),               intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),               intent(in)    :: dnQ(:)
    integer,                    intent(in)    :: nderiv

    !local variable
    TYPE (dnS_t)                  :: PotVal_m1,PotVal_m2,PotVal_Buck
    TYPE (dnS_t)                  :: dnQQ,dnsq,dnX,dnY
    real (kind=Rkind), parameter  :: a0               = 0.52917720835354106_Rkind
    real (kind=Rkind), parameter  :: auTOkcalmol_inv  = 627.51_Rkind

    integer :: n1,n2

    !logical, parameter :: debug=.TRUE.
    logical, parameter :: debug=.FALSE.
    IF (debug .OR. print_level > 0) THEN
      write(out_unitp,*) 'BEGINNING Eval_LinearHBondPot'
      write(out_unitp,*) 'r(:) or QQ,q: ',QML_get_d0_FROM_dnS(dnQ(:))
      write(out_unitp,*) 'nderiv',nderiv
      write(out_unitp,*) 'PubliUnit',QModel%PubliUnit
      flush(out_unitp)
    END IF

    IF (debug .OR. print_level > 1) THEN
      write(out_unitp,*) 'dnQQ'
      CALL QML_Write_dnS(dnQ(1),all_type=.TRUE.)
      write(out_unitp,*) 'dnsq'
      CALL QML_Write_dnS(dnQ(2),all_type=.TRUE.)
      flush(out_unitp)
    END IF

    IF (.NOT. QModel%PubliUnit) THEN
       dnQQ = a0*dnQ(1) ! to convert the bhor into Angstrom.
       dnsq = a0*dnQ(2) ! to convert the bhor into Angstrom.

      IF (debug .OR. print_level > 1) THEN
        write(out_unitp,*) 'dnQQ in Angs'
        CALL QML_Write_dnS(dnQQ,all_type=.TRUE.)
        write(out_unitp,*) 'dnsq in Angs'
        CALL QML_Write_dnS(dnsq,all_type=.TRUE.)
        flush(out_unitp)
      END IF

    ELSE
      dnQQ = dnQ(1)
      dnsq = dnQ(2)
    END IF


    ! new variables for the Morse potentials
    dnX = dnQQ/TWO + dnsq
    dnY = dnQQ/TWO - dnsq

    IF (debug .OR. print_level > 1) THEN
      write(out_unitp,*) 'dnX'
      flush(out_unitp)
      CALL QML_Write_dnS(dnX)
      write(out_unitp,*) 'dnY'
      flush(out_unitp)
      CALL QML_Write_dnS(dnY)
      flush(out_unitp)
    END IF

    PotVal_m1 = dnMorse(dnX,QModel%Morse1)
    IF (debug .OR. print_level > 1) THEN
      write(out_unitp,*) 'PotVal_m1. x:',QML_get_d0_FROM_dnS(dnX)
      CALL QML_Write_dnS(PotVal_m1)
      flush(out_unitp)
    END IF

    PotVal_m2 = dnMorse(dnY,QModel%Morse2)+QModel%Eref2
    IF (debug .OR. print_level > 1) THEN
      write(out_unitp,*) 'PotVal_m2. y:',QML_get_d0_FROM_dnS(dnY)
      CALL QML_Write_dnS(PotVal_m2)
      flush(out_unitp)
    END IF

    PotVal_Buck = dnBuck(dnQQ,QModel%Buck)
    IF (debug .OR. print_level > 1) THEN
      write(out_unitp,*) 'PotVal_Buck. QQ:',QML_get_d0_FROM_dnS(dnQQ)
      CALL QML_Write_dnS(PotVal_Buck)
      flush(out_unitp)
    END IF

    Mat_OF_PotDia(1,1) = PotVal_m1 + PotVal_m2 + PotVal_Buck

    IF (debug .OR. print_level > 1) THEN
      write(out_unitp,*) 'Mat_OF_PotDia(1,1):'
      CALL QML_Write_dnS(Mat_OF_PotDia(1,1))
      flush(out_unitp)
    END IF

    IF (.NOT. QModel%PubliUnit) THEN
      Mat_OF_PotDia(1,1) = Mat_OF_PotDia(1,1)/auTOkcalmol_inv ! to convert the kcal/mol into Hartree
    END IF

    CALL QML_dealloc_dnS(dnX)
    CALL QML_dealloc_dnS(dnY)
    CALL QML_dealloc_dnS(dnQQ)
    CALL QML_dealloc_dnS(dnsq)

    CALL QML_dealloc_dnS(PotVal_m1)
    CALL QML_dealloc_dnS(PotVal_m2)
    CALL QML_dealloc_dnS(PotVal_Buck)

    IF (debug .OR. print_level > 0) THEN
      write(out_unitp,*) 'Mat_OF_PotDia(1,1):'
      CALL QML_Write_dnS(Mat_OF_PotDia(1,1))
      flush(out_unitp)
      write(out_unitp,*) 'END Eval_LinearHBondPot'
      flush(out_unitp)
    END IF


  END SUBROUTINE eval_LinearHBondPot

END MODULE mod_LinearHBondModel
