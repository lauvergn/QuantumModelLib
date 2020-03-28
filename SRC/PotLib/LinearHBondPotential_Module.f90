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

!> @brief Module which makes the initialization, calculation of the Linear H Bond potential (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
MODULE mod_LinearHBondPot
  USE mod_NumParameters
  ! V(Q,q) = Morse1(Q/2+q,param1)+Morse2(Q/2-q,param2)+Eref2+Buckingham(Q)
  !a0=0.52917720835354106d0
  !auTOkcal=627.51
  USE mod_MorsePot
  USE mod_BuckPot

  IMPLICIT NONE

!> @brief Derived type in which the parameters of the linear H-bond potential are set-up.
!> @brief Reference: Eq 3.79 of J. Beutier, thesis.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
  TYPE Param_LinearHBond
     PRIVATE

     TYPE (Param_Morse) :: Morse1
     TYPE (Param_Morse) :: Morse2
     real (kind=Rkind)  :: Eref2  ! energy reference for the second morse (-D*epsi^2)

     TYPE (Param_Buck)  :: Buck

     logical            :: PubliUnit = .FALSE.

     real (kind=Rkind), PUBLIC  :: muQQ= 29156.946380706224_Rkind/TWO  ! reduced mass associated to QQ (O---O)
     real (kind=Rkind), PUBLIC  :: muq = 1837.1526464003414_Rkind      ! reduced mass associated to q (H atom)


  END TYPE Param_LinearHBond
  PRIVATE Read_LinearHBondPot

CONTAINS
!> @brief Subroutine which makes the initialization of the LinearHBond parameters.
!!      Those parameters cannot be modified (PRIVATE).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Para_LinearHBond   TYPE(Para_LinearHBond):   derived type in which the parameters are set-up.
!! @param nio                integer (optional):       file unit to read the parameters.
!! @param read_param         logical (optional):       when it is .TRUE., the parameters are read. Otherwise, they are initialized.
!! @param D,a,req            real (optional):          parameters for the first Morse
!! @param epsi               real (optional):          scaling parameters for the 2d Morse (using parameters of the first Morse)
!! @param Abuck,Bbuck,Cbuck  real (optional):          parameters for the Buckingham potential
!! @param PubliUnit          logical (optional):       when PubliUnit=.TRUE., the units (Angstrom and eV) are used. Default (atomic unit).
  SUBROUTINE Init_LinearHBondPot(Para_LinearHBond,nio,read_param,PubliUnit, &
                                 D,a,req,epsi,Abuck,Bbuck,Cbuck)

    TYPE (Param_LinearHBond),    intent(inout)   :: Para_LinearHBond
    integer,           optional, intent(in)      :: nio
    logical,           optional, intent(in)      :: read_param
    logical,           optional, intent(in)      :: PubliUnit
    real (kind=Rkind), optional, intent(in)      :: D,a,req,Abuck,Bbuck,Cbuck,epsi


    logical           :: read_param_loc
    real (kind=Rkind) :: D_loc,a_loc,req_loc,Abuck_loc,Bbuck_loc,Cbuck_loc,epsi_loc


    IF (present(PubliUnit)) Para_LinearHBond%PubliUnit = PubliUnit

    read_param_loc = .FALSE.
    IF (present(read_param)) read_param_loc = read_param
    IF (read_param_loc .AND. .NOT. present(nio)) THEN
       write(out_unitp,*) ' ERROR in Init_LinearHBondPot'
       write(out_unitp,*) ' read_param = t and The file unit (nio) is not present'
       write(out_unitp,*) ' => impossible to read the input file'
       STOP 'STOP in Init_LinearHBondPot: impossible to read the input file. The file unit (nio) is not present'
    END IF


    ! initalization of the default values
    D_loc     = 60._Rkind
    a_loc     = 2.52_Rkind
    req_loc   = 0.95_Rkind

    epsi_loc  = ONE

    Abuck_loc = 2.32e5_Rkind
    Bbuck_loc = 3.15_Rkind
    Cbuck_loc = 2.31e4_Rkind

    IF (read_param_loc) THEN
      CALL Read_LinearHBondPot(Para_LinearHBond,nio,                    &
                               D_loc,a_loc,req_loc,epsi_loc,            &
                               Abuck_loc,Bbuck_loc,Cbuck_loc)
    ELSE

      IF (present(D))       D_loc     = D
      IF (present(a))       a_loc     = a
      IF (present(req))     req_loc   = req

      IF (present(epsi))    epsi_loc  = epsi

      IF (present(Abuck))   Abuck_loc = Abuck
      IF (present(Bbuck))   Bbuck_loc = Bbuck
      IF (present(Cbuck))   Cbuck_loc = Cbuck

      CALL Init_MorsePot(Para_LinearHBond%Morse1,D=D_loc,            a=a_loc,         req=req_loc)
      CALL Init_MorsePot(Para_LinearHBond%Morse2,D=D_loc*epsi_loc**2,a=a_loc/epsi_loc,req=req_loc)
      CALL Init_BuckPot(Para_LinearHBond%Buck,A=Abuck_loc,B=Bbuck_loc,C=Cbuck_loc)
      Para_LinearHBond%Eref2 = -D_loc*epsi_loc**2
    END IF

    IF (Para_LinearHBond%PubliUnit) THEN
      write(out_unitp,*) 'PubliUnit=.TRUE.,  Q:[Angs,Angs], Energy: [kcal.mol^-1]'
    ELSE
      write(out_unitp,*) 'PubliUnit=.FALSE., Q:[Bohr,Bohr], Energy: [Hartree]'
    END IF

  END SUBROUTINE Init_LinearHBondPot
!> @brief Subroutine wich reads the LinearHBond parameters with a namelist.
!!   This can be called only from the "Init_LinearHBondPot" subroutine.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Para_LinearHBond            TYPE(Param_LinearHBond):    derived type in which the parameters are set-up.
!! @param nio                         integer:                    file unit to read the parameters.
!! @param Dsub,asub,reqsub            real (optional):            parameters for the first Morse
!! @param epsisub                     real (optional):            scaling parameters for the 2d Morse (using parameters of the first Morse)
!! @param Abucksub,Bbucksub,Cbucksub  real (optional):            parameters for the Buckingham potential
  SUBROUTINE Read_LinearHBondPot(Para_LinearHBond,nio,                  &
                                 Dsub,asub,reqsub,epsisub,              &
                                 Abucksub,Bbucksub,Cbucksub)

    TYPE (Param_LinearHBond), intent(inout) :: Para_LinearHBond
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
      write(out_unitp,*) ' ERROR in Read_LinearHBondPot'
      write(out_unitp,*) ' End-of-file or End-of-record'
      write(out_unitp,*) ' The namelist "LinearHBond" is probably absent'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,*)
      STOP ' ERROR in Read_LinearHBondPot'
    ELSE IF (err_read > 0) THEN
      write(out_unitp,*) ' ERROR in Read_LinearHBondPot'
      write(out_unitp,*) ' Some parameter names of the namelist "LinearHBond" are probaly wrong'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,nml=LinearHBond)
      STOP ' ERROR in Read_LinearHBondPot'
    END IF
    write(out_unitp,nml=LinearHBond)

    ! initalization with the read values
    CALL Init_MorsePot(Para_LinearHBond%Morse1,D=D,        a=a,     req=req)
    CALL Init_MorsePot(Para_LinearHBond%Morse2,D=D*epsi**2,a=a/epsi,req=req)
    Para_LinearHBond%Eref2 = -D*epsi**2
    CALL Init_BuckPot(Para_LinearHBond%Buck,A=Abuck,B=Bbuck,C=Cbuck)

  END SUBROUTINE Read_LinearHBondPot
!> @brief Subroutine wich prints the LinearHBond potential parameters from Beutier PhD thesis.
!!
!> @author David Lauvergnat
!! @date 30/07/2019
!!
!! @param nio                integer:                   file unit to print the parameters.
  SUBROUTINE Write0_LinearHBondPot(nio)
    integer, intent(in) :: nio

    write(nio,*) 'LinearHBond parameters:'
    write(nio,*)
    write(nio,*) ' Potential and parameters from J. Beutier PhD thesis (Eq 3.79)'
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
    write(nio,*) 'END LinearHBond parameters'

  END SUBROUTINE Write0_LinearHBondPot
!> @brief Subroutine wich prints the LinearHBond potential parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Para_LinearHBond   TYPE(Param_LinearHBond):   derived type with the Phenol potential parameters.
!! @param nio                integer:                   file unit to print the parameters.
  SUBROUTINE Write_LinearHBondPot(Para_LinearHBond,nio)
    TYPE (Param_LinearHBond), intent(in) :: Para_LinearHBond
    integer, intent(in) :: nio

    write(nio,*) 'LinearHBond current parameters:'
    write(nio,*)
    write(nio,*) 'PubliUnit: ',Para_LinearHBond%PubliUnit
    write(nio,*)
    write(nio,*) '   Morse parameters:   '
    CALL Write_MorsePot(Para_LinearHBond%Morse1,nio)
    CALL Write_MorsePot(Para_LinearHBond%Morse2,nio)
    write(nio,*) '   Buckingham parameters:   '
    CALL Write_BuckPot(Para_LinearHBond%Buck,nio)
    write(nio,*) '  Eref2 = ',Para_LinearHBond%Eref2

    write(nio,*) 'END LinearHBond current parameters'

  END SUBROUTINE Write_LinearHBondPot

!> @brief Subroutine wich calculates the LinearHBond potential with derivatives up to the 2d order if required.
!!
!> @author David Lauvergnat
!! @date 25/05/2018
!!
!! @param PotVal             TYPE(dnMatPot):          derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param Q                  real:                    table of two values for which the potential is calculated (R,theta)
!! @param Para_LinearHBond   TYPE(Param_LinearHBond): derived type with the Morse parameters.
!! @param nderiv             integer:                 it enables to specify up to which derivatives the potential is calculated:
!!                                                    the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
   SUBROUTINE Eval_LinearHBondPot(Mat_OF_PotDia,dnQ,Para_LinearHBond,nderiv)
    USE mod_dnS

    TYPE (Param_LinearHBond),  intent(in)     :: Para_LinearHBond
    TYPE(dnS),               intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE(dnS),               intent(in)     :: dnQ(:) ! QQ and q
    integer,                   intent(in)     :: nderiv

    !local variable
    TYPE(dnS)          :: PotVal_m1,PotVal_m2,PotVal_Buck
    TYPE(dnS)          :: dnQQ,dnsq,dnX,dnY
    real (kind=Rkind), parameter  :: a0      = 0.52917720835354106_Rkind
    real (kind=Rkind), parameter  :: auTOkcalmol_inv  = 627.51_Rkind

    !logical, parameter :: debug=.TRUE.
    logical, parameter :: debug=.FALSE.

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING Eval_LinearHBondPot'
      write(out_unitp,*) 'r(:) or QQ,q: ',get_d0_FROM_dnS(dnQ(1)),get_d0_FROM_dnS(dnQ(2))
      write(out_unitp,*) 'nderiv',nderiv
      flush(out_unitp)
    END IF

    IF (debug) THEN
      write(out_unitp,*) 'dnQQ'
      CALL write_dnS(dnQ(1))
      write(out_unitp,*) 'dnsq'
      CALL write_dnS(dnQ(2))
      flush(out_unitp)
    END IF

    IF (.NOT. Para_LinearHBond%PubliUnit) THEN
       dnQQ = a0*dnQ(1) ! to convert the bhor into Angstrom.
       dnsq = a0*dnQ(2) ! to convert the bhor into Angstrom.

      IF (debug) THEN
        write(out_unitp,*) 'dnQQ in Angs'
        CALL write_dnS(dnQ(1))
        write(out_unitp,*) 'dnsq in Angs'
        CALL write_dnS(dnQ(2))
        flush(out_unitp)
      END IF
    ELSE
      dnQQ = dnQ(1)
      dnsq = dnQ(2)
    END IF


    ! new variables for the Morse potentials
    dnX = dnQQ/TWO + dnsq
    dnY = dnQQ/TWO - dnsq

    IF (debug) THEN
      write(out_unitp,*) 'dnX'
      flush(out_unitp)
      CALL write_dnS(dnX)
      write(out_unitp,*) 'dnY'
      flush(out_unitp)
      CALL write_dnS(dnY)
      flush(out_unitp)
    END IF

    PotVal_m1 = dnMorse(dnX,Para_LinearHBond%Morse1)
    IF (debug) THEN
      write(out_unitp,*) 'PotVal_m1. x:',get_d0_FROM_dnS(dnX)
      CALL Write_dnS(PotVal_m1)
      flush(out_unitp)
    END IF

    PotVal_m2 = dnMorse(dnY,Para_LinearHBond%Morse2)+Para_LinearHBond%Eref2
    IF (debug) THEN
      write(out_unitp,*) 'PotVal_m2. y:',get_d0_FROM_dnS(dnY)
      CALL Write_dnS(PotVal_m2)
      flush(out_unitp)
    END IF

    PotVal_Buck = dnBuck(dnQQ,Para_LinearHBond%Buck)
    IF (debug) THEN
      write(out_unitp,*) 'PotVal_Buck. QQ:',get_d0_FROM_dnS(dnQQ)
      CALL Write_dnS(PotVal_Buck)
      flush(out_unitp)
    END IF

    Mat_OF_PotDia(1,1) = PotVal_m1+PotVal_m2+PotVal_Buck

    IF (debug) THEN
      write(out_unitp,*) 'Mat_OF_PotDia(1,1):'
      CALL write_dnS(Mat_OF_PotDia(1,1))
      flush(out_unitp)
    END IF

    IF (.NOT. Para_LinearHBond%PubliUnit) THEN
      Mat_OF_PotDia(1,1) = Mat_OF_PotDia(1,1)/auTOkcalmol_inv ! to convert the kcal/mol into Hartree
    END IF

    CALL dealloc_dnS(dnX)
    CALL dealloc_dnS(dnY)
    CALL dealloc_dnS(dnQQ)
    CALL dealloc_dnS(dnsq)

    CALL dealloc_dnS(PotVal_m1)
    CALL dealloc_dnS(PotVal_m2)
    CALL dealloc_dnS(PotVal_Buck)

    IF (debug) THEN
      write(out_unitp,*) 'END Eval_LinearHBondPot'
      flush(out_unitp)
    END IF

  END SUBROUTINE Eval_LinearHBondPot

END MODULE mod_LinearHBondPot
