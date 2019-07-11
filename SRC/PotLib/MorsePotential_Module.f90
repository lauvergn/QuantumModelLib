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

!> @brief Module which makes the initialization, calculation of the Morse potential (value, gradient and hessian).
!!
!> @author David Lauvergnat
!> @author Félix MOUHAT
!! @date 10/07/2019
!!
MODULE mod_MorsePot
  USE mod_NumParameters
  IMPLICIT NONE

!> @brief Derived type in which the Morse parameters are set-up.
!> @brief morse(R) = D*(1-exp(-a*(r-Req))**2
!> @brief Default parameters for H-F
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param D              real: Dissociation energy (in Hartree)
!! @param a              real: Scaling parameter (in bohr^-1)
!! @param req            real: Equilibrium distance (in bohr)
!! @param mu             real: Reduced mass of HF (in au)
  TYPE Param_Morse ! V(R) = D*(1-exp(-a*(r-Req))**2
     private
     real (kind=Rkind) :: D   = 0.225_Rkind  !< Dissociation energy for HF (in Hartree)
     real (kind=Rkind) :: a   = 1.1741_Rkind !< Scaling parameter for HF (in bohr^-1)
     real (kind=Rkind) :: req = 1.7329_Rkind !< Equilibrium HF distance (in bohr)
     real (kind=Rkind), PUBLIC :: mu  = 1744.60504565084306291455_Rkind !< Reduced mass of HF (in au)
  END TYPE Param_Morse

  PRIVATE Read_MorsePot


CONTAINS
!> @brief Subroutine which makes the initialization of the Morse parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Para_Morse         TYPE(Param_Morse):   derived type in which the parameters are set-up.
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
!! @param D                  real (optional):    Dissociation energy (in Hartree)
!! @param a                  real (optional):    Scaling parameter (in bohr^-1)
!! @param req                real (optional):    Equilibrium distance (in bohr)
  SUBROUTINE Init_MorsePot(Para_Morse,nio,read_param,D,a,req)
    TYPE (Param_Morse),          intent(inout)   :: Para_Morse
    real (kind=Rkind), optional, intent(in)      :: D,a,req
    integer,           optional, intent(in)      :: nio
    logical,           optional, intent(in)      :: read_param


    logical :: read_param_loc

    read_param_loc = .FALSE.
    IF (present(read_param)) read_param_loc = read_param
    IF (read_param_loc .AND. .NOT. present(nio)) THEN
       write(out_unitp,*) ' ERROR in Init_MorsePot'
       write(out_unitp,*) ' read_param = t and The file unit (nio) is not present'
       write(out_unitp,*) ' => impossible to read the input file'
       STOP 'ERROR in Init_MorsePot: impossible to read the input file. The file unit (nio) is not present'
    END IF

    Para_Morse = Param_Morse(D=0.225_Rkind,a=1.1741_Rkind,Req=1.7329_Rkind) ! default values (HF)

    IF (read_param_loc) THEN
      CALL Read_MorsePot(Para_Morse,nio)
    ELSE
      IF (present(D))   Para_Morse%D = D
      IF (present(a))   Para_Morse%a = a
      IF (present(req)) Para_Morse%req = req
    END IF

  END SUBROUTINE Init_MorsePot

!> @brief Subroutine wich reads the Morse parameters with a namelist.
!!   This can be called only from the "Init_MorsePot" subroutine.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Para_Morse         TYPE(Param_Morse):   derived type in which the parameters are set-up.
!! @param nio                integer          :   file unit to read the parameters.
  SUBROUTINE Read_MorsePot(Para_Morse,nio)
    TYPE (Param_Morse), intent(inout) :: Para_Morse
    integer, intent(in) :: nio

    real (kind=Rkind) :: D,a,req
    integer           :: err_read

    namelist /Morse/ D,a,req

    D   = Para_Morse%D
    a   = Para_Morse%a
    req = Para_Morse%req

    read(nio,nml=Morse,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unitp,*) ' ERROR in Read_MorsePot'
      write(out_unitp,*) ' End-of-file or End-of-record'
      write(out_unitp,*) ' The namelist "Morse" is probably absent'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,*)
      STOP ' ERROR in Read_MorsePot'
    ELSE IF (err_read > 0) THEN
      write(out_unitp,*) ' ERROR in Read_MorsePot'
      write(out_unitp,*) ' Some parameter names of the namelist "Morse" are probaly wrong'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,nml=Morse)
      STOP ' ERROR in Read_MorsePot'
    END IF
    !write(out_unitp,nml=Morse)

    Para_Morse = Param_Morse(D,a,req)

  END SUBROUTINE Read_MorsePot

!> @brief Subroutine wich prints the Morse parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Para_Morse         TYPE(Param_Morse):   derived type with the Morse parameters.
!! @param nio                integer          :   file unit to print the parameters.
  SUBROUTINE Write_MorsePot(Para_Morse,nio)
    TYPE (Param_Morse), intent(in) :: Para_Morse
    integer, intent(in) :: nio

    write(nio,*) 'Morse parameters:'
    write(nio,*) ' For H-F (Default values):'
    write(nio,*) 'Value at: R=2.0 Bohr'
    write(nio,*) 'V        = 0.016304 Hartree'
    write(nio,*) 'gradient = 0.103940'
    write(nio,*) 'hessian  = 0.209272'

    write(nio,*)
    write(nio,*) 'Current parameters:'
    write(nio,*) '    V(R) = D.( 1 - exp(-a.(r-req)) )^2'
    write(nio,*) '  D:   ',Para_Morse%D
    write(nio,*) '  a:   ',Para_Morse%a
    write(nio,*) '  req: ',Para_Morse%req
    write(nio,*) '  mu:  ',Para_Morse%mu
    write(nio,*) 'end Morse parameters'

  END SUBROUTINE Write_MorsePot

!> @brief Subroutine wich calculates the Morse potential with derivatives up to the 2d order if required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE(dnMatPot):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param Para_Morse         TYPE(Param_Morse):   derived type with the Morse parameters.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE Eval_MorsePot(PotVal,r,Para_Morse,nderiv)
    USE mod_dnMatPot
    USE mod_dnSca

    TYPE (Param_Morse), intent(in)    :: Para_Morse
    TYPE(dnMatPot),     intent(inout) :: PotVal
    real (kind=Rkind),  intent(in)    :: r
    integer,            intent(in)    :: nderiv

    !local variables (derived type). They have to be deallocated
    TYPE(dnSca)     :: dnPot,dnR

    !write(out_unitp,*) 'BEGINNING in Eval_MorsePot'
    !flush(out_unitp)

    IF ( Check_NotAlloc_dnMatPot(PotVal,nderiv) ) THEN
      CALL alloc_dnMatPot(PotVal,nsurf=1,ndim=1,nderiv=nderiv)
    END IF

    dnR     = init_dnSca(r,ndim=1,nderiv=nderiv,iQ=1) ! to set up the derivatives

    dnPot = dnMorse(dnR,Para_Morse)

    !transfert the 1D-potential and its derivatives (dnPot) to the matrix form (PotVal)
    CALL sub_dnSca_TO_dnMatPot(dnPot,PotVal,i=1,j=1)


    CALL dealloc_dnSca(dnPot)
    CALL dealloc_dnSca(dnR)

    !write(out_unitp,*) 'Morse PotVal at',r
    !CALL Write_dnMatPot(PotVal)
    !write(out_unitp,*) 'END in Eval_MorsePot'
    !flush(out_unitp)

  END SUBROUTINE Eval_MorsePot

!> @brief Function wich calculates the Morse potential with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @return dnMorse           TYPE(dnSca):           derived type with a value (pot),,if required, its derivatives (gradient (grad) and hessian (hess)).
!! @param dnR                TYPE(dnSca):           derived type with the value of "r" and,if required, its derivatives.
!! @param Para_Morse         TYPE(Param_Morse):   derived type with the Morse parameters.
  FUNCTION dnMorse(dnR,Para_Morse)
    USE mod_dnMatPot
    USE mod_dnSca

    TYPE(dnSca)                          :: dnMorse
    TYPE(dnSca),          intent(in)     :: dnR

    TYPE (Param_Morse), intent(in)    :: Para_Morse

    !local variable
    TYPE(dnSca)     :: dnbeta

    !write(out_unitp,*) 'BEGINNING in dnMorse'
    !write(out_unitp,*) 'dnR'
    !CALL write_dnSca(dnR)

    dnbeta  = exp(-Para_Morse%a*(dnR-Para_Morse%req))
    !write(out_unitp,*) 'dnbeta'
    !CALL write_dnSca(dnbeta)

    dnMorse = Para_Morse%D * (ONE-dnbeta)**2

     CALL dealloc_dnSca(dnbeta)

    !write(out_unitp,*) 'Morse at',get_d0_FROM_dnSca(dnR)
    !CALL Write_dnSca(dnMorse)
    !write(out_unitp,*) 'END in dnMorse'
    !flush(out_unitp)

  END FUNCTION dnMorse

END MODULE mod_MorsePot
