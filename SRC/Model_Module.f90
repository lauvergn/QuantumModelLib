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
!    Copyright 2016 David Lauvergnat
!      with contributions of Félix MOUHAT and Liang LIANG
!
!===========================================================================
!===========================================================================
MODULE mod_Model
!$ USE omp_lib
  USE mod_NumParameters
  USE mod_dnMatPot,       ONLY: dnMatPot,alloc_dnMatPot,dealloc_dnMatPot,Check_NotAlloc_dnMatPot, &
                                get_maxval_OF_dnMatPot,Write_dnMatPot,get_nsurf_FROM_dnMatPot,    &
                                get_ndim_FROM_dnMatPot,dnmatpot2_minus_dnmatpot1,assignment (=)
  USE mod_MorsePot,       ONLY: Param_Morse,Init_MorsePot,Write_MorsePot,Eval_MorsePot
  USE mod_HenonHeilesPot, ONLY: Param_HenonHeiles,Init_HenonHeilesPot,Write_HenonHeilesPot,Eval_HenonHeilesPot
  USE mod_TullyPot,       ONLY: Param_Tully,Init_TullyPot,Write_TullyPot,Eval_TullyPot
  USE mod_LinearHBondPot, ONLY: Param_LinearHBond,Init_LinearHBondPot,Write_LinearHBondPot,Eval_LinearHBondPot
  USE mod_BuckPot,        ONLY: Param_Buck,Init_BuckPot,Write_BuckPot,Eval_BuckPot
  USE mod_PhenolPot,      ONLY: Param_Phenol,Init_PhenolPot,Write_PhenolPot,Eval_PhenolPot
  USE mod_SigmoidPot,     ONLY: Param_Sigmoid,Init_SigmoidPot,Write_SigmoidPot,Eval_SigmoidPot

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Param_Model,Init_Model,Eval_Pot,Write_Model
  PUBLIC :: calc_pot,calc_grad,calc_hess,calc_pot_grad,calc_pot_grad_hess
  PUBLIC :: Check_analytical_numerical_derivatives
  PUBLIC :: Eval_pot_ON_Grid

  real (kind=Rkind)                     :: step = ONETENTH**4

  TYPE Param_Model
    integer :: nsurf       = 1
    integer :: ndim        = 1
    logical :: numeric     = .FALSE.
    logical :: adiabatic   = .TRUE.
    logical :: PubliUnit   = .FALSE. ! when PubliUnit=.TRUE., the units of a reference (publi ...) are used. Default (atomic unit)

    character (len=:), allocatable :: pot_name

    TYPE (dnMatPot)          :: Vec0 ! to get the correct phase of the adiatic couplings

    ! list of potentials ....
    TYPE (Param_Morse)       :: Para_Morse
    TYPE (Param_Buck)        :: Para_Buck
    TYPE (Param_Sigmoid)     :: Para_Sigmoid

    TYPE (Param_LinearHBond) :: Para_LinearHBond
    TYPE (Param_HenonHeiles) :: Para_HenonHeiles
    TYPE (Param_Tully)       :: Para_Tully
    TYPE (Param_Phenol)      :: Para_Phenol

  END TYPE Param_Model


CONTAINS

  SUBROUTINE Read_Model(Para_Model,nio,option1)
  USE mod_Lib
  IMPLICIT NONE

    TYPE (Param_Model), intent(inout) :: Para_Model
    integer, intent(in)               :: nio
    integer, intent(inout)            :: option1

    integer :: ndim,nsurf,nderiv,option
    logical :: adiabatic,numeric,PubliUnit
    character (len=20) :: pot_name

    ! Namelists for input file
    namelist /potential/ ndim,nsurf,pot_name,numeric,adiabatic,option,PubliUnit


    ! Default values defined
    ndim      = 1
    nsurf     = 1
    adiabatic = .true.
    pot_name  = 'morse'
    numeric   = .false.
    PubliUnit = .FALSE.
    option    = -1 ! no option


    write(out_unitp,*) 'Reading input file . . .'
    read(nio,nml=potential)
    !write(out_unitp,nml=potential)

    option1              = option
    Para_Model%ndim      = ndim
    Para_Model%nsurf     = nsurf
    Para_Model%adiabatic = adiabatic
    Para_Model%numeric   = numeric
    Para_Model%pot_name  = strdup(pot_name)
    Para_Model%PubliUnit = PubliUnit

  END SUBROUTINE Read_Model
 
  SUBROUTINE Init_Model(Para_Model,pot_name,ndim,nsurf,                 &
                        read_param,param_file_name,nio_param_file,      &
                        option,PubliUnit)
  USE mod_Lib
  IMPLICIT NONE

    TYPE (Param_Model), intent(inout)      :: Para_Model
    character (len=*), intent(in),optional :: pot_name
    integer, intent(in), optional          :: ndim,nsurf

    logical, intent(in), optional          :: read_param
    integer, intent(in), optional          :: nio_param_file
    character (len=*), intent(in),optional :: param_file_name

    integer, intent(in), optional          :: option
    logical, intent(in), optional          :: PubliUnit

    integer :: option_loc,nio_loc
    logical :: read_param_loc
    character (len=:), allocatable :: param_file_name_loc

    IF (present(option)) THEN
      option_loc = option
    ELSE
      option_loc = -1
    END IF

    IF (present(read_param)) THEN
      read_param_loc = read_param
    ELSE
      read_param_loc = .FALSE.
    END IF
    IF (present(nio_param_file)) THEN
      nio_loc = nio_param_file
    ELSE
      IF (present(param_file_name)) THEN
        IF (len_trim(param_file_name) == 0) THEN
          param_file_name_loc = strdup("input.dat")
          nio_loc = 99
        ELSE
          param_file_name_loc = strdup(param_file_name)
          nio_loc = 99
        END IF
      ELSE
      nio_loc = in_unitp
      END IF
    END IF

    IF (present(pot_name)) THEN
      Para_Model%pot_name  = strdup(pot_name)
    ELSE
      IF (.NOT. read_param_loc) THEN
        write(out_unitp,*) 'ERROR in Init_Model'
        write(out_unitp,*) ' pot_name is not present and read_param=F'
        STOP
      END IF
    END IF

    IF (present(ndim)) THEN
      Para_Model%ndim      = ndim
    ELSE
      Para_Model%ndim      = 1
    END IF
    IF (present(nsurf)) THEN
      Para_Model%nsurf     = nsurf
    ELSE
      Para_Model%nsurf     = 1
    END IF

    IF (present(PubliUnit)) THEN
      Para_Model%PubliUnit      = PubliUnit
    ELSE
      Para_Model%PubliUnit      = .FALSE.
    END IF


    CALL dealloc_dnMatPot(Para_Model%Vec0)

    IF (read_param_loc) THEN

      IF (nio_loc /= in_unitp) THEN
        open(unit=nio_loc,file=param_file_name_loc,status='old',form='formatted')
      END IF
      CALL Read_Model(Para_Model,nio_loc,option_loc)

    END IF

    IF (Para_Model%adiabatic) THEN
      write(out_unitp,*) 'Adiabatic potential . . .'
    ELSE
      write(out_unitp,*) 'Non-adiabatic potential . . .'
    END IF

    IF (Para_Model%numeric) write(out_unitp,*) 'You have decided to perform a numeric checking of the analytic formulas.'


    CALL string_uppercase_TO_lowercase(Para_Model%pot_name)
    write(out_unitp,*) 'Para_Model%pot_name: ',Para_Model%pot_name

    SELECT CASE (Para_Model%pot_name)
    CASE ('morse')
      !! Morse potential: V(R) = D*(1-exp(-a*(r-Req))**2

      Para_Model%ndim      = 1
      Para_Model%nsurf     = 1

      CALL Init_MorsePot(Para_Model%Para_Morse,nio=nio_loc,read_param=read_param_loc)

    CASE ('sigmoid')
      !! sigmoid function: A * 1/2(1+e*tanh((x-B)/C))  remark: e=+/-1
      Para_Model%nsurf     = 1
      Para_Model%ndim      = 1

      CALL Init_SigmoidPot(Para_Model%Para_Sigmoid,nio=nio_loc,read_param=read_param_loc)

    CASE ('buck')
      !! Buckingham potential: V(R) = A*exp(-B*R)-C/R^6

      Para_Model%ndim      = 1
      Para_Model%nsurf     = 1

      CALL Init_BuckPot(Para_Model%Para_Buck,nio=nio_loc,read_param=read_param_loc)

    CASE ('hbond')
      !write(out_unitp,*) 'We are working with Linear H-Bond potential'

      Para_Model%ndim      = 2
      Para_Model%nsurf     = 1

      CALL Init_LinearHBondPot(Para_Model%Para_LinearHBond,               &
                               nio=nio_loc,read_param=read_param_loc,   &
                               PubliUnit=Para_Model%PubliUnit)

    CASE ('henonheiles')

        Para_Model%nsurf     = 1

        CALL Init_HenonHeilesPot(Para_Model%Para_HenonHeiles,ndim=Para_Model%ndim, &
                                 nio=nio_loc,read_param=read_param_loc)

    CASE ('tully')
        !! from Tully, J. Chem. Phys. V93, pp15, 1990
        Para_Model%nsurf     = 2
        Para_Model%ndim      = 1

        write(out_unitp,*) 'option_loc',option_loc
        CALL Init_TullyPot(Para_Model%Para_Tully, option=option_loc,      &
                           nio=nio_loc,read_param=read_param_loc)

    CASE ('phenol')
        !! from Z. Lan, W. Domcke, V. Vallet, A.L. Sobolewski, S. Mahapatra, ...
        !!  J. Chem. Phys. 122 (2005) 224315. doi:10.1063/1.1906218.
        Para_Model%nsurf     = 3
        Para_Model%ndim      = 2

        CALL Init_PhenolPot(Para_Model%Para_Phenol,PubliUnit=Para_Model%PubliUnit)

    CASE DEFAULT
        STOP 'STOP in Init_Model: Other potentials have to be done'
    END SELECT

    IF (read_param_loc .AND. nio_loc /= in_unitp) THEN
       close(nio_loc)
    END IF

    !CALL Write_Model(Para_Model)

  END SUBROUTINE Init_Model

RECURSIVE SUBROUTINE Eval_Pot(Para_Model,Q,PotVal,nderiv,Vec)
  USE mod_diago
  IMPLICIT NONE

    TYPE (Param_Model), intent(inout)         :: Para_Model
    TYPE (dnMatPot), intent(inout)            :: PotVal
    real (kind=Rkind),dimension(:),intent(in) :: Q
    integer, intent(in), optional             :: nderiv
    TYPE (dnMatPot), intent(inout), optional  :: Vec

    ! local variable
    integer                    :: i,j,id,nderiv_loc
    TYPE (dnMatPot)            :: PotVal_dia,Vec_loc

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Eval_Pot'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) ' BEGINNING ',name_sub
      IF (present(nderiv)) write(out_unitp,*) '   nderiv',nderiv
      write(out_unitp,*) '  numeric   ',Para_Model%numeric
      write(out_unitp,*) '  adiabatic ',Para_Model%adiabatic
    END IF

    IF (present(nderiv)) THEN
      nderiv_loc = max(0,nderiv)
      nderiv_loc = min(2,nderiv_loc)
    ELSE
      nderiv_loc = 0
    END IF


    IF ( Check_NotAlloc_dnMatPot(PotVal,nderiv_loc) ) THEN
      CALL alloc_dnMatPot(PotVal,nsurf=Para_Model%nsurf,ndim=Para_Model%ndim,&
                          nderiv=nderiv_loc)
      PotVal = ZERO
      !CALL Write_dnMatPot(PotVal)
      !STOP 'PotVal%dn NOT allocated in Eval_Pot'
    END IF

    IF (Para_Model%numeric .AND. nderiv_loc > 0) THEN
      CALL Eval_Pot_Numeric(Para_Model,Q,PotVal,nderiv_loc)
    ELSE

      SELECT CASE (Para_Model%pot_name)
      CASE ('morse')

        CALL Eval_MorsePot(PotVal,Q(1),Para_Model%Para_Morse,nderiv_loc)

      CASE ('buck')

        CALL Eval_BuckPot(PotVal,Q(1),Para_Model%Para_Buck,nderiv_loc)

      CASE ('sigmoid')

        CALL Eval_SigmoidPot(PotVal,Q(1),Para_Model%Para_Sigmoid,nderiv_loc)

      CASE ('hbond')

        !CALL Eval_LinearHBondPot_old(PotVal,Q,Para_Model%Para_LinearHBond,nderiv_loc)
        CALL Eval_LinearHBondPot(PotVal,Q,Para_Model%Para_LinearHBond,nderiv_loc)

      CASE ('henonheiles')

        CALL Eval_HenonHeilesPot(PotVal,Q,Para_Model%Para_HenonHeiles,nderiv_loc)

      CASE ('tully')

        CALL Eval_TullyPot(PotVal,Q(1),Para_Model%Para_Tully,nderiv_loc)

      CASE ('phenol')

        !CALL Eval_PhenolPot_old(PotVal,Q,Para_Model%Para_Phenol,nderiv_loc)
        CALL Eval_PhenolPot(PotVal,Q,Para_Model%Para_Phenol,nderiv_loc)

      CASE DEFAULT
        STOP 'ERROR in Eval_Pot: Other potentials have to be done'
      END SELECT
    END IF

    IF ( Para_Model%adiabatic .AND. Para_Model%nsurf > 1 .AND. &
        .NOT. (Para_Model%numeric .AND. nderiv_loc > 0) ) THEN
      IF (debug) THEN
        write(out_unitp,*) 'PotVal (dia)'
        CALL Write_dnMatPot(PotVal,6)
        flush(out_unitp)
      END IF

      PotVal_dia = PotVal

      IF (present(Vec)) THEN
        CALL dia_TO_adia(PotVal_dia,PotVal,Vec,Para_Model%Vec0,nderiv)
      ELSE
        CALL dia_TO_adia(PotVal_dia,PotVal,Vec_loc,Para_Model%Vec0,nderiv)
        CALL dealloc_dnMatPot(Vec_loc)
      END IF

      CALL dealloc_dnMatPot(PotVal_dia)


    END IF

    IF (debug) THEN
      IF ( Para_Model%adiabatic) write(out_unitp,*) 'PotVal (adia)'
      IF ( .NOT. Para_Model%adiabatic) write(out_unitp,*) 'PotVal (dia)'
      CALL Write_dnMatPot(PotVal,6)
      write(out_unitp,*) ' END ',name_sub
      flush(out_unitp)
    END IF


  END SUBROUTINE Eval_Pot

  SUBROUTINE Eval_Pot_Numeric(Para_Model,Q,PotVal,nderiv)
  IMPLICIT NONE

    TYPE (Param_Model), intent(inout)         :: Para_Model
    TYPE (dnMatPot), intent(inout)            :: PotVal
    real (kind=Rkind) ,intent(in)             :: Q(:)
    integer, intent(in)                       :: nderiv

    ! local variable
    real (kind=Rkind),dimension(size(Q))  :: Q_loc
    TYPE (dnMatPot)                   :: PotVal_loc0
    integer                           :: i,j


    Q_loc = Q
    CALL alloc_dnMatPot(PotVal_loc0,Para_Model%nsurf,Para_Model%ndim,nderiv=0)

    ! no derivative : PotVal%d0
    PotVal = ZERO
    CALL Eval_Pot(Para_Model,Q,PotVal_loc0,nderiv=0)
    PotVal%d0 = PotVal_loc0%d0

    IF (nderiv >= 1) THEN ! 1st derivatives

        ! Numeric evaluation of forces
      DO i=1,Para_Model%ndim

        Q_loc(i) = Q(i) + step        ! q+dq
        CALL Eval_Pot(Para_Model,Q_loc,PotVal_loc0,nderiv=0) ! Ep(q+dq)
        PotVal%d1(:,:,i) = PotVal_loc0%d0

        IF (nderiv >= 2) THEN
          PotVal%d2(:,:,i,i) = PotVal_loc0%d0
        END IF

        Q_loc(i) = Q(i) - step        ! q-dq
        CALL Eval_Pot(Para_Model,Q_loc,PotVal_loc0,nderiv=0) ! Ep(q-dq)
        PotVal%d1(:,:,i) = (PotVal%d1(:,:,i)-PotVal_loc0%d0)/(TWO*step)

        IF (nderiv >= 2) THEN
          PotVal%d2(:,:,i,i) = (PotVal%d2(:,:,i,i) + PotVal_loc0%d0 - TWO*PotVal%d0)/ &
                                step**2
        END IF

        Q_loc(i) = Q(i)

      END DO
    END IF

    IF (nderiv >= 2) THEN ! 2d derivatives

      DO i=1,Para_Model%ndim
      DO j=i+1,Para_Model%ndim

        Q_loc(i) = Q(i) + step        ! qi+dq
        Q_loc(j) = Q(j) + step        ! qj+dq
        CALL Eval_Pot(Para_Model,Q_loc,PotVal_loc0,nderiv=0)
        PotVal%d2(:,:,j,i) = PotVal_loc0%d0

        Q_loc(i) = Q(i) - step        ! qi-dq
        Q_loc(j) = Q(j) - step        ! qj-dq
        CALL Eval_Pot(Para_Model,Q_loc,PotVal_loc0,nderiv=0)
        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i) + PotVal_loc0%d0

        Q_loc(i) = Q(i) + step        ! qi+dq
        Q_loc(j) = Q(j) - step        ! qj-dq
        CALL Eval_Pot(Para_Model,Q_loc,PotVal_loc0,nderiv=0)
        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i) - PotVal_loc0%d0

        Q_loc(i) = Q(i) - step        ! qi-dq
        Q_loc(j) = Q(j) + step        ! qj+dq
        CALL Eval_Pot(Para_Model,Q_loc,PotVal_loc0,nderiv=0)
        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i) - PotVal_loc0%d0

        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i)/(FOUR*step**2)
        PotVal%d2(:,:,i,j) = PotVal%d2(:,:,j,i)

        Q_loc(i) = Q(i)
        Q_loc(j) = Q(j)
      END DO
      END DO
    END IF

    CALL dealloc_dnMatPot(PotVal_loc0)

  END SUBROUTINE Eval_Pot_Numeric

  SUBROUTINE dia_TO_adia(PotVal_dia,PotVal_adia,Vec,Vec0,nderiv)
    USE mod_diago
    IMPLICIT NONE

    TYPE (dnMatPot), intent(in)               :: PotVal_dia
    TYPE (dnMatPot), intent(inout)            :: PotVal_adia,Vec,Vec0

    integer, intent(in), optional             :: nderiv

    ! local variable
    integer                    :: i,j,id,jd,nderiv_loc,ndim,nsurf
    real (kind=Rkind), allocatable :: Eig(:),tVec(:,:),Vdum(:),Vi(:)
    TYPE (dnMatPot)                :: PotVal_dia_onadia

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='dia_TO_adia'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) ' BEGINNING ',name_sub
      IF (present(nderiv)) write(out_unitp,*) '   nderiv',nderiv
    END IF

    IF (present(nderiv)) THEN
      nderiv_loc = max(0,nderiv)
      nderiv_loc = min(2,nderiv_loc)
    ELSE
      nderiv_loc = 0
    END IF


    IF ( Check_NotAlloc_dnMatPot(PotVal_dia,nderiv_loc) ) THEN
      write(out_unitp,*) ' The diabatic potential MUST be allocated!'
      CALL Write_dnMatPot(PotVal_dia)
      STOP 'PotVal_dia%dn NOT allocated in "dia_TO_adia"'
    END IF
    IF (debug) THEN
      write(out_unitp,*) 'PotVal_dia'
      CALL Write_dnMatPot(PotVal_dia,6)
      flush(out_unitp)
    END IF

    nsurf = get_nsurf_FROM_dnMatPot(PotVal_dia)
    ndim  = get_ndim_FROM_dnMatPot(PotVal_dia)

    IF ( Check_NotAlloc_dnMatPot(PotVal_adia,nderiv_loc) ) THEN
      CALL alloc_dnMatPot(PotVal_adia,nsurf=nsurf,ndim=ndim,nderiv=nderiv_loc)
    END IF
    PotVal_adia = ZERO


    IF ( Check_NotAlloc_dnMatPot(Vec,nderiv_loc) ) THEN
      CALL alloc_dnMatPot(Vec,nsurf=nsurf,ndim=ndim,nderiv=nderiv_loc)
    END IF
    Vec = ZERO

    ! local variables
    CALL alloc_dnMatPot(PotVal_dia_onadia,nsurf=nsurf,ndim=ndim,nderiv=nderiv_loc)
    PotVal_dia_onadia = ZERO


    allocate(Eig(nsurf))
    allocate(tVec(nsurf,nsurf))
    CALL diagonalization(PotVal_dia%d0,Eig,Vec%d0,nsurf)
    IF (Check_NotAlloc_dnMatPot(Vec0,nderiv=0)) THEN
       !$OMP CRITICAL (CRIT_dia_TO_adia)
       IF (debug) write(out_unitp,*) 'init Vec0'
       CALL alloc_dnMatPot(Vec0,nsurf=nsurf,ndim=ndim,nderiv=0)
       Vec0%d0 = Vec%d0
       !$OMP END CRITICAL (CRIT_dia_TO_adia)
    ELSE ! change the phase is required
       IF (debug) write(out_unitp,*) 'Change phase?'
       flush(out_unitp)

       DO i=1,nsurf
         IF (dot_product(Vec0%d0(:,i),Vec%d0(:,i)) < ZERO) Vec%d0(:,i) = -Vec%d0(:,i)
       END DO
    END IF
    tVec = transpose(Vec%d0)

    IF (debug) write(out_unitp,*) 'Eig',Eig

    ! transfortion of PotVal_dia on the adiabatic basis (Vec)
    !PotVal_dia_onadia%d0 = matmul(tVec,matmul(PotVal_dia%d0,Vec%d0))
    DO i=1,nsurf
      PotVal_dia_onadia%d0(i,i) = Eig(i)
    END DO
    IF (nderiv_loc > 0) THEN
      DO id=1,ndim
        PotVal_dia_onadia%d1(:,:,id) = matmul(tVec,matmul(PotVal_dia%d1(:,:,id),Vec%d0))
      END DO
    END IF
    IF (nderiv_loc > 1) THEN
      DO id=1,ndim
      DO jd=1,ndim
        PotVal_dia_onadia%d2(:,:,id,jd) = matmul(tVec,matmul(PotVal_dia%d2(:,:,id,jd),Vec%d0)  )
      END DO
      END DO
    END IF
    deallocate(Eig)
    deallocate(tVec)
    IF (debug) THEN
      write(out_unitp,*) "< Psi_j I dnHdia I Psi_i>"
      CALL Write_dnMatPot(PotVal_dia_onadia,6)
      flush(out_unitp)
    END IF



    ! no derivative
    PotVal_adia%d0 = ZERO
    Vec%d0         = ZERO
    DO i=1,nsurf
      PotVal_adia%d0(i,i) = PotVal_dia_onadia%d0(i,i)
      Vec%d0(i,i)         = ONE
    END DO

    ! 1st order derivatives
    IF (nderiv_loc > 0) THEN

      ! eigenvalue derivatives
      PotVal_adia%d1 = ZERO
      DO id=1,ndim
      DO i=1,nsurf
        PotVal_adia%d1(i,i,id) = PotVal_dia_onadia%d1(i,i,id)
      END DO
      END DO

      ! eigenvector derivatives projected on the eigenvectors
      DO id=1,ndim
      DO i=1,nsurf ! I Psi_i' >
      DO j=1,nsurf ! projection on I Psi_j >
        IF (j /= i) THEN
          Vec%d1(j,i,id) = - PotVal_dia_onadia%d1(j,i,id)/              &
                            ( PotVal_adia%d0(j,j) - PotVal_adia%d0(i,i) )
        ELSE
          Vec%d1(j,i,id) = ZERO
        END IF
      END DO
      END DO
      END DO

    END IF


    ! 2d order derivatives
    IF (nderiv_loc > 1) THEN
      PotVal_adia%d2 = ZERO
      allocate(Vdum(nsurf))
      ! eigenvalue derivatives
      DO id=1,ndim
      DO jd=1,ndim
        DO i=1,nsurf
          Vdum = matmul(PotVal_dia_onadia%d2(:,:,id,jd),Vec%d0(:,i))    + &
                 matmul(PotVal_dia_onadia%d1(:,:,jd),   Vec%d1(:,i,id)) + &
                 matmul(PotVal_dia_onadia%d1(:,:,id),   Vec%d1(:,i,jd))

          !write(out_unitp,*) 'Dum',id,jd,i,':',Vdum
          PotVal_adia%d2(i,i,id,jd) = dot_product(Vec%d0(:,i),Vdum)
        END DO

      END DO
      END DO
      deallocate(Vdum)
    END IF

    CALL alloc_dnMatPot(PotVal_dia_onadia)

    IF (debug) THEN
      write(out_unitp,*) 'PotVal_adia'
      CALL Write_dnMatPot(PotVal_adia,6)

      write(out_unitp,*) 'Vec'
      CALL Write_dnMatPot(Vec,6)
      write(out_unitp,*) ' END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE dia_TO_adia

  SUBROUTINE Write_Model(Para_Model,nio)
  IMPLICIT NONE

    TYPE (Param_Model), intent(in)              :: Para_Model
    integer,            intent(in), optional    :: nio

    integer :: nio_loc

    IF (present(nio)) THEN
      nio_loc = nio
    ELSE
      nio_loc = 6
    END IF

    IF (nio_loc /= 6) THEN
      open(nio_loc,file=trim(adjustl(Para_Model%pot_name))//'.out',form='formatted')
    END IF


    write(nio_loc,*) 'Output file for potential library'
    write(nio_loc,*)
    write(nio_loc,*) 'Potential parameters are written just below'
    write(nio_loc,*)
    write(nio_loc,*) 'nsurf:     ',Para_Model%nsurf
    write(nio_loc,*) 'ndim:      ',Para_Model%ndim
    write(nio_loc,*) 'numeric:   ',Para_Model%numeric
    write(nio_loc,*) 'adiabatic: ',Para_Model%adiabatic
    write(nio_loc,*)

    SELECT CASE (Para_Model%pot_name)
    CASE ('morse')
      CALL Write_MorsePot(Para_Model%Para_Morse,nio=nio_loc)
    CASE ('sigmoid')
      CALL Write_SigmoidPot(Para_Model%Para_Sigmoid,nio=nio_loc)
    CASE ('buck')
      CALL Write_BuckPot(Para_Model%Para_Buck,nio=nio_loc)
   CASE ('hbond')
      CALL Write_LinearHBondPot(Para_Model%Para_LinearHBond,nio=nio_loc)
    CASE ('henonheiles')
        CALL Write_HenonHeilesPot(Para_Model%Para_HenonHeiles,nio=nio_loc)
    CASE ('tully')
        CALL Write_TullyPot(Para_Model%Para_Tully,nio=nio_loc)
    CASE ('phenol')
        CALL Write_PhenolPot(Para_Model%Para_Phenol,nio=nio_loc)
    CASE DEFAULT
        write(nio_loc,*) 'WARNING in Write_Model: Other potentials have to be done'
    END SELECT
 
     IF (nio_loc /= 6) THEN
      close(nio_loc)
    END IF


  END SUBROUTINE Write_Model
  SUBROUTINE Check_analytical_numerical_derivatives(Para_Model,Q,nderiv)
  IMPLICIT NONE

    TYPE (Param_Model), intent(inout)            :: Para_Model
    real (kind=Rkind),dimension(:),intent(in)    :: Q
    integer, intent(in)                          :: nderiv

    TYPE (dnMatPot)           :: PotVal_ana,PotVal_num,PotVal_diff
    logical                   :: numeric_save
    real (kind=Rkind)         :: MaxPot,MaxDiffPot


      CALL alloc_dnMatPot(PotVal_ana,nsurf=Para_Model%nsurf,ndim=Para_Model%ndim,&
                          nderiv=nderiv)

      CALL alloc_dnMatPot(PotVal_num,nsurf=Para_Model%nsurf,ndim=Para_Model%ndim,&
                          nderiv=nderiv)

      numeric_save = Para_Model%numeric


      Para_Model%numeric = .FALSE.
      CALL Eval_Pot(Para_Model,Q,PotVal_ana,nderiv)
      Para_Model%numeric = .TRUE.
      CALL Eval_Pot(Para_Model,Q,PotVal_num,nderiv)
      MaxPot     = get_maxval_OF_dnMatPot(PotVal_ana)


      PotVal_diff = dnMatPot2_MINUS_dnMatPot1(PotVal_num,PotVal_ana)
      MaxDiffPot = get_maxval_OF_dnMatPot(PotVal_diff)

      write(out_unitp,'(a,e9.2)') 'max of the relative Potential diff:',MaxDiffPot/MaxPot
      write(out_unitp,'(a,l9)')   'Potential diff (numer-ana), ZERO?  ',(MaxDiffPot/MaxPot <= step)

      IF (MaxDiffPot/MaxPot > step) THEN
        CALL Write_dnMatPot(PotVal_diff,nio=6)
      END IF

      CALL dealloc_dnMatPot(PotVal_ana)
      CALL dealloc_dnMatPot(PotVal_num)
      CALL dealloc_dnMatPot(PotVal_diff)

      Para_Model%numeric = numeric_save

  END SUBROUTINE Check_analytical_numerical_derivatives

  SUBROUTINE Eval_pot_ON_Grid(Para_Model,Qmin,Qmax,nb_points,nderiv,grid_file)
  IMPLICIT NONE

    TYPE (Param_Model),           intent(inout)   :: Para_Model
    real (kind=Rkind),            intent(in)      :: Qmin(:),Qmax(:)
    integer, optional,            intent(in)      :: nb_points,nderiv
    character (len=*), optional,  intent(in)      :: grid_file

    integer           :: unit_grid_file

    integer           :: i,iq,jq,i1,i2,nb_points_loc,nderiv_loc,ndim_loc,i_Q(Para_Model%ndim)
    real (kind=Rkind) :: dQ(Para_Model%ndim),Q(Para_Model%ndim)

    TYPE (dnMatPot)           :: PotVal,Vec


    IF (size(Qmin) /= Para_Model%ndim .OR. size(Qmax) /= Para_Model%ndim) THEN
       write(out_unitp,*) ' ERROR in Eval_pot_ON_Grid'
       write(out_unitp,*) ' The size of Qmin or Qmax are different from Para_Model%ndim',size(Qmin),size(Qmax),Para_Model%ndim
       write(out_unitp,*) ' => Check the fortran'
       STOP 'ERROR in Eval_pot_ON_Grid: problem with Para_Model%ndim'
    END IF

    IF (present(grid_file)) THEN
      IF (len_trim(grid_file) == 0) THEN
        unit_grid_file = 6
      ELSE
        unit_grid_file = 99
        open(unit=unit_grid_file,file=trim(grid_file) )
      END IF
    ELSE
      unit_grid_file = 6
    END IF

    nb_points_loc = 100
    IF (present(nb_points)) nb_points_loc = nb_points
    nb_points_loc = max(nb_points_loc,2)

    IF (present(nderiv)) THEN
      nderiv_loc = nderiv
    ELSE
      nderiv_loc = 0
    END IF

    dQ       = (Qmax-Qmin) / real(nb_points_loc-1,kind=Rkind)
    ndim_loc = 0
    i_Q      = 0
    DO i=1,Para_Model%ndim
      IF (dQ(i) /= ZERO) THEN
        ndim_loc = ndim_loc + 1
        i_Q(ndim_loc) = i
      END IF
    END DO
    write(out_unitp,*) 'Para_Model%ndim',Para_Model%ndim
    write(out_unitp,*) 'ndim for the grid',ndim_loc
    write(out_unitp,*) 'i_Q',i_Q(1:ndim_loc)


    Q = Qmin

    IF (ndim_loc == 1) THEN
      i1 = i_Q(1)
      DO iq=0,nb_points_loc-1
        Q(i1) = Qmin(i1) + dQ(i1)*real(iq,kind=Rkind)
        IF (Para_Model%nsurf > 1 .AND. Para_Model%adiabatic) THEN
          CALL Eval_Pot(Para_Model,Q,PotVal,nderiv=max(1,nderiv_loc),Vec=Vec)

          IF (nderiv_loc == 0) THEN
            write(unit_grid_file,*) Q(i1),(PotVal%d0(i,i),i=1,Para_Model%nsurf),Vec%d1
          ELSE IF (nderiv_loc == 1) THEN
            write(unit_grid_file,*) Q(i1),(PotVal%d0(i,i),i=1,Para_Model%nsurf),(PotVal%d1(i,i,:),i=1,Para_Model%nsurf)
          ELSE
            write(unit_grid_file,*) Q(i1),(PotVal%d0(i,i),i=1,Para_Model%nsurf),(PotVal%d1(i,i,:),i=1,Para_Model%nsurf), &
                            (PotVal%d2(i,i,:,:),i=1,Para_Model%nsurf)
          END IF

        ELSE
          CALL Eval_Pot(Para_Model,Q,PotVal,nderiv=nderiv_loc)

          IF (nderiv_loc == 0) THEN
            write(unit_grid_file,*) Q(i1),PotVal%d0
          ELSE IF (nderiv_loc == 1) THEN
            write(unit_grid_file,*) Q(i1),PotVal%d0,PotVal%d1
          ELSE
            write(unit_grid_file,*) Q(i1),PotVal%d0,PotVal%d1,PotVal%d2
          END IF

        END IF
      END DO
    ELSE IF (ndim_loc == 2) THEN
      i1 = i_Q(1)
      i2 = i_Q(2)
      DO iq=0,nb_points_loc-1
      DO jq=0,nb_points_loc-1
        Q(i1) = Qmin(i1) + dQ(i1)*real(iq,kind=Rkind)
        Q(i2) = Qmin(i2) + dQ(i2)*real(jq,kind=Rkind)

        IF (Para_Model%nsurf > 1 .AND. Para_Model%adiabatic) THEN
          CALL Eval_Pot(Para_Model,Q,PotVal,nderiv=0,Vec=Vec)
          write(unit_grid_file,*) Q(i1),Q(i2),(PotVal%d0(i,i),i=1,Para_Model%nsurf)
        ELSE
          CALL Eval_Pot(Para_Model,Q,PotVal,nderiv=nderiv_loc)
          write(unit_grid_file,*) Q(i1),Q(i2),PotVal%d0
        END IF
      END DO
      write(unit_grid_file,*)
      END DO


    END IF

    CALL dealloc_dnMatPot(PotVal)
    CALL dealloc_dnMatPot(Vec)

    IF (unit_grid_file /= 6) THEN
      close(unit_grid_file)
    END IF


  END SUBROUTINE Eval_pot_ON_Grid


  SUBROUTINE calc_pot(V,Para_Model,Q)
  IMPLICIT NONE

    TYPE (Param_Model),       intent(inout) :: Para_Model
    real (kind=Rkind),      intent(in)      :: Q(:)
    real (kind=Rkind),      intent(inout)   :: V(:,:) ! it has to be allocated

    TYPE (dnMatPot)                         :: PotVal



    CALL alloc_dnMatPot(PotVal,nsurf=Para_Model%nsurf,ndim=Para_Model%ndim,nderiv=0)

    CALL Eval_Pot(Para_Model,Q,PotVal,nderiv=0)

    V = PotVal%d0

    CALL dealloc_dnMatPot(PotVal)


  END SUBROUTINE calc_pot
  SUBROUTINE calc_pot_grad(V,g,Para_Model,Q)
  IMPLICIT NONE

    TYPE (Param_Model),   intent(inout)     :: Para_Model
    real (kind=Rkind),      intent(in)      :: Q(:)
    real (kind=Rkind),      intent(inout)   :: V(:,:) ! it has to be allocated
    real (kind=Rkind),      intent(inout)   :: g(:,:,:) ! it has to be allocated

    TYPE (dnMatPot)           :: PotVal



    CALL alloc_dnMatPot(PotVal,nsurf=Para_Model%nsurf,ndim=Para_Model%ndim,nderiv=1)

    CALL Eval_Pot(Para_Model,Q,PotVal,nderiv=1)

    V = PotVal%d0
    g = PotVal%d1

    CALL dealloc_dnMatPot(PotVal)


  END SUBROUTINE calc_pot_grad
  SUBROUTINE calc_grad(g,Para_Model,Q)
  IMPLICIT NONE

    TYPE (Param_Model),   intent(inout)     :: Para_Model
    real (kind=Rkind),      intent(in)      :: Q(:)
    real (kind=Rkind),      intent(inout)   :: g(:,:,:) ! it has to be allocated

    TYPE (dnMatPot)           :: PotVal



    CALL alloc_dnMatPot(PotVal,nsurf=Para_Model%nsurf,ndim=Para_Model%ndim,nderiv=1)

    CALL Eval_Pot(Para_Model,Q,PotVal,nderiv=1)

    g = PotVal%d1

    CALL dealloc_dnMatPot(PotVal)


  END SUBROUTINE calc_grad
  SUBROUTINE calc_pot_grad_hess(V,g,h,Para_Model,Q)
  IMPLICIT NONE

    TYPE (Param_Model),   intent(inout)     :: Para_Model
    real (kind=Rkind),      intent(in)      :: Q(:)
    real (kind=Rkind),      intent(inout)   :: V(:,:) ! it has to be allocated
    real (kind=Rkind),      intent(inout)   :: g(:,:,:) ! it has to be allocated
    real (kind=Rkind),      intent(inout)   :: h(:,:,:,:) ! it has to be allocated


    TYPE (dnMatPot)           :: PotVal



    CALL alloc_dnMatPot(PotVal,nsurf=Para_Model%nsurf,ndim=Para_Model%ndim,nderiv=2)

    CALL Eval_Pot(Para_Model,Q,PotVal,nderiv=2)

    V = PotVal%d0
    g = PotVal%d1
    h = PotVal%d2

    CALL dealloc_dnMatPot(PotVal)


  END SUBROUTINE calc_pot_grad_hess
  SUBROUTINE calc_hess(h,Para_Model,Q)
  IMPLICIT NONE

    TYPE (Param_Model),   intent(inout)   :: Para_Model
    real (kind=Rkind),  intent(in)        :: Q(:)
    real (kind=Rkind),  intent(inout)     :: h(:,:,:,:) ! it has to be allocated


    TYPE (dnMatPot)           :: PotVal



    CALL alloc_dnMatPot(PotVal,nsurf=Para_Model%nsurf,ndim=Para_Model%ndim,nderiv=2)

    CALL Eval_Pot(Para_Model,Q,PotVal,nderiv=2)

    h = PotVal%d2

    CALL dealloc_dnMatPot(PotVal)


  END SUBROUTINE calc_hess
END MODULE mod_Model
