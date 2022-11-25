!===========================================================================
!===========================================================================
!This file is part of QuantumModelLib (QML).
!===============================================================================
! MIT License
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
!    Copyright (c) 2022 David Lauvergnat [1]
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

!> @brief Module which makes the initialization, calculation of the Retinal_JPCB2000 potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_Retinal_JPCB2000_m
  USE QMLLib_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the Retinal_JPCB2000 parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_Retinal_JPCB2000_t

   PRIVATE

   !The parameters of the model are in eV:
   !  m=4.84E-4,E1=2.48,W0=3.6,W1=1.09,w=0.19,kappa=0.1,and lambda=0.19.

   real(kind=Rkind) :: m      = ONE/4.84e-4_Rkind
   real(kind=Rkind) :: E1     = 2.48_Rkind
   real(kind=Rkind) :: W0     = 3.6_Rkind
   real(kind=Rkind) :: W1     = 1.09_Rkind
   real(kind=Rkind) :: w      = 0.19_Rkind
   real(kind=Rkind) :: kappa  = 0.1_Rkind
   real(kind=Rkind) :: lambda = 0.19_Rkind
   real(kind=Rkind) :: Cbs    = 0.0_Rkind

   real(kind=Rkind) :: wi_pub(25)       =                                       &
         [   0.0_Rkind,   0.0_Rkind,  792.8_Rkind,  842.8_Rkind,  866.2_Rkind,  &
           882.4_Rkind, 970.3_Rkind,   976.0_Rkind,  997.0_Rkind, 1017.1_Rkind, &
          1089.6_Rkind, 1189.0_Rkind, 1214.7_Rkind, 1238.1_Rkind, 1267.9_Rkind, &
          1317.0_Rkind, 1359.0_Rkind, 1389.0_Rkind, 1428.4_Rkind, 1434.9_Rkind, &
          1451.8_Rkind, 1572.8_Rkind, 1612.1_Rkind, 1629.2_Rkind, 1659.1_Rkind]
   real(kind=Rkind) :: ciwi_inv_pub(25) =                                       &
     [0.0_Rkind,    0.0_Rkind,   0.175_Rkind,  0.2_Rkind,    0.175_Rkind, &
     0.225_Rkind,  0.55_Rkind,  0.3_Rkind,    0.33_Rkind,   0.45_Rkind,  &
     0.125_Rkind,  0.175_Rkind, 0.44_Rkind,   0.5_Rkind,    0.475_Rkind, &
     0.238_Rkind,  0.25_Rkind,  0.25_Rkind,   0.25_Rkind,   0.225_Rkind, &
     0.225_Rkind,  0.25_Rkind,  0.225_Rkind,  0.125_Rkind,  0.225_Rkind]


   real(kind=Rkind), allocatable :: wi(:)
   real(kind=Rkind), allocatable :: ci(:)
   integer,          allocatable :: BathMode_Order(:)

   !integer :: BathMode_Order(25) = [ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,     &
   !                                 14,15,16,17,18,19,20,21,22,23,24,25]

   CONTAINS
    PROCEDURE :: EvalPot_QModel => EvalPot_QML_Retinal_JPCB2000
    PROCEDURE :: Write_QModel    => Write_QML_Retinal_JPCB2000
    PROCEDURE :: Write0_QModel   => Write0_QML_Retinal_JPCB2000
  END TYPE QML_Retinal_JPCB2000_t

  PUBLIC :: QML_Retinal_JPCB2000_t,Init_QML_Retinal_JPCB2000

  CONTAINS
!> @brief Function which makes the initialization of the Retinal_JPCB2000 parameters.
!!
!! @param QModel             TYPE(QML_Retinal_JPCB2000_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_Retinal_JPCB2000(QModel_in,read_param,nio_param_file) RESULT(QModel)
  IMPLICIT NONE

    TYPE (QML_Retinal_JPCB2000_t)                :: QModel ! RESULT

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param

    real (kind=Rkind), parameter :: auTOeV      = 27.211384_Rkind
    real (kind=Rkind), parameter :: auTOcm_inv  = 219474.631443_Rkind
    integer :: i
    character(len=Name_longlen)  :: Temp_string


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_Retinal_JPCB2000'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL Init0_QML_Empty(QModel%QML_Empty_t,QModel_in)

    QModel%nsurf    = 2

    IF (QModel%ndim == 0) THEN
      ! it means that ndim was not present in the CALL Init_Model().
      ! => ndim is set to the default value (2)
      QModel%ndim  = 2
      QModel%pot_name = 'retinal_jpcb2000'
    END IF
    IF (QModel%ndim  > 2) QModel%pot_name = 'retinal_cp2000'

    !The value of QModel%ndim must be in [2:25]
    IF (QModel%ndim < 2 .OR. QModel%ndim > 25) THEN
       write(out_unitp,*) 'Write_QModel'
       CALL QModel%Write_QModel(out_unitp)
       write(out_unitp,*) ' ERROR in ',name_sub
       write(out_unitp,*) ' ndim range MUST be [2:25]. ndim: ',QModel%ndim
       STOP 'ERROR in Init_QML_Retinal_JPCB2000: ndim range MUST be [2:25]'
    END IF
    IF (debug) write(out_unitp,*) 'ndim',QModel%ndim

    IF (allocated(QModel%BathMode_Order)) deallocate(QModel%BathMode_Order)
    IF (allocated(QModel%wi))             deallocate(QModel%wi)
    IF (allocated(QModel%ci))             deallocate(QModel%ci)

    allocate(QModel%BathMode_Order(QModel%ndim))
    allocate(QModel%wi(QModel%ndim))
    allocate(QModel%ci(QModel%ndim))

    IF (read_param) THEN
      QModel%BathMode_Order(:)   = 0
      QModel%BathMode_Order(1:2) = [1,2]
      read(nio_param_file,*) Temp_string,QModel%BathMode_Order(3:QModel%ndim)
      QModel%BathMode_Order(3:QModel%ndim) = QModel%BathMode_Order(3:QModel%ndim) + 2
    ELSE
      QModel%BathMode_Order(:) = [(i,i=1,QModel%ndim)]
    END IF
    write(out_unitp,*) 'BathMode_Order (internal)',QModel%BathMode_Order(3:QModel%ndim)
    write(out_unitp,*) 'BathMode_Order (read    )',QModel%BathMode_Order(3:QModel%ndim)-2


    IF (.NOT. QModel%PubliUnit) THEN
      QModel%m      = QModel%m      * auTOeV   ! 1/m has to be divided by auTOeV
      QModel%E1     = QModel%E1     / auTOeV
      QModel%W0     = QModel%W0     / auTOeV
      QModel%W1     = QModel%W1     / auTOeV
      QModel%w      = QModel%w      / auTOeV
      QModel%kappa  = QModel%kappa  / auTOeV
      QModel%lambda = QModel%lambda / auTOeV

      IF (QModel%ndim > 2) THEN
        QModel%wi(:) = QModel%wi_pub(QModel%BathMode_Order) / auTOcm_inv
      END IF

    ELSE ! if PubliUnit the wi must be convert in eV
      IF (QModel%ndim > 2) THEN
        QModel%wi(:) = QModel%wi_pub(QModel%BathMode_Order) / auTOcm_inv*auTOeV
      END IF
    END IF

    IF (QModel%ndim > 2) THEN
      QModel%ci(:) = QModel%ciwi_inv_pub(QModel%BathMode_Order) * QModel%wi
    END IF

    IF (debug) write(out_unitp,*) 'init Q0 of Retinal_JPCB2000'
    allocate(QModel%Q0(QModel%ndim))
    QModel%Q0 = ZERO

    IF (debug) write(out_unitp,*) 'init d0GGdef of Retinal_JPCB2000'
    CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)

    QModel%d0GGdef(1,1) = ONE/QModel%m
    QModel%d0GGdef(2,2) = QModel%w

    DO i=3,QModel%ndim
      QModel%d0GGdef(i,i) = QModel%wi(i)
    END DO

    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_QML_Retinal_JPCB2000
!> @brief Subroutine wich prints the current QML_Retinal_JPCB2000 parameters.
!!
!! @param QModel            CLASS(QML_Retinal_JPCB2000_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_Retinal_JPCB2000(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_Retinal_JPCB2000_t),   intent(in) :: QModel
    integer,              intent(in) :: nio

    write(nio,*)
    write(nio,*) 'Retinal_JPCB2000 current parameters'
    write(nio,*) '   Parameters of the model from:'
    write(nio,*) '  S. Hahn, G. Stock / Chemical Physics 259 (2000) 297-312'
    write(nio,*) '     doi: 10.1016/S0301-0104(00)00201-9'
    write(nio,*)
    IF (QModel%PubliUnit) THEN
      write(nio,*) ' Unit in eV (as in the publication)'
    ELSE
      write(nio,*) ' Unit in au (atomic units)'
    END IF
    write(nio,*)
    write(nio,*) '   1/m    = ',ONE/QModel%m
    write(nio,*) '   m      = ',QModel%m
    write(nio,*) '   E1     = ',QModel%E1
    write(nio,*) '   W0     = ',QModel%W0
    write(nio,*) '   W1     = ',QModel%W1
    write(nio,*) '   w      = ',QModel%w
    write(nio,*) '   kappa  = ',QModel%kappa
    write(nio,*) '   lambda = ',QModel%lambda
    write(nio,*)
    IF (QModel%ndim > 2) THEN
      write(nio,*) '    wi    =',QModel%wi(3:QModel%ndim)
      write(nio,*) '    ci/wi =',QModel%ci(3:QModel%ndim)/QModel%wi(3:QModel%ndim)
      write(nio,*) '    ci    =',QModel%ci(3:QModel%ndim)
    END IF
    write(nio,*) 'end Retinal_JPCB2000 default parameters'
    write(nio,*)

  END SUBROUTINE Write_QML_Retinal_JPCB2000
!> @brief Subroutine wich prints the default QML_Retinal_JPCB2000 parameters.
!!
!! @param QModel            CLASS(QML_Retinal_JPCB2000_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write0_QML_Retinal_JPCB2000(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_Retinal_JPCB2000_t),   intent(in) :: QModel
    integer,              intent(in) :: nio

    write(nio,*)
    write(nio,*) 'Retinal_JPCB2000 default parameters'
    write(nio,*) '   Parameters (in eV) of the model from:'
    write(nio,*) '  S. Hahn, G. Stock / Chemical Physics 259 (2000) 297-312'
    write(nio,*) '     doi: 10.1016/S0301-0104(00)00201-9'
    write(nio,*)
    write(nio,*) '   1/m    = 4.84e-4'
    write(nio,*) '   E1     = 2.48'
    write(nio,*) '   W0     = 3.6'
    write(nio,*) '   W1     = 1.09'
    write(nio,*) '   w      = 0.19'
    write(nio,*) '   kappa  = 0.1'
    write(nio,*) '   lambda = 0.19'
    write(nio,*)
    write(nio,*) 'end Retinal_JPCB2000 default parameters'
    write(nio,*)

  END SUBROUTINE Write0_QML_Retinal_JPCB2000

!> @brief Subroutine wich calculates the Retinal_JPCB2000 potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(QML_Retinal_JPCB2000_t):   derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_Retinal_JPCB2000(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE ADdnSVM_m
  IMPLICIT NONE

    CLASS(QML_Retinal_JPCB2000_t),   intent(in)    :: QModel
    TYPE (dnS_t),                    intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),                    intent(in)    :: dnQ(:)
    integer,                         intent(in)    :: nderiv

    TYPE (dnS_t) :: VBath
    integer      :: i

    IF (QModel%ndim > 2) THEN
      VBath = HALF*sum( QModel%wi(3:QModel%ndim)*dnQ(3:QModel%ndim)**2 ) + &
              QModel%Cbs * dnQ(2) * sum(dnQ(3:QModel%ndim) )


      Mat_OF_PotDia(1,1) = HALF*QModel%W0*(ONE-cos(dnQ(1))) +                   &
                           HALF*QModel%W*dnQ(2)**2 + VBath

      Mat_OF_PotDia(2,2) = QModel%E1 - HALF*QModel%W1*(ONE-cos(dnQ(1))) +       &
                           HALF*QModel%W*dnQ(2)**2 + QModel%kappa*dnQ(2) +      &
                        VBath + sum( QModel%ci(3:QModel%ndim)*dnQ(3:QModel%ndim) )

      Mat_OF_PotDia(1,2) = QModel%lambda*dnQ(2)
      Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)

      CALL dealloc_dnS(VBath)

    ELSE
      Mat_OF_PotDia(1,1) = HALF*QModel%W0*(ONE-cos(dnQ(1))) +                   &
                           HALF*QModel%W*dnQ(2)**2

      Mat_OF_PotDia(2,2) = QModel%E1 - HALF*QModel%W1*(ONE-cos(dnQ(1))) +       &
                           HALF*QModel%W*dnQ(2)**2 + QModel%kappa*dnQ(2)

      Mat_OF_PotDia(1,2) = QModel%lambda*dnQ(2)
      Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)
    END IF


  END SUBROUTINE EvalPot_QML_Retinal_JPCB2000

END MODULE QML_Retinal_JPCB2000_m
