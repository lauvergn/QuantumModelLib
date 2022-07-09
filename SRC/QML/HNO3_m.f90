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

!> @brief Module which makes the initialization, calculation of the HNO3 potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_HNO3_m
  USE QMLLib_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE


  integer, parameter :: ndim    = 1
  integer, parameter :: max_fit = 8
  integer, parameter :: max_nn  = 10


!> @brief Derived type in which the HNO3 parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_HNO3_t

    PRIVATE

    real (kind=Rkind)  :: F(max_nn,0:max_fit,0:max_fit)  = ZERO
    integer            :: nn(0:ndim,0:max_fit,0:max_fit) = 0
    integer            :: nt(0:max_fit,0:max_fit)        = 0

   CONTAINS
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_HNO3
    PROCEDURE :: Eval_QModel_Func => EvalFunc_QML_HNO3
    PROCEDURE :: Write_QModel     => Write_QML_HNO3
    PROCEDURE :: Write0_QModel    => Write0_QML_HNO3
  END TYPE QML_HNO3_t

  PUBLIC :: QML_HNO3_t,Init_QML_HNO3

  CONTAINS
!> @brief Function which makes the initialization of the HNO3 parameters.
!!
!! @param QModel             TYPE(QML_HNO3_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_HNO3(QModel_in,read_param,nio_param_file) RESULT(QModel)
  IMPLICIT NONE

    TYPE (QML_HNO3_t)                           :: QModel ! RESULT

    TYPE(QML_Empty_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param


    integer :: ii,jj
    logical :: exist
    character (len=:), allocatable  :: FileName

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_HNO3'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL Init0_QML_Empty(QModel%QML_Empty_t,QModel_in)

    IF (QModel%ndim == 0) THEN
      ! it means that ndim was not present in the CALL Init_Model().
      ! => ndim is set to the default value (1)
      QModel%ndim  = 1
    END IF

    !The value of QModel%ndim must be 1 or 9
    IF (QModel%ndim /= 1 .AND. QModel%ndim /= 9) THEN
       write(out_unitp,*) 'Write_QModel'
       CALL QModel%Write_QModel(out_unitp)
       write(out_unitp,*) ' ERROR in ',name_sub
       write(out_unitp,*) ' ndim MUST equal to 1 or 9. ndim: ',QModel%ndim
       STOP 'ERROR in Init_QML_HNO3: ndim MUST equal to 1 or 9'
    END IF

    QModel%nsurf    = 1
    QModel%pot_name = 'hno3'

    QModel%ndimFunc = 1
    QModel%nb_Func  = 72

    ! read the parameters
    jj=0
    DO ii=0,max_fit

      FileName = make_FileName('InternalData/HNO3/inter_' //            &
                                int_TO_char(ii))
      CALL QML_read_para4d_HNO3(QModel%F(:,ii,jj),QModel%nn(:,ii,jj),ndim,       &
                       QModel%nt(ii,jj),max_nn,FileName,exist)
      IF ( .NOT. exist) STOP ' ERROR while reading HNO3 parameters'
    END DO


    !hess(i,j)
    DO jj=1,max_fit
    DO ii=1,max_fit
      FileName = make_FileName('InternalData/HNO3/inter_' //            &
                             int_TO_char(ii) // '_' // int_TO_char(jj) )
      CALL QML_read_para4d_HNO3(QModel%F(:,ii,jj),QModel%nn(:,ii,jj),ndim,       &
                       QModel%nt(ii,jj),max_nn,FileName,exist)
      IF ( .NOT. exist) STOP ' ERROR while reading HNO3 parameters'
    END DO
    END DO


    IF (debug) write(out_unitp,*) 'init Q0 of HNO3' ! for the rigid constraints
    QModel%Q0 = [2.65450_Rkind,2.28_Rkind,2.0_Rkind,Pi/TWO,             &
                 ZERO,ZERO,1.83492_Rkind,1.77157_Rkind,Pi/TWO]

    IF (debug) write(out_unitp,*) 'init d0GGdef of HNO3'
    CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)


    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_QML_HNO3
!> @brief Subroutine wich prints the current QML_HNO3 parameters.
!!
!! @param QModel            CLASS(QML_HNO3_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_HNO3(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_HNO3_t),  intent(in) :: QModel
    integer,              intent(in) :: nio

  END SUBROUTINE Write_QML_HNO3
!> @brief Subroutine wich prints the default QML_HNO3 parameters.
!!
!! @param QModel            CLASS(QML_HNO3_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write0_QML_HNO3(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_HNO3_t),  intent(in) :: QModel
    integer,              intent(in) :: nio

    write(nio,*) 'HNO3 default parameters'
    write(nio,*)
    write(nio,*)
    write(nio,*) 'end HNO3 default parameters'


  END SUBROUTINE Write0_QML_HNO3

!> @brief Subroutine wich calculates the HNO3 potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(QML_HNO3_t):   derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_HNO3(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE ADdnSVM_m
  IMPLICIT NONE

    CLASS(QML_HNO3_t),  intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv



    TYPE (dnS_t)    :: dnDQ(QModel%ndim-1),vh,rot
    integer         :: i1,i2

    rot  = dnQ(1)

    DO i1=1,QModel%ndim-1
      dnDQ(i1) = dnQ(i1+1) - QML_dnvfour_HNO3(rot,i1,0,QModel)
    END DO

    vh = ZERO
    DO i1=1,QModel%ndim-1
    DO i2=1,QModel%ndim-1
      vh = vh + dnDQ(i1)*dnDQ(i2) * QML_dnvfour_HNO3(rot,i1,i2,QModel)
    END DO
    END DO
    Mat_OF_PotDia(1,1) = QML_dnvfour_HNO3(rot,0,0,QModel) + vh * HALF

  END SUBROUTINE EvalPot_QML_HNO3

  SUBROUTINE EvalFunc_QML_HNO3(QModel,Func,dnQ,nderiv)
  USE ADdnSVM_m
  IMPLICIT NONE

    CLASS(QML_HNO3_t),  intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Func(:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv



    TYPE (dnS_t)    :: rot
    integer         :: i1,i2,ifunc

    rot   = dnQ(1)
    ifunc = 0

    DO i1=1,QModel%ndim-1
      ifunc = ifunc + 1
      Func(ifunc) = QML_dnvfour_HNO3(rot,i1,0,QModel)
    END DO

    DO i1=1,QModel%ndim-1
    DO i2=1,QModel%ndim-1
      ifunc = ifunc + 1
      Func(ifunc) = QML_dnvfour_HNO3(rot,i1,i2,QModel)
    END DO
    END DO

  END SUBROUTINE EvalFunc_QML_HNO3

  FUNCTION QML_dnvfour_HNO3(rot,iq,jq,QModel)
  USE ADdnSVM_m
  IMPLICIT NONE

    TYPE (dnS_t)                        :: QML_dnvfour_HNO3

    CLASS(QML_HNO3_t),  intent(in)    :: QModel
    integer,              intent(in)    :: iq,jq
    TYPE (dnS_t),         intent(in)    :: rot

    integer :: i,kl

    IF (iq > max_fit .OR. iq < 0 .OR. jq > max_fit .OR. jq < 0) THEN
      write(out_unitp,*) ' ERROR in QML_dnvfour_HNO3'
      write(out_unitp,*) ' wrong value for iq or jq E [0...8]',iq,jq
      STOP ' ERROR in QML_dnvfour_HNO3: wrong value for iq or jq E [0...8]'
    END IF

    QML_dnvfour_HNO3 = ZERO
    DO i=1,QModel%nn(0,iq,jq)
      IF (QModel%nt(iq,jq) == 41) kl = 2*i-1    ! cos(kl.x)
      IF (QModel%nt(iq,jq) == 42) kl = 4*i-3    ! cos(2.kl.x)
      IF (QModel%nt(iq,jq) == 51) kl = 2*i-0    ! sin(kl.x)
      IF (QModel%nt(iq,jq) == 52) kl = 4*i-0    ! sin(2.kl.x)

      QML_dnvfour_HNO3 = QML_dnvfour_HNO3 + QModel%F(i,iq,jq)*dnFourier(rot,kl)

    END DO

  END FUNCTION QML_dnvfour_HNO3
  SUBROUTINE QML_read_para4d_HNO3(F,n,ndim,nt,max_points,nom1,exist)
  IMPLICIT NONE

   integer,           intent(in)    :: max_points,ndim
   integer,           intent(inout) :: n(0:ndim),nt
   real (kind=Rkind), intent(inout) :: F(max_points)
   character (len=*), intent(in)    :: nom1
   logical,           intent(inout) :: exist

   integer :: no,ios,kl,i

   write(out_unitp,*) 'QML_read_para4d_HNO3: nom1,max_points: ',nom1,max_points


   CALL file_open2(name_file=nom1,iunit=no,lformatted=.TRUE.,           &
                   old=.TRUE.,err_file=ios)
   IF (ios == 0) THEN

     read(no,*) nt
     read(no,*) i ! for nb_fit (not used)

     write(out_unitp,*) 'nom1,nt,ndim: ',nom1,nt,ndim
     read(no,*) n(0:ndim)
     write(out_unitp,*) 'nom1,n ',nom1,n(0:ndim)
     IF (n(0) > max_points) THEN
         write(out_unitp,*) ' ERROR : The number of coefficients (',n(0),') >'
         write(out_unitp,*) '         than max_points (',max_points,')'
         write(out_unitp,*) '         STOP in QML_read_para4d_HNO3'
         STOP 'ERROR in QML_read_para4d_HNO3'
       END IF
       DO kl=1,n(0)
        read(no,*) F(kl)
!       write(out_unitp,*) F(kl)
       END DO
     CLOSE(no)
     exist = .TRUE.
   ELSE
     write(out_unitp,*) 'The file (',nom1,') does not exist !!'
     exist = .FALSE.
   END IF


  END SUBROUTINE QML_read_para4d_HNO3
END MODULE QML_HNO3_m
