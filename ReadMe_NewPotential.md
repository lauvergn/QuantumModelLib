
# Procedure to implement a new potential/model

Let assume the potential/model name is **XXX**.
Three (or four) main steps are required:

1. Creation of a Fortran file, **XXX_m.f90**, containing the potential/model in the SRC/QML directory.
All the required subroutines and the new type (**QML_XXX_t**) are stored in this Fortran file.

2. The Fortran file, **Model_m.f90**, in the SRC diretory has to be modified.

3. The **makefile** has to be modified.

4. It is higly recomended to add a procedure to test the potential in the **SRC/TEST_model.f90** file.

## 1) **XXX_m.f90** file

See the SRC/QML/Test_m.f90 file. This Fortran module contains:

### a) **QML_XXX_t** type

This type is associated to the potential/model: **QML_XXX_t**.
It should contains all parameters of the potential (ideally with private atributes)
Alternatively, you can define the potential data as Fortran parameters directly in the module or in the potential subroutine.
Even, if empty, this type MUST be defined.
The type also contains "type-bound procedures" map to generic procedures:

- **EvalPot_QModel**: a subroutine to evaluates the potential at given coordinate values
- **Write_QModel**:   a subroutine to write informations of the model
- **Cart_TO_Q_QModel**: (optional) a subroutine to transform Cartesian coordinates to the ones used in the potential (internal, normal modes ...)

### b) An initialisation function: **Init_QML_XXX**

This function returns **QModel** of type  **QML_XXX_t** and which contains the data of the new potential.
The argument **QModel_in** is needed for the initialization of **QModel**.
**read_param** and **nio_param_file** are, respectively, a flag to read potential extra parameters (or to change some parameters) and the Fortran unit to read those extra parameters.

Several parameters could be defined in the function:

- **QModel%pot_name**: the name of the potential/model (with the example, Model%pot_name='xxx')
- **QModel%ndim**: the number of degrees of freedom. It can be a fixed value or a default value and/or you can test it with respect of acceptable values.
- **QModel%nsurf**: the number of electronic states (diabatic or adiabatic).
- **QModel%Q0[:]**: some coordinate values, mainly to test the potential.
- **QMode%Gdef[:,:]**: A matrix which contains the metric tensor (covariant components) at a given geometry. The diagonal elements can be view as the inverse of an effective masses (the default is the identity matrix).

Remarks: Special features must be added when the transformation between Cartesian coordinates to the ones used in the potential (Cart_TO_Q=.true.) is required. Then, two other parameters are needed:

- **QModel%ndimCart**: the number of Cartesian degrees of freedom.
- **QModel%ndimQ**: the number of degrees of freedom used in the potential.

Then **QModel%ndim** must be defined as follows:
```fortran
IF (QModel%Cart_TO_Q) THEN
  QModel%ndim       = QModel%ndimCart
ELSE
  QModel%ndim       = QModel%ndimQ
END IF
```

### c) A writing subroutine: **Write_QML_XXX**

This subroutine is bounded to the **QML_XXX_t** type and mapped to the generic procedure **Write_QModel**. The dummy arguments are:

- **QModel** of QML_XXX_t type which will print.
- **nio** is Fortran unit to print **QModel**.

### d) The potential subroutine: **EvalPot_QML_XXX**

This subroutine is bounded to the **QML_XXX_t** type and mapped to the generic procedure **EvalPot_QModel**. The dummy arguments are:

- **QModel** of QML_XXX_t type which contains potential parameters.
- **Mat_OF_PotDia**: a matrix of dnS_t type, to store the potential value (matrix of dnS_t type).
- **dnQ**: a vector of dnS_t type which contains coordinate values (and derivatives).
- **nderiv**: an argument to control the derivative order (normally, it is not used).

  For instance, the following expression is diabatic coupling term between the first and the second diabatic potential and it is proportional to the first coordinate (dnQ(1)):

```fortran
  Mat_OF_PotDia(1,2) =  dnQ(1) * QModel%c
```

  Remark: if required, the intermediate variables MUST be of the dnS_t type.

### e) The Cartesian coordinate transformation subroutine: **Cart_TO_Q_QML_Test**

This subroutine is bounded to the **QML_XXX_t** type and mapped to the generic procedure **Cart_TO_Q_QModel**. The dummy arguments are:

- **QModel** of QML_XXX_t type which contains potential parameters.
- **dnX**: a matrix of dnS_t type which contains the Cartesian coordinates (dnX(3,Nat), where Nat is the number of atoms)
- **dnQ**: a vector of dnS_t type which will contain coordinate values (and derivatives).
- **nderiv**: an argument to control the derivative order (normally, it is not used).

The transformation has to be implemented.

## 2) **Model_m.f90** modifications

Two modifications must be made in the **/SCR/Model_m.f90** file.

Add the use of the new module **QML_XXX_m** in the **Init_Model** subroutine.

```fortran
  USE QML_XXX_m
```

In the select case of the **Init_Model** subroutine, add those Fortran lines:

```fortran
  CASE ('xxx') ! the potential/model name MUST in lowercase letters
    !! xxx-potential
    QModel%QM = Init_QML_XXX(QModel_in,read_param=read_nml,nio_param_file=nio_loc)
```

## 3) **makefile** modifications

To create the source file list and the dependencies run:

```bash
  make dep
```