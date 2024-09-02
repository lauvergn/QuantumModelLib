# QuantumModelLib

 QuantumModelLib or QML* is a free software under the MIT Licence.

  date: 01/09/2024

```
    Copyright (c) 2022 David Lauvergnat [1]
      with contributions of:
        Félix MOUHAT [2]
        Liang LIANG [3]
        Emanuele MARSILI [1,4]
        Evaristo Villaseco Arribas [5]
```
[1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
[2]: Laboratoire PASTEUR, ENS-PSL-Sorbonne Université-CNRS, France
[3]: Maison de la Simulation, CEA-CNRS-Université Paris-Saclay,France
[4]: Durham University, Durham, UK
[5]: Department of Physics, Rutgers University, Newark, New Jersey 07102, USA

\* Originally, it has been developed during the Quantum-Dynamics E-CAM project :
     https://www.e-cam2020.eu/quantum-dynamics

## 1) Installation

   From the QuantumModelLib directory, when make is executated, the **libQMLibFull_XXX_optx_ompy_lapakz.a** must be created (ex: **libQMLibFull_gfortran_opt1_omp1_lapack1**).
   **XXX** is the compiler name and **x**, **y** and **z** are 0/1 when flags are turn off/on. 
   They correspond to OPT (compiler optimzation), OpenMP and Lapack/blas, respectively.

```
   This version works with:
       gfortran 9.0 (linux and macOS)
       ifort/ifx    19
```

## 2) Link the library to your code

When lapack/blas are not linked to the library:

```bash
   gfortran ....   $QuantumModelLib_path/libQMLibFull_XXX_optx_ompy_lapak0.a
```

or with  lapack/blas (linux)

```bash
   gfortran ....   $QuantumModelLib_path/libQMLibFull_XXX_optx_ompy_lapak0.a -llapack -lblas
```
*QuantumModelLib_path* contains the path of the **QuantumModelLib**


## 3) In your Fortan code

 In the following, it shows how to initialize, compute with the driver subroutines (the full library is needed, the Fortran module files are not required)

### 3a1) Initialization of the model (the Potential)

```fortran
  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)
```

```
where
  - ndim       : the number of degree(s) of freedom [integer]
  - nsurf      : the number of electronic surface(s) (adiabatic or diabatic) [integer]
  - pot_name   : the name of the potential or model (phenol, Tully, HenonHeiles ...) [string of characters]
  - adiabatic  : flag (.TRUE. or .FALSE.) [logical]
  - option     : option, to be able to select a model with several options (Tully ...) [integer]
```

he list of available models is given below

Example:
```fortran
  ndim  = 2
  nsurf = 3
  CALL sub_Init_Qmodel(ndim,nsurf,'phenol',.FALSE.,0)
```
It initializes the phenol potential (2D and 3 PES).
=> Computation of the diabatic surface


### 3a2) Initialization of the potential (reading the model)

```fortran
 CALL sub_Read_Qmodel(ndim,nsurf,nio)
```
```
 where
  - ndim       : the number of degree(s) of freedom [integer]
  - nsurf      : the number of electronic surface(s) (adiabatic or diabatic) [integer]
  - nio        : file unit where the namelist is read. It can be the standard unit [integer]
```
Then, the **&potential** namelist is read.
In the following exemple, the 2+1D-retinal model ('Retinal_JPCB2000') is read.

```fortran
  &potential
    pot_name='Retinal_JPCB2000' ! potential surface name
    ndim=3 PubliUnit=f
    adiabatic=t
    Phase_checking=f
     /
```

  It initializes the 2+1D-retinal model (ndim=3).
  For this model, fhe number of electronic surfaces is automatically set up to 2.
    => adiabatic=t      : Computation of the adiabatic surface: 
    => Phase_checking=f : The adiabatic vector phases are not checked between several calculations
    => PubliUnit=f      : The atomic units are used

### 3a3) Initialization (extra)

Some extra parameters can be initialized with specific procedures:

- **Phase_Checking** (default: .FALSE. ) for the adiabatic calculations:
  WARNING: Before the version 23.3, the default was .TRUE.

  When Phase_Checking=.TRUE.:
    The eigenvector phases are checked and changed through overlapp beetween the current and reference eigenvectors.
    The two (and only two) eigenvectors are degenerated a rotation may occur (experimental)

  When Phase_Checking=.FALSE.: The eigenvector phases are not checked

  To change Phase_Checking behavior:

```fortran
  CALL set_Qmodel_Phase_Checking(Phase_Checking)
```


- **Phase_Following** (default: .FALSE. ) for the adiabatic calculations:
  WARNING: Before the version 23.2, the default was .TRUE.
      
  This feature is only relevant if Phase_Checking=.TRUE.
  When Phase_Checking=.TRUE.,  the reference eigenvectors are the ones of the previous evaluation. 
   When Phase_Checking=.FALSE., the reference eigenvectors are the ones of the first evaluation. 

  To change Phase_Following behavior: 
```fortran
  CALL set_Qmodel_Phase_Following(Phase_Following)
```

  - Other parameters can be changed

    or it can be set up while reading the model (see 3a2).

### 3b) Potential energy surface(s), PES, evaluation

```fortran
   CALL sub_Qmodel_V(V,Q)
```

```
  where
    Q(:)   is the ndim coordinates (vector of real(kind=8))
    V(:,:) is PES and a  nsurf x nsurf matrix of real (kind=8)
  Remarks:
    - when adiabatic is set to .TRUE., V(:,:) is a diagonal matrix.
    - Use sub_Qmodel_VG(V,G,Q) to get the potential and the gradient
          G(:,:,:) are real (kind=8) of nsurf x nsurf x ndim
    - Use sub_Qmodel_VGH(V,G,H,Q) to get the potential, the gradient and the hessian
          H(:,:,:,:) are real (kind=8) of nsurf x nsurf x ndim x ndim
    - Use sub_Qmodel_VG_NAC(V,G,NAC,Q) to get the potential, the gradient and the
          non-adiabatic couplings (NAC).
          NAC(:,:,:) are real (kind=8) of nsurf x nsurf x ndim
```

### 3c) Potential energy surface with vibrational adiabatic separation

This feature can be used only when the **model** is read.
Therefore in the initialization with sub_Init_Qmodel **pot_name** must be "read_model".
Then:

(i) The potential must be read as a namelist:
```fortran
  &potential
      pot_name='hbond' ! potential surface name
      Vib_adia=t
      list_act=1
      read_nml=f
      nb_channels=6
  /
```

In this example the 'Hbond' potential is used within the adiabatic separation (**Vib_adia=t**).
**Read_nml=f**, the specific namelist for 'hbond' model is not read.
The number of channels (**nsurf**) is 6
**list_act** is a table which enables to select the active coordinate(s) (here the first one only)
The inactive coordinates are just the remaining coordinates (here 2)

(ii) A basis set must be read for the inactive coordinates:

```fortran
  &basis_nD name='boxAB' nb=64 nq=64 A=-2.5 B=2.5 /
```
Here, it is a 1D basis set (particle-in-a-box), with
  - 64 basis function (nb)
  - 64 grid points (nq)
  - The range of the coordinate is [A,B]

WARNING: It is working only with ONE inactive variable.

The following subroutine enables to get the effective Hamiltonian along $Qact(:)$.
```fortran
  CALL sub_Qmodel_tab_HMatVibAdia(tab_MatH,Q,nb_terms)
```
The table tab_MatH(nsurf,nsurf,nb_terms) contains:
  - Heff: tab_MatH(nsurf,nsurf,1)                 [1 matrix]
  - F2 terms: tab_MatH(nsurf,nsurf,2:...)         [ (nb_act+1)nb_act/2 matrices)]
  - F1 terms: tab_MatH(nsurf,nsurf,...:nb_terms)  [ nb_act matrices ]

### 3d) Get the metric tensor, GGdef

```fortran
  CALL get_Qmodel_GGdef(GGdef)
```

where:
  $GGdef(:,:)$  is the ndim x ndim matrix of real (kind=8)

  Most of the models are associated with a specific kinetic energy operator (KEO)
  Here, it is given through a constant metric tensor, GGdef, so the that the KEO is:

  $\hat{T} = -\frac{1}{2} \sum_{ij} \frac{\partial}{\partial Q_i} GGdef(i,j) \frac{\partial}{\partial Q_j}$


Its diagonal components, $GGdef(i,i)$, can be view as the invers of masses ($1/M_i$)
The volume element, $d\tau$, is:

  $d\tau = dQ_1.dQ_2 ... dQ_{ndim}$


Remark: The metric tensor can be modified:
```fortran
  CALL set_Qmodel_GGdef(GGdef,ndim)
```

where
  - ndim       is the number of degree(s) of freedom it MUST be indentical to the initialized value (with sub_Init_Qmodel)
  - $GGdef(:,:)$  is the new metric tensor a ndim x ndim matrix of real (kind=8)

## 4) Examples and Tests

 With "make all", the libraries are created and several main programs.

### 4a) To test the implementation
 From the main QuantumModelLib directory:

```bash
    make ut
```
 From the Tests directory

```bash
    ./run_test_QML
```

  => Some options are possible, the compiler, OPT, OMP, Lapack

  Or you can run:
```bash
    ./run_tests
```

  => All possible combinations between OPT=0/1, OMP=0/1, Lapack=0/1 will be tested. Beware, this test is long.

### 4b) To run examples
From the main QuantumModelLib directory:

```bash
  ./TEST_driver.x < Tests/DAT_files/Vibadia_HBond.dat > res_driver
```
=> Test sevral models (1 or several surfaces, optimization, Vibration adiabatic separation)

```bash
  ./TEST_VibAdia.x < Tests/DAT_files/Vibadia_HBond.dat > res_VibAdia
```
=> Test a Vibration adiabatic separation model.

```bash
  ./TEST_grid.x > res_grid
```
=> Test the 1D and 2D-cut generation for HenonHeiles potential (it uses subroutines with Fortran modules)

```bash
  ./TEST_OMPloop.x > res_loop
```
=> Test an OpenMP loop ($10^6$) on the HONO model (several seconds)

## 5) Installation + examples with fpm

- You have to edit the fpm.toml and change the QML_path (path of QML directory)
- If you want to remove OpenMP feature, 

- build:
```bash
fpm build
```

Remark: if you get this error "'./././/Ext_lib/QDUtilLib/fpm.toml' could not be found"
Its means that the link between some dependencies are lost. Do:
```bash
  make fpmlink
```

- tests:

```bash
  fpm test --> res
```

The intermediate result files (*.txt grid*) are in the RES_files directory (link to Tests/RES_files).
The tests summaries are in QModel.log file (about 43 tests are performed)

Remark: if you get this error "STOP ERROR in Test_QdnV_FOR_Model: Impossible to open the file"
Its means that the link between RES_files and/or DAT_files are lost. Do:
```bash
  make fpmlink
```

- run examples:

```bash
 fpm run driver --< Tests/DAT_files/Vibadia_HBond.dat --> res_driver
```
=> Test sevral models (1 or several surfaces, optimization, Vibration adiabatic separation)

```bash
 fpm run vibadia --< Tests/DAT_files/Vibadia_HBond.dat --> res_VibAdia
 ```
=> Test a Vibration adiabatic separation model.

```bash
  fpm run grid
```
=> Test the 1D and 2D-cut generation for HenonHeiles potential (it uses subroutines with Fortran modules)

```bash
fpm run omp
```
=> Test an OpenMP loop ($10^6$) on the HONO model (several seconds)


# Model list

