=========================================
=========================================
 QuantumModelLib* is a free software under the MIT Licence.
  date: 25/11/2022

    Copyright (c) 2022 David Lauvergnat [1]
      with contributions of:
        Félix MOUHAT [2]
        Liang LIANG [3]
        Emanuele MARSILI [1,4]

[1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
[2]: Laboratoire PASTEUR, ENS-PSL-Sorbonne Université-CNRS, France
[3]: Maison de la Simulation, CEA-CNRS-Université Paris-Saclay,France
[4]: Durham University, Durham, UK
* Originally, it has been developed during the Quantum-Dynamics E-CAM project :
     https://www.e-cam2020.eu/quantum-dynamics
=========================================
 1) Installation
   From the QuantumModelLib directory, when make is executated, 
   the "libQMLib_XXX_OMPy.a" ("libpot.a") and "libAD_dnSVM_XXX_OMPy.a" must be created.
   XXX is the compiler name and y is 0/1 when OpenMP flag is turn off/on. 

   This version works with:
       gfortran 9.0 (linux and macOS)
       ifort    19
       pgf90    ?
       nagfor   ?

 2) Link "libAD_dnSVM.a" to your code

   gfortran ....   -L$QuantumModelLib_path --LQMLib_XXX_OMPy -lAD_dnSVM_XXX_OMPy

      QuantumModelLib_path contains the path of the QuantumModelLib


 3) In your fortan code
 3a1) Initialization of the model (the Potential)

        CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)

        where
           ndim       : the number of degree(s) of freedom [integer]
           nsurf      : the number of electronic surface(s) (adiabatic or diabatic) [integer]
           pot_name   : the name of the potential or model (phenol, Tully, HenonHeiles ...) [string of characters]
           adiabatic  : flag (.TRUE. or .FALSE.) [logical]
           option     : option, to be able to select a model with several options (Tully ...) [integer]

     The list of available models is given below

     Example:
        ndim  = 2
        nsurf = 3
        CALL sub_Init_Qmodel(ndim,nsurf,'phenol',.FALSE.,0)

        It initializes the phenol potential (2D and 3 PES).
        => Computation of the diabatic surface


        See the list of models below.

 3a2) Initialization of the potential (reading the model)

 CALL sub_Read_Qmodel(ndim,nsurf,nio)

 where
    ndim       : the number of degree(s) of freedom [integer]
    nsurf      : the number of electronic surface(s) (adiabatic or diabatic) [integer]
    nio        : file unit where the namelist is read. It can be the standard unit [integer]

  Then, the &potential namelist is read.
  In the following exemple, the 2+1D-retinal model ('Retinal_JPCB2000') is read.

  &potential
    pot_name='Retinal_JPCB2000' ! potential surface name
    ndim=3 PubliUnit=f
    adiabatic=t
    Phase_checking=f
     /

     It initializes the 2+1D-retinal model (ndim=3).
     For this model, fhe number of electronic surfaces is automatically set up to 2.
     => adiabatic=t      : Computation of the adiabatic surface: 
     => Phase_checking=f : The adiabatic vector phases are not checked between several calculations
     => PubliUnit=f      : The atomic units are used

 3a3) Initialization (extra)
    Some extra parameters can be initialized with specific procedures:

    - Phase_Following (default: .TRUE. ) for the adiabatic calculations:
              CALL set_Qmodel_Phase_Following(Phase_Following)

    - Phase_Checking (default: .TRUE. ) for the adiabatic calculations:
              CALL set_Qmodel_Phase_Checking(Phase_Checking)

    - ...

    or it can be set up while reading the model (see 3a2).

 3b) Potential energy surface(s), PES, evaluation

        CALL sub_Qmodel_V(V,Q)

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

 3c) Potential energy surface with vibrational adiabatic separation
    This feature can be used only when the "model" is read.
    Therefore in the initialization with sub_Init_Qmodel "pot_name" must be "read_model".
    Then:

    (i) The potential must be read as a namelist:

      &potential
          pot_name='hbond' ! potential surface name
          Vib_adia=t
          list_act=1
          read_nml=f
          nb_channels=6
      /

    In this example the 'Hbond' potential is used within the adiabatic separation (Vib_adia=t).
    Read_nml=f, the specific namelist for 'hbond' model is not read.
    The number of channels (nsurf) is 6
    list_act is a table which enables to select the active coordinate(s) (here the first one only)
    The inactive coordinates are just the remaining coordinates (here 2)

    (ii) A basis set must be read for the inactive coordinates:

    &basis_nD name='boxAB' nb=64 nq=64 A=-2.5 B=2.5 /

    Here, it is a 1D basis set (particle-in-a-box), with
      64 basis function (nb)
      64 grid points (nq)
      The range of the coordinate is [A,B]

    WARNING: It is working only with ONE inactive variable.

    The following subroutine enables to get the effective Hamiltonian along Qact(:).
      CALL sub_Qmodel_tab_HMatVibAdia(tab_MatH,Q,nb_terms)
    The table tab_MatH(nsurf,nsurf,nb_terms) contains:
      - Heff: tab_MatH(nsurf,nsurf,1)                 [1 matrix]
      - F2 terms: tab_MatH(nsurf,nsurf,2:...)         [ (nb_act+1)nb_act/2 matrices)]
      - F1 terms: tab_MatH(nsurf,nsurf,...:nb_terms)  [ nb_act matrices ]

 3d) Get the metric tensor, GGdef

        CALL get_Qmodel_GGdef(GGdef)

        where
           GGdef(:,)  is the ndim x ndim matrix of real (kind=8)

     Most of the models are associated with a specific kinetic energy operator (KEO)
     Here, it is given through a constant metric tensor, GGdef, so the that the KEO is:

        T = -1/2 Sum_ij d./dQi GGdef(i,j) d./dQj

     Its diagonal components, GGdef(i,i), can be view as the invers of masses (1/Mi)

     The volume element, dTau, is:

        dTau = dQ1.dQ2 ... dQndim

     Remark: The metric tensor can be modified:
        CALL set_Qmodel_GGdef(GGdef,ndim)

        where
           ndim       is the number of degree(s) of freedom it MUST be indentical to the initialized value (with sub_Init_Qmodel)
           GGdef(:,)  is the new metric tensor a ndim x ndim matrix of real (kind=8)

=========================================
=========================================
=========================================
=========================================
      !! pot_name  = 'Morse'
=========================================
      !! Morse potential: V(R) = D*(1-exp(-a*(r-Req))**2
      !! pot_name  = 'Morse'
      !! ndim      = 1
      !! nsurf     = 1
      !! reduced mass      = 1744.60504565084306291455 au
      !! remark: Default parameters for H-F
=========================================
=========================================
=========================================
=========================================
      !! pot_name  = 'Poly1D'
=========================================
      !! Polynomial potential: V(R) = sum_i coef(i) * (r-Req)**i
      !! pot_name  = 'Poly1D'
      !! ndim      = 1
      !! nsurf     = 1
      !! reduced mass      = 1744.60504565084306291455 au
      !! remark: Default parameters for H-F
=========================================
=========================================
=========================================
=========================================
      !! pot_name  = 'Phenol'
=========================================
      !! Phenol model
      !! pot_name  = 'Phenol'
      !! ndim      = 2 (R=rOH, th=OH-torsion)
      !! nsurf     = 3
      !! Diagonal Metric Tensor(:)      = (/ 0.0005786177, 0.0002550307 /) au
      !! remark:
      !! ref: Z. Lan, W. Domcke, V. Vallet, A.L. Sobolewski, S. Mahapatra, J. Chem. Phys. 122 (2005) 224315. doi:10.1063/1.1906218.
=========================================
=========================================
=========================================
=========================================
      !! pot_name  = 'TwoD'
=========================================
      !! 2D model
      !! pot_name  = 'TwoD'
      !! ndim      = 2 (X,Y)
      !! nsurf     = 2
      !! Reduced masses(:)      = (/ 20000., 6667. /) au
      !! remark: The parameter values have been modified
      !! ref: A. Ferretti, G. Granucci, A. Lami, M. Persico, G. Villani, J. Chem. Phys. 104, 5517 (1996); https://doi.org/10.1063/1.471791
=========================================
=========================================
=========================================
=========================================
      !! pot_name  = 'TwoD_RJDI2014'
=========================================
      !! 2D model
      !! pot_name  = 'TwoD_RJDI2014'
      !! ndim      = 2 (X,Y)
      !! nsurf     = 2
      !! Reduced masses(:)      = [1. , 1.] au
      !! ref:  Ilya G. Ryabinkin, Loïc Joubert-Doriol, and Artur F. Izmaylov, ...
      !!       ... J. Chem. Phys. 140, 214116 (2014); https://doi.org/10.1063/1.4881147
=========================================
=========================================
=========================================
=========================================
      !! pot_name  = 'PSB3'
=========================================
      !! Model for the photo-isomerization of the penta-2,4-dieniminium (PSB3) cation.
      !! pot_name  = 'PSB3'
      !! ndim      = 3
      !! nsurf     = 2
      !! remarks: two options are possible (option = 1,2)
      !! The default is option=1 (ref2).
      !! The parameters for option=2 come from the following reference.
      !! ref1: E. Marsili, M. H. Farag, X. Yang, L. De Vico, and M. Olivucci, JPCA, 123, 1710–1719 (2019).
      !!         https://doi.org/10.1021/acs.jpca.8b10010
      !! ref2: 1 E. Marsili, M. Olivucci, D. Lauvergnat, and F. Agostini, JCTC 16, 6032 (2020).
      !!        https://pubs.acs.org/doi/10.1021/acs.jctc.0c00679
=========================================
=========================================
=========================================
=========================================
      !! pot_name  = 'Retinal_JPCB2000'
=========================================
      !! Model for the photo-isomerization of retinal.
      !! pot_name  = 'Retinal_JPCB2000'
      !! ndim      = 2 or up to 2+23
      !! nsurf     = 2
      !! ref:  S. Hahn, G. Stock / Chemical Physics 259 (2000) 297-312.
      !!              doi: 10.1016/S0301-0104(00)00201-9
=========================================
=========================================
=========================================
=========================================
      !! pot_name  = 'HONO'
=========================================
      !! Model for the HONO.
      !! pot_name  = 'HONO'
      !! ndim      = 6
      !! nsurf     = 1
      !! ref1:  F. Richter, M. Hochlaf, P. Rosmus, F. Gatti, and H.-D. Meyer,
      !!        J. Chem. Phys. 120, 1306 (2004).
      !!       doi: 10.1063/1.1632471
      !! ref2: F. Richter, F. Gatti, C. Léonard, F. Le Quéré, and H.-D. Meyer,
      !!       J. Chem. Phys. 127, 164315 (2007)
      !!       doi: 10.1063/1.2784553
=========================================
=========================================
=========================================
=========================================
      !! pot_name  = 'CH5'
=========================================
      !! H + CH4 -> H-H + CH3 potential
      !!    Quadratic potential along the reaction path'
      !!    Reaction coordinate: R- = 1/2(RCH - RHH)'
      !!    Optimal coordinates along the path at CCSD(T)-F12/cc-pVTZ-F12'
      !!    V0 along the path at CCSD(T)-F12/cc-pVTZ-F12'
      !!    Hessian along the path at MP2/cc-pVDZ'
      !! pot_name  = 'CH5'
      !! ndim      = 12 or 1
      !! nsurf     = 1
      !! option = 4 (default) or 5
=========================================
=========================================
=========================================
=========================================
      !! pot_name  = 'PH4'
=========================================
      !! H + PH3 -> H-H + PH2 potential
      !!    Quadratic potential along the reaction path'
      !!    Reaction coordinate: R- = 1/2(RPH - RHH)'
      !!    Optimal coordinates along the path at MP2/cc-pVTZ'
      !!    V0 along the path at CCSD(T)-F12/cc-pVTZ-F12 (option 4) or ...'
      !!      ... MP2/cc-pVTZ'
      !!    Hessian and gradient along the path at MP2/cc-pVTZ'
      !! pot_name  = 'PH4'
      !! ndim      = 9 or 1
      !! nsurf     = 1
      !! option = 4 (default)
=========================================
=========================================
=========================================
=========================================
      !! pot_name  = 'HOO_DMBE'
=========================================
      !! HOO potential: DMBE IV of Varandas group
      !! pot_name  = 'HOO_DMBE'
      !! ndim      = 3   (R1=dOO,R2=dHO1,R3=dHO2)
      !! nsurf     = 1
      !! ref:    M. R. Pastrana, L. A. M. Quintales, J. Brandão and A. J. C. Varandas'
      !!         JCP, 1990, 94, 8073-8080, doi: 10.1021/j100384a019.
=========================================
=========================================
=========================================
=========================================
      !! pot_name  = 'H3'
=========================================
      !! H3 potential:
      !! pot_name  = 'H3'
      !! option    = 0,1,10,11 (LSTH)
      !! ndim      = 3   (the 3 H-H distances)
      !! nsurf     = 1
      !! refs (option=0):
      !! P. Siegbahn, B. Liu,  J. Chem. Phys. 68, 2457(1978).
      !! D.G. Truhlar and C.J. Horowitz, J. Chem. Phys. 68, 2466 (1978); https://doi.org/10.1063/1.436019
      !! options  0 and 10 : 3D model with IRC functions (1 potential + parameters)
      !! options  0 and  1 : first IRC funtions fitted in polar representation (alpha)
      !! options 10 and 11 : second IRC funtions fitted with the sum and the difference (alpha)
=========================================
=========================================
=========================================
=========================================
      !! pot_name  = 'H2'
=========================================
      !! H2 potential: V(R) = Sum_i a_i * (R-Req)**(i-1)
      !! pot_name  = 'H2'
      !! ndim      = 1
      !! nsurf     = 1
      !! reduced mass      = 1837.1526464003414/2 au
      !! Level: CCSD(T)-F12B/VTZ-F12 (with molpro 2010)
=========================================
=========================================
=========================================
=========================================
      !! pot_name  = 'HCN_Murrell'
=========================================
      !! H3 potential:
      !! pot_name  = 'HCN_Murrell'
      !! option    = 0 (3D-3distances), 1 (3D-Jacobi), 3 IRC+RPH
      !! ndim      = 3
      !! nsurf     = 1
      !! J. N. Murrell, S. Carter and L. O. Halonene, J. Mol. Spectrosc. vo93 p307 1982
      !!  doi: https://doi.org/10.1016/0022-2852(82)90170-9
=========================================
=========================================
=========================================
=========================================
      !! pot_name  = 'H2O'
=========================================
      !! H2O potential:
      !! pot_name  = 'H2O'
      !! option    = 1(default 1)
      !! ndim      = 3
      !! nsurf     = 1
      !! refs: Quadratic model potential for H2O; TIPS force constants taken from:  
      !!       Dang and Pettitt, J. Chem. Phys. 91 (1987)
=========================================
=========================================
=========================================
=========================================
      !! pot_name  = 'ClH2+'
=========================================
      !! ClH2+ potential:
      !! pot_name  = 'ClH2+'
      !! option    = 1 to 6 (default 6)
      !! ndim      = 3
      !! nsurf     = 1
      !! remark, 6 options are possible:
      !!    options = 1,3,5: the coordinates are [angle, R+, R-] (in bohr and radian)
      !!    options = 2,4,6: the coordinates are [R1, R2, angle] (in bohr and radian)
      !!
      !!    options = 1,2: B3LYP/cc-pVTZ 1st version (do not use)
      !!    options = 3,4: B3LYP/cc-pVTZ 2d  version
      !!    options = 5,6: CCSD(T)-F12b/cc-pVTZ-F12
      !! refs: unpublished
=========================================
=========================================
=========================================
=========================================
      !! pot_name  = 'ClH2+_Botschwina'
=========================================
      !! ClH2+ potential:
      !! pot_name  = 'ClH2+_Botschwina'
      !! option    = 1,2 (default 1)
      !! ndim      = 3
      !! nsurf     = 1
      !! remark, 2 options are possible:
      !!    options = 1: CEPA-1
      !!    options = 2: CEPA-1 corrected
      !! ref: Peter Botschwina, J. Chem. Soc., Faraday Trans. 2, 1988, 84(9), 1263-1276'
      !!      DOI: 10.1039/F29888401263
=========================================
=========================================
=========================================
=========================================
=========================================
=========================================
=========================================
=========================================
=========================================
      !! pot_name  = 'buck'
=========================================
      !! Buckingham potential: V(R) = A*exp(-B*R)-C/R^6
      !! pot_name  = 'buck'
      !! ndim      = 1
      !! nsurf     = 1
      !! reduced mass      = 36423.484024390622 au
      !! remark: default parameters for Ar2
      !! ref:  R.A. Buckingham, Proc. R. Soc. A Math. Phys. Eng. Sci. 168 (1938) 264–283. doi:10.1098/rspa.1938.0173
=========================================
=========================================
=========================================
=========================================
      !! pot_name  = 'HBond'
=========================================
      !! LinearHBond potential: Morse1(QQ/2+q,param1)+Morse2(QQ/2-q,param2)+Eref2+Buckingham(QQ)
      !! pot_name  = 'HBond'
      !! ndim      = 2   (QQ,q)
      !! nsurf     = 1
      !! reduced masses      = (/ 29156.946380706224, 1837.1526464003414 /) au
      !! remark:
      !!    A--------------H-----X--------------------B
      !!     <--------------QQ----------------------->
      !!                    <-q->
      !! ref: Dana Codruta Marinica, Marie-Pierre Gaigeot, Daniel Borgis,
      !!    Chemical Physics Letters 423 (2006) 390–394
      !!    DOI: 10.1016/j.cplett.2006.04.007
      !! ref:  Eq 3.79 of J. Beutier, thesis.
      !!
      !! remark: when option=2 is selected, a contribution is added along QQ:
      !!    Dm.exp(-betam.(QQ-QQm)) + Dp.exp(+betap.(QQ-QQp))
      !!
=========================================
=========================================
=========================================
=========================================
      !! pot_name  = '2d_mb'
=========================================
      !! 2D Müller-Brown potential:
      !! pot_name  = '2d_mb'
      !! ndim      = 2   (x,y)
      !! nsurf     = 1
      !! reduced masses      = [1000.,1000.]
      !!
      !! ref:   Klaus Müller and Leo D. Brown, ...
      !!        ... Theoret. Chim. Acta (Berl.) 53, 75-93 (1979)
      !!        https://doi.org/10.1007/BF00547608.
      !!
      !! remark: the option enables one to select among three minima and two TS.
      !!         default (option=1), the first minimum (A in the reference)
      !!
=========================================
=========================================
=========================================
=========================================
      !! pot_name  = 'HenonHeiles'
=========================================
      !! HenonHeiles potential
      !! pot_name  = 'HenonHeiles'
      !! ndim      > 1
      !! nsurf     = 1
      !! reduced masses(:)      = ONE au
      !! ref:  parameters taken from M. Nest, H.-D. Meyer, J. Chem. Phys. 117 (2002) 10499. doi:10.1063/1.1521129
=========================================
=========================================
=========================================
=========================================
      !! pot_name  = 'Tully'
=========================================
      !! Tully potential: three options
      !! pot_name  = 'Tully'
      !! ndim      = 1
      !! nsurf     = 2
      !! reduced mass      = 2000. au
      !! remark: three options are possible (option = 1,2,3)
      !! ref:  Tully, J. Chem. Phys. V93, pp15, 1990
=========================================
=========================================
=========================================
=========================================
      !! pot_name  = '1DSOC_1S1T'
=========================================
      !! Spin Orbit coupling model
      !! pot_name  = '1DSOC_1S1T'
      !! ndim      = 1
      !! nsurf     = 4 or 2
      !! reduced mass      = 20000. au
      !! remarks: 1 singlet and 3 triplet components                           => nsurf     = 4
      !!  or      1 singlet and 1 linear combibation of the triplet components => nsurf     = 2
      !! ref: Giovanni Granucci, Maurizio Persico, and Gloria Spighi, J. Chem. Phys. V137, p22A501 (2012)
=========================================
=========================================
=========================================
=========================================
      !! pot_name  = '1DSOC_2S1T'
=========================================
      !! Spin Orbit coupling model
      !! pot_name  = '1DSOC_2S1T'
      !! ndim      = 1
      !! nsurf     = 4
      !! reduced mass      = 20000. au
      !! remark: 2 singlets and 1 triplet (2 linear combinations of the triplet components are not included) => nsurf     = 4
      !! ref: Giovanni Granucci, Maurizio Persico, and Gloria Spighi, J. Chem. Phys. V137, p22A501 (2012)
=========================================
=========================================
