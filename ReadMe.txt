=========================================
=========================================
 QuantumModelLib* is a free software under LGPL.
  date: 11/04/2020

    Copyright 2016 David Lauvergnat [1]
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
   From the QuantumModelLib directory, execute make
   the "libpot.a" must be created

   This version works with:
       gfortran 8.3, 9.1, 9.2 (macOS)
       ifort    18
       pgf90    17.10-0
       nagfor   7.0 (macOS)

 2) Link "libpot.a" to your code

   gfortran ....   -L$QuantumModelLib_path -lpot

      QuantumModelLib_path contains the path of the QuantumModelLib


 3) In your fortan code
 3a) Initialization of the potential (see the list below)

        CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)

        where
           ndim       is the number of degree(s) of freedom
           nsurf      is the number of electronic surface(s) (adiabatic or diabatic)
           pot_name   is the name of the potential or model (phenol, Tully, HenonHeiles ...)
           adiabatic  is a logical flag (.TRUE. or .FALSE.)
           option     enables to select a model with several options (Tully ...)

     The list of available models is given below

     Example:
        ndim  = 2
        nsurf = 3
        CALL sub_Init_Qmodel(ndim,nsurf,'phenol',.FALSE.,0)

        It initializes the phenol potential (2D and 3 PES).

        Some extra parameters can be initialized with specific procedures:
          - Phase_Following (default: .TRUE. ) for the adiabatic calculations:
              CALL set_Qmodel_Phase_Following(Phase_Following)
          - ...

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
      !! pot_name  = 'Retinal_JPCB2000'
=========================================
      !! Model for the photo-isomerization of retinal.
      !! pot_name  = 'Retinal_JPCB2000'
      !! ndim      = 2
      !! nsurf     = 2
      !! ref:  S. Hahn, G. Stock / Chemical Physics 259 (2000) 297-312.
      !!              doi: 10.1016/S0301-0104(00)00201-9'
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
      !! option    = 1 (LSTH)
      !! ndim      = 3   (the 3 H-H distances)
      !! nsurf     = 1
      !! refs (option=1):
      !! P. Siegbahn, B. Liu,  J. Chem. Phys. 68, 2457(1978).
      !! D.G. Truhlar and C.J. Horowitz, J. Chem. Phys. 68, 2466 (1978); https://doi.org/10.1063/1.436019
=========================================
=========================================
=========================================
=========================================
=========================================
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
      !! ref:  Eq 3.79 of J. Beutier, thesis.
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
