=========================================
=========================================
 QuantumModelLib is a free software under LGPL.
  date: 01/04/2020

    Copyright 2016 David Lauvergnat
      with contributions of Félix MOUHAT and Liang LIANG
=========================================
 1) Installation
   From the QuantumModelLib directory, execute make
   the "libpot.a" must be created

   This version works with:
       gfortran 8.3, 9.1, 9.2 (macOS)
       ifort    17, 18
       pgf90    17.10-0
       nagfor   7.0

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
        CALL sub_Init_Qmodel(2,3,'phenol',.FALSE.,0)

        It initializes the phenol potential (2D and 3 PES).


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


 3c) Get the metric tensor, GGdef

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
           GGdef(:,)  is the new metric tansor a ndim x ndim matrix of real (kind=8)

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
      !! The default is option=1 (unpublished).
      !! The parameters for option=2 come from the following reference.
      !! ref: E. Marsili, M. H. Farag, X. Yang, L. De Vico, and M. Olivucci, JPCA, 123, 1710–1719 (2019). https://doi.org/10.1021/acs.jpca.8b10010
=========================================
=========================================
