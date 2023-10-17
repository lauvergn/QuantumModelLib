USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : real32,real64,real128
IMPLICIT NONE

integer, parameter :: Rk = real64
real (kind=Rk) :: V
real (kind=Rk) :: Q(3)
real (kind=Rk) :: Qsym(3)

Qsym = [ 1.3_Rk,1.2_Rk, -0.5_Rk]
write(*,*) 'Qsym:',Qsym

CALL QML_ClH2p_Qsym_CCSDTF12(V,Qsym)
write(*,*) 'pot (Qsym)',V

Q = [1.7_Rk, 0.7_Rk, 1.3_Rk]
write(*,*) 'Q   :',Q
CALL QML_ClH2p_CCSDTF12(V,Q)
write(*,*) 'pot (Q   )',V
END
  !==================================================================
  !==================================================================
  !  Potential : ClH2+
  !
  !  Abinitio level : CCSD(T)-F12b/cc-pVTZ-F12
  !
  !  Two versions :
  !   1) QML_ClH2p_CCSDTF12(V,Q)
  !      Q(1:3) contains the three coordinates [r1,r2,theta], units [bohr,bohr,radian]
  !      V      is the energy in Hartree
  !      Example, for Q=[1.7, 0.7, 1.3], V = -456.16936282793722
  !   2) QML_ClH2p_Qsym_CCSDTF12(V,Q)
  !      Q(1:3) contains the three coordinates [theta, r+, r-], units [radian,bohr,bohr]
  !             r+ = 1/2(r1+r2) and r- = 1/2(r1-r2)
  !      V      is the energy in Hartree
  !      Example, for Q=[1.3, 1.2, -0.5, 1.3], V = -456.16936282793722
  !
  !  ref: Afansounoudji KMR, Issa R, Sodoga K, Lauvergnat D. 
  !       ChemRxiv. 2023, DOI: 10.26434/chemrxiv-2023-b4hhm
  !==================================================================
  !==================================================================
  SUBROUTINE QML_ClH2p_CCSDTF12(V,Q)
    USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : real32,real64,real128
    IMPLICIT NONE

    integer, parameter :: Rk = real64
    real (kind=Rk),    intent(inout) :: V
    real (kind=Rk),    intent(in)    :: Q(3)

    real (kind=Rk)    :: Qsym(3)

    Qsym(1) = Q(3)
    Qsym(2) = ( Q(1) + Q(2)) * 0.5_Rk
    Qsym(3) = ( Q(1) - Q(2)) * 0.5_Rk

    CALL QML_ClH2p_Qsym_CCSDTF12(V,Qsym)

  END SUBROUTINE QML_ClH2p_CCSDTF12
  SUBROUTINE QML_ClH2p_Qsym_CCSDTF12(V,Qsym)
    USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : real32,real64,real128
    IMPLICIT NONE
  
    integer, parameter :: Rk = real64
    integer, parameter :: ndim = 3

    real (kind=Rk),    intent(inout) :: V
    real (kind=Rk),    intent(in)    :: Qsym(ndim)

    real (kind=Rk), parameter :: Qref(ndim)  = [ &
      1.6459159559216601_Rk, 2.4673412574786124_Rk,0._Rk]

    real (kind=Rk), parameter :: E0 = -460.59052688_Rk

    real (kind=Rk), parameter :: F1(22) = [ &
  2.8367497759041187E-005_Rk,            &
  6.0176972984675098E-006_Rk,            &
 0.28365120303296165_Rk,                 &
 0.28035409261950595_Rk,                 &
  8.3854571559202326E-002_Rk,            &
 -0.27841948015477108_Rk,                 &
 -1.0056865165624051E-002_Rk,            &
 0.17983180018241626_Rk,                 &
 0.18096613130206318_Rk,                 &
 -6.1656237616369422E-003_Rk,            &
 -9.5180546513037839E-002_Rk,            &
  2.1243816535040125E-003_Rk,            &
  4.4652067648935506E-002_Rk,            &
  4.4234367053673118E-002_Rk,            &
 -8.6829390359447808E-003_Rk,            &
 -1.9348848223610559E-002_Rk,            &
  4.8031743074474680E-004_Rk,            &
  8.8885289235987458E-003_Rk,            &
  9.7949638611087610E-003_Rk,            &
  8.1592825661224878E-003_Rk,            &
 -3.7570777670436089E-003_Rk,            &
 -4.9393257715852519E-003_Rk]

  integer, parameter :: ind1(ndim,22) = reshape(   &
                                            [0,1,0,&
                                             1,0,0,&
                                             0,2,0,&
                                             0,0,1,&
                                             2,0,0,&
                                             0,3,0,&
                                             3,0,0,&
                                             0,4,0,&
                                             0,0,2,&
                                             4,0,0,&
                                             0,5,0,&
                                             5,0,0,&
                                             0,6,0,&
                                             0,0,3,&
                                             6,0,0,&
                                             0,7,0,&
                                             7,0,0,&
                                             0,8,0,&
                                             0,0,4,&
                                             8,0,0,&
                                             0,9,0,&
                                             9,0,0],shape=[ndim,22])
  
  real (kind=Rk), parameter :: F2(52) = [ &
  1.0804766937148311E-002_Rk,            &
  7.3203164207750653E-003_Rk,            &
-0.83498706938487566_Rk,            &
 -1.0920616366581173E-002_Rk,            &
 -5.1503804094965404E-002_Rk,            &
  2.6404913597974744E-003_Rk,            &
 -2.8519141629125149E-002_Rk,            &
  1.0775951257136278_Rk,            &
  9.2247277977610849E-003_Rk,            &
  1.2636146891161401E-002_Rk,            &
-0.96075554344298708_Rk,            &
-0.48342782102430981_Rk,            &
 -7.2796742459972164E-003_Rk,            &
  4.5947019212529798E-003_Rk,            &
 -5.1523267597268420E-004_Rk,            &
  1.4154357719871602E-003_Rk,            &
 -5.2588622165617022E-003_Rk,            &
  1.0559706003973223E-002_Rk,            &
 -2.9933865337583966E-003_Rk,            &
  7.4320893409380905E-003_Rk,            &
 0.71555227553533662_Rk,            &
 -2.0852229973068971E-003_Rk,            &
 -2.0023290464512866E-003_Rk,            &
  1.2154077019357550E-002_Rk,            &
 0.73028053963925355_Rk,            &
 -4.2339575097439033E-003_Rk,            &
  1.7350635931079811E-003_Rk,            &
  9.5765267729195185E-003_Rk,            &
-0.36881747572940649_Rk,            &
  4.4057336952482542E-003_Rk,            &
  5.4209837336978645E-003_Rk,            &
-0.60133810347500005_Rk,            &
  6.0900933952191001E-003_Rk,            &
-0.11902943953977349_Rk,            &
  4.0046406882027754E-005_Rk,            &
  2.5955417672450663E-003_Rk,            &
  2.4874659869247816E-002_Rk,            &
  4.3860233465482150E-004_Rk,            &
 -1.3160876734420545E-002_Rk,            &
 -2.0628511197313425E-002_Rk,            &
 -4.0064321891299798E-003_Rk,            &
 -5.9991821921908543E-002_Rk,            &
 -3.5624915901212029E-002_Rk,            &
 -7.7492132483541146E-003_Rk,            &
  2.4890681604352307E-003_Rk,            &
 -6.5472280451834538E-003_Rk,            &
 -4.7324278579810553E-004_Rk,            &
 -5.1564894861429926E-004_Rk,            &
 -4.4426839629724667E-003_Rk,            &
  3.3233926954765514E-002_Rk,            &
  1.8807571632944518E-003_Rk,            &
  1.0667467015363960E-002_Rk]

  integer, parameter :: ind2(ndim,52) = reshape(  &
                                          [1,1,0, &
                                           1,0,1, &
                                           0,1,1, &
                                           1,2,0, &
                                           2,1,0, &
                                           1,3,0, &
                                           2,0,1, &
                                           0,2,1, &
                                           2,2,0, &
                                           3,1,0, &
                                           0,3,1, &
                                           0,1,2, &
                                           1,0,2, &
                                           2,3,0, &
                                           3,0,1, &
                                           1,4,0, &
                                           3,2,0, &
                                           4,1,0, &
                                           2,4,0, &
                                           2,0,2, &
                                           0,4,1, &
                                           1,5,0, &
                                           3,3,0, &
                                           4,0,1, &
                                           0,2,2, &
                                           4,2,0, &
                                           5,1,0, &
                                           1,0,3, &
                                           0,5,1, &
                                           1,6,0, &
                                           3,4,0, &
                                           0,3,2, &
                                           2,5,0, &
                                           0,1,3, &
                                           4,3,0, &
                                           5,0,1, &
                                           3,0,2, &
                                           5,2,0, &
                                           6,1,0, &
                                           4,0,2, &
                                           1,7,0, &
                                           0,4,2, &
                                           0,2,3, &
                                           4,4,0, &
                                           2,0,3, &
                                           2,6,0, &
                                           3,5,0, &
                                           5,3,0, &
                                           6,0,1, &
                                           0,6,1, &
                                           6,2,0, &
                                           7,1,0],shape=[ndim,52])

  real (kind=Rk), parameter :: F3(7) = [ &
  -1.3512518908906826E-002_Rk,            &
   1.1540979493809583E-002_Rk,            &
   2.4854218086488685E-002_Rk,            &
  -5.1957704614613897E-003_Rk,            &
   5.8964351164980380E-003_Rk,            &
  -6.1758760082113694E-003_Rk,            &
  -1.2272880132086163E-002_Rk]

 integer, parameter :: ind3(ndim,7) = reshape(  &
                                         [1,1,1, &
                                          1,2,1, &
                                          2,1,1, &
                                          1,3,1, &
                                          1,1,2, &
                                          2,2,1, &
                                          3,1,1],shape=[ndim,7])

  real (kind=Rk), allocatable :: DQ(:,:)
  real (kind=Rk)              :: Vtemp
  integer                        :: i,j,max_exp

  !write(out_unit,*) ' sub QML_ClH2p_Qsym_CCSDTF12' ; flush(6)
  max_exp = maxval(ind1)
  !write(out_unit,*) ' max_exp',max_exp ; flush(6)

  allocate(DQ(ndim,0:max_exp))
  DQ(:,0) = 1._Rk
  DQ(:,1) = Qsym(:) - Qref(:)
  DQ(3,1) = DQ(3,1)**2 ! because for R-, we have even values
  DO j=2,ubound(DQ,dim=2)
    DQ(:,j) = DQ(:,j-1) * DQ(:,1)
  END DO
  !write(out_unit,*) ' DQ done' ; flush(6)


  V = E0
  !write(out_unit,*) 'V0',V

  DO j=1,size(ind1,dim=2)
    !write(out_unit,*) j,':',ind1(:,j),F1(j) 
    Vtemp = F1(j) 
    DO i=1,ndim
      Vtemp = Vtemp * DQ(i,ind1(i,j))
    END DO
    V = V + Vtemp
  END DO
  
  DO j=1,size(ind2,dim=2)
    Vtemp = F2(j) 
    DO i=1,ndim
      Vtemp = Vtemp * DQ(i,ind2(i,j))
    END DO
    V = V + Vtemp
  END DO

  DO j=1,size(ind3,dim=2)
    Vtemp = F3(j) 
    DO i=1,ndim
      Vtemp = Vtemp * DQ(i,ind3(i,j))
    END DO
    V = V + Vtemp
  END DO
  !write(out_unit,*) ' end QML_ClH2p_Qsym_CCSDTF12' ; flush(6)

  END SUBROUTINE QML_ClH2p_Qsym_CCSDTF12
