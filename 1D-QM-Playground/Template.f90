!> Computes the eigenvalues of a number of states of the particle in the box \$[0,a]\$
!! and compares them with the analytically known result, \$ E_n = (\pi^2n^2/a^2)/2\$ 
!! where \$\hbar=1\$ and \$m=1\$.
program main 

  use, intrinsic :: iso_fortran_env, only : ERROR_UNIT, OUTPUT_UNIT
  use ModuleBSpline

  implicit none

  !.. Bspline parameters and variables
  integer        , parameter   :: BS_NNODS = 51
  ! integer        , parameter   :: BS_NNODS = 101
  ! integer        , parameter   :: BS_NNODS = 201
  ! integer        , parameter   :: BS_NNODS = 501
  integer        , parameter   :: BS_ORDER = 12
  real(kind(1d0)), parameter   :: BS_GRMIN =  0
  ! real(kind(1d0)), parameter   :: BS_GRMAX =  20.d0                 !.. Changed from 1.d0
  ! real(kind(1d0)), parameter   :: BS_GRMAX =  50.d0                 !.. Changed from 1.d0
  ! real(kind(1d0)), parameter   :: BS_GRMAX =  80.d0                 !.. Changed from 1.d0
  real(kind(1d0)), parameter   :: BS_GRMAX =  100.d0                 !.. Changed from 1.d0
  real(kind(1d0)), parameter   :: BS_INTER = BS_GRMAX - BS_GRMIN
  real(kind(1d0))              :: BS_GRID(BS_NNODS)
  type(ClassBSpline)           :: BSpline

  !.. Hamiltonian, overlap, and spectral-decomposition arrays
  real(kind(1d0)), allocatable :: Hmat(:,:), Smat(:,:)
  integer                      :: iNode, nEn, info

  !.. Initializes the BSpline set
  do iNode=1,BS_NNODS
     BS_GRID(iNode) = BS_GRMIN + &
          BS_INTER * dble(iNode-1) / dble(BS_NNODS-1)
  enddo
  call BSpline%Init( &
       BS_NNODS       , &
       BS_ORDER       , &
       BS_GRID        , &
       info           )
  if(info/=0) error stop

  !.. Due to the continuity requirement on the wavefunction,
  !   the first and last BSplines, which do not vanish at the
  !   boundary of the interval, cannot be used
  nEn = BSpline%GetNBsplines() - 2 

  call FillMatrices()

  call GeneralizedEigenproblem()
  
contains

  subroutine FillMatrices()
    integer         :: iBs1, iBs2
    real(kind(1d0)) :: Overlap, KineticEnergy, PotentialEnergy
    real(kind(1d0)) :: parvec(1)
    procedure(D2DFUN), pointer :: fPtrUni, fPtrPot
    fPtrUni => Unity
    fPtrPot => Potential 
    allocate(Smat(BS_ORDER,nEn))
    Smat=0.d0
    allocate(Hmat,source=Smat)
    parvec(1)=0.d0
    do iBs2 = 2, nEn + 1
       do iBs1 = max(2, iBs2-BS_ORDER+1), iBs2
          Overlap         = BSpline%Integral(fPtrUni,iBs1,iBs2)
          KineticEnergy   = BSpline%Integral(fPtrUni,iBs1,iBs2,1,1) / 2.d0
          PotentialEnergy = Bspline%Integral(fPtrPot,iBs1,iBs2,parvec=parvec) 
          Smat(iBs1+BS_ORDER-iBs2,iBs2-1) = Overlap      
          Hmat(iBs1+BS_ORDER-iBs2,iBs2-1) = KineticEnergy + PotentialEnergy
       enddo
    enddo
  end subroutine FillMatrices


  subroutine GeneralizedEigenproblem()
    real(kind(1d0)), allocatable :: Eval(:), Evec(:,:), Work(:)
    real(kind(1d0)), parameter   :: ERROR_THRESHOLD = 1.d-10 
    real(kind(1d0)), parameter   :: PI = acos(-1.d0) 
    integer        , parameter   :: npts=301
    real(kind(1d0))              :: EigenvalueError, x, fun
    integer                      :: iEn, ix, info, uid
    character(len=10)            :: iEnStrn

    !.. Solve the generalized eigenvalue problem for banded matrices
    allocate(work(3*nEn),Eval(nEn),Evec(nEn+2,nEn))
    work=0.d0
    Eval=0.d0
    Evec=0.d0
    call DSBGV( 'V', 'U', nEn, BS_ORDER-1, BS_ORDER-1, &
         Hmat, BS_ORDER, Smat, BS_ORDER, Eval, Evec(2,1), nEn+2, work, info )
    if( info /= 0 )then
       write(ERROR_UNIT,"(a,i0)") "DSBGV info : ", info 
       error stop
    endif
    deallocate(work)

    do iEn = 1, nEn

       !.. Print the eigenvalue on screen
       !  EigenvalueError = Eval(iEn) - (iEn*PI)**2/2.d0 
      !  EigenvalueError = Eval(iEn) + 1.d0 / ( 2.d0 * dble(iEn)**2 )             !.. l = 0
       EigenvalueError = Eval(iEn) + 1.d0 / ( 2.d0 * dble(1+iEn)**2 )           !.. l = 1
       if( EigenvalueError > ERROR_THRESHOLD )exit
       write(OUTPUT_UNIT,"(i4,*(x,e24.16))") iEn, Eval(iEn), EigenvalueError

       !.. Print the eigenfunction on file
       write(iEnStrn,"(i000.3)") iEn
       open(newunit = uid, &
            file    ="EigenState_"//trim(adjustl(iEnStrn)), &
            form    ="formatted", &
            status  ="unknown")
       do ix = 1, npts
          x = BS_GRMIN + (BS_GRMAX - BS_GRMIN)/(npts-1)*(ix-1)
          fun = BSpline%Eval(x,Evec(:,iEn))
          write(uid,"(*(x,e24.16))") x, fun
       enddo
       close(uid)

    enddo
    write(OUTPUT_UNIT,"(a,i0)") "Number of Accurate Eigenvalues : ",iEn-1
  end subroutine GeneralizedEigenproblem

  Pure real(kind(1d0)) function Unity(x,parvec) result(y)
    DoublePrecision, intent(in) :: x
    DoublePrecision, optional, intent(in) :: parvec(*)
    y=1.d0
  end function Unity

  Pure real(kind(1d0)) function Potential(x,parvec) result(res)
    real(kind(1d0)), intent(in) :: x
    real(kind(1d0)), optional, intent(in) :: parvec(*)
    ! res = -1.d0 / dble(x)                                                    !.. l = 0
    res = 1.d0 / dble(x)**2 - 1.d0 / dble(x)                                 !.. l = 1
  end function Potential
  
end program main

