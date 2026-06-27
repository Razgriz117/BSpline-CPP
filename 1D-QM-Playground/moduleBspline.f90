!! ModuleBSplines, Copyright (C) 2020 Luca Argenti, PhD - Some Rights Reserved
!! ModuleBSplines is licensed under a
!! Creative Commons Attribution-ShareAlike 4.0 International License.
!! A copy of the license is available at <http://creativecommons.org/licenses/by-nd/4.0/>.
!!
!! Luca Argenti is Associate Professor of Physics, Optics and Photonics
!! Department of Physics and the College of Optics
!! University of Central Florida
!! 4111 Libra Drive, Orlando, Florida 32816, USA
!! email: luca.argenti@ucf.edu

module ModuleBSpline

  use, intrinsic :: ISO_FORTRAN_ENV

  implicit none

  private

  interface 
     Pure DoublePrecision function D2DFun( x, parvec ) 
       DoublePrecision          , intent(in) :: x
       DoublePrecision, optional, intent(in) :: parvec(*)
     end function D2DFun
  end interface
  public D2DFun

  ! {{{ Private Attributes

  !
  logical                    :: INITIALIZED = .FALSE.
  !
  integer        , parameter :: IOMSG_LENGTH= 100
  !.. Local Gaussian weights
  DoublePrecision, parameter :: PI = 3.14159265358979323846d0
  DoublePrecision, parameter :: GAUSS_POINTS_CONVERGENCE_THRESHOLD = 2.d-15
  integer        , parameter :: NGAUSS = 64
  DoublePrecision            :: Gauss_points(NGAUSS)
  DoublePrecision            :: Gauss_weight(NGAUSS)
  !
  !.. Local Factorial Parameters and variables
  integer, parameter :: DIMFACT = 127
  DoublePrecision    :: fact(0:DIMFACT)
  DoublePrecision    :: tar(0:((DIMFACT+1)*(DIMFACT+2)/2-1))
  DoublePrecision    :: parfactMat(0:DIMFACT,0:DIMFACT)
  !
  !.. B-spline parameters
  integer        , parameter :: MAX_NUMBER_NODES = 10000
  DoublePrecision, parameter :: NOD_THRESHOLD    = 1.d-15
  integer        , parameter :: MAXNR            = 20
  !

  ! }}}

  !> Set of B-splines (spline basis functions)
  !!
  ! {{{ Detailed Description

  !> The spline are defined as piecewise polynomials of order \f$k\f$ 
  !! (maximum degree \f$k-1\f$), \f$\mathcal{C}^\infty\f$ everywhere, except at a fixed set of knots \f$\{t_i\}\f$ 
  !! where they are at least \f$\mathcal{C}^{k_i-2}\f$, with \f$k \geq k_i\geq 1\f$. 
  !! \f$\nu_i=k-k_i+1\f$ is called the knot multiplicity because the same spline basis can be obtained 
  !! from a basis with \f$\sum_i \nu_i\f$ maximally regular knots in the limit where sets of \f$\nu_i\f$ 
  !! contiguous knots coalesce to \f$i\f$-th knots. In other terms, the spline space is specified by an order 
  !! \f$k\f$ and a non-decreasing set of knots. The total dimension of spline space is \f$n+k\f$, but if 
  !! we specialize to the subset which is zero below the lowest and above the highest knots we are left 
  !! with a space of dimension \f$n-k\f$. As a consequence every spline extends at least over \f$k\f$ adjacent 
  !! intervals (\f$k+1\f$ consecutive knots, when counted with their multiplicity). The \f$n-k\f$ independent 
  !! splines restricted to \f$k\f$ intervals, clearly a basis of the spline space, are commonly called B-splines.
  !! Following deBoor\cite{deBoor}, we define the \f$i\f$-th B-spline \f$B_i^k(x)\f$ of order \f$k\f$ and which 
  !! extends from the knot \f$t_i\f$ to the knot \f$t_{i+k}\f$, as follows:
  !! \f{eqnarray}
  !! B_i^1(x)&=&\theta(x-t_i)\,\cdot\,\theta(t_{i+1}-x)\\
  !! B_i^k(x)&=&\frac{x-t_i}{t_{i+k-1}-t_i}\,\,
  !!            B_i^{k-1}(x)+\frac{t_{i+k}-x}{t_{i+k}-t_{i+1}}\,\,
  !!            B_{i+1}^{k-1}(x).
  !! \f}
  !! In the following, unless otherwise stated, we shall refer to standard set of knots where the 
  !! first and last breakpoints are \f$k\f$ times degenerate, while the other nodes are non degenerate:
  !! \f[
  !! t_1=t_2=\ldots=t_k\leq t_{k+1}\leq\ldots\leq t_{n-k}\leq t_{n-k+1}=t_{n-k+2}=\ldots=t_{n}
  !! \f]
  !! The use of B-splines in atomic and molecular physics calculations is reviewed in 
  !! [\cite{rpp.64.1815}](http://iopscience.iop.org/0034-4885/64/12/205/). 
  !! 
  !! B-splines are invariant upon affine transformations: that is, if
  !! \f{equation}
  !! \mathrm{if}\quad x\to x'=a\,x +b,\quad t_i\to t_i'=a\,t_i+b\quad\mathrm{then}
  !! \quad{B_i^{k}}'(x')=B_i^k(x).
  !! \f}
  !! 
  !! It is useful to define the \f$L^2\f$-normalized B-splines as
  !! \f{equation}
  !! \bar{B}_i(x)\equiv B_i(x)/\|B_i\|_{L^2}
  !! \f}
  !! 
  ! }}} 
  type, public :: ClassBSpline
     !
     private
     !> Number of non-degenerate nodes.
     integer :: NNodes
     !
     !> Spline order.
     !! A B-spline of order \f$ k\f$ is a piecewise polynomials
     !! of order \f$k-1\f$ which extends over \f$k\f$ intervals (including
     !! degenerate intervals, i.e., those with zero length).
     !! The quality of the Bspline evaluation is severely compromised
     !! for high order. At order around 35, in the case of just two nodes
     !! (which is the most demanding condition), the numerical error starts
     !! being of the same order of magnitude of the function itself.
     !! For order \f$\sim\f$20, there are occasional points of irregularity
     !! for a reasonable number of intervals (e.g.,\f$\geq 10\f$). Anyhow,
     !! these caveats arise only in cases which already lie outside
     !! the range of applicability of the BSplines as a variational basis, 
     !! which is above order 12 or so.
     integer :: Order
     !
     !> Total number of B-spline. 
     !! \f[
     !! NBsplines = NNodes + Order - 2
     !!\f]
     integer :: NBSplines
     !
     !> grid which contains the node positions. 
     !! - \f$g_{-k+1}\,=\,g_{-k+2}\,=\,\cdots\,=\,g_{-1}\,=\,g_0\f$ 
     !!   is the first \f$ k\f$ times degenerate node.
     !! - \f$g_{n-1}\,=\,g_{n}\,=\,\cdots\,=\,g_{n+k-2}\f$ is the 
     !!   last \f$ k\f$ times degenerate node.
     !! The set of distinct nodes is mapped to the positions [0:n-1] of grid;  
     !! and hence, the upper bound of the i-th BSpline support is grid(i)
     real(kind(1d0)),  allocatable :: grid(:)
     !
     !> Normalization constants.
     !! \f[
     !!    \mathrm{f}_i=\|B_i\|_{L^2}^{-1}=\frac{1}{ \sqrt{ \langle B_i | B_i \rangle } }
     !! \f]
     real(kind(1d0)), private, allocatable :: f(:)
     !
     !
     !> Matrix of polynomial coefficients.
     !! The support $\mathcal{D}_i$ of the \f$i\f$-th B-spline comprises \f$k\f$ consecutive intervals:
     !! \f[
     !!      \mathcal{D}_i = [g_{i-k},g_i] = [g_{i-k+1},g_{i-k+1}] \cup \cdots \cup [g_{i-1}:g_i].
     !! \f]
     !! If we enumerate the individual intervals from \f$0\f$ to \f$k-1\f$, the explicit expression of 
     !! the \f$i\f$-th B-spline in interval $j$ is
     !! \f[
     !!     B_i(x) = \sum_{m=0}^{k-1}
     !!        \left(
     !!                \frac{x-g_{i-k+j}}{g_{i-k+j+1}-g_{i-k+j}}
     !!        \right)^m \,\,  c(m,j,i)
     !! \f]
     real(kind(1d0)), private, allocatable :: c(:,:,:)
     !
     !> Indicates whether ClassBSpline has been initialized or not.
     logical, private :: INITIALIZED=.FALSE.
     !
   contains

     !.. Basics
     generic   :: Init =>  BSplineInit
     procedure :: Free =>  BSplineFree
     procedure :: Save =>  ClassBSpline_Save
     procedure :: Load =>  ClassBSpline_Load

     !.. Accessors
     !> Gets the number of nodes.
     procedure :: GetNNodes         =>  BSplineGetNNodes
     !> Gets the B-splines order.
     procedure :: GetOrder          =>  BSplineGetOrder
     !> Gets the number of B-splines.
     procedure :: GetNBSplines      =>  BSplineGetNBSplines
     !> Gets the position of a requested node.
     procedure :: GetNodePosition   =>  BSplineNodePosition
     !> Gets the normalization factor of a requested B-spline.
     procedure :: GetNormFactor     =>  BSplineNormalizationFactor

     !> Evaluates either a choosen B-spline function, 
     !! or a linear combination of B-splines, both evaluated in certain position. 
     generic   :: Eval              =>  BSplineEval, BSplineFunctionEval

     !> Tabulates in a slected domain either  a choosen B-spline function, 
     !! or a linear combination of B-splines.
     generic   :: Tabulate          =>  BSplineTabulate, BSplineFunctionTabulate

     !> Computes the integral
     !! \f[
     !!    \int_{a}^{b} r^{2}dr
     !!       \frac{d^{n_{1}}Bs_{i}(r)}{dr^{n_{1}}} f(r)
     !!       \frac{d^{n_{2}}Bs_{j}(r)}{dr^{n_{2}}}
     !!\f]
     !! Where \f$f(r)\f$ is a local operator, and if a break point \f$BP\f$ is introduced, 
     !! then the integral is splitted in two parts
     !! \f[
     !!    \int_{a}^{BP} r^{2}dr
     !!       \frac{d^{n_{1}}Bs_{i}(r)}{dr^{n_{1}}} f(r)
     !!       \frac{d^{n_{2}}Bs_{j}(r)}{dr^{n_{2}}} + 
     !!    \int_{BP}^{b} r^{2}dr
     !!       \frac{d^{n_{1}}Bs_{i}(r)}{dr^{n_{1}}} f(r)
     !!       \frac{d^{n_{2}}Bs_{j}(r)}{dr^{n_{2}}}
     !! \f]
     generic   :: Integral =>  BSplineIntegral
     procedure :: BSplineIntegral

     !.. Assignment
     !> Copies the B-spline set from one ClassBSpline to another. 
     generic, public :: ASSIGNMENT(=) => CopyBSplineSet

     !.. Private Interface
     procedure, private :: BSplineInit
     procedure, private :: CopyBSplineSet
     procedure, private :: BSplineEval
     procedure, private :: BSplineFunctionEval
     procedure, private :: BSplineTabulate
     procedure, private :: BSplineFunctionTabulate

  end type ClassBSpline

contains

  subroutine ClassBspline_Save( self, uid )
    class(ClassBspline), intent(in) :: self
    integer            , intent(in)  :: uid
    integer :: i
    write(uid,"(*(x,i0))") &
         self%Nnodes, self%Order, self%NBsplines, &
         size(self%grid), size(self%f), &
         (size(self%c,i),i=1,3) 
    write(uid,"(*(x,e24.16))") self%grid
    write(uid,"(*(x,e24.16))") self%f
    write(uid,"(*(x,e24.16))") self%c
  end subroutine ClassBspline_Save
   
  subroutine ClassBspline_Load( self, uid )
    class(ClassBspline), intent(out) :: self
    integer            , intent(in)  :: uid
    integer :: ng, nf, nc1, nc2, nc3
    call self%free()
    read(uid,*) self%Nnodes, self%Order, self%NBsplines, &
         ng, nf, nc1, nc2, nc3
    allocate(self%grid(ng))
    allocate(self%f   (nf))
    allocate(self%c(nc1,nc2,nc3))
    read(uid,*) self%grid
    read(uid,*) self%f
    read(uid,*) self%c
  end subroutine ClassBspline_Load


  !> Initialize module basic constants
  subroutine Init_Spline_Module()
    call InitGaussPoints
    call InitFactorials
    INITIALIZED=.TRUE.
  end subroutine Init_Spline_Module


  subroutine CheckInitialization(s)
    Class(ClassBSpline), intent(in) :: s
    if(.not.s%INITIALIZED) error stop
  end subroutine CheckInitialization


  !> Initialize the points and weights for 
  !! Gauss integration
  subroutine InitGaussPoints()
    integer         :: i,j
    DoublePrecision :: p1, p2, p3, pp, z, z1
    do i=1,floor((NGAUSS+1)/2.d0)
       z=cos(PI*(i-.25d0)/(NGAUSS+.5d0))
       inna : do
          p1=1.d0
          p2=0.d0
          do j=1,NGAUSS
             p3=p2
             p2=p1
             p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
          enddo
          pp=NGAUSS*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
          if(abs(z-z1)<=GAUSS_POINTS_CONVERGENCE_THRESHOLD)exit inna
       enddo inna
       Gauss_Points(i)=(1.d0-z)/2.d0
       Gauss_Points(NGAUSS+1-i)=(1.d0+z)/2.d0
       Gauss_Weight(i)=1.d0/((1.d0-z*z)*pp*pp)
       Gauss_Weight(NGAUSS+1-i)=1.d0/((1.d0-z*z)*pp*pp)
    enddo
  end subroutine InitGaussPoints


  subroutine InitFactorials()
    integer :: i,j,k,n
    !
    !Initialize Factorials
    fact(0)=1.d0
    do i = 1, DIMFACT
       fact(i)=fact(i-1)*dble(i)
    enddo
    !
    !.. Initialize Binomials
    k=0
    tar=0.d0
    do i = 0,DIMFACT
       do j = 0,i
          tar(k)=fact(i)/fact(i-j)/fact(j)
          k=k+1
       enddo
    enddo
    !
    parfactmat=1.d0
    do n=1,DIMFACT
       do k=1,n
          do i=n-k+1,n
             parfactmat(n,k)=parfactmat(n,k)*dble(i)
          enddo
       enddo
    enddo
    !
  end subroutine InitFactorials

  integer function BSplineGetOrder(s) result(Order)
    class(ClassBSpline), intent(in) :: s
    Order=s%Order
  end function BSplineGetOrder

  integer function BSplineGetNNodes(s) result(NNodes)
    class(ClassBSpline), intent(in) :: s
    NNodes=s%NNodes
  end function BSplineGetNNodes

  integer function BSplineGetNBSplines(s) result(NBSplines)
    class(ClassBSpline), intent(in) :: s
    NBSplines=s%NBSplines
  end function BSplineGetNBSplines

  !> Return the position of a node
  DoublePrecision function BSplineNodePosition(SplineSet,Node) result(Position)
    Class(ClassBSpline), intent(in) :: SplineSet
    integer           , intent(in) :: Node
    if( Node <= 1 )then
       Position = SplineSet%Grid(0)
    elseif( Node >= SplineSet%NNodes )then
       Position = SplineSet%Grid( SplineSet%NNodes-1 )
    else
       Position = SplineSet%Grid( Node-1 )
    endif
    return
  end function BSplineNodePosition

  !> Return the position of a node
  DoublePrecision function BSplineNormalizationFactor(SplineSet,Bs) &
       result(Normalization)
    Class(ClassBSpline), intent(in) :: SplineSet
    integer            , intent(in) :: Bs
    Normalization = 0.d0
    if(Bs<1.or.Bs>SplineSet%NBSplines)return
    Normalization=SplineSet%f(Bs)
  end function BSplineNormalizationFactor


  !> Copies the B-spline set sOrigin to sDestination 
  subroutine CopyBSplineSet(sDestination,sOrigin)
    Class(ClassBSpline), intent(inout):: sDestination
    Class(ClassBSpline), intent(in)   :: sOrigin
    !
    sDestination%NNodes=sOrigin%NNodes
    sDestination%Order=sOrigin%Order
    sDestination%NBSplines=sOrigin%NBSplines
    !
    if(allocated(sDestination%Grid))deallocate(sDestination%Grid)
    allocate(sDestination%Grid(1-sDestination%Order:sDestination%NBSplines))
    sDestination%Grid=sOrigin%Grid
    !
    if(allocated(sDestination%c))deallocate(sDestination%c)
    allocate(sDestination%c(&
         0:sDestination%Order-1 ,&
         0:sDestination%Order-1 ,&
         1:sDestination%NBSplines+sDestination%Order-2))
    sDestination%c=sOrigin%c
    !
    if(allocated(sDestination%f))deallocate(sDestination%f)
    allocate(sDestination%f(sDestination%NBSplines))
    sDestination%f=sOrigin%f
    !
    sDestination%INITIALIZED=.TRUE.
  end subroutine CopyBSplineSet


  subroutine CheckNodesAndOrder(NumberOfNodes,Order,STAT)
    integer, intent(in)  :: NumberOfNodes, Order
    integer, intent(out) :: STAT
    STAT = -1; if ( NumberOfNodes < 2 ) return
    STAT = -2; if (    order      < 1 ) return
    STAT = 0
  end subroutine CheckNodesAndOrder
  

  !> Initialize a BSpline set
  Subroutine BSplineInit(s,NumberOfNodes,order,grid,IOSTAT)
    !
    !> Bspline basis set to be initialized
    Class(ClassBSpline), intent(inout)           :: s
    !> Number of Nodes
    integer            , intent(in)              :: NumberOfNodes
    !> Spline order
    integer            , intent(in)              :: order
    !> Grid of nodes (without repetition of degenerate nodes)
    DoublePrecision    , intent(in)              :: grid(1:NumberOfNodes)
    !> IOSTAT = 0 on successful exit, IOSTAT/=0 otherwise
    !! If IOSTAT is not specified, the programs terminates with
    !! an error.
    integer            , intent(out),   optional :: IOSTAT

    !Local variables
    integer                    :: i, Bs, Status
    procedure(D2DFUN), pointer :: FunPtr

    if(present(IOSTAT)) IOSTAT=0
    if(.not.INITIALIZED)call Init_Spline_Module()

    !.. Check input
    call CheckNodesAndOrder(NumberOfNodes,Order,Status)
    if(Status/=0)then
       if(present(IOSTAT))then
          IOSTAT=Status
          return
       endif
       error stop 
    endif

    s%NNodes    = NumberOfNodes
    s%Order     = order
    s%NBSplines = NumberOfNodes + order - 2

    !.. Check Grid format
    if(  LBOUND(Grid,1) > 1             .or. &
         UBOUND(Grid,1) < NumberOfNodes ) then
       if(present(IOSTAT))then
          IOSTAT=-4
          return
       endif
       error stop
    endif
    do i=2,NumberOfNodes
       if( Grid(i) < Grid(i-1) )then
          if(present(IOSTAT))then
             IOSTAT=1
             return
          else
             error stop
          endif
       endif
    enddo

    !.. Maps the Grid to [0:n-1]; hence, the upper bound
    !   of the i-th BSpline support is g(i)
    if( allocated( s%Grid ) ) deallocate( s%Grid )
    allocate( s%Grid( -(order-1) : NumberOfNodes-1 + order-1 ) )
    s%Grid( 1-order : -1 ) = Grid(1)
    s%Grid( 0:NumberOfNodes-1 ) = Grid(1:NumberOfNodes)
    s%Grid( NumberOfNodes:NumberOfNodes+order-2 ) = Grid(NumberOfNodes)

    !.. Initializes the BSpline coefficients normalized
    !   so that the B-spline set is a partition of unity
    !..
    if( allocated( s%c ) ) deallocate( s%c )
    allocate( s%c( 0:order-1, 0:order-1, 1:s%NBSplines ) )
    s%c = 0.d0
    do Bs=1,s%NBSplines
       call ComputeCoefficientsSingleBSpline(s%Order,s%Grid(Bs-s%Order),s%c(0,0,Bs))
    enddo
    s%INITIALIZED=.TRUE.

    !.. Initialize the normalization factors
    !   ||Bs_i(x)*s%f(i)||_L2=1.d0
    !..
    if( allocated( s%f ) ) deallocate( s%f )
    allocate( s%f( 1 : s%NBSplines ) )
    s%f = 0.d0
    funPtr=>Unity
    do i = 1, s%NBSplines
       s%f(i)=1.d0/sqrt(s%Integral(FunPtr,i,i))
    enddo

  end Subroutine BSplineInit


  Pure DoublePrecision function Unity(x,parvec) result(y)
    DoublePrecision, intent(in) :: x
    DoublePrecision, optional, intent(in) :: parvec(*)
    y=1.d0
  end function Unity


  !> Given the order and a set of order+1 consecutive nodes,
  !> builds the corresponding (uniquely defined) B-spline.
  subroutine ComputeCoefficientsSingleBSpline(order,vec,Coefficients)
    !
    !> Order of the B-spline = number of (possibly degenerate) 
    !> intervals that form the B-spline support = Polynomials degree + 1
    !> (e.g., order=4 means piecewise cubic polynomials, which require
    !> four parameters: \f$ a_0 + a_1 x + a_2 * x^2 + a_3 * x^3 \f$.
    integer        , intent(in)  :: order
    !
    !> vector of order+1 nodes. For Order = 5, for example, vec contains 
    !> the six bounds of the five consecutive intervals that form the 
    !> support of the B-spline, indexed from ONE to SIX.
    !
    !  vec(1)  vec(2)    vec(3)  vec(4) vec(5)    vec(6)
    !    |-------|---------|-------|------|---------|
    DoublePrecision, intent(in)  :: vec(1:order+1)
    !
    !> Matrix of polynomial coefficients in each interval.
    DoublePrecision, intent(out) :: Coefficients(0:order-1,0:order-1)
    !
    !Local variables
    integer         :: k, inter, tempBs
    DoublePrecision :: SpanBSpline
    DoublePrecision :: SpanCurrentInterval
    DoublePrecision :: SpanPriorToInterval
    DoublePrecision :: SpanToTheEnd
    DoublePrecision :: FactorConstant
    DoublePrecision :: FactorMonomial
    DoublePrecision, allocatable :: coe(:,:,:,:)
    !
    !.. Computes the B-spline in each available
    !   sequence of order consecutive intervals
    !
    !
    allocate( coe(0:order-1,0:order-1,order,order) )
    coe=0.d0
    coe(0,0,1,:)=1.d0
    !
    !.. Cycle over the spline order
    !
    do k = 2, order
       !
       !.. Cycle over the (temporary) Bsplines available 
       !   for any given order, which must all be computed
       !
       do tempBs = 1, order - k + 1
          !
          !.. Each B-spline of order k is obtained as the sum of 
          !   the contribution from the two B-splines of order k-1
          !   whose support is comprised in that of the final 
          !   B-spline.
          !
          !.. Contribution of the first B-spline of order k-1
          !   counted only if the B-spline is not degenerate
          !..
          SpanBSpline = vec(tempBs+k-1) - vec(tempBs)
          if( SpanBSpline > NOD_THRESHOLD )then
             do inter=0,k-2
                !
                !.. Constant Component
                SpanPriorToInterval = vec(tempBs+inter) - vec(tempBs)
                FactorConstant = SpanPriorToInterval / SpanBSpline 
                coe(0:k-2,inter,k,tempBs) = coe(0:k-2,inter,k,tempBs) + &
                     coe(0:k-2,inter,k-1,tempBs) * FactorConstant
                !
                !.. Monomial Component
                SpanCurrentInterval = vec(tempBs+inter+1) - vec(tempBs+inter)
                FactorMonomial = SpanCurrentInterval / SpanBSpline
                coe(1:k-1,inter,k,tempBs) = coe(1:k-1,inter,k,tempBs) + &
                     coe(0:k-2,inter,k-1,tempBs) * FactorMonomial
                !
             enddo
          endif
          !
          !.. Contribution of the second B-spline of order k-1
          !   counted only if the B-spline is not degenerate
          !..
          SpanBSpline = vec(tempBs+k) - vec(tempBs+1)
          if( SpanBSpline > NOD_THRESHOLD )then
             do inter=1,k-1
                !
                !.. Constant Component
                SpanToTheEnd = vec(tempBs+k)-vec(tempBs+inter)
                FactorConstant = SpanToTheEnd / SpanBSpline
                coe(0:k-2,inter,k,tempBs) = coe(0:k-2,inter,k,tempBs) + &
                     coe(0:k-2,inter-1,k-1,tempBs+1) * FactorConstant
                !
                !.. Monomial Component
                SpanCurrentInterval = vec(tempBs+inter+1) - vec(tempBs+inter)
                FactorMonomial = SpanCurrentInterval / SpanBSpline
                coe(1:k-1,inter,k,tempBs) = coe(1:k-1,inter,k,tempBs) - &
                     coe(0:k-2,inter-1,k-1,tempBs+1) * FactorMonomial
                !
             enddo
          endif
          !
       enddo
    enddo
    !
    Coefficients=coe(:,:,order,1)
    deallocate(coe)
    !
  end subroutine ComputeCoefficientsSingleBSpline


  !> Deallocates a ClassBSpline variable
  subroutine BSplineFree(s)
    Class(ClassBSpline), intent(inout) :: s
    if(allocated(s%Grid))deallocate(s%Grid)
    if(allocated(s%f))deallocate(s%f)
    if(allocated(s%c))deallocate(s%c)
    s%NNodes=0
    s%Order=0
    s%NBSplines=0
    s%INITIALIZED=.FALSE.
    return
  end subroutine BSplineFree


  !> Returns n if x is in \f$(g_n,g_{n+1}]\f$.
  integer function which_interval(x,s)
    !
    DoublePrecision   , intent(in) :: x
    type(ClassBSpline), intent(in) :: s
    !
    integer, save :: i1 = 0
    integer       :: i2
    !
    if(x>s%Grid(i1).and.x<=s%Grid(i1+1))then
       which_interval=i1
       return
    endif
    if(x>s%Grid(i1+1).and.x<=s%Grid(min(i1+2,s%NNodes+s%Order-1)))then
       which_interval=i1+1
       return
    endif
    i1=0
    i2=0
    !
    if(x<s%Grid(0))then
       which_interval=-1
       return
    elseif(x>s%Grid(s%NNodes-1))then
       which_interval=s%NNodes!-1
       return
    endif
    if(x<s%Grid(i1))i1=0
    if(x>s%Grid(i2))i2=s%NNodes-1
    do while(i2-i1>1)
       if(x>s%Grid((i2+i1)/2))then
          i1=(i1+i2)/2
       else
          i2=(i1+i2)/2
       endif
    enddo
    which_interval=i1
    return
  end function which_interval


  !> Computes the n-th derivative of B-spline Bs, \f$ n=0,1,2,\ldots\f$
  !!
  ! {{{ Detailed Description:

  !! ---------------------
  !! B-splines are positive: expansion coefficients of represented functions are 
  !! akin to the functions themselves. Without rapid oscillations of coefficients, 
  !! the cancellation errors are kept at minimum. 
  !! In particular \cite{deBoor},
  !! if \f$f=\sum_i B_ic_i\f$
  !! \f{equation}
  !! \left|c_i-\frac{M+m}{2}\right|\leq D_{k}\frac{M-m}{2},\quad m
  !!  =\min_{x\in[a,b]} f(x),\,\,M=\max_{x\in[a,b]}f(x)
  !! \f}
  !! where \f$D_k\f$ is a constant that depends only on \f$k\f$ and not on the
  !! particular partition of the \f$[a,b]\f$ interval.
  !! At any point, the evaluation of a linear combination of B-splines 
  !! requires the evaluation of \f$k\f$ basis functions only. 

  ! }}}
  DoublePrecision function BSplineEval(s,x,Bs,n_)
    Class(ClassBSpline), intent(in) :: s
    DoublePrecision    , intent(in) :: x
    integer            , intent(in) :: Bs
    integer, optional  , intent(in) :: n_
    !
    integer         :: i,j,k,n
    DoublePrecision :: r,a,w,OneOverInterval,OneOverIntervalToTheN
    !
    call CheckInitialization(s)
    !
    BSplineEval=0.d0
    if(Bs<1.or.Bs>s%NBSplines)return
    !
    n=0; if(present(n_)) n=n_
    if(n<0) error stop 
    if(n>=s%Order)return
    !
    i=which_interval(x,s)
    if(  i<max(0,Bs-s%Order).or.&
         i>min(s%NNodes-1,Bs-1))return
    !
    j=i-Bs+s%Order
    OneOverInterval=1.d0/(s%Grid(i+1)-s%Grid(i))
    r=(x-s%Grid(i))*OneOverInterval
    !
    a=1.d0
    w=0.d0
    do k=n,s%Order-1
       w=w+a*parfactmat(k,n)*s%c(k,j,Bs)
       a=a*r
    enddo
    !
    OneOverIntervalToTheN=1.d0
    do k=1,n
       OneOverIntervalToTheN=OneOverIntervalToTheN*OneOverInterval
    enddo
    BSplineEval=w*OneOverIntervalToTheN
    !
    return
  end function BSplineEval

  !> Computes the k-th derivative of a function expressed in terms of B-splines.
  DoublePrecision function BSplineFunctionEval(s,x,fc,k,SKIP_FIRST,Bsmin,Bsmax)
    class(ClassBSpline), intent(in) :: s
    DoublePrecision    , intent(in) :: x
    DoublePrecision    , intent(in) :: fc(:)
    integer, optional  , intent(in) :: k
    logical, optional  , intent(in) :: SKIP_FIRST
    integer, optional  , intent(in) :: Bsmin,Bsmax
    integer :: Bsmi, Bsma
    integer :: Bs,i,nskip
    call CheckInitialization(s)
    BSplineFunctionEval=0.d0
    i=which_interval(x,s)
    if(i<0.or.i>s%NNodes-2)return
    Bsmi=1
    nskip=0
    if(present(SKIP_FIRST))then
       if(SKIP_FIRST)then
          Bsmi=2
          nskip=1
       endif
    endif
    if(present(Bsmin))then
       Bsmi = max(Bsmin,Bsmi)
       nskip= max(Bsmi-1,nskip)
    endif
    Bsmi=max(Bsmi,i+1)
    Bsma = min(s%NBSplines,i+s%Order)
    if(present(Bsmax)) Bsma = min(Bsmax,Bsma)
    do Bs=Bsmi,Bsma
       BSplineFunctionEval=BSplineFunctionEval+s%Eval(x,Bs,k)*fc(Bs-nskip)
    enddo
    return
  end function BSplineFunctionEval
  
  !> Computes the  k-th derivatives of a choosen B-spline evaluated in a position array. 
  subroutine BSplineTabulate(s,ndata,xvec,yvec,Bs,n_) 
    Class(ClassBSpline), intent(in) :: s
    integer            , intent(in) :: ndata
    DoublePrecision    , intent(in) :: xvec(1:ndata)
    DoublePrecision    , intent(out):: yvec(1:ndata)
    integer            , intent(in) :: Bs
    integer, optional  , intent(in) :: n_
    !
    integer :: i,n
    n=0;if(present(n_))n=n_
    do i=1,ndata
       yvec(i)=s%Eval(xvec(i),Bs,n)
    enddo
  end subroutine BSplineTabulate

  !> Computes the  k-th derivatives of a function expressed in terms of B-splines and evaluated in an array of positions. 
  subroutine BSplineFunctionTabulate(s,ndata,xvec,yvec,FunVec,n_,SKIP_FIRST) 
    Class(ClassBSpline), intent(in) :: s
    integer            , intent(in) :: ndata
    DoublePrecision    , intent(in) :: xvec(1:ndata)
    DoublePrecision    , intent(out):: yvec(1:ndata)
    DoublePrecision    , intent(in) :: FunVec(:)
    integer, optional  , intent(in) :: n_
    logical, optional  , intent(in) :: SKIP_FIRST
    logical :: SKIP_FIRST_LOC
    !
    integer :: i,n
    SKIP_FIRST_LOC=.FALSE.
    if(present(SKIP_FIRST))then
       SKIP_FIRST_LOC=SKIP_FIRST
    endif
    n=0;if(present(n_))n=n_
    do i=1,ndata
       yvec(i)=s%Eval(xvec(i),FunVec,n,SKIP_FIRST=SKIP_FIRST_LOC)
    enddo
  end subroutine BSplineFunctionTabulate


  !> Computes the integral
  !! \f[
  !!    \int_{a}^{b} \frac{d^{n_{1}}Bs_{i}(r)}{dr^{n_{1}}} f(r)
  !!                 \frac{d^{n_{2}}Bs_{j}(r)}{dr^{n_{2}}} r^{2}dr
  !! \f]
  !! Where \f$f(r)\f$ is a local operator, and if a break point \f$BP\f$ 
  !! is introduced, then the integral is splitted in two parts
  !! \f[
  !!    \int_{a}^{BP} \frac{d^{n_{1}}Bs_{i}(r)}{dr^{n_{1}}} f(r)
  !!                  \frac{d^{n_{2}}Bs_{j}(r)}{dr^{n_{2}}} r^{2}dr + 
  !!    \int_{BP}^{b} \frac{d^{n_{1}}Bs_{i}(r)}{dr^{n_{1}}} f(r)
  !!                  \frac{d^{n_{2}}Bs_{j}(r)}{dr^{n_{2}}}r^{2}dr
  !! \f]
  DoublePrecision function BSplineIntegral( &
       s                 , &
       FunPtr            , & 
       Bs1               , &
       Bs2               , &
       BraDerivativeOrder, &
       KetDerivativeOrder, &
       LowerBound        , &
       UpperBound        , &
       BreakPoint        , &
       parvec            ) &
       result( Integral )
    !
    Class(ClassBSpline)      , intent(in) :: s
    procedure(D2DFun)        , pointer    :: FunPtr
    integer                  , intent(in) :: Bs1
    integer                  , intent(in) :: Bs2
    integer        , optional, intent(in) :: BraDerivativeOrder
    integer        , optional, intent(in) :: KetDerivativeOrder
    DoublePrecision, optional, intent(in) :: LowerBound
    DoublePrecision, optional, intent(in) :: UpperBound
    DoublePrecision, optional, intent(in) :: BreakPoint
    DoublePrecision, optional, intent(in) :: Parvec(*)
    !
    integer         :: n1, n2, iostat
    DoublePrecision :: a, b
    real(kind(1d0)) :: NewBreakPoint
    !
    Integral=0.d0
    call CheckInitialization(s)
    call CheckParameters( iostat )
    if( iostat/=0 )return
    if( present(BreakPoint) )then
       !***
       if ( BreakPoint < a ) then
          NewBreakPoint = a + epsilon(1d0)
       elseif (  BreakPoint > b ) then
          NewBreakPoint = b - epsilon(1d0)
       end if
       !***
       Integral = &
            BSplineDriverIntegral(s,FunPtr,Bs1,Bs2,n1,n2, a, BreakPoint,    Parvec ) + &
            BSplineDriverIntegral(s,FunPtr,Bs1,Bs2,n1,n2, NewBreakPoint, b, Parvec )
    else
       Integral = BSplineDriverIntegral(s,FunPtr,Bs1,Bs2,n1,n2,a,b, Parvec )
    endif
    !
    return
    !
  contains
    !
    subroutine CheckParameters( IOStat )
      integer, intent(out) :: IOStat
      call CheckBSplineIndexes( IOStat ); if( IOStat /= 0 ) return
      call CheckIntegralBounds( IOStat ); if( IOStat /= 0 ) return
      call CheckDerivativeOrder
    end subroutine CheckParameters
    !
    subroutine CheckBSplineIndexes( IOStat )
      integer, intent(out) :: IOStat
      IOStat=1
      if( min(Bs1,Bs2) <  1           ) return      
      if( max(Bs1,Bs2) >  s%NBSplines ) return
      if( abs(Bs1-Bs2) >= s%Order     ) return
      IOStat=0
    end subroutine CheckBSplineIndexes
    !
    subroutine CheckIntegralBounds( IOStat )
      integer, intent(out) :: IOStat
      IOStat=1
      a = s%Grid( max(Bs1,Bs2) - s%Order )
      b = s%Grid( min(Bs1,Bs2) )
      if(present(LowerBound))then
         if( LowerBound >= b )return
         a=max(a,LowerBound)
      endif
      if(present(UpperBound))then
         if( UpperBound <= a )return
         b=min(b,UpperBound)
      endif
      IOStat=0
    end subroutine CheckIntegralBounds
    !
    subroutine CheckDerivativeOrder
      n1=0; if(present(BraDerivativeOrder)) n1=BraDerivativeOrder
      n2=0; if(present(KetDerivativeOrder)) n2=KetDerivativeOrder
    end subroutine CheckDerivativeOrder
    !
  end function BSplineIntegral


  ! {{{ Detailed Description

  !> Compute the integral
  !! \f[
  !!    \int_a^b dr 
  !!    \frac{d^{n_1}B_1(r)}{dr^{n_1}} f(r)
  !!    \frac{d^{n_2}B_2(r)}{dr^{n_2}},
  !! \f]
  !! where \f$f(r)\f$ is a *smooth* function in the
  !! intersection of the supports of the two B-splines.
  !! If the first (last) boundary of the integration 
  !! interval is not specified, then the value a=-oo
  !! (b=+oo) is assumed.
  !!
  !!
  !! Properties of B-spline integrals:
  !! ---------------------------------
  !! Since the support \f$D_i\f$ of the \f$i\f$-th B-spline is formed by
  !! \f$k\f$ consecutive intervals, the integrals between two B-splines and a 
  !! local operator are zero unless their indices differ less than the order \f$k\f$:
  !! \f{equation}
  !! \langle B_i|O|B_j\rangle=\int_{D_i \cap D_j} B_i(x) o(x) B_j(x) dx.
  !! \f}
  !! As a consequence, matrices are sparse and often band-diagonal. 
  !! On uniform grids, matrices of translationally invariant operators, that 
  !! is with kernel \f$o(x,y)=o(x-y)\f$, are Toeplitz:
  !! \f{eqnarray}
  !! \langle B_i|O|B_j\rangle=\langle B_{i+n}|O|B_{j+n}\rangle.
  !! \f}
  !! If an hermitian operator \f$O\f$ is both local and translationally invariant, 
  !! like the identity and the kinetic energy, its matrix on a uniform grid
  !! is Toeplitz and banded, in other terms it is defined by just \f$k\f$ numbers. 
  !! 
  !! With Gauss-type integration technique, the matrix elements of operators 
  !! with polynomial kernels are exact. The error \f$\epsilon\f$ in the approximation
  !! of a \f$\mathcal{C}^k\f$ function \f$f\f$ with B-splines, is bounded by
  !! \f[
  !!     \epsilon\leq \mathrm{const}_k|\mathbf{t}|^k\|D^kf\|,\quad |\mathbf{t}|=
  !!         \max_i(t_{i+1}-t_i),\quad \|f\|=\max_{x\in[a,b]}|f(x)|
  !! \f]
  !!

  ! }}}
  DoublePrecision function BSplineDriverIntegral(s,FunPtr,Bs1,Bs2,n1,n2,a,b,parvec) &
       result( Integral )
    !
    !.. Assumes that the input data have been sanitized
    !
    Class(ClassBSpline), intent(in) :: s
    procedure(D2DFun)  , pointer    :: FunPtr
    integer            , intent(in) :: Bs1, Bs2
    integer            , intent(in) :: n1, n2
    DoublePrecision    , intent(in) :: a, b
    DoublePrecision, optional, intent(in) :: parvec(*)
    !
    integer         :: Interval, IntervalMin, IntervalMax
    DoublePrecision :: LowerBound, UpperBound
    !
    DoublePrecision :: PartialIntegral
    DoublePrecision :: IntervalWidth
    DoublePrecision :: Radius
    integer         :: iGauss
    !
    DoublePrecision :: aPlus, bMinus
    !
    Integral=0.d0
    !
    aPlus  = UpperLimitTo( a )
    bMinus = LowerLimitTo( b )
    IntervalMin = which_interval( aPlus,  s )
    IntervalMax = which_interval( bMinus, s )
    !
    do Interval = IntervalMin, IntervalMax
       LowerBound = max( a, s%Grid( Interval   ) )
       UpperBound = min( b, s%Grid( Interval+1 ) )
       !
       IntervalWidth = UpperBound - LowerBound
       if( IntervalWidth < NOD_THRESHOLD ) cycle
       !
       PartialIntegral = 0.d0
       do iGauss=1,NGAUSS
          !
          Radius = LowerBound + IntervalWidth * Gauss_Points( iGauss )
          PartialIntegral = PartialIntegral + &
               FunPtr(Radius, Parvec) * &
               s%Eval(Radius,Bs1,n1)  * &
               s%Eval(Radius,Bs2,n2)  * &
               Gauss_weight( iGauss )
          !
       enddo
       PartialIntegral = PartialIntegral * IntervalWidth
       !
       Integral = Integral + PartialIntegral
       !
    enddo
    !
  end function BSplineDriverIntegral

  Pure DoublePrecision function UpperLimitTo(x) result(xPlus)
    DoublePrecision, intent(in) :: x
    DoublePrecision :: tin,eps
    eps=epsilon(1.d0)
    tin=tiny(x)*(1.d0+eps)
    xPlus=(x+tin)*(1.d0+eps)
  end function UpperLimitTo

  Pure DoublePrecision function LowerLimitTo(x) result(xMinus)
    DoublePrecision, intent(in) :: x
    DoublePrecision :: tin,eps
    eps=epsilon(1.d0)
    tin=tiny(x)*(1.d0+eps)
    xMinus=(x-tin)*(1.d0-eps)
  end function LowerLimitTo

end module ModuleBSpline
