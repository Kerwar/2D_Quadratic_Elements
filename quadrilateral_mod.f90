module quadrilateral_mod

  use number_types_mod, only: DP 
  use point_mod, only: Point

  implicit none
  
  private
  
  public :: basisF
  
  ! The reference square is the square (-1,-1), (1,-1), (1,1), (-1,1)
  ! We are going to assume Gaussian quadrature of order 2
  real(DP), parameter :: x1 = - 1.0_DP / sqrt(3.0), x2 = 1.0_DP / sqrt(3.0)
  real(DP), parameter :: w1 = 1.0_DP, w2 = 1.0_DP

  abstract interface 
    real(DP) function funct(p) result(res)
      import :: Point, DP
      type(Point) :: p
    end function funct

  end interface

  type, public :: Square
    type(Point) :: a(4)
    real(DP) :: jacobDetInGaussPoint(2,2)
    contains

    procedure :: integral
    procedure :: refToMain 
    procedure :: jacobDet
    ! procedure :: area => area_quadrilateral
    ! procedure :: basis
  end type

  interface Square
    procedure :: newSquare
  end interface

  contains

  real(DP) function integral(self, f) result(res)
    class(square) :: self
    procedure(funct) :: f
    res = f(Point(x1,x1)) * w1 * w1 * self%jacobDetInGaussPoint(1,1) + &
          f(Point(x1,x2)) * w1 * w2 * self%jacobDetInGaussPoint(1,2) + &
          f(Point(x2,x1)) * w2 * w1 * self%jacobDetInGaussPoint(2,1) + &
          f(Point(x2,x2)) * w2 * w2 * self%jacobDetInGaussPoint(2,2) 
    print *, f(Point(x1,x1)), self%jacobDetInGaussPoint(1,1), f(Point(x1,x2)), self%jacobDetInGaussPoint(1,2), &
    f(Point(x2,x1)), self%jacobDetInGaussPoint(2,1), f(Point(x2,x2)), self%jacobDetInGaussPoint(2,2)
  end function integral

  type(Square) function newSquare(a1, a2, a3, a4)
    type(Point) :: a1, a2, a3, a4

    newSquare%a(1) = a1 
    newSquare%a(2) = a2 
    newSquare%a(3) = a3 
    newSquare%a(4) = a4 

    newSquare%jacobDetInGaussPoint(1,1) = newSquare%jacobDet(Point(x1,x1))
    newSquare%jacobDetInGaussPoint(1,2) = newSquare%jacobDet(Point(x1,x2))
    newSquare%jacobDetInGaussPoint(2,1) = newSquare%jacobDet(Point(x2,x1))
    newSquare%jacobDetInGaussPoint(2,2) = newSquare%jacobDet(Point(x2,x2))

  end function newSquare

  type(Point) function refToMain(self, p) result(res)
    class(Square) :: self
    type(Point)   :: p
    real(DP) :: x, y 
    integer  :: i 

    x = sum([(basisF(i,p) * self%a(i)%x, i = 1,4)])
    y = sum([(basisF(i,p) * self%a(i)%y, i = 1,4)])

    res = Point(x,y)

  end function refToMain

  real(DP) function jacobDet(self, p) result(res)
    class(Square) :: self
    type(Point)   :: p
    real(DP)      :: dxdpsi1, dxdpsi2, dydpsi1, dydpsi2
    integer       :: i

    dxdpsi1 = sum([(DX_basisF(i,p) * self%a(i)%x, i = 1,4)])
    !print *, p%x, p%y, [(self%a(i)%x, i = 1,4)]
    dxdpsi2 = sum([(DY_basisF(i,p) * self%a(i)%x, i = 1,4)])

    dydpsi1 = sum([(DX_basisF(i,p) * self%a(i)%y, i = 1,4)])
    dydpsi2 = sum([(DY_basisF(i,p) * self%a(i)%y, i = 1,4)])
    res = abs(dxdpsi1 * dydpsi2 - dxdpsi2 * dydpsi1)
    print *, p%x, p%y, dxdpsi1, dxdpsi2, dydpsi1, dydpsi2, res
  end function jacobDet

  real(DP) function basisF(A, p) result(res)
    integer :: A
    type(Point) :: p

    select case(A)
      case(1)
        res = linear_1D(1, p%x) * linear_1D(1, p%y)
      case(2)
        res = linear_1D(2, p%x) * linear_1D(1, p%y)
      case(3)
        res = linear_1D(2, p%x) * linear_1D(2, p%y)
      case(4)
        res = linear_1D(1, p%x) * linear_1D(2, p%y)
      case DEFAULT
        print *, "basisF basis can't have index", A, " in 2D"
        stop
    end select

  end function basisF

  real(DP) function DX_basisF(A, p) result(res)
    integer :: A
    type(Point) :: p

    select case(A)
      case(1)
        res = D_linear_1D(1, p%x) * linear_1D(1, p%y)
      case(2)
        res = D_linear_1D(2, p%x) * linear_1D(1, p%y)
      case(3)
        res = D_linear_1D(2, p%x) * linear_1D(2, p%y)
      case(4)
        res = D_linear_1D(1, p%x) * linear_1D(2, p%y)
      case DEFAULT
        print *, "basisF basis can't have index", A, " in 2D"
        stop
    end select

  end function DX_basisF

  real(DP) function DY_basisF(A, p) result(res)
    integer :: A
    type(Point) :: p

    select case(A)
      case(1)
        res = linear_1D(1, p%x) * D_linear_1D(1, p%y)
      case(2)
        res = linear_1D(2, p%x) * D_linear_1D(1, p%y)
      case(3)
        res = linear_1D(2, p%x) * D_linear_1D(2, p%y)
      case(4)
        res = linear_1D(1, p%x) * D_linear_1D(2, p%y)
      case DEFAULT
        print *, "basisF basis can't have index", A, " in 2D"
        stop
    end select

  end function DY_basisF

  real(DP) function linear_1D(A, p) result(res)
    integer  :: A
    real(DP) :: p
    
    if (A == 1) then
      res = (1.0_DP - p) / 2.0_DP
    else if (A == 2) then
      res = (1.0_DP + p) / 2.0_DP
    else 
      print *, "Linear Basis can't have index", A 
      stop
    end if
  end function linear_1D

  real(DP) function D_linear_1D(A, p) result(res)
    integer  :: A
    real(DP) :: p
    
    if (A == 1) then
      res = - 1.0_DP / 2.0_DP
    else if (A == 2) then
      res =   1.0_DP / 2.0_DP
    else 
      print *, "Linear Basis Derivative can't have index", A 
      stop
    end if
  end function D_linear_1D

  real(DP) function gauss_quad(f) result(res)

    procedure(funct) :: f
    res = f(Point(x1,x1)) * w1 * w1 + f(Point(x1,x2)) * w1 * w2 + &
    f(Point(x2,x1)) * w2 * w1 + f(Point(x2,x2)) * w2 * w2 
  end function

end module quadrilateral_mod