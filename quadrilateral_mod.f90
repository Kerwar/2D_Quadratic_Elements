module quadrilateral_mod

  use number_types_mod, only: DP 
  use point_mod, only: Point

  implicit none
  
  private
  
  public :: bilinear
  ! The reference square is the square (-1,-1), (1,-1), (1,1), (-1,1)
  ! We are going to assume Gaussian quadrature of order 2
  type, public :: Square
    type(Point) :: a(4)

    contains

    ! procedure :: ref => change_variable
    ! procedure :: area => area_quadrilateral
    ! procedure :: basis
  end type

  contains

! CALCULO DEL AREA DEL TRIANGULO MEDIANTE LA FORMULA DE HERON 

  ! real(DP) function area_quadrilateral(self) result(res)
  !   class(Square) :: self
  !   real(DP) :: edge1, edge2, edge3
  !   real(DP) :: s

  !   edge1 = sqrt((self%a1%x-self%a2%x) ** 2 + (self%a1%y-self%a2%y) ** 2)
  !   edge2 = sqrt((self%a1%x-self%a3%x) ** 2 + (self%a1%y-self%a3%y) ** 2)
  !   edge3 = sqrt((self%a2%x-self%a3%x) ** 2 + (self%a2%y-self%a3%y) ** 2)

  !   s = (edge1 + edge2 + edge3) / 2.0_DP

  !   res = sqrt(s * (s-edge1) * (s-edge2) * (s-edge3))
  ! end function area_triangle

  real(DP) function bilinear(A, p) result(res)
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
        print *, "Bilinear basis can't have index", A, " in 2D"
        stop
    end select

  end function bilinear

  real(DP) function DX_bilinear(A, p) result(res)
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
        print *, "Bilinear basis can't have index", A, " in 2D"
        stop
    end select

  end function DX_bilinear

  real(DP) function DY_bilinear(A, p) result(res)
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
        print *, "Bilinear basis can't have index", A, " in 2D"
        stop
    end select

  end function DY_bilinear

  real(DP) function linear_1D(A, p) result(res)
    integer  :: A
    real(DP) :: p
    
    if (A == 1) then
      res = (1_DP - p) / 2_DP
    else if (A == 2) then
      res = (1_DP + p) / 2_DP
    else 
      print *, "Linear Basis can't have index", A 
      stop
    end if
  end function linear_1D

  real(DP) function D_linear_1D(A, p) result(res)
    integer  :: A
    real(DP) :: p
    
    if (A == 1) then
      res = - 1_DP / 2_DP
    else if (A == 2) then
      res =   1_DP / 2_DP
    else 
      print *, "Linear Basis Derivative can't have index", A 
      stop
    end if
  end function D_linear_1D
end module quadrilateral_mod