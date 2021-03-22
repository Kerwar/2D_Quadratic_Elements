module quadrature_mod
  
  use number_types_mod, only: DP 
  use point_mod, only: Point

  implicit none

  private
  public :: gauss_quad

  real(DP), parameter :: x1 = - 1_DP / sqrt(3.0), x2 = 1_DP / sqrt(3.0)
  real(DP), parameter :: w1 = 1_DP, w2 = 1_DP

  abstract interface 
    real(DP) function funct(p) result(res)
      import :: Point, DP
      type(Point) :: p
    end function funct

  end interface

  contains

  ! Gauss Quadrature from -1 to 1
  real(DP) function gauss_quad(f) result(res)

    procedure(funct) :: f
    res = f(Point(x1,x1)) * w1 * w1 + f(Point(x1,x2)) * w1 * w2 + &
    f(Point(x2,x1)) * w2 * w1 + f(Point(x2,x2)) * w2 * w2 
  end function
end module quadrature_mod 