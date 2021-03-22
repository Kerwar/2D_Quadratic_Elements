program test
  use point_mod, only: Point
  use number_types_mod, only: DP
  use quadrature_mod, only: gauss_quad
  use quadrilateral_mod, only: bilinear

  implicit none 
  abstract interface 
    real(DP) function funct(p) result(res)
      import :: Point, DP
      type(Point) :: p
    end function funct

  end interface
  type(Point) :: po
  procedure(funct) :: testing

  po%x = 1
  po%y = -1
  print *, gauss_quad(testing)
  print *, testing(Point(-1,-1.0)), testing(Point(-1.0/2.0,-1.0)), testing(Point(0,-1.0)), &
  testing(Point(1.0/2.0,-1.0)), testing(Point(1,-1.0))
end program test

real(DP) function testing(p) result(res)
  use point_mod, only: Point
  use number_types_mod, only: DP
  use quadrilateral_mod, only: bilinear

  implicit none
  type(Point) :: p

  res = (1_DP - p%x) / 2_DP * (1_DP - p%y) / 2_DP
end function testing