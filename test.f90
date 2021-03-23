program test
  use point_mod, only: Point
  use number_types_mod, only: DP
    use quadrilateral_mod, only: Square

  implicit none 
  abstract interface 
    real(DP) function funct(p) result(res)
      import :: Point, DP
      type(Point) :: p
    end function funct

  end interface
  type(Square) :: po
  procedure(funct) :: testing

  po = Square(Point(1,2), Point(2,2), Point(2,3), Point(1,3))
  !po = Square(Point(-1,-1), Point(1,-1), Point(1,1), Point(-1,1))
  print *, po%integral(testing)

end program test

real(DP) function testing(p) result(res)
  use point_mod, only: Point
  use number_types_mod, only: DP

  implicit none
  type(Point) :: p

  res = (1.0_DP - (p%x+3.0_DP)/2.0_DP) / 2.0_DP * (1.0_DP - (p%y+5.0_DP)/2.0_DP) / 2.0_DP
end function testing