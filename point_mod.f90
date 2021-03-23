module point_mod

  use number_types_mod, only: DP

  implicit none

  type :: Point
    real(DP) :: x, y 
  end type Point

  interface Point 
    procedure :: newPoint
  end interface

  contains

  type(Point) function newPoint(x,y)
    real(DP) :: x, y
    
    newPoint%x = x
    newPoint%y = y

  end function newPoint
end module point_mod