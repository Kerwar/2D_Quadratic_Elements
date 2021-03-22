module point_mod

  use number_types_mod, only: DP

  implicit none

  type :: Point
    real(DP) :: x, y 
  end type Point
end module point_mod