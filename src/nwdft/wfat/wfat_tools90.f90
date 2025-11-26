elemental function wfat_isnan_e(x) result(w)
  use, intrinsic :: ieee_arithmetic
  implicit none
  double precision, intent(in) :: x
  logical                   :: w
  
  w = ieee_is_nan(x)
end function wfat_isnan_e


      
