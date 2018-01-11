subroutine get_np(buff, np)

implicit none

character*80, intent(in) :: buff
integer, intent(out) :: np
character*5 :: buf

buf=buff(1:5)

select case (buf)
  case ('hexa8')
    np = 8
  case ('quad4')
    np = 4
  case default
    print *,'Element-type ',buf,' not yet defined'
    print *,'Abort'
    stop
end select

end subroutine get_np
