!2017/06/04
! introduce globe_data
subroutine do_whiten(id1,id2)
use globe_data
implicit none
integer i,id1,id2
real(8) f
real(8),dimension(nmax):: s

s=0d0
call smoothf(ss(:,id1,id2),s)         ! smooth the spectrum
do i=1,nk 
   f = dble((i-1))*dom
   if( f .ge. dble(fre(1)) .and. f .le. dble(fre(4)) ) then
      ss(i,id1,id2)=ss(i,id1,id2)/s(i)
   endif
enddo
call taperf(ss(:,id1,id2))
end subroutine
