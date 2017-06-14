subroutine do_water(c,id1,id2)
use globe_data,only : ss,nn
integer :: id1,id2,i
real(8) :: maxa=0.d0,maxb
real(4) :: c
do i=1,nn
   maxa=amax1(maxa,dsqrt(dreal(ss(i,id1,id2))**2+dimag(ss(i,id1,id2))**2))
enddo
do i=1,nn
   amp=dsqrt(dreal(ss(i,id1,id2))**2+dimag(ss(i,id1,id2))**2)
   maxb=amp/amax1(amp,dble(c)*maxa)
   ss(i,id1,id2)=maxb*ss(i,id1,id2)
enddo
end subroutine
