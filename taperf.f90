! Tapering subroutine itself in frequency domain
! 2017/06/04
subroutine taperf(sf)
use globe_data
integer   i,j,id
real(8)   f,d1,d2,dpi,ssa,s(nmax)
complex(8)   sf(nmax)
! ---
dpi = datan(1.0d0)*4.0d0
s=0d0
do i = 1,nk
   f = dble(i-1)*dom
   if(f.le.dble(fre(1))) then
      cycle
   else if(f.le.dble(fre(2))) then
      d1 = dpi/(dble(fre(2))-dble(fre(1)))
      ssa = 1.0d0
      do j = 1,npow
         ssa = ssa*(1.d0-dcos(d1*(dble(fre(1))-f)))/2.0d0
      enddo
      s(i) = ssa
   else if(f.le.dble(fre(3))) then
      s(i) = 1.0d0
   else if(f.le.dble(fre(4))) then
      d2 = dpi/(dble(fre(4))-dble(fre(3)))
      ssa = 1.0d0
      do j = 1,npow
         ssa = ssa*(1.0d0+dcos(d2*(fre(3)-f)))/2.0d0
      enddo
      s(i) = ssa
   endif
enddo
do i=1,nk
   sf(i) = sf(i)*s(i)
enddo
return
end subroutine
