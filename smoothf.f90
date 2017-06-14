! smoothing routine  in frquency domain 
! 2017/06/04
subroutine smoothf(sf,smf)
use globe_data
integer i,j,id
real(8),dimension(nmax) :: smf
real(8) :: summ,f
complex(8),dimension(nmax) :: sf(nmax)
smf = 0.0d0
do i = 1,nk
   f = dble((i-1))*dom
   summ=0d0
   if(f.ge.dble(fre(1)).and.f.le.dble(fre(4)))then
      do j=i-nf,i+nf
         summ=summ+dsqrt(dimag(sf(j))**2+dreal(sf(j))**2)
      enddo
      smf(i) = summ/dble(2.*nf+1.)
   endif
enddo
return
end subroutine
