! smooth the absolute value waveform
! 2017/06/04 --introduce globe_data
subroutine smooth(sigm)
use globe_data
integer(4) :: ii, jj, k
real,dimension(nmax) :: sigm(nmax),temp(nmax)
real:: summ
do ii=1,npts
   k=0;summ=0
   do jj=ii-halfn,ii+halfn
      if (jj.gt.0.and. jj.lt.npts+1) then
         summ=summ+abs(sigm(jj))
         k=k+1
      endif
   enddo
   temp(ii)=summ/k
enddo
sigm(1:npts)=temp(1:npts)
end subroutine
