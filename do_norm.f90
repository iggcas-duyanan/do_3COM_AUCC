! do normalisation using running absolute mean
! 2017/06/04 --introduce globe_data
subroutine do_norm(id)
use globe_data
implicit none
integer :: id
real(4),dimension(nmax):: x1
x1=sig(:,id)
call smooth(x1)      ! get the running absolute average waveform
sig(1:npts,id)=sig(1:npts,id)/x1(1:npts)
end subroutine
