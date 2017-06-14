subroutine analytic(c,id)
use globe_data
integer :: id
complex(8),dimension(nmax) :: s,sf
complex,dimension(nmax)  :: c
s=czero
sf=czero
s(1:npts) = dcmplx(dble(sig(1:npts,id)),0d0)
call dfftw_plan_dft_1d(plan,nn,s,sf,FFTW_FORWARD, FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
sf(1) = sf(1)/2.0d0
sf(nk) = dcmplx(dreal(sf(npts)),0.0d0)
sf(1:nk)=dble(2*dt)*sf(1:nk)
s=czero
sf(nk+1:nn)=czero
call dfftw_plan_dft_1d(plan,nn,sf,s,FFTW_BACKWARD, FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
!c(1:npts)=2.0*real(dreal(s(1:npts)))/nn/dt
c(nk+1:nn)=dcmplx(0d0,0d0)
c(1:nk)=cmplx(sngl(dreal(s(1:nk))),sngl(dimag(s(1:nk))))
end
