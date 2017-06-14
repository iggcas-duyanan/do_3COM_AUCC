! 2017/06/06
subroutine filter(sigm,id)
use globe_data,only : plan,nn,FFTW_FORWARD,FFTW_ESTIMATE,npts,nk,czero,FFTW_BACKWARD,dt,nmax,sig
implicit none
integer::    id
complex(8):: s(nmax),sf(nmax)
real,dimension(nmax) :: sigm
s=czero
sf=czero
s(1:npts) = dcmplx(dble(sig(1:npts,id)),0d0)
call dfftw_plan_dft_1d(plan,nn,s,sf,FFTW_FORWARD, FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
sf(1) = sf(1)/2.0d0
sf(nk) = dcmplx(dreal(sf(npts)),0.0d0)
sf(1:nk)=sf(1:nk)*dble(dt)
call taperf(sf)
s=czero
sf(nk+1:nn)=czero
! make forward FFT for seismogram: sf ==> s
call dfftw_plan_dft_1d(plan,nn,sf,s,FFTW_BACKWARD, FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
sigm(1:npts)=2.0*real(dreal(s(1:npts)))/nn/dt
return
end subroutine
