module globe_data
integer,parameter :: nmax=4000000,nstmax=2000
integer,parameter :: FFTW_ESTIMATE=64,FFTW_MEASURE=1
integer,parameter :: FFTW_FORWARD=-1,FFTW_BACKWARD=1
integer :: nn,nk,npow,nf
integer :: npt2,nsmpl,npts,halfn
integer :: nagc
integer(8) :: plan
real(8) :: dom
real(4),dimension(4) :: fre
real(4),dimension(nmax,3) :: sig
real(4),dimension(nmax,3,3) :: signcc,sigpcc
real(4),dimension(nmax) :: sigt
real(8):: twopi=6.28318530717958,pi=datan2(1.0d0,1.0d0)*4.0d0
real(4):: wu=1.0,dt
complex(8),parameter:: czero=(0.0d0,0.0d0)
complex(8),dimension(nmax,3,3) :: ss
end module
