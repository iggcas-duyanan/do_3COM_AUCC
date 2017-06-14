! 2017/06/05 --introduce globe_data
! npts is the length of sig1 and sig2
subroutine do_ncc(id1,id2)
!use globe_data,only: sig,npt2,npts,sig,nmax,ss
use globe_data
integer id1,id2
complex(8),dimension(nmax) :: s1,s2
complex(8),dimension(nmax) :: ss1,ss2
real(8) aa

s1=czero
s2=czero
ss1=czero
ss2=czero
s1(1:npts)=dcmplx(dble(sig(1:npts,id1)),0d0)
s2(1:npts)=dcmplx(dble(sig(1:npts,id2)),0d0)
call dfftw_plan_dft_1d(plan,npts,s1,ss1,FFTW_FORWARD, FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_1d(plan,npts,s2,ss2,FFTW_FORWARD, FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
aa=-1.0d0
do i=1,nk
   ss(i,id1,id2)=aa*dcmplx(ss1(i))*dcmplx(conjg(ss2(i)))
   aa=-aa
enddo
return
end subroutine
