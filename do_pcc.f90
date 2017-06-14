subroutine do_pcc(id1,id2)
use globe_data
complex,dimension(nmax):: trace,pilot
complex,dimension(nmax):: csig1,csig2
integer ll,id1,id2
call analytic(csig1,id1)
call analytic(csig2,id2)
pilot(1:npts)=csig1(1:npts)/cabs(csig1(1:npts))
trace(1:npts)=csig2(1:npts)/cabs(csig2(1:npts))
call pcc(sigpcc(:,id1,id2),pilot,trace,1-npt2,1+npt2,1,npts,npts,ll,1)
end subroutine
