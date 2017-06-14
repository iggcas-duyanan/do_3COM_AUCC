subroutine do_agc(id1,id2)
use globe_data
real(4) :: sum,temp(nmax),max_val=0
integer ii,jj,i,id1,id2,k
do ii=1,nsmpl
    k=0
    sum=0
    do jj=ii-nagc,ii+nagc
         if (jj-nagc.gt.0 .and. jj+nagc.lt.nsmpl+1) then
              sum=sum+signcc(jj,id1,id2)
              k=k+1
         endif
    enddo
    temp(ii)=sum/k
    if(max_val.lt.temp(ii))max_val=temp(ii)
enddo
do ii=1,npts
    signcc(ii,id1,id2)=signcc(ii,id1,id2)/(temp(ii)+0.001*max_val)
enddo
end subroutine
