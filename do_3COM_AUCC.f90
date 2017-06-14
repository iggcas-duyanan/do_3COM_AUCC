! do auto correlation  junxie
! 2017/06/06-- introduce globe_data and sacio
program do_3COM_AUCC
use sacio
use globe_data
implicit none
type(sac_head) :: sachead(3)
real t1,beg
real f1,f2,c,t_agc
complex(8),dimension(nmax,3,3) :: sss
character (180)command
character (3)com(3),comm
character (20)year_day
character (180)sacf(3)
character (2) net(nstmax)
character (12)sta(nstmax),kst
character (180)dirinn,dirout,output,inpar,sta_list,output_tmp,outdir
integer jday,itemp,nerr
integer tmpn1,tmpn2,doagc
integer multp,dsec,dseg,do_crs
integer begday,endday,ic1,ic2,ic
integer iy,id,ih,ist,nsamp,is,it
integer nzhour,nzmin,nzsec,icheck
integer nft,donorm,do_whitten,n1,nft2
integer year_b,year_e,day_b,day_e,nseg
integer nst,i,nh,j,nnh,nlen,ncom,ib,ie
integer dofilter,dopcc,doccc,iseg,do_n,dowater
logical ext
if(iargc().ne.1)then
   call usage
   stop
endif
call getarg(1,inpar)
inquire(file=inpar,exist=ext)
if(.not.ext)stop 'Hey dude, miss somthing?'
open(10,file=inpar)
read(10,'(a80)')sta_list               ! station list
read(10,*)year_b,day_b,year_e,day_e    ! begin time and end time
read(10,*)dsec,multp,npt2              ! segment length in seconds, overlaping percentage
read(10,*)ncom,comm,do_crs             ! number of component (1 or 3) ! component (if cn eq 1; comm=BHZ else comm=BH) 
read(10,*)f1,f2,dofilter,donorm        ! band pass filter frequency do filter or not, do whittening (before aucc)or not
read(10,*)dopcc,do_whitten             ! 1: output stacked file; 0: output daily cc
read(10,*)dowater,c                    ! do water level normalisation
read(10,*)doagc,t_agc                  ! do agc
read(10,*)t1,npts
read(10,*)dirinn                       ! sac file directory
read(10,*)dirout                       ! output cc file directory
close(10)
if(ncom.ne.1)then
   com(1)=trim(comm)//'Z'
   com(2)=trim(comm)//'N'
   com(3)=trim(comm)//'E'
else
   com(1)=comm
endif
fre(1)=0.8*f1
fre(2)=f1
fre(3)=f2
fre(4)=f2*1.5

open(11,file=sta_list)                                                       ! read in station list
do ist=1,nstmax
   read(11,*,err=13,end=13) net(ist),sta(ist)
enddo
13 close(11)
nst=ist-1

nn=2
npow=1
do while(nn.lt.npts)
   nn=nn*2
   npow=npow+1
enddo
nk = nn/2+1
write(*,'("There are ",i0," station pairs")')nst
write(*,'("Do AUCC from ",i0,"/",i0," to ",i0,"/",i0)')year_b,day_b,year_e,day_e
!if(dofilter.eq.1)write(*,*) "Filter the waveform"
nsmpl=2*npt2+1
dseg=int((1-real(multp)/100.0)*dsec) ! the left points without overlapping
nseg=int((86400-dsec)/dseg)+1         ! number of segments per day.
do is=1,nst                                                                    ! loop over station pair
   write(*,'("Do autocorrelation for station",1x,1a,1x,1a)')trim(net(is)),trim(sta(is))
   write(outdir,'(1a,"/",1a,"_",1a)')trim(dirout),trim(net(is)),trim(sta(is)) ! mkdir for outdir
   write(command,'("mkdir -p ",1a,1x," 2>/dev/null")')trim(outdir)
   do iy=year_b,year_e                                                        ! loop over year
      jday=365
      if(mod(iy,4).eq.0.and.mod(iy,100).ne.0.or.mod(iy,400).eq.0)jday=366
      endday=day_e
      if(iy.ne.year_e)endday=jday
      begday=day_b
      if(iy.ne.year_b)begday=1
      do id=begday,endday                                                  ! loop over day
         write(year_day,'(i0,"_",i3.3)')iy,id
         call system(command)
         do iseg=1,nseg                                                  ! loop over hour segment
            nzhour=(iseg-1)*dseg/3600
            nzmin=mod((iseg-1)*dseg,3600)/60
            nzsec=mod(mod((iseg-1)*dseg,3600),60)
            icheck=0
            do ic=1,ncom                                         ! loop over first component
               write(sacf(ic),'(1a,"/",1a,"/",1a,"_",i2.2,"_",i2.2,"_",i2.2,"_",1a,"_",1a,"_",1a,".SAC")')&
               trim(dirinn),trim(year_day),trim(year_day),nzhour,nzmin,nzsec,trim(net(is)),trim(sta(is)),trim(com(ic))
               inquire(file=sacf(ic),exist=ext)
               icheck=ic
            enddo
            if(icheck.ne.ncom)cycle                   ! check whether all components exist
            do ic=1,ncom                              ! read all data
               call read_sachead(trim(sacf(ic)),sachead(ic),nerr)
               if(ic.eq.1)n1=int(t1/sachead(ic)%delta)+1
               call read_sac(trim(sacf(ic)),sigt(:),sachead(ic),nerr)
               if(nerr.eq.-1)exit
               sig(1:npts,ic)=sigt(n1:n1+npts-1)
            enddo
            if(nerr.eq.-1)cycle                       ! read sac file fails
            dt=sachead(1)%delta
            beg=-npt2*dt         ! begin time of cc
            halfn=int(64/dt)
            dom=dble(1.0/nn/dt)
            nf=int(0.01d0/dom)
            if(doagc.eq.1)nagc=int(t_agc/dt/2)
            if(dofilter.eq.1)then
               do ic=1,ncom
                  call filter(sig(:,ic),ic)
               enddo
            endif
            if(donorm.eq.1)then
               do ic=1,ncom
                  call do_norm(ic)
               enddo
            endif
            sss=czero
            do ic1=1,ncom      ! loop over component one
               ib=1
               ie=ncom
               if(do_crs.ne.1)then
                  ib=ic1
                  ie=ic1
               endif
               do ic2=ib,ie    ! loop over component two
                  call do_ncc(ic1,ic2) ! do normal cross correlation
                  if(dowater.eq.1)call do_water(c,ic1,ic2) ! do water level cross correlation
                  if(do_whitten.eq.1)call do_whiten(ic1,ic2) ! do frequency whittening
                  if(dopcc.eq.1)call do_pcc(ic1,ic2) ! do phase cross correlation
                   
                  call dfftw_plan_dft_1d(plan,nn,ss(:,ic1,ic2),sss(:,ic1,ic2),FFTW_BACKWARD, FFTW_ESTIMATE)
                  call dfftw_execute(plan)
                  call dfftw_destroy_plan(plan)
                  nft2=nn/2-npt2
                  do i=1,nsmpl
                     signcc(2*npt2+2-i,ic1,ic2)=-real(dreal(sss(nft2+i,ic1,ic2)))
                  enddo
                  if(doagc.eq.1)call do_agc(ic1,ic2) ! do average g?? control
                  write(output,'(1a,"/",1a,"_",i2.2,"_",i2.2,"_",i2.2,"_",1a,"_",1a,"_",1a,"_",1a,".SAC")')&
                          trim(outdir),trim(year_day),nzhour,nzmin,nzsec,trim(net(is)),trim(sta(is)),com(ic1),com(ic2)
                  call write_ncf_sac(output,signcc(:,ic1,ic2),sachead(1),sachead(2),nsmpl,1,nerr)
                  if(dopcc.eq.1)then
                     write(output,'(1a,"/pcc_",1a,"_",i2.2,"_",i2.2,"_",i2.2,"_",1a,"_",1a,"_",1a,"_",1a,".SAC")')&
                          trim(outdir),trim(year_day),nzhour,nzmin,nzsec,trim(net(is)),trim(sta(is)),com(ic1),com(ic2)
                     call write_ncf_sac(output,sigpcc(:,ic1,ic2),sachead(1),sachead(2),nsmpl,1,nerr)
                  endif
               enddo
            enddo
         enddo         ! loop over segments
      enddo            ! end loop day
   enddo               ! end loop year
enddo                  ! end loop station pair
end program
