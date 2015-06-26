  subroutine initnet (net_file)
    include 'real8.com'
    dimension inp(5),iop(5)
    include 'network.com'
    character*(*) net_file
!  ------------------------------------------------------------------
!      mat/ he,c,o,ne,mg,si,s,ar,ca,ti,cr,fe,ni/
!
!      data excess/-2.425d0,0.d0,4.737d0,7.043d0,13.932d0,21.492d0
!     2,26.015d0,30.231d0,34.847d0,37.546d0,42.818d0,48.335d0,53.903d0/
!         from netwinv
!      data excess/-2.425d0,0.d0,4.737d0,7.046d0,13.933d0,21.492d0
!     2,26.016d0,30.231d0,34.845d0,37.549d0,42.818d0,48.331d0,53.902d0/
!  ===================================================================
!                                     My Alpha Net
    nt=27
    open(nt,file=net_file,status='old')
    print*,' initnet-nt-file ',nt,net_file
!
    numreac=0
    inpt(:,:)=0
    iout(:,:)=0
    read(nt,*) nmat
    print*,' nmat ',nmat
    do i=1,nmat
       read(nt,*) ii,xmater(i),amater(i),zmater(i),excess(i)
       write(*,7) ii,xmater(i),amater(i),zmater(i),excess(i)
7      format(' init_net-i,x,a,z,e ',i4,a7,10es12.4)
    enddo
10  continue
    read (nt,*) in,(inp(i),i=1,in),io,(iop(i),i=1,io)
!     print*,' in,inp,,io,iop ',in,inp(1:in),io,iop(1:io)
    if(in.gt.0) then
       numreac=numreac+1
       npart(numreac)=in
       inpt(numreac,1:in)=inp(1:in)
       iout(numreac,1:io)=iop(1:io)
       read(nt,11)   (crate(i,numreac),i=1,7)
11     format(4e13.6)
       go to 10
    endif
100 continue
    close(nt)
    print*,' numreac= ',numreac
!     ---------------------------------------------------------------
!      do n=1,numreac
!     print 111,n,(crate(i,n),i=1,7)
! 111     format(i3,7es12.4)
!      enddo
! =================================================================
    if_screen=0
    tmp_nse=6.d9
    print 113,if_screen,tmp_nse
113 format(' initnet : if_screen=',i3,' tmp_nse=',es12.4)
    print*,' end initnet '
    return
  end
