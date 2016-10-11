! This program calculates mean square displacements from disp.out

module calcmsds 

IMPLICIT NONE

INTEGER :: n,i,j,k,nmsdlength,num,nspecies,nstep,nrun,nbin,nconfigs,i1,i2
INTEGER :: mcorrtime,nmsdcalltime
INTEGER, ALLOCATABLE, DIMENSION(:) :: ntype,numspc,normtot
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: norm
DOUBLE PRECISION :: dtime,bin,flag_disp
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xdisp,ydisp,zdisp,z
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: xdispstore,ydispstore,zdispstore
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: xmsd,ymsd,zmsd
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: xds,yds,zds
CHARACTER*20 :: rstfile,filename
CHARACTER*200 :: filein
CHARACTER(len=9) :: indi
LOGICAL :: overflow, restart,readfrominpt_log

contains

    subroutine msdconf

    IMPLICIT NONE

    rstfile='msdrst.dat'
    restart=.false.
    overflow=.false.
    readfrominpt_log=.false.

    indi='123456789'

    !!!!! Get parameters and allocates vectors and matrices !!!!!

    open(10,file='msd.inpt')
    read(10,'(a)') filein
    open(11,file=filein,status='old', form='formatted')

    read(10,*) nspecies
    ALLOCATE(numspc(nspecies))
    ALLOCATE(z(nspecies))
    num=0
    do i=1,nspecies
       read(10,*) numspc(i),z(i)
       num=num+numspc(i)
    end do
!    write(6,*) 'The total number of ions/atoms is:',num

    read(10,*) nmsdlength

    ALLOCATE(xdisp(num))
    ALLOCATE(ydisp(num))
    ALLOCATE(zdisp(num))
    ALLOCATE(xmsd(num,0:nmsdlength))
    ALLOCATE(ymsd(num,0:nmsdlength))
    ALLOCATE(zmsd(num,0:nmsdlength))
    ALLOCATE(xdispstore(num,0:nmsdlength))
    ALLOCATE(ydispstore(num,0:nmsdlength))
    ALLOCATE(zdispstore(num,0:nmsdlength))
    ALLOCATE(xds(0:nmsdlength,nspecies,nspecies))
    ALLOCATE(yds(0:nmsdlength,nspecies,nspecies))
    ALLOCATE(zds(0:nmsdlength,nspecies,nspecies))
    ALLOCATE(norm(num,0:nmsdlength))
    ALLOCATE(normtot(0:nmsdlength))
    ALLOCATE(ntype(num))

    read(10,*) restart
    read(10,*) readfrominpt_log
    if (readfrominpt_log) then
       read(10,*) nmsdcalltime
       read(10,*) dtime
       read(10,*) flag_disp
       read(10,*) nrun
    end if

    do i=1,nspecies
       filename='msd'//indi(i:i)//'.dat'
       open(20+i,file=filename)
       do j=i,nspecies
          filename='msdcollect'//indi(i:i)//indi(j:j)//'.dat'
          open(50+(i*nspecies)+j,file=filename)
       end do
    end do
    open(200,file='work.dat')
    open(201,file='nernst.dat')

    !!!!! Zero arrays !!!!!

    do i=1,num
       xdisp(i)=0.0d0
       ydisp(i)=0.0d0
       zdisp(i)=0.0d0
       do j=1,nmsdlength
          xdispstore(i,j)=0.0d0
          ydispstore(i,j)=0.0d0
          zdispstore(i,j)=0.0d0
       end do
    end do

    do i=1,num
       do j=0,nmsdlength
          xmsd(i,j)=0.0d0
          ymsd(i,j)=0.0d0
          zmsd(i,j)=0.0d0
          norm(i,j)=0
       end do
    end do

    do i=0,nmsdlength
       do j=1,nspecies
          do k=1,nspecies
             xds(i,j,k)=0.0d0
             yds(i,j,k)=0.0d0
             zds(i,j,k)=0.0d0
          end do
       end do
       normtot(i)=0
    end do

    !!!!! If there is a restart file !!!!!

    if (restart) then
       open(12,file=rstfile,status='old',form='formatted')
       do i=0,nmsdlength
          do j=1,num
             read(12,*) xmsd(j,i) 
             read(12,*) ymsd(j,i) 
             read(12,*) zmsd(j,i) 
             read(12,*) norm(j,i) 
          end do
          do j=1,nspecies
             do k=1,nspecies
                read(12,*) xds(i,j,k)
                read(12,*) yds(i,j,k)
                read(12,*) zds(i,j,k)
             end do
          end do
          read(12,*) normtot(i)
       end do
       close(12)
    end if

    !!!!! Types of the different species !!!!!

    if (readfrominpt_log) then
       k=1
       do i=1,nspecies
          do j=1,numspc(i)
             ntype(k)=i
             k=k+1   
          end do
       end do
    !   read(11,*) nbin,bin,nbin
    !   do i=1,num
    !      read(11,*) nbin
    !   end do
    else
       read(11,*) nmsdcalltime,dtime,nrun
       do i=1,num
          read(11,*) ntype(i)
       end do
    end if 
    close(10)

    nstep=0
    mcorrtime=1

    nconfigs=nrun/nmsdcalltime

    do n=1,nconfigs
       nstep=nstep+nmsdcalltime
       do j=1,num
          read(11,*) xdisp(j),ydisp(j),zdisp(j)
!          if(xdisp(j).gt.flag_disp) write(6,*)'warning',xdisp(j),j,nstep
!          if(ydisp(j).gt.flag_disp) write(6,*)'warning',ydisp(j),j,nstep
!          if(zdisp(j).gt.flag_disp) write(6,*)'warning',zdisp(j),j,nstep

       end do
!       write(6,*)'entering msdcalc',nstep
       call msdcalc
    end do
    close(11)

    !!!!! Write out restart file !!!!!

    call msdoutput

    do i=1,nspecies
       close(20+i)
       do j=1,nspecies
          close(50+(i*nspecies)+j)
       end do
    end do
    close(200)
    close(201)

    contains

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!! Subroutine msdcalc !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine msdcalc

        IMPLICIT NONE

        INTEGER :: ipoint
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: xtot,ytot,ztot

        ALLOCATE(xtot(0:nmsdlength,nspecies))
        ALLOCATE(ytot(0:nmsdlength,nspecies))
        ALLOCATE(ztot(0:nmsdlength,nspecies))

        do i=1,num
           xdispstore(i,mcorrtime)=0.0d0
           ydispstore(i,mcorrtime)=0.0d0
           zdispstore(i,mcorrtime)=0.0d0
        end do

        !!!!! Accumulate displacements !!!!!
        do j=1,mcorrtime
           do i=1,num
              xdispstore(i,j)=xdispstore(i,j)+xdisp(i)
              ydispstore(i,j)=ydispstore(i,j)+ydisp(i)
              zdispstore(i,j)=zdispstore(i,j)+zdisp(i)
           end do
        end do

        if (overflow) then
           do j=mcorrtime+1,nmsdlength
              do i=1,num
                 xdispstore(i,j)=xdispstore(i,j)+xdisp(i)
                 ydispstore(i,j)=ydispstore(i,j)+ydisp(i)
                 zdispstore(i,j)=zdispstore(i,j)+zdisp(i)
              end do
           end do
        end if

        !!!!! Zero displacement arrays !!!!!

        do i=1,num
           xdisp(i)=0.0d0
           ydisp(i)=0.0d0
           zdisp(i)=0.0d0
        end do

        do j=1,mcorrtime
           k=mcorrtime-j
           do i=1,nspecies
              xtot(j,i)=0.0d0
              ytot(j,i)=0.0d0
              ztot(j,i)=0.0d0
           end do
           do i=1,num
              ipoint=ntype(i)
              xtot(j,ipoint)=xtot(j,ipoint)+xdispstore(i,j)
              ytot(j,ipoint)=ytot(j,ipoint)+ydispstore(i,j)
              ztot(j,ipoint)=ztot(j,ipoint)+zdispstore(i,j)
              norm(i,k)=norm(i,k)+1
              xmsd(i,k)=xmsd(i,k)+xdispstore(i,j)**2
              ymsd(i,k)=ymsd(i,k)+ydispstore(i,j)**2
              zmsd(i,k)=zmsd(i,k)+zdispstore(i,j)**2
           end do
           normtot(k)=normtot(k)+1
           do i1=1,nspecies
              do i2=i1,nspecies
                 xds(k,i1,i2)=xds(k,i1,i2)+xtot(j,i1)*xtot(j,i2)
                 yds(k,i1,i2)=yds(k,i1,i2)+ytot(j,i1)*ytot(j,i2)
                 zds(k,i1,i2)=zds(k,i1,i2)+ztot(j,i1)*ztot(j,i2)
              end do
           end do
        end do

        if (overflow) then
           do j=mcorrtime+1,nmsdlength
              k=mcorrtime-j+nmsdlength
              do i=1,nspecies
                 xtot(j,i)=0.0d0
                 ytot(j,i)=0.0d0
                 ztot(j,i)=0.0d0
              end do
              do i=1,num
                 ipoint=ntype(i)
                 xtot(j,ipoint)=xtot(j,ipoint)+xdispstore(i,j)
                 ytot(j,ipoint)=ytot(j,ipoint)+ydispstore(i,j)
                 ztot(j,ipoint)=ztot(j,ipoint)+zdispstore(i,j)
                 norm(i,k)=norm(i,k)+1
                 xmsd(i,k)=xmsd(i,k)+xdispstore(i,j)**2
                 ymsd(i,k)=ymsd(i,k)+ydispstore(i,j)**2
                 zmsd(i,k)=zmsd(i,k)+zdispstore(i,j)**2
              end do
              normtot(k)=normtot(k)+1
              do i1=1,nspecies
                 do i2=i1,nspecies
                    xds(k,i1,i2)=xds(k,i1,i2)+xtot(j,i1)*xtot(j,i2)
                    yds(k,i1,i2)=yds(k,i1,i2)+ytot(j,i1)*ytot(j,i2)
                    zds(k,i1,i2)=zds(k,i1,i2)+ztot(j,i1)*ztot(j,i2)
                 end do
              end do
           end do
        end if

        !!!!! Update array counters !!!!!

        if (mod(float(mcorrtime),float(nmsdlength)).eq.0) then
           overflow=.true.
        end if

        mcorrtime=int(mod(float(mcorrtime),float(nmsdlength)))
        mcorrtime=mcorrtime+1

        return

        end subroutine msdcalc

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!! Subroutine msdoutput !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine msdoutput

        IMPLICIT NONE

        INTEGER :: ipoint
        DOUBLE PRECISION :: dum1,dum2,wnernst,work,time,rn
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: somme
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: xtot,ytot,ztot,msd
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: msdcoll

        ALLOCATE(xtot(0:nmsdlength,nspecies))
        ALLOCATE(ytot(0:nmsdlength,nspecies))
        ALLOCATE(ztot(0:nmsdlength,nspecies))
        ALLOCATE(msd(0:nmsdlength,nspecies))
        ALLOCATE(msdcoll(0:nmsdlength,nspecies,nspecies))
        ALLOCATE(somme(nspecies))

        do i=0,nmsdlength
           do i1=1,nspecies
              msd(i,i1)=0.0d0
              do i2=1,nspecies
                 msdcoll(i,i1,i2)=0.0d0
              end do
           end do
        end do

        !!!!! Writing out values for restart file !!!!!

        open(14,file=rstfile,form='formatted')
        do i=0,nmsdlength
           do j=1,num
              write(14,*) xmsd(j,i)
              write(14,*) ymsd(j,i)
              write(14,*) zmsd(j,i)
              write(14,*) norm(j,i)
           end do
           do i1=1,nspecies
              do i2=1,nspecies
                 write(14,*) xds(i,i1,i2)   
                 write(14,*) yds(i,i1,i2)   
                 write(14,*) zds(i,i1,i2)   
              end do
           end do
           write(14,*) normtot(i)
        end do
        close(14)

        !!!!! Average over components and number of molecules !!!!!

        do i=1,num
           do j=0,nmsdlength
              if (norm(i,j).gt.0) then
                 xmsd(i,j)=xmsd(i,j)/float(norm(i,j))
                 ymsd(i,j)=ymsd(i,j)/float(norm(i,j))
                 zmsd(i,j)=zmsd(i,j)/float(norm(i,j))
              end if
           end do
        end do

        do j=0,nmsdlength
           do i1=1,nspecies
              do i2=1,nspecies
                 if (normtot(j).gt.0) then
                    rn=sqrt(float(numspc(i1)*numspc(i2)))*float(normtot(j))
                    xds(j,i1,i2)=xds(j,i1,i2)/rn
                    yds(j,i1,i2)=yds(j,i1,i2)/rn
                    zds(j,i1,i2)=zds(j,i1,i2)/rn
                    msdcoll(j,i1,i2)=xds(j,i1,i2)+yds(j,i1,i2)+zds(j,i1,i2)
                 end if
              end do
           end do
        end do

        do i=0,nmsdlength
           do i1=1,nspecies
              somme(i1)=0.0d0
           end do
           do j=1,num
              ipoint=ntype(j)
              somme(ipoint)=somme(ipoint)+xmsd(j,i)+ymsd(j,i)+zmsd(j,i)
           end do
           do i1=1,nspecies
              if (numspc(i1).gt.0) then
                  msd(i,i1)=somme(i1)/float(numspc(i1))
              end if
           end do
        end do

        dum1=0.0
        dum2=0.0

        do i1=1,nspecies
           write(20+i1,*) dum1,dum2
        end do
        write(200,*) dum1,dum2
        write(201,*) dum1,dum2

        do i1=1,nspecies
           if (numspc(i1).gt.0) then
              do i2=i1,nspecies
                 write(50+(i1*nspecies)+i2,*) dum1,dum2
              end do
           end if
        end do

        do i=0,nmsdlength-1
           time=(dble(i)+1)*dble(nmsdcalltime)*dtime*2.418d-5
           do i1=1,nspecies
              write(20+i1,*) time,msd(i,i1)
           end do
           work=0.0d0
           wnernst=0.0d0
           do i1=1,nspecies
              if (numspc(i1).gt.0) then
                 do i2=i1,nspecies
                    write(50+(i1*nspecies)+i2,*) time,msdcoll(i,i1,i2)
                    msdcoll(i,i1,i2)=msdcoll(i,i1,i2)*z(i1)*z(i2)
                    msdcoll(i,i1,i2)=msdcoll(i,i1,i2)*sqrt(float(numspc(i1)*numspc(i2)))/float(num)
                    if (i1.ne.i2) then
                       work=work+msdcoll(i,i1,i2)
                    end if
                    work=work+msdcoll(i,i1,i2)
                 end do
              end if
              wnernst=wnernst+msd(i,i1)*z(i1)*z(i1)*float(numspc(i1))/float(num)
           end do            
           write(200,*) time,work
           write(201,*) time,wnernst
        end do

!        write(6,*) 
!        write(6,*) '*** Mean squared displacements written out ***' 
!        write(6,*) 

        return

        end subroutine msdoutput

    end subroutine msdconf 

end module calcmsds 
