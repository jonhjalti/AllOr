!AllOr Program for detecting breed of origin.
! Written by Jon H Eiriksson in October 2020
! 
! Needs phased genotypes and pedigree information
! Starts by finding haplotype matches in related animals

!Revised in January 2022
!Simplified and comments added in May 2022
!Version 8 in December 2023

program AllOr

implicit none
	  
integer :: i,j,k,l,m,n,o,p,q                           !loop counters and such
integer :: nosnp,nopbs,nocbs,noped,nopb,nocb           !number of snps, no of pure
integer :: nopbs2                                      !Number of pure breeds squared
integer :: control(3)                                  !Numbers to control the window length
integer :: localped(2,31)                              !For storing pedigree of each individual
integer :: norounds                                    !number of rounds in haplotype comparison
integer :: peddepth, nulocped                          !User controled depth of local pedigree to search for matches

integer,allocatable :: nopure(:)                       !number of reference genotypes
integer,allocatable :: nubrcr(:,:,:)                   !Defenition of crossbred type
integer,allocatable :: pbpoint(:,:),pbanp(:)           !pointer for purebred haplotypes,

character*1,allocatable :: pbhapl(:,:,:)                    !Store the haplotypes of the purebred. (animal, hapl(1/2),snp)
character*1,allocatable :: hapl(:,:)                      ! For keeping the haplotype of the crossbred animal working with each time

integer,allocatable :: numped(:,:)                      !The pedigree expressed as numbers
!integer,allocatable :: genpedpoint(:)                   !For pointing to the genotypes for animals in local pedigree
!integer,allocatable :: pbanbreed(:)                      

integer,allocatable :: roundvis(:,:),nummm(:),comout(:)  !Pointer for rounds in haplotype comparison, number of missing allowed, counting matces in windows
real,allocatable :: apix1(:,:,:),apix2(:,:),apix3(:,:)   !These are used in first round of the assignment to collect info from comparison
integer,allocatable :: breedo(:,:),breedo2(:,:)          !Arrays for the assigned breed for output

real :: high
real :: allowedmm,proport
real :: prod(31)
real :: mult(10)
integer :: loc(5),res(2)
integer :: sbrnow(10),dbrnow(10)
integer :: ndsbr(2,50)
integer :: brinfc
integer :: gapfill                                     !indicates if unassigned gaps should be filled based on neighbouring snps
integer :: mode                                        !code to turn off different parts of the assignment. Mainly for looking for errors. 
integer :: ierr
integer :: nubrn(2)
	  
	  
character*60 :: pedfile,boacodefile(2),haplfile(2)     !Names of infiles
character*25,allocatable :: ped(:,:)                   ! For reading the pedigree file
character*10,allocatable :: pedbrcode(:)               !Used if pedigree should be used for find breed codes
!character*25,allocatable :: pbid(:)                    !ID from the boacode file for purebred
character*25,allocatable :: cbid(:)                    !ID from the boacode file for crossbred
character*25 :: misscode                               !How missing parent is coded
character*25 :: utformat,utformat2
character*25 :: hapform                       !format for reading haplotypes without space (FImpute style)
character*26 :: propformat
character*5 :: pedbrformat
character*10 :: breed(50)
character*10 :: breedc
character*10,allocatable :: breedcb(:)
character :: apixout,fastpopc,possoutr
character*150 :: readline
	  
	  !For defining crosses
integer,allocatable :: crossposs(:,:)
logical,allocatable :: crosslog(:)
logical,allocatable :: matchf(:),matchg(:)
logical,allocatable :: logass(:,:,:,:),komcode(:)
!Addition in 7.2
logical,allocatable :: allcom(:,:,:,:)
!Addition ends

logical :: breedlog(2,8)
logical :: pedhaplog(2,32)
logical :: dummy(2)
logical :: compare
logical :: apout,fastpop,possout
logical :: cycleerror
logical :: pedlog

logical :: step3=.TRUE.
logical :: useped=.TRUE.
logical :: usepop=.TRUE.
logical :: usebc=.TRUE.
logical :: hapread=.FALSE.

real,allocatable :: propanim(:,:)                            !the proportion of alleles of an animal that are traced to breed of origin
real,allocatable :: propsnp(:,:)  	  !for the proportion of alleles at each SNPs that are traced to breed of origin
real :: propall
integer,allocatable :: notasssnp(:,:)                   !For summing up the number of gamets that are not assigned at each SNP position
integer :: notassani(2)                   !Summing number of SNPs not assigned for this individual
integer :: asscount(2,9)

! ---------added dec 2023
logical :: addpure(10)=.FALSE.
integer :: nophap(10)
integer :: maxpure, maxhap
integer :: kr
integer :: locgtpoint(2,31,2,10)                             !Pointing to haplotypes from the local pedigree
integer,allocatable :: genpedpoint(:,:,:),inthap(:)                   !For pointing to the genotypes for animals in local pedigree
character*25,allocatable :: pbid(:,:)                    !ID from the boacode file for purebred
character*25 :: readid
!---------
	  
!=========================================READ CONTROL FILE=================================================================

      !Reading the control file
	  !Allocating too
	  
open(10,file='control.txt')
open(20,file='allor.log')
	  
read(10,*)nosnp
read(10,*)nopbs
write(20,*)nosnp,' markers considered'
write(20,*)nopbs,' pure breeds included'
if(nosnp.gt.999.and.nosnp.lt.10000)write(hapform,'(a1,i4,a3)')'(',nosnp,'a1)'
if(nosnp.ge.10000)write(hapform,'(a1,i5,a3)')'(',nosnp,'a1)'
if(nosnp.le.999)write(hapform,'(a1,i3,a3)')'(',nosnp,'a1)'

nopbs2=nopbs*nopbs
if(nosnp.gt.60000)then
  call errormessage(1,nosnp,0,0.0)
endif
if(nosnp.lt.200)then
  call errormessage(2,nosnp,0,0.0)
endif
if(nopbs.gt.8)then
  write(20,*)nopbs,'are too many breeds '
  call errormessage(1,nopbs,0,0.0)
endif
if(nopbs.lt.2) call errormessage(2,nopbs,0,0.0)


allocate(nopure(nopbs),nubrcr(2,50,nopbs))
allocate(pbpoint(nopbs,2),pbanp(nopbs))

allocate(crossposs(2,nopbs2),crosslog(nopbs2))

nopure=0
nubrcr=0
q=1
!---------------------- addition dec 2023
do i=1,nopbs
  read(10,*)breed(i),nopure(i)
  if(nopure(i).lt.0)then
    addpure(i)=.TRUE.
	nopure(i)=nopure(i)*(-1)
  endif
  write(20,*)breed(i),'maximum of ',nopure(i)*2,' expected haplotypes'
  if(addpure(i))write(20,*)'Haplotypes from crossbred will be added to heplotype library'
enddo


!----------------------

	  
read(10,*)nocbs

p=nocbs+nopbs
write(20,*)nocbs,' different types of crossbreds'
do i=1,nocbs
  read(10,*)breed(i+nopbs),ndsbr(1,i),(nubrcr(1,i,j),j=1,ndsbr(1,i))
  read(10,*)breedc,ndsbr(2,i),(nubrcr(2,i,j),j=1,ndsbr(2,i))
  if(breedc.ne.breed(i+nopbs))then
    write(*,*)'purebred defination error: ',breed(i),breedc
    call errormessage(5,0,0,0.0)
  endif
  write(20,*)breed(nopbs+i)
enddo
	  
read(10,*)noped
read(10,*)peddepth
  
!Read the filenames and number of animals
read(10,'(a)')pedfile
read(10,'(a)')misscode

!nulocped=peddepth*peddepth-1
if(peddepth.eq.1)nulocped=1
if(peddepth.eq.2)nulocped=3
if(peddepth.eq.3)nulocped=7
if(peddepth.eq.4)nulocped=15
if(peddepth.eq.5)nulocped=31

nopb=sum(nopure)
read(10,'(a)')boacodefile(1)
read(10,'(a)')haplfile(1)
read(10,*)nocb
read(10,'(a)')boacodefile(2)
read(10,'(a)')haplfile(2)

read(10,*)control(1),control(2)

	  
read(10,*)allowedmm
read(10,'(i2)')mode
read(10,*)possoutr
	 
	 !translate the mode codes
if(mode.lt.3.and.mode.ge.0)then
  gapfill=mode
elseif(mode.eq.4)then
  gapfill=0
  step3=.FALSE.
elseif(mode.eq.5)then
  gapfill=1
  step3=.FALSE.
elseif(mode.eq.6)then
  gapfill=2
  step3=.FALSE.
elseif(mode.eq.7)then
  gapfill=0
  useped=.FALSE.
elseif(mode.eq.8)then
  gapfill=1
  useped=.FALSE.
elseif(mode.eq.9)then
  gapfill=2
  useped=.FALSE.
elseif(mode.eq.10)then
  gapfill=0
  usepop=.FALSE.
elseif(mode.eq.11)then
  gapfill=1
  usepop=.FALSE.
elseif(mode.eq.12)then
  gapfill=2
  usepop=.FALSE.
elseif(mode.eq.13)then
  gapfill=0
  step3=.FALSE.
  useped=.FALSE.
elseif(mode.eq.14)then
  gapfill=1
  step3=.FALSE.
  useped=.FALSE.
elseif(mode.eq.15)then
  gapfill=2
  step3=.FALSE.
  useped=.FALSE.
elseif(mode.eq.16)then
  gapfill=0
  step3=.FALSE.
  usepop=.FALSE.
elseif(mode.eq.17)then
  gapfill=1
  step3=.FALSE.
  usepop=.FALSE.
elseif(mode.eq.18)then
  gapfill=2
  step3=.FALSE.
  usepop=.FALSE.
elseif(mode.eq.20)then
  usebc=.FALSE.
  gapfill=1
else
  write(*,*)'ERROR, unknown mode: ', mode
  write(20,*)'ERROR, unknown mode: ', mode
  goto 999
endif
	 

possout=.FALSE.

if(possoutr.eq.'y')possout=.TRUE.
	  	 
close(10)
write(20,*)'Reading control file done'	  

            
	 	  
!==============================================READ PED PURE================================================================
!             Read pedigree and phased genotypes of purebred

!allocate(pbhapl(nopb,2,nosnp))             !Allocating an array for the purebred haplotypes	  
!allocate(pbid(nopb))                       !Vector for IDs of purebred reference animals
!allocate(pbanbreed(nopb))      
	
!--------added dec 23
if(any(addpure))then
  maxpure=maxval(nopure)+nocb
else
  maxpure=maxval(nopure)
endif

allocate(pbhapl(2*maxpure,nopbs,nosnp))             !Allocating an array for the purebred haplotypes	  
allocate(pbid(nopb,nopbs))                       !Vector for IDs of purebred reference animals
allocate(hapl(2,nosnp))      !For keeping haplotypes of the current crossbred animal 
kr=nosnp/10
allocate(inthap(kr))
!--------
pbhapl='9'
	
open(11,file=boacodefile(1))
open(12,file=haplfile(1))
	  
ierr=0
i=0
do while (ierr.eq.0)
  read(11,*,iostat=ierr)readline              !Reading boacodefile to check if the number of animals matches expectations
  i=i+1
enddo
i=i-1

if(i.ne.nopb)then
  write(*,*)'Number of purebred in control.txt and ',boacodefile(1), ' do not match',i,nopb
  write(20,*)'Number of purebred in control.txt and ',boacodefile(1), ' do not match',i,nopb
  write(*,*)'Program stops'
  write(20,*)'Program stops'
  goto 999
endif
rewind(11)

read(12,'(a)')readline
if(readline(1:1).eq.' ')hapread=.TRUE.
if(readline(2:2).eq.' ')hapread=.TRUE.
rewind(12)

pbanp=pbpoint(:,1)
nophap=0

do i=1,nopb
  read(11,'(a)')readline
  read(readline,*)readid,breedc
  do j=1,nopbs+1                                            !finding the matchin purebred code  +1 new in 8.1
    l=j                                                   ! l is the number of the breed
    if(breedc.eq.breed(j))exit
  enddo
 
!-------------------------addition dec 23

  if(l.gt.nopbs)then
    write(20,*)'undefined purebreed: ',breedc
	cycle
  endif
  if(hapread)then
    read(12,*)hapl(1,:)
    read(12,*)hapl(2,:)
  else
    read(12,hapform)hapl(1,:)
    read(12,hapform)hapl(2,:)
  endif
  do j=1,2
    m=1
    do k=1,nosnp,10
	  read(hapl(j,k),'(i1)')inthap(m)
	  m=m+1
	enddo
    if(count(inthap.lt.3).gt.kr/2)then     !Check if more than half is missing.
      nophap(l)=nophap(l)+1
	  if(nophap(l).gt.nopure(l)*2)then
        write(*,*)'ERROR, too many purebred ',breed(l)
        goto 999
      endif
      pbhapl(nophap(l),l,:)=hapl(j,:)
	  pbid(nophap(l),l)=readid
    endif
  enddo
  !-----
 
  if(nophap(l).gt.nopure(l)*2)then
    write(*,*)'ERROR, too many purebred ',breed(l)
    goto 999
  endif

enddo
	   
close(11)
close(12)

!------------added dec 2023
do i=1,nopbs
  write(20,*)nophap(i),breed(i),' haplotypes'   
  mult(i)=1.0/real(nophap(i))
enddo
if(any(addpure))then
  maxhap=maxval(nophap)+2*nocb
else
  maxhap=maxval(nophap)
endif

!----------------	  
	  
!Read pedigree
open(11,file=pedfile)

o=3                                                   !This number here means that o first letters in ID are breed codes
write(pedbrformat,'(a2,i1,a1)')'(a',o,')'
ierr=0
i=0

do while (ierr.eq.0)
  read(11,*,iostat=ierr)readline
  i=i+1
enddo
noped=i-1
write(20,*)noped,' lines in pedigree file'
allocate(ped(3,noped),numped(2,noped))
allocate(pedbrcode(noped),genpedpoint(noped,2,nopbs))
genpedpoint=0
numped=0
rewind(11)


do i=1,noped
  if(mod(i,10000).eq.0)write(20,*)i,'read from ped'
  read(11,'(a)')readline
  read(readline,*)ped(:,i)
  read(ped(1,i),pedbrformat)pedbrcode(i)
enddo
write(20,*)noped,'animals read from pedigree'
	 
!Making the pedigree with running numbers for simpler pointing
do i=1,noped
  if(mod(i,10000).eq.0)write(20,*)i,'ped processed'
  if(ped(2,i).ne.misscode)then
    do j=1,noped
      if(ped(1,j).eq.ped(2,i))then
        numped(1,i)=j
        exit
      endif
    enddo
  endif

  if(ped(3,i).ne.misscode)then
    do j=1,noped
      if(ped(1,j).eq.ped(3,i))then
	    numped(2,i)=j
	    exit
	  endif
    enddo
  endif

  !loop for finding the genotypes of purebred animals in pedigree

  !------------------ added dec 23
  
  do l=1,nopbs
    do j=1,nophap(l)
      if(ped(1,i).eq.pbid(j,l))then
	    genpedpoint(i,1,l)=j    
        if(ped(1,i).eq.pbid(j+1,l))genpedpoint(i,2,l)=j+1		                   !check if there is another haplotype
	    exit
	  endif
	enddo
	if(genpedpoint(i,2,l).gt.0)exit
  enddo
  
  !----------------------  
  
  !Þarf að gera einu sinni fyrir hvert kyn og tékka hvort næsta setröð er frá sama. 3víð genpedpoint (nopure, 2 nopbs)
enddo
close(11)
	  
write(20,*)'pedigree processed'

!============================================READ CROSS========================================================
	  
	  !First defing some things
	  
!allocate(cbhapl(2,nosnp))      !For keeping haplotypes of the current crossbred animal
	  
!Calculating how many rounds are needed to cover the chromosome based on number fo markers and specified overlap and window length
!Also allocating
p=(control(1)/5)/control(2)
control(3)=control(1)-p*control(2)


o=2**(nopbs*2)
norounds=(nosnp-control(1))/control(2)+2+2*p
if(mod((nosnp-control(1)),control(2)).eq.0)norounds=norounds-1

allocate(roundvis(3,norounds),nummm(norounds),comout(norounds))
allocate(apix1(2,5,norounds),apix2(nopbs2,norounds))
allocate(apix3(nopbs2,nosnp))
allocate(breedo(2,nosnp),breedo2(2,nosnp))
allocate(logass(2,2,nopbs,nosnp),komcode(o))
!allocate(matchf(norounds),matchm(norounds),matchg(norounds))
allocate(matchf(norounds),matchg(norounds))
allocate(notasssnp(2,nosnp))
allocate(propanim(2,nocb))
allocate(propsnp(2,nosnp))
allocate(cbid(nocb),breedcb(nocb))
	  
allocate(allcom(maxhap,2,nopbs,norounds))
	  
!The roundvis holds the start and end points of the segments considered in each round
!Put to subroutine at some point


roundvis(1,1)=1
roundvis(2,1)=control(3)
j=1
do while(j.le.p)
  j=j+1
  roundvis(1,j)=1
  roundvis(2,j)=roundvis(2,j-1)+control(2)
enddo
	  
do while(roundvis(2,j).le.nosnp)
j=j+1
  roundvis(1,j)=roundvis(1,j-1)+control(2)
  roundvis(2,j)=roundvis(2,j-1)+control(2)
enddo
	  
do i=j,norounds
  roundvis(1,i)=roundvis(1,i-1)+control(2)
  roundvis(2,i)=nosnp
enddo
	  
write(20,*)'Number of cores',norounds

do i=1,norounds
  roundvis(3,i)=roundvis(2,i)-roundvis(1,i)+1    !roundvis 3 holds the number of markers in each round
  nummm(i)=int(roundvis(3,i)*(1.0-allowedmm)+0.5)!Number of allowed mismatces for each round is stored in nummm
enddo

!The values that go into the probability calculations for the pedigree matches. prod(1) is a parent, 2-3 granparents and so on
prod=1.0/32.0
prod(1)=0.5
prod(2:3)=0.25
prod(4:7)=1.0/8.0
prod(8:15)=1.0/16.0

	  
call crossposss(nopbs2,crossposs)	  !Possible breeds 
!	  crosslog=.FALSE.
	  	  
!Making formats for output
if(nosnp.gt.999.and.nosnp.lt.10000)write(utformat,'(a1,i4,a8)')'(',nosnp,'(1x,i1))'
       
if(nosnp.ge.10000)write(utformat,'(a1,i5,a8)')'(',nosnp,'(1x,i1))'
       
if(nosnp.le.999)write(utformat,'(a1,i3,a8)')'(,',nosnp,'(1x,i1))'
       
!Format for the breedoforigin2 output. It needs to be longer for more than 4 breeds
if(nopbs.lt.5)then
  if(nosnp.gt.999.and.nosnp.lt.10000)write(utformat2,'(a1,i4,a8)')'(',nosnp,'(1x,i3))'
  if(nosnp.ge.10000)write(utformat2,'(a1,i5,a8)')'(',nosnp,'(1x,i3))'
  if(nosnp.le.999)write(utformat2,'(a1,i3,a8)')'(,',nosnp,'(1x,i3))'
elseif(nopbs.lt.7)then
  if(nosnp.gt.999.and.nosnp.lt.10000)write(utformat2,'(a1,i4,a8)')'(',nosnp,'(1x,i4))'
  if(nosnp.ge.10000)write(utformat2,'(a1,i5,a8)')'(',nosnp,'(1x,i4))'
  if(nosnp.le.999)write(utformat2,'(a1,i3,a8)')'(,',nosnp,'(1x,i4))'
else
  if(nosnp.gt.999.and.nosnp.lt.10000)write(utformat2,'(a1,i4,a8)')'(',nosnp,'(1x,i5))'
  if(nosnp.ge.10000)write(utformat2,'(a1,i5,a8)')'(',nosnp,'(1x,i5))'
  if(nosnp.le.999)write(utformat2,'(a1,i3,a8)')'(,',nosnp,'(1x,i5))'
 endif
	  
	  
o=nopbs*2+2
if(o.lt.10)then
  write(propformat,'(a16,i1,a8)')'(a25,2(1x,f5.1),',o,'(1x,i6))'
else
  write(propformat,'(a16,i2,a8)')'(a25,2(1x,f5.1),',o,'(1x,i6))'
endif
	  
write(20,*)'Starting reading haplotypes of crossbred'
!reading the animals information
open(11,file=boacodefile(2))
open(12,file=haplfile(2))
open(21,file='breedoforigin.txt')
if(possout)open(27,file='breedoforigin2.txt')
open(23,file='propassigned_animal.txt')
notasssnp=0
komcode=.FALSE.

!=====================================================Here begins the main loop of crossbred animals=============================================

!run through crossbred boacode file
!looking for errors


p=0
ierr=0
do i=1, nocb
  read(11,'(a)',iostat=ierr)readline
  if(ierr.ne.0.and.i.lt.nocb)then
    write(20,*)'Fewer animals in ',boacodefile(2),'than control.txt specifies:',i,'<',nocb
    write(*,*)'Fewer animals in ',boacodefile(2),'than expected:',i,'<',nocb
    write(20,*)'Program terminated'
    write(*,*)'Program terminated'
	goto 999
  endif
  read(readline,*)cbid(i),breedcb(i)
  do j=1,nocbs
    if(breedcb(i).eq.breed(j+nopbs))then
      exit
    endif
    if(j.eq.(nocbs))p=p+1
  enddo
  do j=1,noped
    if(cbid(i).eq.ped(1,j))then            !Finding the crossbred animals in the pedigree
  	  exit
    endif
    if(j.eq.noped)write(20,*)cbid(i),'not in pedigree'
  enddo
enddo
      
write(20,*)nocb,' Crossbred from ',boacodefile(2)
write(20,*)p,'with missing breed code'
write(20,*)'Start reading ',haplfile(2),' and assigning'


ierr=0
do i=1, nocb
  localped=0
  locgtpoint=0
  pedhaplog=.FALSE.
  crosslog=.FALSE.
  breedlog=.FALSE.
  breedo=0
  breedo2=0
  dbrnow=0
  sbrnow=0
  logass=.FALSE.
  cycleerror=.FALSE.
  pedlog=.FALSE.

  write(20,*)'Crossbred animal',i,cbid(i),breedcb(i)

        !Read the crossbred haplotypes
		
  if(hapread)then
    read(12,*,iostat=ierr)hapl(1,:)
    read(12,*,iostat=ierr)hapl(2,:)
  else
    read(12,hapform,iostat=ierr)hapl(1,:)
    read(12,hapform,iostat=ierr)hapl(2,:)
  endif
  if(ierr.ne.0.and.i.lt.nocb)then
    write(*,*)haplfile(2),'had fewer (or shorter) lines than expected:'
	write(*,*)i*2,'<',nocb*2
	write(*,*)'Program stops'
    write(20,*)haplfile(2),'had fewer (or shorter) lines than expected:'
	write(20,*)i*2,'<',nocb*2
	write(20,*)'Program stops'
	goto 999
  endif
  do j=1,nocbs
    if(breedcb(i).eq.breed(j+nopbs))then
      l=j
      nubrn(1)=ndsbr(1,l)
      nubrn(2)=ndsbr(2,l)
      do n=1,nubrn(1)
        sbrnow(n)=nubrcr(1,l,n)
      enddo
      do n=1,nubrn(2)
        dbrnow(n)=nubrcr(2,l,n)
      enddo
      exit
    endif
    if(j.eq.(nocbs))then
      pedlog=.TRUE.
  	write(20,*)'Breedcode of a crossbred not found, will use pedigree'
	write(20,*)'If that was not the intention pleas check the codes in control.txt and ',boacodefile(2)
    endif
  enddo

  call compall(hapl,pbhapl,nopbs,roundvis,allcom,nosnp,nophap,norounds,nummm,maxhap)
  
!Input: hapl of CB, hapl of all PB,roundvis,output,num of snp,num of PB, number of rounds,num of mismatces, length of pbhapl
  
  if(useped)then
    do j=1,noped
      if(cbid(i).eq.ped(1,j))then                          !finding the crossbred animal in the pedigree
        call pedprep(numped,j,nulocped,localped,noped)
        exit
      endif
      if(j.eq.noped)then
        write(20,*)cbid(i),'not in pedigree'
        if(pedlog)then
          pedlog=.FALSE.
          write(20,*)'no breedcode either!' 
          write(20,*)'dont expect much from the results' 
          nubrn(1)=nopbs
          nubrn(2)=nopbs
          do n=1,nopbs
            dbrnow(n)=n
            sbrnow(n)=n
          enddo
        endif
      endif
    enddo
  endif

  if(pedlog)then
    do j=1,noped
      if(cbid(i).eq.ped(1,j))then
        call codefromped(noped,nopbs,localped,pedbrcode,nulocped,sbrnow,dbrnow,nubrn,breed(1:nopbs))
      endif
    enddo
    write(20,*)'sire breeds:',(breed(sbrnow(j)),j=1,nubrn(1))
    write(20,*)'dam breeds:',(breed(dbrnow(j)),j=1,nubrn(2))
  endif
		
  if(.not.usebc)then      !If breedcode is not used all breeds are possible for both dam and sire
    do n=1,nopbs
      dbrnow(n)=n
      sbrnow(n)=n
    enddo
  endif

  do n=1,nubrn(1)
    breedlog(1,sbrnow(n))=.TRUE.
    logass(:,1,sbrnow(n),:)=.TRUE.
  enddo
  do n=1,nubrn(2)
    breedlog(2,dbrnow(n))=.TRUE.
    logass(:,2,dbrnow(n),:)=.TRUE.
  enddo
        
		
!      write(*,*)breedlog
    !make boolean vector indicating if breed combinations are allowed
  do k=1,nopbs2
    do n=1,nubrn(1)
      do o=1,nubrn(2)
        if(sbrnow(n).eq.crossposs(1,k).and.dbrnow(o).eq.crossposs(2,k).or.dbrnow(o).eq.crossposs(1,k).and.&
          sbrnow(n).eq.crossposs(2,k))then
          crosslog(k)=.TRUE.
        endif
  	enddo
    enddo
  enddo


  !-------------------------added dec 23
  if(useped)then
    do j=1,nulocped
	  do l=1,nopbs
        if(genpedpoint(localped(1,j),1,l).gt.0)then
		  locgtpoint(1,j,:,l)=genpedpoint(localped(1,j),:,l)
	      pedhaplog(1,j)=.TRUE.
        endif
        if(genpedpoint(localped(2,j),1,l).gt.0)then
          locgtpoint(2,j,:,l)=genpedpoint(localped(2,j),:,l)
	      pedhaplog(2,j)=.TRUE.
        endif
	  enddo
    enddo
  endif
  
  !---------------------------
  
  comout=0
  apix1=0.0
  apix2=0.0
  
  call getapix1(1,1)
  call getapix1(2,2)

  do j=1,norounds 
    do k=1,nopbs*nopbs
      if(breedlog(1,crossposs(1,k)).and.breedlog(2,crossposs(2,k)))apix2(k,j)=&
           apix1(1,crossposs(1,k),j)+apix1(2,crossposs(2,k),j)
    enddo
  enddo
  !Now for haplotype 1 being the maternal one
  apix1=0.0
  comout=0  

  call getapix1(1,2)
  call getapix1(2,1)
  
  !End NEW
  !Beginning the comparison
  !Starting with first haplotype, finding in pedigree of sire and then in the possible sire breeds.
  

		!Add the probabilities from the second round to apix2
  do j=1,norounds 
    do k=1,nopbs2
      high=apix1(1,crossposs(1,k),j)+apix1(2,crossposs(2,k),j)
      if(breedlog(1,crossposs(2,k)).and.breedlog(2,crossposs(1,k)).and.high.gt.apix2(k,j))apix2(k,j)=high
    enddo
  enddo

  apix3=0.0
  do j=1,nosnp         !Beginning loop for assigning snp round 1
    high=0.00000        !highest so far
    loc=0        !!highest position vector
    dummy=.TRUE.
    do l=1,nopbs2
      do k=1,norounds
        if(j.ge.roundvis(1,k).and.j.le.roundvis(2,k))then
          apix3(l,j)=apix3(l,j)+apix2(l,k)
        endif
      enddo
      if((apix3(l,j)-high).gt.0.0000001)then
        high=apix3(l,j)
        loc(1)=l
        loc(2:5)=0

      elseif(abs(apix3(l,j)-high).lt.0.0000001.and.apix3(l,j).gt.0.0000001)then
        if(loc(2).eq.0)then
          loc(2)=l
        elseif(loc(3).eq.0)then
          loc(3)=l
        elseif(loc(4).eq.0)then
          loc(4)=l
        elseif(loc(5).eq.0)then
          loc(5)=l
        endif
      endif
    enddo
    if(loc(1).gt.0.and.loc(2).eq.0)then
      call brass2(crossposs(:,loc(1)),breedo,logass,j,nosnp,nopbs)
    elseif(loc(2).gt.0)then
      do l=2,5
        if(loc(l).eq.0)exit
        if(crossposs(1,loc(1)).ne.crossposs(1,loc(l)))dummy(1)=.FALSE.
        if(crossposs(2,loc(1)).ne.crossposs(2,loc(l)))dummy(2)=.FALSE.
      enddo
      if(dummy(1).and..not.dummy(2))then
        breedo(:,j)=9
        logass(:,:,:,j)=.FALSE.
        logass(1,1,crossposs(1,loc(1)),j)=.TRUE.
        l=1
        do while(loc(l).gt.0)
          logass(2,1,crossposs(2,loc(l)),j)=.TRUE.
          l=l+1
        enddo
      elseif(dummy(2).and..not.dummy(1))then
        breedo(:,j)=9
        logass(:,:,:,j)=.FALSE.
        logass(2,1,crossposs(2,loc(1)),j)=.TRUE.
        l=1
        do while(loc(l).gt.0)
          logass(1,1,crossposs(1,loc(l)),j)=.TRUE.
          l=l+1
        enddo
      else
        breedo(:,j)=9
      endif
    else
      breedo(:,j)=9
    endif
  enddo
  
  ! fill in if a gap shorter than twice the length of the cores, if the same breed is assigned on both sides


  if(gapfill.gt.0)then
    call gapfi1(nosnp,breedo,logass,control(1),nopbs,1)
    call gapfi1(nosnp,breedo,logass,control(1),nopbs,2)
  endif
  if(gapfill.gt.1)then
    call gapfi1e(nosnp,breedo,logass,control(1),nopbs,1)
    call gapfi1e(nosnp,breedo,logass,control(1),nopbs,2)
  endif

  call logtoint(breedo,logass,nopbs,nosnp)

  if(gapfill.gt.0)then
    call gapfi1(nosnp,breedo,logass,control(1),nopbs,1)
    call gapfi1(nosnp,breedo,logass,control(1),nopbs,2)
  endif
  if(gapfill.gt.1)then
    call gapfi1e(nosnp,breedo,logass,control(1),nopbs,1)
    call gapfi1e(nosnp,breedo,logass,control(1),nopbs,2)
  endif
  
  call logtoint(breedo,logass,nopbs,nosnp)
  !Putting the sire breed on haplotype1 in homozygous snps that are still unassigned
  !Also for maternal if only one possible breed
  !And adding assignment to the other snp if possible
  
  
  if(step3)then
  if(dbrnow(2).eq.0.and.sbrnow(2).eq.0)then
    do j=1,nosnp
      if(breedo(1,j).lt.9.and.breedo(2,j).lt.9)cycle
      if(breedo(1,j).eq.9.and.breedo(2,j).eq.sbrnow(1))then
        breedo(1,j)=dbrnow(1)
      elseif(breedo(2,j).eq.9.and.breedo(1,j).eq.sbrnow(1))then
        breedo(2,j)=dbrnow(1)
      elseif(breedo(1,j).eq.9.and.breedo(2,j).eq.dbrnow(1))then
        breedo(1,j)=sbrnow(1)
      elseif(breedo(2,j).eq.9.and.breedo(1,j).eq.dbrnow(1))then
        breedo(2,j)=sbrnow(1)
      elseif(hapl(1,j).eq.hapl(2,j))then
        res(1)=sbrnow(1)
        res(2)=dbrnow(1)
        call brass2(res,breedo,logass,j,nosnp,nopbs)
      endif
    enddo
  elseif(sbrnow(2).eq.0)then      !if only one sire breed
    do j=1,nosnp
      if(breedo(1,j).lt.9.and.breedo(2,j).lt.9)cycle
      if(breedo(1,j).lt.9.and.breedo(1,j).ne.sbrnow(1))then
        call brass1(sbrnow(1),breedo(:,j),logass(:,:,:,j),nopbs,2)
      elseif(breedo(2,j).lt.9.and.breedo(2,j).ne.sbrnow(1))then
        call brass1(sbrnow(1),breedo(:,j),logass(:,:,:,j),nopbs,1)
      elseif(hapl(1,j).eq.hapl(2,j).and.breedo(1,j).gt.8.and.breedo(2,j).gt.8)then
        call brass1(sbrnow(1),breedo(:,j),logass(:,:,:,j),nopbs,1)
      endif
    enddo
  elseif(dbrnow(2).eq.0)then         !if only one dambreed
    do j=1,nosnp
      if(breedo(1,j).lt.9.and.breedo(2,j).lt.9)cycle
      if(breedo(1,j).lt.9.and.breedo(1,j).ne.dbrnow(1))then
        call brass1(dbrnow(1),breedo(:,j),logass(:,:,:,j),nopbs,2)
      elseif(breedo(2,j).lt.9.and.breedo(2,j).ne.dbrnow(1))then
        call brass1(dbrnow(1),breedo(:,j),logass(:,:,:,j),nopbs,1)
      elseif(hapl(1,j).eq.hapl(2,j).and.breedo(1,j).gt.8.and.breedo(2,j).gt.8)then
        call brass1(dbrnow(1),breedo(:,j),logass(:,:,:,j),nopbs,2)
      endif
    enddo
  endif
  endif
  
  call logtoint(breedo,logass,nopbs,nosnp)
  
  if(possout)call unasscode(logass,nosnp,nopbs,breedo2,komcode)
  
  write(21,utformat)breedo(1,:)
  write(21,utformat)breedo(2,:)
  if(possout)then
    write(27,utformat2)breedo2(1,:)
    write(27,utformat2)breedo2(2,:)
  endif
  
  !Calculate proportion assigned for this animal
  notassani=0
  asscount=0
  
  do j=1,nosnp
    asscount(1,breedo(1,j))=asscount(1,breedo(1,j))+1
    asscount(2,breedo(2,j))=asscount(2,breedo(2,j))+1
  enddo
  
  
  do j=1,nosnp
    if(breedo(1,j).gt.8)notassani(1)=notassani(1)+1
    if(breedo(2,j).gt.8)notassani(2)=notassani(2)+1
    if(breedo(1,j).gt.8)notasssnp(1,j)=notasssnp(1,j)+1
    if(breedo(2,j).gt.8)notasssnp(2,j)=notasssnp(2,j)+1
  enddo
  
  propanim(1,i)=100*(1.00-real(notassani(1))/real(nosnp))
  propanim(2,i)=100*(1.00-real(notassani(2))/real(nosnp))
  
  !---------------------------------Addition December 2023
  !Adding the haplotypes from the crossbred animals to the haplotype library 
  !If more than 98% is assigned in the haplotype and more than half to the same breed
  do k=1,2
    l=maxloc(asscount(k,1:nopbs),dim=1)
    if(addpure(l).and.propanim(k,i).gt.0.98)then
      if(maxval(asscount(k,1:nopbs)).gt.nosnp/2)then       
	    nophap(l)=nophap(l)+1
		pbid(nophap(l),l)=cbid(i)
	    do j=1,nosnp
		  if(breedo(k,j).eq.l)pbhapl(nophap(l),l,j)=hapl(k,j)		    
		enddo
		write(20,*)'Haplotype from crossbred animal ',cbid(i),'added to breed ',breed(l)
		mult(l)=1.0/real(nophap(l))
      endif
    endif
  enddo
   
  !--------------------------------------------
  
  write(23,propformat)cbid(i),propanim(1,i),propanim(2,i),(asscount(1,k),k=1,nopbs)&
  ,asscount(1,9),(asscount(2,l),l=1,nopbs),asscount(2,9)
enddo                                                            !Long animal loop closed
	  
!=============================================================================================================

close(11)
close(12)
close(21)
close(23)

if(possout)close(27)

if(possout)call writecode(nopbs,komcode)

open(22,file='snpsummary.txt')

propall=0

do i=1,nosnp                           !1s loop over SNP to calculate proportion assigned be SNP
  propsnp(1,i)=100-(100*real(notasssnp(1,i))/real(nocb))
  propsnp(2,i)=100-(100*real(notasssnp(2,i))/real(nocb))
  write(22,224)i,propsnp(1,i),propsnp(2,i)
  propall=propall+propsnp(1,i)+propsnp(2,i)
enddo                                  !1s closed
close(22)
propall=propall/2/real(nosnp)
write(20,*)propall
	  
	  
  223 format(a25,2(1x,f5.1),12(1x,i5))
  224 format(i8,2(1x,f5.1))

!==================FINISH=====================================================================================

write(20,*)'Program AllOr finished'
close(20)


	  
  999 stop
  
contains

subroutine getapix1(a,b)         !a is 1 for pateranal 2 for maternal, b is for haplotype 1 and 2
  integer :: a,b
  do l=1,nopbs                                                          !loop over all breeds
    matchf=.FALSE.                                                      !Start having found no matches
    if(.not.breedlog(a,l))cycle                                         ! Checks if breed l is one of the available sire breeds
    if(useped)then
      do j=1,nulocped                                                     !loop over the genotypes that were found in the pedigree
        if(locgtpoint(a,j,1,l).gt.0)then           
          do n=1,norounds
            if(.not.matchf(n).and.allcom(locgtpoint(a,j,1,l),b,l,n))then
			  apix1(b,l,n)=prod(j)           !allcom(maxhap,2,nopbs,norounds)
			  matchf(n)=.TRUE.
			endif
			if(.not.matchf(n).and.locgtpoint(a,j,2,l).gt.0)then  
              if(allcom(locgtpoint(a,j,2,l),b,l,n))then
			    apix1(b,l,n)=prod(j)           !allcom(maxhap,2,nopbs,norounds)
			    matchf(n)=.TRUE.
			  endif			
			endif
	      enddo
	    endif
        if(all(matchf))exit
      enddo
    endif
    if(all(matchf))cycle
    if(usepop)then
	  comout=count(allcom(1:nophap(l),b,l,:),dim=1)  
	  do n=1,norounds
	    if(matchf(n))cycle
		if(real(comout(n))*mult(l).ge.prod(nulocped)/2.0)then
		  apix1(b,l,n)=prod(nulocped)/2.0
		else
          apix1(b,l,n)=real(comout(n))*mult(l)
		endif
      enddo
    endif
  enddo
end subroutine getapix1

end program AllOr
	  
!========================================== Subroutines and functions=========================================
	  
subroutine errormessage(errcode,val1,val2,valr)
  implicit none
  integer errcode,val1,val2
  real valr
!		character*20 text
  write(*,*)'ERROR ',errcode
  if(errcode.eq.1)then
    write(*,*)val1, 'to high'
    stop
  endif
  if(errcode.eq.2)then
    write(*,*)val1, 'to low'
    stop
  endif
  if(errcode.eq.3)then
    write(*,*)val1,val2, 'not matching'
    stop
  endif
  if(errcode.eq.4)then
    write(*,*)'ERROR'
    stop
  endif
end subroutine errormessage
	  
	  
subroutine compall(anhap,pbset,nobr,rc,mats,n,m,l,nmm,mp)    
 !Input: hapl of CB, hapl of all PB,number of breeds,roundvis,output,num of snp,num of PB, number of rounds,num of mismatces
 !Subroutine for finding matching haplotypes in purebred haplotype library for single crossbred animal.
  implicit none
  integer :: n,l,j,mp
  integer :: nobr
  character*1,intent(in) :: anhap(2,n),pbset(mp,nobr,n)
  integer :: i,k,o,c,m(10)
  integer :: nmm(l),rc(3,l)
  logical,intent(out) :: mats(mp,2,nobr,l)          !m virkar ekki her
  logical :: chap(n)
  mats=.FALSE.
  do j=1,nobr         !across the pure breeds
    do i=1,m(j)           !across all pb haplotypes within breeds
	  do k=1,2         !two haplotypes of cb
        chap(:)=anhap(k,:).ne.pbset(i,j,:)     !comparing haplotypes
		do o=1,l   !across rounds
		  c=count(chap(rc(1,o):rc(2,o)))
		  mats(i,k,j,o)=c.le.nmm(o)
		enddo
      enddo
    enddo
  enddo
end subroutine compall
	  
subroutine compan(inn1,inn2,length,nr,rc,mr,outv,mf)
  implicit none
  integer :: length,nr,a
  character*1 :: inn1(length),inn2(2,length)
!  integer*2 :: inn1(length),inn2(2,length)
  integer :: rc(3,nr),mr(nr)
  integer :: outv(nr)
  logical :: mf(nr),compare
  
  outv=0
  do a=1,nr
    if(mf(a))cycle
    if(compare(inn1(rc(1,a):rc(2,a)),inn2(1,rc(1,a):rc(2,a)),rc(3,a),mr(nr)).or.compare(inn1(rc(1,a):rc(2,a)),inn2(2,&
         rc(1,a):rc(2,a)),rc(3,a),mr(nr)))then
      outv(a)=1
      mf(a)=.TRUE.
    endif
  enddo  
!	    return
end subroutine compan
	  
subroutine code(hapinn,haput,linn,lut,num)
  implicit none
  integer  linn,lut,num,a,b,c
  integer  hapinn(linn)
  integer  haput(lut)
  c=0
  haput=0
  do a=1,lut
    do b=1,num
      c=c+1
      haput(a)=haput(a)+2**b*hapinn(c)
      if(c.ge.linn)return
    enddo
  enddo
		
end subroutine code
	  
subroutine crossposss(a,b)
  implicit none
  integer a
  integer b(2,a)
  integer c,d,e,f
  c=int(sqrt(real(a)))
  f=0
  do d=1,c
    do e=1,c
	f=f+1
      b(1,f)=e
      b(2,f)=d
    enddo
  enddo
  return
	  
end subroutine crossposss
	  
subroutine codefromped(a,b,cim,dcv,e,fiv,giv,h,bre)                    !noped, nopbs,locped,pedbrcode,nolocped,output
!Finds out the possible breeds from pedigree information
  implicit none
  integer :: a,b,e,h(2)
  integer :: fiv(10),giv(10)
  integer :: cim(2,e)
  character*10 :: dcv(a),bre(b)
  integer :: i,j,k,l
  integer :: local(32)
  logical :: logbre(b)
  logical :: usebr
  fiv=0
  giv=0
  k=0
  l=0
  logbre=.FALSE.
  usebr=.FALSE.
  !first for sire
  do i=1,e
    if(i.eq.1.or.i.eq.3.or.i.eq.7.or.i.eq.15)k=0
    if(k.eq.0)l=0
    k=k+1
    if(cim(1,i).lt.1)cycle
    do j=1,b
      if(dcv(cim(1,i)).eq.bre(j))then
	    logbre(j)=.TRUE.
        l=l+1
	  endif
    enddo
    if(i.eq.1.or.i.eq.3.or.i.eq.7.or.i.eq.15.or.i.eq.31)then
      if(l.eq.k.or.l.gt.14)usebr=.TRUE.
    endif
  enddo
  j=0
  if(.not.usebr)logbre(1:b)=.TRUE.
  do i=1,b
    if(logbre(i))j=j+1
    if(logbre(i))fiv(j)=i
  enddo
  h(1)=j

  k=0
  l=0
  logbre=.FALSE.
  usebr=.FALSE.
  !first for sire
  do i=1,e
    if(i.eq.1.or.i.eq.3.or.i.eq.7.or.i.eq.15)k=0
    if(k.eq.0)l=0
    k=k+1
    if(cim(2,i).lt.1)cycle
    do j=1,b
      if(dcv(cim(2,i)).eq.bre(j))then
	    logbre(j)=.TRUE.
        l=l+1
	  endif
    enddo
    if(i.eq.1.or.i.eq.3.or.i.eq.7.or.i.eq.15.or.i.eq.31)then
      if(l.eq.k.or.l.gt.14)usebr=.TRUE.
    endif
  enddo
  j=0
  if(.not.usebr)logbre(1:b)=.TRUE.
  do i=1,b
    if(logbre(i))j=j+1
    if(logbre(i))giv(j)=i
  enddo
  h(2)=j
  return
end subroutine codefromped
	  
subroutine logtoint(intbr,logbr,br,snp)
  integer,intent(in) :: br,snp
  logical,intent(in) :: logbr(2,2,br,snp)
  integer :: intbr(2,snp)
  integer :: a,b,c,j
  do j=1,snp
    if(intbr(1,j).gt.8)then
      b=0
  	c=0
      do a=1,br
  	  if(logbr(1,1,a,j).or.logbr(1,2,a,j))then
  	    b=b+1
  		c=a
  	  endif
  	enddo
  	if(b.eq.1)intbr(1,j)=c
    endif
  enddo
  do j=1,snp
    if(intbr(2,j).gt.8)then
      b=0
      c=0
      do a=1,br
        if(logbr(2,1,a,j).or.logbr(2,2,a,j))then
          b=b+1
          c=a
        endif
      enddo
      if(b.eq.1)intbr(2,j)=c
    endif
  enddo
	    
end subroutine logtoint
	  
subroutine gapfi2(snp,intbr,logbr,wl,br)  !Filling in based on neighbouring, both chromosomes
  implicit none
  integer,intent(in) :: snp,wl,br
  integer :: intbr(2,snp)
  integer :: a,b,j
  logical :: logbr(2,2,br,snp)
  integer :: asspar(2,2)
	  
  a=0
  b=0
  asspar=0
  do j=1,snp
    if(intbr(1,j).lt.9)then
      a=j
      asspar(1,:)=intbr(:,j)
    elseif(a.gt.0.and.(j-a).lt.2*wl.and.b.lt.snp)then
      if(b.le.j)then
        do b=j,snp
          if(intbr(1,b).lt.9)then
            asspar(2,:)=intbr(:,b)
            exit
          endif
        enddo
      endif
      if(b.eq.snp)exit
      if((b-a).lt.(2*wl).and.asspar(1,1).eq.asspar(2,1).and.asspar(1,2).eq.asspar(2,2))then
        call brass2(asspar(1,:),intbr,logbr,j,snp,br)
      elseif((b-a).lt.(2*wl).and.asspar(2,1).gt.0)then
        logbr(:,:,:,j)=.FALSE.
        logbr(1,1,asspar(:,1),j)=.TRUE.
        logbr(2,1,asspar(:,2),j)=.TRUE.
      endif
    endif
  enddo
end subroutine gapfi2
	  
subroutine gapfi2e(snp,intbr,logbr,wl,br)  !Filling in based on neighbouring, both chromosomes, at the ends of chromosome
  implicit none
  integer,intent(in) :: snp,wl,br
  integer :: intbr(2,snp)
  integer :: a,b,j
  logical :: logbr(2,2,br,snp)
  integer :: asspar(2)

  b=0
  asspar=0
  if(intbr(1,1).gt.8)then
    do j=1,wl
      if(intbr(1,j).lt.9)then
        asspar(:)=intbr(:,j)
        b=j
        do a=1,b
          call brass2(asspar(:),intbr,logbr,a,snp,br)
        enddo
        exit
      endif
    enddo
  endif
  b=0
  asspar=0
  if(intbr(1,snp).gt.8)then
    do j=snp,snp-wl,-1
      if(intbr(1,j).lt.9)then
        asspar(:)=intbr(:,j)
        b=j
        do a=snp,b,-1
          call brass2(asspar(:),intbr,logbr,a,snp,br)
        enddo
        exit
      endif
    enddo
  endif

end subroutine gapfi2e
	  
subroutine gapfi1(snp,intbr,logbr,wl,br,c)         !neighbouring filling, one chromosome
  implicit none
  integer,intent(in) :: snp,wl,br,c
  integer :: intbr(2,snp)
  integer :: a,b,j
  logical :: logbr(2,2,br,snp)
  integer :: ass(2)

  a=0
  b=0
  ass=0
  do j=1,snp
    if(intbr(c,j).lt.9)then
      a=j
      ass(1)=intbr(c,j)
    elseif(a.gt.0.and.(j-a).lt.2*wl.and.b.lt.snp)then
      if(b.le.j)then
        do b=j,snp
          if(intbr(c,b).lt.9)then
            ass(2)=intbr(c,b)
            exit
          endif
        enddo
      endif
      if(b.eq.snp)exit
      if((b-a).lt.(2*wl).and.ass(1).eq.ass(2))then
        call brass1(ass(1),intbr(:,j),logbr(:,:,:,j),br,c)
      elseif((b-a).lt.(2*wl).and.ass(2).gt.0)then
        logbr(c,:,:,j)=.FALSE.
        logbr(c,1,ass(1),j)=.TRUE.
        logbr(c,1,ass(2),j)=.TRUE.
      endif
    endif
  enddo
end subroutine gapfi1
	  
subroutine gapfi1e(snp,intbr,logbr,wl,br,c)  !Filling in based on neighbouring, both chromosomes, at the ends of chromosome
  implicit none
  integer,intent(in) :: snp,wl,br,c
  integer :: intbr(2,snp)
  integer :: a,b,j
  logical :: logbr(2,2,br,snp)
  integer :: ass

  b=0
  ass=0
  if(intbr(c,1).gt.8)then
    do j=1,wl
      if(intbr(c,j).lt.9)then
        ass=intbr(c,j)
        b=j
        do a=1,b
          call brass1(ass,intbr(:,a),logbr(:,:,:,a),br,c)
        enddo
        exit
      endif
    enddo
  endif
  b=0
  ass=0
  if(intbr(c,snp).gt.8)then
    do j=snp,snp-wl,-1
      if(intbr(c,j).lt.9)then
        ass=intbr(c,j)
        b=j
        do a=snp,b,-1
          call brass1(ass,intbr(:,a),logbr(:,:,:,a),br,c)
        enddo
        exit
      endif
    enddo
  endif
	
  end subroutine gapfi1e
	  
subroutine pedprep(a,b,c,d,e)
  implicit none
  integer :: b,c,e
  integer :: a(2,e)
  integer :: d(2,c)
  integer :: f,g,h
  d(1,1)=a(1,b)
  d(2,1)=a(2,b)
  g=1
  h=1
  do f=2,c
    h=mod(f,2)+1
    if(d(1,g).gt.0)d(1,f)=a(h,d(1,g))
    if(d(2,g).gt.0)d(2,f)=a(h,d(2,g))
    if(h.eq.2)g=g+1
  enddo
  return
end subroutine pedprep
	  
subroutine unasscode(logi,a,b,numb,komid)
  implicit none
  integer :: a,b,c,d,e
  integer :: numb(2,a)
  logical :: logi(2,2,b,a)
  logical :: komid(2**(b+b))
  
  numb=0
  do e=1,2
    do c=1,a
      do d=1,b
        if(logi(e,1,d,c))numb(e,c)=numb(e,c)+2**(d-1)
        if(logi(e,2,d,c))numb(e,c)=numb(e,c)+2**(d+b-1)
      enddo
      komid(numb(e,c))=.TRUE.
    enddo
  enddo
  return
end subroutine unasscode
	  
subroutine writecode(b,komid)
  implicit none
  integer :: a,b,c
  logical :: komid(2**(b+b))
  logical :: logg(b*2)
  real :: prop(b),d(2)
  open(29,file='bo2code.txt')
  do a=1,2**(b+b)
    if(.not.komid(a))cycle
    logg=.FALSE.
    prop=0.0
    d=0.0
    do c=1,b*2
      if(mod(a,2**c).ge.2**(c-1))then
  	  logg(c)=.TRUE.
  	  if(c.le.b)d(1)=d(1)+1.0
  	  if(c.gt.b)d(2)=d(2)+1.0
  	endif
    enddo
    if(d(1).gt.0.001.and.d(2).gt.0.001)then
      d(1)=0.5/d(1)
      d(2)=0.5/d(2)
    elseif(d(1).lt.0.001)then
      d(2)=1.0/d(2)
    elseif(d(2).lt.0.001)then
      d(1)=1.0/d(1)
    endif
    do c=1,b
      if(logg(c))prop(c)=prop(c)+d(1)
  	if(logg(c+b))prop(c)=prop(c)+d(2)
    enddo
    write(29,*)a,(prop(c),c=1,b)
  enddo
  close(29)
  return
end subroutine writecode
	  
subroutine brass1(b,intbr,logbr,f,g)
  implicit none
  integer,intent(in) :: f,g,b
  integer :: d,h,co(2)
  integer :: intbr(2)
  logical :: logbr(2,2,f)
  
  intbr(g)=b
  logbr(g,:,:)=.FALSE.
  logbr(g,1,b)=.TRUE.
          
  h=2
  if(g.eq.2)h=1
  if(intbr(h).gt.8)then
    if(logbr(h,1,b).and.logbr(h,2,b))then
      co=0
  	do d=1,f
  	  if(logbr(h,1,d))co(1)=co(1)+1
  	  if(logbr(h,2,d))co(2)=co(2)+1
  	enddo
  	d=1
  	if(co(1).gt.co(2))d=2
  	logbr(h,d,b)=.FALSE.
    elseif(logbr(h,1,b))then
      logbr(h,1,b)=.FALSE.
    elseif(logbr(h,2,b))then
      logbr(h,2,b)=.FALSE.
    endif
  endif
end subroutine brass1
	  
subroutine brass2(b,intbr,logbr,c,e,f)
  implicit none
  integer :: c,d,e,f
  integer :: b(2),intbr(2,e)
  logical :: logbr(2,2,f,e)
  
  intbr(:,c)=b
  logbr(:,:,:,c)=.FALSE.
  logbr(1,1,b(1),c)=.TRUE.
  logbr(2,1,b(2),c)=.TRUE.
	  
end subroutine brass2
	  
logical function compare(inn1,inn2,length,maxwrong)
  implicit none
  integer :: length,a,b,maxwrong
  character*1 :: inn1(length),inn2(length)
!  integer*2 :: inn1(length),inn2(length)
  b=0
  compare=.TRUE.
  do a=1,length
    if(inn1(a).ne.inn2(a))b=b+1
    if(b.gt.maxwrong)then
      compare=.FALSE.
      exit
    endif
  enddo
end function compare
	  
