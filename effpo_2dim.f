c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c �ʎq�␳���v�Z����T�u���[�`��
c
c �����|�e���V�����̌v�Z(�Q����)
c 2007/1/24, Hara&Hasegawa, Ver.3.2
c
c	Ver1.0	�����[�X�o�[�W����
c	   1.1	���Z�X�\���ɑΉ��iSR�̂݁j
c	   1.2	�����|�e���V�����̌v�Z�̈���@�B�I�ɐݒ�
c	   1.3	���A�k�A�w���ꂼ��Effective Potential�����߂�
c	   1.4	�I�[�~�b�N�ڐG������Effective Potential�������Ă݂�
c	   2.0	��Ǒw�r����Classical�@�ƌ���
c	   2.1	�_�u�����Z�X�\���ɑΉ�
c	   2.2	�ޗ��������ɑΉ�(��MIT�̍\��)
c	   2.3	���E���b�V���̎�舵���̕ύX
c	   3.0	���E�ɋ߂Â��ɏ]���ă�������������iarea_type�ǉ�
c			�\�ʂ�classical�ȃ|�e���V�����Ə��X�Ɉ�v�����邽�߂̗̈挈��)
c	   3.1	�_�u�����Z�X�֘A�̌v�Z���@������ƕύX
c	   3.2	���������C��
c	�@ 3.3.1 exp�̌v�Z���e�[�u�����i�������j
c	   3.4	exp�̌v�Z�̓����l�̏ꏊ���܂Ƃ߂Čv�Z����i�������j
c
c	Input : dx, dz, lhet, iarea, hhm, eg, bktq, lxr1,lzr1,lxr2,lzr2, u, mtemp
c	Output: u_eff1, u_eff2, u_eff3
c
c	�d�l�F�ϕ��͈͂��v�Z�̈���o���ꍇ�̏ꍇ������
c	�ϕ��͈͂��v�Z�̈���o���ꍇ�́A�o���Ƃ���̕\�ʂ̊i�q�_��p����
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	subroutine qeffect(dx, dz, 
     &					lhet, iarea, hhm, eg,bktq,lxr1,lzr1,lxr2,lzr2,
     &					u, u_eff1,u_eff2,u_eff3, mtemp,dltec,ec)
	implicit none
c
c===����===
c---�ϐ��z��p�p�����[�^---
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c---�f�o�C�X�\��---
	real	dx,dz
	integer(2)	lhet(nlayer)
	integer(1)	iarea(nlayer)
c---�̈�ʃp�����[�^---
	real,dimension (nvalley,narea)	:: hhm
	real	eg(nvalley,narea),bktq(ntenum)
	real dhet(0:nlayer),dt
c---�f�o�C�X�����---
	real,dimension (0:nx,0:nz)	:: u				!classical�ȃ|�e���V����
	real,dimension (0:nx,0:nz)	:: u_eff1,u_eff2,u_eff3	!�����|�e���V�����l
	real ec(nvalley,narea)
c----(���Z�X)---
	integer(2),dimension (nrecess)	:: lxr1,lzr1,lxr2,lzr2
c---���M���p�z��---
	integer(2),dimension (0:nx,0:nz) :: mtemp
c----���[�J���ϐ�----
	integer   i,j, ii, n, ix, iz,x,xx,z,zz,ia,cflag
	integer,save ::	  x_exp, x_temp, z_exp, z_temp		!save�Y��
	real	  pi
	parameter(pi  = 3.141592654)
	integer ka,kv
	real,save :: a
	allocatable a(:)
	real(8)	alp
	real(8),save,allocatable		:: cp1,cp2,cp3,ep1,ep2,ep3
	dimension	::	cp1(:,:),cp2(:,:),cp3(:,:),ep1(:,:),ep2(:,:),ep3(:,:)
	integer(1),save,allocatable		::	flag
	dimension	::	flag(:,:)
	integer(1),save,allocatable		::	flag2
	dimension	::	flag2(:,:)
	integer(1),save,allocatable		::	atype
	dimension	::	atype(:,:)
	real(8),save,allocatable	 :: gauss
	dimension :: gauss(:,:,:,:)
	integer(1) atype2,atype3
	real, dimension	(nvalley,narea) :: dltec !delta Ec with Nextnano
c
c----(EP�@�v�Z�̂��߂̉�����)��----
	if (.not. allocated(a)) then	!�ŏ��������s
	  allocate(a(1:narea))		!�ޗ��������ɑΉ�(2006/12/22 Hara)
c
c	  do ia=1,narea		!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  do ia=1,2			!���ً}�G���[��� �v�C��2011/03/21��
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		  a(ia)=sqrt(hhm(1,ia)/bktq(1))/2.0
c
	  enddo
	  x_exp = ifix(maxval(a)*5.0/dx+0.5)	!EP�v�Z�̐ϕ��͈͂�
	  z_exp = ifix(maxval(a)*5.0/dz+0.5)	!���߂�v�Z
c
c----�v�Z�̈挈��-------------
c---	flag=0���v�Z�̈�Aflag=1���v�Z�̈�O-------
	  allocate(flag(0:nx,0:nz))
	  allocate(flag2(0:nx,0:nz))
	  flag=0
	  flag2=0
	  do iz=0,nz
	  do ix=0,nx
		do ii=1,nrecess
		  if(ix.gt.lxr1(ii).and.ix.lt.lxr2(ii)
     &		.and.iz.ge.lzr1(ii).and.iz.lt.lzr2(ii)) then	!�C��2007/1/9
			flag(ix,iz)=1								!�v�Z�̈挈��
		  endif

		  do z=iz-z_exp,iz+z_exp
		  do x=ix-x_exp,ix+x_exp

		    if(x<0)then		
			  if(flag2(ix,iz).eq.0)then
			    flag2(ix,iz)=1
			  elseif(flag2(ix,iz).eq.2)then
			  	flag2(ix,iz)=3
			  endif
			  		
		    elseif(x>nx)then
			  if(flag2(ix,iz).eq.0)then
			    flag2(ix,iz)=1
			  elseif(flag2(ix,iz).eq.2)then
			  	flag2(ix,iz)=3
			  endif
		    endif

			if(ix<lxr1(ii).and.iz>lzr1(ii).and.iz<lzr2(ii)	
     &								.and.x>lxr1(ii))then
			  if(flag2(ix,iz).eq.0)then
			    flag2(ix,iz)=1
			  elseif(flag2(ix,iz).eq.2)then
			  	flag2(ix,iz)=3
			  endif	
			  		  	
			elseif(ix>lxr2(ii).and.iz>lzr1(ii).and.iz<lzr2(ii)	
     &								.and.x<lxr2(ii))then			  
			  if(flag2(ix,iz).eq.0)then
			    flag2(ix,iz)=1
			  elseif(flag2(ix,iz).eq.2)then
			  	flag2(ix,iz)=3
			  endif			  	
			endif
			
		    if(z<0)then
			  if(flag2(ix,iz).eq.0)then
			    flag2(ix,iz)=2
			  elseif(flag2(ix,iz).eq.1)then
			  	flag2(ix,iz)=3
			  endif			  	

		    elseif(z>nz)then
			  if(flag2(ix,iz).eq.0)then
			    flag2(ix,iz)=2
			  elseif(flag2(ix,iz).eq.1)then
			  	flag2(ix,iz)=3
			  endif		
		    endif

			if(iz>lzr2(ii).and.ix>lxr1(ii).and.ix<lxr2(ii)
     &								.and.z<lzr2(ii))then
			  if(flag2(ix,iz).eq.0)then
			    flag2(ix,iz)=2
			  elseif(flag2(ix,iz).eq.1)then
			  	flag2(ix,iz)=3
			  endif		
			endif

		  enddo
		  enddo
		enddo
	  enddo
	  enddo
c
	  allocate(atype(0:nx,0:nz))
c----���E�ɋ߂Â��ɏ]���ă�������������v�Z
c
	  call area_type(lxr1,lzr1,lxr2,lzr2,atype)	
	  atype2=maxval(atype)
c
	  allocate(cp1(0:nx,0:nz),cp2(0:nx,0:nz),cp3(0:nx,0:nz))
	  allocate(ep1(0:nx,0:nz),ep2(0:nx,0:nz),ep3(0:nx,0:nz))

	  allocate(gauss(-x_exp:x_exp,-z_exp:z_exp,narea,0:atype2))				 
	  gauss=0.0
c----EP�v�Z���̎w���v�Z�������e�[�u����-------	
c
c	open(unit=571,file='ep1.txt')
c
	  do atype3=0,atype2
	  if(atype3.eq.1) then
		cycle
	  endif
c
c	  do ia = 1,narea
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  do ia=1,2  !���ً}�G���[��� (�ޗ�2��) �v�C��2011/03/21��
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		alp=dble(a(ia))		
		if(atype3.ne.0) then
		  alp=dble(atype3-1)*dz	!a�̒l��̈�ɂ���ĕω�
		endif			
	    do iz=-z_exp,z_exp
		do ix=-x_exp,x_exp
c
		  gauss(ix,iz,ia,atype3) = dexp(-((dble(ix)*dx/alp)**2
     &						+(dble(iz)*dz/alp)**2)/2.0D0)
     &						*dx*dz
c	write(571,*) ix,iz,ia,atype3,gauss(ix,iz,ia,atype3)
		enddo
  		enddo
	  enddo
	  enddo
	endif
c----��------------------
c
c	write(*,*) "1d"
c
	ep1=0.0D0
	ep2=0.0D0
	ep3=0.0D0
	u_eff1=0.0
	u_eff2=0.0
	u_eff3=0.0

	cp1=u
	cp2=u
	cp3=u

c----�w�e���ޗ����̃�Ec�̍l��-----
c	�I�[�~�b�N�w�̍ޗ�����i���݂�In0.53GaAs)	
c	do i=1,nlayer,2
c	do iz=lhet(i)+1,lhet(i+1)
	do iz=0,nz		!�C��2007/1/9
	  do i=nlayer,1,-1
		if (iz.le.lhet(i)) then
		  ka=iarea(i)
		endif
	  enddo
	do ix=0,nx
	  if (flag(ix,iz).eq.0) then	!�w�e���E�ʂ̃��d���̒���
			cp1(ix,iz) =	cp1(ix,iz) +dltec(1,ka)
			cp2(ix,iz) =	cp2(ix,iz) +dltec(1,ka) -ec(2,ka)
			cp3(ix,iz) =	cp3(ix,iz) +dltec(1,ka) -ec(3,ka)

	  endif
	enddo
	enddo
	
c----EP�v�Z�̐ϕ��v�Z��------
c	�[�������́i�w�e���E�ʂ̈ʒu+�ϕ��͈�)�܂ł��v�Z�̈�
	do iz=0,lhet(nlayer-1)+ifix(a(iarea(nlayer))*5.0/dz+0.5)		! �v�Z�̈�̐ݒ�
	do ix=0,nx

	  if(atype(ix,iz).eq.1) then		!�v�Z�̈�O�̓X�L�b�v
		cycle
	  endif

	  if(flag(ix,iz).eq.0) then		
		do i=nlayer,1,-1
		  if (iz.le.lhet(i)) then
			ka=iarea(i)
		  endif
		enddo

		alp=dble(a(ka))
		atype3 = atype(ix,iz)
		if(atype(ix,iz).ne.0) then
		  alp=dble(atype(ix,iz)-1)*dz	!a�̒l��̈�ɂ���ĕω�
		  atype3 = atype(ix,iz)
		endif

c�@�@�ϕ��͈͂��v�Z�̈���o�Ă��Ȃ��ꍇ�iflag2=0)
	  if(flag2(ix,iz).eq.0)then
		do z=iz-z_exp,iz-1
		  zz=iz-z

		do x=ix-x_exp,ix-1
		  xx=ix-x
	  
c----�ϕ��v�Z����-------
		  ep1(ix,iz)=ep1(ix,iz)+(cp1(x,z)+cp1(ix+xx,z)
     &			+cp1(x,iz+zz)+cp1(ix+xx,iz+zz))			
     &		    *gauss(xx,zz,ka,atype3)
		  ep2(ix,iz)=ep2(ix,iz)+(cp2(x,z)+cp2(ix+xx,z)
     &			+cp2(x,iz+zz)+cp2(ix+xx,iz+zz))		
     &		    *gauss(xx,zz,ka,atype3)
		  ep3(ix,iz)=ep3(ix,iz)+(cp3(x,z)+cp3(ix+xx,z)
     &			+cp3(x,iz+zz)+cp3(ix+xx,iz+zz))
     &		    *gauss(xx,zz,ka,atype3)
		enddo
		enddo
		do x=ix-x_exp,ix-1
		  xx=ix-x
		  ep1(ix,iz)=ep1(ix,iz)+(cp1(x,iz)+cp1(ix+xx,iz)
     &			+cp1(ix,iz-xx)+cp1(ix,iz+xx))
     &			*gauss(xx,0,ka,atype3)
		  ep2(ix,iz)=ep2(ix,iz)+(cp2(x,iz)+cp2(ix+xx,iz)
     &			+cp2(ix,iz-xx)+cp2(ix,iz+xx))
     &			*gauss(xx,0,ka,atype3)
		  ep3(ix,iz)=ep3(ix,iz)+(cp3(x,iz)+cp3(ix+xx,iz)
     &			+cp3(ix,iz-xx)+cp3(ix,iz+xx))
     &			*gauss(xx,0,ka,atype3)
		enddo
		ep1(ix,iz)=ep1(ix,iz)+cp1(ix,iz)*gauss(0,0,ka,atype3)
		ep2(ix,iz)=ep2(ix,iz)+cp2(ix,iz)*gauss(0,0,ka,atype3)
		ep3(ix,iz)=ep3(ix,iz)+cp3(ix,iz)*gauss(0,0,ka,atype3)

c--------------------------------------------------------------
c�@�@�ϕ��͈͂��v�Z�̈悩�牡�������o�Ă���ꍇ�iflag2=1)
	  elseif(flag2(ix,iz).eq.1)then
		do x=ix-x_exp,ix+x_exp

		  cflag=0
		  if(x<0)then		
			x_temp=0
			cflag=1
		  elseif(x>nx)then
			x_temp=nx
			cflag=1	
		  endif
		  do ii=nrecess,1,-1	
			if(ix<lxr1(ii).and.iz>lzr1(ii).and.iz<lzr2(ii)
     &								.and.x>lxr1(ii))then
			  x_temp=lxr1(ii)-1
			  cflag=1	
			  exit
			elseif(ix>lxr2(ii).and.iz>lzr1(ii).and.iz<lzr2(ii)
     &								.and.x<lxr2(ii))then
			  x_temp=lxr2(ii)+1
			  cflag=1	
			  exit
			endif
		  enddo	

		  if(cflag.eq.0)then
			x_temp=x
		  endif

		  xx=ix-x	  

		  ep1(ix,iz)=ep1(ix,iz)+cp1(x_temp,iz)
     &			*gauss(xx,0,ka,atype3)
		  ep2(ix,iz)=ep2(ix,iz)+cp2(x_temp,iz)
     &			*gauss(xx,0,ka,atype3)
		  ep3(ix,iz)=ep3(ix,iz)+cp3(x_temp,iz)
     &			*gauss(xx,0,ka,atype3)	

		do z=iz-z_exp,iz-1
		  zz=iz-z

		  ep1(ix,iz)=ep1(ix,iz)+(cp1(x_temp,z)+cp1(x_temp,iz+zz))
     &		      *gauss(xx,zz,ka,atype3)  
		  ep2(ix,iz)=ep2(ix,iz)+(cp2(x_temp,z)+cp2(x_temp,iz+zz))
     &		      *gauss(xx,zz,ka,atype3)		   
		  ep3(ix,iz)=ep3(ix,iz)+(cp3(x_temp,z)+cp3(x_temp,iz+zz))
     &		      *gauss(xx,zz,ka,atype3)
		enddo
		enddo
c
c------------------------------------------------------------
c�@�@�ϕ��͈͂��v�Z�̈悩��c�������o�Ă���ꍇ�iflag2=2)
	  elseif(flag2(ix,iz).eq.2)then
		do z=iz-z_exp,iz+z_exp

		  cflag=0

		  if(z<0)then
			z_temp=0
			cflag=1	
		  elseif(z>nz)then
			z_temp=nz
			cflag=1
		  endif

		  do ii=nrecess,1,-1
			if(iz>lzr2(ii).and.ix>lxr1(ii).and.ix<lxr2(ii)
     &								.and.z<lzr2(ii))then
			  z_temp=lzr2(ii)+1
			  cflag=1
			  exit
			endif
		  enddo

		  if(cflag.eq.0)then
			z_temp=z
		  endif

		  zz=iz-z

		  ep1(ix,iz)=ep1(ix,iz)+cp1(ix,z_temp)
     &			*gauss(0,zz,ka,atype3)
		  ep2(ix,iz)=ep2(ix,iz)+cp2(ix,z_temp)
     &			*gauss(0,zz,ka,atype3)
		  ep3(ix,iz)=ep3(ix,iz)+cp3(ix,z_temp)
     &			*gauss(0,zz,ka,atype3)	

		do x=ix-x_exp,ix-1
		  xx=ix-x
		  ep1(ix,iz)=ep1(ix,iz)+(cp1(x,z_temp)+cp1(ix+xx,z_temp))
     &		      *gauss(xx,zz,ka,atype3)  
		  ep2(ix,iz)=ep2(ix,iz)+(cp2(x,z_temp)+cp2(ix+xx,z_temp))
     &		      *gauss(xx,zz,ka,atype3)		   
		  ep3(ix,iz)=ep3(ix,iz)+(cp3(x,z_temp)+cp3(ix+xx,z_temp))
     &		      *gauss(xx,zz,ka,atype3)
		enddo
		enddo

c----------------------------------------------------------------
c�@�@�ϕ��͈͂��v�Z�̈悩��c���������o�Ă���ꍇ�iflag2=3)
	  elseif(flag2(ix,iz).eq.3)then
c----x�����̐ϕ��J�n--------
		do x=ix-x_exp,ix+x_exp

c----�ϕ��͈͂��v�Z�̈���o���ꍇ�̏ꍇ����----
c	����:�ϕ��͈͂��v�Z�̈���o���ꍇ�́A�o���Ƃ���̕\�ʂ̊i�q�_��p����

c	cflag=0�Őϕ��͈͂��v�Z�̈��
c	cflag=1�Őϕ��͈͂��v�Z�̈�O
		  cflag=0
		  if(x<0)then		
			x_temp=0
			cflag=1
		  elseif(x>nx)then
			x_temp=nx
			cflag=1	
		  endif
		  do ii=nrecess,1,-1	
			if(ix<lxr1(ii).and.iz>lzr1(ii).and.iz<lzr2(ii)
     &								.and.x>lxr1(ii))then
			  x_temp=lxr1(ii)-1
			  cflag=1	
			  exit
			elseif(ix>lxr2(ii).and.iz>lzr1(ii).and.iz<lzr2(ii)
     &								.and.x<lxr2(ii))then
			  x_temp=lxr2(ii)+1
			  cflag=1	
			  exit
			endif
		  enddo	

		  if(cflag.eq.0)then
			x_temp=x
		  endif

		  xx=ix-x

c----z�����̐ϕ��J�n--------
		do z=iz-z_exp,iz+z_exp

		  cflag=0

		  if(z<0)then
			z_temp=0
			cflag=1	
		  elseif(z>nz)then
			z_temp=nz
			cflag=1
		  endif

		  do ii=nrecess,1,-1
			if(iz>lzr2(ii).and.ix>lxr1(ii).and.ix<lxr2(ii)
     &								.and.z<lzr2(ii))then
			  z_temp=lzr2(ii)+1
			  cflag=1
			  exit
			endif
		  enddo

		  if(cflag.eq.0)then
			z_temp=z
		  endif

		  zz=iz-z
c----�ϕ��v�Z����-------
		  ep1(ix,iz)=ep1(ix,iz)+ 
     &			cp1(x_temp,z_temp)*gauss(xx,zz,ka,atype3)
		  ep2(ix,iz)=ep2(ix,iz)+ 
     &			cp2(x_temp,z_temp)*gauss(xx,zz,ka,atype3)
		  ep3(ix,iz)=ep3(ix,iz)+ 
     &			cp3(x_temp,z_temp)*gauss(xx,zz,ka,atype3)
		enddo
		enddo
	  endif
		ep1(ix,iz)=ep1(ix,iz)/(alp*alp)
		ep2(ix,iz)=ep2(ix,iz)/(alp*alp)
		ep3(ix,iz)=ep3(ix,iz)/(alp*alp)

	  endif
	enddo
	enddo
c----��------

c	ep1�����J�Aep2��L�J�Aep3��X�J�̎����|�e���V�����̒l
	ep1 =ep1 /(2.0*pi)
	ep2 =ep2 /(2.0*pi)
	ep3 =ep3 /(2.0*pi)

	u_eff1=cp1
	u_eff2=cp2
	u_eff3=cp3

	do iz=0,lhet(nlayer-1)+ifix(a(iarea(nlayer))*5.0/dz+0.5)
	do ix=0,nx
	  if(atype(ix,iz).eq.1) then
		u_eff1(ix,iz)= cp1(ix,iz)
		u_eff2(ix,iz)= cp2(ix,iz)
		u_eff3(ix,iz)= cp3(ix,iz)
	  else
		u_eff1(ix,iz)= ep1(ix,iz)
		u_eff2(ix,iz)= ep2(ix,iz)
		u_eff3(ix,iz)= ep3(ix,iz)
	  endif
	enddo
	enddo

c------���t�l�X�U��--J.R. Watling���f��(2012�N�t����)-------------------------
c	���t�l�X�U��  border�Ŕ��˂��l�����邽�߁C�g���l�����ʕ���EP��CP�ɏ㏑������ sato
c	do iz= lhet(nchannel1)-5,lhet(nchannel1)	!�`���l����[���� 5��
c		do ix=0,nx	
c			u_eff1(ix,iz)= cp1(ix,iz)
c			u_eff2(ix,iz)= cp2(ix,iz)
c			u_eff3(ix,iz)= cp3(ix,iz)
c		enddo
c	enddo
c	
c	do iz= lhet(nchannel2)+1,lhet(nchannel2)+7		!�`���l�����[���� 7��
c		do ix=0,nx	
c			u_eff1(ix,iz)= cp1(ix,iz)
c			u_eff2(ix,iz)= cp2(ix,iz)
c			u_eff3(ix,iz)= cp3(ix,iz)
c		enddo
c	enddo
c----------------------------------------------------------------------------

!	ix=101
!	open(unit=11,file='ep.txt')
!	do iz=0,nz
!		write(11,*)iz,cp1(ix,iz), u_eff1(ix,iz)
!	enddo
!	close(11)
!	stop
!
!	iz=32
!	open(unit=11,file='ep.txt')
!	do ix=0,nx
!		write(11,*)ix,cp1(ix,iz), u_eff1(ix,iz)
!	enddo
!	close(11)
!	stop
!
!	open(unit=11,file='ep.txt',recl =1024)
!
!	do iz=0,nz
!		write(11,*)( -u_eff1(ix,iz), ix = 100, 130)
!	enddo
!
!	close(11)
!	stop
c	close(571)
	return
	end
