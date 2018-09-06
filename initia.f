	subroutine initia(
     &				spnum,dx,dz,xmax,zmax,de,dopem,twodeg, 
     &				lnpole1,lnpole2,lhet,iarea,bktq,mtemp,smh,pgm,
     &				p,kp,jpnum,af,n2,E2,i2max,	!120126homma
     &				cxr1,czr1,cxr2,czr2,
c     &				lxr1,lzr1,lxr2,lzr2)
     &				lxr1,lzr1,lxr2,lzr2,czpart2)
c
c === �ϐ���� ===
c	--- ���� ---
c	p(1-6,n) ... ���q���(1-3:k���W(kx,ky,kz)[m^-1],4:�U������[s],5-6:�ʒu(x,y))[m]
c	kp(1-3,n) ... �����J,
c	npmax ... �ő�L��������
c	jpnum ... �f�o�C�X���ɂ��钴���q�̐�
c	spnum ... �����q�̗��q��
c	bktq(T) ... �{���c�}���t�@�N�^�[
c	smh(i) ... i�Ԗڂ̒J�ɂ������2m*e / h
c	nvalley ... �J�̐�
c	d(i) ... �̈�i�ɂ�����h�i�[�s�����Z�x[m^-3]
c	dopem(j) ...j�Ԗڂ�mesh�̃h�[�v�Z�x
c	xmax,zmax ... �̈�̐��@[m]
c	dx,dz ... ���b�V���T�C�Y[m]
c	nx,nz ... �i�q�_�̕����̐�
c	
c	--- ���̑� ---
c	npij ... ���郁�b�V���ɂ����钴���q�̐�
c	e ... �^���G�l���M�[[eV]
c	ak ... �g���x�N�g���̑傫��[m]
c	cf,sf ... �V���p�Ƃɂ�����cos,sin
c	fi ... ���ʊp��
c	i ... �̈�̃J�E���g
c	ix,iz ... ����̈�ɂ�����i�q�_�̃J�E���g
c	n ... �f�o�C�X���ɑ��݂��钴���q���Ɋւ���J�E���g
c
      implicit none
	include 'arraysize.fi'
	real	pi,q
	parameter(pi = 3.1415927,	q = 1.60219e-19)
c
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c
	integer(2) lhet(nlayer)
	integer(1) iarea(nlayer)
	real	twodeg(nlayer)				!2DEG�V�[�g�d�ח�(*1.0e-4 cm^-3)
	integer(2),dimension (npole)	:: lnpole1,lnpole2
	real, dimension (nrecess) :: cxr1,czr1,cxr2,czr2
	integer(2),dimension (nrecess) :: lxr1,lzr1,lxr2,lzr2
	real	dx,dz,xmax,zmax,spnum
	real,	dimension (nenergy)	:: de
	real,	dimension (0:nx,0:nz)	:: 	dopem
	real,	dimension (ntenum)	:: 	bktq
	integer(2),	dimension (0:nx,0:nz)	:: 	mtemp
	real,	dimension (nvalley,narea)	:: smh,af	!120126homma
c	real,	dimension (nvalley,nenergy,narea)	:: pgm
	real,	dimension (nvalley,nenergy,npart)	:: pgm	!07/8/4 �s�����U��
	real,	dimension (6,npmax)	:: p
	integer(1),dimension (3,npmax)	:: kp
	integer jpnum
c
	integer n,npij,m,i,ix,iz,ia,ii	
	real	apij,e,ak,cf,sf,fi
	real	rnd
	integer(2) lz(0:nlayer)
	integer ipart,kpart		!07/8/4 �s�����U��
	real czpart2(npart)		!07/8/4 �s�����U��
	integer ni,i2max							!120126homma
	real ntab,n2(0:i2max),E2(0:i2max)			!120126homma
c
c
c ---- ( ���q�̋�ԕ��z�̐ݒ� ) ----
c

c	open(unit=561,file='test2.bin')

	lz(0)=0
	lz(1:nlayer)=lhet(1:nlayer)
	n=0
	do 30 i=1,nlayer
	do 20 iz=lz(i-1),lz(i)
	do 20 ix=0,nx

	  apij=dopem(ix,iz)*dx*dz/spnum	!���b�V�������蒴���q���i�����j
c	write(*,*) ix,iz,apij
c	  if((ix.ge.lnpole1(2)).and.(ix.le.lnpole2(2)).and.
c     &					(lnpole1(2).ne.lnpole2(2)))	apij=0.0	!�Q�[�g�d�ɉ�
	  if((ix.eq.0).or.(ix.eq.nx)) apij=apij/2.0
	  if((iz.eq.lz(i-1)).or.(iz.eq.lz(i))) apij=apij/2.0

c	-----(���Z�X)----
	  do ii=1,nrecess	!nrecess=1,2
		if(((ix.eq.lxr1(ii)).or.(ix.eq.lxr2(ii)))
     &		 .and.(iz.le.lzr2(ii)))	apij=apij/2
	  enddo

	  npij = ifix(apij+rnd())		!���b�V�������蒴���q���i�����j
	  if(npij.eq.0) goto 20
	  do 10 m=1,npij
	    n=n+1
c	write(*,*) n, npmax
	    if(n.gt.npmax)then
		  write(99,*) '���q���I�[�o�[(initia)'
		  write(*,*) '���q���I�[�o�[(initia)'
		  stop
	    endif
c
c	    do											!120126homma
c	      e = -bktq(mtemp(ix,iz))*alog(rnd())*1.5	!120126homma
c	      if(e.gt.0.0) exit							!120126homma
c	    enddo										!120126homma
c
	    kp(1,n) = 1		!kp(1,n):���q�����ݍ݂�ꏊ�̒JNo.
		kp(2,n) = 1		!kp(2,n):���q�����ݍ݂�ꏊ�̃G�l���M�[No.
		kp(3,n) = i		!kp(3,n):���q�����ݍ݂�ꏊ�̑f��No.
		ia = iarea(i)
c
c----		�A�e�w���Ƃɕs�����Z�x�����߂�------------	!07/8/4 �s�����U��
		do ipart=1,npart
c			�d�ɂ̂Ƃ��͂������΂��H�I(�v�C���I)
			if((ipart.eq.npart-2).or.(ipart.eq.npart-1))cycle 	
			if(p(6,n).le.czpart2(ipart))then
				kpart=ipart 
				exit
			else
			 	kpart=npart
				if(p(6,n).gt.czpart2(npart))then
					write(*,*)'initia�ŕs�����Z�x�G���[1'
					stop
				endif
			endif
		enddo
c----(�t�F���~�E�f�B���b�N���z�ɏ]�������G�l���M�[���z)---- 120126homma
		ntab=0.0d0
		do ni=0,i2max
			ntab=ntab+n2(ni)
			if(rnd().le.ntab)exit
		enddo
		if(ni.gt.i2max)then
			ni=i2max
		endif
		e=E2(ni)
		ak = smh(kp(1,n),ia)*sqrt(e*(1+af(kp(1,n),ia)*e))
c----------------------------------------------------------
c		ak = smh(kp(1,n),ia)*sqrt(e)	!120126homma
	    cf = 1.0-2*rnd()		!cos��
		sf = sqrt(1.0-cf*cf)	!sin��
	    fi = 2.0*pi*rnd()		!��
c
		do while (e.gt.(de(kp(2,n))*nemax)*0.9)
			if(kp(2,n).ge.nenergy) exit
			kp(2,n)	= kp(2,n)+1
		enddo
c
	    p(1,n) = ak*sf*cos(fi)	!kx = ak*sin��*cos��
	    p(2,n) = ak*sf*sin(fi)	!ky = ak*sin��*sin��
	    p(3,n) = ak*cf			!kz = ak*cos��
c	    p(4,n) = -alog(rnd())*pgm(kp(1,n),kp(2,n),ia)
	    p(4,n) = -alog(rnd())*pgm(kp(1,n),kp(2,n),kpart)
	    p(5,n) = dx*(rnd()+float(ix)-0.5)
	    p(6,n) = dz*(rnd()+float(iz)-0.5)
	    if(ix.eq.0) p(5,n) = dx*rnd()*0.5
	    if(ix.eq.nx) p(5,n) = xmax-dx*rnd()*0.5
	    if(iz.eq.0) p(6,n) = dz*rnd()*0.5			!08/3/31 �|��
	    if(iz.eq.nx) p(6,n) = zmax-dz*rnd()*0.5		!08/3/31 �|��
c
c		----(���Z�X)-----
		do ii=1,nrecess
			if((iz.ge.lzr1(ii)).and.(iz.le.lzr2(ii)).and.(ix.eq.lxr1(ii)))
     &			p(5,n) = cxr1(ii)-dx*rnd()*0.5
    			if((iz.ge.lzr1(ii)).and.(iz.le.lzr2(ii)).and.(ix.eq.lxr2(ii)))
     &			p(5,n) = cxr2(ii)+dx*rnd()*0.5
		enddo
	    if(iz.eq.lz(i-1)) p(6,n) = lz(i-1)*dz+dz*rnd()*0.5
	    if(iz.eq.lz(i  )) p(6,n) = lz(i  )*dz-dz*rnd()*0.5
c	    if(iz.eq.nz) p(6,n) = zmax-dz*rnd()*0.5
	
		do ii=1,nrecess
		if((p(5,n).ge.cxr1(ii)).and.(p(5,n).lt.cxr2(ii)))then
		if((p(6,n).ge.czr1(ii)).and.(p(6,n).lt.czr2(ii)))then
		   kp(1,n) = 0		!kp(1,n):���q�����ݍ݂�ꏊ�̒JNo.
		endif
		endif
		enddo

c	write(561,*) n,i,p(5,n),p(6,n)

   10	continue
   20	continue
   30	continue
c
c ---- ( �w�e���E�ʂɂ���ɗ��q�̋�ԕ��z�̐ݒ� ) ----

	do 35 i=1,nlayer-1
c	do 35 i=1,nlayer-1 !11/04/08�C���K�v�ӏ��H
	if(twodeg(i).eq.0.0)cycle
	do 25 ix=0,nx
c	if((ix.ge.lnpole1(2)).and.(ix.le.lnpole2(2)))cycle		! �Q�[�g�d�ɉ�
	iz=lhet(i)
	apij= twodeg(i)*dx/spnum
	if((ix.eq.0).or.(ix.eq.nx)) apij=apij/2.0
	npij = ifix(apij+rnd())
	do 15 m=1,npij
		n=n+1
		if(n.gt.npmax)then
			write(99,*) '���q���I�[�o�[(initia)'
			write(*,*) '���q���I�[�o�[(initia)'
			stop
		endif
c
		kp(1,n) = 1
		kp(2,n) = 1
		kp(3,n) = i+1
		ia = iarea(i)
c
c----		�A�e�w���Ƃɕs�����Z�x�����߂�------------	!07/8/4 �s�����U��
		do ipart=1,npart
c			�d�ɂ̂Ƃ��͂������΂��H�I(�v�C���I)
			if((ipart.eq.npart-2).or.(ipart.eq.npart-1))cycle 	
			if(p(6,n).le.czpart2(ipart))then
				kpart=ipart 
				exit
			else
			 	kpart=npart
				if(p(6,n).gt.czpart2(npart))then
					write(*,*)'initia�ŕs�����Z�x�G���[1'
					stop
				endif
			endif
		enddo
c----(�t�F���~�E�f�B���b�N���z�ɏ]�������G�l���M�[���z)---- 120126homma
		ntab=0.0d0
		do ni=0,i2max
			ntab=ntab+n2(ni)
			if(rnd().le.ntab)exit
		enddo
		if(ni.gt.i2max)then
			ni=i2max
		endif
		e=E2(ni)
		ak = smh(kp(1,n),ia)*sqrt(e*(1+af(kp(1,n),ia)*e))
c----------------------------------------------------------
	    cf = 1.0-2*rnd()		!cos��
		sf = sqrt(1.0-cf*cf)	!sin��
	    fi = 2.0*pi*rnd()		!��
c
		do while (e.gt.(de(kp(2,n))*nemax)*0.9)
			if(kp(2,n).ge.nenergy) exit
			kp(2,n)	= kp(2,n)+1
		enddo
c
	    p(1,n) = ak*sf*cos(fi)	!kx = ak*sin��*cos��
	    p(2,n) = ak*sf*sin(fi)	!ky = ak*sin��*sin��
	    p(3,n) = ak*cf			!kz = ak*cos��
c	    p(4,n) = -alog(rnd())*pgm(kp(1,n),kp(2,n),ia)
	    p(4,n) = -alog(rnd())*pgm(kp(1,n),kp(2,n),kpart)
	    p(5,n) = dx*(rnd()+float(ix)-0.5)
	    p(6,n) = (lhet(i)*dz)+dz*rnd()*0.5

	    if(ix.eq.0) p(5,n) = dx*rnd()*0.5
	    if(ix.eq.nx) p(5,n) = xmax-dx*rnd()*0.5

   15	continue
   25	continue
   35	continue
c
c
	jpnum = n
c	write(*,*) jpnum
	kp(1,jpnum+1:npmax) = 0

c	close(561)
	return
	end
