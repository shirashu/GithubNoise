
	subroutine poisso(dx, dz, lhet, twodeg,iarea,
     &				vb, vi, lnpole1, lnpole2, melpos,
     &				eps, eg, c_ratio,
     &				u, cn, cp, dopem, maceps,
     &				lxr1,lzr1,lxr2,lzr2,dltec)
      implicit none
c
c===�����萔===
	real	qe
	parameter(qe=1.60219e-19)
c
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c
c	nx,nz	... �i�q�̐�
c	npole   	... �d�ɂ̐�
c	nlayer   	... �G���A��
c	nvalley		... �J�̐�
c	lxr1,lzr1	...���Z�X���[�E����[
c	lxr1,lzr1	...���Z�X�E�[�E�E���[
c
c---�f�o�C�X�\��---
	real	dx,dz		!���b�V���T�C�Y[m]
	integer(2),	dimension (nrecess)	:: lxr1,lzr1,lxr2,lzr2
	integer(2)	lhet(nlayer)	!�w�e���E�ʂ̐[��[m]
	integer(1)	iarea(nlayer)
c---2DEG�d�ח�---
	real	twodeg(nlayer)		!2DEG�V�[�g�d�ח�(*1.0e-4 cm^-3)
c---�d��---
	real,	dimension(npole)	:: vb,vi			!vb:�r���h�C���d��,vi:�d��
	integer(2),dimension(npole)	:: lnpole1,lnpole2	!�d�ɂ̕�
	integer(1),dimension(npole)	:: melpos		!
c---�̈�ʃp�����[�^---
	real,	dimension(narea)	:: eps		!�U�d��
	real,	dimension(nvalley,narea):: eg
	real	c_ratio
c---�f�o�C�X�����---
	real,	dimension(0:nx,0:nz)	:: u
	real,	dimension((nx+1)*(nz+1))	:: cn,cp,dopem
c	dopem(j) ...j�Ԗڂ̗̈�̕s�����Z�x
c	
c===poisso�p�ϐ�==========
	real	maceps
c
c====���[�J���ϐ�====
c====poisso�p�ϐ�==========
	integer,save	::	nb1,nb2
	integer,save	::	ierr,m1,nxnz
	integer,save,allocatable	::	ib1(:),ib2(:)
c
	real,save,allocatable	::	bv(:),alp(:),bet(:),p(:),f(:),q(:)

	integer	i,ii,ia,il,ix,iz,intyp
	integer	nxnz2
	integer	lh(0:nlayer)
	real,save	:: vidt(3)
	real	altwodeg
	real,save	:: c_ratio_old
c
c---( �\�ʓd�� )---	
	real,dimension(0:nrecess) :: vsfr 
	real vsf,vsb
	parameter(vsf=0.4, vsb=0.8)
c---( ���ʓd�� )--- 11/08/06��
	real vsbr
	real, dimension	(nvalley,narea):: dltec
c------------------------------------------
c=========================================================================================================
c
	nxnz=(nx+1)*(nz+1)	!�S���b�V����
	m1=nx+1				!nx+1
c
c---( poisson�������̒萔 )----
	if (.not. allocated(p)) then	!�ŏ��������s
		allocate (p(nxnz))
		lh(0)=0;lh(1:nlayer)=lhet(1:nlayer)
		do il = 1, nlayer
			ia = iarea(il)
			p(lh(il-1)*m1+1:lh(il)*m1)	= eps(ia)	!�U�d����p�Ɋi�[
		enddo
		ia = iarea(nlayer)
		p(lh(nlayer)*m1+1:nxnz)	= eps(ia)
		allocate (f(nxnz))	!(0:nx,0:nz))
		allocate (q(nxnz))
		q=0.0
	endif
c
	f=(dopem-cn+cp)*qe		!f,dopem,cn:�S�̔z��
c
c----------------------------------------------------------
c
c=======( ���E����{�������g��} )======================
	if (.not. allocated(bv)) then	!�ŏ��������s
		nxnz2 = (nx+nz)*2	!!!!!!!!!!
		allocate (bv(nxnz2),alp(nxnz2),bet(nxnz2))
		allocate (ib1(nxnz2),ib2(nxnz2))
		vidt=huge(vidt)		!vidt:�S�̔z��
		bv=0;alp=0;bet=0	!������
	endif
c
c	���E�����̐ݒ�E�Đݒ�
c	�i�d���ɕω�������������s�j
c	if(    (vi(1).ne.vidt(1))
c     &   .or.(vi(2).ne.vidt(2))
c     &   .or.(vi(3).ne.vidt(3))
c     &   .or.(c_ratio.ne.c_ratio_old)
c     &   )then
	if(npole.ge.5)then
		vb(4) = vb(1)
		vi(4) = vi(1)
		vb(5) = vb(3)
		vi(5) = vi(3)
	endif
	nb1 = 0;nb2 = 0
c
c---( �\�ʓd�� )--- 2DEG���Z�o 2011/3/25��	
c	Al025 Vg0.0V		9.6263360E+11 cm-2
c	InAs Vg0.0V			6.0841971E+11 cm-2
c
	vsfr(0) = 24.0e15*qe
	vsfr(1) = 24.0e15*qe
	vsfr(2) = 24.0e15*qe
c
c---( ���ʓd�� )--- 2DEG���Z�o 2011/8/6��
c	Al025 Vg0.0V		2.9901177E+11 cm-2
c	InAs Vg0.0V			4.6512364E+11 cm-2
c	InAs AlSb			8.2250168E+22 cm-2	
c
	vsbr	= 1.0e15*qe
c
c---( ���E�����̐ݒ� )---
	do 10 i=1,nxnz			!�S���b�V�����܂�
		ix = (mod(i-1,m1))	!���܂蕪�A���Wx  (m1=nx+1)
		iz = ((i-1)/m1)		!�������������A���Wz
c
c	---�f�o�C�X����[�E�E��[---
		if((i.eq.1).or.(i.eq.m1))then
			if(i.eq. 1)	intyp = 2	!�d�ɍ��[
			if(i.eq.m1)	intyp = 3	!�d�ɉE�[
			do ii=1,npole
				if(lnpole1(ii).ne.lnpole2(ii))then	!���[�E�[�������Ȃ��Ƃ�
				!!-------���Z�X������Γ���Ȃ��@2017/05/22���� ��
					if((melpos(ii).eq.1).and.		!�Q�[�g�̂Ƃ�
     &				  (ix.ge.lnpole1(ii)).and.(ix.le.lnpole2(ii)))then
						goto 1000				!�f�B���N�����E��
				!!-------�@��
					elseif((melpos(ii).eq.intyp).and.			
     &				  (iz.ge.lnpole1(ii)).and.(iz.le.lnpole2(ii)))then
						goto 1000				!�f�B���N�����E��
					endif	
				endif
			enddo
			goto 3000
c
c		elseif	(iz.eq.0)then						!��[
c			do ii=1,npole
c				if((melpos(ii).eq.1).and.(lnpole1(ii).ne.lnpole2(ii)).and.
c    &				(ix.ge.lnpole1(ii)).and.(ix.le.lnpole2(ii)))then
c					goto 1000
c				endif
c			enddo
c			goto 3000
c
c	---�f�o�C�X��[---
		elseif (iz.eq.0)then
			if ((ix.le.minval(lxr1)).or.(ix.ge.maxval(lxr2)))then
				do ii=1,npole
			if((melpos(ii).eq.1).and.(lnpole1(ii).ne.lnpole2(ii)).and.
     &				(ix.ge.lnpole1(ii)).and.(ix.le.lnpole2(ii)))then
						goto 1000
					endif
				enddo
			elseif((ix.gt.minval(lxr1)).and.(ix.lt.maxval(lxr2)))then
				go to 10
			endif 
			goto 3000
c
c	---���Z�X�����[---
		elseif ((iz.gt.lzr1(1)).and.(iz.lt.lzr2(1))
     &						   .and.(ix.eq.lxr1(1)))then
			go to 2000
		elseif ((iz.gt.lzr1(2)).and.(iz.lt.lzr2(2))
     &						   .and.(ix.eq.lxr1(2)))then !2�i�ڃA�E�g
c		write(*,*) 'err'
c-----06/11/9 ��[�@�Q�[�g�ƃ��Z�X�̐ڐG�̂��ߒǉ� ---------
			if (lnpole1(2).eq.lxr1(2))then
				ii=2	!�v���Ӂ@�Q�[�g�d�ɂ��w��
				goto 1000
			endif
c------------------------------------------------------------
			go to 2000
c
c	---���Z�X���E�[---
		elseif ((iz.gt.lzr1(1)).and.(iz.lt.lzr2(1))
     &						   .and.(ix.eq.lxr2(1)))then
			go to 2000
		elseif ((iz.gt.lzr1(2)).and.(iz.lt.lzr2(2))
     &						   .and.(ix.eq.lxr2(2)))then !2�i�ڃA�E�g
c		write(*,*) 'err'
c-----06/11/9 ��[�@�Q�[�g�ƃ��Z�X�̐ڐG�̂��ߒǉ� -----------
			if (lnpole2(2).eq.lxr2(2))then
				ii=2	!�v���Ӂ@�Q�[�g�d�ɂ��w��
				goto 1000
			endif
c-------------------------------------------------------------
			go to 2000
c
c	---���Z�X����[(��i��)---
		elseif ((iz.eq.lzr2(1)).and.(ix.ge.lxr1(1))
     &						   .and.(ix.le.lxr1(2)))then
c-----06/11/9 ��[�@�Q�[�g�ƃ��Z�X�̐ڐG�̂��ߒǉ� -----------2�i�ڃA�E�g
			if ((lnpole1(2).eq.lxr1(2)).and.(ix.eq.lnpole1(2)))then
				ii=2	!�v���Ӂ@�Q�[�g�d�ɂ��w��
				goto 1000
			endif
c-------------------------------------------------------------			
			go to 5000
		elseif ((iz.eq.lzr2(1)).and.(ix.ge.lxr2(2))
     &						   .and.(ix.le.lxr2(1)))then
c-----06/11/9 ��[�@�Q�[�g�ƃ��Z�X�̐ڐG�̂��ߒǉ� -----------2�i�ڃA�E�g
			if ((lnpole2(2).eq.lxr2(2)).and.(ix.eq.lnpole2(2)))then
				ii=2	!�v���Ӂ@�Q�[�g�d�ɂ��w��
				goto 1000
			endif
c-------------------------------------------------------------
			go to 5000
c
c	---���Z�X����[(��i��)---
ccc		elseif ((iz.eq.lzr2(1)).and.(ix.ge.lxr1(1)).and.(ix.le.lxr2(1)))then !���Ƃ��ƂȂ�
		elseif ((iz.eq.lzr2(2)).and.(ix.ge.lxr1(2))
     &						   .and.(ix.le.lxr2(2)))then !2�i�ڃA�E�g
			do ii=1,npole
				if((melpos(ii).eq.1).and.(lnpole1(ii).ne.lnpole2(ii)).and.	
     &				(ix.ge.lnpole1(ii)).and.(ix.le.lnpole2(ii)))then
					goto 1000
				endif
			enddo
			go to 6000 
c
c	---�f�o�C�X���E�[---
		elseif	((ix.eq.0).or.(ix.eq.nx))then
			if(ix.eq. 0)	intyp = 2
			if(ix.eq.nx)	intyp = 3
			do ii=1,npole
				if((melpos(ii).eq.intyp).and.
     &				(iz.ge.lnpole1(ii)).and.(iz.le.lnpole2(ii)))then
					goto 1000
				endif
			enddo
			goto 2000
c
c	---�f�o�C�X���[---
		elseif	(iz.eq.nz)then
			goto 4000
c
c	---�f�o�C�X���g---
		else
			goto 10
		endif
c
c	---���E�����K�p---
c	---Dirichlet�̋��E����---
 1000		continue
		nb1=nb1+1
		ib1(nb1)=i
		bv(nb1)=-vb(ii)+vi(ii)
c
c	---HEMT�p���E---
		if((ix.eq.0).or.(ix.eq.nx))then		!HEMT�̍��[�܂��͉E�[
			if(iz.le.lhet(1))then			!cap�w�Ȃ�
				ia = iarea(1)
c				bv(nb1)=bv(nb1)+0.7*(eg(1,ia)-eg(1,2)) !�o���h�M���b�v�H�ύX�K�v
				bv(nb1)=bv(nb1)-dltec(1,ia)		!2017/12/1 ���
c		write(*,*) iz,ia, bv(nb1)
c	pause
			else
				do il = 2, nlayer
					if((lhet(il-1).lt.iz).and.(iz.le.lhet(il)))then
						ia = iarea(il)
						bv(nb1)=bv(nb1)-dltec(1,ia)
c						if((ia.eq.4))then
c						  bv(nb1)=bv(nb1)+0.42*(eg(1,ia)-eg(1,2))
c						else 
c						  bv(nb1)=bv(nb1)+0.7*(eg(1,ia)-eg(1,2)) !11/06/29�� 
c						endif
c		write(*,*) iz,ia, bv(nb1)
c		pause
					exit
					endif
				enddo
			endif
		endif
		goto 10
c
c	---���R���E����---
 2000		continue
		nb2=nb2+1
		ib2(nb2)=i
		alp(nb2)=0.0
		bet(nb2)=0.0
		goto 10
c
c	---InGaAs�\�ʓd��(��[)---
 3000		continue
		nb2=nb2+1
		ib2(nb2)=i
		alp(nb2)=0.0
		bet(nb2)= vsfr(0)
c		bet(nb2)= 6.20e16*qe
c		bet(nb2)=sqrt(2.0*p(i)*vsf*qe*1.0e+25)
		goto 10
c
c	---InAlAs�\�ʓd��(���Z�X����i��)---
5000		continue
		nb2=nb2+1
		ib2(nb2)=i
		alp(nb2)=0.0
		bet(nb2)= vsfr(1)
c		bet(nb2)= 3.50e16*qe
		goto 10
c
c	---InAlAs�\�ʓd��(���Z�X����i��)---
6000		continue
		nb2=nb2+1
		ib2(nb2)=i
		alp(nb2)=0.0
		bet(nb2)= vsfr(2)
c		bet(nb2)= 4.00e16*qe
		goto 10
c
c	---InAlAs���ʓd��(���[)---
 4000		continue
		nb2=nb2+1
		ib2(nb2)=i
		alp(nb2)=0.0
c		bet(nb2)=sqrt(2.0*p(i)*vsb*qe*1.0e+21)
		bet(nb2) = vsbr
		goto 10
c
c	---���E������K�p���Ȃ�---
   10	continue
	vidt(1:3)=vi(1:3)
	c_ratio_old = c_ratio 
c
c===================================================
c
c---( poisson��������solver��call )----
	call poison(nxnz, m1, dx, dz, 
     &			p, q, f, nb1, ib1, bv, nb2, ib2, alp, bet,
     &			u, ierr, maceps, lxr1, lzr1, lxr2, lzr2)

c	do ix=0,nx
c      write(160,1111)(-u(ix,iz),iz=0,nz)
c      end do 
c1111   format (E15.7,400(',',E15.7))	
	return
	end