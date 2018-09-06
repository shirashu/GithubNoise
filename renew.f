	subroutine renew(
     &                 jpnum,ncon,
     &                 bktq,mtemp,dx,dz,lnpole1,
     &                 lnpole2,melpos,jspsum,vb,
     &                 pgm,smh,p,kp,lhet,iarea,twodeg,balis_flag,
     &				 hhm,af,af4,ecr,de,			!120126homma
     &				 i2max,E2,n2,i3max,E3,n3,	!08/8/6 �|��
     &				x_start,z_start)		!120817sato

	implicit none
c
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c
	integer	jpnum,ncon
c
	real	bktq(ntenum),dx,dz
	integer(2)	mtemp(0:nx,0:nz)
	integer(2)	lnpole1(npole),lnpole2(npole)
	integer(1)	melpos(npole)
	integer(4)	jspsum(npole)
	real	vb(npole)
c
c	real,	dimension (nvalley,nenergy,narea)		:: pgm
	real,	dimension (nvalley,nenergy,npart)		:: pgm	!07/8/4 �s�����U��
	real,	dimension (nvalley,narea)	:: smh
	real	p(6,npmax)
	integer(1)	kp(3,npmax)
	integer(2)	lhet(nlayer)
	integer(1)	iarea(nlayer)
	real	twodeg(nlayer)		!2DEG�V�[�g�d�ח�(*1.0e-4 cm^-3)
c
	integer	jn
	integer	n,nmax
	integer i,j,ni
	integer	intyp,nedge,jpole
	real	hdx,hdz
	character(80) form
	integer(2)	lh
	integer	jcreat
c	integer	jpnum2,jpnum3
	real,allocatable	::	bktqn(:)
	integer,allocatable	::	jnn(:)
c
	integer,allocatable	::	ipn(:,:)
c
c-----�k�ތ���-----
	integer	i2max,i3max
	real E2(0:i2max),n2(0:i2max)
	real E3(0:i3max),n3(0:i3max)
c------------------
c
c----�o���X�e�B�b�N�̌v�Z----------
	integer(1),dimension (0:npmax) ::	balis_flag	!07/11/22 ��[
c-----------------------------------
	real	ecr(7,int(nemax/4))							!120126homma
	real,	dimension (nenergy)	:: de					!120126homma
	real,	dimension (nvalley,narea)	:: hhm,af,af4	!120126homma


c----pass---���ώ��R�s��120817sato
	real	x_start(npmax)
	real	z_start(npmax)
c
c	lh = lhet(maxloc(twodeg,1))	!(�v�C��)
	hdx=dx/2.	!dx ���b�V���̔���
	hdz=dz/2.	!dz ���b�V���̔���
c
	allocate(ipn(0:nx+nz,4))
	ipn = 0		!Fortran90
c
c	--- �d�ׂ̒��������ɍ��킹�ė��q���������� ---
c
c	!$OMP PARALLEL DO SHARED(jspsum,ipn,hdx) PRIVATE(jn,intyp)
	do 40 n=1,jpnum
c		jn = 0

c	---	���q�̑Ή�����d�ɂ̈ʒu�I��  ---
		if(p(6,n).lt.hdz) then						!��[
			intyp = 1
		elseif(p(5,n).lt.hdx) then					!���[
			intyp = 2
		elseif(p(5,n).gt.(dx*float(nx)-hdx)) then	!�E�[
			intyp = 3
		else
	cycle
		endif

c	
c	--- �d�ɋߖT�̃��b�V����I�� ---
		select case (intyp)
		case (1)
			jn=min(nint(p(5,n)/dx),nx)
		case (2,3)
			jn=min(nint(p(6,n)/dz),nz)
c			if((jn.gt.lh-2).and.(jn.lt.lh+2))cycle	! 2DEG���ӂ͏��O
		case default
c			jn=0
			cycle				
		end select
c
c	--- �d�ɂ�I�� ---
		jpole = 0
c
	!SMP$	ASSERT(ITERCNT(1))
		do i=1,npole
			if(vb(i).eq.0.0)then
				if((melpos(i).eq.intyp).and.
     &				(jn.ge.lnpole1(i)).and.(jn.le.lnpole2(i)))then
					if(ipn(jn,intyp).gt.ncon)then
						kp(1,n) = 0
c
c	!$OMP						ATOMIC
						jspsum(i)=jspsum(i)+1
					else
c
c	!$OMP						ATOMIC
						ipn(jn,intyp)=ipn(jn,intyp)+1
					endif
					exit
				endif
			endif
		enddo
   40	continue
c
	nmax = jpnum
c	jpnum2 = jpnum
c	jpnum3 = jpnum
c	do n=1,nmax
c		if(kp(1,n).eq.0)jpnum2=jpnum2-1
c	enddo
	do n=1,nmax
c	--- �������痱�q���������� ---
		do while (kp(1,jpnum).eq.0)
			jpnum=jpnum-1
		enddo
c
		if(jpnum.le.n)	exit
c	--- n�Ԗڂ̗��q�������i�������痱�q���ړ��j ---
c		kp=0�ŏ�����n�Ԗڂ̗��q�̔z��Ɉ�Ԍ��jpnum�̗��q�̃f�[�^��n�Ԗڂ̗��q�̔z��ɓ����D120817sato
		if(kp(1,n).eq.0)then
			p (1:6,n) = p (1:6,jpnum)
			kp(1:3,n) = kp(1:3,jpnum)
			balis_flag(n) = balis_flag(jpnum)	!07/11/22 ��[
			x_start(n) = x_start(jpnum)		!120817sato	
			z_start(n) = z_start(jpnum)		!120817sato	
c			kp(1,jpnum) = 0
			jpnum=jpnum-1
		endif
	enddo
c	do n=1,nmax
c		if(kp(1,n).eq.0)jpnum3=jpnum3-1
c	enddo
c
c	!$OMP BARRIER
c--- ( �d�ɂɂ����闱�q�̐��� ) ---
c	
	allocate(bktqn(ncon*max(nx,nz)),jnn(ncon*max(nx,nz)))	!�v�C��
c	--- �d�ɂ�I�� ---
	do i=1,npole
		if(vb(i).ne.0.0)cycle			!�I�[�~�b�N�d�ɈȊO���O
		if(lnpole1(i).eq.lnpole2(i))cycle		!�d�ɒ��O�Ȃ珜�O
		intyp = melpos(i)				!�d�Ɉʒu
		select case(intyp)
		case(1)	;nedge = nx				!�������d��
		case(2,3);nedge = nz				!�c�����d��
		end select
		n = 0
		do jn=lnpole1(i),lnpole2(i)
c			if(((intyp.eq.2).or.(intyp.eq.3)).and.
c     &			(jn.gt.lh-2).and.(jn.lt.lh+2))cycle	! 2DEG�t�߂͏��O
			if((jn.eq.0).or.(jn.eq.nedge))then		!�f�o�C�X�p
c				ni=(ncon)/2-ipn(jn,intyp)
				cycle
			else
				ni=(ncon)-ipn(jn,intyp)
			endif
c
			if(ni.gt.0)then
				do j=1,ni
					select case(intyp)
					case(1);bktqn(n+j) = bktq(mtemp(jn, 1))	!�d�Ɉʒu ��
					case(2);bktqn(n+j) = bktq(mtemp(1 ,jn))	!��
					case(3);bktqn(n+j) = bktq(mtemp(nx,jn))	!�E
					end select
				enddo
				jnn((n+1):(n+ni))=jn
				n=n+ni
			endif
		enddo
		if(n.le.0)cycle
		if((jpnum+n).gt.npmax) then
			form = "('���q���� npmax=',I12,'�𒴂��܂���(creat)')"
			write(99,form)npmax
			write(*,form)npmax
			stop
		endif
		jcreat = n
c
		call creat2(smh,pgm,lhet,iarea,dx,dz,
     &				intyp,p,kp,jpnum,jcreat,bktqn,jnn,
     &				hhm,af,af4,ecr,de,			!120126homma
     &				i2max,E2,n2,i3max,E3,n3,	!08/8/6 �|��
     &				x_start,z_start)		!120817sato
c
		jspsum(i)=jspsum(i)-jcreat
		jpnum=jpnum+n
	enddo
	deallocate(bktqn,jnn)
	deallocate(ipn)
	return
	end subroutine renew
cc
c==== ( �������v���O���� ) ====
c
c
c--	���q�쐬�T�u���[�`��  ---
	subroutine creat2(smh,pgm,lhet,iarea,dx,dz,
     &				intyp,p,kp,jpnum,jcreat,bktqn,jnn,
     &				hhm,af,af4,ecr,de,				!120126homma
     &				i2max,E2,n2,i3max,E3,n3,		!�|��
     &				x_start,z_start)				!120817sato


	implicit none
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c
	real	smh(nvalley,narea)
c	real	pgm(nvalley,nenergy,narea)
	real	pgm(nvalley,nenergy,npart)	!07/8/4 �s�����U��
	integer(2)	lhet(nlayer)
	integer(1)	iarea(nlayer)
	real	dx,dz
	real	p(6,npmax)
	integer(1)	kp(3,npmax)
	integer	intyp
	integer	jpnum
	integer	jcreat
	real	bktqn(jcreat)
	integer	jnn(jcreat)
c
	real	rnd
	real	akh,akt,fi
c	integer(1)	kv,ken,kl,ka
	integer(1)	kv,ken,kl,ka,kpart	!07/8/4 �s�����U��
	integer(1)	hetero
	integer	n,ni,jn
c	integer	inp
c
c-----�k�ތ���-----
	integer	i,i2max,i3max
	real r1,ntab,e
	real E2(0:i2max),n2(0:i2max)
	real E3(0:i3max),n3(0:i3max)
c------------------
	integer ie2											!120126homma
	real	cf,sf,iz,sk,sq,ecr(7,int(nemax/4))			!120126homma
	real,	dimension (nenergy)	:: de					!120126homma
	real,	dimension (nvalley,narea)	:: hhm,af,af4	!120126homma


c----pass---���ώ��R�s��120817sato
	real	x_start(npmax)
	real	z_start(npmax)

c
	real	pi
	parameter (pi  = 3.14159265358979323846)
c
	kv=1
	ken=1
	ka=1	
	kl=1	!�@
	kpart=1 !kl=1�ł��邱�Ƃ��炱���1�Ƃ��� !07/8/4 �s�����U��
	hetero = -1

	do ni = 1,jcreat				!jcreat
	jn = jnn(ni)
	select case (intyp)
c      case (1);ka(ni) = 1
      case (2,3)
		do while(kl.lt.nlayer)
			if(jn.eq.lhet(kl))then		!���E��
				hetero = nint(rnd())	!0 or 1 random
				kl = kl + hetero
				exit
			elseif(jn.lt.lhet(kl))then
				exit
			endif
			kl=kl+1
		enddo
      end select
	ka = iarea(kl)
c
c-----���q�����G�l���M�[ 08/8/6 �|��-----
c	if((kl.eq.4).or.(kl.eq.6))then	!�`���l���w��In(0.53)Ga(0.47)As
	if(nchannel1.lt.kl.and.kl.le.nchannel2) then !channel 11/04/08�� !����!
		r1 = rnd()
		ntab = 0.0d0
		do i=0,i2max
			ntab = ntab + n2(i)
			if(r1.le.ntab)exit
		enddo
		if(i.gt.i2max)then
			i = i2max
		endif
		e = E2(i)
		akt = smh(kv,ka)*sqrt(e*(1+af(kv,ka)*e))	!120126homma
		akh = akt
c	elseif(kl.eq.5)then		!�`���l���w��InAs ����ύX�ɔ����폜
c		r1 = rnd()
c		ntab = 0.0d0
c		do i=0,i2max
c			ntab = ntab + n2(i)
c			if(r1.le.ntab)exit
c		enddo
c		if(i.gt.i2max)then
cc			i = i2max
c		endif
c		e = E2(i)
c		akt = smh(kv,ka)*sqrt(e)
c		akh = akt
	else
c	20030121�ύX v1.3.1(�Έ�)
c	�������������̎Z�o(���x) *����������͍l�������D
		akt = smh(kv,ka)*sqrt(-bktqn(ni)*alog(rnd()))
c	�������������̎Z�o(�G�l���M�[)
		akh = smh(kv,ka)*sqrt(-bktqn(ni)*alog(rnd()))
	endif
c----------------------------------------
	cf = 1.0-2*rnd()				!120126homma cos��
	sf = sqrt(1.0-cf*cf)			!120126homma sin��
	fi = 2*pi*rnd()					!��
c
c--	�������������̃p�����[�^����  ---
	n = jpnum + ni
	select case (intyp)
      case (1)	!�㑤�i�\�ʁj�d��
		p(3,n) = akt*cf				!120126homma
		p(6,n) = dz*rnd()*0.5
	case (2)	!�����d��
		p(1,n) = akt*sf*cos(fi)		!120126homma
		p(5,n) = dx*rnd()*0.5
      case (3)	!�E���d��
		p(1,n) = akt*sf*cos(fi)		!120126homma
		p(5,n) = dx*(float(nx)-rnd()*0.5)
      end select

c
c--	�������������̃p�����[�^����  ---
	select case (intyp)
      case (1)	!�㑤�i�\�ʁj�d��
		p(1,n) = akh*sf*cos(fi)		!120126homma
		if(jn.eq.0)then
			p(5,n) = dx*rnd()*0.5
		elseif(jn.eq.nx)then
			p(5,n) = dx*(float(nx)-rnd()*0.5)
		else
			p(5,n) = dx*(rnd()+float(jn)-0.5)
		endif
	case (2,3)	!���E���i���ʁj�d��
		p(3,n) = akh*cf				!120126homma
		if(jn.eq.0)then		!��[
			p(6,n) = dz*rnd()*0.5
		elseif(jn.eq.nz)then	!���[
			p(6,n) = dz*(float(nz)-rnd()*0.5)
		elseif(hetero.ne.-1)then	!�w�e���E��
			if(hetero.eq.0)p(6,n) = dz*(float(jn)-rnd()*0.5)
			if(hetero.eq.1)p(6,n) = dz*(float(jn)+rnd()*0.5)
		else						!����ȊO
			p(6,n) = dz*(rnd()+float(jn)-0.5)
		endif
	end select
	n = jpnum + ni
	p(2,n) = akh*sf*sin(fi)			!120126homma �x������
c	p(4,n) = -alog(rnd())*pgm(kv,ken,ka)	!���R�s������
	p(4,n) = -alog(rnd())*pgm(kv,ken,kpart)	!���R�s������	!07/8/4 �s�����U��
	kp(3,n) = kl

		x_start(n) = p(5,n)		!120817sato	
		z_start(n) = p(6,n)		!120817sato

c----(�m�F�p�F�����������q�����̏����G�l���M�[���z)---- 120126homma
		iz = max(min(ifix(p(6,n)/dz+0.5),nz),0)
		sk = p(1,n)*p(1,n)+p(2,n)*p(2,n)+p(3,n)*p(3,n)
		sq = sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*sk)
		e = (sq-1.0)/2/af(kv,ka)
		ie2=max(min((nint((e)/(de(2)*4))+1),(250)),1)
		if(intyp.eq.2)then
			if((iz.eq.22))	!�`���l���̈�����
     &   			ecr(1,ie2) = ecr(1,ie2) + 1
			if((iz.eq.23)) 
     &   			ecr(2,ie2) = ecr(2,ie2) + 1
			if((iz.eq.31)) 
     &   			ecr(3,ie2) = ecr(3,ie2) + 1
			if((iz.eq.34))
     &   			ecr(4,ie2) = ecr(4,ie2) + 1
			ecr(5,ie2) = ecr(5,ie2) + 1
		endif
c------------------------------------------------------
	enddo						!jcreat(���[�v�ϐ�ni)
c
c--	���̑������̃p�����[�^����  ---
	kp(1,jpnum:jpnum+jcreat) = kv
	kp(2,jpnum:jpnum+jcreat) = ken
c
	return
	end subroutine creat2
