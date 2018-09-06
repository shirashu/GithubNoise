	subroutine drift(am_aniso,aff_aniso,fx,fz,dz,dhet,dt,de,pgm,af2,
     &						af4,hm,hhm,
     &				        akx,aky,akz,ts,x,z,t1,tau,
     &						kv,ken,kl,kl2,iflag,ie,ei,sk,
     &						kpart,ka,n)									!07/8/4 �s�����U��


	implicit none
c	
c	�I�����̏�Ԃ�\���t���O�Ƃ��āAiflag��ݒ肷��B
c	iflag=0 : �h���t�g��A�U���C�x���g����
c	iflag=1 : �h���t�g��A���E�ʂɎ������̂Ői������(subroutine border)
c	iflag=2	: �h���t�g��Adt�Ɏ������̂Ŏ��̃��[�v�i���q�j�ɓ���
c
c	��{�I��###�ň͂񂾕����̓G���[�����ł���A���ۂɂ͋N���Ȃ�������
c	���������ꍇ�ɂ̂ݓ��삷��̂ŁA���i�͖������Ă��܂�Ȃ�
c
	include 'arraysize.fi'
	real fx,fz
	real dhet(0:nlayer),dt
	real de(nenergy)
	real pgm(nvalley,nenergy)
	real af2,af4,hm,hhm
c	real,	dimension (nvalley,narea)	:: smh,hhm,hm,af,af2,af4
	real akx,aky,akz,x,z,t1,tau
	real(8)	ts
	integer(1) kv, ken, kl, kl2,ka
	integer iflag
	integer kpart  !07/8/4 �s�����U��
	real, dimension	(nvalley,narea,nvalley)::am_aniso							 !20100624
	real, dimension	(nvalley,narea,nvalley)::aff_aniso							!20100624
c
	real(8)	dkx,dkz
	real	xx,zz,err,dh
c	real(8)	skx,sky
	real(8)	sk,ei
	real	ddt
	real	tauz
	integer	ie
c	real	den
	integer i !,il
	integer n
	real	dz
c
	real(8)	qht,sq,cp,sq21,sq22
	real qh !qh = q/h
	parameter(qh = 1.51925E+15)
	integer klc

c-----------------------------------------------------------------------------------------

c
c



c	--- ���炩�̌�����kv(�JNo.)=0�ł���ꍇ
c		�h���t�g�v�Z�͍s��Ȃ� ---
	if(kv.eq.0) return
c
c#####�т̔���##################################################################
c---	!���q���E�ʏ�ɂ��� �i���C�x���g���������܋N����ꍇ---
c	!iflag=1 �Ƃ��� �h���t�g�͍s��Ȃ�
	if((dhet(kl-1).lt.z).and.(z.lt.dhet(kl)))then
		continue	!�ʏ퓮��
	elseif(z.eq.dhet(kl-1))then
		if(akz.lt.0)then
			kl2 = kl-1
			iflag = 1
			return
		elseif(akz.eq.0)then
			fz = 0.0	!����P�[�X�i��O�����j
			continue
		endif
	elseif(z.eq.dhet(kl))then
		if(akz.gt.0)then
			kl2 = kl+1
			iflag = 1
			return
		elseif(akz.eq.0)then
			fz = 0.0	!����P�[�X (��O����)
			continue
		endif
c	---	���q�p�����[�^�Ɨ��q�ʒu����v���Ȃ��ꍇ(�G���[���)
c	���Z�X�̈�ɂđ��������������G���[���ȉ��ŉ�� 11/05/26 ��
	elseif((kl.eq.1).and.(kpart.eq.2))then
	klc = 0
	do klc=1,3
		if((dhet(klc-1).lt.z).and.(z.lt.dhet(klc))) then
			kl = klc
			exit
		 endif
	enddo
c
c	---	���q�p�����[�^�Ɨ��q�ʒu����v���Ȃ��ꍇ(�G���[����) 
	else
		write(* ,*)'kl�s��(drift)'
		write(99,*)'kl�s��(drift)'
	write(*,*) n,x,z,kl,kpart
	write(*,*) n,dhet(kl-1), dhet(kl)
!SMP$	ASSERT(ITERCNT(1))
c		do il = 1, nlayer
c			if((dhet(il-1).lt.z).and.(z.lt.dhet(il)))then
c				kl = il
c			elseif(z.eq.dhet(il))then
c				kl = il
c				fz = 0.0	!����P�[�X�i��O�����j
c			endif
c		enddo
		kv=0
		return
	endif
c--------------------------------------------------------------
c
c	--- ���炩�̌�����t1>dt�ł���ꍇ
c		�h���t�g�v�Z�͍s��Ȃ��Ŏ��̃��[�v(iflag=2) ---
	ddt = dt*(1.0+epsilon(dt))	!�덷�͈�
	if(t1.lt.dt)then
		continue	!���큨������������
	elseif((t1.ge.dt).and.(t1.lt.ddt))then
		iflag = 2
		return
	elseif(t1.ge.ddt)then
		iflag = 2
		!write( *,*)'t1>dt',t1,dt,ddt
		!write(99,*)'t1>dt',t1,dt,ddt
		return
	endif
c
1000	continue
c
c	--- �т̌v�Z ---
	if(sngl(ts).gt.dt)then
c	--- ���̒P�ʎ��ԓ��ŎU���C�x���g���N���Ȃ��ꍇ
c		��=dt-t1�����h���t�g���Aiflag=2 ---
		tau = dt-t1
		iflag = 2
	else
		tau= sngl(ts-dble(t1))
c	--- ���炩�̌����Ń�=0�ł���ꍇ
c		�h���t�g�v�Z�͍s��Ȃ��ŎU��(iflag=0) ---
		if(tau.gt.0.0)then
			continue	!���큨������������
		elseif(tau.eq.0.0)then
			iflag = 0
			write( *,*)'��=0'
			write(99,*)'��=0'
			return
c	---	tau��0��菬�����ꍇ(�G���[)   ---
		elseif(tau.le.0.0)then
			write(* ,*)'�уG���[(drift)'
			write(99,*)'�уG���[(drift)'
			write(99,*) tau,t1,dt,sngl(ts)
			tau = 0.0
			iflag = 0
			return
		endif
	endif
c###############################################################################
c
c
c#####�h���t�g�^��##############################################################
c	--- ���R�s�����ԃт����h���t�g�v�Z ---
	qht=qh*tau
c
	dkx=-qht*fx		!dkx=��akx: akx�́A���ԃтł̑�����
	dkz=-qht*fz		!dkz=��akz: akz�́A���ԃтł̑�����
c
	if(af4.eq.0.0)then
		sq=1.0
	else
c		sq=1.+af4*hhm*(akx*akx+aky*aky+akz*akz)
c	sq21 = sq
		sq=1.+af4*hhm*
     &	((akx*akx)/am_aniso(kv,ka,1)+(aky*aky)/
     &	 am_aniso(kv,ka,2)+(akz*akz)/am_aniso(kv,ka,3))
     &	 *(am_aniso(kv,ka,1)+am_aniso(kv,ka,2)+
     &	 am_aniso(kv,ka,3))/3.0d0
c	akx
c	sq22 = sq
c	write( *,*) sq21,sq
c		sq=1.+af4*hhm*
c     &	((akx**2.0d0)/am_aniso(kv,ka,1)+(aky**2.0d0)/
c     &	 am_aniso(kv,ka,2)+(akz**2.0d0)/am_aniso(kv,ka,3))
c     &	 *(am_aniso(kv,ka,1)+am_aniso(kv,ka,2)+
c     &	 am_aniso(kv,ka,3))/3.0d0
		if(sq.gt.0.0)then
			sq=1.0/sqrt(sq)
		else
			write(99,*)'sq�̒l���s���ł�(drifttau)'
			stop
		endif
	endif
	cp=hm*tau
	xx = x + cp*(akx+0.5*dkx)*sq
	zz = z + cp*(akz+0.5*dkz)*sq
c###############################################################################
c	check_point
c################�Q�ƃG�l���M�[�e�[�u���͂ݏo���`�F�b�N##########################
	sk	= (akx+dkx)*(akx+dkx)+aky*aky+(akz+dkz)*(akz+dkz)
	if((sk.gt.0.0).and.(sk.le.huge(sk)))then
		if(af2.ne.0.0)then
			ei=(sqrt(1.0+af4*hhm*sk)-1.0)/af2		!eV
		else
			ei=hhm*sk
		endif
	else
		write(99,*)'sk�̒l���s���ł�(drift)',sk,akx+dkx,aky,akz+dkz
		write( *,*)'sk�̒l���s���ł�(drift)',sk,akx+dkx,aky,akz+dkz
	endif
	ie=int(ei/de(ken))+1
c
c=============�G���[����========================================
	if((ei.lt.0.0).and.(ei.ge.huge(ei)))then
		write( *,*)'energy�v�Z�G���[(emcd2)',ei
		write(99,*)'energy�v�Z�G���[(emcd2)',ei
		kv=0
		return
	endif
c===============================================================
c
	if(ie.ge.nemax)then
		if(ken.lt.nenergy)then
c	----	�ċA���� ->�h���t�g�����蒼��---------------
			ts=(ts-dble(t1))/dble(pgm(kv,ken))
     &				*dble(pgm(kv,ken+1))+dble(t1)
			ken = ken+1
			iflag = 0
			goto 1000
c	--------------------------------------------------
		else
			write( *,*)'energy���傫�����܂�',ie  !,n
			write(99,*)'energy���傫�����܂�',ie  !,n
			ie = nemax
		endif
	endif


c#########################################################################
c	�h���t�g�^���O��̈ʒux,xx����̈�AB��ʉ߂��闱�q���J�E���g  120817sato


c#########################################################################
c	--- �h���t�g��ʒu�����E�ʂ��܂����ł��Ȃ���Ίm��E�I�� ---
	if((dhet(kl-1).lt.zz).and.(zz.lt.dhet(kl)))then
		kl2 = kl
		akx = akx+dkx
		akz = akz+dkz
		x = xx
		z = zz
		t1 = t1 + tau
		return
c	--- �h���t�g�オ���E�ʂ��܂����ł����炻�̋��E�ʂ܂Ńh���t�g ---
	elseif(zz.lt.dhet(kl-1))then
		kl2 = kl-1
		dh=dhet(kl-1)		!��Ɋђ�
		continue
	elseif(zz.gt.dhet(kl))then
		kl2 = kl+1
		dh=dhet(kl)		!���Ɋђ�
		continue
c	--- �h���t�g��ɂ��傤�ǊE�ʂɐڐG
c		�h���t�g�Čv�Z�͍s��Ȃ��ŊE�ʃC�x���g(iflag=1) ---
	elseif(zz.eq.dhet(kl-1))then
		kl2 = kl-1
		iflag = 1
		akx = akx+dkx
		akz = akz+dkz
		x = xx
		z = zz
		t1 = t1 + tau
		return
	elseif(zz.eq.dhet(kl))then
		kl2 = kl+1
		iflag = 1
		akx = akx+dkx
		akz = akz+dkz
		x = xx
		z = zz
		t1 = t1 + tau
		return
	endif
c
c	--------------------------------------
	err = dz*0.000005	!1.6e-14
c	call calctau(akx,aky,akz,z,af4,hm,hhm,fz,dh,tauz)
	call calctau(am_aniso,kv,ka,akx,aky,akz,z,af4,hm,
     &				hhm,fz,dh,tauz)

	do i=1,4
		if((tauz.gt.0.0).and.(tauz.le.tau))then
c		---	tauz���펞 ---
			tau = tauz
			ts = dble(t1) + dble(tauz)
c#####�h���t�g�^��##############################################################
c	--- ���R�s�����ԃт����h���t�g�v�Z ---
			qht=qh*tau
c
			dkx=-qht*fx		!dkx=��akx: akx�́A���ԃтł̑�����
			dkz=-qht*fz		!dkz=��akz: akz�́A���ԃтł̑�����
c
			if(af4.eq.0.0)then
				sq=1.0
			else
				sq=1.+af4*hhm*(akx*akx+aky*aky+akz*akz)
	
				sq=1.+af4*hhm*
     &	((akx**2.0d0)/am_aniso(kv,ka,1)+(aky**2.0d0)/
     &	 am_aniso(kv,ka,2)+(akz**2.0d0)/am_aniso(kv,ka,3))
     &	 *(am_aniso(kv,ka,1)+am_aniso(kv,ka,2)+
     &	 am_aniso(kv,ka,3))/3.0d0

				if(sq.gt.0.0)then
					sq=1.0/sqrt(sq)
				else
					write(99,*)'sq�̒l���s���ł�(drifttau)'
					stop
				endif
			endif
			cp=hm*tau
			xx = x + cp*(akx+0.5*dkx)*sq
			zz = z + cp*(akz+0.5*dkz)*sq
c###############################################################################
c
			if((dh-err.le.zz).and.(zz.le.dh+err))then
c				---	zz�̈ʒu���덷�͈͓�(zz�}(0..err) = dh) ---
				akx = akx+dkx
				akz = akz+dkz
				x = xx
				z = dh
				t1 = t1 + tau
				iflag = 1
	return
			elseif((dhet(kl-1).lt.zz).and.(zz.lt.dhet(kl)))then
c			elseif(kl.eq.kl2)then
c				---	zz�̈ʒu���i���O ---
				akx=akx+dkx
				akz=akz+dkz
				x = xx
				z = zz
				t1 = t1 + tau
				call calctau(am_aniso,kv,ka,akx,aky,akz,z,af4,hm,
     &				hhm,fz,dh,tauz)
		cycle
			else
c				---	zz�̈ʒu���i�����Ă��� ---
				tauz = tauz*0.9
		cycle
			endif
		elseif((tauz.eq.0.0).and.(z.eq.dh))then
c		---	tauz�ُ펞 ---
			write(*,*)'tauz=0,z=dh'
	return
		else
c******************** �G���[���� **********************************
			write( *,*)'tauz�G���[(drift)' !,sq21,sq22,tauz,tau
			write(99,*)'tauz�G���[(drift)'
			write(99,*) tauz,tau,t1,fz
			write(99,'(4(I3))') iflag,kl,kv,i
			write(99,*) akx,aky,akz,sngl(ts),x,z,xx,zz
			if(abs(z-dh).lt.abs(zz-dh))then
				tau = 0.0
				z = dh
				iflag = 1
	return
			else
				akx=akx+dkx
				akz=akz+dkz
				x = xx
				z = dh
				t1 = t1 + tau
				iflag = 1
	return
			endif
		endif
	enddo
c		---	���[�v���Ɏ���������Ȃ������ꍇ ---
	if(zz.ne.dh)then
		write(* ,*)'�������܂���(drift)'
		write(99,*)'�������܂���(drift)'
		write(99,*) tau,t1,fz
		write(99,'(4(I3))') iflag,kl,kv,i
		write(99,*) akx,aky,akz,sngl(ts),x,z,xx,zz
		z = dh
		iflag = 1
	return
	endif
c
	write( *,*)'�_���G���[(drift)'
	write(99,*)'�_���G���[(drift)'
c	**********************************************************
	stop
	end subroutine drift
c
c
c
c
      subroutine calctau(am_aniso,kv,ka,akx,aky,akz,z,af4,
     &	hm,hhm,fz,dh,tauz)
c	���q�����E��dh�ɐڐG���鎞��tauz���v�Z����B
c	hm*akz*t - Fz*(q/2m*)*(t**2) -(dh-z)sq = 0 
c	�̕������̉� tauz1,tauz2(tauz1<tauz2)�𓱂��A
c	tauz1>0�̂Ƃ�tauz=tauz1�Atauz1<0�̂Ƃ�tauz=tauz1�Ƃ���B
c	�d���ł�tauz1=tauz2=tauz�A
c	���Ȃ��̏ꍇtauz=���ƂȂ�B
c
	implicit none
	

	real, dimension	(3,5,3)::am_aniso							 !20100624
	integer(1) kv,ka															!20100624

c	real	z
	real	akx,aky,akz,z
c	integer	iflag
	real	af4,hm,hhm
	real	fz
	real	dh
	real	tauz
c
	real	tauz1,tauz2
	real(8)	akz2
	real	akk
	real(8)	sq,sq2
	real(8)	hqfz
	real(8)	b24ac
	real(8)	ac4
c==========
	real	hq
	parameter(hq = 6.58218E-16)
c
	akz2=akz*akz
	akk =akx*akx+aky*aky+akz2
	akk = ((akx**2.0d0)/am_aniso(kv,ka,1)+(aky**2.0d0)/
     & am_aniso(kv,ka,2)+(akz**2.0d0)/am_aniso(kv,ka,3))
     & *(am_aniso(kv,ka,1)+am_aniso(kv,ka,2)+
     & am_aniso(kv,ka,3))/3.0d0

	if(af4.eq.0.0)then
		sq=1.0
	else
		sq	=1.+af4*hhm*akk
		if(sq.ge.0.0)then
			sq=sqrt(sq)
		else
			write(99,*)'sq�̒l���s���ł�(calctau)'
			stop
		endif
	endif
c
	if(fz.eq.0.0)then
c		�ѓ��̍���������̂ŁA���`�Ȏ��ŉ���
		tauz = (dh-z)/akz/hm*sq
	else
		hqfz = hq/fz
		ac4  = dble(dh-z)*fz/hhm*sq
		if(dh.eq.z)then
			tauz1 = min(0.0,sngl(hqfz*(2.0*akz)))
			tauz2 = max(0.0,sngl(hqfz*(2.0*akz)))
			if(tauz1.gt.0.0) tauz = tauz1
			if(tauz1.le.0.0) tauz = tauz2
		elseif(abs(ac4*8e5).lt.abs(akz2))then
c			fz���������ƌ덷���傫���Ȃ�̂Ő��`�ߎ��ŉ���
			tauz = (dh-z)/akz/hm*sq
		else
c			������ Fz*(q/2m*)*(t**2)-hm*akz*t+(dh-z)sq = 0 ������
c			t = qFz/h*(akz�}sqrt(akz**2-(dh-z)*Fz/hhm*sq)
			b24ac = (akz2-ac4)
			if(b24ac.gt.0.0)then
c				�������̉� tauz1,tauz2
				sq2   = sqrt(b24ac)
				tauz1 = min(hqfz*(akz-sq2),hqfz*(akz+sq2))
				tauz2 = max(hqfz*(akz-sq2),hqfz*(akz+sq2))
				if(tauz1.gt.0.0) tauz = tauz1
				if(tauz1.le.0.0) tauz = tauz2
			elseif(b24ac.eq.0.0)then
c				�������̉� �� �d��
				tauz = hqfz*(akz)
			else
c				�������̉��Ȃ� tauz����
				tauz = huge(tauz)
			endif
		endif
	endif
c
	return
	end
