c
c-----���q�����E�ʂŔ��˂܂��͐i������C�x���g-----
	subroutine border_roughness1(
     &					hhm,af,af2,af4,eg,
     &					akx,aky,akz,kv,kl,kl2,ka,iarea,lhet,iz,	!2006/12/09 Hara
     &					swk_rou,ken,ix,de,					!121029sato
     &					roughness1_countx,roughness1_counte)
	implicit none
c
	include 'arraysize.fi'
	common	/arraydata/nx,nz					!121029sato
	integer	nx,nz
	real,	dimension (nvalley,narea)	:: hhm,af,af2,af4,eg
	real	akx,aky,akz
	integer(1)	kv,kl,kl2,ka,iz	!2006/12/09 Hara
	integer(1)	iarea(nlayer)
	integer(2)	lhet(nlayer)	!2006/12/09 Hara
c
c	---	���[�J���ϐ�	---
c	skx,sky,skz:�e�����g���̓����i�[����ϐ�
c	vb:��Ǎ������i�[����ϐ�
c	ex,ey,ez:�e�����^���G�l���M�[���i�[����ϐ�
	real(8)	skx,sky,skz
	real	vb
	real(8)	ex,ey,ez
	integer(1)	ka2

c------���t�l�X�U��--R.P. Joshi���f��(2012�N�H)-------------------------!121029sato
	real	swk_rou(nemax,nvalley,nenergy,narea)		!���t�l�X�U�����[�g
	real	swkd_rou				!					!���t�l�X�U�����[�g���`���
	integer,dimension (0:nx,narea,2) :: roughness1_countx
	integer,dimension (nemax,narea,2) :: roughness1_counte
	real	de(nenergy)
	real	den_2deg
	real(8)	ei_2deg
	real(8)	sk_2deg
	integer	ie_2deg
	integer(1)	ken
	integer		ix
	real	cs,sn
	real	rnd		!����
c

c=======================================================================================!121029sato
	if((iz.eq.lhet(nchannel1)).or.(iz.eq.lhet(nchannel2))) then		!�`���l���̃w�e���E�ʂɂ���Ƃ�
c	open(unit=1299,file='zztriela.txt')	!�m�F�p!

c------���t�l�X�U��--2DEG�G�l���M�[-------------------------

	sk_2deg	= akx*akx+aky*aky			!2DEG�̉^���ʂ���G�l���M�[�����߂�
	if((sk_2deg.ge.0.0).and.(sk_2deg.le.huge(sk_2deg)))then
		if(af2(kv,ka).ne.0.0)then
			ei_2deg=(sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*sk_2deg)-1.0)/af2(kv,ka)		!eV
		else
			ei_2deg=hhm(kv,ka)*sk_2deg
		endif
	else
		write(99,*)'sk_2deg�̒l���s���ł�',sk_2deg,akx,aky
		write( *,*)'sk_2deg�̒l���s���ł�',sk_2deg,akx,aky
	endif
									!ei_2deg	!2deg�G�l���M�[���� 
	ie_2deg = int(ei_2deg/de(ken))+1			!2deg�G�l���M�[������(���U�l) 
	den_2deg=ei_2deg/de(ken)-float(ie_2deg-1)	!de_2deg:ie_2deg��ei_2deg�̍��i�Y���j

c------���t�l�X�U��--�U�����[�g�̐��`���-------------------------
						
	if(ie_2deg.eq.1) then						
		swkd_rou=swk_rou(ie_2deg,kv,ken,ka)	
											
	else									
		swkd_rou=den_2deg*swk_rou(ie_2deg,kv,ken,ka)	+
     &					(1-den_2deg)*swk_rou(ie_2deg-1,kv,ken,ka)	
													
		if(swkd_rou.gt.1.0)then		!�G���[����	
			write( *,*)'swkd_rou��1.0�𒴂��܂���',swkd_rou
			write(99,*)'swkd_rou��1.0�𒴂��܂���',swkd_rou
			swkd_rou = swkd_rou/swkd_rou	!�ċK�i�� swkd_rou=1					
		endif													
	endif
c------���t�l�X�U��--�U�����[�g�Q��-------------------------

	if(rnd().lt.swkd_rou)then		!���t�l�X�U�����[�g�������Ƃ�
		cs  = 1.0-2.0*rnd()
		if(abs(cs).ge.1.0)then
		  sn = 0.0
		else
		  sn  = sqrt(1.0-cs*cs)
		endif

		akx = sqrt(sk_2deg)*cs		!xy���ʓ�����
		aky = sqrt(sk_2deg)*sn		!xy���ʓ�����		
		
c------���t�l�X�U��--�J�E���^-------------------------
		
		if(iz.eq.lhet(nchannel1)) then		!�`���l���@��w�e���E��
			roughness1_countx(ix,ka,1) =roughness1_countx(ix,ka,1) +1
			roughness1_counte(ie_2deg,ka,1) 
     &							  =roughness1_counte(ie_2deg,ka,1) +1
		endif	
		if(iz.eq.lhet(nchannel2)) then		!�`���l���@���w�e���E��
			roughness1_countx(ix,ka,2) =roughness1_countx(ix,ka,2) +1
			roughness1_counte(ie_2deg,ka,2) 
     &							  =roughness1_counte(ie_2deg,ka,2) +1
		endif


	endif		!���t�l�X�U�����[�g�������Ƃ�

	endif		!�`���l���̃w�e���E�ʂɂ���Ƃ�

c=========================================================================================!121029sato
c
c	--- z�����̉^���G�l���M�[�v�Z ---
	if(af4(kv,ka).ne.0.0)then
		skz = dble(akz)*dble(akz)
		ez=(sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*skz)-1.0)/af2(kv,ka)		!eV
	else
		ez=hhm(kv,ka)*akz*akz
	endif
c	---------------------------------
c
	ka2 = iarea(kl2)
c	vb =0.7*(eg(kv,ka2)-eg(kv,ka))		!Vb:��Ǎ���[eV]
c	if(iz.lt.(nlayer-2)-5) then	!�h���t�g�O�ʒu��classical�̈�̂Ƃ�(2006/08/08����)!2006/12/09 Hara
c	if(iz.lt.(nlayer-4)-4) then	!�h���t�g�O�ʒu��classical�̈�̂Ƃ�(2006/08/08����)!2006/12/22 Hara
c		vb =0.7*(eg(kv,ka2)-eg(kv,ka))		!Vb:��Ǎ���[eV]
c	else						!�h���t�g�O�ʒu��EP�̈�̂Ƃ�
		vb = 0.0			!�Ȃ� 0 eV ?? 2011/03/23��
c	endif
c
	if(ez.gt.vb)then					!���Vb��藱�q�̃G�l���M�[�� ������
c	===	��ǂ����z����ꍇ�i�i���j	===
		ez = ez-vb			!z�����̃G�l���M�[�����ǃG�l���M�[���Ђ�
		skx = dble(akx)*dble(akx)
		sky = dble(aky)*dble(aky)
		if(af4(kv,ka).ne.0.0)then
			ex=(sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*skx)-1.0)/af2(kv,ka)		!eV
			ey=(sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*sky)-1.0)/af2(kv,ka)		!eV
		else
			ex=hhm(kv,ka)*skx		!�����O�̎��͒P���v�Z
			ey=hhm(kv,ka)*sky
		endif
c
		kl = kl2
		ka = ka2
c		---	���߂��G�l���M�[����g������i�����O�ł��v�Z�\�j
		akx=sign(sngl(sqrt(ex*(1.0+af(kv,ka2)*ex)/hhm(kv,ka2))),akx)
		aky=sign(sngl(sqrt(ey*(1.0+af(kv,ka2)*ey)/hhm(kv,ka2))),aky)
		akz=sign(sngl(sqrt(ez*(1.0+af(kv,ka2)*ez)/hhm(kv,ka2))),akz)
c		---	sign(a,b) = a�̐�Βl��b�̕������|�����l
c
	else
c	===	��ǂ����z���Ȃ��ꍇ�i���ˁj	===
		akz=-akz
	endif
c
      return
	end