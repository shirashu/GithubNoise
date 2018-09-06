c
c-----���q�����E�ʂŔ��˂܂��͐i������C�x���g-----
	subroutine border(
     &					hhm,af,af2,af4,eg,
     &					akx,aky,akz,kv,kl,kl2,ka,iarea,lhet,iz)	!2006/12/09 Hara
	implicit none
c
	include 'arraysize.fi'
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
c
c=================================================================================================
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