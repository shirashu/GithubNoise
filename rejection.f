	subroutine rejection(hhm,hm,af2,af4,ec,p,kp,iarea,bkq,
     &					 efermi,ix,iz,
     &					 akx,aky,akz,kv,kl,rflag,
c     &					 adkx,adky,adkz,iak)
     &					 adkx,adky,adkz,iak,avsumtel,epA,u,eg)		!08/8/6 �|��

	implicit none
c===����===
c---�ϐ��z��p�p�����[�^---
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c---��{�p�����[�^---
	real(8) bkq
	real rnd
c---�f�o�C�X�\��---
	integer(1)	iarea(nlayer)
      real	eg(nvalley,narea)
c---���q���---
	integer(1),dimension (3,npmax)	:: kp
	real	p(6,npmax)
c---�̈�ʃp�����[�^---
	real,	dimension (nvalley,narea) :: hhm,hm,af2,af4
	real,	dimension (nvalley,narea) :: ec
      real,	dimension (0:nx,0:nz)	:: epA !101220
	real,	dimension (0:nx,0:nz)	:: u					!classical�ȃ|�e���V����
c----(���b�V��&�J��fermi�G�l���M�[,�d�q���x)---
	real efermi(0:nx,0:nz,nvalley)
	real avsumtel(0:nx,0:nz,nvalley)
c----(dt�O�̃��b�V�����J�ʃh���t�g�g������)---	
	real adkz(0:nx,0:nz,0:nvalley)
	real adkx(0:nx,0:nz,0:nvalley)
	real adky(0:nx,0:nz,0:nvalley)
c----(���b�V��&�J�ʃt�F���~-�f�B���b�N���z�֐�)---
c	real fd(0:nx,0:nz,nvalley)
	real fd		!yama fd��3�����z�񂩂�ϐ��ɕύX 08/1/2
	real fbkt
	real akx,aky,akz,x,z
	real(8) sk,eir,sq
	real iak
	real fdexp
c===���[�J���ϐ�===
	integer n,rflag
	integer	ix,iz
	integer(1)	kv,ken,kl,ka

	sk=0.0
	sq=0.0
	eir=0.0
	rflag=0
	fbkt=0.0
	iak=0.0
	fd=0.0

	ka = iarea(kl)
	iak=akx+aky+akz
c	�h���t�g����drikx���������G�l���M�[�̎Z�o

	sk=(abs(akx)-abs(adkx(ix,iz,kv)))*(abs(akx)-abs(adkx(ix,iz,kv)))+
     &   (abs(aky)-abs(adky(ix,iz,kv)))*(abs(aky)-abs(adky(ix,iz,kv)))+
     &   (abs(akz)-abs(adkz(ix,iz,kv)))*(abs(akz)-abs(adkz(ix,iz,kv)))

		if(af4(kv,ka).ne.0.0)then
			sq = sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*sk)
c			eir=(sq-1.0)/af2(kv,ka)+ec(kv,ka)
			eir=(sq-1.0)/af2(kv,ka)
c     &			        -epA(ix,iz)+u(ix,iz)-0.7*(eg(1,ka)-eg(1,2))			!101220
		else
c			eir=hhm(kv,ka)*sk+ec(kv,ka)
			eir=hhm(kv,ka)*sk
c     &                   -epA(ix,iz)+u(ix,iz)-0.7*(eg(1,ka)-eg(1,2))			!101220
		endif

		fbkt=bkq*avsumtel(ix,iz,kv)			!08/5/15 �|��
		fbkt=1.0d0/fbkt

		fdexp = (eir-efermi(ix,iz,kv))*fbkt

		if(fdexp.gt.10.0) then
			rflag = 1
			goto 100
		endif

c		fd(ix,iz,kv)=exp((eir-efermi(ix,iz,kv))*fbkt)+1.0	!Fermi_Dirac Distribution
		fd = exp((eir-efermi(ix,iz,kv))*fbkt)+1.0d0			!Fermi_Dirac Distribution

c		fd(ix,iz,kv)=1.0d0/fd(ix,iz,kv)
		fd = 1.0d0/fd
c
c		if ((1.0d0-fd(ix,iz,kv)).lt.rnd())then
		if ((1.0d0-fd).lt.rnd())then
			rflag = 1
		endif

  100 continue
c
	return
	end