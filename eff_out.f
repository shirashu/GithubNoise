c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c Effective Potential �̏o�͂��s���T�u���[�`��
c
c �����|�e���V�����̏o��
c 2007/1/9, Hara, Ver.1.1
c
c	Input : ict, u_eff
c	Output: NOTHING

c	����F���ώ����|�e���V�������o��	
c	�o�͐�Feff_potential.txt,eff_potential2.txt�ɏo��
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	subroutine eff_out(ict, u_eff1,u_eff2,u_eff3,epA,epB,epC,
     &					epA2,epB2)	!120126homma
	implicit none
c
c===����===
c---�ϐ��z��p�p�����[�^---
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c--- �V�~�����[�V�������� ---
	common /outp/joutc,jouti,jouta,jsave,jheat,jcput,joutpi,sw
	integer	joutc,jouti,jouta,jsave,jheat,jcput,joutpi,sw(8)
c---��{�p�����[�^---
	integer	ict
c---�f�o�C�X�����--- 120126homma
	real,	dimension (0:nx,0:nz)	:: u_eff1,u_eff2,u_eff3
	real,	dimension (0:nx,0:nz)	:: epA,epB,epC,epA2,epB2
c----���[�J���ϐ�---- 120126homma
	integer	i,j
	integer,save	::	ncount
	real(8),save,allocatable 		:: avu_eff1(:,:)
	real(8),save,allocatable 		:: avu_eff2(:,:)
	real(8),save,allocatable		:: avu_eff3(:,:)
	integer ix,iz
c
c=========================================================================================================
c
c----(���A���^�C���lepA,epB,epC)--- !120126homma
	if (.not. allocated(avu_eff1)) then		!�ŏ��������s
		allocate(avu_eff1(0:nx,0:nz))
		allocate(avu_eff2(0:nx,0:nz))
		allocate(avu_eff3(0:nx,0:nz))
		avu_eff1 = 0.0
		avu_eff2 = 0.0
		avu_eff3 = 0.0
		ncount	= 0
	endif
	avu_eff1 = avu_eff1 + u_eff1
	avu_eff2 = avu_eff2 + u_eff2
	avu_eff3 = avu_eff3 + u_eff3
	ncount	= ncount + 1
	epA = u_eff1
	epB = u_eff2
	epC = u_eff3
c----(���ϒlepA2,epB2)------------- !120126homma
	if(modulo(ict,jouta).eq.(jouta-1))then
		avu_eff1 = avu_eff1/ncount
		avu_eff2 = avu_eff2/ncount
		avu_eff3 = avu_eff3/ncount
		rewind 72; rewind 76; rewind 77; rewind 78
		do j=0,nz
		do i=0,nx
			write(72,*) sngl(avu_eff1(i,j))	!eff_potential2.txt �o��
			write(76,*) sngl(avu_eff2(i,j))	!eff_epB2.txt
			write(77,*) sngl(epA(i,j))		!eff_epA.txt
			write(78,*) sngl(epB(i,j))		!eff_epB.txt
		enddo
		enddo
		epA2 = avu_eff1
		epB2 = avu_eff2
		avu_eff1 = 0.0
		avu_eff2 = 0.0
		avu_eff3 = 0.0
		ncount = 0
	endif
c----------------------------------
c	
	return
	end
