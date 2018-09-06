c-----����Ef�����肷��T�u���[�`��(���J�̂�)-----
	subroutine fermi_calc(efermi,ka3,
     &					  tel1,tel2,tel3,efermi1,efermi2,efermi3,
     &					  avsumconc,avsumtel,avsumtei1,ntab1,etab1,
     &                    avsumteiA,avsumconcA,lhet,am,aff)		!�|�ݕύX

	implicit none
c	===����===
c---�ϐ��z��p�p�����[�^---
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c---��{�p�����[�^---
	integer	ix,iz,iv,ia,i1,ind,flag3 !���������
	integer(1)	ka,ka3(0:nx,0:nz,nvalley)
	real, dimension (nvalley,narea)::am	
	real(8) aff(nvalley,narea)					!120201
c---�Z�x�EEf�ETel---
	real efermi(0:nx,0:nz,nvalley)
	real,dimension(:),save,allocatable	:: ntabb
c
	real,	dimension(0:nx,0:nz)	:: tel1,efermi1
	real,	dimension(0:nx,0:nz)	:: tel2,efermi2
	real,	dimension(0:nx,0:nz)	:: tel3,efermi3
c
c-----08/8/6 �|�ݒǉ�-----
	real avsumconc(0:nx,0:nz,nvalley)
	real avsumconcA(0:nx)          !101221�d�q�Z�x�`���l������2
	real avsumtel(0:nx,0:nz,nvalley)
	real(8) avsumtei1(0:nx,0:nz,nvalley)
	real(8) avsumteiA(0:nx)                 !101221�d�q�G�l���M�[�`���l������2
	real ini_fermi,ini_tel		!Fermi_Level,Tel�̏����l
c
	real ntab1(0:300000,5)		!�Z�x(�ԍ��A�ޗ�)
	real etab1(0:300000,5)		!�G�l���M�[(�ԍ��A�ޗ�)
	integer(2)	lhet(nlayer)
c-------------------------
c
	if(.not. allocated(ntabb))then
		allocate(ntabb(0:5))
		efermi  = -1.0595		!���q�̂��Ȃ��Ƃ����ka=4��X�JEg/2
		efermi1 = -1.0595		!���q�̂��Ȃ��Ƃ����ka=4��X�JEg/2
		efermi2 = -1.0595		!���q�̂��Ȃ��Ƃ����ka=4��X�JEg/2
		efermi3 = -1.0595		!���q�̂��Ȃ��Ƃ����ka=4��X�JEg/2
	endif
c
c-----���J�̂ݍl�����Ă��邽��Ec=0�ƍl���A��Ԗ��x�֐���(E-Ec)**(1/2)��Ec=0
c	�S�Ă̒J���l������ꍇ�AEc=0�ȊO���l����K�v����
	iv=1				!���J�̂ݍl��

c-----------101221�d�q�G�l���M�[�C�Z�x�`���l������------------
c-----�v���O�����Ɋ܂߂�F�V�[�g���ρ@�v���O��������O���F���b�V������ !120126homma
c	do iz = lhet(nchannel1)+1,lhet(nchannel2)-1
c	     do ix=0,nx
c	          avsumconc(ix,iz,1) =avsumconcA(ix)
c	          avsumtei1(ix,iz,1) =avsumteiA(ix)
c		 enddo
c      enddo
c-------------------------------------------------------------

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	do iz=38,46			!Channel�w��Ef��Tel�Z�o InAs
	do iz = lhet(nchannel1)+1,lhet(nchannel2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		do ix=1,nx-1
c
	ka = ka3(ix,iz,iv)		!�f��No.��1�`4
c
c-----���q�̂��Ȃ��Ƃ���̃t�F���~���x���Ɠd�q���x-----
      if(ka.eq.0)then
	goto 1201
      endif
c
	if(avsumconc(ix,iz,iv).eq.0.0)then
		efermi(ix,iz,iv) = -1.0595  !ka=4��X�J��Eg��1/2�̒l(=minEg)
		avsumtel(ix,iz,iv)=0.0		!���q�̂��Ȃ����b�V���͓d�q���x0[K]
		goto 1201
	endif
c------------------------------------------------------
c

c
c	write(*,*) 'ka',ka
	if(ka.eq.1)then			!InAs	11/04/21��
		goto 12
	elseif(ka.eq.2)then		!In(0.52)Al(0.48)As	 11/04/21��
		goto 13
	elseif(ka.eq.3)then		!!InAs	�@11/04/21��
		goto 14
	elseif(ka.eq.4)then		!!InAs
		goto 15
	elseif(ka.eq.5)then     !!InAs
		goto 16
	else
	continue
		write(99,*)'�K����ޗ����݂���܂���'
		stop
	endif
ccccccccccccccccccccccccccccccccccccccccccccccccc
c	�ȍ~�ݒ���@��������Ȃ�!!2011/03/22�@��
cccccccccccccccccccccccccccccccccccccccccccccccc
c
c-----InAs (old_In(0.52)Al(0.48)As)-----11/04/21��
c   12 continue
c		write(99,*)'InAlAs�ɓ���܂���'
c		stop
   12 continue    
	if((avsumconc(ix,iz,iv).gt.0.0).and.
     &(avsumconc(ix,iz,iv).le.1e18))then
		avsumconc(ix,iz,iv) = 1e18
	elseif(avsumconc(ix,iz,iv).ge.1e25)then
		avsumconc(ix,iz,iv) = 1e25
	endif
c
	do i1=0,100000,10
		if(avsumconc(ix,iz,iv).le.ntab1(i1,1))then
			avsumconc(ix,iz,iv) = ntab1(i1,1)
			exit
		endif
	enddo
c
	if(avsumtei1(ix,iz,iv).le.etab1(i1,1))then
		avsumtei1(ix,iz,iv) = etab1(i1,1)
	endif
c
c-----Ef�̏����l�̐ݒ�-----
	ini_fermi = -0.30
	ini_tel   = 300
	if((avsumconc(ix,iz,iv).ge.5e24)
     &.and.(avsumconc(ix,iz,iv).le.1e25))then
		ini_fermi = 0.20
	elseif((avsumconc(ix,iz,iv).ge.1e24)
     &.and.(avsumconc(ix,iz,iv).lt.5e24))then
		ini_fermi = 0.10
	elseif((avsumconc(ix,iz,iv).ge.5e23)
     &.and.(avsumconc(ix,iz,iv).lt.1e24))then
		ini_fermi = 0.03
	endif
c
	goto 1230
c
c-----Al(0.52)In(0.48)As (old_In(0.53)Ga(0.47)As)----- 11/04/21��
   13 continue
	if((avsumconc(ix,iz,iv).gt.0.0).and.
     &(avsumconc(ix,iz,iv).le.1e18))then
		avsumconc(ix,iz,iv) = 1e18
	elseif(avsumconc(ix,iz,iv).ge.1e25)then
		avsumconc(ix,iz,iv) = 1e25
	endif
c
	do i1=0,100000,10
		if(avsumconc(ix,iz,iv).le.ntab1(i1,2))then
			avsumconc(ix,iz,iv) = ntab1(i1,2)
			exit
		endif
	enddo
c
	if(avsumtei1(ix,iz,iv).le.etab1(i1,2))then
		avsumtei1(ix,iz,iv) = etab1(i1,2)
	endif
c
c-----Ef�̏����l�̐ݒ�-----
	ini_fermi = -0.30  ! -0.05
	ini_tel   = 300
	if((avsumconc(ix,iz,iv).ge.5e24)
     &.and.(avsumconc(ix,iz,iv).le.1e25))then
		ini_fermi = 0.20
	elseif((avsumconc(ix,iz,iv).ge.1e24)
     &.and.(avsumconc(ix,iz,iv).lt.5e24))then
		ini_fermi = 0.10
	elseif((avsumconc(ix,iz,iv).ge.5e23)
     &.and.(avsumconc(ix,iz,iv).lt.1e24))then
		ini_fermi = 0.03
	endif
c
	goto 1230
c
c-----InAs (old_InAs�@)-----
   14 continue
	if((avsumconc(ix,iz,iv).gt.0.0).and.
     &(avsumconc(ix,iz,iv).le.1e18))then
		avsumconc(ix,iz,iv) = 1e18
	elseif(avsumconc(ix,iz,iv).ge.1e25)then
		avsumconc(ix,iz,iv) = 1e25
	endif
c
c-----���ϔZ�x�ƕ��σG�l���M�[�̕␳-----
	do i1=0,100000,10
		if(avsumconc(ix,iz,iv).le.ntab1(i1,3))then
			avsumconc(ix,iz,iv) = ntab1(i1,3)
			exit
		endif
	enddo
c
	if(avsumtei1(ix,iz,iv).le.etab1(i1,3))then
		avsumtei1(ix,iz,iv) = etab1(i1,3)
	endif
c
c-----Ef�̏����l�̐ݒ�-----
	ini_fermi = -0.30
	ini_tel   = 300
	if((avsumconc(ix,iz,iv).ge.5e24)
     &.and.(avsumconc(ix,iz,iv).le.1e25))then
		ini_fermi = 0.30
	elseif((avsumconc(ix,iz,iv).ge.1e24)
     &.and.(avsumconc(ix,iz,iv).lt.5e24))then
		ini_fermi = 0.20
	elseif((avsumconc(ix,iz,iv).ge.5e23)
     &.and.(avsumconc(ix,iz,iv).lt.1e24))then
		ini_fermi = 0.05
	endif
c
	goto 1230
c
c-----InP-----
   15 continue
		write(99,*)'InP�ɓ���܂���'
		stop
c
c
c-----channel�ޗ�-----
   16 continue
	if((avsumconc(ix,iz,iv).gt.0.0).and.
     &(avsumconc(ix,iz,iv).le.1e18))then
		avsumconc(ix,iz,iv) = 1e18
	elseif(avsumconc(ix,iz,iv).ge.1e25)then
		avsumconc(ix,iz,iv) = 1e25
	endif
c
c-----���ϔZ�x�ƕ��σG�l���M�[�̕␳-----
	do i1=0,100000,10
		if(avsumconc(ix,iz,iv).le.ntab1(i1,5))then
			avsumconc(ix,iz,iv) = ntab1(i1,5)
			exit
		endif
	enddo
c
	if(avsumtei1(ix,iz,iv).le.etab1(i1,5))then
		avsumtei1(ix,iz,iv) = etab1(i1,5)
	endif
c
c-----Ef�̏����l�̐ݒ�-----
	ini_fermi = -0.30
	ini_tel   = 300
	if((avsumconc(ix,iz,iv).ge.5e24)
     &.and.(avsumconc(ix,iz,iv).le.1e25))then
		ini_fermi = 0.20
	elseif((avsumconc(ix,iz,iv).ge.1e24)
     &.and.(avsumconc(ix,iz,iv).lt.5e24))then
		ini_fermi = 0.10
	elseif((avsumconc(ix,iz,iv).ge.5e23)
     &.and.(avsumconc(ix,iz,iv).lt.1e24))then
		ini_fermi = 0.05
	endif
c
	goto 1230



c-----���ω����ꂽ�Z�x�ƃG�l���M�[����Ef��Tel�����߂� �|��-----
 1230 continue
	call newton_raphson(ka,ini_fermi,ini_tel,avsumconc,avsumtel,
     &                    avsumtei1,efermi,ix,iz,iv,ind,flag3,am,aff) !��������� �ύX
c
 1201 continue
c
c-----��������Ȃ������Ƃ�-----
	flag3 = 0	!����������̃t���O
	if(ind.ne.0)then
		ini_fermi = ini_fermi + 0.001   !��������� �ύX
		if((ini_fermi.lt.0.85).and.(ini_tel.lt.2500))then
			goto 1230
		endif
		if((ini_fermi.gt.0.85).and.(ini_tel.lt.2500))then
			ini_fermi = -0.30
			ini_tel   = ini_tel + 100
			goto 1230
		endif
	endif
c
c-----NaN�v�Z���ꂽ�Ƃ�-----
	if(isnan(efermi(ix,iz,iv)))then
		ini_fermi = ini_fermi + 0.001   !��������� �ύX
		if((ini_fermi.lt.0.85).and.(ini_tel.lt.2500))then
			goto 1230
		endif
		if((ini_fermi.gt.0.85).and.(ini_tel.lt.2500))then
			ini_fermi = -0.30
			ini_tel   = ini_tel + 100
			goto 1230
		endif
	endif
c
c-----�G���[�o��-----
	if(ind.ne.0)then
		write(3,*)'N=',avsumconc(ix,iz,iv),'E=',avsumtei1(ix,iz,iv)
		write(3,*)'ix=',ix,'iz=',iz,'ka=',ka,'ind=',ind
		write(3,*)'ini_fermi=',ini_fermi
		write(3,*)'Tel=',avsumtel(ix,iz,iv),'Ef=',efermi(ix,iz,iv)
	endif
c
	if(isnan(efermi(ix,iz,iv)))then
		write(3,*)'N=',avsumconc(ix,iz,iv),'E=',avsumtei1(ix,iz,iv)
		write(3,*)'ix=',ix,'iz=',iz,'ka=',ka,'ind=',ind
		write(3,*)'ini_fermi=',ini_fermi
		write(3,*)'Tel=',avsumtel(ix,iz,iv),'Ef=',efermi(ix,iz,iv)
	endif
c--------------------
		enddo
	enddo
c
	return
	end