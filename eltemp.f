c-----�e�Z����Ԃ̕��ϔZ�x�ƕ��σG�l���M�[���v�Z����T�u���[�`��-----
	subroutine eltemp(adkx,adky,adkz,
     &				  jpnum,spnum,p,kp,
     &				  af2,af4,hhm,ec,
     &				  dx,dz,xmax,zmax,iarea,lhet,
     &				  tel,mconc,cn,
c     &				  ka3)
     &				  ka3,tel1,ict,avsumtel,avsumtel1,		!�|�ݕύX
     &				  avsumconc,avsumtei1,sflag,			!�|�ݕύX
     &				  efermi,efermi1,avsumconc1,avsumtei11,tel2,
     &                  epA,epB,epC,u,eg,						!120126homma
     &                  avsumtei1_1,avsumconc_1,avsumteiA,avsumconcA)	!�|�ݕύX

	implicit none
c===����===
c---�ϐ��z��p�p�����[�^---
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c---��{�p�����[�^---
	real(8) pi,q,h,bk
	parameter(pi = 3.141592654, q  = 1.60219e-19)
	parameter(h  = 1.05459e-34, bk = 1.38066e-23)
	real	spnum
	integer	jpnum
c---�f�o�C�X�\��---
	real	dx,dz,xmax,zmax
	integer(1)	iarea(nlayer)
	integer(2)  lhet(nlayer) !�ǉ��@11/06/22�@��
c---�̈�ʃp�����[�^---
	real,	dimension (nvalley,narea)	:: hhm,af2,af4
	real,	dimension (nvalley,narea)	:: ec
c----(���Z�X)---
	integer ii
	integer(2),dimension (nrecess) :: lxr1,lzr1,lxr2,lzr2
c---���q���---
	real,	dimension   (6,npmax)	:: p
	integer(1),dimension (3,npmax)	:: kp
c===���[�J���ϐ�===
	real pdx,pdz
	integer(1)	kv,ken,kl,ka
	real akx,aky,akz
	integer n,eflag
	integer	ix,iz,il,iv
c---kd�Ɋւ���p�����[�^---	
	real adkx(0:nx,0:nz,0:nvalley)
	real adky(0:nx,0:nz,0:nvalley)
	real adkz(0:nx,0:nz,0:nvalley)
	real tkx,tky,tkz
c----(���b�V��&�J�ʓd�q�Z�x)---
	real,	dimension (0:nx,0:nz)	:: cn
	real cnmesh(0:nx,0:nz)
	real kvmesh(0:nx,0:nz,nvalley)
	real kvmesh2(0:nx)  !101221
	real mconc(0:nx,0:nz,nvalley)
	real,	dimension (0:nx,0:nz)	:: epA, epB, epC ,ka4
c
c-----08/8/6 �|��-----
	real avsumtel(0:nx,0:nz,nvalley)
	real avsumconc(0:nx,0:nz,nvalley)
	real avsumconc_1(0:nx,0:nz,nvalley)  !101221�d�q�Z�x�`���l������
	real avsumconcA(0:nx)                !101221�d�q�Z�x�`���l������2
	real(8) avsumtei1(0:nx,0:nz,nvalley)
      real(8) avsumtei1_1(0:nx,0:nz,nvalley)  !101221�d�q�G�l���M�[�`���l������
	real(8) avsumteiA(0:nx)                 !101221�d�q�G�l���M�[�`���l������2
	real efermi(0:nx,0:nz,nvalley)
c
	real,save,allocatable :: sumtel
	dimension :: sumtel(:,:,:)
	real,save,allocatable :: sumconc
	dimension :: sumconc(:,:,:)
	real(8),save,allocatable :: sumtei1
	real(8),save,allocatable :: sumtei2  !101221
	dimension :: sumtei1(:,:,:)
	dimension :: sumtei2(:)
	real,	dimension (0:nx,0:nz) :: efermi1		!���J�̃t�F���~���x��
	real,	dimension(0:nx,0:nz) :: avsumtel1		!���J�̓d�q���x
	real,	dimension(0:nx,0:nz) :: avsumconc1		!���J�̓d�q�Z�x
	real(8),	dimension(0:nx,0:nz) :: avsumtei11	!���J�̓d�q�G�l���M�[
c
	integer	ict,sflag,ict2
	real,	dimension(0:nx,0:nz) :: tel1
	real,	dimension(0:nx,0:nz) :: tel2
c------------------
	real,	dimension (0:nx,0:nz)	:: u	!120126homma
      real	eg(nvalley,narea)				!120126homma
c
	real tel(0:nx,0:nz,nvalley)
	real(8) tei1(0:nx,0:nz,nvalley)
	real(8) tei2(0:nx)  !101221
c
	real(8) ei1,sk,sq
	real bkq23
	integer(1) ka3(0:nx,0:nz,nvalley)
	real(8),save,allocatable :: ei2(:,:,:)

	if (.not. allocated(ei2))then				!�ŏ��������s
		allocate(ei2(0:nx,0:nz,nvalley))
		allocate(sumtel(0:nx,0:nz,nvalley))		!08/8/6 �|��
		allocate(sumconc(0:nx,0:nz,nvalley))	!08/8/6 �|��
		allocate(sumtei1(0:nx,0:nz,nvalley))	!08/8/6 �|��
		allocate(sumtei2(0:nx))	!101221
		ei2 = 0.0
		tel = 0.0
		ka3 = 0.0
		mconc = 0.0
c-----08/8/6 �|��-----
		avsumtel = 0.0
		avsumconc = 0.0
		sumconc = 0.0
		sumtel = 0.0
		sumtei1= 0.0
	    sumtei2= 0.0
		tel1 = 0.0
		avsumtel1 = 0.0
		avsumconc1 = 0.0
		avsumtei11 = 0.0
c--------------------
	endif

	tel2=0.0					!08/10/10 �|��
	kvmesh=0.0 ; cnmesh=0.0		!�K�v
	tei1 = 0.0					!�K�v

	tkx = 0.0 ; tky = 0.0 ; tkz = 0.0
	sk = 0.0

	bkq23 = 2.0/(3.0*bk/q)
	pdx = 1.0/dx
	pdz = 1.0/dz

	do n=1,jpnum
		if(kp(1,n).eq.0)cycle

		ei1 = 0.0
		kv	= kp(1,n)	!�JNo.
		kl  = kp(3,n)	!�wNo.
		ka	= iarea(kl)	!�f��No.

		ix = min(nx,max(0,nint(p(5,n)*pdx)))
		iz = min(nz,max(0,nint(p(6,n)*pdz)))
		ka3(ix,iz,kv) = ka 

		akx	= p(1,n)
		aky	= p(2,n)
		akz	= p(3,n)

		tkx = abs(p(1,n)) - abs(adkx(ix,iz,kv))
		tky = abs(p(2,n)) - abs(adky(ix,iz,kv))
		tkz = abs(p(3,n)) - abs(adkz(ix,iz,kv))

		sk = tkx*tkx + tky*tky + tkz*tkz

		if(af4(kv,ka).ne.0.0)then
			sq = sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*sk)
c			ei1=(sq-1.0)/af2(kv,ka)+ec(kv,ka)
			ei1=(sq-1.0)/af2(kv,ka)					!08/11/4 �|��
c     &			-epA(ix,iz)+u(ix,iz)-0.7*(eg(1,ka)-eg(1,2))		!120126homma
		else
c			ei1=hhm(kv,ka)*sk+ec(kv,ka)
			ei1=hhm(kv,ka)*sk						!08/11/4 �|��
c     &			-epA(ix,iz)+u(ix,iz)-0.7*(eg(1,ka)-eg(1,2))		!120126homma
		endif

c	------------�J��L���W�v-------------
		kvmesh(ix,iz,kv) = kvmesh(ix,iz,kv)+1.0		! ���b�V�����J�ʗ��q��
		cnmesh(ix,iz) = cnmesh(ix,iz)+1.0			! ���b�V�������q���a

c	-----�G�l���M�[Ev(k-kd(r))�̑��a-----
		tei1(ix,iz,kv)=tei1(ix,iz,kv)+ei1
		ei2(ix,iz,kv)= ei1

c     ------�G�l���M�[�`���l������-------
      if(kv.eq.1)then
	  if((lhet(nchannel1)+1.le.iz).and.
     &	(iz.le.lhet(nchannel2))) then !channel 11/06/22��
              tei2(ix)=tei2(ix)+ei1
			kvmesh2(ix)=kvmesh2(ix)+1.0
	  endif
	endif

	enddo

	do iv=1,nvalley
		do iz=0,nz
			do ix=0,nx
c	���Z�X���ւ̏���
				eflag=0
				do ii = 1,nrecess
					if((ix.gt.lxr1(ii)).and.(ix.lt.lxr2(ii))
     &					.and.(iz.ge.lzr1(ii)).and.(iz.lt.lzr2(ii)))then
					tei1(ix,iz,iv)=0.0
					mconc(ix,iz,kv) = 0.0
					eflag=1
					exit
					endif
				enddo

			if(eflag.eq.1)cycle
c	���Z�X�O�ł̏���
			if(kvmesh(ix,iz,iv).eq.0.0)then
				tei1(ix,iz,iv)=0.0
				cycle
			endif

			if(cnmesh(ix,iz).ne.0.0) then

				tei1(ix,iz,iv)=tei1(ix,iz,iv)/kvmesh(ix,iz,iv) !�G�l���M�[���q������

				kvmesh(ix,iz,iv) = kvmesh(ix,iz,iv) / cnmesh(ix,iz)		!��L��

				tel(ix,iz,iv) = bkq23 * tei1(ix,iz,iv)			!�J�ʓd�q���x�Z�o(Boltzmann)
				mconc(ix,iz,iv) = cn(ix,iz) * kvmesh(ix,iz,iv)	!�J�ʓd�q�Z�x�Z�o

			endif

c	-----08/8/6 �|�� �Z�x,�G�l���M�[�̃X�e�b�v���ω�-----
		sumconc(ix,iz,iv) = sumconc(ix,iz,iv) + mconc(ix,iz,iv)
		sumtei1(ix,iz,iv) = sumtei1(ix,iz,iv) + tei1(ix,iz,iv)
c	----------------------------------------------

			enddo
		enddo
	enddo

c----------�d�q�G�l���M�[���ω��Z-------
      do ix=0,nx
        tei2(ix)=tei2(ix)/kvmesh2(ix)
	  sumtei2(ix)=sumtei2(ix)+tei2(ix)
      enddo
c
c
c-----08/6/3 �|�� �Z�x,�G�l���M�[�̃X�e�b�v���ω�-----
c	sflag=0�̎��ɂ͂��̂܂܂̒l,sflag=1�̎��ɂ̓X�e�b�v���ς��s��
c	�e�Z����Ԃɂ�����d�q�Z�x�Ɠd�q���x���X�V����
	sflag=0
	ict2=0
	if((mod(abs(ict),sumcnt).eq.0))then
		sflag=1
		do iv=1,nvalley
			do iz=0,nz
				do ix=0,nx
				avsumconc(ix,iz,iv) = sumconc(ix,iz,iv) / sumcnt
				avsumtei1(ix,iz,iv) = sumtei1(ix,iz,iv) / sumcnt
				enddo
			enddo
		enddo

	    avsumconc_1 = avsumconc
          avsumtei1_1 = avsumtei1

		sumconc = 0.0
		sumtei1 = 0.0


          avsumconcA = 0.0   !avsumconcA������101215

	   do iz = lhet(nchannel1)+1,lhet(nchannel2) !�C��11/06/22�� 
	     do ix=0,nx
	          avsumconcA(ix)=avsumconcA(ix)+avsumconc(ix,iz,1)
		 enddo
		     ict2=ict2+1
         enddo  

	    do ix=0,nx
	          avsumconcA(ix)=avsumconcA(ix)/ict2
                avsumteiA(ix) =sumtei2(ix)/sumcnt
		enddo
	          sumtei2 = 0.0
	endif

c-----------------------------------------------------
c
c-----08/8/6 �|�� �o�͊֌W-----
	do iz=0,nz
		do ix=0,nx
			avsumtel1(ix,iz)  = avsumtel(ix,iz,1)	!output�ŗ��p ���J�̂�
			efermi1(ix,iz)    = efermi(ix,iz,1)		!output�ŗ��p ���J�̂�
			avsumconc1(ix,iz) = avsumconc_1(ix,iz,1)	!output�ŗ��p ���J�̂�
			avsumtei11(ix,iz) = avsumtei1_1(ix,iz,1)	!output�ŗ��p ���J�̂�
			tel2(ix,iz)		  = tel(ix,iz,2)		!output�ŗ��p L�J�̓d�q���x
		enddo
	enddo

	return
	end