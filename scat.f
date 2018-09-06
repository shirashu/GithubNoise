c
c===( �U���ߒ��̌v�Z )===
c				call scat(
c     &					 smh(1,ka), af(1,ka),
c     &				     swk(1,1,1,ken,mtp,ka), escat(1,1,ka),
c     &			         iarg(1,ka), iband(1,1,ka),
c     &		             eps(ka), bktq(mtp), dopem(ix,iz),
c     &			         ap(1),ap(2),ap(3), kv,
c     &	                 sk, skx, sky, ei, ef, ie, den,
c     &					 hescat(1,1,ka),hef,hcss,jp,iscat)
      subroutine scat(
     &                am_aniso,aff_aniso,hole_am_aniso,hole_aff_aniso,
     &				smh, af,
     &				iiflag,iiix,iiiz,	!I.I.�␳�p�ǉ� 11/08��
     &                swk, escat, iarg, iband,
     &                eps, bktq, dopem,
     &                akx, aky, akz, kv, kl,
     &                sk, ei, ef, ie, den,
     &                hescat, hef, hcss, jp, iscat,
     &				dx,dz,jpnum,kp,p,pgm,n,tava,
     &				nava,cnava,x,z,count_flag,six,siscat,skv,			!10/05/07 �ǉ�
     &				cxpole1,cxpole2,balis_scat,balis_flag2,balis_n,		!07/03/15 !11/20
     &				ccs,ccs2,cncs,cncs2,allback_scat,back_scat,			!08/1/28
     &				n_scat_p,n_scat_n,dn3,								!09/2/19 �|��
     &				lxr1,lxr2,ec,eg,hm,ka,
     &			     hiXL,u_eff1,u_eff2,u_eff3,ix,iz,Eth_table,ken,u,			!08/2/19 �|�� �U���Ȃ��p
     &				lhet,de,					!120330	����
     &				scatpoint,xx_goal,zz_goal,x_start,z_start,		!120817sato
     &					II_S,dltec)			!120921sato
c
c%%%%%%%%%% ���� 03/04/11 �ύX %%%%%%%%%%
c
c === �ϐ���� ===
c	---	common ---	
c	nvalley ... �J�̍ő吔
c	nscat ... �l������U���@�\�̑���
c	nemax ... �G�l���M�[�X�e�b�v�̍ő�l
c
c	--- ���� ---
c	*** �o���h�\���Ɋւ���p�����[�^ ***
c	smh(iv) ... iv�Ԗڂ̒J�� ��2m*/h
c	af ... ��������o���h�p�����[�^��
c
c	*** �f�o�C�X�����p�����[�^ ***
c	eps		... �ÓI�U�d����
c	bktq	... �{���c�}���t�@�N�^�[Kb*T*q
c	dopem	... �s�����Z�x
c
c	*** �U���Ɋւ���p�����[�^ ***
c	swk(iscat,ie,iv) ... �K�i�����ꂽ�U�����[�g
c	escat(iscat,iv)	... �U�����ɓd�q������(����)�G�l���M�[
c	iarg(iscat) 	... �I��Ԃɂ�����d�q�̕��������߂�p�����[�^
c	iband(iscat,iv)	... �U���ɂ��J�ڐ�̒J
c		(iscat:�U���@�\,ie:�U���O�̃G�l���M�[,iv:�U���O�̒J)
c
c	*** ���q�Ɋւ������ ***
c	akx,aky,akz ...	���q�̔g��
c	kv ... �����J
c	ei ... �U���O�̃G�l���M�[
c	ef ... �U����̃G�l���M�[
c	ie ... �U���O�̃G�l���M�[�ide�ŗʎq���ρj
c	den... �ʎq���덷�ide*(ie+den)=ei�j
c
c	*** ���M�Ɋւ���p�����[�^ ***
c	hescat(iscat,iv) ... �U�����N�����Ƃ��̔��M�G�l���M�[
c	hef 	... �U���ɂ�锭�M�G�l���M�[
c	hcss	... ���M�v�Z�X�C�b�`�i=1:�v�Z�j
c	jp		... �U���O�J
c	iscat	... �������U��No.�i0:�U���Ȃ��j
c
	implicit none
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
	real,	dimension (nvalley)	:: smh,af
	real(8),	dimension (nscat,nemax,nvalley):: swk
	real,	dimension (nscat,nvalley)	:: escat
	integer(1),dimension (nscat)	:: iarg
	integer(1),dimension (nscat,nvalley)	:: iband
	real	dn3								!09/2/19 �|��
	real	eps,bktq
	real	dopem
	real	akx,aky,akz
	integer(1)	kv, kl
	real(8)	sk
	real(8)	ei
	real	ef, den
	integer ie
	real,	dimension (nscat,nvalley) :: hescat
	real	hef
	integer hcss
	integer jp,iscat








c---(�Փ˓d���ǉ�)---nishino_model
	real,	dimension (2,nemax,nvalley)::	swk_temp								!2010NH
	real,	dimension (nscat,nemax,nvalley)::	temp_allswk							!2010NH
	real,	dimension (nscat,nemax,nvalley)::	temp_allswkd						!2010NH
	real(8),	dimension (2,nemax,nvalley)::	swk_ii								!2010NH
	real(8),	dimension (2,nemax,nvalley)::	swk_ii_x							!2010NH
	real(8),	dimension (2,nemax,nvalley)::	swk_ii_y							!2010NH
	real(8),	dimension (2,nemax,nvalley)::	swk_ii_z							!2010NH
	real(8),	dimension (2,nemax,nvalley)::	swk_ii_xyz							!2010NH
	real(8)	eth4,eth4x,eth4y,eth4z,eth4xyz,alen,dtie								!2010NH
	real dx,dz																		!2010NH
	integer ix,iz,iinscat															!2010NH
c	integer	nx,nz																	!2010NH
	integer,dimension (0:nx,0:nz,0:nvalley) :: nava									!2010NH
	integer,dimension (0:nx,0:nvalley) :: cnava										!2010NH
	real, dimension	(nvalley,narea,nvalley)::am_aniso								!2010NH
	real, dimension	(nvalley,narea,nvalley)::aff_aniso								!2010NH
	real, dimension	(nvalley,narea,nvalley)::hole_am_aniso							!2010NH
	real, dimension	(nvalley,narea,nvalley)::hole_aff_aniso							!2010NH
	real, dimension	(narea,nvalley)::hiXL											!2010NH
	real, dimension (nx,nz)	::	u_eff1,u_eff2,u_eff3,u											!2010NH
	real, dimension	(4,20000,2) :: Eth_table										!2010NH
	real(8), dimension (10) :: test_scat_for_ii										!2010NH
	real(8)		akxyz,d1,d2,d3,d4													!2010NH
	real(8)		tmd1,tmd2,tmd3,tmd4													!2010NH
c---( �Փ˓d��_�^����-�G�l���M�[�ۑ��^ nishino_model 2009�N1��18��)�̕ϐ�---			!2010NH
	real(8) iihhm(3), iihm(3), iipaf2(3), iiaf4(3), iieg							!2010NH
      real(8) iiam(3),iiec(3),iiaf(3)													!2010NH
      real(8) iieps,iispf,iiep														!2010NH
      real(8) iismh(3)														!2010NH
	real(8)	iiegmin(nvalley,narea)
	real(8) ii_effc_pot																!2010NH
	real(8) ii_normalize															!2010NH
	real(8) ii_eth, ii_normal, ii_scat_ratio										!2010NH
      integer iv,ethx,ethy,ethz														!2010NH
	integer(1) ka																	!2010NH
	integer(1) ken																	!2010NH
!	real(8) pi,q,h,bk,ep0,am0														!2010NH
	real(8) am0																		!2010NH
!	parameter(pi  = 3.141592654, q   = 1.60219e-19)									!2010NH
	parameter(am0 = 9.10953e-31)													!2010NH
!	real(8)	ei,sk																	!2010NH
!�ǉ���																				!2010NH
      real(4) iidt																	!2010NH
c	real(4) highestX,highestL														!2010NH
      real(8) skvg,skvlx																!2010NH
      real(8) skinit,skfinal,skhole													!2010NH
      real(8) einit,efinal,ehole														!2010NH
      real(8) amhole,affhole,af4hole,paf2hole,hhmhole									!2010NH
      real(8) skraitio																!2010NH
!      real(8) akx,aky,akz															!2010NH
      real(8) evlg(3)																	!2010NH
	real,	dimension (nvalley,narea)	:: hm										!2010NH
	real,	dimension (nvalley,narea) :: ec											!2010NH
	real	eg(nvalley,narea)														!2010NH
	real	dltec(nvalley,narea)			!2017/12/1 ��؋M��

	integer(2)	lhet(nlayer)			!120330sato			!�w�e���E�ʈʒu
	real	de(nenergy)			!120330sato		�G�l���M�[�X�e�b�v
	
	real II_S(narea)		!120921sato


ccccccc���W�F�N�V�����E�Փ˓d���␳�p 11/08��
	integer iiflag
	integer iiix,iiiz
ccccccccccccccccccccccccccccccccccccccccccccc

c	real,save,allocatable	::	u_eff1	!2006/12/09 Hara
c	dimension	::	u_eff1(:,:)
c	allocate u_eff1(0:nx,0:nz)

	real    tava
	integer   n,m1
c
	integer jpnum
	real    ,dimension   (6,npmax)	:: p
	real    ,dimension   (nvalley,nenergy) :: pgm
	integer(1)  kp(3,npmax)
	real	x,z		!07/03/15
c-----�U���J�E���g-------------------!10/05/07
	integer	  count_flag,six,siscat,skv
c----�o���X�e�B�b�N�̌v�Z----------
	real	cxpole1(npole),cxpole2(npole)		!07/11/22
	integer(1)	balis_flag2				!07/11/22
	integer		balis_scat,balis_n					!07/11/22
c
	integer,dimension (0:nx,nscat,nvalley) :: n_scat_p !(+)yama071223
	integer,dimension (0:nx,nscat,nvalley) :: n_scat_n !(-)yama071223
c-----�U���p�̏W�v-------------------!08/1/21
	real(8),dimension (0:nx) ::	ccs,cncs
	real(8),dimension (0:nx,nscat) ::	ccs2,cncs2
c-----����U���̏W�v-------------------!08/1/28
	real(8),dimension (0:nx,0:nscat) ::	allback_scat
	real(8),dimension (0:nx,0:nscat) ::	back_scat
	integer back_scat_flag
c
c----(���Z�X)---
	integer(2),dimension (nrecess) :: lxr1,lxr2
c
c-----���[�J���ϐ�-----------------------
c	swkd(iscat) ... ���`�⊮�����U�����[�g
	real,	dimension (nscat)::	swkd
	real,	dimension (nscat)::	tempswkd
	real	swkdn
	real	ki,kf,skk,sf,cf,sb,qd2,r2,cb,f,fai,sn,cs,r1
	real(8)	sq
	real	x1,x2,x3,a11,a12,a13,a21,a22,a23,a32,a33
c
	real	pi,q,h
      parameter(pi=3.14159, q=1.60219e-19, h  = 1.05459e-34)
	real	rnd

c-----120817sato
	integer	scatpoint
	real	xx_goal
	real	zz_goal
	real	x_start(npmax)		!120817sato	�Փ˓d�����q�����p
	real	z_start(npmax)
c
	hef = 0.0
	iscat = 0
	ef = ei
	back_scat_flag = 0	!08/1/28
c
	jp = kv
c
	r1 = rnd()
c
ccccc���W�F�N�V�����E�Փ˓d�� �␳�p ������11/08��
	iiix = 0
	iiiz = 0
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	--- �������̂��߁A���炩�ɎU���������Ȃ��Ƃ��ɂ�return ---
	if(r1.gt.swk(nscat,ie,kv))then
		if(ie.eq.1)then
			return
		else
			if(r1.gt.swk(nscat,ie-1,kv))then
				return
			endif
		endif
	endif
c	----------------------------------------------------------
c
c	-----�`���l���w��񃊃Z�X�������͎U���Ȃ�-----
	if(scatsw.eq.0)then							!09/2/19 �|�� �U��ON/OFF�X�C�b�`
		ix = nint(x/dx)
c		if((nlayer-4.lt.kl).and.(kl.le.nlayer-1))then	!�`���l����
		if(nchannel1.lt.kl.and.kl.le.nchannel2) then !channel 11/04/08��
			if((ix.ge.lxr1(2)).and.(ix.lt.lxr2(2))) then
				return
			endif
		endif
	endif
c	------------------------------------------------
c


c---( �Փ˓d�� )--!20100624 nishino_model start!!!!!!!!!!!!!							!2010NH
c	�Փ˓d����Eth�ˑ��̎U�����[�g�̃e�[�u����										!2010NH
																					!2010NH
	akxyz = sqrt(akx**2.0d0 + aky**2.0d0 + akz**2.0d0)								!2010NH
																					!2010NH
	tmd1 = sqrt((abs(akx)-abs(akxyz))**2.0d0+aky**2.0d0+akz**2.0d0)					!2010NH
	tmd2 = sqrt(akx**2.0d0+(abs(aky)-abs(akxyz))**2.0d0+akz**2.0d0)					!2010NH
	tmd3 = sqrt(akx**2.0d0+aky**2.0d0+(abs(akz)-abs(akxyz))**2.0d0)					!2010NH
	tmd4 = sqrt( ( abs(akx)-(abs(akxyz)/sqrt(3.0d0)) )**2.0d0+						!2010NH
     &	( abs(aky)-(abs(akxyz)/sqrt(3.0d0)) )**2.0d0+								!2010NH
     &    ( abs(akz)-(abs(akxyz)/sqrt(3.0d0)) )**2.0d0)								!2010NH	
																					!2010NH
																					!2010NH
	if(tmd1==0.or.tmd2==0.or.tmd3==0.or.tmd4==0)then								!2010NH
	tmd1 = 1.0d0																	!2010NH
	tmd2 = 1.0d0																	!2010NH
	tmd3 = 1.0d0																	!2010NH
	tmd4 = 1.0d0																	!2010NH
	end if																			!2010NH
																					!2010NH
																					!2010NH
	d1 = 1.0d0/tmd1																	!2010NH
	d2 = 1.0d0/tmd2																	!2010NH
	d3 = 1.0d0/tmd3																	!2010NH
	d4 = 1.0d0/tmd4																	!2010NH
																					!2010NH
!	write(*,*) 'scat'																!2010NH
!	if(ka==3.or.ka==5)then	!�R���|�W�b�g�`���l���p									!2010NH
c	if(ka==3)then			!�V���O���`���l���p
!	if(ka==1)then			!�V���O���`���l���p			!11/07/04��
																					!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
																					!2010NH	
c	if(ka==3)then		!!!															!2010NH
	if(ka==1)then			!11/07/04��											
																					!2010NH
c	if(iz >= 55.or.ix > 230.or.ix == 0.or.iz == 0)then								!2010NH
c		if(iz >= 50.or.ix > 500.or.ix == 0.or.iz == 0)then
	 if(iz >= (lhet(nchannel2)+5).or.ix > nx.or.ix == 0.or.iz == 0)then		!120413sato				
			!�Փ˓d���v�Z�͈�		�`���l�����E�ʂ���+5(���q�̍L���蕪)
			!						�f�o�C�X�傫���O			������		 !120330sato

		  iiegmin(kv,1) = eg(kv,1) 												!2010NH
c	 else																			!2010NH	
c		  iiegmin = eg(kv,ka) - ec(kv,ka) 												!2010NH
c    &		+ (u(ix,iz)-0.7*(eg(1,ka)-eg(1,2)) ) - u_eff1(ix,iz)						!2010NH	
	 else
		if(kv.eq.1)then																		!2010NH	
 		  iiegmin(1,1)=eg(1,1) 
     &	   +(u(ix,iz)+dltec(1,1))-u_eff1(ix,iz)						!2010NH	�@c 2017/12/1 ���
		elseif(kv.eq.2)then
			iiegmin(2,1)=eg(2,1) 
     &	   +(u(ix,iz)+dltec(1,1)-ec(2,1))-u_eff2(ix,iz)
		elseif(kv.eq.3)then
			iiegmin(3,1)=eg(3,1) 
     &	   +(u(ix,iz)+dltec(1,1)-ec(3,1))-u_eff3(ix,iz)
		endif
	 end if																		!2010NH	
																					!2010NH
																					!2010NH
																					!2010NH
       ethx=int(iiegmin(kv,ka)*10000)	!eg 0.170														!2010NH
	 if(ethx<=18000)then																!2010NH
        eth4x = Eth_table(2,ethx,1)														!2010NH
	  eth4y = Eth_table(2,ethx,1)														!2010NH
	  eth4z = Eth_table(3,ethx,1)														!2010NH
	  eth4xyz = Eth_table(4,ethx,1)													!2010NH
	 else																			!2010NH
 	  eth4x = iiegmin(kv,ka)*1.2d0															!2010NH
 	  eth4y = iiegmin(kv,ka)*1.2d0															!2010NH
 	  eth4z = iiegmin(kv,ka)*1.2d0															!2010NH
 	  eth4xyz = iiegmin(kv,ka)*1.2d0															!2010NH
	 end if																			!2010NH
	  if(iiegmin(kv,ka)<=0)then																!2010NH
 		  eth4x = 1000																	!2010NH
 		  eth4y = 1000																	!2010NH
 		  eth4z = 1000																	!2010NH
 		  eth4xyz = 1000																	!2010NH																					
	  end if																			!2010NH
																					!2010NH
																					!2010NH																					
																					!2010NH																					
c	  if(ka==5)then		!!!!														!2010NH
c																					!2010NH
c																					!2010NH
c		if(iz >= 55.or.ix > 230.or.ix == 0.or.iz == 0)then								!2010NH
c		  iiegmin = eg(kv,ka) - ec(kv,ka)													!2010NH
c		else																			!2010NH	
c		  iiegmin = eg(kv,ka) - ec(kv,ka) 												!2010NH
c     &		+ (u(ix,iz)-0.7*(eg(1,ka)-eg(1,2)) ) - u_eff1(ix,iz)						!2010NH	
c		end if																			!2010NH
c																					!2010NH
c																					!2010NH
c																					!2010NH
c		ethx=int(iiegmin*10000)															!2010NH
c		if(ethx<=18000)then																!2010NH
c		  eth4x = Eth_table(2,ethx,2)														!2010NH
c		  eth4y = Eth_table(2,ethx,2)														!2010NH
c		  eth4z = Eth_table(3,ethx,2)														!2010NH
c		  eth4xyz = Eth_table(4,ethx,2)													!2010NH
c		else																			!2010NH
c 		  eth4x = iiegmin*1.2d0															!2010NH
c		  eth4y = iiegmin*1.2d0															!2010NH
c 		  eth4z = iiegmin*1.2d0															!2010NH
c 		  eth4xyz = iiegmin*1.2d0															!2010NH
c		end if																			!2010NH
c		if(iiegmin<=0)then																!2010NH
c 		  eth4x = 1000																	!2010NH
c 		  eth4y = 1000																	!2010NH
c 		  eth4z = 1000																	!2010NH
c 		  eth4xyz = 1000																	!2010NH																					
c		end if																			!2010NH
c																					!2010NH
c																					!2010NH																					
c	  end if				!!!!														!2010NH																					
																					!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
c	  dtie = 0.000698567																!2010NH
	  dtie = de(ken)		!120330sato	�U�����[�g�G�l���M�[�X�e�b�v

																					!2010NH
	  alen=(ei+iiec(kv))																!2010NH
																					!2010NH
	  if (alen.gt.eth4x) then															!2010NH
		swk_ii_x(1,ie,kv) = 															!2010NH
     &	  II_S(ka)*((abs(alen-eth4x)/eth4x)**2.0d0)									!2010NH
		swk_ii_x(1,ie-1,kv) = 														!2010NH
     &	  II_S(ka)*((abs(alen-dtie-eth4x)/eth4x)**2.0d0)								!2010NH
	  else 																			!2010NH
		swk_ii_x(1,ie,kv)=0.0
		if(ie.gt.1)then															!2010NH
			swk_ii_x(1,ie-1,kv)=0.0
		else
c			write(*,*) ie-1,kv
		endif
																!2010NH
	  end if																			!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
	if (alen.gt.eth4y) then															!2010NH
	  swk_ii_y(1,ie,kv) = 															!2010NH
     &	  II_S(ka)*((abs(alen-eth4y)/eth4y)**2.0d0)									!2010NH
		swk_ii_y(1,ie-1,kv) = 														!2010NH
     &	  II_S(ka)*((abs(alen-dtie-eth4y)/eth4y)**2.0d0)								!2010NH
	else 																			!2010NH
	  swk_ii_y(1,ie,kv)=0.0	
	  if(ie.gt.1)then														!2010NH
	      swk_ii_y(1,ie-1,kv)=0.0
	  else
c	      write(*,*) ie-1,kv
	  end if														!2010NH
	end if																			!2010NH
																					!2010NH
																					!2010NH
	if (alen.gt.eth4z) then															!2010NH
	  swk_ii_z(1,ie,kv) = 															!2010NH
     &	  II_S(ka)*((abs(alen-eth4z)/eth4z)**2.0d0)									!2010NH
		swk_ii_z(1,ie-1,kv) = 														!2010NH
     &	  II_S(ka)*((abs(alen-dtie-eth4z)/eth4z)**2.0d0)								!2010NH
	else 																			!2010NH
	  swk_ii_z(1,ie,kv)=0.0	
	  if(ie.gt.1)then														!2010NH
	      swk_ii_y(1,ie-1,kv)=0.0
	  else
c	      write(*,*) ie-1,kv
	  end if													!2010NH														!2010NH
	end if																			!2010NH
																					!2010NH
																					!2010NH
	if (alen.gt.eth4xyz) then														!2010NH
	  swk_ii_xyz(1,ie,kv) = 														!2010NH
     &	  II_S(ka)*((abs(alen-eth4xyz)/eth4xyz)**2.0d0)								!2010NH
		swk_ii_xyz(1,ie-1,kv) = 													!2010NH
     &	  II_S(ka)*((abs(alen-dtie-eth4xyz)/eth4xyz)**2.0d0)							!2010NH
	else 																			!2010NH
	  swk_ii_xyz(1,ie,kv)=0.0
	  if(ie.gt.1)then														!2010NH
	      swk_ii_y(1,ie-1,kv)=0.0
	  else
c	      write(*,*) ie-1,kv
	  end if														!2010NH													!2010NH
	end if																			!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
      swk_ii(1,ie,kv) =																!2010NH
     & (swk_ii_x(1,ie,kv)*d1 + 														!2010NH
     &  swk_ii_y(1,ie,kv)*d2 + 														!2010NH
     &  swk_ii_z(1,ie,kv)*d3 +														!2010NH
     &  swk_ii_xyz(1,ie,kv)*d4														!2010NH
     &     )/(d1+d2+d3+d4)															!2010NH
																					!2010NH
																					!2010NH	
      if(ie.gt.1)then
	swk_ii(1,ie-1,kv) = 															!2010NH
     & (swk_ii_x(1,ie-1,kv)*d1 + 														!2010NH
     &  swk_ii_y(1,ie-1,kv)*d2 + 														!2010NH
     &  swk_ii_z(1,ie-1,kv)*d3 +														!2010NH
     &  swk_ii_xyz(1,ie-1,kv)*d4														!2010NH
     &     )/(d1+d2+d3+d4)
      end if															!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
	do iscat=1,nscat																!2010NH
	temp_allswkd(iscat,ie,kv) = 													!2010NH
     &		swk(iscat,ie,kv)
      if(ie.gt.1)then														!2010NH
	temp_allswkd(iscat,ie-1,kv) = 													!2010NH
     &		swk(iscat,ie-1,kv)														!2010NH
	endif
	end do																			!2010NH
																					!2010NH
																					!2010NH
	temp_allswkd(12,ie,kv) 															!2010NH
     &				= swk(11,ie,kv) + swk_ii(1,ie,kv)*pgm(kv,ken)					!2010NH
	if(ie.gt.1)then
	temp_allswkd(12,ie-1,kv) 														!2010NH
     &			= swk(11,ie-1,kv) + swk_ii(1,ie-1,kv)*pgm(kv,ken)					!2010NH
	endif
																					!2010NH
																					!2010NH
	if(temp_allswkd(12,ie,kv) >1.0d0)then											!2010NH
	write(*,*) 'swk_rate_over'
	write(99,*)	'swk_rate_over'										!2010NH
	write( *,*) 'scat12',temp_allswkd(12,ie,kv),swk(11,ie,kv)						!2010NH
     	write( *,*)	swk_ii(1,ie,kv),pgm(kv,ken)											!2010NH
	write( *,*) 'Eth',eth4x,eth4z,eth4xyz,iiegmin									!2010NH
	write( *,*) 'iiegmin',iiegmin,eg(kv,ka)-ec(kv,ka) 								!2010NH
	temp_allswkd(12,ie,kv) = swk(12,ie,kv)
	if(ie.gt.1)then																			!2010NH
	temp_allswkd(12,ie-1,kv) = swk(12,ie-1,kv)										!2010NH
	end if
																					!2010NH
	end if																			!2010NH
																					!2010NH
																					!2010NH	
!	if(swk_ii(1,ie,kv).ne.0)then													!2010NH
!	swk_ii(1,ie,kv) = swk_ii(1,ie,kv)												!2010NH
!	write(*,*) swk_ii(1,ie,kv)														!2010NH
!	end if																			!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
c	---	�U�����[�g�̐��`�⊮	I.I.��p ---										!2010NH
	if(ie.eq.1) then																!2010NH
		swkd(1:nscat)=temp_allswkd(1:nscat,ie,kv)									!2010NH
																					!2010NH	
	else																			!2010NH
																					!2010NH
		swkd(1:nscat)=den*temp_allswkd(1:nscat,ie,kv)+								!2010NH
     &					(1-den)*temp_allswkd(1:nscat,ie-1,kv)						!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
		if(swkd(nscat).gt.1.0)then		!�G���[����									!2010NH
c			write( *,*)'swkd��1.0�𒴂��܂���',swkd,iscat,ie,den,kv					!2010NH
c     			write( *,*)	swk_ii(1,ie,kv)*pgm(kv,ken),pgm(kv,ken)						!2010NH
c			write(99,*)'swkd��1.0�𒴂��܂���',swkd,iscat,ie,den,kv					!2010NH
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	InAs�̏ꍇ�C�U�����[�gswk�̍ő�lgm=maxval(swk)���Փ˓d���̔������Ȃ��G�l���M�[�l
c														(getswk�Q��)
c	���̂��߁C�Փ˓d���̃��[�g���v�Z������(swkd)�ƍő�lgm�𒴂��Ă��܂��ꍇ������D	  
c	
c	�ȈՑ΍�F�U�����[�g1�𒴂����l�ōċK�i�� (���s�Q��)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 120413sato
			swkdn = swkd(nscat)														!2010NH
			swkd(1:nscat)=swkd(1:nscat)/swkdn		!�ċK�i��						!2010NH
		endif																		!2010NH
	endif																			!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
	else					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!					!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
																					!2010NH	
	iiegmin = eg(kv,ka) - ec(kv,ka)													!2010NH
																					!2010NH
																					!2010NH
c	---	�U�����[�g�̐��`�⊮	 ---												!2010NH
	if(ie.eq.1) then																!2010NH
		swkd(1:nscat)=swk(1:nscat,ie,kv)											!2010NH
																					!2010NH
	else																			!2010NH
		swkd(1:nscat)=den*swk(1:nscat,ie,kv)+										!2010NH
     &					(1-den)*swk(1:nscat,ie-1,kv)								!2010NH
																					!2010NH
		if(swkd(nscat).gt.1.0)then		!�G���[����									!2010NH
c			write( *,*)'swkd��1.0�𒴂��܂���',swkd,iscat,ie,den,kv					!2010NH
c			write(99,*)'swkd��1.0�𒴂��܂���',swkd,iscat,ie,den,kv					!2010NH
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	InAs�̏ꍇ�C�U�����[�gswk�̍ő�lgm=maxval(swk)���Փ˓d���̔������Ȃ��G�l���M�[�l
c														(getswk�Q��)
c	���̂��߁C�Փ˓d���̃��[�g���v�Z������(swkd)�ƍő�lgm�𒴂��Ă��܂��ꍇ������D	  
c	
c	�ȈՑ΍�F�U�����[�g1�𒴂����l�ōċK�i�� (���s�Q��)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 120413sato

			swkdn = swkd(nscat)														!2010NH
			swkd(1:nscat)=swkd(1:nscat)/swkdn		!�ċK�i��						!2010NH
		endif																		!2010NH
	endif																			!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
	endif					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!					!2010NH




















c	--------------------------------
c
c	--- �U�����[�g�Ɨ�������U��No.(iscat)���� ---
!SMP$	ASSERT(ITERCNT(1))
	do iscat=1,nscat
		if(r1.le.swkd(iscat))then
			goto 1000
		endif
	enddo
c	----------------------------------------------
c
c	--- ���[�v�𔲂��遨�U���Ȃ�(����I��)  ---
	hef = 0.0
	iscat = 0
	ef = ei
	scatpoint = 0	!120817sato
	return
c
c
c	### �U���m��(iscat)��̏�Ԍ��� ####
 1000	continue
c
	xx_goal = x	!120817sato
	zz_goal = z	!120817sato
	scatpoint = 1
	 
c	---(�`���l�����U���p�x)---
c	ix = nint(p(5,n)/dx)
c	iz = nint(p(6,n)/dz)
	ix = nint(x/dx)	!07/03/15
	iz = nint(z/dz)	!07/03/15
c	if((iscat.ne.0).and.(nlayer-4.lt.kl.and.kl.le.nlayer-1)) then	!channel 2006/12/22
	if((iscat.ne.0)
     & .and.(nchannel1.lt.kl.and.kl.le.nchannel2))then !channel 11/04/08��			
			count_flag = 1			!09/12/11 �U���J�E���g
			six = ix
			siscat = iscat
			skv = kv
c-----08/8/6 �|��-----
		if(akx.ge.0.0) then
			n_scat_p(ix,iscat,kv)=n_scat_p(ix,iscat,kv)+1
		else
			n_scat_n(ix,iscat,kv)=n_scat_n(ix,iscat,kv)+1
		endif
c---------------------
	endif
c -----------(�o���X�e�B�b�N�̌v�Z-------- 07/11/20 ��[
	if((iscat.ne.0).and.
     &	((cxpole1(2).lt.x).and.(cxpole2(2).gt.x))) then
c	iscat=0�ł͂Ȃ��A�d�ɉ��̗��q�Ȃ�B
		balis_scat = balis_scat + 1
		if(balis_flag2.ne.2)then
			balis_n = balis_n+1
			balis_flag2 = 2
		endif
	endif
c------------------------------------------------------
c-----����U���̏W�v-------------------!08/1/28
	if(akx.gt.0)then
		back_scat_flag = 1
	endif
c--------------------------------------------
c
c	---	( ���M�G�l���M�[ )	  ---
	if(hcss.eq.1) hef = hescat(iscat,kv)
c	---
c
c		!������sk,ei�͈����Ŏ󂯎�����l
	ef=ei+escat(iscat,kv)	!�U���㗱�q�G�l���M�[
	if(ef.lt.0.0)then
c		write( *,*)'�U����G�l���M�[�l�����ł�(scat)',sq
c		write(99,*)'�U����G�l���M�[�l�����ł�(scat)',sq
	hef = 0.0
	iscat = 0
	ef = ei

	return		!�U���Ȃ�
	endif
c
	kv = iband(iscat,kv)	!�J���U���ŕω�
c
c---( �U����̔g���̐�Βlkf���� )---
	if(af(kv).eq.0.0)then
		kf = smh(kv)*sqrt(ef)
	else
		sq = ef*(1.+af(kv)*ef)
		if(sq.ge.0.0)then
			kf = smh(kv)*sqrt(sq)
		else
			write( *,*)'sq�̒l���s���ł�(scat)',sq
			write(99,*)'sq�̒l���s���ł�(scat)',sq
			stop
		endif
	endif
c
	ki = sqrt(sk)				!�U���O�̔g���̐�Βl
c
c---( �U���̎�ނ��ƂɎU�����Ԍ��� )---
	select case (iarg(iscat))
c
c---( �������U��(��L�ɐ����w�t�H�m���U��,�����t�H�m���U��,�����U��) )---
	case(1)
		cs  = 1.0-2.0*rnd()
c
		if(abs(cs) .ge. 1.0)then
		  sn = 0.0
		else
		  sn  = sqrt(1.0-cs*cs)
		endif
c
		fai = 2.0*pi*rnd()
		akx  = kf*cs
		aky  = kf*sn*cos(fai)
		akz  = kf*sn*sin(fai)

c-----�U���p�����߂�------------- !07/1/21 ��������emcd2��
c		if((iscat.ne.0).and.(nlayer-4.lt.kl.and.kl.le.nlayer-1)) then	!�`���l���� 2006/12/22	
	if((iscat.ne.0)
     & .and.(nchannel1.lt.kl.and.kl.le.nchannel2))then !channel 11/04/08��
			cncs(ix) = cncs(ix) + 1.0
			ccs(ix) = ccs(ix) + cs
			cncs2(ix,iscat) = cncs2(ix,iscat) + 1.0
			ccs2(ix,iscat) = ccs2(ix,iscat) + cs
	endif
c-------------------------------
c-----����U�����J�E���g------------- !07/1/28
c	if((iscat.ne.0).and.(nlayer-4.lt.kl.and.kl.le.nlayer-1)) then  !�`���l���� 2006/12/22
c	if((iscat.ne.0).and.(kl.gt.nlayer-3).and.(kl.le.nlayer-2)) then	 !�`���l���� 11/03/23��	
	if((iscat.ne.0)
     & .and.(nchannel1.lt.kl.and.kl.le.nchannel2))then !channel 11/04/08��
			allback_scat(ix,0) = allback_scat(ix,0) +1.0	
			allback_scat(ix,iscat) = allback_scat(ix,iscat) +1.0
		if((akx.lt.0).and.(back_scat_flag.eq.1))then
			back_scat(ix,0) = back_scat(ix,0) + 1.0
			back_scat(ix,iscat) = back_scat(ix,iscat) + 1.0
		endif
	endif
c--------------------------------
	return
c
c---( �L�ɐ����w�t�H�m���U�� )---
	case(2)
	  f = 2.0*ki*kf/(ki-kf)/(ki-kf)
	  if (f.le.0.0) then
	    write(99,*)'f���s���ł�(scat)',f,ki,kf
	    write( *,*)'f���s���ł�(scat)',f,ki,kf
	    return
	  endif
	  cb  = (1.0 + f - (1.0+2.0*f)**rnd())/f
	goto 20
c
c---( �s�����U�� )---
	case(3)
	  r2=rnd()
c	  qd2=q*dopem/eps/bktq
	  qd2=q*dn3/eps/bktq		!09/2/19 �|��
	  cb=1.0-r2/(0.5+(1.0-r2)*2*sk/qd2)
	goto 20



























c
c---( �Փ˓d�� )--!20100624 nishino_model start!!!!!!!!!!!!!							!2010NH
	case(4)																			!2010NH										
																					!2010NH																					
c	goto 720																		!2010NH											
																					!2010NH	
c	write (*,*) '�Փ˓d���������I'													!2010NH									
c																					!2010NH										
c	---( �I��Ԃ̈ʒu���� )---														!2010NH										
																					!2010NH	
!---( �֐��ѕ� )---																	!2010NH									
c	iiegmin = 0.354																	!2010NH									
																					!2010NH									
!	if(iz >= 55.or.ix > 230.or.ix == 0.or.iz == 0)then								!2010NH	
!		ii_effc_pot = 0d0															!2010NH								
!	else																			!2010NH								
!		ii_effc_pot = u_eff1(ix,iz)													!2010NH								
!	end if																			!2010NH								
!	iiegmin = eg(kv,ka) - ec(kv,ka) 												!2010NH
!    &	+ (u(ix,iz)-0.7*(eg(1,ka)-eg(1,2)) ) - u_eff1(ix,iz)						!2010NH								
!																					!2010NH
!	write(*,*) ei,ii_effc_pot,iiegmin												!2010NH
																					!2010NH
	if (ei <= iiegmin(kv,ka))then															!2010NH
c		if(ka==3.and.kv==1)then	!11/07/26��
		if(ka==1.and.kv==1)then								!2010NH
c		write(*,*) 'ii_reject1'														!2010NH
		!�R�����g�A�E�g InAs�΍�
		endif																		!2010NH

	goto 780																		!2010NH
	end if																			!2010NH
																					!2010NH
	iiegmin = eg(kv,ka) - ec(kv,ka)													!2010NH
																					!2010NH
																					!2010NH
!---( �d�q�̗L������ )---																!2010NH
																					!2010NH
      iiam(1) =																		!2010NH
     & (am_aniso(1,ka,1)*d1 + 														!2010NH
     &  am_aniso(1,ka,2)*d2 + 														!2010NH
     &  am_aniso(1,ka,3)*d3 															!2010NH
     &     )/(d1+d2+d3)																!2010NH
																					!2010NH
      iiam(2) =																		!2010NH
     & (am_aniso(2,ka,1)*d1 + 														!2010NH
     &  am_aniso(2,ka,2)*d2 + 														!2010NH
     &  am_aniso(2,ka,3)*d3 															!2010NH
     &     )/(d1+d2+d3)																!2010NH
																					!2010NH
      iiam(3) =																		!2010NH
     & (am_aniso(3,ka,1)*d1 + 														!2010NH
     &  am_aniso(3,ka,2)*d2 + 														!2010NH
     &  am_aniso(3,ka,3)*d3 															!2010NH
     &     )/(d1+d2+d3)																!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
!---( hole�̗L������ )---																!2010NH
																					!2010NH
      amhole =																		!2010NH
     & (hole_am_aniso(1,ka,1)*d1 + 													!2010NH
     &  hole_am_aniso(1,ka,2)*d2 + 													!2010NH
     &  hole_am_aniso(1,ka,3)*d3 														!2010NH
     &     )/(d1+d2+d3)																!2010NH
																					!2010NH
																					!2010NH
																					!2010NH
!---( �`���ђ�̃G�l���M�[ )---														!2010NH
c	iiec(1)	= 0.00		!���J														!2010NH
c	iiec(2)	= 1.0661		!L�J													!2010NH
c	iiec(3)	= 1.5895		!X�J													!2010NH
																					!2010NH																					
!---( �`���ђ�̃G�l���M�[ )---														!2010NH
	iiec(1)	= 0.00		!���J														!2010NH
	iiec(2)	= ec(2,ka)		!L�J										!2010NH
	iiec(3)	= ec(3,ka)		!X�J										!2010NH
																					!2010NH																					
!---( ����������p�����[�^�� )---														!2010NH
																					!2010NH
      iiaf(1) =																		!2010NH
     & (aff_aniso(1,ka,1)*d1 + 														!2010NH
     &  aff_aniso(1,ka,2)*d2 + 														!2010NH
     &  aff_aniso(1,ka,3)*d3 															!2010NH
     &     )/(d1+d2+d3)																!2010NH
																					!2010NH
      iiaf(2) =																		!2010NH
     & (aff_aniso(2,ka,1)*d1 + 														!2010NH
     &  aff_aniso(2,ka,2)*d2 + 														!2010NH
     &  aff_aniso(2,ka,3)*d3 															!2010NH
     &     )/(d1+d2+d3)																!2010NH
																					!2010NH
      iiaf(3) =																		!2010NH
     & (aff_aniso(3,ka,1)*d1 + 														!2010NH
     &  aff_aniso(3,ka,2)*d2 + 														!2010NH
     &  aff_aniso(3,ka,3)*d3 															!2010NH
     &     )/(d1+d2+d3)																!2010NH
																					!2010NH
																					!2010NH
!---( hole ����������p�����[�^��)---													!2010NH
																					!2010NH
      affhole =																		!2010NH
     & (hole_aff_aniso(1,ka,1)*d1 + 													!2010NH
     &  hole_aff_aniso(1,ka,2)*d2 + 													!2010NH
     &  hole_aff_aniso(1,ka,3)*d3 													!2010NH
     &     )/(d1+d2+d3)																!2010NH
																					!2010NH
      do iv=1,nvalley																		!2010NH
																					!2010NH
      iipaf2(iv)	= sngl(0.5D0/dble(iiaf(iv)))										!2010NH
	iiaf4(iv)	= 4.0*iiaf(iv)														!2010NH
	iismh(iv)	= sqrt(2.0*iiam(iv))*sqrt(q)/h										!2010NH
	iihhm(iv)	= h/iiam(iv)/q*h/2.													!2010NH
	iihm(iv)	= h/iiam(iv)   														!2010NH
																					!2010NH	
	end do 																			!2010NH
																					!2010NH	
	!---( hole�̗L������ )---														!2010NH
c	amhole	= 0.57162d0*am0	!���J													!2010NH
																					!2010NH
!---( ����������p�����[�^�� )---														!2010NH
	affhole =(1.0d0-amhole/am0)**2/iiegmin(kv,ka)            !���J		 					!2010NH
      paf2hole=sngl(0.5D0/dble(affhole))												!2010NH
      af4hole =4.0*affhole															!2010NH
      hhmhole	= h/amhole/q*h/2.														!2010NH
																					!2010NH
																					!2010NH
      iv=kv 																			!2010NH
 																					!2010NH   
																					!2010NH
      einit=ei																		!2010NH
																					!2010NH
 																					!2010NH   
      if(iv==1)then																	!2010NH
																					!2010NH
      skinit=sk       																!2010NH
																					!2010NH
      else																			!2010NH
 																					!2010NH   
c      evlg(2)=1.47d0																	!2010NH
c      evlg(3)=2.50d0																	!2010NH
      evlg(2)=hiXL(ka,2)																!2010NH
      evlg(3)=hiXL(ka,3)																!2010NH
      evlg(1)=evlg(iv)																!2010NH
 																					!2010NH   
      skvg=((2.0d0*iiaf(1)*evlg(1)+1.0d0)**2.0d0-1.0d0)								!2010NH
     &	/(iiaf4(1)*iihhm(1))														!2010NH
																					!2010NH     
      skvlx=((2.0d0*iiaf(iv)*(evlg(iv)-iiec(iv))+1.0d0)**2.0d0-1.0d0)								!121031sato
     &	/(iiaf4(iv)*iihhm(iv))														!2010NH
 																					!2010NH       
      skinit=skvg+skvlx-sk															!2010NH
																					!2010NH    
      end if																			!2010NH
																					!2010NH
      iidt=sk/100																		!2010NH
 																					!2010NH       
      do skfinal=0,skinit,iidt														!2010NH
																					!2010NH
        skhole=skinit-2.0d0*skfinal													!2010NH
 																					!2010NH       
      efinal=(sqrt(abs(1.0+iiaf4(1)*iihhm(1)*skfinal))-1.0d0)*iipaf2(1)				!2010NH
      ehole=(sqrt(abs(1.0+af4hole*hhmhole*skhole))-1.0d0)								!2010NH
     &	*paf2hole+iiegmin(kv,ka)															!2010NH
																					!2010NH        
 !       write(6,*) efinal,ehole														!2010NH
 !       write(6,*) skfinal,skhole													!2010NH
 																					!2010NH           
        if(einit<=(2.0d0*efinal+ehole))then											!2010NH
  																					!2010NH      
            skraitio=skfinal/skinit													!2010NH
            sk=skfinal																!2010NH
 																					!2010NH           
            akx=akx*skraitio															!2010NH
            aky=aky*skraitio							!121031sato
            akz=akz*skraitio							!121031sato
            kv=1																		!2010NH
																					!2010NH
             goto 720         														!2010NH
 																					!2010NH   
        end if																		!2010NH
																					!2010NH
      end do																			!2010NH
																					!2010NH	
	write(*,*) 'ii_reject2'															!2010NH

	goto 770																		!2010NH
c	akx=0																			!2010NH
c     aky=0																			!2010NH
c     akz=0																			!2010NH
																					!2010NH    
720   continue																		!2010NH
c																					!2010NH
c																					!2010NH
c	---(���b�V���ʏՓ˓d�������p�x�J�E���g)---										!2010NH
	nava(ix,iz,kv) = nava(ix,iz,kv)+1												!2010NH
	nava(ix,iz,0) = nava(ix,iz,0)+1													!2010NH
c	if(nlayer-4.lt.kl.and.kl.le.nlayer-1) then	!channel 2006/12/22	
	if(nchannel1.lt.kl.and.kl.le.nchannel2)then !channel 11/04/08��
		cnava(ix,kv) = cnava(ix,kv)+1												!2010NH
		cnava(ix,0) = cnava(ix,0)+1													!2010NH
	endif
c
cccccc�Փ˓d���J�E���g�␳�p�@11/08��
	iiix = ix
	iiiz = iz
ccccccccccccccccccccccccccccccccccccc													!2010NH
c																					!2010NH
c	---( �V�������q�̐��� )---														!2010NH
	if(jpnum.ge.npmax) write(*,*)'���q�����������܂�-'								!2010NH
	m1=jpnum+1  !m�Ԗڂ̗��q�ɂ��̎�����											!2010NH
c																					!2010NH
	kp(1,m1)=1																		!2010NH
	kp(2,m1)=1																		!2010NH
	kp(3,m1)=kp(3,n)																!2010NH
c																					!2010NH
	p(1,m1) =akx																	!2010NH
	p(2,m1) =aky																	!2010NH
	p(3,m1) =akz																	!2010NH
c																					!2010NH
cc	taba=p(4,n)  !�A�o�����V�F���N���鎞��											!2010NH
	p(4,m1)=tava-log(rnd())*pgm(kp(1,m1),kp(2,m1))									!2010NH
c																					!2010NH
	p(5,m1)=p(5,n)																	!2010NH
	p(6,m1)=p(6,n)																	!2010NH
c																					!2010NH
	x_start(m1) = p(5,m1)			!1120817sato
	z_start(m1) = p(6,m1)

	
	jpnum=jpnum+1  !���q����
c
cccccc���W�F�N�V�����E�Փ˓d���@�d������m�F�p 11/08�� 
c		write(*,*) 'scat',n,ix,iz,jpnum
c		write(*,*) 'scat',nava(ix,iz,0),cnava(ix,0)
cccccc
cccccc���W�F�N�V�����E�Փ˓d���@�d����� I.I.�t���O���� 11/08��
	iiflag = 1														!2010NH
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc																						!2010NH
c
770   continue																		!2010NH
780	continue																		!2010NH
	return																			!2010NH
c																					!2010NH
	end select																		!2010NH
																					!2010NH
c---( �Փ˓d�� )--!20100624 nishino_model finish!!!!!!!!!!!!!							!2010NH



































c
c---( ���ʊp�̌��� )---
   20	continue
	if((-1.0.le.cb).and.(cb.le.1.0))then
	  sb = sqrt(1.0-cb*cb)
	else
	  write(99,*)'�U���pcb���s���ł�(scat)',cb
	  write( *,*)'�U���pcb���s���ł�(scat)',cb
        sb = 0.0
	endif
c
	fai = 2.0*pi*rnd()
	cf  = cos(fai)
	sf  = sin(fai)

c-----�U���p�����߂�------------- !07/1/21
c		if((iscat.ne.0).and.(nlayer-4.lt.kl.and.kl.le.nlayer-1)) then	!channel 2006/12/22		
	if((iscat.ne.0)
     & .and.(nchannel1.lt.kl.and.kl.le.nchannel2))then !channel 11/04/08��
			cncs(ix) = cncs(ix) + 1.0
			ccs(ix) = ccs(ix) + cb
			cncs2(ix,iscat) = cncs2(ix,iscat) + 1.0
			ccs2(ix,iscat) = ccs2(ix,iscat) + cb
		endif
c-------------------------------
c
c---( ���W���ϊ� )---
	skk = akx*akx+aky*aky
c
	if(skk .ne. 0.0)then
c
	  skk = sqrt(skk)
c
	  a11 = aky/skk
	  a12 = akx*akz/skk/ki
	  a13 = akx/ki
	  a21 =-akx/skk
	  a22 = aky*akz/skk/ki
	  a23 = aky/ki
	  a32 =-skk/ki
	  a33 = akz/ki
c
	else
c
	  a11 = 1.0
	  a12 = 0.0
	  a13 = 0.0
	  a21 = 0.0
	  a22 = 1.0
	  a23 = 0.0
	  a32 = 0.0
	  a33 = 1.0
c
	endif
c
	x1  = kf*sb*cf
	x2  = kf*sb*sf
	x3  = kf*cb
c
	akx = a11*x1+a12*x2+a13*x3
	aky = a21*x1+a22*x2+a23*x3
	akz =        a32*x2+a33*x3
c
c-----����U�����J�E���g------------- !07/1/28
c	if((iscat.ne.0).and.(nlayer-4.lt.kl.and.kl.le.nlayer-1)) then !channel 2006/12/22
	if((iscat.ne.0)
     & .and.(nchannel1.lt.kl.and.kl.le.nchannel2))then !channel 11/04/08��
		allback_scat(ix,0) = allback_scat(ix,0) +1.0	
		allback_scat(ix,iscat) = allback_scat(ix,iscat) +1.0
		if((akx.lt.0).and.(back_scat_flag.eq.1))then
			back_scat(ix,0) = back_scat(ix,0) + 1.0
			back_scat(ix,iscat) = back_scat(ix,iscat) + 1.0
		endif
	endif

	if (ka>5)then
	write(*,*) 'scat ka error ka kv akxyz',ka,kv,akx,aky,akz
	ka = 5
c	read(*,*) 
	end if
	if (kv>3)then
	write(*,*) 'scat ka error ka kv akxyz',ka,kv,akx,aky,akz
	kv = 3
c	read(*,*) 
	end if


c--------------------------------
      return
c
	end