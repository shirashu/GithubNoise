	subroutine emcd2(
     &                 am_aniso,aff_aniso,hole_am_aniso,hole_aff_aniso,
     &				 dt,de,jpnum,
     &                 dx,dz,xmax,zmax,lhet,iarea,twodeg,dopem,
     &                 cxpole1,cxpole2,melpos,jspsum,
     &                 smh,hm,hhm,af,af2,af4,eps,eg,bktq,
     &                 swk,pgm,escat,iarg,iband,
     &                 p,kp,u,cn,hef_mesh,
     &                 mtemp,hescat,hef_scat,hcss,
     &				 cxr1,czr1,cxr2,czr2,lxr1,lzr1,lxr2,lzr2,czpart2,	!07/8/4 �s�����U��
     &				 nava,cnava,n_scat,count_scat,ict,ep_flag,	!Effective Potential�p�ǉ��@Hara 2006/12/09
     &				 balis_flag,balis_scat,balis_n,balis_all,
     &				 ccs,ccs2,cncs,cncs2,allback_scat,back_scat,	!08/1/21
     &				 ec,adkx,adky,adkz,efermi,n_scat_p,n_scat_n,	!08/8/6 �|��
     &			     sscnt,avsumconc,avsumtel,dn3,hiXL,
     &                 epA,epB,epC,epA2,epB2,Eth_table,	!09/2/19 �|�� !120126homma
     &            xxx,vvv,basho_reflection,basho_roughness,
     &            split,delta,lambda,count_reflection,
     &            count_roughness,average,
     &			pass,pass_r,x_start,xx_goal,scatpoint,			!120817sato
     &			x_mean_free_path_sum,x_mean_free_path_count,
     &			z_start,zz_goal,mean_free_path_sum,
     &			  mean_free_path_sum2,bscat_xflag,fix_u,			!15/1/2takahashi
     &				II_S,swk_rou,roughness1_countx,roughness1_counte,		!120921sato	!121029sato		
     &			dltec)

	implicit none
c===����===

c---�ϐ��z��p�p�����[�^---
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c---��{�p�����[�^---
	integer	ict		!2006/12/09 Hara
	real	dt
	integer	jpnum
	integer	ep_flag		!2006/12/09 Hara
c---�f�o�C�X�\��---
	real	dx,dz,xmax,zmax
	integer(2)	lhet(nlayer)
	integer(1)	iarea(nlayer)
	real	twodeg(nlayer)
	real	dopem(0:nx,0:nz)
	real czpart2(npart)	!07/8/4 �s�����U��
	integer ipart,kpart	!07/8/4 �s�����U��
c---�d��---
	real	cxpole1(npole),cxpole2(npole),melpos(npole)
	integer(4)	jspsum(npole)
c---�̈�ʃp�����[�^---
	real,	dimension (nvalley,narea)	:: smh,hhm,hm,af,af2,af4
	real	eps(narea),eg(nvalley,narea),bktq(ntenum)
	real	dltec(nvalley,narea) !band offsets by nextnano
c---�U���p�����[�^---
	real	de(nenergy)
c	real,	dimension (nscat,nemax,nvalley,nenergy,ntenum,narea)::	swk
c	real,	dimension (nvalley,nenergy,narea)	:: pgm
	real,	dimension (nscat,nemax,nvalley,nenergy,ntenum,npart)::swk	!07/8/4 �s�����U��
	real,	dimension (nvalley,nenergy,npart)	:: pgm	!07/8/4 �s�����U��
	real,	dimension (nscat,nvalley,narea)	:: escat
	integer(1),dimension (nscat,nvalley,narea)	:: iband
	integer(1),dimension (nscat,narea)	:: iarg
	real,	dimension (npart)	:: dn3					!09/2/19 �|��
c---���q���---
	real	,dimension   (6,npmax)	:: p
	integer(1),dimension (3,npmax)	:: kp
	real	,dimension   (6)	:: sp				 !121009			
	integer(1),dimension (6)	:: skp				 !121009				
c---�f�o�C�X�����---
	real,	dimension (0:nx,0:nz)	:: u,cn,hef_mesh
	real,	dimension (0:nx,0:nz)	:: epA,epB,epC,epA2,epB2	!120126homma
	real,save,allocatable	::	u_eff1,u_eff2,u_eff3	!2006/12/09 Hara
	dimension	::	u_eff1(:,:),u_eff2(:,:),u_eff3(:,:)
c---���M���p�z��---
	integer(2),dimension (0:nx,0:nz) :: mtemp
	real,	dimension (nscat,nvalley,narea) :: hescat
	real hef
	real,	dimension (0:nx,nscat,nvalley) :: hef_scat
	integer hcss
c----(���Z�X)---
	integer ii
	real, dimension (nrecess) :: cxr1,czr1,cxr2,czr2
	integer(2),dimension (nrecess) :: lxr1,lzr1,lxr2,lzr2
c---(�Փ˓d��)---
	integer,dimension (0:nx,0:nz,0:nvalley) :: nava
	integer,dimension (0:nx,0:nvalley) :: cnava
	integer,dimension (0:nx,nscat,nvalley) :: n_scat
	integer,dimension (0:nx) :: count_scat

	real, dimension	(nvalley,narea,nvalley)::am_aniso		!20100624
	real, dimension	(nvalley,narea,nvalley)::aff_aniso		!20100624
	real, dimension	(nvalley,narea,nvalley)::hole_am_aniso	!20100624
	real, dimension	(nvalley,narea,nvalley)::hole_aff_aniso	!20100624
	real, dimension	(narea,nvalley)::hiXL
c	real(4) highestX,highestL
	real, dimension	(4,20000,2) :: Eth_table
	real II_S(narea)		!120921sato

c----�o���X�e�B�b�N�̌v�Z----------
	integer(1),dimension (0:npmax) ::	balis_flag		!07/11/22
	integer(1) balis_flag2
	integer	balis_scat,balis_n,balis_all
c-----�U���p�̏W�v-------------------!08/1/21
	real(8),dimension (0:nx) ::	ccs,cncs
	real(8),dimension (0:nx,nscat) ::	ccs2,cncs2	
c-----����U���̏W�v-------------------!08/1/28
	real(8),dimension (0:nx,0:nscat) ::	allback_scat
	real(8),dimension (0:nx,0:nscat) ::	back_scat	
c-----�k�ތ���-----
	real(8) q,bk,bkq
	parameter(q  = 1.60219e-19)
	parameter(bk = 1.38066e-23)
	integer,dimension(0:nx,nscat,nvalley) :: n_scat_p !(+)yama071223
	integer,dimension(0:nx,nscat,nvalley) :: n_scat_n !(-)yama071223
	real drikx(0:nx,0:nz,0:nvalley)
	real driky(0:nx,0:nz,0:nvalley)
	real drikz(0:nx,0:nz,0:nvalley)
	real adkz(0:nx,0:nz,0:nvalley)
	real adkx(0:nx,0:nz,0:nvalley)
	real adky(0:nx,0:nz,0:nvalley)
	real subdkx,subdky,subdkz
	integer	ix2,iz2
	integer(1) kv2
	integer ix3
c
	real dri(0:nx,0:nz,0:nvalley)
	real,	dimension (nvalley,narea) :: ec
	real iak,jak
	integer siflag,sfflag,smpf
	integer,dimension (0:nx,nscat,4) :: sscnt
	integer back_scat_flag,fix_u		!14/12/29takahashi
	integer bscat_xflag !14/12/29takahashi
c
	real efermi(0:nx,0:nz,nvalley)
	real avsumconc(0:nx,0:nz,nvalley)	!08/8/6 �|��
	real avsumtel(0:nx,0:nz,nvalley)	!08/8/6 �|��
	integer rflag,swrej
	integer pflag
c
c----���[�J���ϐ�----
	real akx,aky,akz,x,z,t1,tau
	real(8) ts
c	real bkx,bky,bkz,bx,bz,bt1,btau
c	real(8) bp4
	integer jp,iscat,mtp
	real,save,allocatable	::	dhet
	dimension	::	dhet(:)
c
	real	fx,fz
	integer n,iflag
	real	rnd
	real	ef
	integer(1)	kv,ken,kl,ka,kl2,ka2
	integer	ie
	real(8)	sk,ei,sq
	real	den
	integer	ix,iz,iv	!08/8/6 �|��
	character(80) form
	real	pdx,pdz
	integer	nstat,nend
	integer iii
c-----�U���J�E���g-------------------!10/05/07
	integer	  count_flag,six,siscat,skv
c
ccccc���W�F�N�V�����E�Փ˓d�� �␳�p�@11/08��
	integer iiflag		!�ǉ��@��
	integer iiix,iiiz
cccccccccccccccccccccccccccccccccccccccccccc

c------���t�l�X�U��--J.R. Watling���f��(2012�N�t����)-------------------------
      real delta,lambda,average
	integer split,count_reflection,count_roughness,
     &        vvv,xxx 
	integer,dimension (0:nx,2) :: basho_roughness,basho_reflection

c-----------------------------------------------------------------------------------------
c-----	�̈�ʉ߂��闱�q���J�E���g120817sato
	real(8),dimension (10) :: pass,pass_r		!�J�E���^�Cr�̓��W�F�N�V�����΍�
	real	x_start(npmax)			!drift���n�߂�ʒu��ۑ�����z��
	real	x_start_1				!���[�v������ϐ�
	real	x_start_r				!rejection�΍�
	real	xx_goal					!drift���I���ʒu��ۑ�����ϐ�
	integer	scatpoint			!�U���D��U������
	real	x_mean_free_path_sum
	integer x_mean_free_path_count

	real	z_start(npmax)			!drift���n�߂�ʒu��ۑ�����z��
	real	z_start_1				!���[�v������ϐ�
	real	z_start_r				!rejection�΍�
	real	zz_goal					!drift���I���ʒu��ۑ�����ϐ�
	real	mean_free_path_sum
	real,dimension (0:nx,3) :: 	mean_free_path_sum2		!���ώ��R�s��x���W	 1:�J�E���^,2:x�������ώ��R�s��,3:���ώ��R�s��
	integer xcenter				!x_start��xx_goal�̒��S�_

c------���t�l�X�U��--R.P. Joshi���f��(2012�N�H)-------------------------
	real	swk_rou(nemax,nvalley,nenergy,narea)		!121029sato
	integer,dimension (0:nx,narea,2) :: roughness1_countx
	integer,dimension (nemax,narea,2) :: roughness1_counte

c=========================================================================================================
c
c
	if (.not. allocated(dhet)) then	!�ŏ��������s
		swrej=5 !yama���̏����l��0�ȊO�̒l�Ȃ�OK
		allocate (dhet(0:nlayer))
		dhet(0) = -huge(dhet(0))
		dhet(1:nlayer-1) = zmax*float(lhet(1:nlayer-1))/nz	 !�w�e���E�ʈʒu
		dhet(nlayer) = huge(dhet(nlayer))
c
		allocate (u_eff1(0:nx,0:nz),u_eff2(0:nx,0:nz),u_eff3(0:nx,0:nz))	!2006/12/09 Hara
		u_eff1=0.0
		u_eff2=0.0
		u_eff3=0.0
c----------------------�o���X�e�B�b�N�̌v�Z 07/11/22 ��[--------------
		balis_scat = 0	
		balis_all = 0	
		balis_n = 0	
		cncs = 0.0;cncs2 = 0.0
		ccs = 0.0;ccs2 = 0.0 
		allback_scat = 0.0; back_scat = 0.0
c
		do n=1,jpnum
			if((cxpole1(2).lt.p(5,n)).and.(cxpole2(2).gt.p(5,n)))then	
				balis_flag(n) = 1
				balis_all=balis_all + 1
c					if(balis_flag.eq.0)then
c						balis_n=balis_n + 1
c					endif
c			elseif((cxpole1(2).lt.p(5,n)).and.(cxpole2(2).gt.p(5,n)))then
			else
				balis_flag(n) = 0 
			endif 
		enddo
c--------------------------------------------------
c
		siflag = 0 ; sfflag = 0
		smpf = 0
	endif
c
c-----�k�ތ���-----
	drikx=0.0
	driky=0.0
	drikz=0.0
c	adkx=0.0;adky=0.0;adkz=0.0
	dri=0.0
	pflag=0
	subdkx=0.0;subdky=0.0;subdkz=0.0
	ix2=0; iz2=0; kv2=0
c	flag_d=0
c	flag_s=0
c	flag1=0
	iak=0
	bkq=bk/q
c	akcnt=0.0
c
	pdx = 1.0/dx
	pdz = 1.0/dz
c
	if (ep_flag.eq.1) then	
c					!jpot�ɍ��킹�Ď����|�e���V�������v�Z(2006/05/20)2006/12/09
c
c		�ǉ�  �����|�e���V�����v�Z�@2005/12/02 Hara
		call qeffect(dx, dz, lhet, iarea, 
     &					hhm, eg, bktq, lxr1,lzr1,lxr2,lzr2,
     &					u, u_eff1,u_eff2,u_eff3, mtemp,dltec,ec)
c
		call eff_out(ict, u_eff1,u_eff2,u_eff3,epA,epB,epC,epA2,epB2)	!120126homma	
c					!�����|�e���V�����̏o��
c
	endif
c
	nstat	= 1
	nend	= jpnum
c
	if(ict.eq.-1)then
		open(unit=560,file='ryuushi.dat')
	endif
c
c	----(�Փ˓d��)----
 1000	continue
c
	do n=nstat,nend		!1~nend�܂ł̗��q�J��Ԃ�
c
c-----�k�ތ���-----
c	flag_d=0
c	flag_s=0
	iak=0
c	nak=0
	pflag = 0
c 2001 continue			!121009
	count_flag = 0		!10/05/07 �U���J�E���g	!121009
		
c		iiflag = 0		!���W�F�N�V�����E�Փ˓d�����p�t���O11/08	 !�R�����g�A�E�g121009

	jak=0.0
c------------------
c
c		n�ڂ̗��q�́E�E�E
		kv	= kp(1,n)	!�JNo.
		ken	= kp(2,n)	!�G�l���M�[�e�[�u��No.
		kl  = kp(3,n)	!�wNo
		ka	= iarea(kl)	!�f��No�D
		kl2 = kl		!kl2:���q���i��(����)����ꏊ�̑w
		ka2	= ka		!ka2:���q���i��(����)����ꏊ�̑f��
c
c		kpn	= kp(n)-1
c		ka	= kpn/(nvalley*nenergy)+1		!�JNo.
c		ken	= mod(kpn,nenergy)/nvalley+1	!�G�l���M�[�e�[�u��No.
c		kv	= mod(kpn,nvalley)+1			!�f��No�D
c
		akx	= p(1,n)
		aky	= p(2,n)
		akz	= p(3,n)
		ts	= dble(p(4,n))
		x	= p(5,n)
		z	= p(6,n)
		t1	= 0.0
c	check_point!!
		x_start_1 = x_start(n)		!120817sato
		x_start_r = x_start_1
		z_start_1 = z_start(n)		!120817sato
		z_start_r = z_start_1
c
	if(ict.eq.-1)then
		write(560,*) x*1e9,z*1e9,n,ict
	endif
c
c		tau = tau
c
c---	�@�d�q������w�̕s�����Z�x�𒲂ׂ� ----	  !07/8/4 �s�����U��
c
		do ipart=1,npart
			if((ipart.eq.npart-2).or.(ipart.eq.npart-1))cycle  !�d�Ɉʒu
			if(p(6,n).le.czpart2(ipart))then
				kpart=ipart
				exit	
			else
				kpart=npart
				if(p(6,n).gt.czpart2(npart))then
					write(*,*)'emcd�ŕs�����Z�x�G���[1'
					stop
				endif
			endif
		enddo
c--------------------------------------------------
c
		iflag = 0
		balis_flag2 = balis_flag(n)		!�o���X�e�B�b�N�̌v�Z 07/11/22 ��[
c
c
	loop : do while(kp(1,n).ne.0)
c
c
c************���q�h���t�g����**************************************************
c#####�d�E���x�̎擾############################################################
c		call field(x,z,pdx,pdz,u,fx,fz,dhet,kl,
c     &				cxr1,czr1,cxr2,czr2,lxr1,lzr1,lxr2,lzr2)
c
c	�����|�e���V��������d�E���x���擾����	  2005/12/02 Hara �ύX2006/12/09
		if (kv.eq.1) then
			call field(x,z,pdx,pdz,u_eff1,fx,fz,dhet,kl,
     &				cxr1,czr1,cxr2,czr2,lxr1,lzr1,lxr2,lzr2)
		elseif (kv.eq.2) then
			call field(x,z,pdx,pdz,u_eff2,fx,fz,dhet,kl,
     &				cxr1,czr1,cxr2,czr2,lxr1,lzr1,lxr2,lzr2)
		elseif (kv.eq.3) then
			call field(x,z,pdx,pdz,u_eff3,fx,fz,dhet,kl,
     &				cxr1,czr1,cxr2,czr2,lxr1,lzr1,lxr2,lzr2)
		endif
c###############################################################################
c
c
		call drift(
c    &				        fx, fz, dz, dhet, dt, de, pgm(1,1,ka),
     &				        am_aniso,aff_aniso,fx, fz, dz, dhet, dt, 
     &						de, pgm(1,1,kpart),	!07/8/4 �s�����U��
     &				        af2(kv,ka), af4(kv,ka),
     &				        hm(kv,ka), hhm(kv,ka),
     &				        akx,aky,akz,ts,x,z,t1,tau,
     &			            kv,ken,kl,kl2,iflag,ie,ei,sk,
     &						kpart,ka,n)									!07/8/4 �s�����U��


c-----�k�ތ��� drift����(k-kd)�̊e���b�V���E�J�ł̑��a-----
		ix = min(nx,max(0,nint(x*pdx)))
		iz = min(nz,max(0,nint(z*pdz)))
c
			subdkx = akx-p(1,n)
			subdky = aky-p(2,n)
			subdkz = akz-p(3,n)
c			subdkx = akx
c			subdky = aky
c			subdkz = akz
			ix2 = ix
			iz2 = iz
			kv2 = kv
c----------------------------------------------------------
c
		call surf(akx,akz,x,z,kv,
     &		           jspsum,cxpole1,cxpole2,melpos,xmax,zmax,
     &				   cxr1,czr1,cxr2,czr2) !2011/3/25��				
c
c-------------------�o���X�e�B�b�N�̌v�Z 07/11/22 ��[--------------------
		if((cxpole1(2).lt.x).and.(cxpole2(2).gt.x))then
 			if(balis_flag2.eq.0)then
				balis_flag2 = 1
				balis_all=balis_all + 1
			endif
c		elseif((cxpole1(2).lt.x).and.(cxpole2(2).gt.x))then
		else
			balis_flag2 = 0 
		endif 		
c----------------------------------------------------
c
c************�C�x���g����******************************************************
c	----���q�����ł����ꍇ----
		if(kv.eq.0)then
	exit loop
		endif
c
c-----������Ԃ̏����i�[-----	!121009
c	p(1-6,n) ... ���q���(1-3:k���W(kx,ky,kz)[m^-1],4:�U������[s],5-6:�ʒu(x,y)[m])
c	kp(1-6,n) ... �����J
		sp = 0	;	skp = 0			
		sp(1)  = akx
		sp(2)  = aky
		sp(3)  = akz
		sp(4)  = ts	
		sp(5)  = x
		sp(6)  = z
		skp(1) = kv		!�JNo.
		skp(2) = ken		!�G�l���M�[�e�[�u��No.
		skp(3) = kl		!�wNo
c------------------------------
 2001 continue

		iiflag = 0		!���W�F�N�V�����E�Փ˓d�����p�t���O11/08
		
c-----�I��Ԃ�������Ԃ̏��ɖ߂�----- !121009
c	p(1-6,n) ... ���q���(1-3:k���W(kx,ky,kz)[m^-1],4:�U������[s],5-6:�ʒu(x,y)[m])
c	kp(1-6,n) ... �����J
	if(rflag.eq.1)then	!reje1
		akx	= sp(1)
		aky	= sp(2)
		akz	= sp(3)
		ts	= sp(4)
		x	= sp(5)
		z	= sp(6)
		kv	= skp(1)
		ken	= skp(2)
		kl	= skp(3)
		iscat = 0			
c		siscat = 0
	endif
c-----------------------------
c	----�т�dt�ɒB�����ꍇ(iflag=2)----
c	----�U���C�x���g�����������ꍇ(iflag=0)----
		if((iflag.eq.0).or.(iflag.eq.2)) then
c
c
c			---	���� ---
			if(iflag.eq.2)then
	exit loop
			endif
c
			if(iflag.eq.0)then
c	----�U���C�x���g�����������ꍇ(iflag=0)----
c
				den=ei/de(ken)-float(ie-1)		!de:ie��ei�̍��i�Y���j
				if((den.gt.1.0).or.(den.lt.0.0))then
					form="('den�̒l���s���ł�(scat) den=',f,'ei=',e)"
					write(* ,form)den,ei
					write(99,form)den,ei
					
c	 check_point
				endif

				ix = min(nx,max(0,nint(x*pdx)))
				iz = min(nz,max(0,nint(z*pdz)))
				mtp = mtemp(ix,iz)		!���q�̂��郁�b�V���̉��x
c
c-----yama�ǉ�-----
c	�U���C�x���g�����������ꍇ�A���ł��J�E���g!
c	���� +������kx�������q	 siflag=2�@
c	���� -������kx�������q	 siflag=1
			ix3 = ix
			if(akx.gt.0.0) then
				siflag = 2
			else
				siflag = 1
			endif
c------------------
c	 check_point	error message
c
c	if((bscat_xflag.eq.1).and.(p(5,n).ge.(cxpole2(2)+65e-9)))then	!���q���h���C���̈�ɂ�����
	if((bscat_xflag.eq.1).and.(p(5,n).ge.(cxpole2(2))))then	!���q���h���C���̈�ɂ�����
	if(ei.ge.0.1)then
	goto 20
	endif
	endif
		call scat(
     &				 am_aniso,aff_aniso,hole_am_aniso,
     &				 hole_aff_aniso,smh(1,ka),
     &				 af(1,ka),
     &				 iiflag,iiix,iiiz,	!!I.I.�␳�p�ǉ� 11/08��
c     &			     swk(1,1,1,ken,mtp,ka), escat(1,1,ka),
     &			     swk(1,1,1,ken,mtp,kpart), escat(1,1,ka),	!07/8/4 �s�����U��
     &		         iarg(1,ka), iband(1,1,ka),
     &	             eps(ka), bktq(mtp), dopem(ix,iz),
     &		         akx, aky, akz, kv, kl,
     &                 sk, ei, ef, ie, den,
     &				 hescat(1,1,ka),hef,hcss,jp,iscat,
c     &				 dx,dz,jpnum,kp,p,pgm(1,1,ka),n,t1,
     &				 dx,dz,jpnum,kp,p,pgm(1,1,kpart),n,t1,	!07/8/4 �s�����U��
     &				 nava,cnava,x,z,count_flag,six,siscat,skv,			!10/05/07 �U���J�E���g
     &				 cxpole1,cxpole2,
     &				 balis_scat,balis_flag2,balis_n,		!07/03/15
     &				 ccs,ccs2,cncs,cncs2,allback_scat,back_scat,		!08/1/28
     &				 n_scat_p,n_scat_n,dn3(kpart),						!09/2/19 �|��
     &				lxr1,lxr2,ec,eg,hm(kv,ka),ka,
     &			     hiXL,u_eff1,u_eff2,u_eff3,ix,iz,Eth_table,ken,u,			!09/2/19 �|��
     &				lhet,de,											!120330	����
     &				scatpoint,xx_goal,zz_goal,x_start,z_start,II_S,		!120817sato
     &				dltec)
c
c-----yama�ǉ�-----
c	�U���C�x���g�����������ꍇ�A�w�e�����̂�!!�J�E���g
c	�U���� +������kx�������q	 sfflag=3
c	�U���� -������kx�������q	 sfflag=1
20		if((iscat.ne.0).and.
     &		(nchannel1.lt.kl.and.kl.le.nchannel2))then !channel2011/05/25��
c     &				(nlayer-4.lt.kl).and.(kl.le.nlayer-1)) then  !�`���l����
c     &				(kl.gt.nlayer-3).and.(kl.le.nlayer-2)) then  !�`���l����
				if(akx.gt.0.0) then
					sfflag = 3 
				else
					sfflag = 1
				endif
c
				smpf = siflag + sfflag -1	!smpf = 1����4
				sscnt(ix3,iscat,smpf)=sscnt(ix3,iscat,smpf)+1	!�eix�̎U���������ssflag
				siflag = 0 ; sfflag = 0
				smpf = 0
c
c-----kx�� ����(+) ��(+) >> smpf = 4
c-----kx�� ����(-) ��(+) >> smpf = 3
c-----kx�� ����(+) ��(-) >> smpf = 2
c-----kx�� ����(-) ��(+) >> smpf = 1
c
		endif
c------------------
c
				if(	(hef.ne.0.0) .and.
     &				(hcss.eq.1 ) .and.
     &				(kv .ne.0  )) then
c
					call charge_heat(
     &						 x,z,pdx,pdz,hef,hef_mesh,hef_scat,
     &						 jp,iscat)
				endif


c				ts = dble(t1)-dble(log(rnd())*pgm(kv,ken,ka))
				ts = dble(t1)-dble(log(rnd())*pgm(kv,ken,kpart))	!07/8/4 �s�����U��			
c=============�G���[����========================================
				if(sngl(ts).lt.t1)then
					write( *,*)'�ј_���G���[2(emcd2-�U����)'
					write(99,*)'�ј_���G���[2(emcd2-�U����)'
					kv=0
	exit loop
				endif
c===============================================================
c
			endif
c
c	----�w�e����ǏՓ˃C�x���g������(iflag=1)----
		elseif(iflag.eq.1)then
c		----(���Z�X)----
		  do ii = 1,nrecess	
			if((x.gt.cxr1(ii)).and.(x.lt.cxr2(ii)).and.(z.eq.czr2(ii)))then
			  akz=-akz
c
				if ((ii.eq.2)
     &			  .and.(czr2(1).eq.czr2(2))) then !�ً}�����كP�[�X11/04/08��
				  akz=-akz
				endif
c
			  iflag = 0
			endif
		  enddo     

		  if(iflag.eq.1)then
				call border(hhm,af,af2,af4,eg,
     &				akx,aky,akz,kv,kl,kl2,ka,iarea,
     &				lhet,min(nz,max(0,nint(z*pdz))))

c
c------���t�l�X�U��--R.P. Joshi���f��(2012�N�H)-------------------------
c				call border_roughness1(hhm,af,af2,af4,eg,
c     &				akx,aky,akz,kv,kl,kl2,ka,iarea,
c     &				lhet,min(nz,max(0,nint(z*pdz))),
c     &				swk_rou,ken,min(nx,max(0,nint(x*pdx))),de,		!121029sato
c     &				roughness1_countx,roughness1_counte)
c---------------------------------------------------------------------------

c------���t�l�X�U��--J.R. Watling���f��(2012�N�t����)-------------------------
c				call border_roughness2(hhm,af,af2,af4,eg,
c    &				akx,aky,akz,kv,kl,kl2,ka,iarea,
c     &				lhet,min(nz,max(0,nint(z*pdz))),
c     &				min(nx,max(0,nint(x*pdx))),
c     &            xxx,vvv,basho_reflection,basho_roughness,
c     &            split,delta,lambda,count_reflection,
c     &            count_roughness,average,
c     &			epA,u)			!���t�l�X�U���pep���
c-----------------------------------------------------------------------------

c										     ! �h���t�g�O�ʒu2006/12/09Hara
		  endif
c			ts = dble(t1)-dble(log(rnd())*pgm(kv,ken,ka))
			ts = dble(t1)-dble(log(rnd())*pgm(kv,ken,kpart)) !07/8/4 �s�����U��
c=============�G���[����========================================
			if(sngl(ts).lt.t1)then
				write( *,*)'�ј_���G���[3(emcd2-��ǌ�)'
				write(99,*)'�ј_���G���[3(emcd2-��ǌ�)'
				kv=0
	exit loop
			endif
c===============================================================
c
c	----������ł��Ȃ��ꍇ�i�G���[�j----
		else
			write( *,*)'iflag�G���[(emcd2)',iflag
			write(99,*)'iflag�G���[(emcd2)',iflag
			kv=0
	exit loop
		endif
c#########################################################################
c	�h���t�g�^���O��̈ʒux,xx����̈�AB��ʉ߂��闱�q���J�E���g  120817sato

	call path(dx,z,dhet,cxpole1,cxpole2,
     &			pass,pass_r,x_start_1,x_start_r,xx_goal,scatpoint,			!120817sato
     &			x_mean_free_path_sum,x_mean_free_path_count,
     &			z_start_1,z_start_r,zz_goal,mean_free_path_sum,
     &			mean_free_path_sum2,xcenter)

c
c	enddo loop		 !121009
c


c-----Rejection-Technique-----
	ix = min(nx,max(0,nint(x*pdx)))
	iz = min(nz,max(0,nint(z*pdz)))	
	
c	write(*,*) nx,x,pdx
c	write(*,*) nz,z,pdz
c	write(*,*) ix,iz,nx-1

	jak=akx+aky+akz

c	����Tel��Ef���ςȒl��������Rejection�ɓ���Ȃ�
	if(ict.gt.strej)then		!ict > strej=-19998
		if((rejcnt.gt.0).and.(kv.eq.1))then			!�|�� �ύX rejcnt=20,���J

			if((ix.gt.1).and.(ix.lt.(nx-1)).and.	! 1< ix <309
c
cccc!Rejection�`���l���w�K�p!!����!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     &	(lhet(nchannel1)+1.le.iz).and.
     &	(iz.le.lhet(nchannel2))) then !channel Rejection 11/05/25��
c     &	(iz.ge.lhet(3)+1).and.(iz.le.lhet(4)-1))then !�`���l��4�w�ڌŒ�
c     &			(iz.ge.38).and.(iz.le.46))then		!�|�� �`���l���w�̂�Rejection
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
	if((avsumtel(ix,iz,kv).ge.250).and.				!�K�p�d�q���x�͈�
     &(avsumtel(ix,iz,kv).le.3000))then

c	-----08/5/15 �|��-----1�x�͂����Ă����܂��B
c			if((iak.eq.jak).and.(pflag.lt.rejcnt))then		!�O��Ɠ����G�l���M�[�Ȃ�2001��
c				pflag = pflag + 1
c				go to 2001
c			else
c	----------------------	

			call rejection(
     &						hhm,hm,af2,af4,ec,p,kp,iarea,bkq,
     &						efermi,ix,iz,
     &						akx,aky,akz,kv,kl,rflag,
     &						adkx,adky,adkz,iak,avsumtel,epA,u,eg)	!�|�ݕύX

				if((rflag.eq.1).and.(pflag.le.rejcnt)) then
					pflag = pflag + 1
c
cccccc���W�F�N�V�����E�Փ˓d�� �d����� ���q�E�񐔕␳11/08��	
					if(iiflag.eq.1) then		!!!
					  jpnum = jpnum -1		!!!
					  nava(iiix,iiiz,kv) = nava(iiix,iiiz,kv)-1
				  	  nava(iiix,iiiz,0) = nava(iiix,iiiz,0)-1						
					  cnava(iiix,kv) = cnava(iiix,kv)-1
					  cnava(iiix,0) = cnava(iiix,0)-1
cccccc
cccccc���W�F�N�V�����E�Փ˓d���@�d����� �m�F�p 11/08��
c					  write(*,*) 'emcd',n,ix,iz,jpnum,pflag,ict
c					  write(*,*) 'emcd',nava(iiix,iiiz,0),cnava(iiix,0)
c					  pause
cccccc
					endif						!!!  
cccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc120817sato

					if(scatpoint.ge.1)then		!�U�����N�������i��U��(scatpoint=0)�ɂȂ�Ȃ��������j
					if((dhet(nchannel1).lt.z).and.(z.le.dhet(nchannel2)))then	!channel���Ɍ���

					pass = pass - pass_r		!120817sato rejection�񐔕����炷
					pass_r = 0.0				!120817sato	������
					x_start_1 = x_start_r		!�O�̈ʒu�ɖ߂�
					z_start_1 = z_start_r		!�O�̈ʒu�ɖ߂�
					x_mean_free_path_sum = x_mean_free_path_sum 
     &										-abs(xx_goal-x_start_1)
					mean_free_path_sum = mean_free_path_sum
     &					-sqrt((xx_goal-x_start_1)**2+(zz_goal-z_start_1)**2)

					x_mean_free_path_count = x_mean_free_path_count -1			!�J�E���^
					
			mean_free_path_sum2(xcenter,1)
     &			= mean_free_path_sum2(xcenter,1) -1.0							!�J�E���^
			mean_free_path_sum2(xcenter,2)
     &			= mean_free_path_sum2(xcenter,2)-abs(xx_goal-x_start_1)
			mean_free_path_sum2(xcenter,3)
     &			= mean_free_path_sum2(xcenter,3)
     &			-sqrt((xx_goal-x_start_1)**2+(zz_goal-z_start_1)**2)

					endif	 !channel���Ɍ���
					endif	!scatpoint
	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

					go to 2001
				endif

	

			endif
		endif
	endif
	endif

	enddo loop		!121009

		drikx(ix,iz,kv)=drikx(ix,iz,kv) + subdkx
		driky(ix,iz,kv)=driky(ix,iz,kv) + subdky
		drikz(ix,iz,kv)=drikz(ix,iz,kv) + subdkz
		dri(ix2,iz2,kv2) = dri(ix2,iz2,kv2)+1.0
c-----------------------------
c
c	---(�`���l�����U���p�x)---	!10/05/07
	if(count_flag.eq.1) then	
		if((skv.ne.0).and.(siscat.ne.0)) then	!121009sato	
		n_scat(six,siscat,skv) = n_scat(six,siscat,skv) + 1
		count_scat(six) = count_scat(six) + 1
		endif
	endif
c	--------------------------
		p(1,n) = akx
		p(2,n) = aky
		p(3,n) = akz
		p(4,n) = sngl(ts)
		p(5,n) = x
		p(6,n) = z
		kp(1,n)	= kv	!�JNo.
		kp(2,n)	= ken	!�G�l���M�[�e�[�u��No.
		kp(3,n)	= kl	!�wNo�D
c		kp(n)	= ka *(nvalley*nenergy)		!ka:�f��No�D
c     &			+ ken* nvalley				!ken:�G�l���M�[�e�[�u��No.
c     &			+ kv						!kv:�JNo.

		balis_flag(n) = balis_flag2		!�o���X�e�B�b�N�̌v�Z 07/11/22 ��[

		scatpoint = 0	!���̗��q�ɂȂ�O�ɏ�����120817sato
		pass_r = 0.0	!120817sato 
	
		x_start(n) = x_start_1
		z_start(n) = z_start_1

	enddo	!1~nend�܂ł̗��q�J��Ԃ�
c
c	----(�Փ˓d��)----
	nstat = n
	nend = jpnum
	if(n.lt.(jpnum+1))goto 1000
c
	p(4,1:jpnum)=p(4,1:jpnum)-dt
c
c-----�k�ތ��� �h���t�g�g�������̊e���b�V���E�J�ɂ����闱�q����-----
c	���݂̓��J�����l�� iv=1�̂�
c	do iv=1,nvalley
		do iz=0,nz
			do ix=0,nx
				if(dri(ix,iz,1).eq.0.0)cycle
					adkx(ix,iz,1) = drikx(ix,iz,1) / dri(ix,iz,1)
					adky(ix,iz,1) = driky(ix,iz,1) / dri(ix,iz,1)
					adkz(ix,iz,1) = drikz(ix,iz,1) / dri(ix,iz,1)
			enddo
		enddo
c	enddo
c-------------------------------------------------------------------
	if(ict.eq.-1)then
		close(560)
	endif
c

	return
	end