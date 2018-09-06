	subroutine sim_dt2(
     &              am_aniso,aff_aniso,hole_am_aniso,hole_aff_aniso,
     &			  dt,spnum,istp,ict,jpnum,
     &              dx,dz,xmax,zmax,lhet,iarea,dopem,twodeg,c_ratio,
     &              cxpart1,cxpart2,lxpart1,lxpart2,
     &              vb,vi,cur,cxpole1,cxpole2,lnpole1,lnpole2,
     &              melpos,jspsum,ncon,
     &              smh,hhm,hm,af,af2,af4,eps,eg,ec,bktq,
     &              de,swk,pgm,escat,iarg,iband,
     &              p,kp,
     &              u,cn,cp,cloud,hef_mesh,
     &              maceps,hescat,mtemp,
     &			  hef_scat, hcss, count, jpot, sstat,
     & 			  cxr1,czr1,cxr2,czr2,lxr1,lzr1,lxr2,lzr2,
c     &			  nava,cnava,n_scat,count_scat)
c     &			  nava,cnava,n_scat,count_scat,czpart2)
     &			  nava,cnava,n_scat,count_scat,czpart2,balis_flag,		!07/8/4 �s�����U��
     &			  n_scat_p,n_scat_n,	!08/8/6 �|��
     &			  efermi,avsumtel,avsumconc,i2max,E2,n2,i3max,E3,n3,	!08/8/6 �|��
     &        	  avsumconc1,avsumtei11,ntab1,etab1,dn3,hiXL,
     &              epA,epB,epC,epA2,epB2,ecr,							!120126homma
     &              avsumconc_1,avsumconcA,Eth_table,				!09/2/19 �|��
     &			  am,aff,II_S,
     &			  basho_reflection,basho_roughness,
     &		      pass,pass_r,x_start,xx_goal,scatpoint,
     &			  x_mean_free_path_sum,x_mean_free_path_count,			!120817sato	
     &			  z_start,zz_goal,mean_free_path_sum,
     &			  mean_free_path_sum2,swk_rou,
     &			  roughness1_countx,roughness1_counte,
c############circuit 100110(IKEDA)##############################################
     &			  IDS1_stack,IG1_stack,ISS1_stack,i_or_c,jc_on)	

c%%%%%%%%%% ���� 03/04/11 �ύX %%%%%%%%%%
c
	USE DFPORT		!cput�擾�p
	implicit none
	include 'arraysize.fi'
c############circuit 070720##############################################
	include 'Circuit.fi'
c


c---�ϐ��z��p�p�����[�^---
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c
c---��{�p�����[�^---
	real	dt,spnum
	integer	istp,ict
	integer	jpnum
	integer jpot
	common /step/jinit,jstat,jstep,jdisp
	integer	jinit,jstat,jstep,jdisp
c---�f�o�C�X�\��---
	real	dx,dz,xmax,zmax
	integer(2)	lhet(nlayer)
	integer(1),dimension (nlayer)	::iarea
	real,	dimension (0:nx,0:nz)	:: dopem
	real,	dimension (npart)	:: cxpart1, cxpart2
	integer(2),dimension (npart)	::lxpart1,lxpart2
	real czpart2(npart)	!07/8/4 �s�����U��

c---!120201-----------------------
	real, dimension (nvalley,narea)::am	
	real(8) aff(nvalley,narea)					!120201
c---�d��---
	real,	dimension (npole)	:: vb,vi,cur,cxpole1,cxpole2
	integer(2),dimension (npole)	:: lnpole1,lnpole2
	integer(1),dimension (npole)	:: melpos
	integer(4)	jspsum(npole)
	integer	ncon
c---�̈�ʃp�����[�^---
	real,	dimension (nvalley,narea)	:: smh,hhm,hm,af,af2,af4
	real	eps(narea),eg(nvalley,narea),bktq(ntenum)
	real	ec(nvalley,narea)
	real	twodeg(nlayer)				!2DEG�V�[�g�d�ח�(*1.0e-4 cm^-3)
	real	c_ratio
c---�U���p�����[�^---
	real,	dimension (nenergy)		:: de
c	real,	dimension (nscat,nemax,nvalley,nenergy,ntenum,narea):: swk
c	real,	dimension (nvalley,nenergy,narea)	:: pgm
	real,	dimension (nscat,nemax,nvalley,nenergy,ntenum,npart)::	swk	!07/8/4 �s�����U��
	real,	dimension (nvalley,nenergy,npart)	:: pgm	!07/8/4 �s�����U��
	real,	dimension (nscat,nvalley,narea)	:: escat
	integer(1),dimension (nscat,narea)	:: iarg
	integer(1),dimension (nscat,nvalley,narea)	:: iband
	real,	dimension (npart)	::	dn3									!09/2/19 �|��
c---���q���---
	real	p(6,npmax)
	integer(1),dimension (3,npmax)	:: kp
c---�f�o�C�X�����---
	real,	dimension (0:nx,0:nz)	:: u,cn,cp,hef_mesh
	real,	dimension (0:nx,0:nz)	:: epA,epB,epC,epA2,epB2			!120126homma
c---poisso�p�z��---
	real	maceps
c---���o�͗p---
	integer sstat
c---���Z�X----------
	real,dimension (nrecess)	:: cxr1,czr1,cxr2,czr2
	integer(2),	dimension (nrecess)	:: lxr1,lzr1,lxr2,lzr2
c---�Փ˓d��----
	integer,dimension (0:nx,0:nz,0:nvalley) :: nava
	integer,dimension (0:nx,0:nvalley) :: cnava
	integer,dimension (0:nx,nscat,nvalley) :: n_scat
	integer,dimension (0:nx) :: count_scat
c########circuit 140720###################################
	real IG1_stack(-1:jtp00),IDS1_stack(-1:jtp00),ISS1_stack(-1:jtp00)
	integer i_or_c,jc_on
c
	integer,dimension (0:nx,nscat,nvalley) :: n_scat_p !(+)yama071223
	integer,dimension (0:nx,nscat,nvalley) :: n_scat_n !(-)yama071223
	real, dimension	(nvalley,narea,nvalley)::am_aniso							 !20100624
	real, dimension	(nvalley,narea,nvalley)::aff_aniso							!20100624
	real, dimension	(nvalley,narea,nvalley)::hole_am_aniso						!20100624
	real, dimension	(nvalley,narea,nvalley)::hole_aff_aniso							!20100624
c	real(4) highestX,highestL
	real, dimension	(narea,nvalley)::hiXL
	real, dimension	(4,20000,2) :: Eth_table
	real II_S(narea)		!120921sato

c-----�k�ތ���-----
c	real(8) bkq
	real adkz(0:nx,0:nz,0:nvalley)
	real adkx(0:nx,0:nz,0:nvalley)
	real adky(0:nx,0:nz,0:nvalley)
c
	real mconc(0:nx,0:nz,nvalley)
	real tel(0:nx,0:nz,nvalley),efermi(0:nx,0:nz,nvalley)
c
c-----08/8/6 �|��-----
	real avsumtel(0:nx,0:nz,nvalley)
	real avsumconc(0:nx,0:nz,nvalley)
	real avsumconc_1(0:nx,0:nz,nvalley)  !101221�d�q�Z�x�`���l������
      real avsumconcA(0:nx)                !101221�d�q�Z�x�`���l������2
	real(8) avsumtei1(0:nx,0:nz,nvalley)
	real(8) avsumtei1_1(0:nx,0:nz,nvalley) !101221�d�q�G�l���M�[�`���l������
	real(8) avsumteiA(0:nx)                !101221�d�q�G�l���M�[�`���l������2
	integer	sflag,i2max,i3max
	real,	dimension(0:nx,0:nz) ::	avsumtel1
	real,	dimension(0:nx,0:nz) ::	avsumconc1
	real(8),	dimension(0:nx,0:nz) ::	avsumtei11
c
	real E2(0:i2max),n2(0:i2max)
	real E3(0:i3max),n3(0:i3max)
c
	real ntab1(0:300000,5)		!�Z�x(�ԍ��A�ޗ�)
	real etab1(0:300000,5)		!�G�l���M�[(�ԍ��A�ޗ�)
c
	real,	dimension (0:nx,0:nz)	:: tel1,efermi1
	real,	dimension (0:nx,0:nz)	:: tel2,efermi2
	real,	dimension (0:nx,0:nz)	:: tel3,efermi3
	integer(1) ka3(0:nx,0:nz,nvalley)
	integer,dimension (0:nx,nscat,4) :: sscnt
	integer bscat_xflag,fix_u
c------------------
c
c%%%%%%%%%% ���� 03/04/11 �ǉ� %%%%%%%%%%
c---���M���p�z��---
	integer(2),dimension((nx+1)*(nz+1))::mtemp
	real,	dimension(nscat,nvalley,narea) :: hescat
	real,	dimension(0:nx,nvalley,nscat)::hef_scat
	integer hcss	!�M�v�Z�t���O�A0:�v�Z���Ȃ��A1:�v�Z����
	integer	count	!�M�v�Z�J�E���g
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c	--- �Z���d�׉_ ---
	real	cloud(ngx*2+1,ngz*2+1,ngn)
c----------�o���X�e�B�b�N�̌v�Z---------------
	integer(1),dimension (0:npmax) ::	balis_flag			!07/11/20
	integer		balis_scat,balis_n,balis_all
c-----�U���p�̏W�v-------------------!08/1/21
	real(8),dimension (0:nx) ::	ccs,cncs
	real(8),dimension (0:nx,nscat) ::	ccs2,cncs2
c-----����U���̏W�v-------------------!08/1/28
	real(8),dimension (0:nx,0:nscat) ::	allback_scat
	real(8),dimension (0:nx,0:nscat) ::	back_scat
	real	ecr(7,int(nemax/4))		!120126homma
c
c---�����ϐ�---
	integer i	!,l,n
	real	t,ps
	real	ddt
	real	ts
c	real(8)	sk,ei
	integer icn
	integer	hh, mm, ss, s
	integer	hhr,mmr,ssr,sr,crate,cmax
	integer ibord,jpot2
c	integer	nstep,ipp
	character(255) form
c
c	integer n1
c	integer,dimension(:),allocatable::ix,iz
	REAL(4) ta(2)
	real(8),save,allocatable	:: cput(:)
c	REAL(8)	RTC
c

c------���t�l�X�U��--J.R. Watling���f��(2012�N�t����)-------------------------
      real delta,lambda,average
	integer split,count_reflection,count_roughness,
     &        vvv,xxx 
	integer,dimension (0:nx,2) :: basho_roughness,basho_reflection

c-------------------------------------------------------------------------
c-----	�̈�ʉ߂��闱�q���J�E���g120817sato
	real(8),dimension (10) :: pass,pass_r		!�J�E���^�Cr�̓��W�F�N�V�����΍�
	real	x_start(npmax)
	real	xx_goal
	real	z_start(npmax)
	real	zz_goal
	integer	scatpoint			!�U���D��U������
	real	x_mean_free_path_sum		!���ώ��R�s��
	real	mean_free_path_sum		!���ώ��R�s��
	integer x_mean_free_path_count
	real,dimension (0:nx,3) :: 	mean_free_path_sum2		!���ώ��R�s��x���W	 1:�J�E���^,2:x�������ώ��R�s��,3:���ώ��R�s��

c------���t�l�X�U��--R.P. Joshi���f��(2012�N�H)-------------------------
	real	swk_rou(nemax,nvalley,nenergy,narea)
	integer,dimension (0:nx,narea,2) :: roughness1_countx
	integer,dimension (nemax,narea,2) :: roughness1_counte

c
	real, dimension	(nvalley,narea):: dltec
	call cpu_time(ts)		!CPU���Ԏ擾
	s =	nint(ts)			!�����^�̊J�n��������̌o�ߎ��ԁi�b�j �ɕύX
	CALL SYSTEM_CLOCK(COUNT=sr,COUNT_RATE=crate,COUNT_MAX=cmax)
	if(sr.gt.sstat)then
		sr = (sr-sstat)/crate
	else
		sr = (sr-sstat+cmax)/crate
	endif
c
c
c--	�����ԃf�B�X�v���C�o��	--
	if((jdisp.eq.0).or.(modulo(ict,jdisp).eq.0))then
		t=dt*float(ict)						!�V�~�����[�V�������Ԑݒ�
		ps = t*1.0e12
c
		hh = mod(s/3600,1000)	!hour   = sec/(1h/sec=3600) (0 �` 999)
		mm = mod(s/  60,  60)	!min	(0 �` 59)
		ss = mod(s     ,  60)	!second	(0 �` 59)
c
		hhr = mod(sr/3600,1000)	!hour   = sec/(1h/sec=3600) (0 �` 999)
		mmr = mod(sr/  60,  60)	!min	(0 �` 59)
		ssr = mod(sr     ,  60)	!second	(0 �` 59)
c
		form = "(F8.3,'ps  ict=',I8,4X,'jpnum=',I8,2X,
     &			I3,2(':',I2),',',I3,2(':',I2),F7.2)"
		write(*,form)ps,ict,jpnum,hh,mm,ss,hhr,mmr,ssr,c_ratio
		write(3,form)ps,ict,jpnum,hh,mm,ss,hhr,mmr,ssr,c_ratio
		form = '(2(I8),X,I3,2('':'',I2),F7.2)'
		write(13,form)ict,jpnum,hh,mm,ss,c_ratio	!jpnum.txt
		rewind	40;	write(40)	npmax,xmax,zmax,spnum
		rewind	41;	write(41)	jpnum,p,kp
	endif
c
c ----�V�~�����[�V������----
c--	jpot�ɂ����emcd�𕡐���񂷂悤�ɏo���� --
c--	�ʏ�jpot=1�Ŏg�p���� --
	cn = 0		!cn:�S�̔z��
	icn = 0
	ddt = dt
	dt=dt/float(jpot)
	i=istp		!�G���[���p
c
c---	!�z���`�E�ŏ��������s ----
	if (.not. allocated(cput)) then
		allocate (cput(8))				!08/8/6 �|�ݕύX
		cput = dtime(ta)
		cput = 0.0
		adkx=0.0;adky=0.0;adkz=0.0		!08/8/6 �|��
		sscnt = 0						!08/8/6 �|��
		efermi = -1.0595				!08/8/6 �|��
		tel = 0.0						!08/8/6 �|��
	endif
c---------------------------------
c
	jpot2 = jpot
	do i=1,jpot
c
c----		(���q�̋����̃V�~�����[�g )----
		cput(8) = cput(8)+dtime(ta)		!08/8/6 �|�ݕύX
		if(jdisp.eq.0)write(*,*) 'emcd'
c
		call emcd2(
     &                 am_aniso,aff_aniso,hole_am_aniso,hole_aff_aniso,
     &				 dt,de,jpnum,
     &                 dx,dz,xmax,zmax,lhet,iarea,twodeg,dopem,
     &                 cxpole1,cxpole2,melpos,jspsum,
     &                 smh,hm,hhm,af,af2,af4,eps,eg,bktq,
     &                 swk,pgm,escat,iarg,iband,
     &                 p,kp,u,cn,hef_mesh,
     &                 mtemp,hescat,hef_scat,hcss,
     & 			     cxr1,czr1,cxr2,czr2,lxr1,lzr1,lxr2,lzr2,czpart2,	!07/8/4 �s�����U��
     &				 nava,cnava,n_scat,count_scat,ict,i,
     &				 balis_flag,balis_scat,balis_n,balis_all,	!2006/12/09 Hara
     &				 ccs,ccs2,cncs,cncs2,allback_scat,back_scat,		!08/1/28
     &				 ec,adkx,adky,adkz,efermi,n_scat_p,n_scat_n,
     &				 sscnt,avsumconc,avsumtel,dn3,hiXL,
     &                 epA,epB,epC,epA2,epB2,Eth_table,				!09/2/19 �|�� 120126homma
     &            xxx,vvv,basho_reflection,basho_roughness,
     &            split,delta,lambda,count_reflection,
     &            count_roughness,average,
     &			pass,pass_r,x_start,xx_goal,scatpoint,
     &			x_mean_free_path_sum,x_mean_free_path_count,	!120817sato
     &			z_start,zz_goal,mean_free_path_sum,
     &			  mean_free_path_sum2,bscat_xflag,fix_u,			!15/1/2takahashi
     &				II_S,swk_rou,roughness1_countx,roughness1_counte,			!120921sato	!121029sato	
     &            dltec)
		cput(1) = cput(1)+dtime(ta)
c
c---		(�d�ɕt�߂ł̗��q�̑��� )---
		if(jdisp.eq.0)write(*,*) 'renew'
c
		call renew(
     &                jpnum,ncon,
     &                bktq,mtemp,dx,dz,lnpole1,lnpole2,melpos,jspsum,vb,
     &                pgm,smh,p,kp,lhet,iarea,twodeg,balis_flag,
     &				hhm,af,af4,ecr,de,				!120126homma
     &				i2max,E2,n2,i3max,E3,n3,		!�|��
     &				x_start,z_start)			!120817sato
		cput(2) = cput(2)+dtime(ta)
c
c
c	enddo		!C
c	jpot2 = 1	!C
c
c---	(���b�V��������̗��q�Z�x�̃J�E���g)---
		if(jdisp.eq.0)write(*,*) 'charge'
		ibord= lhet(nlayer-1)+2
		call charge(jpnum,dx,dz,spnum,p,kp,cn,icn,jpot2,
     &                ibord,cloud,lxr1,lzr1,lxr2,lzr2)
		cput(3) = cput(3)+dtime(ta)
c
	enddo		!C
c
	dt = ddt
	if(hcss.eq.1 )count = count + 1
c	count = count + 1
c
c---	(���q�Z�x����|�e���V�����̓��o(�|�A�\��������) )---
	if(jdisp.eq.0)write(*,*) 'poisso'
	maceps = epsilon(maceps)**2  !*10.0
	if(modulo(ict,50).eq.(50-1))then		
c		maceps = 2.0E-13 !7		!�v�Z�@�C�v�V����
	else
c		maceps = 2.0E-13 !7		!�v�Z�@�C�v�V����
	endif
	if(fix_u.ne.1)then
	call poisso(dx, dz, lhet, twodeg,iarea,
     &			vb,vi,lnpole1,lnpole2,melpos,eps,eg,c_ratio,
     &			u,cn,cp,dopem,maceps,
     &			lxr1,lzr1,lxr2,lzr2,dltec)
	endif
	cput(4) = cput(4)+dtime(ta)
c
	call eltemp(adkx,adky,adkz,
     &			jpnum,spnum,p,kp,
     &			af2,af4,hhm,ec,
     &			dx,dz,xmax,zmax,iarea,lhet,
     &			tel,mconc,cn,
c     &			ka3)
     &			ka3,tel1,ict,avsumtel,avsumtel1,		!�|�ݕύX
     &			avsumconc,avsumtei1,sflag,				!�|�ݕύX
     &			efermi,efermi1,avsumconc1,avsumtei11,tel2,
     &            epA,epB,epC,u,eg,						!120126homma
     &            avsumtei1_1,avsumconc_1,avsumteiA,avsumconcA)	!�|�ݕύX

	cput(5) = cput(5)+dtime(ta)

c-----�|�� sflag�̔���-----
	if(sflag.eq.1)then	!�X�e�b�v���ς�����Ƃ�����
		call fermi_calc(efermi,ka3,
     &					tel1,tel2,tel3,efermi1,efermi2,efermi3,
     &					avsumconc,avsumtel,avsumtei1,ntab1,etab1,
     &                    avsumteiA,avsumconcA,lhet,am,aff)	!�|�ݕύX
	endif
c--------------------------
	cput(6) = cput(6)+dtime(ta)
c
c---	(�o��)	  ----
	if(jdisp.eq.0)write(*,*) 'output'
	call output(
     &			dt,spnum,de,jpnum,ict,
     &			dx,dz,xmax,zmax,lhet,iarea,
     &			cxpart1,cxpart2,lxpart1,lxpart2,
     &			cxpole1,cxpole2,lnpole1,lnpole2,jspsum,melpos,
     &			hhm,hm,af2,af4,eps,ec,p,kp,u,cn,hef_mesh,
     &			cur,hef_scat,count,hcss,cput,
     &			nava,cnava,n_scat,count_scat,
     &			balis_scat,balis_n,balis_all,
     &			ccs,ccs2,cncs,cncs2,allback_scat,back_scat,		!08/1/28
     &			tel1,tel2,tel3,efermi1,efermi2,efermi3,			!�|�ݒǉ�
     &			n_scat_p,n_scat_n,sscnt,						!�|�ݒǉ�
     &			avsumtel1,avsumconc1,avsumtei11,
     &            epA,epB,epC,epA2,epB2,eg,ecr,			!120126homma
     &            avsumteiA,avsumconcA,				!�|�ݒǉ�
     &            xxx,vvv,basho_reflection,basho_roughness,
     &            split,delta,lambda,count_reflection,
     &            count_roughness,average,
     &			pass,x_mean_free_path_sum,x_mean_free_path_count,			!120817sato	
     &			mean_free_path_sum,mean_free_path_sum2,
     &			roughness1_countx,roughness1_counte,bscat_xflag,fix_u,		!15/1/2takahashi
c############circuit 2014/12/01(takahashi)###################################
     &			IDS1_stack,IG1_stack,ISS1_stack,istp,i_or_c,jc_on)

	cput(7) = cput(7)+dtime(ta)
c
	if(jdisp.eq.0)write(*,*) 'other'
c
c############################
c ���V�ǉ�20030603	
c
c	�폜(ver.4.9.10)
c	
    	return
	end
