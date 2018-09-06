c################################################################
c		�d���h�炬�̃f�o�C�X�����z�Z�o�v���O����
c
c 
c 2014/12/20,�����@�@�i�}�g��̍��삳��̓d���h�炬�Z�o���@�j
c################################################################
c
      subroutine current_fluctuation(jpnum,dx,dz,xmax,zmax,
     &fix_u,ec,cur,dt,cxpole1,cxpole2,
     & p,kp,hhm,hm,af4,af2,iarea,adkx,adky,adkz,spnum,cff,bscat_xflag)																
c===����===
c---�ϐ��z��p�p�����[�^---
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c---��{�p�����[�^---
	real(8) q,h,bk,pi
	real spnum,dt
	parameter(h  = 1.05459e-34, q  = 1.60219e-19,		!h��h(ber)�ł���P�ʂ�[Js]
     &		 bk = 1.38066e-23, pi = 3.1415927)			!bk�̒P�ʂ�[JK-1]
	real cur(3)
	real di_node,dis_node
c
	integer	jpnum
	integer cff,bscat_xflag,fix_u
c---�f�o�C�X�\��---
	real	dx,dz,xmax,zmax
	integer(1)	iarea(nlayer)
c---�d��---
	real	cxpole1(npole),cxpole2(npole)
c---�̈�ʃp�����[�^---
	real,	dimension (nvalley,narea)	:: hm,hhm,af2,af4	!hhm=h(ber)^2/(2mq)
	real,	dimension (nvalley,narea) :: ec					!hm = h(ber)/m*
c---���q���---
	real,	dimension   (6,npmax)	:: p
	integer(1),dimension (3,npmax)	:: kp
c---kd�Ɋւ���p�����[�^---	
	real adkx(0:nx,0:nz,0:nvalley)
	real adky(0:nx,0:nz,0:nvalley)
	real adkz(0:nx,0:nz,0:nvalley)
	real akx(0:nx,0:nvalley)
	real aky(0:nx,0:nvalley)
	real akz(0:nx,0:nvalley)
c
	real tkx,tky,tkz
c--- �V�~�����[�V�������� ---
c	common /outp/jouti,joutp2,joutp,jheat,jout2d,joute,joute2,jcput,sw
c	integer	jouti,joutp2,joutp,jheat,jout2d,joute,joute2,jcput,sw(9)
c	integer	hcss
c	
c--- ���[�J���ϐ� ---						
	integer(8),save :: jt				!do���[�v�J�E���^
	integer,save :: idx, ixmax			!�G����͗p���b�V����
	integer ix,iz,iix	 						!���q�̑��݂��郁�b�V��
	real(8),save :: sv,sv_v			!���b�V�����̗��q��
	allocatable sv(:),sv_v(:,:)			
	real(8),save :: sn,sn_v,ssn,ssn_v					!�J�ʃ��b�V�������q��
	allocatable sn(:),sn_v(:,:),ssn(:),ssn_v(:,:)		
	real(8),save :: atv1,atv1_v	,atvth_v		
	allocatable atv1(:),atv1_v(:,:),atvth_v(:,:)		
	real(8),save :: atn,atn_v			
	allocatable atn(:),atn_v(:,:)		
	real(8),save :: dvi,dvi2,dvi_v,dvi2_v	
	allocatable dvi(:),dvi2(:),dvi_v(:,:),dvi2_v(:,:)
	real(8),save :: dv1,dv1_v,dv2,dv2_v
	allocatable dv1(:),dv1_v(:,:),dv2(:),dv2_v(:,:)
	real(8),save :: sdv2,sdv2_v
	allocatable sdv2(:),sdv2_v(:,:)
	real(8),save :: sdvi,sdvi2,sdvi_v,sdvi2_v,sdv1,sdv1_v	
	allocatable sdvi(:),sdvi2(:),sdvi_v(:,:),sdvi2_v(:,:),
     &			sdv1(:),sdv1_v(:,:)			
	real(8),save :: dn,dn2,dn_v,dn2_v
	allocatable dn(:),dn2(:),dn_v(:,:),dn2_v(:,:)
	real(8),save :: sdn,sdn2,sdn_v,sdn2_v
	allocatable sdn(:),sdn2(:),sdn_v(:,:),sdn2_v(:,:)
	real(8),save :: avi,avi_v
	allocatable avi(:),avi_v(:,:)
	real(8),save :: di2,di2_v
	allocatable di2(:),di2_v(:,:)
	real(8),save :: di_s,di_t,di_st,di2_stc
	allocatable di_s(:),di_t(:),di_st(:),di2_stc(:)
	real(8),save :: di_sv,di_tv,di_stv,di2_stcv		
	allocatable di_sv(:,:),di_tv(:,:),di_stv(:,:),di2_stcv(:,:)
	real(8),save :: tel_v,stel_v,stel					!�d�q���x
	allocatable tel_v(:,:),stel_v(:,:),stel(:)
	real(8),save :: Sv_1,Si_1,Sv_2,Si_2,dv0,dv00						!Gonzaletz�̌v�Z�ŗp����ϐ�
	allocatable Sv_1(:),Si_1(:),Sv_2(:),Si_2(:),dv0(:),dv00(:)
	real(8),save :: var_x,atvth							!���z�֐��̕��U
	allocatable var_x(:),atvth(:)
c	real(8),save :: rv1,rv2,rv3,rv4,rv5
c	allocatable rv1(:),rv2(:),rv3(:),rv4(:),rv5(:)
c	real(8),save :: rv1_v,rv2_v,rv3_v,rv4_v,rv5_v
c	allocatable rv1_v(:,:),rv2_v(:,:),rv3_v(:,:),rv4_v(:,:),rv5_v(:,:)
	real(8) vth,ei,ei1,sid,sid2
c	
	integer(8),save :: k,m					!���ԕ��ς̃J�E���^
c
	integer kv,ka,kl					
	real sk,sk1,sk2,gk,v1,L1,L2								!v1�F�Q���x�ix�����jL1,2:Lgeff�͈̔͑I��
	real(8),save :: a,b,c,vave_all,n_all,w,freq,var			!vave_all�̓`���l���̈�S�̂̃h���t�g���x
c	
c-------------------���b�V���̐ݒ�------------------------------
	if(.not. allocated(atv1))then
		jt = 0
		idx = 1									!�܂Ƃ߂郁�b�V���̐�  
		ixmax = ifix((xmax-float(idx-1)*dx)/(float(idx)*dx))
		k = 0	
		m = 0	
		a = (q*spnum / xmax)**2					!spnum=75000
c---------------simpson method �̐ϕ��͈͂̐ݒ�-----------------
c		b = 0.0
c		c = ixmax*dx
		di_node = 0.0
		dis_node = 0.0		
c--------------Gonzaletz�̕��@�ɂ��v�Z�ɂ�����ݒ�------------
		freq = 1.0e11 							!100GHz
c---------------------------------------------------------------				
c		allocate(rv1(-1:nx),rv1_v(-1:nx,nvalley))
		allocate(sv(-1:nx),sv_v(-1:nx,nvalley))
		allocate(sn(-1:nx),sn_v(-1:nx,nvalley))
		allocate(atv1(-1:nx),atv1_v(-1:nx,nvalley))
		allocate(atn(-1:nx),atn_v(-1:nx,nvalley))
c----------------------------�[���N���A-------------------------------------
c	    rv1 = 0.0		!�v�Z�덷�ۏ�(subroutine sigma)
c	    rv1_v = 0.0		!�o�O����i��񗎂�or���̃T�u���[�`���Ƃ̑��d��`�j
	    sv = 0.0		!�h���t�g���x�̍��v�l�i�S�X�e�b�v�j
	    sv_v = 0.0
	    sn = 0.0		!���q�����v�l
	    sn_v = 0.0   
	    atv1 = 0.0		!�h���t�g���x�̎��ԕ���
	    atv1_v = 0.0    
		atn = 0.0		!���q���̎��ԕ���
		atn_v = 0.0
		vave_all = 0.0
		n_all = 0.0
		var = 0.0
	endif
c
c######Ramo-Shockley�̃T���v�����O�̈�̐ݒ�(�����Q�[�g������R�w���Œ�`)######
c-----------------------���U�[�o�[�̕��ϗ��q���̒��o----------------------------
	if((jt.ge.7500).and.(jt.le.8000)) then
		do n=1,jpnum
		ix = min(nx,max(0,nint(p(5,n)/(idx*dx))))
		atn(ix) = atn(ix) + 1
		enddo
	k = k +1
	endif
c
	if(jt.eq.8000)then
		do ix=1,ixmax
		atn(ix) = atn(ix) / k
		enddo
		k = 0
c
		do ix=1,ixmax
			if((ix.ge.150).and.(ix.le.170))then
			n_all = n_all + atn(ix)
			k = k +1
			endif		
		enddo
		n_all = n_all /k						!���U�[�o�[�̕��ϗ��q��
		k = 0
c--�����Q�[�g���̐ݒ�i���q�������U�[�o�̕��ϒl���10���ȉ��i��j�ɂȂ����Ƃ���)--
			do ix=170,ixmax
				if(atn(ix).le.(0.80*n_all))then
				L1 = ix
				exit
				endif
			enddo
c
			do ix=245,ixmax
				if(atn(ix).ge.(0.80*n_all))then
				L2 = ix
				exit
				endif
			enddo
c-----------Ramo-Shockley�̃T���v�����O�̈�̐ݒ�(�C�ӂɒ�`)----------------
c				L1 =
c				L2 =
c----------------------------------------------------------------------------
		n_all = 0.0
		atn = 0.0
	endif
c###########�e���b�V�����Ƃ̃h���t�g���x�C���q���̎��ԕ���##################
c---------------���b�V�����Ƃ̃h���t�g���x,���q���̍��v---------------------
	if((jt.ge.8000).and.(jt.le.9000)) then
		if(jt.eq.8000) then
			bscat_xflag = 0					!����U���}���t���O	
			fix_u = 0						!�|�e���V�����Œ�t���O
			if(fix_u.eq.1) write(*,*)'fixed potential'
			if(bscat_xflag.eq.1) write(*,*)'prohibite back scattering'
		endif		
	do n=1,jpnum
	kv = kp(1,n)	!�J���
	kl = kp(3,n)	!�wNo�D
	ka = iarea(kl)	!�wNo.����f��No.��
      ix = min(nx,max(0,nint(p(5,n)/(idx*dx))))					!���q�̃��b�V���ʒu
	sk1 = p(1,n)*p(1,n)+p(2,n)*p(2,n)+p(3,n)*p(3,n)				!�k�ތ��ʕK�v�H
      v1 = p(1,n)*hm(kv,ka)/sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*sk1)	!�Q���x�i�h���t�g�����j
	sv(ix) = sv(ix) + v1										!�e�Z�N�V�������Ƃ̃h���t�g���x�̘a
		sv_v(ix,kv) = sv_v(ix,kv) + v1								!�J��
c     call sigma(v1, sv(ix), rv1(ix))								!�h���t�g�����̘a
c     call sigma(v1, sv_v(ix,kv), rv1_v(ix,kv))					!�h���t�g�����̘a(�J�ʁj
		
      sn(ix) = sn(ix)+1.0											!���b�V�����̗��q��
	sn_v(ix,kv) = sn_v(ix,kv)+1.0									!�J�� 
c---------�T���v�����O�̈�S�̂ł̃h���t�g���x&�L�����A���̘a
	if((ix.ge.L1).and.(ix.le.L2))then					!Ramo-Shockley�̃T���v�����O�̈�
	vave_all = vave_all + v1									
	n_all =n_all + 1.0
	endif															
	end do
	di_node	= di_node + cur(3)									!�[�q�d���̎��ԃX�e�b�v�a
	k = k + 1
	endif
c-------------------���b�V�����Ƃł̕��σh���t�g���x,���q��(���ԕ���)------------------
	if((jt.ge.9001).and.(jt.le.10000)) then	
		if(jt.eq.9001)then
		do ix=1,ixmax
			if(sn(ix).ne.0)then
				atv1(ix) = sv(ix)/sn(ix)		   		!�e���b�V���ł̑��x�̕��ϒl
				atn(ix)   = sn(ix) / k	
			else
				atv1(ix) = 0.0
				atn(ix) = 0.0
			endif
		enddo
		di_node	= di_node / k							!�[�q�d���̎��ԕ��ϒl
c--------------------------------------�J��--------------------------------------------
		do kv=1,nvalley
			do ix=1,ixmax
				if(sn_v(ix,kv).ne.0)then
					atv1_v(ix,kv) = sv_v(ix,kv)/sn_v(ix,kv)	
					atn_v(ix,kv) = sn_v(ix,kv) / k
c
				else
					atv1_v(ix,kv) = 0.0
					atn_v(ix,kv) = 0.0
c
				endif
			enddo
		enddo
c
		vave_all = vave_all / n_all
		n_all = n_all / k								!�����Q�[�g�����̕��ϗ��q��
c
c		allocate(rv1(-1:nx),rv1_v(-1:nx,nvalley))
c		allocate(rv2(-1:nx),rv2_v(-1:nx,nvalley))
c		allocate(rv3(-1:nx),rv3_v(-1:nx,nvalley))
c		allocate(rv4(-1:nx),rv4_v(-1:nx,nvalley))
c		allocate(rv5(-1:nx),rv5_v(-1:nx,nvalley))
c
		allocate(ssn(-1:nx),ssn_v(-1:nx,nvalley))	
		allocate(dvi(-1:nx),dvi2(-1:nx),
     &		dvi_v(-1:nx,nvalley),dvi2_v(-1:nx,nvalley))	
		allocate (dv1(-1:nx),dv1_v(-1:nx,nvalley))
		allocate (dv2(-1:nx),dv2_v(-1:nx,nvalley))
		allocate(sdvi(-1:nx),sdvi2(-1:nx),
     &		sdvi_v(-1:nx,nvalley),sdvi2_v(-1:nx,nvalley))	
		allocate (sdv1(-1:nx),sdv1_v(-1:nx,nvalley))
		allocate (sdv2(-1:nx),sdv2_v(-1:nx,nvalley))
		allocate(dn(-1:nx),dn2(1:nx),avi(-1:nx),avi_v(-1:nx,nvalley),
     &		dn_v(-1:nx,nvalley),dn2_v(-1:nx,nvalley))	
		allocate(sdn(-1:nx),sdn2(-1:nx),
     &		sdn_v(-1:nx,nvalley),sdn2_v(-1:nx,nvalley))
		allocate(di2(-1:nx),di2_v(-1:nx,nvalley))
		allocate(di_s(-1:nx),di_t(-1:nx),
     &				di_st(-1:nx),di2_stc(-1:nx))
		allocate(di_sv(-1:nx,nvalley),di_tv(-1:nx,nvalley),
     &									di_stv(-1:nx,nvalley))
		allocate(di2_stcv(-1:nx,nvalley))
		allocate(tel_v(-1:nx,nvalley),stel_v(-1:nx,nvalley))
		allocate(atvth_v(-1:nx,nvalley),stel(-1:nx))
		allocate(Sv_1(-1:nx),Si_1(-1:nx),Sv_2(-1:nx),atvth(-1:nx),
     &			Si_2(-1:nx),dv0(-1:nx),dv00(-1:nx),var_x(-1:nx))

c			
c		rv1 = 0.0 ;	rv1_v = 0.0	;rv2 = 0.0	;rv2_v = 0.0	
c		rv3 = 0.0 ;	rv3_v = 0.0	;rv4 = 0.0	;rv4_v = 0.0
c		rv5 = 0.0 ;	rv5_v = 0.0		
	    ssn=0.0 ; ssn_v=0.0 ; avi_v=0.0 ; sn=0.0 ; sn_v=0.0 ;sv_v=0.0
		sdvi=0.0 ; sdvi2=0.0 ; sdvi_v=0.0 ; sdvi2_v=0.0 ; avi=0.0
		sdv1=0.0 ; sdv1_v=0.0 ; sdv2=0.0 ; sdv2_v=0.0 ; stel_v=0.0
		stel_v=0.0 ; sdn=0.0 ; sdn2=0.0 ; sdn_v=0.0 ; sdn2_v=0.0
		Sv_1=0.0 ; Si_1=0.0 ; Sv_2=0.0 ; Si_2=0.0 ; dv0=0.0 ; dv00=0.0
		var_x=0.0 ; atvth=0.0
	   endif
c##########################�d���h�炬��I�O�Q###############################
c		�M���x�Ǝ��ԕ��σh���t�g���x�̍�����d���h�炬���Z�o
c		(�h���t�g���x�̗h�炬���Z�o�j
c##########################################################################
c
c---------------------------�M���xvth�̌v�Z--------------------------------
	do n=1,jpnum
		kv = kp(1,n)
		kl = kp(3,n)
		ka = iarea(kl)	!�wNo.����f��No.��
	    ix = min(nx,max(0,nint(p(5,n)/(idx*dx))))
		iz = min(nz,max(0,nint(p(6,n)/dz)))				!���b�V����eltemp�Ɠ���
		iix = min(nx,max(0,nint(p(5,n)/dx)))
c
		akx = 0.0
		aky = 0.0
		akz = 0.0 
c
		tkx = abs(p(1,n)) - abs(adkx(ix,iz,kv))
		tky = abs(p(2,n)) - abs(adky(ix,iz,kv))
		tkz = abs(p(3,n)) - abs(adkz(ix,iz,kv))	
c		
c		sk = tkx*tkx
		sk = p(1,n)*p(1,n)
		sk1 = p(1,n)*p(1,n)+p(2,n)*p(2,n)+p(3,n)*p(3,n)				!�k�ތ��ʕK�v�H
c		sk2 = tkx*tkx + tky*tky + tkz*tkz
	    v1 = p(1,n)*hm(kv,ka)/sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*sk1)	!�Q���x�i�h���t�g�����j
c
	    sn(ix) = sn(ix) + 1.0									!���b�V�����̗��q��
		sn_v(ix,kv) = sn_v(ix,kv) + 1.0						!���b�V�����̗��q��(�J��) 
		ssn(ix) = ssn(ix) + 1.0								!�S�X�e�b�v�̃��b�V���ʑ����q��
		ssn_v(ix,kv) = ssn_v(ix,kv) + 1.0					
c	   
		if(af4(kv,ka).ne.0.0)then									!����������l��
			sq = sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*sk)
			sq1 = sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*sk1)				!�k�ލl������
			ei1=(sq1-1.0)/af2(kv,ka)+ec(kv,ka)
			ei=(sq-1.0)/af2(kv,ka)+ec(kv,ka)						!�d�q�̃G�l���M�[(Ekx�̂�)
c     &			-epA(ix,iz)+u(ix,iz)-0.7*(eg(1,ka)-eg(1,2))	
		else
			ei1=hhm(kv,ka)*sk1+ec(kv,ka)
			ei=hhm(kv,ka)*sk+ec(kv,ka)						
c     &			-epA(ix,iz)+u(ix,iz)-0.7*(eg(1,ka)-eg(1,2))		
		endif
c		x�����̃G�l���M�[���璊�o�����M���x�i1/2*m*v(x)^2=Ex�j
		vth = sqrt((2.0*q*hm(kv,ka)*ei) / h)				!x���������̔M���x
		atvth(ix) =atvth(ix) + vth
c
c		�g�[�^���̃G�l���M�[����̔M���x
c		vth = sqrt((2.0*q*hm(kv,ka)*ei1) / h)			
c		vth = sqrt(sk/sk1)*vth)				!x���������̔M���x(1��ƃZ�b�g)
c
c		�g�[�^���̓d�q���x
		stel_v(ix,kv) = stel_v(ix,kv) + (2.0/3.0)*q*ei1/bk	!�d�q���x�̘a
		stel(ix) = stel(ix) + (2.0/3.0)*q*ei1/bk			!�d�q���x�̘a
		var_x(ix) = var_x(ix) + sqrt((2.0/3.0)*q*ei1*hm(kv,ka)/h)	!���z�֐��̕��U((KbTe/m*)^1/2)
c		�i���U���M���xvth�Ƃ݂Ȃ���j
c
c		x�����̓d�q���x
c		stel_v(ix,kv) = stel_v(ix,kv) + q*ei/bk			!�d�q���x�̘a
c		stel(ix) = stel(ix) + q*ei/bk					!�d�q���x�̘a

c
c---------------�e���b�V�����Ƒ��x�̗h�炬�̘a(���q�ɑ΂���a------------
c		�M���x�Ƀh���t�g�����������Ă���Ƃ���
c		���σh���t�g���x�̍l������2�ʂ�ł���i���b�V�����j
		if(p(1,n).ge.0)then								!�M���x�𐳕��ŏꍇ����
		dvi(ix)   = vth  - atv1(ix)						!���x�̕΍�
		dvi_v(ix,kv)  = vth  - atv1_v(ix,kv)			!�e�J����
		sv_v(ix,kv) = sv_v(ix,kv) + vth 
		else 
		dvi(ix)   = -vth  - atv1(ix)					!���x�̕΍�
		dvi_v(ix,kv)  = -vth - atv1_v(ix,kv)			!�e�J����
		sv_v(ix,kv) = sv_v(ix,kv) - vth					!���_�I�ɂ�atv1�Ɠ������Ȃ�
		endif
c			�iRamo-Shockley�̒藝�̃T���v���̈�S�́j
c		if(p(1,n).ge.0)then								!�M���x�𐳕��ŏꍇ����
c		dvi(ix)   = vth - vave_all						!���x�̕΍�
c		dvi_v(ix,kv)  = vth - vave_all					!�e�J����
c		sv_v(ix,kv) = sv_v(ix,kv) + vth 
c		else 
c		dvi(ix)   = -vth - vave_all						!���x�̕΍�
c		dvi_v(ix,kv)  = -vth - vave_all					!�e�J����
c		sv_v(ix,kv) = sv_v(ix,kv) - vth					!���_�I�ɂ�atv1�Ɠ������Ȃ�
c		endif
c
		dvi2(ix)  = dvi(ix)**2.0							!���x�̗h�炬
		dvi2_v(ix,kv) = dvi_v(ix,kv)**2.0
c					�@�i���b�V�����j
		dv1(ix) = v1 - atv1(ix)							!�h���t�g���x�̗h�炬
		dv2(ix) = dv1(ix)**2.0
		dv1_v(ix,kv) = v1 - atv1_v(ix,kv)
		dv2_v(ix,kv) = dv1_v(ix,kv)**2.0
c			�iRamo-Shockley�̒藝�̃T���v���̈�S�́j
c		dv1(ix) = v1 - vave_all							!�h���t�g���x�̗h�炬
c		dv2(ix) = dv1(ix)**2.0
c		dv1_v(ix,kv) = v1 - vave_all
c		dv2_v(ix,kv) = dv1_v(ix,kv)**2.0	
c		
c		�M���x�l��
		avi(ix) = avi(ix) + dvi2(ix)				!�d���h�炬�i���I����͎��Ԃɑ΂��Ă̘a�j
		avi_v(ix,kv) = avi_v(ix,kv) + dvi2_v(ix,kv)
c
c		�h���t�g�����̂�
c		avi(ix) = avi(ix) + dv2(ix)					!�d���h�炬�i���I����͎��Ԃɑ΂��Ă̘a�j
c		avi_v(ix,kv) = avi_v(ix,kv) + dv2_v(ix,kv)	
c
		sdvi(ix) = sdvi(ix) +  dvi(ix) 
		sdvi2(ix) = sdvi2(ix) + dvi2(ix) 
		sdvi_v(ix,kv) = sdvi_v(ix,kv) + dvi_v(ix,kv) 
		sdvi2_v(ix,kv) = sdvi2_v(ix,kv)+ dvi2_v(ix,kv)
c
		sdv1(ix) = sdv1(ix) + dv1(ix) 
		sdv1_v(ix,kv) = sdv1_v(ix,kv) + dvi_v(ix,kv) 
		sdv2(ix) = sdv2(ix) + dv2(ix) 
		sdv2_v(ix,kv) = sdv2_v(ix,kv) + dv2_v(ix,kv)
c
c------------------�v�Z�덷�ۏ�p�i�o�O����j--------------------------------
c		call sigma(dvi(ix), sdvi(ix), rv1(ix))					!���x�̗h�炬�̘a
c		call sigma(dvi2(ix), sdvi2(ix), rv2(ix))
c	    call sigma(dvi_v(ix,kv), sdvi_v(ix,kv), rv1_v(ix,kv))
c	    call sigma(dvi2_v(ix,kv), sdvi2_v(ix,kv), rv2_v(ix,kv))	
c
c		call sigma(dv1(ix), sdv1(ix), rv3(ix))					!�h���t�g���x�̗h�炬�̘a
c		call sigma(dv1_v(ix,kv), sdv1_v(ix,kv), rv3_v(ix,kv))	!���x�̗h�炬�̘a
	enddo	
c------------------�e���b�V�����Ƃ̗��q���̗h�炬�̘a�i���Ԃɑ΂��āj-----------
	do ix=1,ixmax
		dn(ix) = sn(ix) - atn(ix) 
		dn2(ix) = dn(ix)**2.0 
		dn_v(ix,1:3) = sn_v(ix,1:3) - atn_v(ix,1:3)				!�e�J����
		dn2_v(ix,1:3) = dn_v(ix,1:3)**2.0
c
		sdn(ix) = sdn(ix) + dn(ix)
		sdn2(ix) = sdn2(ix) + dn2(ix)
		sdn_v(ix,1:3) = sdn_v(ix,1:3) + dn_v(ix,1:3)
		sdn2_v(ix,1:3) = sdn2_v(ix,1:3) + dn2_v(ix,1:3)
	enddo
c
c	Gonzaletz�̕��@�ɂ��Sv�ASin�̌v�Z
	do ix=1,ixmax
		if(m.eq.0)then
		dv0(ix) = sdv1(ix) / sn(ix)	
		endif
		dv1(ix) = sdv1(ix) / sn(ix)								!�Z�N�V�������Ƃ̑��x�΍��i�A���T���u�����ρj
		w =2.0 * pi * freq
		Sv_1(ix) = Sv_1(ix) + 4.0*dv0(ix) * dv1(ix) *cos(w*dt*m) * dt   !�_�Ō��Ă�
		Si_1(ix)=Si_1(ix)+a/(dx*idx)*xmax**2.0*sn(ix)*spnum*Sv_1(ix)	!�_�Ō��Ă�
		Sid = Sid + Si_1(ix) * dx*idx									
c		Cv �� n �̕␳�ɂ��v�Z
		if(m.eq.0)then
		dv00(ix) = sdv1(ix) / sn(ix) * sqrt(sn(ix))
		endif
		dv1(ix) = sdv1(ix) / sn(ix)	* sqrt(sn(ix))				
		Sv_2(ix) = Sv_2(ix) + 4.0*dv00(ix)*dv1(ix)*cos(w*dt*m)*dt 
		Si_2(ix)=Si_2(ix)+a/(dx*idx)*xmax**2*sn(ix)*spnum*Sv_1(ix)
		Sid2 = Sid2 + Si_2(ix) * dx
	enddo
	if(jt.eq.9000)then
		Sid = Sid / m
		Sid2 = Sid2 / m	
	endif
c------------------------------------------------------------------------
c		call sigma(dn(ix), sdn(ix), rv4(ix))					!���x�̗h�炬�̘a
c		call sigma(dn2(ix), sdn2(ix), rv5(ix))
c	    call sigma(dn_v(ix,kv), sdn_v(ix,kv), rv4_v(ix,kv))				
c	    call sigma(dn2_v(ix,kv), sdn2_v(ix,kv), rv5_v(ix,kv))
c
		sn = 0.0		
	    sn_v = 0.0 
c
		dis_node = dis_node + (cur(3) - di_node)**2.0				!�[�q�d���̗h�炬�̐ώZ�l
		m = m + 1
	endif
c----------------------�e���b�V�����Ƒ��x�̗h�炬�̘a�i���Ԃɑ΂��āj-----------
	if(jt.eq.9002)then
	do ix=1,ixmax
		if(ssn(ix).ne.0)then
			dvi(ix) = sdvi(ix) / ssn(ix) 
			dvi2(ix) = sdvi2(ix) / ssn(ix) 
			dv1(ix) = sdv1(ix) / ssn(ix)
			dv2(ix) = sdv2(ix) / ssn(ix)
			stel(ix) = stel(ix) / ssn(ix)						!�e���b�V���ʓd�q���x���z
			var_x(ix) = var_x(ix) / ssn(ix)						!�e���b�V���ʂ̕��z�֐��̕��U
			atvth(ix) = atvth(ix) / ssn(ix)
		else
			dvi(ix)  = 0.0
			dvi2(ix) = 0.0
			dv1(ix)  = 0.0
			dv2(ix)  = 0.0
		endif
	enddo
c
		do ix=1,ixmax
			if(ssn_v(ix,kv).ne.0)then
				dvi_v(ix,1:3)  = sdvi_v(ix,1:3) / ssn_v(ix,1:3)
				dvi2_v(ix,1:3) = sdvi2_v(ix,1:3) / ssn_v(ix,1:3)
				dv1_v(ix,1:3) = sdv1_v(ix,1:3) / ssn_v(ix,1:3)
				dv2_v(ix,1:3) = sdv2_v(ix,1:3) / ssn_v(ix,1:3)
c
				tel_v(ix,1:3) = stel_v(ix,1:3) / ssn_v(ix,1:3)	!�e���b�V���ʓd�q���x���z
				atvth_v(ix,1:3) = sv_v(ix,1:3) / ssn_v(ix,1:3)	!�e���b�V�����Ƃ̃L�����A�̑��x�̕���
c
			else
				dvi_v(ix,1:3)  = 0.0
				dvi2_v(ix,1:3) = 0.0
				dv1_v(ix,1:3) = 0.0
				dv2_v(ix,1:3) = 0.0
c
				tel_v(ix,1:3) = 0.0
				atvth_v(ix,1:3) = 0.0
			endif
		enddo
c---------------�e���b�V�����Ƃ̓d���̗h�炬�i���ԕ��ρj------------------
	do ix=1,ixmax			
c
		avi(ix) = avi(ix)/m 				!�d���h�炬
		avi_v(ix,1:3) = avi_v(ix,1:3)/m 
c
c------------------�e���b�V�����Ƃ̒����q���̗h�炬�i���ԕ��ρj---------------
c
		dn(ix) = sdn(ix)/m 							
		dn2(ix) = sdn2(ix)/m 
		dn_v(ix,1:3) = sdn_v(ix,1:3)/m 			
		dn2_v(ix,1:3) = sdn2_v(ix,1:3)/m 
c--------------------�e���b�V���ɂ�����d���h�炬�̎Z�o---------------------
		di2(ix) = a * avi(ix) / spnum							!(a=q/l)^2
		di2_v(ix,1:3) = a * avi_v(ix,1:3)/ spnum				!�J��(spnum�̉e����l��)
c##################thermal,shot noise�̕���(�h���t�g���x�̗h�炬)##################
	di_s(ix)  = a * dn2(ix) * (atv1(ix)**2.0)	/ spnum						!shot noise
	di_t(ix)  = a * dv2(ix) * (atn(ix)**2.0)	/ spnum						!thermal noise�@?spnum�̈������悭�킩��Ȃ�
	di_st(ix) = abs(a * 2.0*atn(ix)*dn(ix)*atv1(ix)*dv1(ix))/spnum		!cross term
	di2_stc(ix) = di_s(ix) + di_t(ix) +	di_st(ix)						!�h���t�g���x�̗h�炬
c																		(��̒l�ƈ�v���邩�̊m�F)
	di_sv(ix,1:3)  = a * dn2_v(ix,1:3) * (atv1_v(ix,1:3)**2.0)/ spnum		!�J��	
	di_tv(ix,1:3)  = a * dvi2_v(ix,1:3) * (atn_v(ix,1:3)**2.0)/ spnum			
	di_stv(ix,1:3) = a * abs(2.0*atn_v(ix,1:3)*dn_v(ix,1:3)*
     &                   		atv1_v(ix,1:3)*dvi_v(ix,1:3))/ spnum	
	di2_stcv(ix,1:3) = di_sv(ix,1:3)+di_tv(ix,1:3)+di_stv(ix,1:3)
c##################thermal,shot noise�̕���(�M���x�̗h�炬)##################
c	di_s(ix)  = a * dn2(ix) * (atvth(ix)**2.0)						!shot noise
c	di_t(ix)  = a * dv2(ix) * ((atn(ix)*spnum)**2.0)					!thermal noise
c	di_st(ix) = abs(a * 2.0*atn(ix)*spnum*dn(ix)*atvth(ix)*dv1(ix))	!cross term
c	di2_stc(ix) = di_s(ix) + di_t(ix) +	di_st(ix)					!�h���t�g���x�̗h�炬
c																	(��̒l�ƈ�v���邩�̊m�F)
c	di_sv(ix,1:3)  = a * dn2_v(ix,1:3) * (atvth_v(ix,1:3)**2.0)			!�J��	
c	di_tv(ix,1:3)  = a * dvi2_v(ix,1:3) * ((atn_v(ix,1:3)*spnum)**2.0)			
c	di_stv(ix,1:3) = a * abs(2*atn_v(ix,1:3)*spnum*dn_v(ix,1:3)*
c     &                   		atvth_v(ix,1:3)*dvi_v(ix,1:3))	
c	di2_stcv(ix,1:3) = di_sv(ix,1:3)+di_tv(ix,1:3)+di_stv(ix,1:3)
c
c--------------------�e���b�V���̗��q(�d�q)���̗h�炬--------------------------
     	dn(ix) = sdn(ix) * spnum							
	dn2(ix) = sdn2(ix) * spnum**2.0
	dn_v(ix,1:3) = sdn_v(ix,1:3) * spnum			
	dn2_v(ix,1:3) = sdn2_v(ix,1:3) * spnum**2.0
c
	atn(ix) = atn(ix) * spnum
c----------------------------���z�֐��̕��U�̘a---------------------------------
	if((ix.ge.L1).and.(ix.le.L2))then		!Ramo-Shockley�̃T���v�����O�̈�
	var = var + var_x(ix)
	endif
	enddo	
	var = var/(L2-L1)								!���ς̕��z�֐����U�l
c-------------------�[�q�d���̗h�炬�i���ԕ��ρj--------------------------------
		dis_node	= dis_node/m/spnum
c##############################�o�͊֘A#########################################
	do kv=1,nvalley
		do ix=1,ixmax
c
5	format('�[�q�d���̂�炬�i���U�j = ',e15.7)
10	format((I12,E12.3))
100	format((I12,E12.3,E12.3,E12.3,E12.3,E12.3))
150	format((I12,E12.3,E12.3))
170	format((I12,E12.3,E12.3,E12.3))
200	format((I12,E12.3,E12.3,E12.3,E12.3))
300	format(4(A12))
400	format(6(A12))
500	format(5(A12))
1000	format('���z���狁�߂���炬= ',e15.7)
1001	format('<��v> = ',e15.7)
1002	format('<vd> = ',e15.7)
1003	format('<Nch> = ',e15.7)
1004	format('�[�q: <��Id/Id> = ',e15.7)
1005	format('���S�Ɍ��ߎ� : <��Id/Id> = ',e15.7)
1006	format('�����q�� = ',e15.7)
c
	if(kv.eq.1)then
		if(ix.eq.1)then
		write(900,400) 'ix','dvi(x)','dvi2(x)','dv1(x)','dv2(x)','g'
		write(901,500) 'ix','atvth(x)','atv1(x)','vave','g'
		write(902,400) 'ix','dn(x)','dn2(x)','dn2_v(x)','atn(x)','g'			
		write(903,300) 'ix','di2(x)','di2_v(x)','g'
		write(904,400) 'ix','di_s','di_t','di_st','di2_stc','g'
		write(905,400) 'ix','di_sv','di_tv','di_stv','di2_stcv','g'
		write(906,400) 'ix','tel(x)','var(x)','var','g'
		endif
c
	write(900,200) ix,dvi(ix),dvi2(ix),dv1(ix),dv2(ix)
	write(901,170) ix,atvth(ix),atv1(ix),vave_all
c	write(901,200) ix,dv1(ix),dv2(ix),dv2_v(ix,kv),atv1(ix)
	write(902,200) ix,dn(ix),dn2(ix),dn2_v(ix,kv),atn(ix)				
	write(903,150) ix,di2(ix),di2_v(ix,kv)
	write(904,200) ix,di_s(ix),di_t(ix),di_st(ix),di2_stc(ix)
	write(905,200) ix,di_sv(ix,kv),di_tv(ix,kv),  
     &                	di_stv(ix,kv),di2_stcv(ix,kv)
	write(906,170) ix,stel(ix),var_x(ix),var
	endif
c
	if(kv.eq.2)then
		if(ix.eq.1)then
		write(900,400) 'ix','dvi(x)','dvi2(x)','dv1(x)','dv2(x)','l'
		write(901,500) 'ix','atvth(x)','atv1(x)','vave','l'
		write(902,400) 'ix','dn(x)','dn2(x)','dn2_v(x)','atn(x)','l'			
		write(903,300) 'ix','di2(x)','di2_v(x)','l'
		write(904,400) 'ix','di_s','di_t','di_st','di2_stc','l'
		write(905,400) 'ix','di_sv','di_tv','di_stv','di2_stcv','l'
		write(906,400) 'ix','tel(x)','var(x)','var','l'
		endif
c
	write(900,200) ix,dvi(ix),dvi2(ix),dv1(ix),dv2(ix)
	write(901,170) ix,atvth(ix),atv1(ix),vave_all
c	write(901,200) ix,dv1(ix),dv2(ix),dv2_v(ix,kv),atv1(ix)
	write(902,200) ix,dn(ix),dn2(ix),dn2_v(ix,kv),atn(ix)				
	write(903,150) ix,di2(ix),di2_v(ix,kv)
	write(904,200) ix,di_s(ix),di_t(ix),di_st(ix),di2_stc(ix)
	write(905,200) ix,di_sv(ix,kv),di_tv(ix,kv),  
     &                	di_stv(ix,kv),di2_stcv(ix,kv)
	write(906,170) ix,stel(ix),var_x(ix),var
	endif
c
	if(kv.eq.3)then
		if(ix.eq.1)then
		write(900,400) 'ix','dvi(x)','dvi2(x)','dv1(x)','dv2(x)','x'
		write(901,500) 'ix','atvth(x)','atv1(x)','vave','x'
		write(902,400) 'ix','dn(x)','dn2(x)','dn2_v(x)','atn(x)','x'			
		write(903,300) 'ix','di2(x)','di2_v(x)','x'
		write(904,400) 'ix','di_s','di_t','di_st','di2_stc','x'
		write(905,400) 'ix','di_sv','di_tv','di_stv','di2_stcv','x'
		write(906,400) 'ix','tel(x)','var(x)','var','x'
		endif
c
	write(900,200) ix,dvi(ix),dvi2(ix),dv1(ix),dv2(ix)
	write(901,170) ix,atvth(ix),atv1(ix),vave_all
c	write(901,200) ix,dv1(ix),dv2(ix),dv2_v(ix,kv),atv1(ix)
	write(902,200) ix,dn(ix),dn2(ix),dn2_v(ix,kv),atn(ix)				
	write(903,150) ix,di2(ix),di2_v(ix,kv)
	write(904,200) ix,di_s(ix),di_t(ix),di_st(ix),di2_stc(ix)
	write(905,200) ix,di_sv(ix,kv),di_tv(ix,kv),  
     &                	di_stv(ix,kv),di2_stcv(ix,kv)
	write(906,170) ix,stel(ix),var_x(ix),var
c
	cff=1	!�I���t���O
	endif
c
		enddo
	enddo
c
c----------------------Sv,Sin�o��------------------------
	do ix=1,ixmax
	if(ix.eq.1)then
	write(908,300) 'ix','Sv','Sin','Sid'
	write(909,300) 'ix','Sv*','Sin*','Sid*'
	endif
	write(908,170) ix*idx,Sv_1(ix),Si_1(ix),Sid
	write(909,170) ix*idx,Sv_2(ix),Si_2(ix),Sid2
	enddo
c
c	deallocate(rv1_v,rv2,rv2_v,rv3,rv3_v,rv4,rv4_v,rv5,rv5_v)
c	deallocate(sn,sn_v,dvi,dvi2,dvi_v,dvi2_v,dv1,dv1_v)
c	deallocate(sdvi,sdvi2,sdvi_v,sdvi2_v,sdv1,sdv1_v)
c	deallocate(dn,dn2,dn_v,dn2_v,sdn,sdn2,sdn_v,sdn2_v)
c	deallocate(di2,di2_v,di_s,di_t,di_st,di2_stc)
c	deallocate(di_sv,di_tv,di_stv,di2_stcv)
c
c---------�[�q�d���̗h�炬-----------
	write(907,5)dis_node
c--���삳��̕��@�̓d���̗h�炬�̘a--
	do ix=1,ixmax
	di2(-1) = di2(-1) + di2(ix)
	enddo
	write(907,1000) di2(-1)
	write(907,1001) var
	write(907,1002) vave_all
	write(907,1003) n_all
	write(907,1004) sqrt(dis_node)/di_node
	write(907,1005) var/vave_all/sqrt(n_all)/sqrt(spnum)
	write(907,1006) spnum
c	call simpson(ixmax,di2,idx)
	write(*,*)"L1=",L1,"L2=",L2
	stop
	endif
	jt = jt + 1
c
	return
	end subroutine