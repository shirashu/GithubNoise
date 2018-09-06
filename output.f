	subroutine output(
     &			dt,spnum,de,jpnum,ict,
     &			dx,dz,xmax,zmax,lhet,iarea,
     &			cxpart1,cxpart2,lxpart1,lxpart2,
     &			cxpole1,cxpole2,lnpole1,lnpole2,jspsum,melpos,
     &			hhm,hm,af2,af4,eps,ec,p,kp,u,cn,hef_mesh,
     &			cur,hef_scat,count,hcss,cput,
     &			nava,cnava,n_scat,count_scat,
     &			balis_scat,balis_n,balis_all,
     &			ccs,ccs2,cncs,cncs2,allback_scat,back_scat,		!08/1/28
     &			tel1,tel2,tel3,efermi1,efermi2,efermi3,		!�|�ݒǉ�
     &			n_scat_p,n_scat_n,sscnt,					!�|�ݒǉ�
     &			avsumtel1,avsumconc1,avsumtei11,
     &            epA,epB,epC,epA2,epB2,eg,ecr,		!120126homma
     &            avsumteiA,avsumconcA,			!�|�ݒǉ�
     &            xxx,vvv,basho_reflection,basho_roughness,
     &            split,delta,lambda,count_reflection,
     &            count_roughness,average,
     &			pass,x_mean_free_path_sum,x_mean_free_path_count,			!120817sato		
     &			mean_free_path_sum,mean_free_path_sum2,
     &			roughness1_countx,roughness1_counte,bscat_xflag,fix_u,		!15/1/2takahashi
c############circuit 2014/12/01(takahashi)###################################
     &			IDS1_stack,IG1_stack,ISS1_stack,istp,i_or_c,jc_on)		

	implicit none
c
	include 'arraysize.fi'
c############circuit 2014/12/01(takahashi)###################################
	include 'Circuit.fi'
	real q
	parameter(q   = 1.60219e-19)
c
c---�ϐ��z��p�p�����[�^---
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c
c--- �V�~�����[�V�������� ---
	common /outp/joutc,jouti,jouta,jsave,jheat,jcput,joutpi,sw
	integer	joutc,jouti,jouta,jsave,jheat,jcput,joutpi,sw(8)
	integer	hcss
c
	common /step/jinit,jstat,jstep,jdisp
	integer	jinit,jstat,jstep,jdisp

c---��{�p�����[�^---
	real	dt,spnum,de(nenergy)
	integer	jpnum,ict
c---�f�o�C�X�\��---
	real	dx,dz,xmax,zmax
	integer(2),dimension (nlayer)	::lhet
	integer(1),dimension (nlayer)	::iarea
	real,	dimension (npart)	:: cxpart1, cxpart2
	integer(2),dimension (npart)	::lxpart1,lxpart2
c---�d��---
	real,	dimension (npole)	:: cxpole1,cxpole2
	integer(2),dimension (npole)	:: lnpole1,lnpole2
	integer(4),dimension (npole)	:: jspsum
	integer(1),dimension (npole)	:: melpos
c---�̈�ʃp�����[�^---
	real,	dimension (nvalley,narea)	:: hhm,hm,af2,af4
	real,	dimension (narea)	:: eps
	real,	dimension (nvalley,narea) :: ec,eg
c---���q���---
	real	p(6,npmax)
	integer(1),dimension (3,npmax)	:: kp
c---�f�o�C�X�����---
	real,	dimension (0:nx,0:nz)		:: u,cn
	real,	dimension ((nx+1)*(nz+1))	:: hef_mesh
	real,	dimension (nx+1,nscat,nvalley)	::hef_scat
	real,	dimension (0:nx,0:nz)	:: epA,epB,epC,epA2,epB2	!120126homma
c---���o�͗p�ϐ�---
	real	cur(npole)
	real(8)	cput(6)
c---�Փ˓d��---
	integer,dimension (0:nx,0:nz,0:nvalley) :: nava
	integer,dimension (0:nx,0:nvalley) :: cnava
	integer,dimension (0:nx,nscat,nvalley) :: n_scat
	integer,dimension (0:nx) :: count_scat
c-----�k�ތ���-----
	integer,dimension (0:nx,nscat,nvalley) :: n_scat_p !(+)yama071223
	integer,dimension (0:nx,nscat,nvalley) :: n_scat_n !(-)yama071223
	real,	dimension (0:nx,0:nz)	:: tel1,efermi1
	real,	dimension (0:nx,0:nz)	:: tel2,efermi2
	real,	dimension (0:nx,0:nz)	:: tel3,efermi3
c
	real,	dimension(0:nx,0:nz)	::	avsumtel1
	real,	dimension(0:nx,0:nz)	::	avsumconc1
	real avsumconcA(0:nx)          !101221�d�q�Z�x�`���l������2
	real(8),	dimension(0:nx,0:nz)	::	avsumtei11
	real(8) avsumteiA(0:nx)                 !101221�d�q�G�l���M�[�`���l������2
	real,save,allocatable 	:: av_tel1(:,:)
	real,save,allocatable 	:: av_conc1(:,:)
	real(8),save,allocatable 	:: av_tei11(:,:)
c
c---kd�Ɋւ���p�����[�^---	
	real adkx(0:nx,0:nz,0:nvalley)			!2014/12/01(takahashi)
	real adky(0:nx,0:nz,0:nvalley)
	real adkz(0:nx,0:nz,0:nvalley)
c
c-----yama�ǉ� 071020-----
	real,save,allocatable 	:: avef1(:,:),avtel1(:,:)
	real,save,allocatable 	:: avef2(:,:),avtel2(:,:)
	real,save,allocatable 	:: avef3(:,:),avtel3(:,:)
c
	integer ie2
	integer nemax4
	integer cff,bscat_xflag,fix_u					!�d���ӂ���v�Z�t���O(takahashi)
	real ecnt(10,1,nvalley,250)
	real ecr(7,int(nemax/4))		!120126homma
c
	real(8),save,allocatable :: cn_scat_p(:,:,:) !071223yama
	real(8),save,allocatable :: cn_scat_n(:,:,:) !071223yama
c
	integer,dimension (0:nx,nscat,4) :: sscnt
c
c----�o���X�e�B�b�N�̌v�Z----------
	integer	balis_scat,balis_n,balis_all
c-----�U���p�̏W�v-------------------!08/1/21
	real(8),dimension (0:nx) ::	ccs,cncs
	real(8),dimension (0:nx,nscat) ::	ccs2,cncs2
c-----����U���̏W�v-------------------!08/1/28
	real(8),dimension (0:nx,0:nscat) ::	allback_scat
	real(8),dimension (0:nx,0:nscat) ::	back_scat			
c##########circuit 2014/12/01(takahashi)###############################
	real IG1_stack(-1:jtp00),IDS1_stack(-1:jtp00),ISS1_stack(-1:jtp00)
	integer	istp,i_or_c,jc_on
c---���[�J���ϐ�---
	real,	dimension (:,:,:),allocatable::hef_fin	!���M�̃R�����g�A�E�g����
	real,	dimension (:),allocatable::hef_allscat	!���݂͎g���Ă��Ȃ�
	integer	count
c	integer	i,j,k,n
	integer i,n					!j,k�͌��ݎg���Ă��Ȃ�
	integer	ix,iz,ixiz			!ixiz�͌��݂͎g���Ă��Ȃ�
	integer	nxnz				!���M�֌W�̔z�񐔁B�������A�R�����g�A�E�g�̏ꏊ
	integer iscat
	real(8)	sk,ei
	integer	kv,ken,kl,ka,ie
	integer,save :: ncount			!save�Y��@
	real	ddt
	character(80) form
	real(8)	cpusum
	real	buff
c	
	real(8),save,allocatable 	:: avu(:,:),avcn(:,:)
	real(8),save,allocatable 	:: avu2(:),avcn2(:)
	real(8),save,allocatable	:: efield(:)
	real(8)	gk,sq,vave
	real(8),save,allocatable	:: v(:)
	real(8),save,allocatable 	:: vnmesh(:,:,:)
	integer(8),save,allocatable :: nvmesh(:,:)
	real(8),allocatable			:: vnout(:,:,:)
	real(8),save,allocatable 	:: vnx(:,:),cvnx(:,:)
	integer(8),save,allocatable :: nvx(:),cnvx(:)
	real(8),allocatable			:: vnxout(:,:),cvnxout(:,:)
	real(8),save,allocatable 	:: enmesh(:,:)
	integer(8),save,allocatable :: kvmesh(:,:,:)
	integer(8),save,allocatable :: cnmesh(:,:)
	real(8),save,allocatable 	:: enx(:),cenx(:)
	integer(8),save,allocatable :: kvx(:,:),ckvx(:,:)
	integer(8),save,allocatable :: cnx(:),ccnx(:)
	real(8),allocatable			:: dkvx(:,:),cdkvx(:,:)
	real(8),allocatable			:: kvratio(:,:,:)
	real(8),allocatable			:: dkvmesh(:,:,:)
	real(8),allocatable			:: dcnmesh(:,:)
	integer(8),save,allocatable :: iener(:,:)	!�G�l���M�[�ʗ��q��
	real(8),allocatable			:: dener(:,:)
	real(8),save	::	count2					!save�Y��
	real(8),save,allocatable :: vnava(:,:)
	real(8),save,allocatable :: cn_scat(:,:,:)

	real(8),save,allocatable 	:: en(:,:),cen(:,:)	!���b�V���ʕ��σG�l���M�[
	integer(8),save,allocatable :: dn(:),cdn(:)	!���b�V���ʗ��q��

	real(8),save,allocatable	:: sk2(:),sq2(:),ei2(:)
	real(8),save,allocatable	:: en2(:,:),cen2(:,:)
	real(8), dimension (0:nx) :: cn_nx

	real(8),save,allocatable 	:: cvnxp(:),cvnxm(:)	!07/11/14 ��[
	integer(8),save,allocatable :: cnvxp(:),cnvxm(:)		!07/11/14 ��[
	real(8),allocatable			:: cvnxoutp(:),cvnxoutm(:)
	real(8),save	::	count3								!save�Y��
	integer	dist_ix			!07/11/22 ��[
	integer,save,allocatable	:: dist_kxp(:,:),dist_kxm(:,:)	!07/11/22 ��[
	real,save	::	dkx											!07/11/22 ��[
	real balis_freq,balis_per									!07/11/22 ��[

	real(8),save,allocatable 	:: ccs_all(:)	!08/1/21 ��[
	real(8),save,allocatable 	:: ccs_scat(:,:)	!08/1/21 ��[
c
c	sw(1)...�d��
c	sw(2)...�G��
c	sw(3)...�|�e���V�����E�d�q�Z�x,2-DEG�w��ԁE�o�C�i���[
c	sw(4)...�G�l���M�[�E�J��L��
c	sw(5)...�d�q���x
c	sw(6)...�Փ˓d���p�x�E�U���p�x
c	sw(7)...���M��
c	sw(8)...cpu����
c

c------���t�l�X�U��--J.R. Watling���f��(2012�N�t����)-------------------------
      real delta,lambda,average
	integer split,count_reflection,count_roughness,
     &        vvv,xxx 
	integer,dimension (0:nx,2) :: basho_roughness,basho_reflection
      integer yyy

c-----	�̈�ʉ߂��闱�q���J�E���g120817sato
	real(8),dimension (10) :: pass
	real	x_mean_free_path_sum
	real	mean_free_path_sum
	integer x_mean_free_path_count
	real,dimension (0:nx,3) :: 	mean_free_path_sum2		!���ώ��R�s��x���W	 1:�J�E���^,2:x�������ώ��R�s��,3:���ώ��R�s��

c------���t�l�X�U��--R.P. Joshi���f��(2012�N�H)-------------------------
	integer,dimension (0:nx,narea,2) :: roughness1_countx			!x�����ɑ΂��郉�t�l�X�U���̐�(�`���l���㉺)
	integer,dimension (nemax,narea,2) :: roughness1_counte		!�G�l���M�[�ɑ΂��郉�t�l�X�U���̐�(�`���l���㉺)	

c ----( �d�� )----
	if(sw(1).gt.0)then
	if((jouti.eq.1).or.
     &	(modulo(ict,joutc).eq.(joutc-1)))then
		ddt = dt*jouti
		call gauss(
     &				ddt,spnum,jpnum,
     &				dx,dz,xmax,lhet,iarea,
     &				cxpart1,cxpart2,lxpart1,lxpart2,
     &				cxpole1,cxpole2,lnpole1,lnpole2,melpos,
     &				hm,hhm,af4,eps,p,kp,u,cur)
		form = '(I8,3('','',E11.4))'
		write(7 ,form) ict,(cur(i),i=3,1,-1)	!current.txt
		form = '(3('' '',E11.4))'
		write(45 ,form) (cur(i),i=3,1,-1)		!Y,Sparameter�pcurrent_t.txt
		form = '(I8,5('','',I8))'
		write(11,form) ict,(jspsum(i),i=1,npole)	!jspsum.txt
		if(ict.ge.0)then						!Vc�͖{�񂵂���o�͂��J�n���Ă���
c			�t�@�C�� �� unit=31 ... 'current_d.txt'
c			�t�@�C�� �� unit=32 ... 'current_s.txt'
			write(31,*)cur(3)		!�h���C���d��
			write(32,*)cur(1)		!�\�[�X�d��
			write(100,*) cur(3),'	',cur(2),'	',cur(1)     !vc.txt(DGS)
		endif
c###################circuit 2014/12/01(takahashi)#########################
	if(((jc_on.eq.1).or.(jc_on.eq.2).or.
     &             (jc_on.eq.3).or.(jc_on.eq.4)))then
		if(istp.gt.0)then
			IDS1_stack(istp) = cur(3)		!*Wg	drain	
			IG1_stack(istp)  = cur(2)		!*Wg	gate
			ISS1_stack(istp) = cur(1)		!*Wg	source
		endif
	endif
!	�ǉ�(2005/09/30 Hara)
! ----( �d���G���̃f�o�C�X�����z�v�Z )----
	if((sw(2).gt.0).and.(cff.ne.1))then
		call current_fluctuation(jpnum,dx,dz,xmax,zmax,
     &fix_u,ec,cur,dt,cxpole1,cxpole2,
     & p,kp,hhm,hm,af4,af2,iarea,adkx,adky,adkz,spnum,cff,bscat_xflag)	
	endif
c
	endif
	endif
c
c	----( �|�e���V�����E�d�q�Z�x�o�� )----
	if(sw(3).gt.0)then
		if (.not. allocated(avu)) then	!�ŏ��������s
			allocate(avu(0:nx,0:nz),avcn(0:nx,0:nz),
c-----�k�ތ���-----
     &				avef1(0:nx,0:nz),avtel1(0:nx,0:nz),
     &				avef2(0:nx,0:nz),avtel2(0:nx,0:nz),
     &				avef3(0:nx,0:nz),avtel3(0:nx,0:nz) )
c------------------
c-----08/8/6 �|��-----
	allocate(av_tel1(0:nx,0:nz))
	allocate(av_conc1(0:nx,0:nz))
	allocate(av_tei11(0:nx,0:nz))
	av_tel1 = 0.0
	av_conc1 = 0.0
	av_tei11 = 0.0
c---------------------
			avu  = 0.0
			avcn = 0.0
			ncount = 0
c-----�k�ތ���-----
			avef1 = 0.0; avef2 = 0.0; avef3 = 0.0
			avtel1 = 0.0; avtel2 = 0.0; avtel3 = 0.0
c------------------
		endif
		avu  = avu+u
		avcn = avcn+cn
		ncount= ncount+1
c-----�k�ތ���-----
		avef1 = avef1 + efermi1
		avef2 = avef2 + efermi2
		avef3 = avef3 + efermi3
		avtel1 = avtel1 + tel1
		avtel2 = avtel2 + tel2
		avtel3 = avtel3 + tel3
		av_tel1 = av_tel1 + avsumtel1	!08/8/6 �|�� ���J
		av_conc1=av_conc1 + avsumconc1	!08/8/6 �|�� ���J
		av_tei11=av_tei11 + avsumtei11	!08/8/6 �|�� ���J
c	-------------
c
c	----( �o�߃|�e���V�����E�d�q�Z�x�o�� )----
		if(modulo(ict,joutpi).eq.(joutpi-1))then
			!unit=25: potential.txt
			rewind 25 !; write(25,*)nx,nz,ict
			write(25,'(f12.7)') sngl(u)
			!unit=26: density.txt
			rewind 26 !; write(26,*)nx,nz,ict
			write(26,'(e15.7)') sngl(cn)
c	----( 2DEG�����|�e���V�����E�d�q�Z�x�o�� )----
			if(modulo(ict,(jouta)).eq.0)then
				rewind	27
				rewind	28
			endif
			write(27,'(I7,x)')ict
			write(28,'(I7,x)')ict
c			iz = lhet(nlayer-1)	!�v�C��
			iz = lhet(nchannel2)	
			!unit=27: potential2d.txt
c			write(27,'(f12.7)') (sngl(u(0:nx,iz)))	!07/03/14
			!unit=28: density2d.txt
c			write(28,'(e15.7)') (sngl(cn(0:nx,iz)))
		endif
c	---( ���ϗ��q���x, �|�e���V����, �d�E���x, �񎟌��d�q�Z�x�o�� )----
		if(modulo(ict,jouta).eq.(jouta-1))then
			! ----( �|�e���V�������Ϗo�� )----
			avu    = avu/ncount
			!unit=38: potential_ave.txt
			write(38,*)nx,nz,ict
			write(38,'(f16.7)')avu
c
			! ----( �d�E���x���Ϗo�� )----
			allocate(avu2(0:nx))
			allocate(efield(1:nx))
			avu2=0.0
			efield=0.0
c			do iz=lhet(nlayer-4),lhet(nlayer-1)	!channel 2006/12/22 3-6 <- nlayer 7
			do iz=lhet(nchannel1),lhet(nchannel2)	!channel 2011/04/08
			do ix=0,nx
			avu2(ix)=avu2(ix)+avu(ix,iz)
			enddo
			enddo
			do ix=1,nx
			efield(ix) = avu2(ix)-avu2(ix-1)
			enddo
			!unit=65: efield.txt
			rewind 65
			write(65,*)nx,dx,ict
c			write(65,'(e16.7)')efield/(lhet(nlayer-1)-lhet(nlayer-4)+1)/dx !channel 2006/12/22
	write(65,'(e16.7)')efield/(lhet(nchannel2)-lhet(nchannel1)+1)/dx !channel 2011/04/08
			deallocate(avu2,efield)
c
			! ----( ���q�Z�x���Ϗo�� )----
			avcn    = avcn/ncount
			!unit=39: cn_ave.txt
			write(39,*)nx,nz,ict
			write(39,'(e16.7)')avcn
c
			! ----( �񎟌����q�Z�x���z���Ϗo�� )----
			allocate(avcn2(0:nx))
c			do iz=lhet(nlayer-4),lhet(nlayer-1)	!channel 2006/12/22
			do iz=lhet(nchannel1),lhet(nchannel2) !channel 2011/04/08
			do ix=0,nx
			avcn2(ix)=avcn2(ix)+avcn(ix,iz)
			enddo
			enddo
			!unit=66: avcn2.txt
			rewind 66
			write(66,*)nx,dx,ict
			write(66,'(e16.7)')avcn2*dz/1.0e4	!cm-2
			deallocate(avcn2)
c
			!unit=9: state.txt
			do iz = 0, nz
			do ix = 0, nx
				write(9,*) sngl(avu(ix,iz)),',',sngl(avcn(ix,iz))
			enddo
			enddo
c
c-----�k�ތ���-----
	if(ict.gt.0)then
c-----(�t�F���~���x�����Ϗo��)-----
			avef1 = avef1 / ncount
			avef2 = avef2 / ncount
			avef3 = avef3 / ncount
			write(130,*)nx,nz,ict
			write(130,'(f16.7)')avef1
			write(131,*)nx,nz,ict
			write(131,'(f16.7)')avef2
			write(132,*)nx,nz,ict
			write(132,'(f16.7)')avef3
c
c-----(�d�q���x���Ϗo��)-----
			avtel1 = avtel1 / ncount
			avtel2 = avtel2 / ncount
			avtel3 = avtel3 / ncount
			write(133,*)nx,nz,ict
			write(133,'(f16.7)')avtel1
			write(134,*)nx,nz,ict
			write(134,'(f16.7)')avtel2
			write(135,*)nx,nz,ict
			write(135,'(f16.7)')avtel3
c-------------------
c-----08/8/6 �|��-----
			av_tel1 = av_tel1 / ncount
			av_conc1 = av_conc1 / ncount
			av_tei11 = av_tei11 / ncount
			write(150,*)nx,nz,ict			!avsumtel1.txt
			write(150,'(f16.7)')av_tel1		!avsumtel1.txt
			write(151,*)nx,nz,ict			!avsumconc1.txt
			write(151,'(e16.7)')av_conc1	!avsumconc1.txt
			write(152,*)nx,nz,ict			!avsumtei11.txt
			write(152,'(f16.7)')av_tei11	!avsumtei11.txt
	endif
c	----------------
c
			avu  = 0.0;	avcn = 0.0
			ncount= 0
c
c-----�k�ތ���-----
			avef1 = 0.0; avtel1 = 0.0
			avef2 = 0.0; avtel2 = 0.0
			avef3 = 0.0; avtel3 = 0.0
			av_tel1 = 0.0; av_conc1 = 0.0; av_tei11 = 0.0	!08/8/6 �|��
c-------------------
c
			if(modulo(ict,jsave).eq.(jsave-1))then
				rewind	40;	write(40)	npmax,xmax,zmax,spnum
				rewind	41;	write(41)	jpnum,p,kp
			endif
		endif
! �ǉ�--- �ϑ����O�̃f�[�^(05/09/30 Hara)
		if(ict.eq.-1)then
		!density_bef.txt �o��
		write(70,'(e15.7)') cn
		endif
	endif
c
c ----( �G�l���M�[�W�v )----
	if(sw(4).gt.0)then
	if(modulo(ict,jouti).eq.(jouti-1))then
		if (.not. allocated(enmesh)) then
			allocate(enmesh(0:nx,0:nz))
			allocate(kvmesh(0:nx,0:nz,nvalley))
			allocate(cnmesh(0:nx,0:nz))
			allocate(enx(0:nx),cenx(0:nx))
			allocate(kvx(0:nx,nvalley),ckvx(0:nx,nvalley))
			allocate(cnx(0:nx),ccnx(0:nx))
			allocate(iener(nenergy,nemax))
			allocate(en(0:nx,0:2),cen(0:nx,0:2))
			allocate(dn(0:nx),cdn(0:nx))
			allocate(sk2(1:3),sq2(1:3),ei2(1:3))	!07/2/20
			allocate(en2(0:nx,1:3),cen2(0:nx,1:3))	!07/2/20	
			enmesh = 0.0; kvmesh = 0; cnmesh = 0
			count2 =0.0
			iener  = 0
			enx = 0.0	;cenx = 0.0
			kvx = 0		;ckvx = 0
			cnx = 0		;ccnx = 0
			en	= 0.0	;cen = 0.0
			dn = 0		;cdn = 0
			en2=0.0     ;cen2=0.0
			ecnt = 0.0	;nemax4=nemax/4			!yama�ǉ�07/10/22
		endif

		rewind	170;rewind	175
c		;rewind	171;rewind	172;rewind	173;rewind	174
c		rewind  176;rewind  177;rewind  178;rewind  179;
		rewind	180;rewind	181;rewind	182
		rewind	183;rewind	184;rewind  185
		rewind  186;rewind  187
		rewind  188;rewind  189

		do n=1,jpnum
			if(kp(1,n).eq.0)cycle
			kv	= kp(1,n)
			ken	= kp(2,n)
			kl	= kp(3,n)
			ka = iarea(kl)
			sk = p(1,n)*p(1,n)+p(2,n)*p(2,n)+p(3,n)*p(3,n)
			sk2(1:3) = p(1:3,n)*p(1:3,n)		!07/2/20
			if(af4(kv,ka).ne.0.0)then
				sq = sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*sk)
				ei=(sq-1.0)/af2(kv,ka)+ec(kv,ka)
				sq2(1:3) = sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*sk2(1:3))		!07/2/20	
   				ei2(1:3) = (sq2(1:3)-1.0)/af2(kv,ka) !07/2/20
			else
				ei=hhm(kv,ka)*sk+ec(kv,ka)
				ei2(1:3) = hhm(kv,ka)*sk2(1:3)	!07/2/20	
			endif
c
			ix = max(min(ifix(p(5,n)/dx+0.5),nx),0)
			iz = max(min(ifix(p(6,n)/dz+0.5),nz),0)

c----------�Ō�̃`���l���̓d�q�̃G�l���M�[���z--------2010/12/21hisa	 ���b�V����1nm
              if(modulo(ict,jouta).eq.(jouta-1))then
	           if(kv.eq.1) then
c                     if((iz.ge.38).and.(iz.le.47))		then
					if((lhet(nchannel1)+1.le.iz).and. !channel 11/05/26��
     &			 		(iz.le.lhet(nchannel2))) then 		
	                      write(170,*) p(5,n),ei-epA2(ix,iz),n !iz_ene.txt
					!---------- 120126homma
						if(modulo(n,2).eq.0)then
							write(173,*)p(5,n),ei-epA(ix,iz),n
c						endif
c						if(modulo(n,3).eq.0)then
							write(174,*)p(5,n),ei-epA2(ix,iz),n
						endif
					!---------- 
					endif
					if((iz.ge.22).and.(iz.le.23))		then
						write(171,*)p(5,n),ei-epA2(ix,iz),n	
					elseif((iz.ge.1).and.(iz.le.29))		then
						write(172,*)p(5,n),ei-epA2(ix,iz),n
 					endif
c                   
	               if(ix.eq.110)		then	!�v�C��11/05/26�� ix=50��120126homma ix=110
                         write(180,*) p(6,n),ei-epA2(ix,iz),n
	               elseif(ix.eq.210)		then
                         write(181,*) p(6,n),ei-epA2(ix,iz),n		
	               elseif(ix.eq.250)		then
                         write(182,*) p(6,n),ei-epA2(ix,iz),n
	               elseif(ix.eq.290)	then
                         write(183,*) p(6,n),ei-epA2(ix,iz),n
	               elseif(ix.eq.450)	then
                         write(184,*) p(6,n),ei-epA2(ix,iz),n
	               endif
c
	           elseif((kv.eq.2).or.(kv.eq.3))		then
c                     if((iz.ge.38).and.(iz.le.47))		then
					if((lhet(nchannel1)+1.le.iz).and. !channel 11/05/26��
     &					(iz.le.lhet(nchannel2))) then 
	                      write(175,*) p(5,n),ei-ec(kv,ka)-epB2(ix,iz),n !iz_ene_L.txt 120126homma
						!---------- 120126homma
						if(modulo(n,2).eq.0)then
							write(176,*)p(5,n),ei-ec(kv,ka)-epB2(ix,iz),n
						endif
						if(modulo(n,3).eq.0)then
							write(177,*)p(5,n),ei-ec(kv,ka)-epB2(ix,iz),n
						endif
						!---------- 
					endif
c
	               if(ix.eq.110)		then	!�v�C��11/05/26����120126homma ix=110
                         write(185,*) p(6,n),ei-ec(kv,ka)-epB2(ix,iz),n	!120126homma
	               elseif(ix.eq.210)		then
                         write(186,*) p(6,n),ei-ec(kv,ka)-epB2(ix,iz),n	!120126homma
	               elseif(ix.eq.250)		then
                         write(187,*) p(6,n),ei-ec(kv,ka)-epB2(ix,iz),n	!120126homma
	               elseif(ix.eq.290)	then
                         write(188,*) p(6,n),ei-ec(kv,ka)-epB2(ix,iz),n	!120126homma
	               elseif(ix.eq.450)	then
                         write(189,*) p(6,n),ei-ec(kv,ka)-epB2(ix,iz),n	!120126homma
	               endif
	           endif
	       endif
c------------------------------------------------------

			! --- ���q�G�l���M�[���z�W�v
			enmesh(ix,iz)=enmesh(ix,iz)+ei		!���b�V�����G�l���M�[���a
			enx(ix)=enx(ix)+ei				! ���b�V�����G�l���M�[���a

			! --- �J��L���W�v
			kvmesh(ix,iz,kv)=kvmesh(ix,iz,kv)+1		!���b�V�����J�ʗ��q��
			kvx(ix,kv)=kvx(ix,kv)+1				!�J�ʗ��q���c�������a

			! --- ���b�V�������q���J�E���g
			cnmesh(ix,iz)=cnmesh(ix,iz)+1	!���b�V�������q���a
			cnx(ix)=cnx(ix)+1				! ���q���c�������a
c
c-----�k�ތ���-----
c-----�G�l���M�[�ʗ��q���J�E���g-----
c			ie2=max(min((nint(ei/(de(2)*4))+1),(nemax4)),1)
	if(kv.eq.1)then		!120126homma
		ie2=max(min((nint((ei-epA(ix,iz)+
     &	              (u(ix,iz)-0.7*(eg(1,ka)-eg(1,2)))) !!11/06/29��
     &                    /(de(2)*4))+1),(nemax4)),1)	!101220
	elseif(kv.eq.2)then		!120126homma
		ie2=max(min((nint((ei-ec(kv,ka)-epB(ix,iz)+
     &	              (u(ix,iz)-0.7*(eg(1,ka)-eg(1,2))))
     &                    /(de(2)*4))+1),(nemax4)),1)
	elseif(kv.eq.3)then		!120126homma
		ie2=max(min((nint((ei-ec(kv,ka)-epC(ix,iz)+
     &	              (u(ix,iz)-0.7*(eg(1,ka)-eg(1,2))))
     &                    /(de(2)*4))+1),(nemax4)),1)
     	endif
c		if(nlayer-4.lt.kl.and.kl.le.nlayer-1) then  �`���l�����ł݂�̂�����H
			if((ix.eq.10).and.(iz.eq.36))   ! iz=30�֕ύX(�ύX�Oiz=42)
     &   					ecnt(1,1,kv,ie2) = ecnt(1,1,kv,ie2) + 1
			if((ix.eq.100).and.(iz.eq.36)) 
     &   					ecnt(2,1,kv,ie2) = ecnt(2,1,kv,ie2) + 1
			if((ix.eq.210).and.(iz.eq.36)) 
     &   					ecnt(3,1,kv,ie2) = ecnt(3,1,kv,ie2) + 1
			if((ix.eq.215).and.(iz.eq.36)) 
     &   					ecnt(4,1,kv,ie2) = ecnt(4,1,kv,ie2) + 1
			if((ix.eq.250).and.(iz.eq.36)) 
     &   					ecnt(5,1,kv,ie2) = ecnt(5,1,kv,ie2) + 1
			if((ix.eq.285).and.(iz.eq.36)) 
     &   					ecnt(6,1,kv,ie2) = ecnt(6,1,kv,ie2) + 1
			if((ix.eq.290).and.(iz.eq.36)) 
     &   					ecnt(7,1,kv,ie2) = ecnt(7,1,kv,ie2) + 1
			if((ix.eq.400).and.(iz.eq.36)) 
     &   					ecnt(8,1,kv,ie2) = ecnt(8,1,kv,ie2) + 1
			if((ix.eq.490).and.(iz.eq.36)) 
     &   					ecnt(9,1,kv,ie2) = ecnt(9,1,kv,ie2) + 1
c			if((ix.eq.300).and.(iz.eq.36)) 
c     &   					ecnt(10,1,kv,ie2) = ecnt(10,1,kv,ie2) + 1
c
			! --- �G�l���M�[�e�[�u�������z�W�v
			ie=max(min((nint(ei/de(ken))+1),nemax),1)
			iener(ken,ie)=iener(ken,ie)+1	

			!---�`���l����---!!!!!!!!!!!!!!!!!!!!!!!!
c			if(nlayer-4.lt.kl.and.kl.le.nlayer-1) then	!channel 2006/12/22
			if(nchannel1.lt.kl.and.kl.le.nchannel2) then !channel 2011/04/08
				cenx(ix)=cenx(ix)+ei			!���b�V�����G�l���M�[���a
				ckvx(ix,kv)=ckvx(ix,kv)+1		!���b�V�����J�ʗ��q��
				ccnx(ix)=ccnx(ix)+1			!���b�V�������q���a
			endif

c			---�ꎟ�����σG�l���M�[enx---
			en(ix,1) = en(ix,1) + ei-ec(kv,ka)
			en(ix,2) = en(ix,2) + ec(kv,ka)
			en(ix,0) = en(ix,0) + ei
			dn(ix) = dn(ix) + 1
	        en2(ix,1) = en2(ix,1) +ei2(1)
			en2(ix,2) = en2(ix,2) +ei2(2)
			en2(ix,3) = en2(ix,3) +ei2(3)
c			---�ꎟ���`���l�������σG�l���M�[cenx---
c			if(nlayer-4.lt.kl.and.kl.le.nlayer-1)then	!channel 2007/2/2
			if(nchannel1.lt.kl.and.kl.le.nchannel2) then
			cen(ix,1) = cen(ix,1) + ei-ec(kv,ka)
			cen(ix,2) = cen(ix,2) + ec(kv,ka)
			cen(ix,0) = cen(ix,0) + ei
			cdn(ix) = cdn(ix) + 1
		    cen2(ix,1) = cen2(ix,1) +ei2(1)
			cen2(ix,2) = cen2(ix,2) +ei2(2)
			cen2(ix,3) = cen2(ix,3) +ei2(3)
			endif
c
		end do
c
		count2 = count2 + 1.0

		allocate(dener(nenergy,nemax))
		dener = dble(iener)/count2
		rewind 20
		do ie=1,nemax
			write(20,"(10(f10.3))")(dener(i,ie),i=1,nenergy)
		enddo
		deallocate(dener)
		iener = 0
c
c---( �G�l���M�[,�J��L���o�� )---
		if(modulo(ict,jouta).eq.(jouta-1))then
			forall(ix=0:nx,ccnx(ix).ne.0)			!���Ōv�Z����B���ɕʂ̗��q������
			cn_nx(ix)=dble(ccnx(ix))/count2
			end forall
c		---�ꎟ�����σG�l���M�[�o��enx---
			forall(ix=0:nx,dn(ix).ne.0)
				en(ix,0:2) = en(ix,0:2)/dble(dn(ix))	!07/1/22�@��[				
			end forall
			rewind 68
			do ix=0,nx
				write(68,'(3(F8.4))')en(ix,0:2)
			enddo
			forall(ix=0:nx,dn(ix).ne.0)	!channel 2007/2/2
				en2(ix,1:3) = en2(ix,1:3)/dble(dn(ix)) !07/1/22�@��[			
			end forall
			rewind 88
			do ix=0,nx
				write(88,'(3(F8.4))')en2(ix,1:3)
			enddo
c		---�ꎟ���`���l�������σG�l���M�[�o��cenx---
			forall(ix=0:nx,cdn(ix).ne.0)
				cen(ix,0:2) = cen(ix,0:2)/dble(cdn(ix))		!07/1/22�@��[
			end forall
			rewind 69
			do ix=0,nx
				write(69,'(3(F8.4))')cen(ix,0:2)
			enddo
			forall(ix=0:nx,cdn(ix).ne.0)	!channel 2007/2/2	!07/1/22�@��[
				cen2(ix,1:3) = cen2(ix,1:3)/dble(cdn(ix))		!07/1/22�@��[	
			end forall
			rewind 89
			do ix=0,nx
				write(89,'(3(F8.4))')cen2(ix,1:3)
			enddo
			en = 0.0; cen = 0.0
			dn = 0;   cdn = 0
			en2	= 0.0;	cen2 = 0.0
c			deallocate(en,dn)
c			deallocate(cen,cdn)
c
c---( �ꎟ���G�l���M�[,�J��L���o�� )---
			allocate(dkvx(0:nx,1:nvalley),cdkvx(0:nx,1:nvalley))
			dkvx = kvx
			cdkvx = ckvx
			!unit=49: dn_zave.txt
			rewind 49
			write(49,*)nx,ict
			write(49,'(e12.5)') dkvx/count2
			!unit=59: cdn_zave.txt
			rewind 59
			write(59,*)nx,ict
			write(59,'(e12.5)') cdkvx/count2
			forall(ix=0:nx,cnx(ix).ne.0)
				enx(ix) = enx(ix)/dble(cnx(ix))
			end forall
			do i = 1,nvalley
				forall(ix=0:nx,cnx(ix).ne.0)
					dkvx(ix,i)  =dkvx(ix,i)/dble(cnx(ix))
				end forall
			enddo
			forall(ix=0:nx,ccnx(ix).ne.0)
				cenx(ix)    =cenx(ix)  /dble(ccnx(ix))
			end forall
			do i = 1,nvalley
				forall(ix=0:nx,ccnx(ix).ne.0)
					cdkvx(ix,i)  =cdkvx(ix,i)/dble(ccnx(ix))
				end forall
			enddo
			!unit=46: kv_zave.txt
			rewind 46
			write(46,*)nx,nvalley,ict
			write(46,'(f10.7)')dkvx
			!unit=56: ckv_zave.txt
			rewind 56
			write(56,*)nx,nvalley,ict
			write(56,'(f10.7)')cdkvx
			!unit=47: en_zave.txt
			rewind 47
			write(47,*)nx,dx,ict
			write(47,'(e12.5)')enx
			!unit=57: cen_zave.txt
			rewind 57
			write(57,*)nx,dx,ict
			write(57,'(e12.5)')cenx
			enx=0.0 ;kvx=0 ;cnx=0
			cenx=0.0;ckvx=0;ccnx=0
			deallocate(dkvx,cdkvx)
c			deallocate(enx,cenx)
c			deallocate(kvx,ckvx)
c			deallocate(cnx,ccnx)

c---( �񎟌��G�l���M�[,�J��L���o�� )---
			allocate(dkvmesh(0:nx,0:nz,nvalley),dcnmesh(0:nx,0:nz))
			allocate(kvratio(0:nx,0:nz,nvalley))
			dkvmesh = kvmesh
			dcnmesh = cnmesh
			kvratio = 0.0
			forall(ix=0:nx,iz=0:nz,cnmesh(ix,iz).ne.0)
				enmesh(ix,iz)=enmesh(ix,iz)/dcnmesh(ix,iz)
				kvratio(ix,iz,1:nvalley)=dkvmesh(ix,iz,1:nvalley)
     &					/dcnmesh(ix,iz)
				dkvmesh(ix,iz,1:nvalley)=dkvmesh(ix,iz,1:nvalley)/count2
			end forall
			rewind	21;write(21,*)nx,nz,ict;write(21,'(e12.5)')enmesh 	!enmesh.txt
			rewind	23;write(23,*)nx,nz,ict;write(23,'(f10.7)')kvratio	!kvmesh.txt
			rewind	24;write(24,*)nx,nz,nvalley,ict;write(24,'(e12.5)')dkvmesh	!dnmesh.txt
			enmesh=0.0;vnmesh=0.0;kvmesh=0.0
			deallocate(dkvmesh)
			deallocate(kvratio)
c			deallocate(enmesh,kvmesh,cnmesh)

			! ----( ���q�Z�x���z���Ϗo�� )----
			dcnmesh = dcnmesh*spnum/dx/dz/count2
			!unit=37: density_ave.txt
			write(37,*)nx,nz,ict
			write(37,'(e16.7)')dcnmesh
			cnmesh = 0; count2 = 0;
			deallocate(dcnmesh)
c
c-----�k�ތ���-----
c-----�G�l���M�[�ʗ��q���o��-----
			rewind	120;rewind	121;rewind	122;rewind	123;rewind	124
			rewind	125;rewind	126;rewind	127;rewind	128
c															;rewind	129
			rewind	153		!120126homma
	do ix=1,10
	do kv=1,3
	do ie=1,nemax4
c			forall(ix=1:10,kv=1:3,ie2=1:nemax4)
			if(count2.ne.0)	 ecnt(ix,1,kv,ie2)=ecnt(ix,1,kv,ie2)/count2
c			endforall
	enddo
	enddo
	enddo
			do ie2 = 1,nemax4
				write(120,*) ecnt(1,1,1,ie2),ecnt(1,1,2,ie2),ecnt(1,1,3,ie2)
				write(121,*) ecnt(2,1,1,ie2),ecnt(2,1,2,ie2),ecnt(2,1,3,ie2)
				write(122,*) ecnt(3,1,1,ie2),ecnt(3,1,2,ie2),ecnt(3,1,3,ie2)
				write(123,*) ecnt(4,1,1,ie2),ecnt(4,1,2,ie2),ecnt(4,1,3,ie2)
				write(124,*) ecnt(5,1,1,ie2),ecnt(5,1,2,ie2),ecnt(5,1,3,ie2)
				write(125,*) ecnt(6,1,1,ie2),ecnt(6,1,2,ie2),ecnt(6,1,3,ie2)
				write(126,*) ecnt(7,1,1,ie2),ecnt(7,1,2,ie2),ecnt(7,1,3,ie2)
				write(127,*) ecnt(8,1,1,ie2),ecnt(8,1,2,ie2),ecnt(8,1,3,ie2)
				write(128,*) ecnt(9,1,1,ie2),ecnt(9,1,2,ie2),ecnt(9,1,3,ie2)
c				write(129,*) ecnt(10,1,1,ie2),ecnt(10,1,2,ie2),ecnt(10,1,3,ie2)
	write(153,*) ecr(1,ie2),ecr(2,ie2),ecr(3,ie2),ecr(4,ie2),ecr(5,ie2)	!120126homma
			enddo
			ecnt=0.0
c-------------------

		endif
	endif
	endif
c
c----( ���x���z�o�� )----
	if(sw(5).gt.0)then
	if(modulo(ict,jouti).eq.(jouti-1))then
		if (.not. allocated(v)) then
c			v_mesh...���b�V���ʑ��x�a, sv_mesh...���b�V���ʕ��ϑ��x
c			1...x����, 2...y����, 3...z����, 0...��Βl
			allocate(v(1:3))
 			allocate(vnmesh(0:3,0:nx,0:nz))
			allocate(nvmesh(0:nx,0:nz))
			allocate(vnx(0:3,0:nx),cvnx(0:3,0:nx))
			allocate(nvx(0:nx),cnvx(0:nx))
			allocate(cvnxp(0:nx),cvnxm(0:nx))		!07/11/14 ��[
			allocate(cnvxp(0:nx),cnvxm(0:nx))		!07/11/14 ��[
			allocate(dist_kxp(0:nx,0:500),dist_kxm(0:nx,0:500))
c	-----08/3/24 �|�� �\����傫������Ƃ��ɂ͏��c������c���Ƃ�-----
c			allocate(dist_kxp(0:nx,0:2000),dist_kxm(0:nx,0:2000))
			vnmesh = 0.0
			nvmesh = 0
			vnx = 0.0	;cvnx = 0.0
			nvx = 0		;cnvx = 0
			cvnxp = 0.0	;cvnxm = 0.0				!07/11/14 ��[
			cnvxp = 0	;cnvxm = 0					!07/11/14 ��[
			dkx = 5.0e-8							!2e7�ŋK�i��
			dist_kxp = 0;dist_kxm = 0	
			count3 =0.0								!08/8/6 �|��
		endif
c
		do n=1,jpnum
			sk = p(1,n)*p(1,n)+p(2,n)*p(2,n)+p(3,n)*p(3,n)
			kv = kp(1,n)
			kl = kp(3,n)
			ka = iarea(kl)
			gk = hhm(kv,ka)*sk
			v(1:3)  = p(1:3,n)*hm(kv,ka)/sqrt(1.0+af4(kv,ka)*gk)
			vave = sqrt(v(1)**2+v(2)**2+v(3)**2)
			ix = nint(p(5,n)/dx)
			iz = nint(p(6,n)/dz)
c			---���b�V���ʕ��ϑ��x---
			vnmesh(0,ix,iz)=vnmesh(0,ix,iz)+vave
			vnmesh(1:3,ix,iz)=vnmesh(1:3,ix,iz)+v(1:3)
c			---�ꎟ�����ϑ��x---
			vnx(0,ix)=vnx(0,ix)+vave
			vnx(1:3,ix)=vnx(1:3,ix)+v(1:3)
			! --- ���b�V���ʑ��x�W�v���J�E���g
			nvmesh(ix,iz)=nvmesh(ix,iz)+1
			nvx(ix)=nvx(ix)+1
			!---�`���l����---
c			if(nlayer-4.lt.kl.and.kl.le.nlayer-1) then
			if(nchannel1.lt.kl.and.kl.le.nchannel2) then !channel 11/04/07��
				cvnx(0,ix)=cvnx(0,ix)+vave
				cvnx(1:3,ix)=cvnx(1:3,ix)+v(1:3)
				cnvx(ix)=cnvx(ix)+1

				dist_ix = nint(p(1,n)*dkx)	!07/11/22 ��[
c---------------------------------------------------------- !07/11/14 ��[
c
c-----cvn xp(x):�v���X�����̑��x�������a
c-----cnv xp(x):�v���X�����̑��x�������q�����a
c-----cvn xm(x):�}�C�i�X�����̑��x�������a
c-----cnv xm(x):�}�C�i�X�����̑��x�������q�����a
c-----cvnxoutp(ix):�v���X�������x�o��
c-----cvnxoutm(ix):�}�C�i�X�������x�o��


				if(v(1).gt.0) then
					if(dist_ix.ge.100) dist_ix=100	!yama080105
					cvnxp(ix) = cvnxp(ix)+v(1)
					cnvxp(ix)=cnvxp(ix)+1
					dist_kxp(ix,dist_ix) = dist_kxp(ix,dist_ix) + 1	!07/11/22 ��[
				elseif(v(1).lt.0) then
					if(dist_ix.le.-100) dist_ix=-100	!yama080105
					cvnxm(ix) = cvnxm(ix)+v(1)
					cnvxm(ix)=cnvxm(ix)+1
					dist_ix = -dist_ix						!07/11/22 ��[
					dist_kxm(ix,dist_ix) = dist_kxm(ix,dist_ix) + 1
				endif
c----------------------------------------------------------
			endif
		enddo
		count3 = count3 + 1.0
c	
c	  ---�d�q���x�o��---
		if(modulo(ict,jouta).eq.(jouta-1))then
			allocate(vnxout(0:nx,0:3),cvnxout(0:nx,0:3))
			allocate(cvnxoutp(0:nx),cvnxoutm(0:nx))
			do ix=0,nx
c				if(nvx(ix).ne.0)then
				if(cnvx(ix).ne.0)then	!08/1/22	��[
					vnxout(ix,0:3)=vnx(0:3,ix)/dble(nvx(ix))
					cvnxout(ix,0:3)=cvnx(0:3,ix)/dble(cnvx(ix))
c					cvnxoutp(ix)=cvnxp(ix)/dble(cnvxp(ix))
c					cvnxoutm(ix)=cvnxm(ix)/dble(cnvxm(ix))
				else
					vnxout(ix,0:3)=0.0
					cvnxout(ix,0:3)=0.0
c					cvnxoutp(ix)=0.0
c					cvnxoutm(ix)=0.0
				endif

				if(cnvxp(ix).ne.0)then
					cvnxoutp(ix)=cvnxp(ix)/dble(cnvxp(ix))
				else
					cvnxoutp(ix)=0.0
				endif

				if(cnvxm(ix).ne.0)then
					cvnxoutm(ix)=cvnxm(ix)/dble(cnvxm(ix))
				else
					cvnxoutm(ix)=0.0
				endif
			enddo
			!unit=48: vn_zave.txt
			rewind 48
			write(48,*)nx,dx,ict
			write(48,'(e14.6)')vnxout
			!unit=58: cvn_zave.txt
			rewind 58
			write(58,*)nx,dx,ict
			write(58,'(e14.6)')cvnxout
			!unit=108: cvn_zave2.txt
			rewind 108
			write(108,*)nx,dx,ict
			do ix=0,nx
				write(108,'(2(e14.6))')cvnxoutp(ix),cvnxoutm(ix)
			enddo
			!unit=109: cdn_zave2.txt
			rewind 109
			write(109,*)nx,dx,ict
			do ix=0,nx
				write(109,'(2(e12.5))') cnvxp(ix)/count3,cnvxm(ix)/count3
			enddo
			!unit=111: kx_distribution.txt
			rewind 111
			do dist_ix = 0,100
				do ix=0,nx
					write(111,'(2(f10.5))') dble(dist_kxp(ix,dist_ix))/count3
     &							,dble(dist_kxm(ix,dist_ix))/count3
				enddo
			enddo

c			deallocate(vnx,cvnx)
c			deallocate(nvx,cnvx)
			deallocate(vnxout,cvnxout)
			deallocate(cvnxoutp,cvnxoutm)

			allocate(vnout(0:nx,0:nz,0:3))
			do iz=0,nz
			do ix=0,nx
				if(nvmesh(ix,iz).ne.0)then
					vnout(ix,iz,0:3)=vnmesh(0:3,ix,iz)/nvmesh(ix,iz)
				else
					vnout(ix,iz,0:3)=0.0
				endif
			enddo
			enddo
			rewind	22;write(22,*)nx,nz,ict;write(22,'(e14.6)')vnout	!vnmesh.txt
c			deallocate(v,vnmesh,nvmesh)
			deallocate(vnout)
			vnmesh = 0.0	!���Z�b�g
			nvmesh = 0
			vnx = 0.0	;cvnx = 0.0
			nvx = 0		;cnvx = 0
			cvnxp = 0.0	;cvnxm = 0.0
			cnvxp = 0		;cnvxm = 0
			count3 =0.0
			dist_kxp = 0;dist_kxm = 0

		endif
	endif
	endif
c
c----( �Փ˓d���p�x )----
c	nava...���b�V���ʏՓ˓d���p�x, vnava...�J�ʏՓ˓d���p�x
	if(sw(6).gt.0)then
	if(modulo(ict,jouta).eq.(jouta-1))then
	  if (.not. allocated(vnava)) then
 		allocate(vnava(0:nx,1:nvalley))
 		allocate(cn_scat(0:nx,1:nscat,nvalley))
		allocate(ccs_all(0:nx))
		allocate(ccs_scat(0:nx,1:nscat))
c-----yama�ǉ�-----
	 	allocate(cn_scat_p(0:nx,1:nscat,nvalley))
	 	allocate(cn_scat_n(0:nx,1:nscat,nvalley))
		vnava = 0.0
		cn_scat = 0.0
		ccs_all = 0.0
		ccs_scat = 0.0
c-----yama�ǉ�-----
		cn_scat_p = 0.0
		cn_scat_n = 0.0
	endif
		!---���b�V���ʏՓ˓d���p�x---
		!unit=61: nava.txt
		rewind 61
		write(61,'(I9)') nava(0:nx,0:nz,0)
		nava = 0
		!---�ꎟ���Փ˓d���p�x---
		!unit=60: ava_zave.txt
		rewind 60
		write(60,*)nx,dx,ict
		write(60,'(e12.4)') cnava(0:nx,0)*spnum/dx/dz/jouta/dt	!times/s/m3
		!unit=62: cnave.txt
		rewind 62
		write(62,*)nx,dx,ict
		write(62,'(I9)') cnava(0:nx,0) !times
		!---�J�ʏՓ˓d���p�x---
		do kv = 1, nvalley
			forall(ix=0:nx,cnava(ix,0).ne.0)
				vnava(ix,kv) = dble(cnava(ix,kv))/dble(cnava(ix,0))*100
			end forall
		enddo
		!unit=63: vnava.txt
		rewind 63
		write(63,*)nx,nvalley,ict
		write(63,'(f8.4)')vnava
		cnava=0
		vnava=0

		!----( �`���������U���p�x�o�� )----
		!---�ꎟ���U���p�x---
			do kv = 1, nvalley
				do iscat = 1, nscat
				  forall(ix=0:nx,count_scat(ix).ne.0)
c				  cn_scat(ix,iscat,kv) = dble(n_scat(ix,iscat,kv))
c   &								/dble(count_scat(ix))*spnum/dx/dz/jouta/dt
				  cn_scat(ix,iscat,kv) = dble(n_scat(ix,iscat,kv))
     &								*spnum/dx/dz/jouta/dt
c-----08/8/6 �|��-----
				  cn_scat_p(ix,iscat,kv) = dble(n_scat_p(ix,iscat,kv))
     &								*spnum/dx/dz/jouta/dt
				  cn_scat_n(ix,iscat,kv) = dble(n_scat_n(ix,iscat,kv))
     &							*spnum/dx/dz/jouta/dt
c---------------------
				  end forall
				enddo
			enddo
c		endif

		!unit=64: cnscat.txt
		rewind 64
		write(64,*)nx,dx,nscat,nvalley,ict
c-----�k�ތ���-----
		rewind 136
		rewind 137
		write(136,*)nx,dx,nscat,nvalley,ict
		write(137,*)nx,dx,nscat,nvalley,ict
c------------------
		form = '(e12.4, 2(I6))'
		do kv = 1, nvalley
		do iscat = 1, nscat
		do ix=0, nx
				write(64,form) cn_scat(ix,iscat,kv),n_scat(ix,iscat,kv),
     &			           count_scat(ix)
c-----�k�ތ���-----
				write(136,form) cn_scat_p(ix,iscat,kv),n_scat_p(ix,iscat,kv),
     &							count_scat(ix)
				write(137,form) cn_scat_n(ix,iscat,kv),n_scat_n(ix,iscat,kv),
     &							count_scat(ix)
c------------------
		enddo
		enddo
		enddo
c
c-----yama�ǉ�-----
		rewind 138
		rewind 139
		rewind 140
		rewind 141
		rewind 142
		rewind 143
		rewind 144
		rewind 145
		write(138,*)nx,dx,nscat,nvalley+1,ict
		write(139,*)nx,dx,nscat,nvalley+1,ict
		write(140,*)nx,dx,nscat,nvalley+1,ict
		write(141,*)nx,dx,nscat,nvalley+1,ict
		write(142,*)nx,dx,nscat,nvalley+1,ict
		write(143,*)nx,dx,nscat,nvalley+1,ict
		write(144,*)nx,dx,nscat,nvalley+1,ict
		write(145,*)nx,dx,nscat,nvalley+1,ict
		do iscat = 1, nscat
		do ix=0, nx
			write(138,*) float(sscnt(ix,iscat,1))/float(jouta)
			write(139,*) float(sscnt(ix,iscat,2))/float(jouta)
			write(140,*) float(sscnt(ix,iscat,3))/float(jouta)
			write(141,*) float(sscnt(ix,iscat,4))/float(jouta)
			write(142,*) sscnt(ix,iscat,1)
			write(143,*) sscnt(ix,iscat,2)
			write(144,*) sscnt(ix,iscat,3)
			write(145,*) sscnt(ix,iscat,4)
		enddo
		enddo
c
c----�o���X�e�B�b�N�̌v�Z----------
		balis_freq = float(balis_scat)/float(balis_all)
		balis_per =	float(balis_all-balis_n)/float(balis_all)*100
		
		rewind 110
		write(110,'(2(e12.5))') balis_freq,balis_per

		balis_all = 0; balis_scat = 0
		balis_n = 0; balis_freq = 0
		balis_per = 0
c--------------------------
c-----�U���p�̌v�Z---------
		do ix = 0, nx
			if(cncs(ix).ne.0)then
				ccs_all(ix)=ccs(ix)/cncs(ix)
			else
				ccs_all(ix)=0.0
			endif
			do iscat=1, nscat
				if(cncs2(ix,iscat).ne.0)then
					ccs_scat(ix,iscat)= ccs2(ix,iscat)/cncs2(ix,iscat)
				else
					ccs_scat(ix,iscat) = 0.0
				endif
			enddo
		enddo

		rewind 112
		rewind 113
		do ix=0, nx
				write(112,*) ccs_all(ix)
		enddo
		write(113,*)nx,dx,nscat,nvalley,ict
		do iscat = 1, nscat
		do ix=0, nx
				write(113,*) ccs_scat(ix,iscat)
		enddo
		enddo

		ccs_all = 0.0; ccs_scat = 0.0
		ccs = 0.0;ccs2 = 0.0
		cncs = 0.0;cncs2 = 0.0
c--------------------------
c-----����U���̃J�E���g--------
		rewind 114
		rewind 115
		rewind 116
		do ix=0, nx
				write(114,*) allback_scat(ix,0),back_scat(ix,0)
		enddo
		write(115,*)nx,dx,nscat,nvalley,ict
		write(116,*)nx,dx,nscat,nvalley,ict
		do iscat = 1, nscat
		do ix=0, nx
				write(115,*) allback_scat(ix,iscat)
				write(116,*) back_scat(ix,iscat)
		enddo
		enddo
		allback_scat = 0.0; back_scat = 0.0
c-------------------------------------
		n_scat = 0;count_scat=0
		cn_scat= 0;
		cn_nx=0
		n_scat_p = 0					!08/8/6 �|��
		n_scat_n = 0					!08/8/6 �|��
	endif
	endif
c
c%%%%%%%%%% ���M���֌W�o�� %%%%%%%%%%
c	�t�@�C�� �� unit=30 ... 'heat_generation.txt'
c	�t�@�C�� �� unit=33 ... 'heat_generation_mesh.txt'
c	�t�@�C�� �� unit=36 ... 'heat_generation_all.txt'
	if(sw(7).gt.0)then
	if((modulo(ict,jheat).eq.(jheat-1)).and.(hcss.eq.1))then

			buff=q*spnum/(dx*dz*dt*float(count))
			hef_mesh=hef_mesh*buff	!hef_mesh:�S�̔z��
			hef_scat=hef_scat*buff	!hef_scat:�S�̔z��
			rewind	30
			write(30,'(E15.7)') hef_mesh		!heat_generation.txt  �o��
			write(36,'(E15.7)') hef_mesh		!heat_generation_all.txt  �o��

			hef_mesh = 0.0		!�S�̔z��w��

			rewind	33
			write(33,'(E15.7)') hef_scat		!heat_generation_scat.txt  �o��

			hef_scat = 0.0		!�S�̔z��w��

c			nxnz = (nx+1)*(nz+1)
c			allocate (hef_fin(nxnz,nscat,nvalley))
c			hef_fin=hef_scat*buff	!hef_fin,hef_scat:�S�̔z��

!SMP$	ASSERT (ITERCNT(70))	!�X�p�R���p�œK���⏕����
c			do k=1,nvalley
c			do j=1,nscat
c				hef_fin( 0,0:nz,j,k)=hef_fin( 0,0:nz,j,k)*2.0
c				hef_fin(nx,0:nz,j,k)=hef_fin(nx,0:nz,j,k)*2.0
c				hef_fin(0:nx, 0,j,k)=hef_fin(0:nx, 0,j,k)*2.0
c				hef_fin(0:nx,nz,j,k)=hef_fin(0:nx,nz,j,k)*2.0
c			enddo
c			enddo

c			hef_scat = 0.0		!�S�̔z��w��

c			rewind	33			!heat_generation_scat.txt �o��
c			write(33,'(E15.7)') hef_fin
c			deallocate (hef_fin)
c
			count = 0
	endif
	endif
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c------���t�l�X�U��--J.R. Watling���f��(2012�N�t����)-------------------------
      if(modulo(ict,jouta).eq.(jouta-1)) then

	    open(1919,file='roughness_data.txt')
          write(1919,*) '�p�x������',split
          write(1919,*) '���ʍ���',delta
		write(1919,*) '���ʑ��֒�',lambda
	    write(1919,*) '���ː�',count_reflection
          write(1919,*) '�U����',count_roughness
          write(1919,*) '���ϊp�x�ω���',average

          open(0721,file='reflection_count_data.txt')
	    do vvv=0,nx
	        write(0721,*) (basho_reflection(vvv,yyy),yyy=1,2)
	    end do

	    open(6969,file='roughness_count_data.txt')
	    do xxx=0,nx
	        write(6969,*) (basho_roughness(xxx,yyy),yyy=1,2)
          end do

	end if

c-----	�̈�ʉ߂��闱�q���J�E���g120817sato
c-----pass(1)A�O����̈�ʉ�,(2)A�O����̈��,(3)�̈������B�O,(4)�̈������A�O,
c-----pass(5)�̈���A��,(6)�̈�A����B�O,(7)�̈�A����A�O	 (10)�̈���A��(��)
c-----pass(8)B�O����̈�ʉ�,(9)B�O����̈��

      if(modulo(ict,jouta).eq.(jouta-1)) then

		x_mean_free_path_sum = x_mean_free_path_sum 
     &							/ x_mean_free_path_count
		mean_free_path_sum = mean_free_path_sum 
     &							/ x_mean_free_path_count

		write(852,*) 'A�O����̈�ʉ�',',',pass(1)
		write(852,*) 'A�O����̈��',',',pass(2)
		write(852,*) '�̈������B�O',',',pass(3)
		write(852,*) '�̈������A�O',',',pass(4)
		write(852,*) '�̈���A��',',',pass(5)
		write(852,*) '�̈�A����B�O',',',pass(6)
		write(852,*) '�̈�A����A�O',',',pass(7)
		write(852,*) 'B�O����̈�ʉ�',',',pass(8)
		write(852,*) 'B�O����̈��',',',pass(9)
		write(852,*) '���ώ��R�s��x',',',x_mean_free_path_sum
		write(852,*) '���ώ��R�s��',',',mean_free_path_sum
		write(852,*) ''

	    do xxx=0,nx
			mean_free_path_sum2(xxx,2)
     &			=mean_free_path_sum2(xxx,2)/mean_free_path_sum2(xxx,1)
			mean_free_path_sum2(xxx,3)
     &			=mean_free_path_sum2(xxx,3)/mean_free_path_sum2(xxx,1)
			
			write(853,*)xxx,',',mean_free_path_sum2(xxx,1),','
     &					,mean_free_path_sum2(xxx,2),','
     &					,mean_free_path_sum2(xxx,3)
		enddo
		write(853,*)''
		
		x_mean_free_path_sum = 0.0	!������
		mean_free_path_sum = 0.0	!������
		x_mean_free_path_count = 0
		mean_free_path_sum2 = 0.0	!������
	endif

c------���t�l�X�U��--R.P. Joshi���f��(2012�N�H)-------------------------
      if(modulo(ict,jouta).eq.(jouta-1))then
		do ix=0,nx
		write(1210,*)ix,',',roughness1_countx(ix,1:narea,1)		!�`���l����w�e���E��x�������t�l�X��
		write(1211,*)ix,',',roughness1_countx(ix,1:narea,2)		!�`���l�����w�e���E��x�������t�l�X��
		enddo
		do ie=1,nemax
		write(1212,*)ix,',',roughness1_counte(ie,1:narea,1)		!�`���l����w�e���E�ʃG�l���M�[���t�l�X��
		write(1213,*)ix,',',roughness1_counte(ie,1:narea,2)		!�`���l�����w�e���E�ʃG�l���M�[���t�l�X��
		enddo
		 
		write(1210,*)''
		write(1211,*)''
		write(1212,*)''
		write(1213,*)''
		
		roughness1_countx = 0	!������
		roughness1_counte = 0	!������
	endif

c ----( cpu���Ԕz���o�� )----
	if(sw(8).gt.0)then
	if(modulo(ict,jcput).eq.(jcput-1))then
		rewind	35
		cpusum = sum(cput)
		write(35,"('emcd   = ',f7.3,'%')") cput(1)/cpusum*100
		write(35,"('renew  = ',f7.3,'%')") cput(2)/cpusum*100
		write(35,"('charge = ',f7.3,'%')") cput(3)/cpusum*100
		write(35,"('poisso = ',f7.3,'%')") cput(4)/cpusum*100
		write(35,"('output = ',f7.3,'%')") cput(5)/cpusum*100
		write(35,"('other  = ',f7.3,'%')") cput(6)/cpusum*100
		write(35,"('cpusum = ',f15.1,'s')") cpusum
	endif
	endif

	end subroutine output