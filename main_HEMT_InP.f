c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c InP-HEMT��EMC�V�~�����[�V����
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	implicit none
c	�ϐ���`���[���i�����^�j
c	n... -> �z���`�p�ϐ��E�ő�l(��`��͒萔)
c	m... -> mesh�ʓ���Ԋi�[�ϐ�(��`��͒萔)
c	l... -> �i�q���W�p�i�[�ϐ�(��`��͒萔)
c	j... -> ��؂莞�ԁE���[�v�p�ő�l(�قڒ萔)
c	i... -> ���[�v�p�C���N�������g�ϐ�
c	k... -> �p�����[�^�i�[�ϐ�
c
	real qe
	parameter(qe   = 1.60219e-19)
c
c	--- �z���`�p���� ---
	include 'arraysize.fi'
c############circuit 1030 takahashi##############################################
	include 'Circuit.fi'
c
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c	--- �V�~�����[�V�������� ---
	real	dt					!dt:�P�ʎ���[s]
	real	spnum				!spnum:�����q��	epp
	integer	istp,ict			!�X�e�b�v
	real,	dimension (:),allocatable	:: de	!�G�l���M�[�����P��	  [eV]
	integer np1							!�P���b�V��������̗��q�� [��]
	integer	jpot
c	--- ���q��� ---
	integer	jpnum	!�����q��
	real,	dimension (:,:),allocatable	:: p	!���q�p�����[�^(kx,ky,kz,��,x,z)
	integer(1),dimension (:,:),allocatable	:: kp	!�JNo.,�G�l���M�[No,�G���ANo
c	--- �s�����Z�x�ɂ��̈敪�� ---
c	dx,dz:���b�V�����Axmax,zmax:�f�o�C�X��[m]
c	cxpart1,czpart1,cxpart2,czpart2		�G���A�̍��W(coordinate)[m]
c	lxpart1,lzpart1,lxpart2,lzpart2		�G���A�̃��b�V��No.(lattice)
	real	dx,dz,xmax,zmax		
	real,	dimension (:),allocatable	:: cxpart1,czpart1,cxpart2,czpart2
	integer(2),dimension (:),allocatable	:: lxpart1,lzpart1,lxpart2,lzpart2
	real,	dimension (:),allocatable	:: dconc		!doping concentration[m-3]
	real,dimension(:),allocatable::	dn3					!09/2/19 �|��

c
c	---���Z�X�̈�---
c	cxrecess1,czrecess1,cxrecess2,czrecess2		���Z�X�G���A�̍��W(coordinate)[m]
c	lxrecess1,lzrecess1,lxrecess2,lzrecess2		�G���A�̃��b�V��No.(lattice)
	real, dimension (:),allocatable :: cxrecess1,czrecess1
	real, dimension (:),allocatable :: cxrecess2,czrecess2
	integer(2),dimension (:),allocatable :: lxrecess1,lzrecess1
	integer(2),dimension (:),allocatable :: lxrecess2,lzrecess2
c
c	---�Փ˓d��---
c	nava:���b�V���ʏՓ˓d���p�x	
	integer,allocatable :: nava(:,:,:),cnava(:,:)
	integer,allocatable :: n_scat(:,:,:),count_scat(:)
c
c	--- �ޗ��ʃp�����[�^ ---
c	!twodeg:2DEG�V�[�g�d�ז��x[m^-3],eps:�U�d����s,iarea:�f�o�C�X�ޗ�
	real,	dimension (:),allocatable	:: twodeg,eps
	integer(1),dimension(:),allocatable :: iarea	
	integer(2),dimension (:),allocatable	:: lhet		!�w�e���ڍ��̈ʒu
c	--- �d�ɔz�u�Ɠd�ׂ̈ړ� ---
	real	lg					!�Q�[�g��[m-3]
	integer	ncon				!�d�Ƀ��b�V���̂���ׂ����q�� [��]
	real,	dimension	(:),allocatable	:: cxpole1,cxpole2	!�e�d�ɂ̍��W[m]
	integer(2),dimension(:),allocatable	:: lnpole1,lnpole2	!�e�d�ɂ̍��W[mesh]
	integer(1),dimension(:),allocatable	:: melpos 		!�e�d�Ɉʒu(1:��A2:���A3:�E)intype
	real,	dimension	(:),allocatable	:: vb,vi		!vb:�r���h�C���d��,vi:�d�ɓd��
	real,	dimension	(:),allocatable	:: bias,dvi		!�d�ɓd������p�ϐ�
	integer(2),dimension(:),allocatable	:: idvi			!vi=bias+(dvi*idvi)
	real,	dimension	(:),allocatable	:: cur			!�d�ɕʓd���l
	integer(4),dimension(:),allocatable	:: jspsum		!�e�d�ɂɏ����������q�� jqn
	real	c_ratio
c	--- �o���h�\���p�����[�^ ---
	real,	dimension (:,:),allocatable	:: smh,hhm,hm,af,af2,af4	!param�Q��
	real,	dimension (:,:),allocatable :: eg		!�G�l���M�[�M���b�v
	real,	dimension (:,:),allocatable	:: ec		!�R���_�N�V�����o���h�G�l���M�[
	real, dimension	(nvalley,narea):: dltec
c	--- �U�����[�g ---
	real,	dimension (:,:,:,:,:,:),allocatable	:: swk		!�K�i���U�����[�g
	real,	dimension (:,:,:),allocatable		:: pgm		!�U�����[�g�ő�l�̋t��
	real,	dimension (:,:,:),allocatable		:: escat	!�U���G�l���M�[
	real,	dimension (:,:),allocatable			:: am			!120201 ���ω������L������
	integer(1),dimension (:,:),allocatable		:: iarg		!�U�����
	integer(1),dimension (:,:,:),allocatable	:: iband	!�U����J�ڒJ
c	real	dn1				!�s�����Z�x
c	--- ���b�V��������p�����[�^ ---
	real,	dimension (:),allocatable	:: u		!�|�e���V����
	real,   dimension (:),allocatable	:: epA,epB,epC,epA2,epB2				!120126homma
	real,	dimension (:),allocatable	:: cn		!�d�ז��x(�d�q)
	real,	dimension (:,:),allocatable	:: cp		!�d�ז��x(�z�[��)
	real,	dimension (:),allocatable	:: hef_mesh	!���M�ʕ��z
	real,	dimension (:),allocatable	:: dopem	!mesh�ʕs�����Z�x dtype
	integer(1),dimension (:),allocatable	:: marea	!mesh�ʃG���ANo. iatype
c	--- ���M�֌W�p�����[�^ ---
	real	btmp,dtmp	!�U�����[�g�v�Z���x�Abtmp:�x�[�X���x�Adtmp:���s
	integer(2),dimension(:),allocatable		:: mtemp	!���b�V���ʉ��x(T=btmp+mtemp*dtmp)
	real,	dimension (:),allocatable		:: bktq	![Kb*T*q],Kb:�{���c�}���萔
	real,	dimension (:,:,:),allocatable	:: hescat	!�U�������o�G�l���M�[
	real,	dimension (:,:,:),allocatable	:: hef_scat	!mesh�ʔ��M��
	integer hcss	!�M�v�Z�t���O�A0:�v�Z���Ȃ��A1:�v�Z����
	integer	count
c	--- ���o�͗p ---
	integer sstat
c	--- �Z���d�׉_ ---
	real,allocatable :: cloud(:)	!�񎟌��K�E�X���z
	integer	ibord						!�K�E�X���z������[��
c
c-----�k�ތ���-----
	real,	dimension (:,:,:),allocatable	:: avsumtel
	real,	dimension (:,:,:),allocatable	:: avsumconc
	real,   dimension (:,:,:),allocatable   :: avsumconc_1  !101221�d�q�Z�x�`���l������
	real,   dimension (:),allocatable       :: avsumconcA   !101221�d�q�Z�x�`���l������2
	real,	dimension (:,:,:),allocatable	:: efermi
	
c
	real,	dimension (:,:),allocatable	:: avsumconc1
	real(8),	dimension (:,:),allocatable	:: avsumtei11
c
	real,	dimension(:),allocatable :: E2
	real,	dimension(:),allocatable :: n2
	real,	dimension(:),allocatable :: E3
	real,	dimension(:),allocatable :: n3
c
	real ntab1(0:300000,5)		!�Z�x(�ԍ��A�ޗ�)
	real etab1(0:300000,5)		!�G�l���M�[(�ԍ��A�ޗ�)
c
	real de2,e2max,de3,e3max
	integer i2max,i3max
c
	integer,allocatable :: n_scat_p(:,:,:)	!(+)yama071223
	integer,allocatable :: n_scat_n(:,:,:)	!(+)yama071223
c------------------
	real,	dimension(:,:),allocatable :: ecr		!121019sato
c-----Rejection�̔����������-----
	real(8) aff(nvalley,narea)
c
c----------�o���X�e�B�b�N�̌v�Z---------------
	integer(1),dimension (0:npmax) ::	balis_flag			!07/11/22
c
c	--- poisso�p�ϐ� ---
	real	maceps				!�v�Z�@�C�v�V����
c----�Փ˓d��-----
	real,	dimension (:,:,:),allocatable		:: am_aniso						!20100624
	real,	dimension (:,:,:),allocatable		:: aff_aniso					!20100624
	real,	dimension (:,:,:),allocatable		:: hole_am_aniso				!20100624
	real,	dimension (:,:,:),allocatable		:: hole_aff_aniso				!20100624
	real,	dimension	(:,:),allocatable		:: hiXL
	real,	dimension	(:,:,:),allocatable		:: Eth_table					!20100624
c	real(4) highestX,highestL
	real, dimension(:),allocatable ::  II_S	!120921sato

c------���t�l�X�U��--J.R. Watling���f��(2012�N�t����)-------------------------
	integer,dimension (:,:),allocatable :: basho_roughness
	integer,dimension (:,:),allocatable :: basho_reflection
c-----	�̈�ʉ߂��闱�q���J�E���g120817sato
	real(8),dimension (10) :: pass,pass_r		!�J�E���^�Cr�̓��W�F�N�V�����΍�
	real,	dimension (:),allocatable	:: x_start
	real	xx_goal
	real,	dimension (:),allocatable	:: z_start
	real	zz_goal
	integer	scatpoint			!�U���D��U������
	real	x_mean_free_path_sum
	real	mean_free_path_sum
	integer x_mean_free_path_count
	real,dimension (:,:),allocatable :: mean_free_path_sum2	!���ώ��R�s��x���W	 1:�J�E���^,2:x�������ώ��R�s��,3:���ώ��R�s��
c------���t�l�X�U��--R.P. Joshi���f��(2012�N�H)-------------------------
	integer ia,iv,ie,ien
	real,	dimension (:,:,:,:),allocatable	::swk_rou	
	integer,dimension (:,:,:),allocatable :: roughness1_countx
	integer,dimension (:,:,:),allocatable :: roughness1_counte

c
c	--- ���[�v���� ---
	common /step/jinit,jstat,jstep,jdisp
	integer	jinit,jstat,jstep,jdisp
	common /outp/joutc,jouti,jouta,jsave,jheat,jcput,joutpi,sw
	integer	joutc,jouti,jouta,jsave,jheat,jcput,joutpi,sw(8)
	integer	nig,nic
c	integer iii
c==========================
c############circuit 141030 takahashi##############################################
c�@MC�Ɖ�H�V�~�����[�V�����̃C���^�t�F�[�X��VG1,VD1�AVS1�����IG1,�c
c##################################################################################
	integer icount,iic,iic2,i_out
	real pi,VDC,VAC,f
	parameter(pi  = 3.141592)
	integer iiostp2,iiostp	!�P�����̊Ԃɉ��_�o�͂��邩
	integer jc_on	!1:�}�b�`���O��H,2:�񐶕����̉�H, 3:Source�̂݊񐶕����̉�H,4:OFF

	real VG1(-1:jtp00),VD1(-1:jtp00),VS1(-1:jtp00),VO(-1:jtp00)				!��H�o�́i�����͂̓d�ɂւ̈���d���j
	real V1(-1:jtp00),IG1(-1:jtp00),IDS1(-1:jtp00),ID1(-1:jtp00)
	real ISS1(-1:jtp00),IO(-1:jtp00)								!CT�Ԃ̕��ϓd��(MC�o��)
	real Wg,tk
	integer dtc,vave,iave,vave00						
	real CT												!��H�̃X�e�b�v���ԕ�
	real t
	real ID1_ave,IG1_ave,IS1_ave						
	integer icstp
	real IG1_stack(-1:jtp00),IDS1_stack(-1:jtp00),ISS1_stack(-1:jtp00)	!MC�o�͊i�[�p�z��
	integer i_or_c	
	integer	cflag						!��H�p�T�u���[�`�����ʃt���O				
c######### 1 #############################################################
c---Ex�p-------------------
	real V_in,V_out,Vd_add,Vg_add,tdcount
	integer idcount
c==========================
c%%%%%%%%%���[�J���ϐ�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	integer	npmaxr
	real	xmaxr,zmaxr,spnumr
c	--- ���[�v���� ---
	integer i1,i2
	integer	ni1,ni2
	integer	k1,k2
	integer	icn
c	--- �o�͗p�ϐ� ---
	character(80) form,form_d,form_g
	integer date_time(8)
	real	t1,t2
c
	integer	nxnz,mx,i
	integer	ix,iz
	integer	jx,jz
c
	real cpp
c
	CALL SYSTEM_CLOCK(COUNT=sstat)
c
	call init_rand
c
c---( �t�@�C���̃I�[�v�� )---
c	open(unit=0,)	!�W���G���[�o�́i��ʁj
	open(unit=3,file='disp.txt')
	open(unit=4,file='data.txt')
c	open(unit=5,)	!������ԂŃL�[�{�[�h����
c	open(unit=6,)	!������Ԃŉ�ʏo��
	open(unit=7,file='current.txt')
c	open(unit=8,file='read')
	open(unit=9,file='state.txt')
	open(unit=10,file='outdata.txt')
	open(unit=11,file='jspsum.txt')
	open(unit=13,file='jpnum.txt')
c	open(unit=15,file='swkpara.txt')
c	open(unit=16,file='swk.txt')

	open(unit=20,file='entable.txt')
	open(unit=21,file='enmesh.txt')
	open(unit=22,file='vnmesh.txt')
	open(unit=23,file='kvmesh.txt')
	open(unit=24,file='dnmesh.txt')
	open(unit=25,file='potential.txt')
	open(unit=26,file='density.txt')
	open(unit=27,file='potential2d.txt')
	open(unit=28,file='density2d.txt')

	open(unit=30,file='heat_generation.txt')
	open(unit=31,file='current_d.txt')
	open(unit=32,file='current_s.txt')
	open(unit=33,file='heat_generation_scat.txt')
	open(unit=35,file='cput.txt')
	open(unit=36,file='heat_generation_all.txt')
	open(unit=37,file='density_ave.txt')
	open(unit=38,file='potential_ave.txt')
	open(unit=39,file='cn_ave.txt')
	open(unit=40,file='pinit.bin',form='UNFORMATTED')
	open(unit=41,file='pdata.bin',form='UNFORMATTED')
	open(unit=45,file='current_t.txt')
	open(unit=46,file='kv_zave.txt')
	open(unit=47,file='en_zave.txt')
	open(unit=48,file='vn_zave.txt')
	open(unit=49,file='dn_zave.txt')

	open(unit=56,file='ckv_zave.txt')
	open(unit=57,file='cen_zave.txt')
	open(unit=58,file='cvn_zave.txt')
	open(unit=59,file='cdn_zave.txt')
	open(unit=60,file='ava_zave.txt')

	open(unit=61,file='nava.txt')
	open(unit=62,file='cnava.txt')
	open(unit=63,file='vnava.txt')
	open(unit=64,file='cnscat.txt')

	open(unit=65,file='efield.txt')
	open(unit=66,file='avcn2.txt')

c	open(unit=37,file='potential.bin',form='UNFORMATTED')
c	open(unit=38,file='density.bin',form='UNFORMATTED')

	open(unit=68,file='enx.txt')
	open(unit=69,file='cenx.txt')

	open(unit=70,file='density_bef.txt')
	open(unit=72,file='eff_potential2.txt')
	open(unit=75,file='hole.txt')
	open(unit=76,file='eff_epB2.txt.txt')		!120126homma
	open(unit=77,file='eff_epA.txt')			!120126homma
	open(unit=78,file='eff_epB.txt')			!120126homma

	open(unit=88,file='enx2.txt')
	open(unit=89,file='cenx2.txt')

	open(unit=99,file='error.txt')

c
!	noise�p�Ɏg���t�@�C��(100�`)   �ǉ�    Hara 2005/09/30
	open(unit=100,file='vc.txt')
	open(unit=101,file='vcvel.txt')
	open(unit=102,file='vcdis.txt')			! ���b�V�����̑��x�̂ӂ��
	open(unit=103,file='vcdis(2).txt')
	open(unit=104,file='part_num.txt')
c
	open(unit=108,file='cvn_zave2.txt')
	open(unit=109,file='cdn_zave2.txt')
	open(unit=110,file='balistic.txt')
	open(unit=111,file='kx_distribution.txt')
	open(unit=112,file='ccs_all.txt')
	open(unit=113,file='ccs_scat.txt')
	open(unit=114,file='all_back_scat.txt')
	open(unit=115,file='allback_scat.txt')
	open(unit=116,file='allback_scat2.txt')									
c
c-----�k�ތ��ʊ֘A�̏o�̓t�@�C��-----
	open(unit=120,file='av_enct10.txt')
	open(unit=121,file='av_enct50.txt')
	open(unit=122,file='av_enct210.txt')
	open(unit=123,file='av_enct215.txt')
	open(unit=124,file='av_enct250.txt')
	open(unit=125,file='av_enct285.txt')
	open(unit=126,file='av_enct290.txt')
	open(unit=127,file='av_enct400.txt')
	open(unit=128,file='av_enct490.txt')
c	open(unit=129,file='av_enct200.txt')
c
	open(unit=130,file='avef1.txt')	
	open(unit=131,file='avef2.txt')
	open(unit=132,file='avef3.txt')	
c
	open(unit=133,file='avtel1.txt')
	open(unit=134,file='avtel2.txt')
	open(unit=135,file='avtel3.txt')
c
	open(unit=136,file='cnscat_p.txt')
	open(unit=137,file='cnscat_n.txt')
c
	open(unit=138,file='asscnt1.txt')
	open(unit=139,file='asscnt2.txt')
	open(unit=140,file='asscnt3.txt')
	open(unit=141,file='asscnt4.txt')
	open(unit=142,file='sscnt1.txt')
	open(unit=143,file='sscnt2.txt')
	open(unit=144,file='sscnt3.txt')
	open(unit=145,file='sscnt4.txt')
c
	open(unit=150,file='avsumtel1.txt')		!08/8/6 �|��
	open(unit=151,file='avsumconc1.txt')	!08/8/6 �|�� ���J�̓d�q�Z�x(�Z������)
	open(unit=152,file='avsumtei11.txt')	!08/8/6 �|�� ���J�̓d�q�G�l���M�[(�Z������)
	open(unit=153,file='renew_ec.txt')		!120126homma
c
	open(unit=170,file='iz_ene.txt')
	open(unit=171,file='iz23_ene.txt')
	open(unit=172,file='iz1-29_ene.txt')
	open(unit=173,file='iz2_ene.txt')		!120126homma
	open(unit=174,file='iz3_ene.txt')		!120126homma
	open(unit=175,file='iz_ene_L.txt')
	open(unit=176,file='iz2_ene_L.txt')		!120126homma
	open(unit=177,file='iz3_ene_L.txt')		!120126homma
c	open(unit=178,file='iz43_ene_L.txt')
c	open(unit=179,file='iz45_ene_L.txt')
c
	open(unit=180,file='ix50_ene.txt')
	open(unit=181,file='ix210_ene.txt')
	open(unit=182,file='ix250_ene.txt')
	open(unit=183,file='ix290_ene.txt')
	open(unit=184,file='ix450_ene.txt')
	open(unit=185,file='ix50_ene_L.txt')
	open(unit=186,file='ix210_ene_L.txt')
	open(unit=187,file='ix250_ene_L.txt')
	open(unit=188,file='ix290_ene_L.txt')
	open(unit=189,file='ix450_ene_L.txt')
c 
c############circuit 141030##############################################
	open(unit=221,file='power.csv')
      open(unit=222,file='Device.csv')
      open(unit=223,file='in_v.csv')
      open(unit=224,file='in_i.csv')
	open(unit=225,file='out_v.csv')
	open(unit=226,file='V2_Vo.csv')
	open(unit=227,file='I1_IR.csv')
	open(unit=228,file='out_i.csv')
	open(unit=229,file='Ex_v.csv')
	open(unit=230,file='Ex_i.csv')
c-----nari--------------------------------------------------------------
	open(unit=300,file='current_ypara.txt')
	open(unit=401,file='Ex_current_ypara_NoWg.txt')		!txt
	open(unit=402,file='Ex_current_ypara_kakeruWg.txt')  !txt
c---------------------------------------------------------------------
	open(unit=500,file='Ypara.txt')
	open(unit=501,file='Spara.S11.S21.txt')
	open(unit=502,file='Spara.S12.S22.txt')
	open(unit=503,file='Gausian_Wave.txt')
c	open(unit=190,file='pass_ict1.txt')
c	open(unit=191,file='pass_ict2.txt')
c      open(unit=160,file='plsdmon2.csv')
	open(unit=560,file='ryuushi.dat')
c	open(unit=630,file='IIeief.csv')
c------------------------------------
	
	open(unit=852,file='pass.txt')		!120817sato
	open(unit=853,file='mean_free_path.csv')		!120817sato
c
c	noise�p�Ɏg���t�@�C��(900�`)   �ǉ�    Takahashi 2014/12/23
	open(unit=900,file='dis_dv.txt')			!���J��
	open(unit=901,file='dis_v.txt')			!�h���t�g���x�̗h�炬
	open(unit=902,file='dn2.txt')			! ���b�V�����̑��x�̂ӂ��
	open(unit=903,file='di2.txt')			!
	open(unit=904,file='di2(s_t_c).txt')	!
	open(unit=905,file='di2_valley(s_t_c).txt')
	open(unit=906,file='�d�q���x��.txt')
	open(unit=907,file='���S�Ɍ��̒藝.txt')
	open(unit=908,file='�S���U���X�̃m�C�Y��͖@.txt')
	open(unit=909,file='�S���U���X�̃m�C�Y��͖@�i�␳�j.txt')
c====	input	===========================================================
c---( �z��f�[�^�T�C�Y�ǂݏo�� )---
c	open(unit=8,file='arraysize.data')
c	call arraysize
c	close(8)
c
c---( �f�o�C�X��������ǂݏo�� )----
	allocate(bias(npole),dvi(npole),de(nenergy))
	open(unit=8,file='SimCond.data')
	call condition(dt,np1,bias,dvi,nig,nic,de,jpot)
	close(8)
c
c#####circuit 141030##�f�o�C�X��������ǂݏo��#######################
	open(unit=301,file='CirDevCond.data')
	read(301,*) Wg		!�Q�[�g��
	read(301,*) dtc		!��H�X�e�b�v
	read(301,*) vave	!Id,Ig,Is�𕽋ς���X�e�b�v��
	read(301,*) vave00
	read(301,*) jc_on
	close(301)
c
c---( �f�o�C�X�\���ǂݏo�� )---
	allocate(cxpart1(npart),czpart1(npart))
	allocate(cxpart2(npart),czpart2(npart))
	allocate(lxpart1(npart),lzpart1(npart))
	allocate(lxpart2(npart),lzpart2(npart))
c---(���Z�X)---
	allocate(cxrecess1(nrecess),czrecess1(nrecess))
	allocate(cxrecess2(nrecess),czrecess2(nrecess))
	allocate(lxrecess1(nrecess),lzrecess1(nrecess))
	allocate(lxrecess2(nrecess),lzrecess2(nrecess))
c
	allocate(dconc(npart))
	allocate(lhet(nlayer),iarea(nlayer),twodeg(nlayer))
	allocate(cxpole1(npole),cxpole2(npole))
	allocate(lnpole1(npole),lnpole2(npole))
	allocate(melpos(npole))
	allocate(vb(npole))
	open(unit=8,file='DevCond.data')
	call device(np1,dx,dz,
     &				cxpart1,czpart1,cxpart2,czpart2,
     &                lxpart1,lzpart1,lxpart2,lzpart2,
     &				cxrecess1,czrecess1,cxrecess2,czrecess2,
     &                lxrecess1,lzrecess1,lxrecess2,lzrecess2,
     &				dconc,xmax,zmax,
     &				lhet,iarea,twodeg,
     &				cxpole1,cxpole2,lnpole1,lnpole2,vb,
     &				melpos,lg,spnum,ncon)
	close(8)
c
c
c	---( �f�o�C�X�\���̐ݒ� )---
	nxnz = (nx+1)*(nz+1)
	mx=nx+1
c
c	---( �s�����Z�x�̐ݒ� )---
	allocate(dopem(nxnz),marea(nxnz))
	do i=1,npart
		do jz=lzpart1(i),lzpart2(i)
		do jx=lxpart1(i),lxpart2(i)
			dopem(jx+mx*jz+1) = dconc(i)
			marea(jx+mx*jz+1) = i
		enddo
		enddo
	enddo
c
c	---(���Z�X)---	
	do i=1,nrecess
		do jz=lzrecess1(i),lzrecess2(i)-1
		do jx=lxrecess1(i)+1,lxrecess2(i)-1
			dopem(jx+mx*jz+1) = 0
			marea(jx+mx*jz+1) = 0
		enddo
		enddo
	enddo
c
c---(�z�[���Z�x)---
	allocate(cp(0:nx,0:nz))
	cp(0:nx,0:nz) = 0.0
	open(unit=75,file='hole.txt')
	do iz = 0,nz
		do ix = lxpart1(1),lxpart2(1)
			read(75,*,ERR=800,END=800) cpp
			cp(ix,iz)=cpp
		enddo
c
c	---	n+�̈�̃f�[�^��[   ---
		do ix = 0,lxpart1(1)-1
			cp(ix,iz)=cp(lxpart1(1),iz)
		enddo
		do ix = lxpart2(1)+1,nx
			cp(ix,iz)=cp(lxpart2(1),iz)
		enddo
	enddo
	close(75)
	form = "(x,'�z�[���Z�x���z���m�F���܂����B�ǂݍ��݂܂�')"
	write(*,form)
	goto 900
800	write(*,*) 'hole.txt���ǂݍ��߂܂���'
	write(99,*) 'hole.txt���ǂݍ��߂܂���'
900	continue
c
c	open(unit=75,file='hole.txt')		! �z�[���Z�x����
c	do i = 1, nxnz
c		read(75,*,ERR=800,END=800) cp(i)
c	enddo
c	close(75)
c	form = "(x,'�z�[���Z�x���z���m�F���܂����B�ǂݍ��݂܂�')"
c	write(*,form)
c	goto 900
c800	continue
c	write(*,*) 'hole.txt���ǂݍ��߂܂���'
c	write(99,*) 'hole.txt���ǂݍ��߂܂���'
c	cp(1:nxnz) = 0.0
c900	continue
c	cp(1:nxnz) = cp(1:nxnz)*spnum/dx/dz
c
c---( �f�o�C�X���x�����ǂݏo�� )----
	allocate(mtemp(nxnz))
	open(unit=8,file='latt_temp_fin.txt')
	call temperature(lxpart1(1),lxpart2(1),mtemp,btmp,dtmp)
	close(8)
c
c---( ���ʏo�͏����ǂݏo�� )----
	open(unit=8,file='output.data')
	call outinit(dt)
	close(8)
c
c--	�J�n�����\�� ---
	call date_and_time (values = date_time)	!�����擾
	form = "(x,'�J�n����',2x,I2,'��',I2,'��',I2,2(':',I2.2))"
	write(*,form) date_time(2:3),date_time(5:7)
	write(4,form) date_time(2:3),date_time(5:7)	!���t�����\��
	write(3,form) date_time(2:3),date_time(5:7)
c
c---( �����萔, �ޗ��萔����я��p�����[�^ )---
	write(*,*) 'param �v�Z���E�E�E'
c
	allocate(nava(0:nx,0:nz,0:nvalley),cnava(0:nx,0:nvalley))
	allocate(n_scat(0:nx,nscat,nvalley),count_scat(0:nx))
	allocate(bktq(ntenum))
	allocate(smh(nvalley,narea),hhm(nvalley,narea),hm(nvalley,narea))
	allocate(af(nvalley,narea),af2(nvalley,narea),af4(nvalley,narea))
	allocate(eg(nvalley,narea),ec(nvalley,narea))
	allocate(am(nvalley,narea))
	allocate(am_aniso(nvalley,narea,nvalley))									!20100624
	allocate(aff_aniso(nvalley,narea,nvalley))									!20100624
	allocate(hole_am_aniso(nvalley,narea,nvalley))								!20100624
	allocate(hole_aff_aniso(nvalley,narea,nvalley))	
	allocate(hiXL(narea,nvalley))							!20100624
	allocate(Eth_table(4,20000,2))												!20100624
c	allocate(pgm(nvalley,nenergy,narea),escat(nscat,nvalley,narea))
c	allocate(swk(nscat,nemax,nvalley,nenergy,ntenum,narea))
	allocate(pgm(nvalley,nenergy,npart),escat(nscat,nvalley,narea))	!07/8/4 �s�����U��
	allocate(swk(nscat,nemax,nvalley,nenergy,ntenum,npart))	!07/8/4 �s�����U��
	allocate(iarg(nscat,narea),iband(nscat,nvalley,narea))	
	allocate(hescat(nscat,nvalley,narea))
	allocate(eps(narea))
	allocate(dn3(npart))									!09/2/19 �|��
	allocate(II_S(narea))
c
	nava  = 0;	cnava = 0
	n_scat =0;	count_scat=0
	pgm = 0.0
c
	call param2(
     &			  de,dx,dz,btmp,dtmp,dconc,smh,hhm,hm,am_aniso,am,
     &			  aff_aniso,hole_am_aniso,hole_aff_aniso,af,af2,
     &			  af4,eg,ec,eps,bktq, 
     &			  swk,pgm,escat,iband,iarg,
     &			  hescat,dn3,hiXL,Eth_table,aff,dltec,								!09/2/19 �|��
     &			  II_S)				!120201
c
c---( dx,dz,dt���Ó����ǂ������� )---
	call exam(dx,dz,dt,minval(dconc),
     &		  eps(1),bktq(ntenum),hm(1,1),jpot,pgm)
c
c-----�I�[�~�b�N�d�ɂ̔����G�l���M�[�Ɋւ���v�Z-----
c-----In(0.53)Ga(0.47)As-----
	de2=0.002d0			!de[eV]		!�G�l���M�[���ݕ�
	e2max=0.8d0			!Emax[eV]	!�G�l���M�[�ő�l
	i2max=e2max/de2
	allocate(E2(0:i2max))
	allocate(n2(0:i2max))
c
c-----InAs-----
	de3=0.002d0			!de[eV]		!�G�l���M�[���ݕ�
	e3max=0.8d0			!Emax[eV]	!�G�l���M�[�ő�l
	i3max=e3max/de3
	allocate (E3(0:i3max))
	allocate (n3(0:i3max))
	call ohmic_theory(de2,de3,i2max,i3max,E2,E3,n2,n3,aff,am)
c----------------------------------------------------
c
c---( ��{�f�[�^�����o�� )---
	write(4,*) 'dt=',dt
	write(4,*) 'jstep=',jstep
	write(4,*) 'jinit=',jinit
	write(4,*) 'np1=',np1
	write(4,*) 'spnum=',spnum
	write(4,*) 'nx=',nx,' nz=',nz
	write(4,*) 'dx=',dx,' dz=',dz
	write(4,*) 'lg=',lg		!cxpole2(2)-cxpole1(2)
	write(4,*) 'btmp=',btmp,'��T=',dtmp
	write(4,*) 'ntenum=',ntenum
c
c--	�����I�ɃQ�[�g�E�h���C���̏㏸���������� ---
	allocate(idvi(npole))
	form	= "(x,A12,'���ω����܂���B')"
	if(nig.ge.nic)then
		k1	= 3
		k2	= 2
		if(dvi(k2).eq.0.0)then
			write(* ,form)'�Q�[�g�d��'
			write(3 ,form)'�Q�[�g�d��'
			write(99,form)'�Q�[�g�d��'
		endif
		ni1 = nic-1
		ni2 = nig-1
	else
		k1	= 2
		k2	= 3
		if(dvi(k2).eq.0.0)then
			write(* ,form)'�h���C���d��'
			write(3 ,form)'�h���C���d��'
			write(99,form)'�h���C���d��'
		endif
		ni1 = nig-1
		ni2 = nic-1
	endif
c	----------------------------------------------
c
c--	�f�[�^�����v���O�����p�o�� ---
      open(unit=10,file='outdata.txt')
	write(10,*) nx,nz				!���b�V����
	write(10,*) jstep,jinit,jouta	!�e�X�e�b�v��
	write(10,*) nvalley,nscat		!�J���A�U����
	write(10,*) ni2+1,ni1+1			!�������[�v���A�O�����[�v��
	write(10,*) dvi(k1),bias(k1)		!�O�����[�v�x�[�X�d���A�ω��d��
	write(10,*) dvi(k2),bias(k2)		!�������[�v�x�[�X�d���A�ω��d��
	write(10,*) lhet(1),0.7*(eg(1,2)-eg(1,3))
	write(10,*) dt
	write(10,*) btmp,dtmp,ntenum		!btmp:�x�[�X���x�Adtmp:���s
	close(10)
c
c---( �K�E�X���z�ɏ]�����d�׉_�̌v�Z )---
	allocate(cloud((ngx*2+1)*(ngz*2+1)*ngn))
	if(jdisp.eq.0)write(*,*) '�K�E�X�d�׉_�v�Z'
	call setcloud(dx,dz,cloud)
c
c
c	-------------------------------------------------
c=====���q�������z===========================================================
	allocate(cur(npole),jspsum(npole))
	allocate(vi(npole))
	count = 0
	jspsum = 0
	idvi	= 0	
c--	�ȉ��̏����̎��ɂ�t=0�ŃX�e�b�v�d����������	  ---
	if((nig.eq.1).and.(nic.le.2))then		!�X�e�b�v�o�͗p
		idvi(k1)	= 1
		idvi(k2)	= 1
	endif
	hcss = sw(7) !0		!�M�v�Z�t���O�A0:�v�Z���Ȃ��A1:�v�Z����
	allocate(hef_scat(0:nx,nvalley,nscat))
	nxnz = (nx+1)*(nz+1)	!�S���b�V����
	allocate(hef_mesh(nxnz))
	if(hcss.eq.1)then		!�X�e�b�v�o�͗p
		hef_scat = 0.0		!�S�̔z��
		hef_mesh = 0.0		!�S�̔z��
	endif
c
	vi		= bias		!�S�̔z��
c---( ���������̐ݒ� )---
	if(jdisp.eq.0)write(*,*) '���q�����z�u'
	allocate(p(6,npmax),kp(3,npmax))
c ---- ( ���q�̋�ԕ��z�̐ݒ� ) ----
c	goto 100
	read(40,ERR=100,END=100)	npmaxr,xmaxr,zmaxr,spnumr !pinit.bin�͂ǂ�����?
	form = "(x,'���q���z�t�@�C�����m�F���܂����B�ǂݍ��݂܂�')"
	write(*,form)
	write(99,form)
	if(npmax.ne.npmaxr)then
		form= "(x,'npmax���Ⴂ�܂��Bnpmax��',
     &				I7,'��',I7,'�ɂ��Ă�������')"
c		write(* ,form)npmax,npmaxr		!�G���[���o�� 11/03/19 ��
c		write(99,form)npmax,npmaxr		!�G���[���o�� 11/03/19 ��
		write(*,form)
		write(*,*)npmax,npmaxr
		write(99,form)
		write(99,*)npmax,npmaxr

		goto 200
	endif
	if(spnum.ne.spnumr)then
		form= "(x,'spnum���Ⴂ�܂��B')"
c		write(*,form)spnum,spnumr		!�G���[���o�� 11/03/19 ��
c		write(99,form)spnum,spnumr		!�G���[���o�� 11/03/19 ��
		write(*,form)
		write(*,*)spnum,spnumr
		write(99,form)
		write(99,*)spnum,spnumr

		goto 200
	endif
	read(41,ERR=200)	jpnum,p,kp
	goto 400
c
c---	�G���[����	 ---
200	continue
	form = "(x,'���q���z�t�@�C�����ǂݍ��߂܂���ł����B')"
	write(*,form)
	write(99,form)
	stop
100	continue
	call initia(
     &				spnum,dx,dz,xmax,zmax,de,dopem,twodeg,
     &				lnpole1,lnpole2,lhet,iarea,bktq,mtemp,smh,pgm,
     &				p,kp,jpnum,af,n2,E2,i2max,	!120126homma
     &				cxrecess1,czrecess1,cxrecess2,czrecess2,
     &				lxrecess1,lzrecess1,lxrecess2,lzrecess2,
     &				czpart2)	!07/8/4 �s�����U��
400	continue
c
c---(	�d�ɕt�߂ł̗��q�̑��� )---
	if(jdisp.eq.0)write(*,*) 'renew'
c----------�o���X�e�B�b�N�̌v�Z---------------
	allocate(ecr(7,int(nemax/4)))	!121019sato
	balis_flag = 0
	ecr = 0.0		!120126homma
c---------------------------------------------
c----pass������ 120817sato
	allocate(x_start(npmax))
	allocate(z_start(npmax))
	allocate(mean_free_path_sum2(0:nx,3))
	pass = 0.0
	pass_r = 0.0
	x_start = 0.0
	xx_goal = 0.0
	z_start = 0.0
	zz_goal = 0.0
	scatpoint = 0	!������	   =1:�U��,=0:��U��
	x_mean_free_path_sum = 0.0
	x_mean_free_path_count = 0
	mean_free_path_sum = 0.0
	mean_free_path_sum2 = 0.0	
c--------------------------------

	call renew(
     &			jpnum,ncon,
     &			bktq,mtemp,dx,dz,lnpole1,
     &			lnpole2,melpos,jspsum,vb,
     &			pgm,smh,p,kp,lhet,iarea,twodeg,balis_flag,
     &			hhm,af,af4,ecr,de,			!120126homma
     &			i2max,E2,n2,i3max,E3,n3,	!08/8/6 �|��
     &			x_start,z_start)		!120817sato

c
c---(	���b�V��������̗��q�Z�x�̃J�E���g )---
	if(jdisp.eq.0)write(*,*) '1st charge'
	nxnz = (nx+1)*(nz+1)
	allocate(u(nxnz),cn(nxnz))
	allocate(epA(nxnz),epB(nxnz),epC(nxnz),epA2(nxnz),epB2(nxnz))	!120126homma


	!c_ratio = 0.6


	cn = 0.0;icn=0;u=0.0
	epA=0.0;epB=0.0;epC=0.0;epA2=0.0;epB2=0.0	!120126homma
	ibord= lhet(nlayer-1)+2
c		ibord= lhet(nlayer-1)+2 !11/04/08�C���K�v�ӏ��H
	call charge(jpnum,dx,dz,spnum,p,kp,cn,icn,1,ibord,cloud,
     &			lxrecess1,lzrecess1,lxrecess2,lzrecess2)
c
	rewind	26;	write(26,'(e15.7)') cn	!density.txt �o��
c	rewind	38;	write(38) 2,nx,nz,cn		!density.bin �o��
c
c-----�k�ތ���ON/OFF�̕\��-----
	if(rejcnt.eq.0)then
		write(*,*)'�p�E���̔r����OFF�Ōv�Z���܂�'
		write(3,*)'�p�E���̔r����OFF�Ōv�Z���܂�'
	elseif(rejcnt.gt.0)then
		write(*,*)'�p�E���̔r����ON Rejection������',rejcnt,'��Ōv�Z���܂�'
		write(3,*)'�p�E���̔r����ON Rejection������',rejcnt,'��Ōv�Z���܂�'
	else
		write(*,*)'rejcnt���s���ł�'
		write(3,*)'rejcnt���s���ł�'
		stop
	endif
c------------------------------
c---(	���q�Z�x����|�e���V�����̓��o(�|�A�\��������) )---
!	if(jdisp.eq.0)
	write(*,*) '1st poisso'
	maceps = epsilon(maceps) !�v�Z�@�C�v�V����
	maceps = 1.0E-12 !7		!�v�Z�@�C�v�V����
	call poisso(dx, dz, lhet, twodeg,iarea,
     &			vb,vi,lnpole1,lnpole2,melpos,
     &			eps,eg,c_ratio,
     &			u,cn,cp,dopem, maceps,
     &			lxrecess1,lzrecess1,lxrecess2,lzrecess2,dltec)
c
c
!	if(jdisp.eq.0)
	write(*,*) 'output u&cn'
	rewind	25;	write(25,'(f12.7)') u	!potential.txt �o��
c	rewind	37;	write(37) 2,nx,nz,u 		!potential.bin �o��
c
c-----�k�ތ���-----
	ntab1=0.0
	etab1=0.0
	call fermi_calc_table(ntab1,etab1,aff,am)
c
	allocate(efermi(0:nx,0:nz,nvalley))
	allocate(avsumtel(0:nx,0:nz,nvalley))
	allocate(avsumconc(0:nx,0:nz,nvalley))
      allocate(avsumconc_1(0:nx,0:nz,nvalley)) !101221 �d�q�Z�x�`���l������
	allocate(avsumconcA(0:nx))               !101221 �d�q�Z�x�`���l������
c
	allocate(avsumconc1(0:nx,0:nz))
	allocate(avsumtei11(0:nx,0:nz))
c
	efermi=0.0
	avsumtel=0.0
	avsumconc=0.0
	avsumconc_1=0.0
	avsumconcA=0.0
c
	avsumconc1=0.0
	avsumtei11=0.0
c
	allocate(n_scat_p(0:nx,nscat,nvalley),n_scat_n(0:nx,nscat,nvalley))
	n_scat_p =0
	n_scat_n =0
c------------------
c------���t�l�X�U��--J.R. Watling���f��(2012�N�t����)-------------------------
	allocate(basho_roughness(0:nx,2))
	allocate(basho_reflection(0:nx,2))
	basho_roughness = 0
	basho_reflection = 0
c
c------���t�l�X�U��--R.P. Joshi���f��(2012�N�H)-------------------------
	allocate(swk_rou(nemax,nvalley,nenergy,narea))		!121029sato
	allocate(roughness1_countx(0:nx,narea,2))
	allocate(roughness1_counte(nemax,narea,2))
	open(unit=1209,file='swk_roughness.txt')	!���t�l�X�U�����[�g
c	���t�l�X�U���J�E���^
	open(unit=1210,file='roughness_x1.csv')	!�`���l����w�e���E��x�������t�l�X��
	open(unit=1211,file='roughness_x2.csv')	!�`���l�����w�e���E��x�������t�l�X��
	open(unit=1212,file='roughness_e1.csv')	!�`���l����w�e���E�ʃG�l���M�[���t�l�X��
	open(unit=1213,file='roughness_e2.csv')	!�`���l�����w�e���E�ʃG�l���M�[���t�l�X��
c	open(unit=1245,file='swk_roughness_main.csv')	!���t�l�X�U�����[�g(�m�F�p)

	do ia = 1, narea
	do iv = 1, nvalley
	do ie=1,nemax
	do ien=1,nenergy
		read(1209,*)swk_rou(ie,iv,ien,ia)		!���t�l�X�U�����[�g�K�i����_R.P.Joshi	
	end do	!ien
	end do	!ie
	end do	!iv
	end do !ia


	roughness1_countx=0		!������
	roughness1_counte=0		!������
c	
c	iv = 1
c	do ie=1,nemax
c	ien=1
c		write(1245,*)swk_rou(ie,iv,ien,1:2)		!���t�l�X�U�����[�g�K�i���o��_R.P.Joshi	
c	!ien
c	end do	!ie
c	!iv
c	!ia
c	write(*,*)'���t�l�X�U�����[�g�ǂݎ��'


c====================================================================================
c---( �����q�����e�J�����v�Z )----
c
c	(�o�̓t�H�[�}�b�g�ݒ�)
	form_d	= "(x,'�h���C���d��',F8.3,'V ',F7.1,'ps-> ',F7.1,'ps')"
	form_g	= "(x,'�Q�[�g�d��  ',F8.3,'V ',F7.1,'ps-> ',F7.1,'ps')"
c
c
	write(*,*) 'jpnum=',jpnum
	write(3,*) 'jpnum=',jpnum
	write(4,*) 'jpnum=',jpnum
c
c
	ict = -jinit
c	write(*,*) -jinit,jstat-1
c######circuit 100110(IKEDA)###### �������� ###########################################
      if(jc_on.eq.1)then
	open(unit=304,file='InitiaCir.data')
		read(304,*) VG1(-1)		
		read(304,*) VD1(-1)		
		read(304,*) VS1(-1)
		read(304,*) vi(1)	!Is	
		read(304,*) vi(2)	!IG1	
		read(304,*) vi(3)	!Id
	close(304)
	endif

      if(jc_on.eq.2)then
	open(unit=304,file='ExInitiaCir.data')
		read(304,*) VG1(-1)		
		read(304,*) VD1(-1)		
		read(304,*) VS1(-1)
		read(304,*) vi(1)	!Vs	
		read(304,*) vi(2)	!VG	
		read(304,*) vi(3)	!Vd
		read(304,*) V_in	!Vd2�̏����l
		read(304,*) V_out	!Vd2�̏����l
		read(304,*) Vd_add	!Vd2�̑�����
		read(304,*) Vg_add	!Vg2�̑�����
		read(304,*) tdcount	!Ex��H�̈��Ԃ�����̎���
	close(304)
	endif

      if(jc_on.eq.3)then
	open(unit=304,file='S_ExInitiaCir.data')
		read(304,*) VG1(-1)		
		read(304,*) VD1(-1)		
		read(304,*) VS1(-1)
		read(304,*) vi(1)	!Vs	
		read(304,*) vi(2)	!VG	
		read(304,*) vi(3)	!Vd
		read(304,*) V_in	!Vd1�̏����lnnn
		read(304,*) V_out	!Vg1�̏����l	
		read(304,*) Vd_add	!Vd1�̑�����
		read(304,*) Vg_add	!Vg1�̑�����
		read(304,*) tdcount	!S_Ex��H�̈��Ԃ�����̎���
	close(304)
	endif
c********(�i�s�g�Z�o)**************************************************
      if(jc_on.eq.4)then
	open(unit=304,file='ExInitiaCir.data')
		read(304,*) VG1(0)		
		read(304,*) VD1(0)		
		read(304,*) VS1(0)
		read(304,*) vi(1)	!Vs	
		read(304,*) vi(2)	!VG	
		read(304,*) vi(3)	!Vd
		read(304,*) V_in	!Vd2�̏����l
		read(304,*) V_out	!Vd2�̏����l
		read(304,*) Vd_add	!Vd2�̑�����
		read(304,*) Vg_add	!Vg2�̑�����
		read(304,*) tdcount	!Ex��H�̈��Ԃ�����̎���
	close(304)
	endif
c***********************************************************************
	idcount=anint(tdcount/dt)
c###############################################################
c
c---( ����Ԃ֎��������� )---
	if(-jinit.le.jstat-1)then
		write(*,*) '�������E�E�E'
		t1 = dt*float(-jinit);	t2 = dt*float(jstat-1)
		t1 = t1 * 1.0e12;		t2 = t2 * 1.0e12	!sec -> ps
		write(*,form_d) vi(3)
		write(4,form_d) vi(3)
		write(3,form_d) vi(3)
		write(*,form_g) vi(2),t1,t2
		write(4,form_g) vi(2),t1,t2
		write(3,form_g) vi(2),t1,t2
		istp = 0
c
c###################circuit 2014/12/01(takahashi)### �������� ##############
		i_or_c=0
		do iave=-1,jtp00
			IDS1_stack(iave)=0.0		
			IG1_stack(iave)=0.0
			ISS1_stack(iave)=0.0
		enddo	
c###################circuit 2014/12/01(takahashi)###########################
		do ict=-jinit,jstat-1
c	�`�F�b�N�p
c		write(190,*) 'ict',ict,'jpnum',jpnum
			call sim_dt2(
     &              am_aniso,aff_aniso,hole_am_aniso,hole_aff_aniso,
     &			  dt,spnum,istp,ict,jpnum,
     &              dx,dz,xmax,zmax,lhet,iarea,dopem,twodeg,c_ratio,
     &              cxpart1,cxpart2,lxpart1,lxpart2,
     &              vb,vi,cur,cxpole1,cxpole2,lnpole1,lnpole2,
     &              melpos,jspsum,ncon,
     &              smh,hhm,hm,af,af2,af4,eps,eg,ec,bktq,
     &              de,swk,pgm,escat,iarg,iband,
     &              p,kp,u,cn,cp,cloud,hef_mesh,
     &              maceps,hescat,mtemp,
     &		      hef_scat,hcss,count,jpot,sstat,
     &			  cxrecess1,czrecess1,cxrecess2,czrecess2,
     &			  lxrecess1,lzrecess1,lxrecess2,lzrecess2,
c     &			  nava,cnava,n_scat,count_scat)
c     &			  nava,cnava,n_scat,count_scat,czpart2)
     &			  nava,cnava,n_scat,count_scat,czpart2,balis_flag,		!07/8/4 �s�����U��
     &			  n_scat_p,n_scat_n,									!08/8/6 �|��
     &			  efermi,avsumtel,avsumconc,i2max,E2,n2,i3max,E3,n3,	!08/8/6 �|��
     &		      avsumconc1,avsumtei11,ntab1,etab1,dn3,hiXL,
     &              epA,epB,epC,epA2,epB2,ecr,					!120126homma
     &              avsumconc_1,avsumconcA,Eth_table,			!09/2/19 �|��
     &			  am,aff,II_S,
     &			  basho_reflection,basho_roughness,	
     &			  pass,pass_r,x_start,xx_goal,scatpoint,
     &			  x_mean_free_path_sum,x_mean_free_path_count,			!120817sato
     &			  z_start,zz_goal,mean_free_path_sum,
     &			  mean_free_path_sum2,swk_rou,							!121029sato	
     &			  roughness1_countx,roughness1_counte,
c############circuit 100110(IKEDA)##############################################
     &			  IDS1_stack,IG1_stack,ISS1_stack,i_or_c,jc_on)	
c
		istp = istp+1
		enddo
	endif
c
c
	rewind	40;	write(40)	npmax,xmax,zmax,spnum
	rewind	41;	write(41)	jpnum,p,kp
c
	hcss = sw(7)	!�M�v�Z�t���O�A0:�v�Z���Ȃ��A1:�v�Z����
c
c
	do i1=0,ni1
		vi(k1) = bias(k1) + dvi(k1)*idvi(k1)	! �h���C��(�Q�[�g)�d�����Đݒ�
		if(k1.eq.3)then
			write(*,form_d) vi(3)	!�h���C���d���E��ʏo��
			write(4,form_d) vi(3)	!�h���C���d���E�t�@�C���o��
			write(3,form_d) vi(3)	!�h���C���d���E�t�@�C���o��
		else
			write(*,form_g) vi(2)	!�Q�[�g�d���E��ʏo��
			write(4,form_g) vi(2)	!�Q�[�g�d���E�t�@�C���o��
			write(3,form_g) vi(2)	!�Q�[�g�d���E�t�@�C���o��
		endif
	do i2=0,ni2
		vi(k2) = bias(k2) + dvi(k2)*idvi(k2)	! �Q�[�g(�h���C��)�d�����Đݒ�
		jspsum = 0								!jspsum:�d�ɂɓ��������q��
		t1 = dt*float(ict);	t2 = dt*float(ict+jstep-1)
		t1 = t1 * 1.0e12;	t2 = t2 * 1.0e12			!sec -> ps
		if(k2.eq.2)then
			write(*,form_g) vi(2),t1,t2	!�Q�[�g�d���E��ʏo��
			write(4,form_g) vi(2),t1,t2	!�Q�[�g�d���E�t�@�C���o��
			write(3,form_g) vi(2),t1,t2	!�Q�[�g�d���E�t�@�C���o��
		else
			write(*,form_d) vi(3),t1,t2	!�h���C���d���E��ʏo��
			write(4,form_d) vi(3),t1,t2	!�h���C���d���E�t�@�C���o��
			write(3,form_d) vi(3),t1,t2	!�h���C���d���E�t�@�C���o��
		endif
c		istp  	: �P��ԓ��ł̃X�e�b�v����
c		ict	: �S�̂�ʂ����X�e�b�v����
c		istp �� ict 
c##########circuit 2014/12/01(takahashi)######################################
		icstp=0
		i_or_c=1
		cflag=0				!��H�p�T�u���[�`�����ʃt���O
c
	if(jc_on.ne.5)then
		if(istp.gt.0)then
			ID1_ave = 0.0		
			IG1_ave = 0.0
			IS1_ave = 0.0
			do iave=istp-vave00,istp
				ID1_ave = ID1_ave + IDS1_stack(iave)
				IG1_ave = IG1_ave + IG1_stack(iave)
				IS1_ave = IS1_ave + ISS1_stack(iave)
			enddo
			ID1_ave=ID1_ave/vave00 		
			IG1_ave=IG1_ave/vave00 
			IS1_ave=IS1_ave/vave00 
		write(*,*) '�����d�����xId(0)=',ID1_ave

c
			IDS1(0) = ID1_ave*Wg		
			IG1(0)  = IG1_ave*Wg	
			ISS1(0) = -IS1_ave*Wg
		endif
c
		do iave=-1,jtp00
			IDS1_stack(iave)=0.0	
			IG1_stack(iave) =0.0
			ISS1_stack(iave)=0.0
		enddo
	endif
c#################################################################
		do istp=0,jstep-1
c	write(*,*) 'istp',istp
c		write(191,*) 'ict',ict,'jpnum',jpnum
c		nflag=1					!121102_nishida	�{�񂵃t���O
c#######################takahasi 2014 11/25########################
c	 pre_circuit���g���Ĕ񕽍t��Ԃł̉�H����������������A
c�@  	 pulse_circuit�ŕ��t��Ԃ̉�H������������
c##################################################################
			if(cflag.eq.0)then
				if((mod(istp-1,dtc).eq.0).or.(istp.eq.0)
     &				.and.(jc_on.ne.5))then
c	
				t=istp*dt
				CT=dtc*dt
c
				call pre_circuit(IG1,IDS1,ISS1,ID1,VO,VG1,V1,VD1,VS1,CT,t,
     &			icstp,Vd_add,Vg_add,idcount,istp,V_in,V_out,Wg,jc_on)
c			
c----------bias�̕ύX---------------------------------------------
c
				vi(3)=VD1(icstp)
				vi(2)=VG1(icstp)
				vi(1)=VS1(icstp)
c
				icstp = icstp + 1
c------------------------------------------------------------------
				endif
				cflag=2
			endif
c
	if(cflag.eq.2)then
c
			call sim_dt2(
     &              am_aniso,aff_aniso,hole_am_aniso,hole_aff_aniso,
     &			  dt,spnum,istp,ict,jpnum,
     &              dx,dz,xmax,zmax,lhet,iarea,dopem,twodeg,c_ratio,
     &              cxpart1,cxpart2,lxpart1,lxpart2,
     &              vb,vi,cur,cxpole1,cxpole2,lnpole1,lnpole2,
     &              melpos,jspsum,ncon,
     &              smh,hhm,hm,af,af2,af4,eps,eg,ec,bktq,
     &              de,swk,pgm,escat,iarg,iband,
     &              p,kp,u,cn,cp,cloud,hef_mesh,
     &              maceps,hescat,mtemp,
     &		      hef_scat,hcss,count,jpot,sstat,
     &			  cxrecess1,czrecess1,cxrecess2,czrecess2,
     &			  lxrecess1,lzrecess1,lxrecess2,lzrecess2,
c     &			  nava,cnava,n_scat,count_scat)
c     &			  nava,cnava,n_scat,count_scat,czpart2)
     &			  nava,cnava,n_scat,count_scat,czpart2,balis_flag,		!07/8/4 �s�����U��
     &			  n_scat_p,n_scat_n,									!08/8/6 �|��
     &			  efermi,avsumtel,avsumconc,i2max,E2,n2,i3max,E3,n3,	!08/8/6 �|��
     &		      avsumconc1,avsumtei11,ntab1,etab1,dn3,hiXL,
     &              epA,epB,epC,epA2,epB2,ecr,					!120126homma
     &              avsumconc_1,avsumconcA,Eth_table,				!09/2/19 �|��
     &			  am,aff,II_S,
     &			  basho_reflection,basho_roughness,	
     &			  pass,pass_r,x_start,xx_goal,scatpoint,
     &			  x_mean_free_path_sum,x_mean_free_path_count,			!120817sato
     &			  z_start,zz_goal,mean_free_path_sum,
     &			  mean_free_path_sum2,swk_rou,							!121029sato	
     &			  roughness1_countx,roughness1_counte,
c############circuit 100110(IKEDA)##############################################
     &			  IDS1_stack,IG1_stack,ISS1_stack,i_or_c,jc_on)	

c############circuit 100110(IKEDA)##############################################
c------��H�V�~�����[�V����----------
c
				t=istp*dt	
c
		if((mod(istp,dtc).eq.0).and.(istp.ne.0)
     &					.and.(jc_on.ne.5))then		
				if(istp.gt.0)then			
				ID1_ave = 0.0		
				IG1_ave = 0.0
				IS1_ave = 0.0
				do iave=istp-vave,istp
				ID1_ave = ID1_ave + IDS1_stack(iave)	
				IG1_ave = IG1_ave + IG1_stack(iave)
				IS1_ave = IS1_ave + ISS1_stack(iave)
				enddo
				ID1_ave=ID1_ave/vave 		
				IG1_ave=IG1_ave/vave 
				IS1_ave=IS1_ave/vave 
c
c
				IDS1(icstp) = ID1_ave*Wg		
				IG1(icstp)  = IG1_ave*Wg
				ISS1(icstp) = -IS1_ave*Wg
c
				cflag=1
				endif
		endif
	endif
c
c	if(jc_on.eq.1)then
c		call CIRCUIT_G(IG1, IDS1, ISS1, VG1, VD1, VS1, CT, t, icstp)	
c	endif
c	if(jc_on.eq.2)then
c		call Ex_CIRCUIT_G(IG1, IDS1, ISS1, VG1, VD1, VS1, CT, t,icstp,
c     &				Vd_add,Vg_add,idcount,istp,V_in,V_out,Wg,jc_on)	
c	endif
c	if(jc_on.eq.3)then
c		call S_Ex_CIRCUIT_G(IG1, IDS1, ISS1, VG1, VD1, VS1, CT, t,icstp,
c     &						Vd_add,Vg_add,idcount,istp,V_in,V_out,Wg)	
c	endif
c********(�i�s�g�Z�o)************************************************************
	if((jc_on.eq.4).and.(cflag.eq.1))then
		call pulse_circuit(IG1,IDS1,ISS1,ID1,VO,VG1,V1,VD1,VS1,CT,t,
     &			icstp,Vd_add,Vg_add,idcount,istp,V_in,V_out,Wg,jc_on)
		cflag=0
c********************************************************************************
c		
		if(mod(icstp,1000).eq.0)then
		 write(*,*) VG1(icstp),VD1(icstp),IG1(icstp),IDS1(icstp),icstp
		endif
	endif
c
c############circuit 14 11/20(takahashi)##############################################
	
			ict=ict+1		! ----ict �C���N�������g----
		enddo
		idvi(k2) = idvi(k2) + 1		!���[�v�����̓d����dvi��������
		!c_ratio = c_ratio + 0.0
	enddo
		idvi(k2) = idvi(k2) - ni2 - 1		!���[�v�����̓d����������
		idvi(k1) = idvi(k1) + 1		!���[�v�O���̓d����dvi��������
	enddo
!
!	file close
	close(100)
	close(101)
	close(102)
	close(103)
	close(104)
c	close(190)
c	close(191)
	close(560)
c
	stop
	end