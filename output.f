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
     &			tel1,tel2,tel3,efermi1,efermi2,efermi3,		!竹岸追加
     &			n_scat_p,n_scat_n,sscnt,					!竹岸追加
     &			avsumtel1,avsumconc1,avsumtei11,
     &            epA,epB,epC,epA2,epB2,eg,ecr,		!120126homma
     &            avsumteiA,avsumconcA,			!竹岸追加
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
c---変数配列用パラメータ---
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c
c--- シミュレーション条件 ---
	common /outp/joutc,jouti,jouta,jsave,jheat,jcput,joutpi,sw
	integer	joutc,jouti,jouta,jsave,jheat,jcput,joutpi,sw(8)
	integer	hcss
c
	common /step/jinit,jstat,jstep,jdisp
	integer	jinit,jstat,jstep,jdisp

c---基本パラメータ---
	real	dt,spnum,de(nenergy)
	integer	jpnum,ict
c---デバイス構造---
	real	dx,dz,xmax,zmax
	integer(2),dimension (nlayer)	::lhet
	integer(1),dimension (nlayer)	::iarea
	real,	dimension (npart)	:: cxpart1, cxpart2
	integer(2),dimension (npart)	::lxpart1,lxpart2
c---電極---
	real,	dimension (npole)	:: cxpole1,cxpole2
	integer(2),dimension (npole)	:: lnpole1,lnpole2
	integer(4),dimension (npole)	:: jspsum
	integer(1),dimension (npole)	:: melpos
c---領域別パラメータ---
	real,	dimension (nvalley,narea)	:: hhm,hm,af2,af4
	real,	dimension (narea)	:: eps
	real,	dimension (nvalley,narea) :: ec,eg
c---粒子状態---
	real	p(6,npmax)
	integer(1),dimension (3,npmax)	:: kp
c---デバイス内状態---
	real,	dimension (0:nx,0:nz)		:: u,cn
	real,	dimension ((nx+1)*(nz+1))	:: hef_mesh
	real,	dimension (nx+1,nscat,nvalley)	::hef_scat
	real,	dimension (0:nx,0:nz)	:: epA,epB,epC,epA2,epB2	!120126homma
c---入出力用変数---
	real	cur(npole)
	real(8)	cput(6)
c---衝突電離---
	integer,dimension (0:nx,0:nz,0:nvalley) :: nava
	integer,dimension (0:nx,0:nvalley) :: cnava
	integer,dimension (0:nx,nscat,nvalley) :: n_scat
	integer,dimension (0:nx) :: count_scat
c-----縮退効果-----
	integer,dimension (0:nx,nscat,nvalley) :: n_scat_p !(+)yama071223
	integer,dimension (0:nx,nscat,nvalley) :: n_scat_n !(-)yama071223
	real,	dimension (0:nx,0:nz)	:: tel1,efermi1
	real,	dimension (0:nx,0:nz)	:: tel2,efermi2
	real,	dimension (0:nx,0:nz)	:: tel3,efermi3
c
	real,	dimension(0:nx,0:nz)	::	avsumtel1
	real,	dimension(0:nx,0:nz)	::	avsumconc1
	real avsumconcA(0:nx)          !101221電子濃度チャネル平均2
	real(8),	dimension(0:nx,0:nz)	::	avsumtei11
	real(8) avsumteiA(0:nx)                 !101221電子エネルギーチャネル平均2
	real,save,allocatable 	:: av_tel1(:,:)
	real,save,allocatable 	:: av_conc1(:,:)
	real(8),save,allocatable 	:: av_tei11(:,:)
c
c---kdに関するパラメータ---	
	real adkx(0:nx,0:nz,0:nvalley)			!2014/12/01(takahashi)
	real adky(0:nx,0:nz,0:nvalley)
	real adkz(0:nx,0:nz,0:nvalley)
c
c-----yama追加 071020-----
	real,save,allocatable 	:: avef1(:,:),avtel1(:,:)
	real,save,allocatable 	:: avef2(:,:),avtel2(:,:)
	real,save,allocatable 	:: avef3(:,:),avtel3(:,:)
c
	integer ie2
	integer nemax4
	integer cff,bscat_xflag,fix_u					!電流ふらつき計算フラグ(takahashi)
	real ecnt(10,1,nvalley,250)
	real ecr(7,int(nemax/4))		!120126homma
c
	real(8),save,allocatable :: cn_scat_p(:,:,:) !071223yama
	real(8),save,allocatable :: cn_scat_n(:,:,:) !071223yama
c
	integer,dimension (0:nx,nscat,4) :: sscnt
c
c----バリスティックの計算----------
	integer	balis_scat,balis_n,balis_all
c-----散乱角の集計-------------------!08/1/21
	real(8),dimension (0:nx) ::	ccs,cncs
	real(8),dimension (0:nx,nscat) ::	ccs2,cncs2
c-----後方散乱の集計-------------------!08/1/28
	real(8),dimension (0:nx,0:nscat) ::	allback_scat
	real(8),dimension (0:nx,0:nscat) ::	back_scat			
c##########circuit 2014/12/01(takahashi)###############################
	real IG1_stack(-1:jtp00),IDS1_stack(-1:jtp00),ISS1_stack(-1:jtp00)
	integer	istp,i_or_c,jc_on
c---ローカル変数---
	real,	dimension (:,:,:),allocatable::hef_fin	!発熱のコメントアウト部分
	real,	dimension (:),allocatable::hef_allscat	!現在は使っていない
	integer	count
c	integer	i,j,k,n
	integer i,n					!j,kは現在使っていない
	integer	ix,iz,ixiz			!ixizは現在は使っていない
	integer	nxnz				!発熱関係の配列数。ただし、コメントアウトの場所
	integer iscat
	real(8)	sk,ei
	integer	kv,ken,kl,ka,ie
	integer,save :: ncount			!save忘れ　
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
	integer(8),save,allocatable :: iener(:,:)	!エネルギー別粒子数
	real(8),allocatable			:: dener(:,:)
	real(8),save	::	count2					!save忘れ
	real(8),save,allocatable :: vnava(:,:)
	real(8),save,allocatable :: cn_scat(:,:,:)

	real(8),save,allocatable 	:: en(:,:),cen(:,:)	!メッシュ別平均エネルギー
	integer(8),save,allocatable :: dn(:),cdn(:)	!メッシュ別粒子数

	real(8),save,allocatable	:: sk2(:),sq2(:),ei2(:)
	real(8),save,allocatable	:: en2(:,:),cen2(:,:)
	real(8), dimension (0:nx) :: cn_nx

	real(8),save,allocatable 	:: cvnxp(:),cvnxm(:)	!07/11/14 川端
	integer(8),save,allocatable :: cnvxp(:),cnvxm(:)		!07/11/14 川端
	real(8),allocatable			:: cvnxoutp(:),cvnxoutm(:)
	real(8),save	::	count3								!save忘れ
	integer	dist_ix			!07/11/22 川端
	integer,save,allocatable	:: dist_kxp(:,:),dist_kxm(:,:)	!07/11/22 川端
	real,save	::	dkx											!07/11/22 川端
	real balis_freq,balis_per									!07/11/22 川端

	real(8),save,allocatable 	:: ccs_all(:)	!08/1/21 川端
	real(8),save,allocatable 	:: ccs_scat(:,:)	!08/1/21 川端
c
c	sw(1)...電流
c	sw(2)...雑音
c	sw(3)...ポテンシャル・電子濃度,2-DEG層状態・バイナリー
c	sw(4)...エネルギー・谷占有率
c	sw(5)...電子速度
c	sw(6)...衝突電離頻度・散乱頻度
c	sw(7)...発熱率
c	sw(8)...cpu時間
c

c------ラフネス散乱--J.R. Watlingモデル(2012年春応物)-------------------------
      real delta,lambda,average
	integer split,count_reflection,count_roughness,
     &        vvv,xxx 
	integer,dimension (0:nx,2) :: basho_roughness,basho_reflection
      integer yyy

c-----	領域通過する粒子をカウント120817sato
	real(8),dimension (10) :: pass
	real	x_mean_free_path_sum
	real	mean_free_path_sum
	integer x_mean_free_path_count
	real,dimension (0:nx,3) :: 	mean_free_path_sum2		!平均自由行程x座標	 1:カウンタ,2:x方向平均自由行程,3:平均自由行程

c------ラフネス散乱--R.P. Joshiモデル(2012年秋)-------------------------
	integer,dimension (0:nx,narea,2) :: roughness1_countx			!x方向に対するラフネス散乱の数(チャネル上下)
	integer,dimension (nemax,narea,2) :: roughness1_counte		!エネルギーに対するラフネス散乱の数(チャネル上下)	

c ----( 電流 )----
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
		write(45 ,form) (cur(i),i=3,1,-1)		!Y,Sparameter用current_t.txt
		form = '(I8,5('','',I8))'
		write(11,form) ict,(jspsum(i),i=1,npole)	!jspsum.txt
		if(ict.ge.0)then						!Vcは本回しから出力を開始している
c			ファイル → unit=31 ... 'current_d.txt'
c			ファイル → unit=32 ... 'current_s.txt'
			write(31,*)cur(3)		!ドレイン電流
			write(32,*)cur(1)		!ソース電流
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
!	追加(2005/09/30 Hara)
! ----( 電流雑音のデバイス内分布計算 )----
	if((sw(2).gt.0).and.(cff.ne.1))then
		call current_fluctuation(jpnum,dx,dz,xmax,zmax,
     &fix_u,ec,cur,dt,cxpole1,cxpole2,
     & p,kp,hhm,hm,af4,af2,iarea,adkx,adky,adkz,spnum,cff,bscat_xflag)	
	endif
c
	endif
	endif
c
c	----( ポテンシャル・電子濃度出力 )----
	if(sw(3).gt.0)then
		if (.not. allocated(avu)) then	!最初だけ実行
			allocate(avu(0:nx,0:nz),avcn(0:nx,0:nz),
c-----縮退効果-----
     &				avef1(0:nx,0:nz),avtel1(0:nx,0:nz),
     &				avef2(0:nx,0:nz),avtel2(0:nx,0:nz),
     &				avef3(0:nx,0:nz),avtel3(0:nx,0:nz) )
c------------------
c-----08/8/6 竹岸-----
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
c-----縮退効果-----
			avef1 = 0.0; avef2 = 0.0; avef3 = 0.0
			avtel1 = 0.0; avtel2 = 0.0; avtel3 = 0.0
c------------------
		endif
		avu  = avu+u
		avcn = avcn+cn
		ncount= ncount+1
c-----縮退効果-----
		avef1 = avef1 + efermi1
		avef2 = avef2 + efermi2
		avef3 = avef3 + efermi3
		avtel1 = avtel1 + tel1
		avtel2 = avtel2 + tel2
		avtel3 = avtel3 + tel3
		av_tel1 = av_tel1 + avsumtel1	!08/8/6 竹岸 Γ谷
		av_conc1=av_conc1 + avsumconc1	!08/8/6 竹岸 Γ谷
		av_tei11=av_tei11 + avsumtei11	!08/8/6 竹岸 Γ谷
c	-------------
c
c	----( 経過ポテンシャル・電子濃度出力 )----
		if(modulo(ict,joutpi).eq.(joutpi-1))then
			!unit=25: potential.txt
			rewind 25 !; write(25,*)nx,nz,ict
			write(25,'(f12.7)') sngl(u)
			!unit=26: density.txt
			rewind 26 !; write(26,*)nx,nz,ict
			write(26,'(e15.7)') sngl(cn)
c	----( 2DEG部分ポテンシャル・電子濃度出力 )----
			if(modulo(ict,(jouta)).eq.0)then
				rewind	27
				rewind	28
			endif
			write(27,'(I7,x)')ict
			write(28,'(I7,x)')ict
c			iz = lhet(nlayer-1)	!要修正
			iz = lhet(nchannel2)	
			!unit=27: potential2d.txt
c			write(27,'(f12.7)') (sngl(u(0:nx,iz)))	!07/03/14
			!unit=28: density2d.txt
c			write(28,'(e15.7)') (sngl(cn(0:nx,iz)))
		endif
c	---( 平均粒子密度, ポテンシャル, 電界強度, 二次元電子濃度出力 )----
		if(modulo(ict,jouta).eq.(jouta-1))then
			! ----( ポテンシャル平均出力 )----
			avu    = avu/ncount
			!unit=38: potential_ave.txt
			write(38,*)nx,nz,ict
			write(38,'(f16.7)')avu
c
			! ----( 電界強度平均出力 )----
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
			! ----( 粒子濃度平均出力 )----
			avcn    = avcn/ncount
			!unit=39: cn_ave.txt
			write(39,*)nx,nz,ict
			write(39,'(e16.7)')avcn
c
			! ----( 二次元粒子濃度分布平均出力 )----
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
c-----縮退効果-----
	if(ict.gt.0)then
c-----(フェルミレベル平均出力)-----
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
c-----(電子温度平均出力)-----
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
c-----08/8/6 竹岸-----
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
c-----縮退効果-----
			avef1 = 0.0; avtel1 = 0.0
			avef2 = 0.0; avtel2 = 0.0
			avef3 = 0.0; avtel3 = 0.0
			av_tel1 = 0.0; av_conc1 = 0.0; av_tei11 = 0.0	!08/8/6 竹岸
c-------------------
c
			if(modulo(ict,jsave).eq.(jsave-1))then
				rewind	40;	write(40)	npmax,xmax,zmax,spnum
				rewind	41;	write(41)	jpnum,p,kp
			endif
		endif
! 追加--- 観測直前のデータ(05/09/30 Hara)
		if(ict.eq.-1)then
		!density_bef.txt 出力
		write(70,'(e15.7)') cn
		endif
	endif
c
c ----( エネルギー集計 )----
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
			ecnt = 0.0	;nemax4=nemax/4			!yama追加07/10/22
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

c----------最後のチャネルの電子のエネルギー分布--------2010/12/21hisa	 メッシュ幅1nm
              if(modulo(ict,jouta).eq.(jouta-1))then
	           if(kv.eq.1) then
c                     if((iz.ge.38).and.(iz.le.47))		then
					if((lhet(nchannel1)+1.le.iz).and. !channel 11/05/26原
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
	               if(ix.eq.110)		then	!要修正11/05/26原 ix=50→120126homma ix=110
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
					if((lhet(nchannel1)+1.le.iz).and. !channel 11/05/26原
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
	               if(ix.eq.110)		then	!要修正11/05/26原→120126homma ix=110
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

			! --- 粒子エネルギー分布集計
			enmesh(ix,iz)=enmesh(ix,iz)+ei		!メッシュ内エネルギー総和
			enx(ix)=enx(ix)+ei				! メッシュ内エネルギー総和

			! --- 谷占有率集計
			kvmesh(ix,iz,kv)=kvmesh(ix,iz,kv)+1		!メッシュ内谷別粒子数
			kvx(ix,kv)=kvx(ix,kv)+1				!谷別粒子数縦方向総和

			! --- メッシュ内粒子数カウント
			cnmesh(ix,iz)=cnmesh(ix,iz)+1	!メッシュ内粒子総和
			cnx(ix)=cnx(ix)+1				! 粒子数縦方向総和
c
c-----縮退効果-----
c-----エネルギー別粒子数カウント-----
c			ie2=max(min((nint(ei/(de(2)*4))+1),(nemax4)),1)
	if(kv.eq.1)then		!120126homma
		ie2=max(min((nint((ei-epA(ix,iz)+
     &	              (u(ix,iz)-0.7*(eg(1,ka)-eg(1,2)))) !!11/06/29原
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
c		if(nlayer-4.lt.kl.and.kl.le.nlayer-1) then  チャネル内でみるのもあり？
			if((ix.eq.10).and.(iz.eq.36))   ! iz=30へ変更(変更前iz=42)
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
			! --- エネルギーテーブル内分布集計
			ie=max(min((nint(ei/de(ken))+1),nemax),1)
			iener(ken,ie)=iener(ken,ie)+1	

			!---チャネル内---!!!!!!!!!!!!!!!!!!!!!!!!
c			if(nlayer-4.lt.kl.and.kl.le.nlayer-1) then	!channel 2006/12/22
			if(nchannel1.lt.kl.and.kl.le.nchannel2) then !channel 2011/04/08
				cenx(ix)=cenx(ix)+ei			!メッシュ内エネルギー総和
				ckvx(ix,kv)=ckvx(ix,kv)+1		!メッシュ内谷別粒子数
				ccnx(ix)=ccnx(ix)+1			!メッシュ内粒子総和
			endif

c			---一次元平均エネルギーenx---
			en(ix,1) = en(ix,1) + ei-ec(kv,ka)
			en(ix,2) = en(ix,2) + ec(kv,ka)
			en(ix,0) = en(ix,0) + ei
			dn(ix) = dn(ix) + 1
	        en2(ix,1) = en2(ix,1) +ei2(1)
			en2(ix,2) = en2(ix,2) +ei2(2)
			en2(ix,3) = en2(ix,3) +ei2(3)
c			---一次元チャネル内平均エネルギーcenx---
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
c---( エネルギー,谷占有率出力 )---
		if(modulo(ict,jouta).eq.(jouta-1))then
			forall(ix=0:nx,ccnx(ix).ne.0)			!他で計算する。他に別の粒子数から
			cn_nx(ix)=dble(ccnx(ix))/count2
			end forall
c		---一次元平均エネルギー出力enx---
			forall(ix=0:nx,dn(ix).ne.0)
				en(ix,0:2) = en(ix,0:2)/dble(dn(ix))	!07/1/22　川端				
			end forall
			rewind 68
			do ix=0,nx
				write(68,'(3(F8.4))')en(ix,0:2)
			enddo
			forall(ix=0:nx,dn(ix).ne.0)	!channel 2007/2/2
				en2(ix,1:3) = en2(ix,1:3)/dble(dn(ix)) !07/1/22　川端			
			end forall
			rewind 88
			do ix=0,nx
				write(88,'(3(F8.4))')en2(ix,1:3)
			enddo
c		---一次元チャネル内平均エネルギー出力cenx---
			forall(ix=0:nx,cdn(ix).ne.0)
				cen(ix,0:2) = cen(ix,0:2)/dble(cdn(ix))		!07/1/22　川端
			end forall
			rewind 69
			do ix=0,nx
				write(69,'(3(F8.4))')cen(ix,0:2)
			enddo
			forall(ix=0:nx,cdn(ix).ne.0)	!channel 2007/2/2	!07/1/22　川端
				cen2(ix,1:3) = cen2(ix,1:3)/dble(cdn(ix))		!07/1/22　川端	
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
c---( 一次元エネルギー,谷占有率出力 )---
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

c---( 二次元エネルギー,谷占有率出力 )---
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

			! ----( 粒子濃度分布平均出力 )----
			dcnmesh = dcnmesh*spnum/dx/dz/count2
			!unit=37: density_ave.txt
			write(37,*)nx,nz,ict
			write(37,'(e16.7)')dcnmesh
			cnmesh = 0; count2 = 0;
			deallocate(dcnmesh)
c
c-----縮退効果-----
c-----エネルギー別粒子数出力-----
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
c----( 速度分布出力 )----
	if(sw(5).gt.0)then
	if(modulo(ict,jouti).eq.(jouti-1))then
		if (.not. allocated(v)) then
c			v_mesh...メッシュ別速度和, sv_mesh...メッシュ別平均速度
c			1...x成分, 2...y成分, 3...z成分, 0...絶対値
			allocate(v(1:3))
 			allocate(vnmesh(0:3,0:nx,0:nz))
			allocate(nvmesh(0:nx,0:nz))
			allocate(vnx(0:3,0:nx),cvnx(0:3,0:nx))
			allocate(nvx(0:nx),cnvx(0:nx))
			allocate(cvnxp(0:nx),cvnxm(0:nx))		!07/11/14 川端
			allocate(cnvxp(0:nx),cnvxm(0:nx))		!07/11/14 川端
			allocate(dist_kxp(0:nx,0:500),dist_kxm(0:nx,0:500))
c	-----08/3/24 竹岸 構造を大きくするときには上をcし下のcをとる-----
c			allocate(dist_kxp(0:nx,0:2000),dist_kxm(0:nx,0:2000))
			vnmesh = 0.0
			nvmesh = 0
			vnx = 0.0	;cvnx = 0.0
			nvx = 0		;cnvx = 0
			cvnxp = 0.0	;cvnxm = 0.0				!07/11/14 川端
			cnvxp = 0	;cnvxm = 0					!07/11/14 川端
			dkx = 5.0e-8							!2e7で規格化
			dist_kxp = 0;dist_kxm = 0	
			count3 =0.0								!08/8/6 竹岸
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
c			---メッシュ別平均速度---
			vnmesh(0,ix,iz)=vnmesh(0,ix,iz)+vave
			vnmesh(1:3,ix,iz)=vnmesh(1:3,ix,iz)+v(1:3)
c			---一次元平均速度---
			vnx(0,ix)=vnx(0,ix)+vave
			vnx(1:3,ix)=vnx(1:3,ix)+v(1:3)
			! --- メッシュ別速度集計数カウント
			nvmesh(ix,iz)=nvmesh(ix,iz)+1
			nvx(ix)=nvx(ix)+1
			!---チャネル内---
c			if(nlayer-4.lt.kl.and.kl.le.nlayer-1) then
			if(nchannel1.lt.kl.and.kl.le.nchannel2) then !channel 11/04/07原
				cvnx(0,ix)=cvnx(0,ix)+vave
				cvnx(1:3,ix)=cvnx(1:3,ix)+v(1:3)
				cnvx(ix)=cnvx(ix)+1

				dist_ix = nint(p(1,n)*dkx)	!07/11/22 川端
c---------------------------------------------------------- !07/11/14 川端
c
c-----cvn xp(x):プラス方向の速度成分総和
c-----cnv xp(x):プラス方向の速度を持つ粒子数総和
c-----cvn xm(x):マイナス方向の速度成分総和
c-----cnv xm(x):マイナス方向の速度を持つ粒子数総和
c-----cvnxoutp(ix):プラス方向速度出力
c-----cvnxoutm(ix):マイナス方向速度出力


				if(v(1).gt.0) then
					if(dist_ix.ge.100) dist_ix=100	!yama080105
					cvnxp(ix) = cvnxp(ix)+v(1)
					cnvxp(ix)=cnvxp(ix)+1
					dist_kxp(ix,dist_ix) = dist_kxp(ix,dist_ix) + 1	!07/11/22 川端
				elseif(v(1).lt.0) then
					if(dist_ix.le.-100) dist_ix=-100	!yama080105
					cvnxm(ix) = cvnxm(ix)+v(1)
					cnvxm(ix)=cnvxm(ix)+1
					dist_ix = -dist_ix						!07/11/22 川端
					dist_kxm(ix,dist_ix) = dist_kxm(ix,dist_ix) + 1
				endif
c----------------------------------------------------------
			endif
		enddo
		count3 = count3 + 1.0
c	
c	  ---電子速度出力---
		if(modulo(ict,jouta).eq.(jouta-1))then
			allocate(vnxout(0:nx,0:3),cvnxout(0:nx,0:3))
			allocate(cvnxoutp(0:nx),cvnxoutm(0:nx))
			do ix=0,nx
c				if(nvx(ix).ne.0)then
				if(cnvx(ix).ne.0)then	!08/1/22	川端
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
			vnmesh = 0.0	!リセット
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
c----( 衝突電離頻度 )----
c	nava...メッシュ別衝突電離頻度, vnava...谷別衝突電離頻度
	if(sw(6).gt.0)then
	if(modulo(ict,jouta).eq.(jouta-1))then
	  if (.not. allocated(vnava)) then
 		allocate(vnava(0:nx,1:nvalley))
 		allocate(cn_scat(0:nx,1:nscat,nvalley))
		allocate(ccs_all(0:nx))
		allocate(ccs_scat(0:nx,1:nscat))
c-----yama追加-----
	 	allocate(cn_scat_p(0:nx,1:nscat,nvalley))
	 	allocate(cn_scat_n(0:nx,1:nscat,nvalley))
		vnava = 0.0
		cn_scat = 0.0
		ccs_all = 0.0
		ccs_scat = 0.0
c-----yama追加-----
		cn_scat_p = 0.0
		cn_scat_n = 0.0
	endif
		!---メッシュ別衝突電離頻度---
		!unit=61: nava.txt
		rewind 61
		write(61,'(I9)') nava(0:nx,0:nz,0)
		nava = 0
		!---一次元衝突電離頻度---
		!unit=60: ava_zave.txt
		rewind 60
		write(60,*)nx,dx,ict
		write(60,'(e12.4)') cnava(0:nx,0)*spnum/dx/dz/jouta/dt	!times/s/m3
		!unit=62: cnave.txt
		rewind 62
		write(62,*)nx,dx,ict
		write(62,'(I9)') cnava(0:nx,0) !times
		!---谷別衝突電離頻度---
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

		!----( チャンル内散乱頻度出力 )----
		!---一次元散乱頻度---
			do kv = 1, nvalley
				do iscat = 1, nscat
				  forall(ix=0:nx,count_scat(ix).ne.0)
c				  cn_scat(ix,iscat,kv) = dble(n_scat(ix,iscat,kv))
c   &								/dble(count_scat(ix))*spnum/dx/dz/jouta/dt
				  cn_scat(ix,iscat,kv) = dble(n_scat(ix,iscat,kv))
     &								*spnum/dx/dz/jouta/dt
c-----08/8/6 竹岸-----
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
c-----縮退効果-----
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
c-----縮退効果-----
				write(136,form) cn_scat_p(ix,iscat,kv),n_scat_p(ix,iscat,kv),
     &							count_scat(ix)
				write(137,form) cn_scat_n(ix,iscat,kv),n_scat_n(ix,iscat,kv),
     &							count_scat(ix)
c------------------
		enddo
		enddo
		enddo
c
c-----yama追加-----
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
c----バリスティックの計算----------
		balis_freq = float(balis_scat)/float(balis_all)
		balis_per =	float(balis_all-balis_n)/float(balis_all)*100
		
		rewind 110
		write(110,'(2(e12.5))') balis_freq,balis_per

		balis_all = 0; balis_scat = 0
		balis_n = 0; balis_freq = 0
		balis_per = 0
c--------------------------
c-----散乱角の計算---------
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
c-----後方散乱のカウント--------
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
		n_scat_p = 0					!08/8/6 竹岸
		n_scat_n = 0					!08/8/6 竹岸
	endif
	endif
c
c%%%%%%%%%% 発熱率関係出力 %%%%%%%%%%
c	ファイル → unit=30 ... 'heat_generation.txt'
c	ファイル → unit=33 ... 'heat_generation_mesh.txt'
c	ファイル → unit=36 ... 'heat_generation_all.txt'
	if(sw(7).gt.0)then
	if((modulo(ict,jheat).eq.(jheat-1)).and.(hcss.eq.1))then

			buff=q*spnum/(dx*dz*dt*float(count))
			hef_mesh=hef_mesh*buff	!hef_mesh:全体配列
			hef_scat=hef_scat*buff	!hef_scat:全体配列
			rewind	30
			write(30,'(E15.7)') hef_mesh		!heat_generation.txt  出力
			write(36,'(E15.7)') hef_mesh		!heat_generation_all.txt  出力

			hef_mesh = 0.0		!全体配列指定

			rewind	33
			write(33,'(E15.7)') hef_scat		!heat_generation_scat.txt  出力

			hef_scat = 0.0		!全体配列指定

c			nxnz = (nx+1)*(nz+1)
c			allocate (hef_fin(nxnz,nscat,nvalley))
c			hef_fin=hef_scat*buff	!hef_fin,hef_scat:全体配列

!SMP$	ASSERT (ITERCNT(70))	!スパコン用最適化補助命令
c			do k=1,nvalley
c			do j=1,nscat
c				hef_fin( 0,0:nz,j,k)=hef_fin( 0,0:nz,j,k)*2.0
c				hef_fin(nx,0:nz,j,k)=hef_fin(nx,0:nz,j,k)*2.0
c				hef_fin(0:nx, 0,j,k)=hef_fin(0:nx, 0,j,k)*2.0
c				hef_fin(0:nx,nz,j,k)=hef_fin(0:nx,nz,j,k)*2.0
c			enddo
c			enddo

c			hef_scat = 0.0		!全体配列指定

c			rewind	33			!heat_generation_scat.txt 出力
c			write(33,'(E15.7)') hef_fin
c			deallocate (hef_fin)
c
			count = 0
	endif
	endif
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c------ラフネス散乱--J.R. Watlingモデル(2012年春応物)-------------------------
      if(modulo(ict,jouta).eq.(jouta-1)) then

	    open(1919,file='roughness_data.txt')
          write(1919,*) '角度分割数',split
          write(1919,*) '凹凸高さ',delta
		write(1919,*) '凹凸相関長',lambda
	    write(1919,*) '反射数',count_reflection
          write(1919,*) '散乱数',count_roughness
          write(1919,*) '平均角度変化量',average

          open(0721,file='reflection_count_data.txt')
	    do vvv=0,nx
	        write(0721,*) (basho_reflection(vvv,yyy),yyy=1,2)
	    end do

	    open(6969,file='roughness_count_data.txt')
	    do xxx=0,nx
	        write(6969,*) (basho_roughness(xxx,yyy),yyy=1,2)
          end do

	end if

c-----	領域通過する粒子をカウント120817sato
c-----pass(1)A外から領域通過,(2)A外から領域内,(3)領域内からB外,(4)領域内からA外,
c-----pass(5)領域内連続,(6)領域連続後B外,(7)領域連続後A外	 (10)領域内連続(仮)
c-----pass(8)B外から領域通過,(9)B外から領域内

      if(modulo(ict,jouta).eq.(jouta-1)) then

		x_mean_free_path_sum = x_mean_free_path_sum 
     &							/ x_mean_free_path_count
		mean_free_path_sum = mean_free_path_sum 
     &							/ x_mean_free_path_count

		write(852,*) 'A外から領域通過',',',pass(1)
		write(852,*) 'A外から領域内',',',pass(2)
		write(852,*) '領域内からB外',',',pass(3)
		write(852,*) '領域内からA外',',',pass(4)
		write(852,*) '領域内連続',',',pass(5)
		write(852,*) '領域連続後B外',',',pass(6)
		write(852,*) '領域連続後A外',',',pass(7)
		write(852,*) 'B外から領域通過',',',pass(8)
		write(852,*) 'B外から領域内',',',pass(9)
		write(852,*) '平均自由行程x',',',x_mean_free_path_sum
		write(852,*) '平均自由行程',',',mean_free_path_sum
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
		
		x_mean_free_path_sum = 0.0	!初期化
		mean_free_path_sum = 0.0	!初期化
		x_mean_free_path_count = 0
		mean_free_path_sum2 = 0.0	!初期化
	endif

c------ラフネス散乱--R.P. Joshiモデル(2012年秋)-------------------------
      if(modulo(ict,jouta).eq.(jouta-1))then
		do ix=0,nx
		write(1210,*)ix,',',roughness1_countx(ix,1:narea,1)		!チャネル上ヘテロ界面x方向ラフネス回数
		write(1211,*)ix,',',roughness1_countx(ix,1:narea,2)		!チャネル下ヘテロ界面x方向ラフネス回数
		enddo
		do ie=1,nemax
		write(1212,*)ix,',',roughness1_counte(ie,1:narea,1)		!チャネル上ヘテロ界面エネルギーラフネス回数
		write(1213,*)ix,',',roughness1_counte(ie,1:narea,2)		!チャネル下ヘテロ界面エネルギーラフネス回数
		enddo
		 
		write(1210,*)''
		write(1211,*)''
		write(1212,*)''
		write(1213,*)''
		
		roughness1_countx = 0	!初期化
		roughness1_counte = 0	!初期化
	endif

c ----( cpu時間配分出力 )----
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