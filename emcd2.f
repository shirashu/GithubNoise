	subroutine emcd2(
     &                 am_aniso,aff_aniso,hole_am_aniso,hole_aff_aniso,
     &				 dt,de,jpnum,
     &                 dx,dz,xmax,zmax,lhet,iarea,twodeg,dopem,
     &                 cxpole1,cxpole2,melpos,jspsum,
     &                 smh,hm,hhm,af,af2,af4,eps,eg,bktq,
     &                 swk,pgm,escat,iarg,iband,
     &                 p,kp,u,cn,hef_mesh,
     &                 mtemp,hescat,hef_scat,hcss,
     &				 cxr1,czr1,cxr2,czr2,lxr1,lzr1,lxr2,lzr2,czpart2,	!07/8/4 不純物散乱
     &				 nava,cnava,n_scat,count_scat,ict,ep_flag,	!Effective Potential用追加　Hara 2006/12/09
     &				 balis_flag,balis_scat,balis_n,balis_all,
     &				 ccs,ccs2,cncs,cncs2,allback_scat,back_scat,	!08/1/21
     &				 ec,adkx,adky,adkz,efermi,n_scat_p,n_scat_n,	!08/8/6 竹岸
     &			     sscnt,avsumconc,avsumtel,dn3,hiXL,
     &                 epA,epB,epC,epA2,epB2,Eth_table,	!09/2/19 竹岸 !120126homma
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
c===引数===

c---変数配列用パラメータ---
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c---基本パラメータ---
	integer	ict		!2006/12/09 Hara
	real	dt
	integer	jpnum
	integer	ep_flag		!2006/12/09 Hara
c---デバイス構造---
	real	dx,dz,xmax,zmax
	integer(2)	lhet(nlayer)
	integer(1)	iarea(nlayer)
	real	twodeg(nlayer)
	real	dopem(0:nx,0:nz)
	real czpart2(npart)	!07/8/4 不純物散乱
	integer ipart,kpart	!07/8/4 不純物散乱
c---電極---
	real	cxpole1(npole),cxpole2(npole),melpos(npole)
	integer(4)	jspsum(npole)
c---領域別パラメータ---
	real,	dimension (nvalley,narea)	:: smh,hhm,hm,af,af2,af4
	real	eps(narea),eg(nvalley,narea),bktq(ntenum)
	real	dltec(nvalley,narea) !band offsets by nextnano
c---散乱パラメータ---
	real	de(nenergy)
c	real,	dimension (nscat,nemax,nvalley,nenergy,ntenum,narea)::	swk
c	real,	dimension (nvalley,nenergy,narea)	:: pgm
	real,	dimension (nscat,nemax,nvalley,nenergy,ntenum,npart)::swk	!07/8/4 不純物散乱
	real,	dimension (nvalley,nenergy,npart)	:: pgm	!07/8/4 不純物散乱
	real,	dimension (nscat,nvalley,narea)	:: escat
	integer(1),dimension (nscat,nvalley,narea)	:: iband
	integer(1),dimension (nscat,narea)	:: iarg
	real,	dimension (npart)	:: dn3					!09/2/19 竹岸
c---粒子状態---
	real	,dimension   (6,npmax)	:: p
	integer(1),dimension (3,npmax)	:: kp
	real	,dimension   (6)	:: sp				 !121009			
	integer(1),dimension (6)	:: skp				 !121009				
c---デバイス内状態---
	real,	dimension (0:nx,0:nz)	:: u,cn,hef_mesh
	real,	dimension (0:nx,0:nz)	:: epA,epB,epC,epA2,epB2	!120126homma
	real,save,allocatable	::	u_eff1,u_eff2,u_eff3	!2006/12/09 Hara
	dimension	::	u_eff1(:,:),u_eff2(:,:),u_eff3(:,:)
c---発熱率用配列---
	integer(2),dimension (0:nx,0:nz) :: mtemp
	real,	dimension (nscat,nvalley,narea) :: hescat
	real hef
	real,	dimension (0:nx,nscat,nvalley) :: hef_scat
	integer hcss
c----(リセス)---
	integer ii
	real, dimension (nrecess) :: cxr1,czr1,cxr2,czr2
	integer(2),dimension (nrecess) :: lxr1,lzr1,lxr2,lzr2
c---(衝突電離)---
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

c----バリスティックの計算----------
	integer(1),dimension (0:npmax) ::	balis_flag		!07/11/22
	integer(1) balis_flag2
	integer	balis_scat,balis_n,balis_all
c-----散乱角の集計-------------------!08/1/21
	real(8),dimension (0:nx) ::	ccs,cncs
	real(8),dimension (0:nx,nscat) ::	ccs2,cncs2	
c-----後方散乱の集計-------------------!08/1/28
	real(8),dimension (0:nx,0:nscat) ::	allback_scat
	real(8),dimension (0:nx,0:nscat) ::	back_scat	
c-----縮退効果-----
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
	real avsumconc(0:nx,0:nz,nvalley)	!08/8/6 竹岸
	real avsumtel(0:nx,0:nz,nvalley)	!08/8/6 竹岸
	integer rflag,swrej
	integer pflag
c
c----ローカル変数----
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
	integer	ix,iz,iv	!08/8/6 竹岸
	character(80) form
	real	pdx,pdz
	integer	nstat,nend
	integer iii
c-----散乱カウント-------------------!10/05/07
	integer	  count_flag,six,siscat,skv
c
cccccリジェクション・衝突電離 補正用　11/08原
	integer iiflag		!追加　原
	integer iiix,iiiz
cccccccccccccccccccccccccccccccccccccccccccc

c------ラフネス散乱--J.R. Watlingモデル(2012年春応物)-------------------------
      real delta,lambda,average
	integer split,count_reflection,count_roughness,
     &        vvv,xxx 
	integer,dimension (0:nx,2) :: basho_roughness,basho_reflection

c-----------------------------------------------------------------------------------------
c-----	領域通過する粒子をカウント120817sato
	real(8),dimension (10) :: pass,pass_r		!カウンタ，rはリジェクション対策
	real	x_start(npmax)			!driftを始める位置を保存する配列
	real	x_start_1				!ループさせる変数
	real	x_start_r				!rejection対策
	real	xx_goal					!driftが終わる位置を保存する変数
	integer	scatpoint			!散乱．空散乱判定
	real	x_mean_free_path_sum
	integer x_mean_free_path_count

	real	z_start(npmax)			!driftを始める位置を保存する配列
	real	z_start_1				!ループさせる変数
	real	z_start_r				!rejection対策
	real	zz_goal					!driftが終わる位置を保存する変数
	real	mean_free_path_sum
	real,dimension (0:nx,3) :: 	mean_free_path_sum2		!平均自由行程x座標	 1:カウンタ,2:x方向平均自由行程,3:平均自由行程
	integer xcenter				!x_startとxx_goalの中心点

c------ラフネス散乱--R.P. Joshiモデル(2012年秋)-------------------------
	real	swk_rou(nemax,nvalley,nenergy,narea)		!121029sato
	integer,dimension (0:nx,narea,2) :: roughness1_countx
	integer,dimension (nemax,narea,2) :: roughness1_counte

c=========================================================================================================
c
c
	if (.not. allocated(dhet)) then	!最初だけ実行
		swrej=5 !yama仮の初期値→0以外の値ならOK
		allocate (dhet(0:nlayer))
		dhet(0) = -huge(dhet(0))
		dhet(1:nlayer-1) = zmax*float(lhet(1:nlayer-1))/nz	 !ヘテロ界面位置
		dhet(nlayer) = huge(dhet(nlayer))
c
		allocate (u_eff1(0:nx,0:nz),u_eff2(0:nx,0:nz),u_eff3(0:nx,0:nz))	!2006/12/09 Hara
		u_eff1=0.0
		u_eff2=0.0
		u_eff3=0.0
c----------------------バリスティックの計算 07/11/22 川端--------------
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
c-----縮退効果-----
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
c					!jpotに合わせて実効ポテンシャルを計算(2006/05/20)2006/12/09
c
c		追加  実効ポテンシャル計算　2005/12/02 Hara
		call qeffect(dx, dz, lhet, iarea, 
     &					hhm, eg, bktq, lxr1,lzr1,lxr2,lzr2,
     &					u, u_eff1,u_eff2,u_eff3, mtemp,dltec,ec)
c
		call eff_out(ict, u_eff1,u_eff2,u_eff3,epA,epB,epC,epA2,epB2)	!120126homma	
c					!実効ポテンシャルの出力
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
c	----(衝突電離)----
 1000	continue
c
	do n=nstat,nend		!1~nendまでの粒子繰り返し
c
c-----縮退効果-----
c	flag_d=0
c	flag_s=0
	iak=0
c	nak=0
	pflag = 0
c 2001 continue			!121009
	count_flag = 0		!10/05/07 散乱カウント	!121009
		
c		iiflag = 0		!リジェクション・衝突電離回避用フラグ11/08	 !コメントアウト121009

	jak=0.0
c------------------
c
c		n個目の粒子の・・・
		kv	= kp(1,n)	!谷No.
		ken	= kp(2,n)	!エネルギーテーブルNo.
		kl  = kp(3,n)	!層No
		ka	= iarea(kl)	!素材No．
		kl2 = kl		!kl2:粒子が進入(反射)する場所の層
		ka2	= ka		!ka2:粒子が進入(反射)する場所の素材
c
c		kpn	= kp(n)-1
c		ka	= kpn/(nvalley*nenergy)+1		!谷No.
c		ken	= mod(kpn,nenergy)/nvalley+1	!エネルギーテーブルNo.
c		kv	= mod(kpn,nvalley)+1			!素材No．
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
c---	①電子がいる層の不純物濃度を調べる ----	  !07/8/4 不純物散乱
c
		do ipart=1,npart
			if((ipart.eq.npart-2).or.(ipart.eq.npart-1))cycle  !電極位置
			if(p(6,n).le.czpart2(ipart))then
				kpart=ipart
				exit	
			else
				kpart=npart
				if(p(6,n).gt.czpart2(npart))then
					write(*,*)'emcdで不純物濃度エラー1'
					stop
				endif
			endif
		enddo
c--------------------------------------------------
c
		iflag = 0
		balis_flag2 = balis_flag(n)		!バリスティックの計算 07/11/22 川端
c
c
	loop : do while(kp(1,n).ne.0)
c
c
c************粒子ドリフト部分**************************************************
c#####電界強度の取得############################################################
c		call field(x,z,pdx,pdz,u,fx,fz,dhet,kl,
c     &				cxr1,czr1,cxr2,czr2,lxr1,lzr1,lxr2,lzr2)
c
c	実効ポテンシャルから電界強度を取得する	  2005/12/02 Hara 変更2006/12/09
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
     &						de, pgm(1,1,kpart),	!07/8/4 不純物散乱
     &				        af2(kv,ka), af4(kv,ka),
     &				        hm(kv,ka), hhm(kv,ka),
     &				        akx,aky,akz,ts,x,z,t1,tau,
     &			            kv,ken,kl,kl2,iflag,ie,ei,sk,
     &						kpart,ka,n)									!07/8/4 不純物散乱


c-----縮退効果 drift成分(k-kd)の各メッシュ・谷での総和-----
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
     &				   cxr1,czr1,cxr2,czr2) !2011/3/25原				
c
c-------------------バリスティックの計算 07/11/22 川端--------------------
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
c************イベント部分******************************************************
c	----粒子が消滅した場合----
		if(kv.eq.0)then
	exit loop
		endif
c
c-----初期状態の情報を格納-----	!121009
c	p(1-6,n) ... 粒子状態(1-3:k座標(kx,ky,kz)[m^-1],4:散乱時刻[s],5-6:位置(x,y)[m])
c	kp(1-6,n) ... 所属谷
		sp = 0	;	skp = 0			
		sp(1)  = akx
		sp(2)  = aky
		sp(3)  = akz
		sp(4)  = ts	
		sp(5)  = x
		sp(6)  = z
		skp(1) = kv		!谷No.
		skp(2) = ken		!エネルギーテーブルNo.
		skp(3) = kl		!層No
c------------------------------
 2001 continue

		iiflag = 0		!リジェクション・衝突電離回避用フラグ11/08
		
c-----終状態を初期状態の情報に戻す----- !121009
c	p(1-6,n) ... 粒子状態(1-3:k座標(kx,ky,kz)[m^-1],4:散乱時刻[s],5-6:位置(x,y)[m])
c	kp(1-6,n) ... 所属谷
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
c	----τがdtに達した場合(iflag=2)----
c	----散乱イベントが発生した場合(iflag=0)----
		if((iflag.eq.0).or.(iflag.eq.2)) then
c
c
c			---	分岐 ---
			if(iflag.eq.2)then
	exit loop
			endif
c
			if(iflag.eq.0)then
c	----散乱イベントが発生した場合(iflag=0)----
c
				den=ei/de(ken)-float(ie-1)		!de:ieとeiの差（ズレ）
				if((den.gt.1.0).or.(den.lt.0.0))then
					form="('denの値が不正です(scat) den=',f,'ei=',e)"
					write(* ,form)den,ei
					write(99,form)den,ei
					
c	 check_point
				endif

				ix = min(nx,max(0,nint(x*pdx)))
				iz = min(nz,max(0,nint(z*pdz)))
				mtp = mtemp(ix,iz)		!粒子のいるメッシュの温度
c
c-----yama追加-----
c	散乱イベントが発生した場合、いつでもカウント!
c	初め +方向のkxを持つ粒子	 siflag=2　
c	初め -方向のkxを持つ粒子	 siflag=1
			ix3 = ix
			if(akx.gt.0.0) then
				siflag = 2
			else
				siflag = 1
			endif
c------------------
c	 check_point	error message
c
c	if((bscat_xflag.eq.1).and.(p(5,n).ge.(cxpole2(2)+65e-9)))then	!粒子がドレイン領域にいたら
	if((bscat_xflag.eq.1).and.(p(5,n).ge.(cxpole2(2))))then	!粒子がドレイン領域にいたら
	if(ei.ge.0.1)then
	goto 20
	endif
	endif
		call scat(
     &				 am_aniso,aff_aniso,hole_am_aniso,
     &				 hole_aff_aniso,smh(1,ka),
     &				 af(1,ka),
     &				 iiflag,iiix,iiiz,	!!I.I.補正用追加 11/08原
c     &			     swk(1,1,1,ken,mtp,ka), escat(1,1,ka),
     &			     swk(1,1,1,ken,mtp,kpart), escat(1,1,ka),	!07/8/4 不純物散乱
     &		         iarg(1,ka), iband(1,1,ka),
     &	             eps(ka), bktq(mtp), dopem(ix,iz),
     &		         akx, aky, akz, kv, kl,
     &                 sk, ei, ef, ie, den,
     &				 hescat(1,1,ka),hef,hcss,jp,iscat,
c     &				 dx,dz,jpnum,kp,p,pgm(1,1,ka),n,t1,
     &				 dx,dz,jpnum,kp,p,pgm(1,1,kpart),n,t1,	!07/8/4 不純物散乱
     &				 nava,cnava,x,z,count_flag,six,siscat,skv,			!10/05/07 散乱カウント
     &				 cxpole1,cxpole2,
     &				 balis_scat,balis_flag2,balis_n,		!07/03/15
     &				 ccs,ccs2,cncs,cncs2,allback_scat,back_scat,		!08/1/28
     &				 n_scat_p,n_scat_n,dn3(kpart),						!09/2/19 竹岸
     &				lxr1,lxr2,ec,eg,hm(kv,ka),ka,
     &			     hiXL,u_eff1,u_eff2,u_eff3,ix,iz,Eth_table,ken,u,			!09/2/19 竹岸
     &				lhet,de,											!120330	佐藤
     &				scatpoint,xx_goal,zz_goal,x_start,z_start,II_S,		!120817sato
     &				dltec)
c
c-----yama追加-----
c	散乱イベントが発生した場合、ヘテロ内のみ!!カウント
c	散乱後 +方向のkxを持つ粒子	 sfflag=3
c	散乱後 -方向のkxを持つ粒子	 sfflag=1
20		if((iscat.ne.0).and.
     &		(nchannel1.lt.kl.and.kl.le.nchannel2))then !channel2011/05/25原
c     &				(nlayer-4.lt.kl).and.(kl.le.nlayer-1)) then  !チャネル内
c     &				(kl.gt.nlayer-3).and.(kl.le.nlayer-2)) then  !チャネル内
				if(akx.gt.0.0) then
					sfflag = 3 
				else
					sfflag = 1
				endif
c
				smpf = siflag + sfflag -1	!smpf = 1から4
				sscnt(ix3,iscat,smpf)=sscnt(ix3,iscat,smpf)+1	!各ixの散乱方向種類ssflag
				siflag = 0 ; sfflag = 0
				smpf = 0
c
c-----kxが 初め(+) 後(+) >> smpf = 4
c-----kxが 初め(-) 後(+) >> smpf = 3
c-----kxが 初め(+) 後(-) >> smpf = 2
c-----kxが 初め(-) 後(+) >> smpf = 1
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
				ts = dble(t1)-dble(log(rnd())*pgm(kv,ken,kpart))	!07/8/4 不純物散乱			
c=============エラー処理========================================
				if(sngl(ts).lt.t1)then
					write( *,*)'τ論理エラー2(emcd2-散乱後)'
					write(99,*)'τ論理エラー2(emcd2-散乱後)'
					kv=0
	exit loop
				endif
c===============================================================
c
			endif
c
c	----ヘテロ障壁衝突イベントが発生(iflag=1)----
		elseif(iflag.eq.1)then
c		----(リセス)----
		  do ii = 1,nrecess	
			if((x.gt.cxr1(ii)).and.(x.lt.cxr2(ii)).and.(z.eq.czr2(ii)))then
			  akz=-akz
c
				if ((ii.eq.2)
     &			  .and.(czr2(1).eq.czr2(2))) then !緊急回避特異ケース11/04/08原
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
c------ラフネス散乱--R.P. Joshiモデル(2012年秋)-------------------------
c				call border_roughness1(hhm,af,af2,af4,eg,
c     &				akx,aky,akz,kv,kl,kl2,ka,iarea,
c     &				lhet,min(nz,max(0,nint(z*pdz))),
c     &				swk_rou,ken,min(nx,max(0,nint(x*pdx))),de,		!121029sato
c     &				roughness1_countx,roughness1_counte)
c---------------------------------------------------------------------------

c------ラフネス散乱--J.R. Watlingモデル(2012年春応物)-------------------------
c				call border_roughness2(hhm,af,af2,af4,eg,
c    &				akx,aky,akz,kv,kl,kl2,ka,iarea,
c     &				lhet,min(nz,max(0,nint(z*pdz))),
c     &				min(nx,max(0,nint(x*pdx))),
c     &            xxx,vvv,basho_reflection,basho_roughness,
c     &            split,delta,lambda,count_reflection,
c     &            count_roughness,average,
c     &			epA,u)			!ラフネス散乱用ep障壁
c-----------------------------------------------------------------------------

c										     ! ドリフト前位置2006/12/09Hara
		  endif
c			ts = dble(t1)-dble(log(rnd())*pgm(kv,ken,ka))
			ts = dble(t1)-dble(log(rnd())*pgm(kv,ken,kpart)) !07/8/4 不純物散乱
c=============エラー処理========================================
			if(sngl(ts).lt.t1)then
				write( *,*)'τ論理エラー3(emcd2-障壁後)'
				write(99,*)'τ論理エラー3(emcd2-障壁後)'
				kv=0
	exit loop
			endif
c===============================================================
c
c	----いずれでもない場合（エラー）----
		else
			write( *,*)'iflagエラー(emcd2)',iflag
			write(99,*)'iflagエラー(emcd2)',iflag
			kv=0
	exit loop
		endif
c#########################################################################
c	ドリフト運動前後の位置x,xxから領域ABを通過する粒子をカウント  120817sato

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

c	もしTelやEfが変な値だったらRejectionに入らない
	if(ict.gt.strej)then		!ict > strej=-19998
		if((rejcnt.gt.0).and.(kv.eq.1))then			!竹岸 変更 rejcnt=20,Γ谷

			if((ix.gt.1).and.(ix.lt.(nx-1)).and.	! 1< ix <309
c
cccc!Rejectionチャネル層適用!!注意!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     &	(lhet(nchannel1)+1.le.iz).and.
     &	(iz.le.lhet(nchannel2))) then !channel Rejection 11/05/25原
c     &	(iz.ge.lhet(3)+1).and.(iz.le.lhet(4)-1))then !チャネル4層目固定
c     &			(iz.ge.38).and.(iz.le.46))then		!竹岸 チャネル層のみRejection
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
	if((avsumtel(ix,iz,kv).ge.250).and.				!適用電子温度範囲
     &(avsumtel(ix,iz,kv).le.3000))then

c	-----08/5/15 竹岸-----1度はずしておきます。
c			if((iak.eq.jak).and.(pflag.lt.rejcnt))then		!前回と同じエネルギーなら2001へ
c				pflag = pflag + 1
c				go to 2001
c			else
c	----------------------	

			call rejection(
     &						hhm,hm,af2,af4,ec,p,kp,iarea,bkq,
     &						efermi,ix,iz,
     &						akx,aky,akz,kv,kl,rflag,
     &						adkx,adky,adkz,iak,avsumtel,epA,u,eg)	!竹岸変更

				if((rflag.eq.1).and.(pflag.le.rejcnt)) then
					pflag = pflag + 1
c
ccccccリジェクション・衝突電離 重複回避 粒子・回数補正11/08原	
					if(iiflag.eq.1) then		!!!
					  jpnum = jpnum -1		!!!
					  nava(iiix,iiiz,kv) = nava(iiix,iiiz,kv)-1
				  	  nava(iiix,iiiz,0) = nava(iiix,iiiz,0)-1						
					  cnava(iiix,kv) = cnava(iiix,kv)-1
					  cnava(iiix,0) = cnava(iiix,0)-1
cccccc
ccccccリジェクション・衝突電離　重複回避 確認用 11/08原
c					  write(*,*) 'emcd',n,ix,iz,jpnum,pflag,ict
c					  write(*,*) 'emcd',nava(iiix,iiiz,0),cnava(iiix,0)
c					  pause
cccccc
					endif						!!!  
cccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc120817sato

					if(scatpoint.ge.1)then		!散乱が起きた時（空散乱(scatpoint=0)にならなかった時）
					if((dhet(nchannel1).lt.z).and.(z.le.dhet(nchannel2)))then	!channel内に限定

					pass = pass - pass_r		!120817sato rejection回数分減らす
					pass_r = 0.0				!120817sato	初期化
					x_start_1 = x_start_r		!前の位置に戻す
					z_start_1 = z_start_r		!前の位置に戻す
					x_mean_free_path_sum = x_mean_free_path_sum 
     &										-abs(xx_goal-x_start_1)
					mean_free_path_sum = mean_free_path_sum
     &					-sqrt((xx_goal-x_start_1)**2+(zz_goal-z_start_1)**2)

					x_mean_free_path_count = x_mean_free_path_count -1			!カウンタ
					
			mean_free_path_sum2(xcenter,1)
     &			= mean_free_path_sum2(xcenter,1) -1.0							!カウンタ
			mean_free_path_sum2(xcenter,2)
     &			= mean_free_path_sum2(xcenter,2)-abs(xx_goal-x_start_1)
			mean_free_path_sum2(xcenter,3)
     &			= mean_free_path_sum2(xcenter,3)
     &			-sqrt((xx_goal-x_start_1)**2+(zz_goal-z_start_1)**2)

					endif	 !channel内に限定
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
c	---(チャネル内散乱頻度)---	!10/05/07
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
		kp(1,n)	= kv	!谷No.
		kp(2,n)	= ken	!エネルギーテーブルNo.
		kp(3,n)	= kl	!層No．
c		kp(n)	= ka *(nvalley*nenergy)		!ka:素材No．
c     &			+ ken* nvalley				!ken:エネルギーテーブルNo.
c     &			+ kv						!kv:谷No.

		balis_flag(n) = balis_flag2		!バリスティックの計算 07/11/22 川端

		scatpoint = 0	!次の粒子になる前に初期化120817sato
		pass_r = 0.0	!120817sato 
	
		x_start(n) = x_start_1
		z_start(n) = z_start_1

	enddo	!1~nendまでの粒子繰り返し
c
c	----(衝突電離)----
	nstat = n
	nend = jpnum
	if(n.lt.(jpnum+1))goto 1000
c
	p(4,1:jpnum)=p(4,1:jpnum)-dt
c
c-----縮退効果 ドリフト波数成分の各メッシュ・谷における粒子平均-----
c	現在はΓ谷だけ考慮 iv=1のみ
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