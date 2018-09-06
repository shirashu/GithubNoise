c################################################################
c		電流揺らぎのデバイス内分布算出プログラム
c
c 
c 2014/12/20,高橋　　（筑波大の佐野さんの電流揺らぎ算出方法）
c################################################################
c
      subroutine current_fluctuation(jpnum,dx,dz,xmax,zmax,
     &fix_u,ec,cur,dt,cxpole1,cxpole2,
     & p,kp,hhm,hm,af4,af2,iarea,adkx,adky,adkz,spnum,cff,bscat_xflag)																
c===引数===
c---変数配列用パラメータ---
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c---基本パラメータ---
	real(8) q,h,bk,pi
	real spnum,dt
	parameter(h  = 1.05459e-34, q  = 1.60219e-19,		!hはh(ber)であり単位は[Js]
     &		 bk = 1.38066e-23, pi = 3.1415927)			!bkの単位は[JK-1]
	real cur(3)
	real di_node,dis_node
c
	integer	jpnum
	integer cff,bscat_xflag,fix_u
c---デバイス構造---
	real	dx,dz,xmax,zmax
	integer(1)	iarea(nlayer)
c---電極---
	real	cxpole1(npole),cxpole2(npole)
c---領域別パラメータ---
	real,	dimension (nvalley,narea)	:: hm,hhm,af2,af4	!hhm=h(ber)^2/(2mq)
	real,	dimension (nvalley,narea) :: ec					!hm = h(ber)/m*
c---粒子状態---
	real,	dimension   (6,npmax)	:: p
	integer(1),dimension (3,npmax)	:: kp
c---kdに関するパラメータ---	
	real adkx(0:nx,0:nz,0:nvalley)
	real adky(0:nx,0:nz,0:nvalley)
	real adkz(0:nx,0:nz,0:nvalley)
	real akx(0:nx,0:nvalley)
	real aky(0:nx,0:nvalley)
	real akz(0:nx,0:nvalley)
c
	real tkx,tky,tkz
c--- シミュレーション条件 ---
c	common /outp/jouti,joutp2,joutp,jheat,jout2d,joute,joute2,jcput,sw
c	integer	jouti,joutp2,joutp,jheat,jout2d,joute,joute2,jcput,sw(9)
c	integer	hcss
c	
c--- ローカル変数 ---						
	integer(8),save :: jt				!doループカウンタ
	integer,save :: idx, ixmax			!雑音解析用メッシュ数
	integer ix,iz,iix	 						!粒子の存在するメッシュ
	real(8),save :: sv,sv_v			!メッシュ内の粒子数
	allocatable sv(:),sv_v(:,:)			
	real(8),save :: sn,sn_v,ssn,ssn_v					!谷別メッシュ内粒子数
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
	real(8),save :: tel_v,stel_v,stel					!電子温度
	allocatable tel_v(:,:),stel_v(:,:),stel(:)
	real(8),save :: Sv_1,Si_1,Sv_2,Si_2,dv0,dv00						!Gonzaletzの計算で用いる変数
	allocatable Sv_1(:),Si_1(:),Sv_2(:),Si_2(:),dv0(:),dv00(:)
	real(8),save :: var_x,atvth							!分布関数の分散
	allocatable var_x(:),atvth(:)
c	real(8),save :: rv1,rv2,rv3,rv4,rv5
c	allocatable rv1(:),rv2(:),rv3(:),rv4(:),rv5(:)
c	real(8),save :: rv1_v,rv2_v,rv3_v,rv4_v,rv5_v
c	allocatable rv1_v(:,:),rv2_v(:,:),rv3_v(:,:),rv4_v(:,:),rv5_v(:,:)
	real(8) vth,ei,ei1,sid,sid2
c	
	integer(8),save :: k,m					!時間平均のカウンタ
c
	integer kv,ka,kl					
	real sk,sk1,sk2,gk,v1,L1,L2								!v1：群速度（x成分）L1,2:Lgeffの範囲選択
	real(8),save :: a,b,c,vave_all,n_all,w,freq,var			!vave_allはチャネル領域全体のドリフト速度
c	
c-------------------メッシュの設定------------------------------
	if(.not. allocated(atv1))then
		jt = 0
		idx = 1									!まとめるメッシュの数  
		ixmax = ifix((xmax-float(idx-1)*dx)/(float(idx)*dx))
		k = 0	
		m = 0	
		a = (q*spnum / xmax)**2					!spnum=75000
c---------------simpson method の積分範囲の設定-----------------
c		b = 0.0
c		c = ixmax*dx
		di_node = 0.0
		dis_node = 0.0		
c--------------Gonzaletzの方法による計算における設定------------
		freq = 1.0e11 							!100GHz
c---------------------------------------------------------------				
c		allocate(rv1(-1:nx),rv1_v(-1:nx,nvalley))
		allocate(sv(-1:nx),sv_v(-1:nx,nvalley))
		allocate(sn(-1:nx),sn_v(-1:nx,nvalley))
		allocate(atv1(-1:nx),atv1_v(-1:nx,nvalley))
		allocate(atn(-1:nx),atn_v(-1:nx,nvalley))
c----------------------------ゼロクリア-------------------------------------
c	    rv1 = 0.0		!計算誤差保障(subroutine sigma)
c	    rv1_v = 0.0		!バグあり（情報落ちor他のサブルーチンとの多重定義）
	    sv = 0.0		!ドリフト速度の合計値（全ステップ）
	    sv_v = 0.0
	    sn = 0.0		!粒子数合計値
	    sn_v = 0.0   
	    atv1 = 0.0		!ドリフト速度の時間平均
	    atv1_v = 0.0    
		atn = 0.0		!粒子数の時間平均
		atn_v = 0.0
		vave_all = 0.0
		n_all = 0.0
		var = 0.0
	endif
c
c######Ramo-Shockleyのサンプリング領域の設定(実効ゲート長を空乏層幅で定義)######
c-----------------------リザーバーの平均粒子数の抽出----------------------------
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
		n_all = n_all /k						!リザーバーの平均粒子数
		k = 0
c--実効ゲート長の設定（粒子数がリザーバの平均値より10％以下（上）になったところ)--
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
c-----------Ramo-Shockleyのサンプリング領域の設定(任意に定義)----------------
c				L1 =
c				L2 =
c----------------------------------------------------------------------------
		n_all = 0.0
		atn = 0.0
	endif
c###########各メッシュごとのドリフト速度，粒子数の時間平均##################
c---------------メッシュごとのドリフト速度,粒子数の合計---------------------
	if((jt.ge.8000).and.(jt.le.9000)) then
		if(jt.eq.8000) then
			bscat_xflag = 0					!後方散乱抑制フラグ	
			fix_u = 0						!ポテンシャル固定フラグ
			if(fix_u.eq.1) write(*,*)'fixed potential'
			if(bscat_xflag.eq.1) write(*,*)'prohibite back scattering'
		endif		
	do n=1,jpnum
	kv = kp(1,n)	!谷情報
	kl = kp(3,n)	!層No．
	ka = iarea(kl)	!層No.から素材No.へ
      ix = min(nx,max(0,nint(p(5,n)/(idx*dx))))					!粒子のメッシュ位置
	sk1 = p(1,n)*p(1,n)+p(2,n)*p(2,n)+p(3,n)*p(3,n)				!縮退効果必要？
      v1 = p(1,n)*hm(kv,ka)/sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*sk1)	!群速度（ドリフト成分）
	sv(ix) = sv(ix) + v1										!各セクションごとのドリフト速度の和
		sv_v(ix,kv) = sv_v(ix,kv) + v1								!谷別
c     call sigma(v1, sv(ix), rv1(ix))								!ドリフト成分の和
c     call sigma(v1, sv_v(ix,kv), rv1_v(ix,kv))					!ドリフト成分の和(谷別）
		
      sn(ix) = sn(ix)+1.0											!メッシュ内の粒子数
	sn_v(ix,kv) = sn_v(ix,kv)+1.0									!谷別 
c---------サンプリング領域全体でのドリフト速度&キャリア数の和
	if((ix.ge.L1).and.(ix.le.L2))then					!Ramo-Shockleyのサンプリング領域
	vave_all = vave_all + v1									
	n_all =n_all + 1.0
	endif															
	end do
	di_node	= di_node + cur(3)									!端子電流の時間ステップ和
	k = k + 1
	endif
c-------------------メッシュごとでの平均ドリフト速度,粒子数(時間平均)------------------
	if((jt.ge.9001).and.(jt.le.10000)) then	
		if(jt.eq.9001)then
		do ix=1,ixmax
			if(sn(ix).ne.0)then
				atv1(ix) = sv(ix)/sn(ix)		   		!各メッシュでの速度の平均値
				atn(ix)   = sn(ix) / k	
			else
				atv1(ix) = 0.0
				atn(ix) = 0.0
			endif
		enddo
		di_node	= di_node / k							!端子電流の時間平均値
c--------------------------------------谷別--------------------------------------------
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
		n_all = n_all / k								!実効ゲート長内の平均粒子数
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
c##########################電流揺らぎΔI＾２###############################
c		熱速度と時間平均ドリフト速度の差から電流揺らぎを算出
c		(ドリフト速度の揺らぎも算出）
c##########################################################################
c
c---------------------------熱速度vthの計算--------------------------------
	do n=1,jpnum
		kv = kp(1,n)
		kl = kp(3,n)
		ka = iarea(kl)	!層No.から素材No.へ
	    ix = min(nx,max(0,nint(p(5,n)/(idx*dx))))
		iz = min(nz,max(0,nint(p(6,n)/dz)))				!メッシュはeltempと同期
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
		sk1 = p(1,n)*p(1,n)+p(2,n)*p(2,n)+p(3,n)*p(3,n)				!縮退効果必要？
c		sk2 = tkx*tkx + tky*tky + tkz*tkz
	    v1 = p(1,n)*hm(kv,ka)/sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*sk1)	!群速度（ドリフト成分）
c
	    sn(ix) = sn(ix) + 1.0									!メッシュ内の粒子数
		sn_v(ix,kv) = sn_v(ix,kv) + 1.0						!メッシュ内の粒子数(谷別) 
		ssn(ix) = ssn(ix) + 1.0								!全ステップのメッシュ別総粒子数
		ssn_v(ix,kv) = ssn_v(ix,kv) + 1.0					
c	   
		if(af4(kv,ka).ne.0.0)then									!非放物線性考慮
			sq = sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*sk)
			sq1 = sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*sk1)				!縮退考慮せず
			ei1=(sq1-1.0)/af2(kv,ka)+ec(kv,ka)
			ei=(sq-1.0)/af2(kv,ka)+ec(kv,ka)						!電子のエネルギー(Ekxのみ)
c     &			-epA(ix,iz)+u(ix,iz)-0.7*(eg(1,ka)-eg(1,2))	
		else
			ei1=hhm(kv,ka)*sk1+ec(kv,ka)
			ei=hhm(kv,ka)*sk+ec(kv,ka)						
c     &			-epA(ix,iz)+u(ix,iz)-0.7*(eg(1,ka)-eg(1,2))		
		endif
c		x成分のエネルギーから抽出した熱速度（1/2*m*v(x)^2=Ex）
		vth = sqrt((2.0*q*hm(kv,ka)*ei) / h)				!x方向成分の熱速度
		atvth(ix) =atvth(ix) + vth
c
c		トータルのエネルギーからの熱速度
c		vth = sqrt((2.0*q*hm(kv,ka)*ei1) / h)			
c		vth = sqrt(sk/sk1)*vth)				!x方向成分の熱速度(1つ上とセット)
c
c		トータルの電子温度
		stel_v(ix,kv) = stel_v(ix,kv) + (2.0/3.0)*q*ei1/bk	!電子温度の和
		stel(ix) = stel(ix) + (2.0/3.0)*q*ei1/bk			!電子温度の和
		var_x(ix) = var_x(ix) + sqrt((2.0/3.0)*q*ei1*hm(kv,ka)/h)	!分布関数の分散((KbTe/m*)^1/2)
c		（分散≒熱速度vthとみなせる）
c
c		x成分の電子温度
c		stel_v(ix,kv) = stel_v(ix,kv) + q*ei/bk			!電子温度の和
c		stel(ix) = stel(ix) + q*ei/bk					!電子温度の和

c
c---------------各メッシュごと速度の揺らぎの和(粒子に対する和------------
c		熱速度にドリフト成分が入っているとした
c		平均ドリフト速度の考え方が2通りできる（メッシュ毎）
		if(p(1,n).ge.0)then								!熱速度を正負で場合分け
		dvi(ix)   = vth  - atv1(ix)						!速度の偏差
		dvi_v(ix,kv)  = vth  - atv1_v(ix,kv)			!各谷ごと
		sv_v(ix,kv) = sv_v(ix,kv) + vth 
		else 
		dvi(ix)   = -vth  - atv1(ix)					!速度の偏差
		dvi_v(ix,kv)  = -vth - atv1_v(ix,kv)			!各谷ごと
		sv_v(ix,kv) = sv_v(ix,kv) - vth					!理論的にはatv1と等しくなる
		endif
c			（Ramo-Shockleyの定理のサンプル領域全体）
c		if(p(1,n).ge.0)then								!熱速度を正負で場合分け
c		dvi(ix)   = vth - vave_all						!速度の偏差
c		dvi_v(ix,kv)  = vth - vave_all					!各谷ごと
c		sv_v(ix,kv) = sv_v(ix,kv) + vth 
c		else 
c		dvi(ix)   = -vth - vave_all						!速度の偏差
c		dvi_v(ix,kv)  = -vth - vave_all					!各谷ごと
c		sv_v(ix,kv) = sv_v(ix,kv) - vth					!理論的にはatv1と等しくなる
c		endif
c
		dvi2(ix)  = dvi(ix)**2.0							!速度の揺らぎ
		dvi2_v(ix,kv) = dvi_v(ix,kv)**2.0
c					　（メッシュ毎）
		dv1(ix) = v1 - atv1(ix)							!ドリフト速度の揺らぎ
		dv2(ix) = dv1(ix)**2.0
		dv1_v(ix,kv) = v1 - atv1_v(ix,kv)
		dv2_v(ix,kv) = dv1_v(ix,kv)**2.0
c			（Ramo-Shockleyの定理のサンプル領域全体）
c		dv1(ix) = v1 - vave_all							!ドリフト速度の揺らぎ
c		dv2(ix) = dv1(ix)**2.0
c		dv1_v(ix,kv) = v1 - vave_all
c		dv2_v(ix,kv) = dv1_v(ix,kv)**2.0	
c		
c		熱速度考慮
		avi(ix) = avi(ix) + dvi2(ix)				!電流揺らぎ（注！これは時間に対しての和）
		avi_v(ix,kv) = avi_v(ix,kv) + dvi2_v(ix,kv)
c
c		ドリフト成分のみ
c		avi(ix) = avi(ix) + dv2(ix)					!電流揺らぎ（注！これは時間に対しての和）
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
c------------------計算誤差保障用（バグあり）--------------------------------
c		call sigma(dvi(ix), sdvi(ix), rv1(ix))					!速度の揺らぎの和
c		call sigma(dvi2(ix), sdvi2(ix), rv2(ix))
c	    call sigma(dvi_v(ix,kv), sdvi_v(ix,kv), rv1_v(ix,kv))
c	    call sigma(dvi2_v(ix,kv), sdvi2_v(ix,kv), rv2_v(ix,kv))	
c
c		call sigma(dv1(ix), sdv1(ix), rv3(ix))					!ドリフト速度の揺らぎの和
c		call sigma(dv1_v(ix,kv), sdv1_v(ix,kv), rv3_v(ix,kv))	!速度の揺らぎの和
	enddo	
c------------------各メッシュごとの粒子数の揺らぎの和（時間に対して）-----------
	do ix=1,ixmax
		dn(ix) = sn(ix) - atn(ix) 
		dn2(ix) = dn(ix)**2.0 
		dn_v(ix,1:3) = sn_v(ix,1:3) - atn_v(ix,1:3)				!各谷ごと
		dn2_v(ix,1:3) = dn_v(ix,1:3)**2.0
c
		sdn(ix) = sdn(ix) + dn(ix)
		sdn2(ix) = sdn2(ix) + dn2(ix)
		sdn_v(ix,1:3) = sdn_v(ix,1:3) + dn_v(ix,1:3)
		sdn2_v(ix,1:3) = sdn2_v(ix,1:3) + dn2_v(ix,1:3)
	enddo
c
c	Gonzaletzの方法によるSv、Sinの計算
	do ix=1,ixmax
		if(m.eq.0)then
		dv0(ix) = sdv1(ix) / sn(ix)	
		endif
		dv1(ix) = sdv1(ix) / sn(ix)								!セクションごとの速度偏差（アンサンブル平均）
		w =2.0 * pi * freq
		Sv_1(ix) = Sv_1(ix) + 4.0*dv0(ix) * dv1(ix) *cos(w*dt*m) * dt   !点で見てる
		Si_1(ix)=Si_1(ix)+a/(dx*idx)*xmax**2.0*sn(ix)*spnum*Sv_1(ix)	!点で見てる
		Sid = Sid + Si_1(ix) * dx*idx									
c		Cv ∝ n の補正による計算
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
c		call sigma(dn(ix), sdn(ix), rv4(ix))					!速度の揺らぎの和
c		call sigma(dn2(ix), sdn2(ix), rv5(ix))
c	    call sigma(dn_v(ix,kv), sdn_v(ix,kv), rv4_v(ix,kv))				
c	    call sigma(dn2_v(ix,kv), sdn2_v(ix,kv), rv5_v(ix,kv))
c
		sn = 0.0		
	    sn_v = 0.0 
c
		dis_node = dis_node + (cur(3) - di_node)**2.0				!端子電流の揺らぎの積算値
		m = m + 1
	endif
c----------------------各メッシュごと速度の揺らぎの和（時間に対して）-----------
	if(jt.eq.9002)then
	do ix=1,ixmax
		if(ssn(ix).ne.0)then
			dvi(ix) = sdvi(ix) / ssn(ix) 
			dvi2(ix) = sdvi2(ix) / ssn(ix) 
			dv1(ix) = sdv1(ix) / ssn(ix)
			dv2(ix) = sdv2(ix) / ssn(ix)
			stel(ix) = stel(ix) / ssn(ix)						!各メッシュ別電子温度分布
			var_x(ix) = var_x(ix) / ssn(ix)						!各メッシュ別の分布関数の分散
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
				tel_v(ix,1:3) = stel_v(ix,1:3) / ssn_v(ix,1:3)	!各メッシュ別電子温度分布
				atvth_v(ix,1:3) = sv_v(ix,1:3) / ssn_v(ix,1:3)	!各メッシュごとのキャリアの速度の平均
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
c---------------各メッシュごとの電流の揺らぎ（時間平均）------------------
	do ix=1,ixmax			
c
		avi(ix) = avi(ix)/m 				!電流揺らぎ
		avi_v(ix,1:3) = avi_v(ix,1:3)/m 
c
c------------------各メッシュごとの超粒子数の揺らぎ（時間平均）---------------
c
		dn(ix) = sdn(ix)/m 							
		dn2(ix) = sdn2(ix)/m 
		dn_v(ix,1:3) = sdn_v(ix,1:3)/m 			
		dn2_v(ix,1:3) = sdn2_v(ix,1:3)/m 
c--------------------各メッシュにおける電流揺らぎの算出---------------------
		di2(ix) = a * avi(ix) / spnum							!(a=q/l)^2
		di2_v(ix,1:3) = a * avi_v(ix,1:3)/ spnum				!谷別(spnumの影響非考慮)
c##################thermal,shot noiseの分離(ドリフト速度の揺らぎ)##################
	di_s(ix)  = a * dn2(ix) * (atv1(ix)**2.0)	/ spnum						!shot noise
	di_t(ix)  = a * dv2(ix) * (atn(ix)**2.0)	/ spnum						!thermal noise　?spnumの扱いがよくわからない
	di_st(ix) = abs(a * 2.0*atn(ix)*dn(ix)*atv1(ix)*dv1(ix))/spnum		!cross term
	di2_stc(ix) = di_s(ix) + di_t(ix) +	di_st(ix)						!ドリフト速度の揺らぎ
c																		(上の値と一致するかの確認)
	di_sv(ix,1:3)  = a * dn2_v(ix,1:3) * (atv1_v(ix,1:3)**2.0)/ spnum		!谷別	
	di_tv(ix,1:3)  = a * dvi2_v(ix,1:3) * (atn_v(ix,1:3)**2.0)/ spnum			
	di_stv(ix,1:3) = a * abs(2.0*atn_v(ix,1:3)*dn_v(ix,1:3)*
     &                   		atv1_v(ix,1:3)*dvi_v(ix,1:3))/ spnum	
	di2_stcv(ix,1:3) = di_sv(ix,1:3)+di_tv(ix,1:3)+di_stv(ix,1:3)
c##################thermal,shot noiseの分離(熱速度の揺らぎ)##################
c	di_s(ix)  = a * dn2(ix) * (atvth(ix)**2.0)						!shot noise
c	di_t(ix)  = a * dv2(ix) * ((atn(ix)*spnum)**2.0)					!thermal noise
c	di_st(ix) = abs(a * 2.0*atn(ix)*spnum*dn(ix)*atvth(ix)*dv1(ix))	!cross term
c	di2_stc(ix) = di_s(ix) + di_t(ix) +	di_st(ix)					!ドリフト速度の揺らぎ
c																	(上の値と一致するかの確認)
c	di_sv(ix,1:3)  = a * dn2_v(ix,1:3) * (atvth_v(ix,1:3)**2.0)			!谷別	
c	di_tv(ix,1:3)  = a * dvi2_v(ix,1:3) * ((atn_v(ix,1:3)*spnum)**2.0)			
c	di_stv(ix,1:3) = a * abs(2*atn_v(ix,1:3)*spnum*dn_v(ix,1:3)*
c     &                   		atvth_v(ix,1:3)*dvi_v(ix,1:3))	
c	di2_stcv(ix,1:3) = di_sv(ix,1:3)+di_tv(ix,1:3)+di_stv(ix,1:3)
c
c--------------------各メッシュの粒子(電子)数の揺らぎ--------------------------
     	dn(ix) = sdn(ix) * spnum							
	dn2(ix) = sdn2(ix) * spnum**2.0
	dn_v(ix,1:3) = sdn_v(ix,1:3) * spnum			
	dn2_v(ix,1:3) = sdn2_v(ix,1:3) * spnum**2.0
c
	atn(ix) = atn(ix) * spnum
c----------------------------分布関数の分散の和---------------------------------
	if((ix.ge.L1).and.(ix.le.L2))then		!Ramo-Shockleyのサンプリング領域
	var = var + var_x(ix)
	endif
	enddo	
	var = var/(L2-L1)								!平均の分布関数分散値
c-------------------端子電流の揺らぎ（時間平均）--------------------------------
		dis_node	= dis_node/m/spnum
c##############################出力関連#########################################
	do kv=1,nvalley
		do ix=1,ixmax
c
5	format('端子電流のゆらぎ（分散） = ',e15.7)
10	format((I12,E12.3))
100	format((I12,E12.3,E12.3,E12.3,E12.3,E12.3))
150	format((I12,E12.3,E12.3))
170	format((I12,E12.3,E12.3,E12.3))
200	format((I12,E12.3,E12.3,E12.3,E12.3))
300	format(4(A12))
400	format(6(A12))
500	format(5(A12))
1000	format('分布から求めたゆらぎ= ',e15.7)
1001	format('<σv> = ',e15.7)
1002	format('<vd> = ',e15.7)
1003	format('<Nch> = ',e15.7)
1004	format('端子: <δId/Id> = ',e15.7)
1005	format('中心極限近似 : <δId/Id> = ',e15.7)
1006	format('超粒子数 = ',e15.7)
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
	cff=1	!終了フラグ
	endif
c
		enddo
	enddo
c
c----------------------Sv,Sin出力------------------------
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
c---------端子電流の揺らぎ-----------
	write(907,5)dis_node
c--佐野さんの方法の電流の揺らぎの和--
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