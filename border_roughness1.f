c
c-----粒子が境界面で反射または進入するイベント-----
	subroutine border_roughness1(
     &					hhm,af,af2,af4,eg,
     &					akx,aky,akz,kv,kl,kl2,ka,iarea,lhet,iz,	!2006/12/09 Hara
     &					swk_rou,ken,ix,de,					!121029sato
     &					roughness1_countx,roughness1_counte)
	implicit none
c
	include 'arraysize.fi'
	common	/arraydata/nx,nz					!121029sato
	integer	nx,nz
	real,	dimension (nvalley,narea)	:: hhm,af,af2,af4,eg
	real	akx,aky,akz
	integer(1)	kv,kl,kl2,ka,iz	!2006/12/09 Hara
	integer(1)	iarea(nlayer)
	integer(2)	lhet(nlayer)	!2006/12/09 Hara
c
c	---	ローカル変数	---
c	skx,sky,skz:各方向波数の二乗を格納する変数
c	vb:障壁高さを格納する変数
c	ex,ey,ez:各方向運動エネルギーを格納する変数
	real(8)	skx,sky,skz
	real	vb
	real(8)	ex,ey,ez
	integer(1)	ka2

c------ラフネス散乱--R.P. Joshiモデル(2012年秋)-------------------------!121029sato
	real	swk_rou(nemax,nvalley,nenergy,narea)		!ラフネス散乱レート
	real	swkd_rou				!					!ラフネス散乱レート線形補間
	integer,dimension (0:nx,narea,2) :: roughness1_countx
	integer,dimension (nemax,narea,2) :: roughness1_counte
	real	de(nenergy)
	real	den_2deg
	real(8)	ei_2deg
	real(8)	sk_2deg
	integer	ie_2deg
	integer(1)	ken
	integer		ix
	real	cs,sn
	real	rnd		!乱数
c

c=======================================================================================!121029sato
	if((iz.eq.lhet(nchannel1)).or.(iz.eq.lhet(nchannel2))) then		!チャネルのヘテロ界面にあるとき
c	open(unit=1299,file='zztriela.txt')	!確認用!

c------ラフネス散乱--2DEGエネルギー-------------------------

	sk_2deg	= akx*akx+aky*aky			!2DEGの運動量からエネルギーを求める
	if((sk_2deg.ge.0.0).and.(sk_2deg.le.huge(sk_2deg)))then
		if(af2(kv,ka).ne.0.0)then
			ei_2deg=(sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*sk_2deg)-1.0)/af2(kv,ka)		!eV
		else
			ei_2deg=hhm(kv,ka)*sk_2deg
		endif
	else
		write(99,*)'sk_2degの値が不正です',sk_2deg,akx,aky
		write( *,*)'sk_2degの値が不正です',sk_2deg,akx,aky
	endif
									!ei_2deg	!2degエネルギー実数 
	ie_2deg = int(ei_2deg/de(ken))+1			!2degエネルギー整数化(離散値) 
	den_2deg=ei_2deg/de(ken)-float(ie_2deg-1)	!de_2deg:ie_2degとei_2degの差（ズレ）

c------ラフネス散乱--散乱レートの線形補間-------------------------
						
	if(ie_2deg.eq.1) then						
		swkd_rou=swk_rou(ie_2deg,kv,ken,ka)	
											
	else									
		swkd_rou=den_2deg*swk_rou(ie_2deg,kv,ken,ka)	+
     &					(1-den_2deg)*swk_rou(ie_2deg-1,kv,ken,ka)	
													
		if(swkd_rou.gt.1.0)then		!エラー処理	
			write( *,*)'swkd_rouが1.0を超えました',swkd_rou
			write(99,*)'swkd_rouが1.0を超えました',swkd_rou
			swkd_rou = swkd_rou/swkd_rou	!再規格化 swkd_rou=1					
		endif													
	endif
c------ラフネス散乱--散乱レート参照-------------------------

	if(rnd().lt.swkd_rou)then		!ラフネス散乱レート満たすとき
		cs  = 1.0-2.0*rnd()
		if(abs(cs).ge.1.0)then
		  sn = 0.0
		else
		  sn  = sqrt(1.0-cs*cs)
		endif

		akx = sqrt(sk_2deg)*cs		!xy平面等方性
		aky = sqrt(sk_2deg)*sn		!xy平面等方性		
		
c------ラフネス散乱--カウンタ-------------------------
		
		if(iz.eq.lhet(nchannel1)) then		!チャネル　上ヘテロ界面
			roughness1_countx(ix,ka,1) =roughness1_countx(ix,ka,1) +1
			roughness1_counte(ie_2deg,ka,1) 
     &							  =roughness1_counte(ie_2deg,ka,1) +1
		endif	
		if(iz.eq.lhet(nchannel2)) then		!チャネル　下ヘテロ界面
			roughness1_countx(ix,ka,2) =roughness1_countx(ix,ka,2) +1
			roughness1_counte(ie_2deg,ka,2) 
     &							  =roughness1_counte(ie_2deg,ka,2) +1
		endif


	endif		!ラフネス散乱レート満たすとき

	endif		!チャネルのヘテロ界面にあるとき

c=========================================================================================!121029sato
c
c	--- z方向の運動エネルギー計算 ---
	if(af4(kv,ka).ne.0.0)then
		skz = dble(akz)*dble(akz)
		ez=(sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*skz)-1.0)/af2(kv,ka)		!eV
	else
		ez=hhm(kv,ka)*akz*akz
	endif
c	---------------------------------
c
	ka2 = iarea(kl2)
c	vb =0.7*(eg(kv,ka2)-eg(kv,ka))		!Vb:障壁高さ[eV]
c	if(iz.lt.(nlayer-2)-5) then	!ドリフト前位置がclassical領域のとき(2006/08/08改良)!2006/12/09 Hara
c	if(iz.lt.(nlayer-4)-4) then	!ドリフト前位置がclassical領域のとき(2006/08/08改良)!2006/12/22 Hara
c		vb =0.7*(eg(kv,ka2)-eg(kv,ka))		!Vb:障壁高さ[eV]
c	else						!ドリフト前位置がEP領域のとき
		vb = 0.0			!なぜ 0 eV ?? 2011/03/23原
c	endif
c
	if(ez.gt.vb)then					!障壁Vbより粒子のエネルギーが 高いか
c	===	障壁を乗り越える場合（進入）	===
		ez = ez-vb			!z方向のエネルギーから障壁エネルギーをひく
		skx = dble(akx)*dble(akx)
		sky = dble(aky)*dble(aky)
		if(af4(kv,ka).ne.0.0)then
			ex=(sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*skx)-1.0)/af2(kv,ka)		!eV
			ey=(sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*sky)-1.0)/af2(kv,ka)		!eV
		else
			ex=hhm(kv,ka)*skx		!αが０の時は単純計算
			ey=hhm(kv,ka)*sky
		endif
c
		kl = kl2
		ka = ka2
c		---	求めたエネルギーから波数決定（α＝０でも計算可能）
		akx=sign(sngl(sqrt(ex*(1.0+af(kv,ka2)*ex)/hhm(kv,ka2))),akx)
		aky=sign(sngl(sqrt(ey*(1.0+af(kv,ka2)*ey)/hhm(kv,ka2))),aky)
		akz=sign(sngl(sqrt(ez*(1.0+af(kv,ka2)*ez)/hhm(kv,ka2))),akz)
c		---	sign(a,b) = aの絶対値にbの符号を掛けた値
c
	else
c	===	障壁を乗り越えない場合（反射）	===
		akz=-akz
	endif
c
      return
	end