c
c-----粒子が境界面で反射または進入するイベント-----
	subroutine border(
     &					hhm,af,af2,af4,eg,
     &					akx,aky,akz,kv,kl,kl2,ka,iarea,lhet,iz)	!2006/12/09 Hara
	implicit none
c
	include 'arraysize.fi'
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
c
c=================================================================================================
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