	subroutine drift(am_aniso,aff_aniso,fx,fz,dz,dhet,dt,de,pgm,af2,
     &						af4,hm,hhm,
     &				        akx,aky,akz,ts,x,z,t1,tau,
     &						kv,ken,kl,kl2,iflag,ie,ei,sk,
     &						kpart,ka,n)									!07/8/4 不純物散乱


	implicit none
c	
c	終了時の状態を表すフラグとして、iflagを設定する。
c	iflag=0 : ドリフト後、散乱イベント発生
c	iflag=1 : ドリフト後、境界面に至ったので進入判定(subroutine border)
c	iflag=2	: ドリフト後、dtに至ったので次のループ（粒子）に入る
c
c	基本的に###で囲んだ部分はエラー処理であり、実際には起きない条件が
c	発生した場合にのみ動作するので、普段は無視してかまわない
c
	include 'arraysize.fi'
	real fx,fz
	real dhet(0:nlayer),dt
	real de(nenergy)
	real pgm(nvalley,nenergy)
	real af2,af4,hm,hhm
c	real,	dimension (nvalley,narea)	:: smh,hhm,hm,af,af2,af4
	real akx,aky,akz,x,z,t1,tau
	real(8)	ts
	integer(1) kv, ken, kl, kl2,ka
	integer iflag
	integer kpart  !07/8/4 不純物散乱
	real, dimension	(nvalley,narea,nvalley)::am_aniso							 !20100624
	real, dimension	(nvalley,narea,nvalley)::aff_aniso							!20100624
c
	real(8)	dkx,dkz
	real	xx,zz,err,dh
c	real(8)	skx,sky
	real(8)	sk,ei
	real	ddt
	real	tauz
	integer	ie
c	real	den
	integer i !,il
	integer n
	real	dz
c
	real(8)	qht,sq,cp,sq21,sq22
	real qh !qh = q/h
	parameter(qh = 1.51925E+15)
	integer klc

c-----------------------------------------------------------------------------------------

c
c



c	--- 何らかの原因でkv(谷No.)=0である場合
c		ドリフト計算は行わない ---
	if(kv.eq.0) return
c
c#####τの判定##################################################################
c---	!粒子が界面上にいて 進入イベントがすぐさま起きる場合---
c	!iflag=1 として ドリフトは行わない
	if((dhet(kl-1).lt.z).and.(z.lt.dhet(kl)))then
		continue	!通常動作
	elseif(z.eq.dhet(kl-1))then
		if(akz.lt.0)then
			kl2 = kl-1
			iflag = 1
			return
		elseif(akz.eq.0)then
			fz = 0.0	!特殊ケース（例外処理）
			continue
		endif
	elseif(z.eq.dhet(kl))then
		if(akz.gt.0)then
			kl2 = kl+1
			iflag = 1
			return
		elseif(akz.eq.0)then
			fz = 0.0	!特殊ケース (例外処理)
			continue
		endif
c	---	粒子パラメータと粒子位置が一致しない場合(エラー回避)
c	リセス領域にて多数発生する特殊エラーを以下で回避 11/05/26 原
	elseif((kl.eq.1).and.(kpart.eq.2))then
	klc = 0
	do klc=1,3
		if((dhet(klc-1).lt.z).and.(z.lt.dhet(klc))) then
			kl = klc
			exit
		 endif
	enddo
c
c	---	粒子パラメータと粒子位置が一致しない場合(エラー処理) 
	else
		write(* ,*)'kl不正(drift)'
		write(99,*)'kl不正(drift)'
	write(*,*) n,x,z,kl,kpart
	write(*,*) n,dhet(kl-1), dhet(kl)
!SMP$	ASSERT(ITERCNT(1))
c		do il = 1, nlayer
c			if((dhet(il-1).lt.z).and.(z.lt.dhet(il)))then
c				kl = il
c			elseif(z.eq.dhet(il))then
c				kl = il
c				fz = 0.0	!特殊ケース（例外処理）
c			endif
c		enddo
		kv=0
		return
	endif
c--------------------------------------------------------------
c
c	--- 何らかの原因でt1>dtである場合
c		ドリフト計算は行わないで次のループ(iflag=2) ---
	ddt = dt*(1.0+epsilon(dt))	!誤差範囲
	if(t1.lt.dt)then
		continue	!正常→何もせず続く
	elseif((t1.ge.dt).and.(t1.lt.ddt))then
		iflag = 2
		return
	elseif(t1.ge.ddt)then
		iflag = 2
		!write( *,*)'t1>dt',t1,dt,ddt
		!write(99,*)'t1>dt',t1,dt,ddt
		return
	endif
c
1000	continue
c
c	--- τの計算 ---
	if(sngl(ts).gt.dt)then
c	--- この単位時間内で散乱イベントが起きない場合
c		τ=dt-t1だけドリフトし、iflag=2 ---
		tau = dt-t1
		iflag = 2
	else
		tau= sngl(ts-dble(t1))
c	--- 何らかの原因でτ=0である場合
c		ドリフト計算は行わないで散乱(iflag=0) ---
		if(tau.gt.0.0)then
			continue	!正常→何もせず続く
		elseif(tau.eq.0.0)then
			iflag = 0
			write( *,*)'τ=0'
			write(99,*)'τ=0'
			return
c	---	tauが0より小さい場合(エラー)   ---
		elseif(tau.le.0.0)then
			write(* ,*)'τエラー(drift)'
			write(99,*)'τエラー(drift)'
			write(99,*) tau,t1,dt,sngl(ts)
			tau = 0.0
			iflag = 0
			return
		endif
	endif
c###############################################################################
c
c
c#####ドリフト運動##############################################################
c	--- 自由行程時間τだけドリフト計算 ---
	qht=qh*tau
c
	dkx=-qht*fx		!dkx=Δakx: akxの、時間τでの増加量
	dkz=-qht*fz		!dkz=Δakz: akzの、時間τでの増加量
c
	if(af4.eq.0.0)then
		sq=1.0
	else
c		sq=1.+af4*hhm*(akx*akx+aky*aky+akz*akz)
c	sq21 = sq
		sq=1.+af4*hhm*
     &	((akx*akx)/am_aniso(kv,ka,1)+(aky*aky)/
     &	 am_aniso(kv,ka,2)+(akz*akz)/am_aniso(kv,ka,3))
     &	 *(am_aniso(kv,ka,1)+am_aniso(kv,ka,2)+
     &	 am_aniso(kv,ka,3))/3.0d0
c	akx
c	sq22 = sq
c	write( *,*) sq21,sq
c		sq=1.+af4*hhm*
c     &	((akx**2.0d0)/am_aniso(kv,ka,1)+(aky**2.0d0)/
c     &	 am_aniso(kv,ka,2)+(akz**2.0d0)/am_aniso(kv,ka,3))
c     &	 *(am_aniso(kv,ka,1)+am_aniso(kv,ka,2)+
c     &	 am_aniso(kv,ka,3))/3.0d0
		if(sq.gt.0.0)then
			sq=1.0/sqrt(sq)
		else
			write(99,*)'sqの値が不正です(drifttau)'
			stop
		endif
	endif
	cp=hm*tau
	xx = x + cp*(akx+0.5*dkx)*sq
	zz = z + cp*(akz+0.5*dkz)*sq
c###############################################################################
c	check_point
c################参照エネルギーテーブルはみ出しチェック##########################
	sk	= (akx+dkx)*(akx+dkx)+aky*aky+(akz+dkz)*(akz+dkz)
	if((sk.gt.0.0).and.(sk.le.huge(sk)))then
		if(af2.ne.0.0)then
			ei=(sqrt(1.0+af4*hhm*sk)-1.0)/af2		!eV
		else
			ei=hhm*sk
		endif
	else
		write(99,*)'skの値が不正です(drift)',sk,akx+dkx,aky,akz+dkz
		write( *,*)'skの値が不正です(drift)',sk,akx+dkx,aky,akz+dkz
	endif
	ie=int(ei/de(ken))+1
c
c=============エラー処理========================================
	if((ei.lt.0.0).and.(ei.ge.huge(ei)))then
		write( *,*)'energy計算エラー(emcd2)',ei
		write(99,*)'energy計算エラー(emcd2)',ei
		kv=0
		return
	endif
c===============================================================
c
	if(ie.ge.nemax)then
		if(ken.lt.nenergy)then
c	----	再帰処理 ->ドリフトからやり直し---------------
			ts=(ts-dble(t1))/dble(pgm(kv,ken))
     &				*dble(pgm(kv,ken+1))+dble(t1)
			ken = ken+1
			iflag = 0
			goto 1000
c	--------------------------------------------------
		else
			write( *,*)'energyが大きすぎます',ie  !,n
			write(99,*)'energyが大きすぎます',ie  !,n
			ie = nemax
		endif
	endif


c#########################################################################
c	ドリフト運動前後の位置x,xxから領域ABを通過する粒子をカウント  120817sato


c#########################################################################
c	--- ドリフト後位置が境界面をまたいでいなければ確定・終了 ---
	if((dhet(kl-1).lt.zz).and.(zz.lt.dhet(kl)))then
		kl2 = kl
		akx = akx+dkx
		akz = akz+dkz
		x = xx
		z = zz
		t1 = t1 + tau
		return
c	--- ドリフト後が境界面をまたいでいたらその境界面までドリフト ---
	elseif(zz.lt.dhet(kl-1))then
		kl2 = kl-1
		dh=dhet(kl-1)		!上に貫通
		continue
	elseif(zz.gt.dhet(kl))then
		kl2 = kl+1
		dh=dhet(kl)		!下に貫通
		continue
c	--- ドリフト後にちょうど界面に接触
c		ドリフト再計算は行わないで界面イベント(iflag=1) ---
	elseif(zz.eq.dhet(kl-1))then
		kl2 = kl-1
		iflag = 1
		akx = akx+dkx
		akz = akz+dkz
		x = xx
		z = zz
		t1 = t1 + tau
		return
	elseif(zz.eq.dhet(kl))then
		kl2 = kl+1
		iflag = 1
		akx = akx+dkx
		akz = akz+dkz
		x = xx
		z = zz
		t1 = t1 + tau
		return
	endif
c
c	--------------------------------------
	err = dz*0.000005	!1.6e-14
c	call calctau(akx,aky,akz,z,af4,hm,hhm,fz,dh,tauz)
	call calctau(am_aniso,kv,ka,akx,aky,akz,z,af4,hm,
     &				hhm,fz,dh,tauz)

	do i=1,4
		if((tauz.gt.0.0).and.(tauz.le.tau))then
c		---	tauz正常時 ---
			tau = tauz
			ts = dble(t1) + dble(tauz)
c#####ドリフト運動##############################################################
c	--- 自由行程時間τだけドリフト計算 ---
			qht=qh*tau
c
			dkx=-qht*fx		!dkx=Δakx: akxの、時間τでの増加量
			dkz=-qht*fz		!dkz=Δakz: akzの、時間τでの増加量
c
			if(af4.eq.0.0)then
				sq=1.0
			else
				sq=1.+af4*hhm*(akx*akx+aky*aky+akz*akz)
	
				sq=1.+af4*hhm*
     &	((akx**2.0d0)/am_aniso(kv,ka,1)+(aky**2.0d0)/
     &	 am_aniso(kv,ka,2)+(akz**2.0d0)/am_aniso(kv,ka,3))
     &	 *(am_aniso(kv,ka,1)+am_aniso(kv,ka,2)+
     &	 am_aniso(kv,ka,3))/3.0d0

				if(sq.gt.0.0)then
					sq=1.0/sqrt(sq)
				else
					write(99,*)'sqの値が不正です(drifttau)'
					stop
				endif
			endif
			cp=hm*tau
			xx = x + cp*(akx+0.5*dkx)*sq
			zz = z + cp*(akz+0.5*dkz)*sq
c###############################################################################
c
			if((dh-err.le.zz).and.(zz.le.dh+err))then
c				---	zzの位置が誤差範囲内(zz±(0..err) = dh) ---
				akx = akx+dkx
				akz = akz+dkz
				x = xx
				z = dh
				t1 = t1 + tau
				iflag = 1
	return
			elseif((dhet(kl-1).lt.zz).and.(zz.lt.dhet(kl)))then
c			elseif(kl.eq.kl2)then
c				---	zzの位置が進入前 ---
				akx=akx+dkx
				akz=akz+dkz
				x = xx
				z = zz
				t1 = t1 + tau
				call calctau(am_aniso,kv,ka,akx,aky,akz,z,af4,hm,
     &				hhm,fz,dh,tauz)
		cycle
			else
c				---	zzの位置が進入している ---
				tauz = tauz*0.9
		cycle
			endif
		elseif((tauz.eq.0.0).and.(z.eq.dh))then
c		---	tauz異常時 ---
			write(*,*)'tauz=0,z=dh'
	return
		else
c******************** エラー処理 **********************************
			write( *,*)'tauzエラー(drift)' !,sq21,sq22,tauz,tau
			write(99,*)'tauzエラー(drift)'
			write(99,*) tauz,tau,t1,fz
			write(99,'(4(I3))') iflag,kl,kv,i
			write(99,*) akx,aky,akz,sngl(ts),x,z,xx,zz
			if(abs(z-dh).lt.abs(zz-dh))then
				tau = 0.0
				z = dh
				iflag = 1
	return
			else
				akx=akx+dkx
				akz=akz+dkz
				x = xx
				z = dh
				t1 = t1 + tau
				iflag = 1
	return
			endif
		endif
	enddo
c		---	ループ中に収束しきれなかった場合 ---
	if(zz.ne.dh)then
		write(* ,*)'収束しません(drift)'
		write(99,*)'収束しません(drift)'
		write(99,*) tau,t1,fz
		write(99,'(4(I3))') iflag,kl,kv,i
		write(99,*) akx,aky,akz,sngl(ts),x,z,xx,zz
		z = dh
		iflag = 1
	return
	endif
c
	write( *,*)'論理エラー(drift)'
	write(99,*)'論理エラー(drift)'
c	**********************************************************
	stop
	end subroutine drift
c
c
c
c
      subroutine calctau(am_aniso,kv,ka,akx,aky,akz,z,af4,
     &	hm,hhm,fz,dh,tauz)
c	粒子が境界面dhに接触する時間tauzを計算する。
c	hm*akz*t - Fz*(q/2m*)*(t**2) -(dh-z)sq = 0 
c	の方程式の解 tauz1,tauz2(tauz1<tauz2)を導き、
c	tauz1>0のときtauz=tauz1、tauz1<0のときtauz=tauz1とする。
c	重解ではtauz1=tauz2=tauz、
c	解なしの場合tauz=∞となる。
c
	implicit none
	

	real, dimension	(3,5,3)::am_aniso							 !20100624
	integer(1) kv,ka															!20100624

c	real	z
	real	akx,aky,akz,z
c	integer	iflag
	real	af4,hm,hhm
	real	fz
	real	dh
	real	tauz
c
	real	tauz1,tauz2
	real(8)	akz2
	real	akk
	real(8)	sq,sq2
	real(8)	hqfz
	real(8)	b24ac
	real(8)	ac4
c==========
	real	hq
	parameter(hq = 6.58218E-16)
c
	akz2=akz*akz
	akk =akx*akx+aky*aky+akz2
	akk = ((akx**2.0d0)/am_aniso(kv,ka,1)+(aky**2.0d0)/
     & am_aniso(kv,ka,2)+(akz**2.0d0)/am_aniso(kv,ka,3))
     & *(am_aniso(kv,ka,1)+am_aniso(kv,ka,2)+
     & am_aniso(kv,ka,3))/3.0d0

	if(af4.eq.0.0)then
		sq=1.0
	else
		sq	=1.+af4*hhm*akk
		if(sq.ge.0.0)then
			sq=sqrt(sq)
		else
			write(99,*)'sqの値が不正です(calctau)'
			stop
		endif
	endif
c
	if(fz.eq.0.0)then
c		τ二乗の項が消えるので、線形な式で解く
		tauz = (dh-z)/akz/hm*sq
	else
		hqfz = hq/fz
		ac4  = dble(dh-z)*fz/hhm*sq
		if(dh.eq.z)then
			tauz1 = min(0.0,sngl(hqfz*(2.0*akz)))
			tauz2 = max(0.0,sngl(hqfz*(2.0*akz)))
			if(tauz1.gt.0.0) tauz = tauz1
			if(tauz1.le.0.0) tauz = tauz2
		elseif(abs(ac4*8e5).lt.abs(akz2))then
c			fzが小さいと誤差が大きくなるので線形近似で解く
			tauz = (dh-z)/akz/hm*sq
		else
c			方程式 Fz*(q/2m*)*(t**2)-hm*akz*t+(dh-z)sq = 0 を解く
c			t = qFz/h*(akz±sqrt(akz**2-(dh-z)*Fz/hhm*sq)
			b24ac = (akz2-ac4)
			if(b24ac.gt.0.0)then
c				方程式の解 tauz1,tauz2
				sq2   = sqrt(b24ac)
				tauz1 = min(hqfz*(akz-sq2),hqfz*(akz+sq2))
				tauz2 = max(hqfz*(akz-sq2),hqfz*(akz+sq2))
				if(tauz1.gt.0.0) tauz = tauz1
				if(tauz1.le.0.0) tauz = tauz2
			elseif(b24ac.eq.0.0)then
c				方程式の解 → 重解
				tauz = hqfz*(akz)
			else
c				方程式の解なし tauz→∞
				tauz = huge(tauz)
			endif
		endif
	endif
c
	return
	end
