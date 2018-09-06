	subroutine condition(dt,np1,bias,dvi,nig,nic,de,jpot)

	include 'arraysize.fi'
	common /step/jinit,jstat,jstep,jdisp
	integer	jinit,jstat,jstep,jdisp

	real	dt
	integer	np1
	real,dimension(npole):: bias,dvi
	integer nig,nic
	real	de(nenergy)
	integer	jpot

	real	dem(nenergy)
	real	dinit,dstat,dstep
	integer i
	character(80) form

c	---(シミュレーション条件)---
	read(8,*) dt		!2.0e-15	!１ステップの実時間
	read(8,*) dinit		!50.0e-12	!初期化に要する時間ステップ
	read(8,*) dstat		!0固定		!観測スタート時間ステップ
	read(8,*) dstep		!200.0e-12	!状態あたりのステップ数
	read(8,*) jpot		!1			!ポテンシャル計算比
	read(8,*) np1		!400		!１メッシュあたりの粒子数
	read(8,*) (bias(i), i = 1, 3)	!初期値sgd
	read(8,*) (dvi(i), i = 1, 3)	!増加量
	read(8,*) nig		!5
	read(8,*) nic		!5
	read(8,*) (dem(i), i = 1, nenergy)	!demax

	de	= dem/float(nemax)	!全体配列

	jinit = nint(dinit/dt)
	jstat = nint(dstat/dt)
	jstep = nint(dstep/dt)
	jpot  = max(1,jpot)

	return
	end
