	subroutine outinit(dt)
c
      implicit none
	common /step/jinit,jstat,jstep,jdisp
	common /outp/joutc,jouti,jouta,jsave,jheat,jcput,joutpi,sw
c
	integer	jinit,jstat,jstep,jdisp
	integer	joutc,jouti,jouta,jsave,jheat,jcput,joutpi,sw(8)
c
	real	ddisp,doutc,douti,douta,dsave,dheat,dcput,doutpi
	real	dt
c
	integer i
c
	integer jdata
c
c k:出力スイッチ	電流, 雑音, ポテンシャル, エネルギー, 電子速度, 衝突電離, 発熱率,cpu時間
c	output.data読み込み
	read(8,*) (sw(i),	i = 1, 8)
	read(8,*) ddisp			!画面出力ステップ
	read(8,*) doutc			!電流出力ステップ
	read(8,*) douti			!経過計算出力ステップ
	read(8,*) douta			!平均状態出力ステップ
	read(8,*) dsave			!バイナリー出力ステップ
	read(8,*) dheat			!発熱率出力ステップ
	read(8,*) dcput			!cpu時間出力ステップ
	read(8,*) doutpi			!経過計算出力ステップ	

	jdisp	= max(0,nint(ddisp/dt))
c	write(*,*)  jdisp
	joutc	= jdata(doutc ,dt)
	jouti	= jdata(douti ,dt)
	jouta	= jdata(douta ,dt)
	jsave	= jdata(dsave ,dt)
	jheat	= jdata(dheat ,dt)
	jcput	= jdata(dcput ,dt)
	joutpi	= jdata(doutpi ,dt)
	end subroutine outinit


	integer function jdata(ddata,dt)
	real ddata,dt

	if(ddata.gt.0.0)then
		jdata	  = max(nint(ddata /dt),1)
	elseif(ddata.lt.0.0)then
		jdata	  = -1
	endif

	end function jdata