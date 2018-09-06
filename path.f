	subroutine path(dx,z,dhet,cxpole1,cxpole2,
     &			pass,pass_r,x_start_1,x_start_r,xx_goal,scatpoint,			!120817sato
     &			x_mean_free_path_sum,x_mean_free_path_count,
     &			z_start_1,z_start_r,zz_goal,mean_free_path_sum,
     &			mean_free_path_sum2,xcenter)

c---変数配列用パラメータ---
	include 'arraysize.fi'
	common	/arraydata/nx,nz
	integer	nx,nz

c---基本パラメータ---
	real	dx
	real	z		!07/03/15
	real	dhet(0:nlayer)						!ヘテロ界面(m)
	real	cxpole1(npole),cxpole2(npole)		!ゲート電極(2)

c-----	領域通過する粒子をカウント120817sato
	real(8),dimension (10) :: pass,pass_r		!カウンタ，rはリジェクション対策
	real	x_start_1				!ループさせる変数
	real	x_start_r				!rejection対策
	real	xx_goal					!driftが終わる位置を保存する変数
	integer	scatpoint			!散乱．空散乱判定
	real	x_mean_free_path_sum
	integer x_mean_free_path_count

	real	z_start_1				!ループさせる変数
	real	z_start_r				!rejection対策
	real	zz_goal					!driftが終わる位置を保存する変数
	real	mean_free_path_sum
	real,dimension (0:nx,3) :: 	mean_free_path_sum2		!平均自由行程x座標	 1:カウンタ,2:x方向平均自由行程,3:平均自由行程
	integer xcenter				!x_startとxx_goalの中心点
		
c-----	ローカル変数  -------
	real	regionA,regionB						!指定する領域

c==================================================================================

	regionA = cxpole1(2)		!!とりあえずゲート直下を指定
	regionB = cxpole2(2)

c==================================================================================

c	ドリフト開始位置が初期状態ならxx_goalを入れて平均自由行程をカウントしない
	if(x_start_1.eq.0.0)then
		x_start_1 = xx_goal
		x_start_r = x_start_1
		scatpoint = 0			!平均自由行程をカウントするループを飛ばすため
	endif
	if(z_start_1.eq.0.0)then
		z_start_1 = zz_goal
		z_start_r = z_start_1
		scatpoint = 0			!平均自由行程をカウントするループを飛ばすため
	endif

c-----pass(1)A外から領域通過,(2)A外から領域内,(3)領域内からB外,(4)領域内からA外,
c-----pass(5)領域内連続,(6)領域連続後B外,(7)領域連続後A外	 (10)領域内連続(仮)
c-----pass(8)B外から領域通過,(9)B外から領域内
					
		!	平均自由行程  !

	if(scatpoint.ge.1)then		!散乱が起きた時（空散乱(scatpoint=0)にならなかった時）

	if((dhet(nchannel1).lt.z).and.(z.lt.dhet(nchannel2)))then	!channel内に限定

	if((x_start_1.lt.regionA).and.(xx_goal.gt.regionB)) then		!x |AB| xx
		pass(1) = pass(1) +1.0
		pass_r(1) = pass_r(1) +1.0	
c		write(*,*)'pass(1)',pass(1)
c		write(852,*)'pass(1)',pass(1)
c		write(*,*)'pass_r(1)',pass_r(1)
c		write(852,*)'pass_r(1)',pass_r(1)
	endif
	if((x_start_1.lt.regionA).and.(xx_goal.gt.regionA)
     &			.and.(xx_goal.lt.regionB)) then		!x |A xx B|
		pass(2) = pass(2) +1.0
		pass_r(2) = pass_r(2) +1.0	
c		write(*,*)'pass(2)',pass(2)
c		write(852,*)'pass(2)',pass(2)
c		write(*,*)'pass(2)_r',pass_r(2)	
c		write(852,*)'pass_r(2)',pass_r(2)
	
	endif	 
	if((x_start_1.gt.regionA).and.(x_start_1.lt.regionB)
     &			.and.(xx_goal.gt.regionB)) then		!|A x B| xx
			pass(3) = pass(3) +1.0
			pass_r(3) = pass_r(3) +1.0
		if(pass(10).gt.0)then
			pass(6) = pass(6) +pass(10)
			pass_r(6) = pass_r(6) +pass_r(10)
			pass(10) = 0.0
			pass_r(10) = 0.0
		endif
	endif
	if((xx_goal.lt.regionA).and.(x_start_1.gt.regionA)
     &			.and.(x_start_1.lt.regionB)) then		!xx |A x B|
			pass(4) = pass(4) +1.0
			pass_r(4) = pass_r(4) +1.0
		if(pass(10).gt.0)then
			pass(7) = pass(7) +pass(10)
			pass_r(7) = pass_r(7) +pass_r(10)
			pass(10) = 0.0
			pass_r(10) = 0.0
		endif
	endif	
	if((x_start_1.gt.regionA).and.(x_start_1.lt.regionB)
     &		.and.(xx_goal.gt.regionA).and.(xx_goal.lt.regionB)) then		!|A x xx B|
		pass(5) = pass(5) +1.0
		pass_r(5) = pass_r(5) +1.0
		pass(10) = pass(10) +1.0
		pass_r(10) = pass_r(10) +1.0	
	endif	
	
	if((xx_goal.lt.regionA).and.(x_start_1.gt.regionB)) then		!xx |AB| x
		pass(8) = pass(8) +1.0
		pass_r(8) = pass_r(8) +1.0	
	endif
	if((x_start_1.gt.regionB).and.(xx_goal.gt.regionA)
     &			.and.(xx_goal.lt.regionB)) then		!|A xx B| x
		pass(9) = pass(9) +1.0
		pass_r(9) = pass_r(9) +1.0	
	endif
	
		x_mean_free_path_count = x_mean_free_path_count +1			!カウンタ
		x_mean_free_path_sum = x_mean_free_path_sum+abs(xx_goal-x_start_1)		!平均自由行程の合計をカウンタで割る(output)
		mean_free_path_sum = mean_free_path_sum
     &		+sqrt((xx_goal-x_start_1)**2+(zz_goal-z_start_1)**2)

		xcenter = nint((x_start_1+xx_goal)/dx/2)		!driftする中心点(整数型)
		mean_free_path_sum2(xcenter,1) 
     &			= mean_free_path_sum2(xcenter,1) +1.0		!カウンタ
		mean_free_path_sum2(xcenter,2) 	
     &			= mean_free_path_sum2(xcenter,2) +abs(xx_goal-x_start_1)
		mean_free_path_sum2(xcenter,3) = mean_free_path_sum2(xcenter,3) 
     &		+sqrt((xx_goal-x_start_1)**2+(zz_goal-z_start_1)**2)

	endif	 !channel内に限定

	x_start_1 = xx_goal		!散乱した箇所を次のドリフト開始位置にする
	z_start_1 = zz_goal		!散乱した箇所を次のドリフト開始位置にする
		
	endif	!scatpoint

	return
	end