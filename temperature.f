	subroutine temperature(mxpart1,mxpart2,ntemp,btmp,dtmp)
c---	デバイス内温度分布をファイルから読み込みテーブルにするサブルーチン
c
	implicit none
c === 変数解説 ===
c ---	input ---
c	ntenum ... 温度分布分割数
c	mxpart1,mxpart2 ... 	各エリアの座標[mesh]
c ---	output ---
c	btmp ... ベース温度(デバイス内の最低温度)
c	dtmp ... ntemp1当たりの温度上昇分ΔT
c	ntemp(ix,iz) ... 各メッシュのbtmpからの温度上昇分テーブル
c
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c
	integer(2) mxpart1,mxpart2
	integer(2) ntemp(0:nx,0:nz)		!温度分布テーブル
	real	btmp,dtmp
c
c###	ローカル変数  ###
	integer	ix,iz
	real	tmin,tmax
	real,allocatable ::temp(:,:)
	real	tem
	character(80) form
c
c---( デバイス構造 )---
	allocate(temp(0:nx,0:nz))
c	ntemp = 1
	ntenum = ntemax
	btmp = 300.0		!ベース温度[K]
	dtmp = 0.2			!温度刻みΔT
c
	tmax = 0.0
	tmin = huge(tmin)
c---	ファイルから温度分布読み込み(n+領域は除く)   ---
c---	エラー時goto 100,ファイルend時goto 200	---
	do iz = 0,nz
	do ix = mxpart1,mxpart2
		read(8,'(f15.7)',ERR=100,END=200)tem
		tmin = min(tmin,tem)
		tmax = max(tmax,tem)
		temp(ix,iz)=tem
	enddo
c
c---	n+領域のデータ補充   ---
	do ix = mxpart2+1,nx
		temp(ix,iz)=temp(mxpart2,iz)
	enddo
	do ix = 0,mxpart1-1
		temp(ix,iz)=temp(mxpart1,iz)
	enddo
	enddo
c
	btmp = tmin		!ベース温度[K]
	dtmp = (tmax-tmin)/(ntenum)	!0.2			!温度刻みΔT
c
	if(tmax.eq.tmin)then
		form= '(''温度変化がみられません'')'
		write(*,form)
		write(99,form)
		goto 400
		return
	elseif((tmax.lt.tmin).or.(tmin.lt.0.0))then
		form= '(''温度設定がおかしくなっています。'')'
		write(*,form)
		write(99,form)
		stop
	endif
c
c---	テーブル変換	 ---
	do ix = 0,nx
	do iz = 0,nz
		ntemp(ix,iz) = ifix((temp(ix,iz)-btmp)/dtmp)+1
		ntemp(ix,iz) = min(int2(ntenum),max(int2(1),ntemp(ix,iz)))
	enddo
	enddo
	deallocate(temp)
	return
c
c---	エラー処理	 ---
200	continue
100	continue
	form = "(x,'温度分布ファイルが読み込めませんでした。')"
	write(*,form)
	write(99,form)
	if(tmax.ge.tmin) btmp = tmin
c
400	continue
	dtmp = 0.0			!温度刻みΔT
	ntemp = 1			!ntemp:全体配列
	ntenum = 1
c
	form = "(x,'全て',F9.4,'Kで計算します。')"
	write(*,form)btmp
	write(99,form)btmp
c
	return
c	
	end
