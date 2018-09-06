c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c 量子補正を計算するサブルーチン
c
c 実効ポテンシャルの領域分類(２次元)
c 2007/1/23, Hara, Ver.1.0
c 解説：表面に近づいたときの、αの値を決めるための準備(EPVer3.0に対応)
c
c	atype=0が通常の領域、atype=1が表面の格子点
c	atype=2が表面から１つ離れた格子点
c	atype=3が表面から2つ離れた格子点
c
c	Ver1.0	リリースバージョン
c
c	Input : lxr1,lzr1,lxr2,lzr2
c	Output: atype
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	subroutine area_type(lxr1,lzr1,lxr2,lzr2,atype)
	implicit none
c
c===引数===
c---変数配列用パラメータ---
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c----(リセス)---
	integer(2),dimension (nrecess) :: lxr1,lzr1,lxr2,lzr2
c
	integer(1),	dimension (0:nx,0:nz)	:: atype
c	
c----(ローカル変数)----
	integer ix,iz,ii
c
	atype = 0
c
	do iz=0,nz
	  atype(2,iz)=3
	  atype(nx-2,iz)=3
	enddo
	do ix=0,nx
	  atype(ix,2)=3
	  atype(ix,nz-2)=3
	enddo
c	
	do ii=1,nrecess
	  do iz=lzr1(ii),lzr2(ii)
		atype(lxr1(ii)-2,iz)=3
		atype(lxr2(ii)+2,iz)=3
	  enddo
	  do ix=lxr1(ii),lxr2(ii)
		atype(ix,lzr2(ii)+2)=3
	  enddo
	enddo
c
	do iz=0,nz
	  atype(1,iz)=2
	  atype(nx-1,iz)=2
	enddo
	do ix=0,nx
	  atype(ix,1)=2
	  atype(ix,nz-1)=2
	enddo
c	
	do ii=1,nrecess
	  do iz=lzr1(ii),lzr2(ii)
		atype(lxr1(ii)-1,iz)=2
		atype(lxr2(ii)+1,iz)=2
	  enddo
	  do ix=lxr1(ii),lxr2(ii)
		atype(ix,lzr2(ii)+1)=2
	  enddo
	enddo
c
	do iz=0,nz
	  atype(0,iz)=1
	  atype(nx,iz)=1
	enddo
	do ix=0,nx
	  atype(ix,0)=1
	  atype(ix,nz)=1
	enddo
c
	do ii=1,nrecess
	  do iz=lzr1(ii),lzr2(ii)
		atype(lxr1(ii),iz)=1
		atype(lxr2(ii),iz)=1
	  enddo
	  do ix=lxr1(ii),lxr2(ii)
		atype(ix,lzr2(ii))=1
	  enddo
	enddo
c
	return
	end
