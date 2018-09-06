c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c Effective Potential の出力を行うサブルーチン
c
c 実効ポテンシャルの出力
c 2007/1/9, Hara, Ver.1.1
c
c	Input : ict, u_eff
c	Output: NOTHING

c	解説：平均実効ポテンシャルを出力	
c	出力先：eff_potential.txt,eff_potential2.txtに出力
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	subroutine eff_out(ict, u_eff1,u_eff2,u_eff3,epA,epB,epC,
     &					epA2,epB2)	!120126homma
	implicit none
c
c===引数===
c---変数配列用パラメータ---
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c--- シミュレーション条件 ---
	common /outp/joutc,jouti,jouta,jsave,jheat,jcput,joutpi,sw
	integer	joutc,jouti,jouta,jsave,jheat,jcput,joutpi,sw(8)
c---基本パラメータ---
	integer	ict
c---デバイス内状態--- 120126homma
	real,	dimension (0:nx,0:nz)	:: u_eff1,u_eff2,u_eff3
	real,	dimension (0:nx,0:nz)	:: epA,epB,epC,epA2,epB2
c----ローカル変数---- 120126homma
	integer	i,j
	integer,save	::	ncount
	real(8),save,allocatable 		:: avu_eff1(:,:)
	real(8),save,allocatable 		:: avu_eff2(:,:)
	real(8),save,allocatable		:: avu_eff3(:,:)
	integer ix,iz
c
c=========================================================================================================
c
c----(リアルタイム値epA,epB,epC)--- !120126homma
	if (.not. allocated(avu_eff1)) then		!最初だけ実行
		allocate(avu_eff1(0:nx,0:nz))
		allocate(avu_eff2(0:nx,0:nz))
		allocate(avu_eff3(0:nx,0:nz))
		avu_eff1 = 0.0
		avu_eff2 = 0.0
		avu_eff3 = 0.0
		ncount	= 0
	endif
	avu_eff1 = avu_eff1 + u_eff1
	avu_eff2 = avu_eff2 + u_eff2
	avu_eff3 = avu_eff3 + u_eff3
	ncount	= ncount + 1
	epA = u_eff1
	epB = u_eff2
	epC = u_eff3
c----(平均値epA2,epB2)------------- !120126homma
	if(modulo(ict,jouta).eq.(jouta-1))then
		avu_eff1 = avu_eff1/ncount
		avu_eff2 = avu_eff2/ncount
		avu_eff3 = avu_eff3/ncount
		rewind 72; rewind 76; rewind 77; rewind 78
		do j=0,nz
		do i=0,nx
			write(72,*) sngl(avu_eff1(i,j))	!eff_potential2.txt 出力
			write(76,*) sngl(avu_eff2(i,j))	!eff_epB2.txt
			write(77,*) sngl(epA(i,j))		!eff_epA.txt
			write(78,*) sngl(epB(i,j))		!eff_epB.txt
		enddo
		enddo
		epA2 = avu_eff1
		epB2 = avu_eff2
		avu_eff1 = 0.0
		avu_eff2 = 0.0
		avu_eff3 = 0.0
		ncount = 0
	endif
c----------------------------------
c	
	return
	end
