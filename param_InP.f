	subroutine param2(
     &			  de,dx,dz,btmp,dtmp,dconc,smh,hhm,hm,am_aniso,am,
     &			  aff_aniso,hole_am_aniso,hole_aff_aniso,af,af2,af4,
     &			  eg,ec,eps,bktq, 
     &			  swk,pgm,escat,iband,iarg,
     &			  hescat,dn3,hiXL,Eth_table,aff,dltec,					!09/2/19 竹岸
     &			  II_S)		!120921sato
c
c===( 物理定数, 材料定数および諸パラメータ )====
c
c === 変数解説 ===
c ---	input ---
c	nvalley ... 谷の数
c	nenergy ... エネルギーでの分割
c	nemax ... エネルギーステップの最大値
c	nscat ... 考慮する散乱機構の総数
c	de(iv) ... swkのkのエネルギー計算におけるエネルギーステップ[eV]
c	dx,dz ... メッシュサイズ[m]
c	temp ... 温度[K]
c	dn1 ... 不純物散乱レートの計算に代表値として用いる濃度[m^-3]
c
c ---	output ---
c	smh(iv) ... iv番目の谷の √2m*/h
c	hhm(iv) ... iv番目の谷の h^2/(2m*q)
c	hm(iv) ... iv番目の谷の h/m*^2
c	af,af2,af4 ... 非放物線バンドパラメータα,2α,4α
c	eg ... 禁制帯幅Eg
c	eps ... 静的誘電率
c	bktq ... ボルツマンファクター KbTq
c	swk(iv,j,iscat) ... 規格化された散乱レート(iv:谷,j:散乱機構,iscat:散乱前のエネルギー)
c	gm ... 散乱レートを規格化するΓ[1/s]
c	pgm ... (1/Γ) [s]
c	escat(iscat,iv) ... 散乱時に電子が得る(失う)エネルギー
c	iarg(iscat,iv) ... 終状態における電子の方向を決めるパラメータ
c	iarg(iv,j)= 1:等方的散乱, 2:有極性フォノン散乱, 3:不純物散乱
c	iband(iv,j) ... 散乱による遷移先の谷
c	hwo(iv) ... iv谷の有極性光学フォノンのエネルギー
c	hw(iv,j)	= バンド端散乱のフォノンエネルギーhwij（i→j）
c
	implicit none
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum

c	real	dconc(npole)
	real	dconc(npart)	!07/2/9
	real	dx,dz,btmp,dtmp
	integer ia
	integer ipart		!07/8/4 不純物散乱
	real, dimension	(nvalley,narea)::smh,hhm,hm,af,af2,af4,eg,dltec
c	real, dimension (nvalley,narea)::am									!120201
	real, dimension	(nvalley,narea) :: ec 
	real	eps(narea),bktq(ntenum),de(nenergy)
c	real	swk(nscat,nemax,nvalley,nenergy,ntenum,narea)
c	real	pgm(nvalley,nenergy,narea),escat(nscat,nvalley,narea)
	real	swk(nscat,nemax,nvalley,nenergy,ntenum,npart)	!07/8/4 不純物散乱
	real	pgm(nvalley,nenergy,npart),escat(nscat,nvalley,narea) !07/8/4 不純物散乱    
	integer(1),dimension (nscat,narea)	:: iarg
	integer(1),dimension (nscat,nvalley,narea)	:: iband
c
	real,	dimension (nscat,nvalley,narea) :: hescat
c
c	----------------------------------------
	real	qd21,bimp,temp
	real,dimension(:),allocatable:: z,dos,poe,poa,aco
	real,dimension(:,:),allocatable:: ope,opa
c
	real,	dimension (:,:),allocatable		:: hwo		!光学フォノンエネルギー
	real,	dimension (:,:,:),allocatable	:: hw		!光学フォノンバンド端エネルギー
c	----------------------------------------
c	real,dimension(:,:),allocatable:: am
c	real,dimension(:,:,:),allocatable:: am_aniso
	real,dimension(:,:,:),allocatable:: gm
	real	dn1
c	real,dimension(:),allocatable::	dn3		!07/8/4 不純物散乱
	real	dn3(npart)						!09/2/19 竹岸
	real(8) egmin(narea),epf(narea),ep(narea),rou(narea)
	real(8) aff(nvalley,narea)
c	real	hwnp,qeps
	real(8) sv(narea),cl(narea)
c	real(8)	qh
	real am(nvalley,narea)
	real(8), dimension(:),allocatable :: wo,no
	real(8), dimension(:,:),allocatable :: w,n
	real(8), dimension(:,:),allocatable	:: da
	real(8), dimension(:,:,:),allocatable :: d
c	real(8)	dop
	integer iv,jv,para_num1,para_num2,para_num3,para_num4,para_num5
	integer material_num1,material_num2,material_num3,
     &	material_num4,material_num5,material_num6,material_num7
      integer  material_num8,material_num9,material_num10
c
	real(8) pi,q,h,bk,ep0,am0
	parameter(pi  = 3.141592654, q   = 1.60219e-19)
	parameter(h   = 1.05459e-34, bk  = 1.38066e-23)
	parameter(ep0 = 8.85419e-12, am0 = 9.10953e-31)
c
	real(8)	ei,sei
	integer ie,iscat,ien,itp
	integer nemin
c
c---(合金散乱)---
	real(8),dimension(:),allocatable :: alloy	
	real,dimension(:),allocatable :: inx,cx,lc,va,ea
c
c---(衝突電離)---
	real, dimension	(nvalley,narea,nvalley)::am_aniso,aff_aniso						!20100624
	real, dimension	(nvalley,narea,nvalley)::hole_am_aniso							!20100624	
	real, dimension	(nvalley,narea,nvalley)::hole_aff_aniso							!20100624
c	real, dimension	(nvalley,narea,nvalley)::temp_am_aniso							!20100624
c	real, dimension	(nvalley,narea,nvalley)::temp_aff_aniso							!20100624
c	real, dimension	(nvalley,narea,nvalley)::temp_hole_am_aniso						!20100624
c	real, dimension	(nvalley,narea,nvalley)::temp_hole_aff_aniso					!20100624
c	real, dimension	(narea,nvalley)::highestXL					!20100624
c	real, dimension	(narea,nvalley)::hiXL					!20100624
c	real, dimension	(5,3) :: hiXL					!20100624
	real, dimension (narea,nvalley) :: hiXL
	real, dimension	(5,34) :: temp_read					!20100624
	real, dimension	(5,34) :: temp_material 
	real, dimension	(4,20000,2) :: Eth_table											!20100624 
	real, dimension(:),allocatable ::  eth,a,b
	real II_S(narea)		!120921sato






	real temp_param1(12)
	real temp_param2(12)
	integer i,j,k


c----(ループ用変数100826)---
	integer roop_num,roop_num2

c
c---------------------------------------------------------
	if(nvalley.ne.3)then
		write(*,*)'nvalleyの値が不正です'
		write(99,*)'nvalleyの値が不正です'
		stop
	endif
c
	allocate (z(nvalley),dos(nvalley))
	allocate (poe(nvalley),poa(nvalley),aco(nvalley))
	allocate (ope(nvalley,nvalley),opa(nvalley,nvalley))
c	allocate (am(nvalley,narea))
c	allocate (am_aniso(nvalley,narea,nvalley))
	allocate (wo(nvalley),no(nvalley))
	allocate (da(nvalley,narea))
	allocate (w(nvalley,nvalley))
	allocate (n(nvalley,nvalley))
	allocate (d(nvalley,nvalley,narea))
	allocate (hwo(nvalley,narea))		!光学フォノンエネルギー
	allocate (hw(nvalley,nvalley,narea))	!光学フォノンバンド端エネルギー
	allocate (alloy(nvalley))
	allocate (inx(narea),cx(narea),lc(narea),va(narea),ea(narea))
	allocate (eth(narea),a(narea),b(narea))
	swk	= 0.0		!Fortran90で可能な効率的表現
c
c---( 谷の数 )---
	z(1) = 1.0d0
	z(2) = 4.0d0
	z(3) = 3.0d0

	temp_material = 0d0

c	テキストからパラメータの読み込み
	open(unit=199,file='param1.txt')	!InSb
	open(unit=200,file='param2.txt')	!Al0.15In0.85Sb
	open(unit=201,file='param3.txt')	!Al0.25In0.75Sb
	open(unit=202,file='param4.txt')	!AlSb
	open(unit=203,file='param5.txt')

	open(unit=204,file='material1.txt')
	open(unit=205,file='material2.txt')
	open(unit=206,file='material3.txt')
	open(unit=207,file='material4.txt')
	open(unit=208,file='set_paramaterial.txt')
	open(unit=209,file='material5.txt')

c	write(*,*) npart

c	Eg依存Ethテーブル読み取り（X,Y,Z軸別）	
	open(unit=210,file='Eth_table.txt')

	do i=1,20000
	read(210,*) Eth_table(1,i,1), Eth_table(2,i,1), 
     &				Eth_table(3,i,1), Eth_table(4,i,1)	
	end do

	open(unit=211,file='Eth_table2.txt')

	do i=1,20000
	read(211,*) Eth_table(1,i,2), Eth_table(2,i,2), 
     &				Eth_table(3,i,2), Eth_table(4,i,2)	
	end do

c	set_paramaterial
	read(208,*) para_num1
	read(208,*) para_num2
	read(208,*)	para_num3
	read(208,*) para_num4
	read(208,*) para_num5
	read(208,*)	material_num1
	read(208,*) material_num2
	read(208,*) material_num3
	read(208,*)	material_num4
	read(208,*) material_num5
	read(208,*) material_num6
	read(208,*)	material_num7
	read(208,*)	material_num8
	read(208,*)	material_num9
	read(208,*)	material_num10

c	para_num1 = 1
c	para_num2 = 2
c	para_num3 = 3
c	para_num4 = 4
c	para_num5 = 5

c	material_num1 = 1
c	material_num2 = 2
c	material_num3 = 2
c	material_num4 = 3
c	material_num5 = 2
c	material_num6 = 2
c	material_num7 = 4
c	material_num8 = 4
c	material_num9 = 2
c	material_num10 = 3


c	material1_InSb (bak_In(0.52)Al(0.48)As)
	read(204,*) temp_material(1,1)
	read(204,*) temp_material(1,2)
	read(204,*) temp_material(1,3)
	read(204,*) temp_material(1,4)
	read(204,*) temp_material(1,5)
	read(204,*) temp_material(1,6)
	read(204,*) temp_material(1,7)
	read(204,*) temp_material(1,8)
	read(204,*) temp_material(1,9)
	read(204,*) temp_material(1,10)
	read(204,*) temp_material(1,11)
	read(204,*) temp_material(1,12)
	read(204,*) temp_material(1,13)
	read(204,*) temp_material(1,14)
	read(204,*) temp_material(1,15)
	read(204,*) temp_material(1,16)
	read(204,*) temp_material(1,17)
	read(204,*) temp_material(1,18)
	read(204,*) temp_material(1,19)
	read(204,*) temp_material(1,20)
	read(204,*) temp_material(1,21)
	read(204,*) temp_material(1,22)
	read(204,*) temp_material(1,23)
	read(204,*) temp_material(1,24)
	read(204,*) temp_material(1,25)
	read(204,*) temp_material(1,26)

c	material2_Al0.15In0.85Sb (bak_InAs)
	read(205,*) temp_material(2,1)
	read(205,*) temp_material(2,2)
	read(205,*) temp_material(2,3)
	read(205,*) temp_material(2,4)
	read(205,*) temp_material(2,5)
	read(205,*) temp_material(2,6)
	read(205,*) temp_material(2,7)
	read(205,*) temp_material(2,8)
	read(205,*) temp_material(2,9)
	read(205,*) temp_material(2,10)
	read(205,*) temp_material(2,11)
	read(205,*) temp_material(2,12)
	read(205,*) temp_material(2,13)
	read(205,*) temp_material(2,14)
	read(205,*) temp_material(2,15)
	read(205,*) temp_material(2,16)
	read(205,*) temp_material(2,17)
	read(205,*) temp_material(2,18)
	read(205,*) temp_material(2,19)
	read(205,*) temp_material(2,20)
	read(205,*) temp_material(2,21)
	read(205,*) temp_material(2,22)
	read(205,*) temp_material(2,23)
	read(205,*) temp_material(2,24)
	read(205,*) temp_material(2,25)
	read(205,*) temp_material(2,26)

c	material3_Al0.25In0.75Sb (bak_GaAs)
	read(206,*) temp_material(3,1)
	read(206,*) temp_material(3,2)
	read(206,*) temp_material(3,3)
	read(206,*) temp_material(3,4)
	read(206,*) temp_material(3,5)
	read(206,*) temp_material(3,6)
	read(206,*) temp_material(3,7)
	read(206,*) temp_material(3,8)
	read(206,*) temp_material(3,9)
	read(206,*) temp_material(3,10)
	read(206,*) temp_material(3,11)
	read(206,*) temp_material(3,12)
	read(206,*) temp_material(3,13)
	read(206,*) temp_material(3,14)
	read(206,*) temp_material(3,15)
	read(206,*) temp_material(3,16)
	read(206,*) temp_material(3,17)
	read(206,*) temp_material(3,18)
	read(206,*) temp_material(3,19)
	read(206,*) temp_material(3,20)
	read(206,*) temp_material(3,21)
	read(206,*) temp_material(3,22)
	read(206,*) temp_material(3,23)
	read(206,*) temp_material(3,24)
	read(206,*) temp_material(3,25)
	read(206,*) temp_material(3,26)

c	material4_AlSb (bak_InP)
	read(207,*) temp_material(4,1)
	read(207,*) temp_material(4,2)
	read(207,*) temp_material(4,3)
	read(207,*) temp_material(4,4)
	read(207,*) temp_material(4,5)
	read(207,*) temp_material(4,6)
	read(207,*) temp_material(4,7)
	read(207,*) temp_material(4,8)
	read(207,*) temp_material(4,9)
	read(207,*) temp_material(4,10)
	read(207,*) temp_material(4,11)
	read(207,*) temp_material(4,12)
	read(207,*) temp_material(4,13)
	read(207,*) temp_material(4,14)
	read(207,*) temp_material(4,15)
	read(207,*) temp_material(4,16)
	read(207,*) temp_material(4,17)
	read(207,*) temp_material(4,18)
	read(207,*) temp_material(4,19)
	read(207,*) temp_material(4,20)
	read(207,*) temp_material(4,21)
	read(207,*) temp_material(4,22)
	read(207,*) temp_material(4,23)
	read(207,*) temp_material(4,24)
	read(207,*) temp_material(4,25)
	read(207,*) temp_material(4,26)

c	material5_InSb
	read(209,*) temp_material(5,1)
	read(209,*) temp_material(5,2)
	read(209,*) temp_material(5,3)
	read(209,*) temp_material(5,4)
	read(209,*) temp_material(5,5)
	read(209,*) temp_material(5,6)
	read(209,*) temp_material(5,7)
	read(209,*) temp_material(5,8)
	read(209,*) temp_material(5,9)
	read(209,*) temp_material(5,10)
	read(209,*) temp_material(5,11)
	read(209,*) temp_material(5,12)
	read(209,*) temp_material(5,13)
	read(209,*) temp_material(5,14)
	read(209,*) temp_material(5,15)
	read(209,*) temp_material(5,16)
	read(209,*) temp_material(5,17)
	read(209,*) temp_material(5,18)
	read(209,*) temp_material(5,19)
	read(209,*) temp_material(5,20)
	read(209,*) temp_material(5,21)
	read(209,*) temp_material(5,22)
	read(209,*) temp_material(5,23)
	read(209,*) temp_material(5,24)
	read(209,*) temp_material(5,25)
	read(209,*) temp_material(5,26)

c
c	param1
c---(組成比)---
	read(199,*) temp_read(1,1)
c---( 伝導帯底のエネルギー )---
	read(199,*) temp_read(1,2)	!Ｌ谷
	read(199,*) temp_read(1,3)	!Ｘ谷
	read(199,*) temp_read(1,4)
c---( 電子の有効質量 )---
	read(199,*) temp_read(1,5)	!Γ谷(||)
	read(199,*) temp_read(1,6)	!Γ谷(⊥)
	read(199,*) temp_read(1,7)	!Ｌ谷(|)
	read(199,*) temp_read(1,8)	!Ｌ谷(t1)
	read(199,*) temp_read(1,9)	!Ｌ谷(t2)
	read(199,*) temp_read(1,10)	!X谷(|)
	read(199,*) temp_read(1,11)	!X谷(t1)
	read(199,*) temp_read(1,12)	!X谷(t2)
	read(199,*) temp_read(1,13)	!Z谷(|)
	read(199,*) temp_read(1,14)	!Z谷(t1)
	read(199,*) temp_read(1,15)	!Z谷(t2)
c---(平均した有効質量を仮に読み込む)---
	read(199,*) temp_read(1,16)		!m*(Γ)
	read(199,*) temp_read(1,17)		!m*(L)
	read(199,*) temp_read(1,18)		!(m*(X) =
	read(199,*) temp_read(1,19)		!m*(Z) =
	read(199,*) temp_read(1,20)		!m*(X,Z))
c---( 非放物線性パラメータα )---
	read(199,*) temp_read(1,21)	!Γ谷
	read(199,*) temp_read(1,22)	!Ｌ谷
	read(199,*) temp_read(1,23)
	read(199,*) temp_read(1,24)	!Ｘ谷
	read(199,*) temp_read(1,25)
	read(199,*) temp_read(1,26)
c---( 禁制帯幅 )---
	read(199,*) temp_read(1,27)
	read(199,*) temp_read(1,28)	!JAP94(2003)4096
c---( 格子定数,体積 )---
	read(199,*) temp_read(1,29)
c---( ホールの有効質量 )---
	read(199,*) temp_read(1,30)
	read(199,*) temp_read(1,31)
c---( ホールの非放物線性パラメータα )---
	read(199,*) temp_read(1,32)
c---( X谷、L谷の最大値_I.I.用)---
	read(199,*) temp_read(1,33)
	read(199,*) temp_read(1,34)
	read(199,*) temp_read(1,35)
	read(199,*) temp_read(1,36)
	read(199,*) temp_read(1,37)

c	param2
c---(組成比)---
	read(200,*) temp_read(2,1)
c---( 伝導帯底のエネルギー )---
	read(200,*) temp_read(2,2)	!Ｌ谷
	read(200,*) temp_read(2,3)	!Ｘ谷
	read(200,*) temp_read(2,4)
c---( 電子の有効質量 )---
	read(200,*) temp_read(2,5)	!Γ谷(||)
	read(200,*) temp_read(2,6)	!Γ谷(⊥)
	read(200,*) temp_read(2,7)	!Ｌ谷(|)
	read(200,*) temp_read(2,8)	!Ｌ谷(t1)
	read(200,*) temp_read(2,9)	!Ｌ谷(t2)
	read(200,*) temp_read(2,10)	!X谷(|)
	read(200,*) temp_read(2,11)	!X谷(t1)
	read(200,*) temp_read(2,12)	!X谷(t2)
	read(200,*) temp_read(2,13)	!Z谷(|)
	read(200,*) temp_read(2,14)	!Z谷(t1)
	read(200,*) temp_read(2,15)	!Z谷(t2)
c---(平均した有効質量を仮に読み込む)---
	read(200,*) temp_read(2,16)		!m*(Γ)
	read(200,*) temp_read(2,17)		!m*(L)
	read(200,*) temp_read(2,18)		!(m*(X) =
	read(200,*) temp_read(2,19)		!m*(Z) =
	read(200,*) temp_read(2,20)		!m*(X,Z))
c---( 非放物線性パラメータα )---
	read(200,*) temp_read(2,21)	!Γ谷
	read(200,*) temp_read(2,22)	!Ｌ谷
	read(200,*) temp_read(2,23)
	read(200,*) temp_read(2,24)	!Ｘ谷
	read(200,*) temp_read(2,25)
	read(200,*) temp_read(2,26)
c---( 禁制帯幅 )---
	read(200,*) temp_read(2,27)
	read(200,*) temp_read(2,28)	!JAP94(2003)4096
c---( 格子定数,体積 )---
	read(200,*) temp_read(2,29)
c---( ホールの有効質量 )---
	read(200,*) temp_read(2,30)
	read(200,*) temp_read(2,31)
c---( ホールの非放物線性パラメータα )---
	read(200,*) temp_read(2,32)
c---( X谷、L谷の最大値_I.I.用)---
	read(200,*) temp_read(2,33)
	read(200,*) temp_read(2,34)

c	param3
c---(組成比)---
	read(201,*) temp_read(3,1)
c---( 伝導帯底のエネルギー )---
	read(201,*) temp_read(3,2)	!Ｌ谷
	read(201,*) temp_read(3,3)	!Ｘ谷
	read(201,*) temp_read(3,4)
c---( 電子の有効質量 )---
	read(201,*) temp_read(3,5)	!Γ谷(||)
	read(201,*) temp_read(3,6)	!Γ谷(⊥)
	read(201,*) temp_read(3,7)	!Ｌ谷(|)
	read(201,*) temp_read(3,8)	!Ｌ谷(t1)
	read(201,*) temp_read(3,9)	!Ｌ谷(t2)
	read(201,*) temp_read(3,10)	!X谷(|)
	read(201,*) temp_read(3,11)	!X谷(t1)
	read(201,*) temp_read(3,12)	!X谷(t2)
	read(201,*) temp_read(3,13)	!Z谷(|)
	read(201,*) temp_read(3,14)	!Z谷(t1)
	read(201,*) temp_read(3,15)	!Z谷(t2)
c---(平均した有効質量を仮に読み込む)---
	read(201,*) temp_read(3,16)		!m*(Γ)
	read(201,*) temp_read(3,17)		!m*(L)
	read(201,*) temp_read(3,18)		!(m*(X) =
	read(201,*) temp_read(3,19)		!m*(Z) =
	read(201,*) temp_read(3,20)		!m*(X,Z))
c---( 非放物線性パラメータα )---
	read(201,*) temp_read(3,21)	!Γ谷
	read(201,*) temp_read(3,22)	!Ｌ谷
	read(201,*) temp_read(3,23)
	read(201,*) temp_read(3,24)	!Ｘ谷
	read(201,*) temp_read(3,25)
	read(201,*) temp_read(3,26)
c---( 禁制帯幅 )---
	read(201,*) temp_read(3,27)
	read(201,*) temp_read(3,28)	!JAP94(2003)4096
c---( 格子定数,体積 )---
	read(201,*) temp_read(3,29)
c---( ホールの有効質量 )---
	read(201,*) temp_read(3,30)
	read(201,*) temp_read(3,31)
c---( ホールの非放物線性パラメータα )---
	read(201,*) temp_read(3,32)
c---( X谷、L谷の最大値_I.I.用)---
	read(201,*) temp_read(3,33)
	read(201,*) temp_read(3,34)

c	param4
c---(組成比)---
	read(202,*) temp_read(4,1)
c---( 伝導帯底のエネルギー )---
	read(202,*) temp_read(4,2)	!Ｌ谷
	read(202,*) temp_read(4,3)	!Ｘ谷
	read(202,*) temp_read(4,4)
c---( 電子の有効質量 )---
	read(202,*) temp_read(4,5)	!Γ谷(||)
	read(202,*) temp_read(4,6)	!Γ谷(⊥)
	read(202,*) temp_read(4,7)	!Ｌ谷(|)
	read(202,*) temp_read(4,8)	!Ｌ谷(t1)
	read(202,*) temp_read(4,9)	!Ｌ谷(t2)
	read(202,*) temp_read(4,10)	!X谷(|)
	read(202,*) temp_read(4,11)	!X谷(t1)
	read(202,*) temp_read(4,12)	!X谷(t2)
	read(202,*) temp_read(4,13)	!Z谷(|)
	read(202,*) temp_read(4,14)	!Z谷(t1)
	read(202,*) temp_read(4,15)	!Z谷(t2)
c---(平均した有効質量を仮に読み込む)---
	read(202,*) temp_read(4,16)		!m*(Γ)
	read(202,*) temp_read(4,17)		!m*(L)
	read(202,*) temp_read(4,18)		!(m*(X) =
	read(202,*) temp_read(4,19)		!m*(Z) =
	read(202,*) temp_read(4,20)		!m*(X,Z))
c---( 非放物線性パラメータα )---
	read(202,*) temp_read(4,21)	!Γ谷
	read(202,*) temp_read(4,22)	!Ｌ谷
	read(202,*) temp_read(4,23)
	read(202,*) temp_read(4,24)	!Ｘ谷
	read(202,*) temp_read(4,25)
	read(202,*) temp_read(4,26)
c---( 禁制帯幅 )---
	read(202,*) temp_read(4,27)
	read(202,*) temp_read(4,28)	!JAP94(2003)4096
c---( 格子定数,体積 )---
	read(202,*) temp_read(4,29)
c---( ホールの有効質量 )---
	read(202,*) temp_read(4,30)
	read(202,*) temp_read(4,31)
c---( ホールの非放物線性パラメータα )---
	read(202,*) temp_read(4,32)
c---( X谷、L谷の最大値_I.I.用)---
	read(202,*) temp_read(4,33)
	read(202,*) temp_read(4,34)

c	param5
c---(組成比)---
	read(203,*) temp_read(5,1)
c---( 伝導帯底のエネルギー )---
	read(203,*) temp_read(5,2)	!Ｌ谷
	read(203,*) temp_read(5,3)	!Ｘ谷
	read(203,*) temp_read(5,4)
c---( 電子の有効質量 )---
	read(203,*) temp_read(5,5)	!Γ谷(||)
	read(203,*) temp_read(5,6)	!Γ谷(⊥)
	read(203,*) temp_read(5,7)	!Ｌ谷(|)
	read(203,*) temp_read(5,8)	!Ｌ谷(t1)
	read(203,*) temp_read(5,9)	!Ｌ谷(t2)
	read(203,*) temp_read(5,10)	!X谷(|)
	read(203,*) temp_read(5,11)	!X谷(t1)
	read(203,*) temp_read(5,12)	!X谷(t2)
	read(203,*) temp_read(5,13)	!Z谷(|)
	read(203,*) temp_read(5,14)	!Z谷(t1)
	read(203,*) temp_read(5,15)	!Z谷(t2)
c---(平均した有効質量を仮に読み込む)---
	read(203,*) temp_read(5,16)		!m*(Γ)
	read(203,*) temp_read(5,17)		!m*(L)
	read(203,*) temp_read(5,18)		!(m*(X) =
	read(203,*) temp_read(5,19)		!m*(Z) =
	read(203,*) temp_read(5,20)		!m*(X,Z))
c---( 非放物線性パラメータα )---
	read(203,*) temp_read(5,21)	!Γ谷
	read(203,*) temp_read(5,22)	!Ｌ谷
	read(203,*) temp_read(5,23)
	read(203,*) temp_read(5,24)	!Ｘ谷
	read(203,*) temp_read(5,25)
	read(203,*) temp_read(5,26)
c---( 禁制帯幅 )---
	read(203,*) temp_read(5,27)
c---( 衝突電離スレショルド)---
	read(203,*) temp_read(5,28)	!JAP94(2003)4096
c---( 格子定数,体積 )---
	read(203,*) temp_read(5,29)
c---( ホールの有効質量 )---
	read(203,*) temp_read(5,30)
	read(203,*) temp_read(5,31)
c---( ホールの非放物線性パラメータα )---
	read(203,*) temp_read(5,32)
c---( X谷、L谷の最大値_I.I.用)---
	read(203,*) temp_read(5,33)
	read(203,*) temp_read(5,34)
	read(203,*) temp_read(5,35)
	read(203,*) temp_read(5,36)
	read(203,*) temp_read(5,37)

c---------------------InSb(old_In(0.52)Al(0.48)As )------------
c
c---( In組成比 )---
c	inx(1) = temp_read(1,1)		!0.52
	inx(1) = temp_read(para_num1,1)		!0.52
	cx(1) = inx(1)*(1-inx(1))
c
c---( 格子定数,体積 )---
c	lc(1) = temp_read(1,29)*1.0e-10		!5.8687e-10
	lc(1) = temp_read(para_num1,29)*1.0e-10		!5.8687e-10
	va(1) = lc(1)**3/4.0

c---( 禁制帯幅 )---
	egmin(1) = temp_read(1,27)	!1.457
c
c---( フォノン散乱の諸パラメータ )---
	rou(1) = temp_material(material_num1,1)*
     &	inx(1)+temp_material(material_num2,1)*(1-inx(1))
    	!rou(1) = temp_material(material_num1,1) !4878	!半導体の比重(kg/m^3)
	sv(1)  = temp_material(material_num1,2)*
     &	inx(1)+temp_material(material_num2,2)*(1-inx(1))	!修正11/07/25原
     	!sv(1)  = temp_material(material_num1,2) !4679	!半導体中の音速(m/s)
c
c---( 誘電率 )---
	eps(1)	= (temp_material(material_num1,3)*
     &	inx(1)+temp_material(material_num2,3)*(1-inx(1)))*ep0
      !eps(1)	= temp_material(material_num1,3)*ep0 !12.42*ep0		!半導体の誘電率εs
	epf(1)  = (temp_material(material_num1,4)*
     &	inx(1)+temp_material(material_num2,4)*(1-inx(1)))*ep0
      !epf(1)	= temp_material(material_num1,4)*ep0 !10.23*ep0		!光学的誘電率ε∞
	ep(1)   = 1.0/(1.0/epf(1)-1.0/eps(1))
c
c---( 電子の有効質量 )---
c	am(1,1)	= temp_read(1,16)*am0		!0.083*am0	!Γ谷
c	am(2,1)	= temp_read(1,17)*am0		!0.304*am0	!Ｌ谷
c	am(3,1)	= temp_read(1,18)*am0		!0.496*am0	!Ｘ谷

	am(1,1)	= temp_read(para_num1,16)*am0		!0.083*am0	!Γ谷
	am(2,1)	= temp_read(para_num1,17)*am0		!0.304*am0	!Ｌ谷
	am(3,1)	= temp_read(para_num1,18)*am0		!0.496*am0	!Ｘ谷
	
c	am_aniso(1,1,1) = temp_read(1,5)*am0	!0.083*am0
c	am_aniso(1,1,2) = temp_read(1,5)*am0	!0.083*am0
c	am_aniso(1,1,3) = temp_read(1,6)*am0	!0.083*am0
c	am_aniso(2,1,1) = temp_read(1,7)*am0	!0.304*am0
c	am_aniso(2,1,2) = temp_read(1,7)*am0	!0.304*am0
c	am_aniso(2,1,3) = temp_read(1,7)*am0	!0.304*am0
c	am_aniso(3,1,1) = temp_read(1,10)*am0	!0.304*am0
c	am_aniso(3,1,2) = temp_read(1,10)*am0	!0.304*am0
c	am_aniso(3,1,3) = temp_read(1,10)*am0	!0.304*am0

	am_aniso(1,1,1) = temp_read(para_num1,5)*am0	!0.083*am0
	am_aniso(1,1,2) = temp_read(para_num1,5)*am0	!0.083*am0
	am_aniso(1,1,3) = temp_read(para_num1,6)*am0	!0.083*am0
	am_aniso(2,1,1) = temp_read(para_num1,7)*am0	!0.304*am0
	am_aniso(2,1,2) = temp_read(para_num1,7)*am0	!0.304*am0
	am_aniso(2,1,3) = temp_read(para_num1,7)*am0	!0.304*am0
	am_aniso(3,1,1) = temp_read(para_num1,10)*am0	!0.304*am0
	am_aniso(3,1,2) = temp_read(para_num1,10)*am0	!0.304*am0
	am_aniso(3,1,3) = temp_read(para_num1,10)*am0	!0.304*am0
c

c	hole_am_aniso(1,1,1) = temp_read(1,30)*am0	!0.57162*am0
c	hole_am_aniso(1,1,2) = temp_read(1,30)*am0	!0.57162*am0
c	hole_am_aniso(1,1,3) = temp_read(1,31)*am0	!0.57162*am0

	hole_am_aniso(1,1,1) = temp_read(para_num1,30)*am0	!0.57162*am0
	hole_am_aniso(1,1,2) = temp_read(para_num1,30)*am0	!0.57162*am0
	hole_am_aniso(1,1,3) = temp_read(para_num1,31)*am0	!0.57162*am0

c	hole_aff_aniso(1,1,1) = temp_read(1,32)	!55.7159
c	hole_aff_aniso(1,1,2) = temp_read(1,32)	!55.7159
c	hole_aff_aniso(1,1,3) = temp_read(1,32)	!55.7159

	hole_aff_aniso(1,1,1) = temp_read(para_num1,32)	!55.7159
	hole_aff_aniso(1,1,2) = temp_read(para_num1,32)	!55.7159
	hole_aff_aniso(1,1,3) = temp_read(para_num1,32)	!55.7159

c---( 非放物線性パラメータα )---
c	aff(1,1) = temp_read(1,21)		!0.543	!Γ谷
c	aff(2,1) = temp_read(1,22)		!0.415	!Ｌ谷
c	aff(3,1) = temp_read(1,24)		!0.204	!Ｘ谷

c	aff_aniso(1,1,1) =  temp_read(1,21)		!0.543
c	aff_aniso(1,1,2) =  temp_read(1,21)		!0.543
c	aff_aniso(1,1,3) =  temp_read(1,21)		!0.543
c	aff_aniso(2,1,1) =  temp_read(1,22)		!0.415
c	aff_aniso(2,1,2) =  temp_read(1,22)		!0.415
c	aff_aniso(2,1,3) =  temp_read(1,22)		!0.415
c	aff_aniso(3,1,1) =  temp_read(1,24)		!0.204
c	aff_aniso(3,1,2) =  temp_read(1,24)		!0.204
c	aff_aniso(3,1,3) =  temp_read(1,24)		!0.204

	aff(1,1) = temp_read(para_num1,21)		!0.543	!Γ谷
	aff(2,1) = temp_read(para_num1,22)		!0.415	!Ｌ谷
	aff(3,1) = temp_read(para_num1,24)		!0.204	!Ｘ谷

	aff_aniso(1,1,1) =  temp_read(para_num1,21)		!0.543
	aff_aniso(1,1,2) =  temp_read(para_num1,21)		!0.543
	aff_aniso(1,1,3) =  temp_read(para_num1,21)		!0.543
	aff_aniso(2,1,1) =  temp_read(para_num1,22)		!0.415
	aff_aniso(2,1,2) =  temp_read(para_num1,22)		!0.415
	aff_aniso(2,1,3) =  temp_read(para_num1,22)		!0.415
	aff_aniso(3,1,1) =  temp_read(para_num1,24)		!0.204
	aff_aniso(3,1,2) =  temp_read(para_num1,24)		!0.204
	aff_aniso(3,1,3) =  temp_read(para_num1,24)		!0.204

c
c---( 伝導帯底のエネルギー )---
c	ec(1,1)	= 0.0		!Γ谷
c	ec(2,1)	= temp_read(1,2)			!0.50		!Ｌ谷
c	ec(3,1)	= temp_read(1,3)			!0.60		!Ｘ谷

	ec(1,1)	= 0.0		!Γ谷
	ec(2,1)	= temp_read(para_num1,2)			!0.50		!Ｌ谷
	ec(3,1)	= temp_read(para_num1,3)			!0.60		!Ｘ谷
c
c---( 音響フォノンの変位ポテンシャル (da(i) = Ξd))---
	da(1,1)	= (temp_material(material_num1,5)*inx(1)+
     &	temp_material(material_num2,5)*(1-inx(1)))*q
c      da(1,1)	= temp_material(material_num1,5)*q !5.93*q	!Γ谷
	da(2,1)	= (temp_material(material_num1,6)*inx(1)+
     &	temp_material(material_num2,6)*(1-inx(1)))*q
c      da(2,1)	= temp_material(material_num1,6)*q !7.23*q	!Ｌ谷
	da(3,1)	= (temp_material(material_num1,7)*inx(1)+
     &	temp_material(material_num2,7)*(1-inx(1)))*q
c      da(3,1)	= temp_material(material_num1,7)*q !9.02*q	!Ｘ谷

c
c---( 光学フォノンの変位ポテンシャル(eV/m) (d(i,j) = Ｄij))---	 
	d(1,1,1) = temp_material(material_num1,8)*1.0e10*q !0.0			!ΓtoΓ
	d(2,1,1) = (temp_material(material_num1,9)*inx(1)+
     &	temp_material(material_num2,9)*(1-inx(1)))*1.0e10*q
c      d(2,1,1) = temp_material(material_num1,9)*1.0e10*q !5.25e10*q	!ΓtoＬ
	d(3,1,1) = (temp_material(material_num1,10)*inx(1)+
     &	temp_material(material_num2,10)*(1-inx(1)))*1.0e10*q
c      d(3,1,1) = temp_material(material_num1,10)*1.0e10*q !3.82e10*q	!ΓtoＸ
	d(2,2,1) = (temp_material(material_num1,11)*inx(1)+
     &	temp_material(material_num2,11)*(1-inx(1)))*1.0e10*q
c      d(2,2,1) = temp_material(material_num1,11)*1.0e10*q !6.55e10*q	!ＬtoＬ
	d(3,2,1) = (temp_material(material_num1,12)*inx(1)+
     &	temp_material(material_num2,12)*(1-inx(1)))*1.0e10*q
c      d(3,2,1) = temp_material(material_num1,12)*1.0e10*q !8.60e10*q	!ＬtoＸ
	d(3,3,1) = (temp_material(material_num1,13)*inx(1)+
     &	temp_material(material_num2,13)*(1-inx(1)))*1.0e10*q
c      d(3,3,1) = temp_material(material_num1,13)*1.0e10*q !5.72e10*q	!ＸtoＸ
c
c---(	有極性光学フォノンのエネルギー )---
	hwo(1:3,1) = temp_material(material_num1,14)*
     &	inx(1)+temp_material(material_num2,14)*(1-inx(1))
c      hwo(1:3,1) = temp_material(material_num1,14) !0.0404
c	hwo(1,1) = temp_material(material_num1,15) !0.0
c	hwo(2,1) = temp_material(material_num1,16) !0.0404
c	hwo(3,1) = temp_material(material_num1,17) !0.0
c
c---(	バンド端散乱のフォノンエネルギー(eV/m) (hw(i,j) = hwij))---
	hw(1,1,1) = 0.0			
	!hw(1,1,1) = temp_material(material_num1,18) !0.0		!ΓtoΓ
	hw(2,1,1) = (temp_material(material_num1,19)*inx(1)+
     &	temp_material(material_num2,19)*(1-inx(1)))
      !hw(2,1,1) = temp_material(material_num1,19) !0.0293	!ΓtoＬ
	hw(3,1,1) = (temp_material(material_num1,20)*inx(1)+
     &	temp_material(material_num2,20)*(1-inx(1)))
      !hw(3,1,1) = temp_material(material_num1,20) !0.0293	!ΓtoＸ
	hw(2,2,1) = (temp_material(material_num1,21)*inx(1)+
     &	temp_material(material_num2,21)*(1-inx(1)))
      !hw(2,2,1) = temp_material(material_num1,21) !0.0303	!ＬtoＬ
	hw(3,2,1) = (temp_material(material_num1,22)*inx(1)+
     &	temp_material(material_num2,22)*(1-inx(1)))
      !hw(3,2,1) = temp_material(material_num1,22) !0.0322	!ＬtoＸ
	hw(3,3,1) = (temp_material(material_num1,23)*inx(1)+
     &	temp_material(material_num2,23)*(1-inx(1)))
      !hw(3,3,1) = temp_material(material_num1,23) !0.0268	!ＸtoＸ
c
c

c
c---(合金散乱)---
c	ea = 0.0		!合金散乱非考慮
c	ea(1) = temp_material(material_num1,24) !0.47
c	ea(1) = (temp_material(material_num1,24)*inx(1)/0.52+
	ea(1) = (temp_material(material_num1,24)*inx(1)+
     &	temp_material(material_num2,24)*(1-inx(1)))
c	write(*,*) temp_material(material_num1,24)
c	write(*,*) temp_material(material_num2,24)
c	write(*,*) ea(1)
c
c---( X谷、L谷の最大値_I.I.用)---
c	hiXL(1,2) = temp_read(2,33)
c	hiXL(1,3) = temp_read(3,34)
	hiXL(1,2) = temp_read(para_num1,33)
	hiXL(1,3) = temp_read(para_num1,34)
c	highestX(1,1) = temp_read(1,33)
c	highestL(1,2) = temp_read(1,34)

c---( 衝突電離 )---
c	eth = 1000		!衝突電離非考慮		
c	eth(1) = temp_read(1,28)		!0.86		!JAP94(2003)4096
	eth(1) = temp_read(1,28)		!0.86		!JAP94(2003)4096
	a(1)	= temp_material(material_num1,25) !1.0e12
	II_S(1)	= a(1)		
	b(1)	= temp_material(material_num1,26) !2.0
c
c	eth = 0.80		!JAP94(2003)4096
c	a	= 2.0e12
c	b	= 2.0

c	eth = 0.75		!JAP96(2004)5650
c	a	= 2.0e10
c	b	= 2.5
c
c---------------------Al0.15In0.85Sb(old_In(0.53)Ga(0.47)As )-------------
c
c
c---( In組成比 )---		
	inx(2) = temp_read(para_num2,1)		!0.53
	cx(2) = inx(2)*(1-inx(2))
c
c---( 格子定数,体積 )---
	lc(2) = temp_read(2,29)*1.0e-10
	 !(temp_read(2,29)*inx(2)+5.65325*(1-inx(2)))*1.0e-10 !temp_read(2,29)*1.0e-10		!5.8687e-10
	va(2) = lc(2)**3/4.0

c---( 禁制帯幅 )---
	egmin(2) = temp_read(2,27)
	 !temp_read(para_num2,27)*inx(2)+1.424*(1-inx(2)) !temp_read(2,27)	!0.675
c
c---( フォノン散乱の諸パラメータ )---
	rou(2) = temp_material(material_num3,1)*
     &	inx(2)+temp_material(material_num4,1)*(1-inx(2))
     	!5469*inx(2)+5310*(1-inx(2))	!半導体の比重(kg/m^3)
	sv(2)  = temp_material(material_num3,2)*
     &	inx(2)+temp_material(material_num4,2)*(1-inx(2))	!修正11/07/25原
     	!4742*inx(2)+5240*(1-inx(2))	!半導体中の音速(m/s)
c
c---( 誘電率 )---
	eps(2)	= (temp_material(material_num3,3)*
     &	inx(2)+temp_material(material_num4,3)*(1-inx(2)))*ep0 !<-ミス!!
      !(13.88*inx(2)+12.90*(1-inx(2)))*ep0 !13.88*ep0		!半導体の誘電率εs
	epf(2)  = (temp_material(material_num3,4)*
     &	inx(2)+temp_material(material_num4,4)*(1-inx(2)))*ep0
      !(11.34*inx(2)+10.89*(1-inx(2)))*ep0 !11.34*ep0		!光学的誘電率ε∞
	ep(2)   = 1.0/(1.0/epf(2)-1.0/eps(2))
c
c---( 電子の有効質量 )---
	am(1,2)	= temp_read(para_num2,16)*am0		!0.04591*am0	!Γ谷
	am(2,2)	= temp_read(para_num2,17)*am0		!m0.17174*am0	!Ｌ谷
	am(3,2)	= temp_read(para_num2,18)*am0		!0.35054*am0	!Ｘ谷
	
	am_aniso(1,2,1) = temp_read(para_num2,5)*am0	!0.04591*am0
	am_aniso(1,2,2) = temp_read(para_num2,5)*am0	!0.04591*am0
	am_aniso(1,2,3) = temp_read(para_num2,6)*am0	!0.04591*am0
	am_aniso(2,2,1) = temp_read(para_num2,7)*am0	!0.17174*am0
	am_aniso(2,2,2) = temp_read(para_num2,7)*am0	!0.17174*am0
	am_aniso(2,2,3) = temp_read(para_num2,7)*am0	!0.17174*am0
	am_aniso(3,2,1) = temp_read(para_num2,10)*am0	!0.35054*am0
	am_aniso(3,2,2) = temp_read(para_num2,10)*am0	!0.35054*am0
	am_aniso(3,2,3) = temp_read(para_num2,10)*am0	!0.35054*am0

c
	hole_am_aniso(1,2,1) = temp_read(para_num2,30)*am0	!0.57162*am0
	hole_am_aniso(1,2,2) = temp_read(para_num2,30)*am0	!0.57162*am0
	hole_am_aniso(1,2,3) = temp_read(para_num2,31)*am0	!0.57162*am0

	hole_aff_aniso(1,2,1) = temp_read(para_num2,32)	!55.7159
	hole_aff_aniso(1,2,2) = temp_read(para_num2,32)	!55.7159
	hole_aff_aniso(1,2,3) = temp_read(para_num2,32)	!55.7159

c---( 非放物線性パラメータα )---
	aff(1,2) = temp_read(para_num2,21)		!1.450	!Γ谷
	aff(2,2) = temp_read(para_num2,22)		!0.466	!Ｌ谷
	aff(3,2) = temp_read(para_num2,24)		!0.133	!Ｘ谷

	aff_aniso(1,2,1) =  temp_read(para_num2,21)		!1.450
	aff_aniso(1,2,2) =  temp_read(para_num2,21)		!1.450
	aff_aniso(1,2,3) =  temp_read(para_num2,21)		!1.450
	aff_aniso(2,2,1) =  temp_read(para_num2,22)		!0.466
	aff_aniso(2,2,2) =  temp_read(para_num2,22)		!0.466
	aff_aniso(2,2,3) =  temp_read(para_num2,22)		!0.466
	aff_aniso(3,2,1) =  temp_read(para_num2,24)		!0.133
	aff_aniso(3,2,2) =  temp_read(para_num2,24)		!0.133
	aff_aniso(3,2,3) =  temp_read(para_num2,24)		!0.133
c
c---( 伝導帯底のエネルギー )---
	ec(1,2)	= 0.000		!Γ谷
	ec(2,2)	= temp_read(para_num2,2)			!0.738		!Ｌ谷
	ec(3,2)	= temp_read(para_num2,3)			!1.079		!Ｘ谷
c
c---( 音響フォノンの変位ポテンシャル (da(i) = Ξd))---
	da(1,2)	= (temp_material(material_num3,5)*inx(2)+
     &	temp_material(material_num4,5)*(1-inx(2)))*q
      !(5.887*inx(2)+7.00*(1-inx(2)))*q !5.887*q	!Γ谷
	da(2,2)	= (temp_material(material_num3,6)*inx(2)+
     &	temp_material(material_num4,6)*(1-inx(2)))*q
      !(10.80*inx(2)+9.20*(1-inx(2)))*q !10.80*q	!Ｌ谷
	da(3,2)	= (temp_material(material_num3,7)*inx(2)+
     &	temp_material(material_num4,7)*(1-inx(2)))*q
      !(9.657*inx(2)+9.27*(1-inx(2)))*q !9.657*q	!Ｘ谷
c
c---( 光学フォノンの変位ポテンシャル(eV/m) (d(i,j) = Ｄij))---	 
	d(1,1,2) = 0.0			!ΓtoΓ
	d(2,1,2) = (temp_material(material_num3,9)*inx(2)+
     &	temp_material(material_num4,9)*(1-inx(2)))*1.0e10*q
      !(7.83*inx(2)+5.25*(1-inx(2)))*1.0e10*q !7.83e10*q	!ΓtoＬ
	d(3,1,2) = (temp_material(material_num3,10)*inx(2)+
     &	temp_material(material_num4,10)*(1-inx(2)))*1.0e10*q
      !(11.32*inx(2)+5.48*(1-inx(2)))*1.0e10*q !11.32e10*q	!ΓtoＸ
	d(2,2,2) = (temp_material(material_num3,11)*inx(2)+
     &	temp_material(material_num4,11)*(1-inx(2)))*1.0e10*q
      !(6.40*inx(2)+5.94*(1-inx(2)))*1.0e10*q !6.40e10*q	!ＬtoＬ
	d(3,2,2) = (temp_material(material_num3,12)*inx(2)+
     &	temp_material(material_num4,12)*(1-inx(2)))*1.0e10*q
      !(6.80*inx(2)+5.01*(1-inx(2)))*1.0e10*q !6.80e10*q	!ＬtoＸ
	d(3,3,2) = (temp_material(material_num3,13)*inx(2)+
     &	temp_material(material_num4,13)*(1-inx(2)))*1.0e10*q
      !(8.54*inx(2)+2.99*(1-inx(2)))*1.0e10*q !8.54e10*q	!ＸtoＸ
c
c---(	有極性光学フォノンのエネルギー )---
	hwo(1:3,2) = temp_material(material_num3,14)*
     &	inx(2)+temp_material(material_num4,14)*(1-inx(2))
      !0.0328*inx(2)+0.03536*(1-inx(2)) !0.0328
c	hwo(1,2) = 0.0
c	hwo(2,2) = 0.0328
c	hwo(3,2) = 0.0
c
c---(	バンド端散乱のフォノンエネルギー(eV/m) (hw(i,j) = hwij))---
	hw(1,1,2) = 0.0			!ΓtoΓ
	hw(2,1,2) = (temp_material(material_num3,19)*inx(2)+
     &	temp_material(material_num4,19)*(1-inx(2)))
      !(25.4*inx(2)+22.69*(1-inx(2)))*1.0e-3 !0.0254		!ΓtoＬ
	hw(3,1,2) = (temp_material(material_num3,20)*inx(2)+
     &	temp_material(material_num4,20)*(1-inx(2)))
      !0.0257		!ΓtoＸ
	hw(2,2,2) = (temp_material(material_num3,21)*inx(2)+
     &	temp_material(material_num4,21)*(1-inx(2)))
      !(24.8*inx(2)+24.97*(1-inx(2)))*1.0e-3 !0.0248		!ＬtoＬ
	hw(3,2,2) = (temp_material(material_num3,22)*inx(2)+
     &	temp_material(material_num4,22)*(1-inx(2)))
      !(30.2*inx(2)+21.85*(1-inx(2)))*1.0e-3 !0.0302		!ＬtoＸ
	hw(3,3,2) = (temp_material(material_num3,23)*inx(2)+
     &	temp_material(material_num4,23)*(1-inx(2)))
      !(28.4*inx(2)+24.31*(1-inx(2)))*1.0e-3 !0.0284		!ＸtoＸ

c
c---(合金散乱)---
c	ea = 0.0		!合金散乱非考慮
c	ea(2) = temp_material(material_num3,24) !1.5		!08/11/10 竹岸
c	ea(2) = (temp_material(material_num3,24)*inx(2)/0.47+
c     &	temp_material(material_num4,24)*(1-inx(2))/0.47)
	ea(2) = (temp_material(material_num3,24)*inx(2)+
     &	temp_material(material_num4,24)*(1-inx(2)))	!非考慮
c	write(*,*) 1-inx(2)
c	write(*,*) 'ea(2)',ea(2)
c
c---( X谷、L谷の最大値_I.I.用)---
	hiXL(2,2) = temp_read(para_num2,33)
	hiXL(2,3) = temp_read(para_num2,34)

c---( 衝突電離 )---
c	eth = 1000		!衝突電離非考慮		
	eth(2) = temp_read(para_num2,28)		!0.808793		!JAP94(2003)4096
	a(2)	= temp_material(material_num3,25) !1.0e12
	II_S(2)	= a(2)		
	b(2)	= temp_material(material_num3,26) !2.0
c
c	eth = 0.80		!JAP94(2003)4096
c	a	= 2.0e12
c	b	= 2.0

c	eth = 0.75		!JAP96(2004)5650
c	a	= 2.0e10
c	b	= 2.5
c
c---------------------Al0.25In0.75Sb(old_InSb�@チャネル)-----------------
c---( In組成比 )---
	inx(3) = temp_read(para_num3,1)		!0.25		
	cx(3) = inx(3)*(1-inx(3))
c
c---( 格子定数,体積 )---
	lc(3) =temp_read(para_num3,29)*1.0e-10	
	 !temp_read(3,29)*1.0e-10		!5.872*1e-10       !格子定数
	va(3) = lc(3)**3/4.0
c
c---( 禁制帯幅 )---
	egmin(3) = temp_read(para_num3,27)	!0.472448  !InASのものを直ち
c
c---( フォノン散乱の諸パラメータ )---
	rou(3) = temp_material(material_num5,1)*
     &	inx(3)+temp_material(material_num6,1)*(1-inx(3))
     	!rou(3) = 5680*inx(3)+5310*(1-inx(3)) !5680	!半導体の比重(kg/m^3)
	sv(3)  = temp_material(material_num5,2)*
     &	inx(3)+temp_material(material_num6,2)*(1-inx(3))	!修正11/07/25原
     	!sv(3)  = 4280*inx(3)+5240*(1-inx(3)) !4280	!半導体中の音速(m/s)
c
c---( 誘電率 )---
	eps(3)	= (temp_material(material_num5,3)*
     &	inx(3)+temp_material(material_num6,3)*(1-inx(3)))*ep0
      !eps(3)	= (15.1*inx(3)+12.90*(1-inx(3)))*ep0 !15.1*ep0		!半導体の誘電率εs
	epf(3)  = (temp_material(material_num5,4)*
     &	inx(3)+temp_material(material_num6,4)*(1-inx(3)))*ep0
      !epf(3)  = (12.3*inx(3)+10.89*(1-inx(3)))*ep0 !12.3*ep0		!光学的誘電率ε∞
	ep(3)   = 1.0/(1.0/epf(3)-1.0/eps(3))
c
c---( 電子の有効質量 )---
	am(1,3)	= temp_read(para_num3,16)*am0		!0.042653035*am0	!Γ谷
	am(2,3)	= temp_read(para_num3,17)*am0		!0.1915490*am0	!Ｌ谷
	am(3,3)	= temp_read(para_num3,18)*am0		!0.3691148*am0	!Ｘ谷

	am_aniso(1,3,1) = temp_read(para_num3,5)*am0	!0.042653035*am0
	am_aniso(1,3,2) = temp_read(para_num3,5)*am0	!0.042653035*am0
	am_aniso(1,3,3) = temp_read(para_num3,6)*am0	!0.042653035*am0
	am_aniso(2,3,1) = temp_read(para_num3,7)*am0	!0.1915490*am0
	am_aniso(2,3,2) = temp_read(para_num3,7)*am0	!0.1915490*am0
	am_aniso(2,3,3) = temp_read(para_num3,7)*am0	!0.1915490*am0
	am_aniso(3,3,1) = temp_read(para_num3,10)*am0	!0.3691148*am0
	am_aniso(3,3,2) = temp_read(para_num3,10)*am0	!0.3691148*am0
	am_aniso(3,3,3) = temp_read(para_num3,10)*am0	!0.3691148*am0
c
	hole_am_aniso(1,3,1) = temp_read(para_num3,30)*am0	!0.57162*am0
	hole_am_aniso(1,3,2) = temp_read(para_num3,30)*am0	!0.57162*am0
	hole_am_aniso(1,3,3) = temp_read(para_num3,31)*am0	!0.57162*am0

	hole_aff_aniso(1,3,1) = temp_read(para_num3,32)	!55.7159
	hole_aff_aniso(1,3,2) = temp_read(para_num3,32)	!55.7159
	hole_aff_aniso(1,3,3) = temp_read(para_num3,32)	!55.7159

c---( 非放物線性パラメータα )---
	aff(1,3) = temp_read(para_num3,21)		!1.601	!Γ谷
	aff(2,3) = temp_read(para_num3,22)		!0.293	!Ｌ谷
	aff(3,3) = temp_read(para_num3,24)		!0.138	!Ｘ谷

	aff_aniso(1,3,1) =  temp_read(para_num3,21)		!1.601
	aff_aniso(1,3,2) =  temp_read(para_num3,21)		!1.601
	aff_aniso(1,3,3) =  temp_read(para_num3,21)		!1.601
	aff_aniso(2,3,1) =  temp_read(para_num3,22)		!0.293
	aff_aniso(2,3,2) =  temp_read(para_num3,22)		!0.293
	aff_aniso(2,3,3) =  temp_read(para_num3,22)		!0.293
	aff_aniso(3,3,1) =  temp_read(para_num3,24)		!0.138
	aff_aniso(3,3,2) =  temp_read(para_num3,24)		!0.138
	aff_aniso(3,3,3) =  temp_read(para_num3,24)		!0.138
c
c---( 伝導帯底のエネルギー )---
	ec(1,3)	= 0.000		!Γ谷
	ec(2,3)	= temp_read(para_num3,2)			!0.928	!Ｌ谷
	ec(3,3)	= temp_read(para_num3,3)			!1.318	!Ｘ谷
c
c---( 音響フォノンの変位ポテンシャル (da(i) = Ξd))---
	da(1,3)	= (temp_material(material_num5,5)*inx(3)+
     &	temp_material(material_num6,5)*(1-inx(3)))*q
      !da(1,3)	= (5.93*inx(3)+7.00*(1-inx(3)))*q !5.93*q	!Γ谷
	da(2,3)	= (temp_material(material_num5,6)*inx(3)+
     &	temp_material(material_num6,6)*(1-inx(3)))*q
      !da(2,3)	= (7.23*inx(3)+9.20*(1-inx(3)))*q !7.23*q	!L谷
	da(3,3)	= (temp_material(material_num5,7)*inx(2)+
     &	temp_material(material_num6,7)*(1-inx(2)))*q
      !da(3,3)	= (9.02*inx(3)+9.27*(1-inx(3)))*q !9.02*q	!Ｘ谷
c
c---( 光学フォノンの変位ポテンシャル(eV/m) (d(i,j) = Ｄij))---	 
	d(1,1,3) = 0.0										!ΓtoΓ
	d(2,1,3) = (temp_material(material_num5,9)*inx(3)+
     &	temp_material(material_num6,9)*(1-inx(3)))*1.0e10*q
      !d(2,1,3) = ((5.59*inx(3)+5.25*(1-inx(3)))*1.0e10*q !5.59*1.0e10*q	!ΓtoＬ
	d(3,1,3) = (temp_material(material_num5,10)*inx(3)+
     &	temp_material(material_num6,10)*(1-inx(3)))*1.0e10*q
      !d(3,1,3) = (6.35*inx(3)+5.48*(1-inx(3)))*1.0e10*q !6.35*1.0e10*q	!ΓtoＸ
	d(2,2,3) = (temp_material(material_num5,11)*inx(3)+
     &	temp_material(material_num6,11)*(1-inx(3)))*1.0e10*q
      !d(2,2,3) = (6.35*inx(3)+5.94*(1-inx(3)))*1.0e10*q !6.35*1.0e10*q	!ＬtoＬ
	d(3,2,3) = (temp_material(material_num5,12)*inx(3)+
     &	temp_material(material_num6,12)*(1-inx(3)))*1.0e10*q
      !d(3,2,3) = (5.59*inx(3)+5.01*(1-inx(3)))*1.0e10*q !5.59*1.0e10*q	!ＬtoＸ
	d(3,3,3) = (temp_material(material_num5,13)*inx(3)+
     &	temp_material(material_num6,13)*(1-inx(3)))*1.0e10*q
      !d(3,3,3) = (3.36*inx(3)+2.99*(1-inx(3)))*1.0e10*q !3.36*1.0e10*q	!ＸtoＸ
c
c---(	有極性光学フォノンのエネルギー )---
	hwo(1:3,3) = temp_material(material_num5,14)*
     &	inx(3)+temp_material(material_num6,14)*(1-inx(3))
      !hwo(1:3,3) = 0.0302*inx(3)+0.03536*(1-inx(3))
c	hwo(1,2) = 0.0
c	hwo(2,2) = 0.0328
c	hwo(3,2) = 0.0
c
c---(	バンド端散乱のフォノンエネルギー(eV/m) (hw(i,j) = hwij))---
	hw(1,1,3) = 0.0										!ΓtoΓ
	hw(2,1,3) = (temp_material(material_num5,19)*inx(3)+
     &	temp_material(material_num6,19)*(1-inx(3)))
      !hw(2,1,3) = (17.45*inx(3)+22.69*(1-inx(3)))*1.0e-3 !17.45*1.0e-3	!ΓtoＬ
	hw(3,1,3) = (temp_material(material_num5,20)*inx(3)+
     &	temp_material(material_num6,20)*(1-inx(3)))
      !hw(3,1,3) = (19.23*inx(3)+23.45*(1-inx(3)))*1.0e-3 !19.23*1.0e-3	!ΓtoＸ
	hw(2,2,3) = (temp_material(material_num5,21)*inx(3)+
     &	temp_material(material_num6,21)*(1-inx(3)))
      !hw(2,2,3) = (19.23*inx(3)+24.97*(1-inx(3)))*1.0e-3 !19.23*1.0e-3	!ＬtoＬ
	hw(3,2,3) = (temp_material(material_num5,22)*inx(3)+
     &	temp_material(material_num6,22)*(1-inx(3)))
      !hw(3,2,3) = (17.45*inx(3)+21.85*(1-inx(3)))*1.0e-3 !17.45*1.0e-3	!ＬtoＸ
	hw(3,3,3) = (temp_material(material_num5,23)*inx(3)+
     &	temp_material(material_num6,23)*(1-inx(3)))
      !hw(3,3,3) = (19.26*inx(3)+24.31*(1-inx(3)))*1.0e-3 !19.26*1.0e-3	!ＸtoＸ
c
c---(合金散乱)---
c	ea = 0.0		!合金散乱非考慮
c	ea(3) = temp_material(material_num5,24) !ea(3) = 0.0		!08/11/10 竹岸 InAsなので0.0
c	ea(3) = (temp_material(material_num5,24)*inx(3)/0.47+
c     &	temp_material(material_num6,24)*(1-inx(3))/0.47)
	ea(3) = (temp_material(material_num5,24)*inx(3)+	!11/07/25原 0.47修正
     &	temp_material(material_num6,24)*(1-inx(3)))
c	write(*,*) material_num5
c	write(*,*) material_num6
c	write(*,*) temp_material(material_num5,24)
c	write(*,*) temp_material(material_num6,24)
c	write(*,*) 'ea(3)',ea(3)

c
c---( X谷、L谷の最大値_I.I.用)---
	hiXL(3,2) = temp_read(para_num3,33)
	hiXL(3,3) = temp_read(para_num3,34)

c---( 衝突電離 )---
c	eth = 1000		!衝突電離非考慮
	eth(3) = temp_read(para_num3,28)		!0.5179486
	a(3)	= temp_material(material_num5,25) !a(3)   = 1.0e12
	II_S(3)	= a(3)
	b(3)	= temp_material(material_num5,26) !b(3)   = 2.0
c
c	eth = 0.80		!JAP94(2003)4096
c	a	= 2.0e12
c	b	= 2.0

c	eth = 0.75		!JAP96(2004)5650
c	a	= 2.0e10
c	b	= 2.5
c
c---------------------AlSb(old_InP)---------------------
c
c---( In組成比 )---
	inx(4) = temp_read(para_num4,1)		!1
	cx(4) = inx(4)*(1-inx(4))
c---( 禁制帯幅 )---
	egmin(4) = temp_read(para_num4,27)	!1.344
c---( 格子定数,体積 )---
	lc(4) = temp_read(para_num4,29)*1.0e-10		!5.8687e-10
	va(4) = lc(4)**3/4.0
c
c
c---( フォノン散乱の諸パラメータ )---
	rou(4) = temp_material(material_num7,1)*
     &	inx(4)+temp_material(material_num8,1)*(1-inx(4))
     	!rou(4) = 4810	!半導体の比重(kg/m^3)
	sv(4)  = temp_material(material_num7,2)*
     &	inx(4)+temp_material(material_num8,2)*(1-inx(4))	!修正11/07/25原
     	!sv(4)  = 5130	!半導体中の音速(m/s)
c
c---( 誘電率 )---
	eps(4)	= (temp_material(material_num7,3)*
     &	inx(4)+temp_material(material_num8,3)*(1-inx(4)))*ep0
      !eps(4)	= 12.56*ep0		!半導体の誘電率εs
	epf(4)  = (temp_material(material_num7,4)*
     &	inx(4)+temp_material(material_num8,4)*(1-inx(4)))*ep0
      !epf(4)	= 9.61*ep0		!光学的誘電率ε∞
	ep(4)   = 1.0/(1.0/epf(4)-1.0/eps(4))
c
c---( 電子の有効質量 )---
	am(1,4)	= temp_read(para_num4,16)*am0		!0.08*am0	!Γ谷
	am(2,4)	= temp_read(para_num4,17)*am0		!0.25*am0	!Ｌ谷
	am(3,4)	= temp_read(para_num4,18)*am0		!0.325*am0	!Ｘ谷

	am_aniso(1,4,1) = temp_read(para_num4,5)*am0	!0.042653035*am0
	am_aniso(1,4,2) = temp_read(para_num4,5)*am0	!0.042653035*am0
	am_aniso(1,4,3) = temp_read(para_num4,6)*am0	!0.042653035*am0
	am_aniso(2,4,1) = temp_read(para_num4,7)*am0	!0.1915490*am0
	am_aniso(2,4,2) = temp_read(para_num4,7)*am0	!0.1915490*am0
	am_aniso(2,4,3) = temp_read(para_num4,7)*am0	!0.1915490*am0
	am_aniso(3,4,1) = temp_read(para_num4,10)*am0	!0.3691148*am0
	am_aniso(3,4,2) = temp_read(para_num4,10)*am0	!0.3691148*am0
	am_aniso(3,4,3) = temp_read(para_num4,10)*am0	!0.3691148*am0

c
	hole_am_aniso(1,4,1) = temp_read(para_num4,30)*am0	!0.57162*am0
	hole_am_aniso(1,4,2) = temp_read(para_num4,30)*am0	!0.57162*am0
	hole_am_aniso(1,4,3) = temp_read(para_num4,31)*am0	!0.57162*am0

	hole_aff_aniso(1,4,1) = temp_read(para_num4,32)	!55.7159
	hole_aff_aniso(1,4,2) = temp_read(para_num4,32)	!55.7159
	hole_aff_aniso(1,4,3) = temp_read(para_num4,32)	!55.7159

c---( 非放物線性パラメータα )---
	aff(1,4) = temp_read(para_num4,21)		!0.83	!Γ谷
	aff(2,4) = temp_read(para_num4,22)		!0.23	!Ｌ谷
	aff(3,4) = temp_read(para_num4,24)		!0.38	!Ｘ谷

	aff_aniso(1,4,1) =  temp_read(para_num4,21)		!0.83
	aff_aniso(1,4,2) =  temp_read(para_num4,21)		!0.83
	aff_aniso(1,4,3) =  temp_read(para_num4,21)		!0.83
	aff_aniso(2,4,1) =  temp_read(para_num4,22)		!0.23
	aff_aniso(2,4,2) =  temp_read(para_num4,22)		!0.23
	aff_aniso(2,4,3) =  temp_read(para_num4,22)		!0.23
	aff_aniso(3,4,1) =  temp_read(para_num4,24)		!0.38
	aff_aniso(3,4,2) =  temp_read(para_num4,24)		!0.38
	aff_aniso(3,4,3) =  temp_read(para_num4,24)		!0.38
c
c---( 伝導帯底のエネルギー )---
	ec(1,4)	= 0.0		!Γ谷
	ec(2,4)	= temp_read(para_num4,2)			!0.540		!Ｌ谷
	ec(3,4)	= temp_read(para_num4,3)			!0.775		!Ｘ谷
c
c---( 音響フォノンの変位ポテンシャル (da(i) = Ξd))---
	da(1,4)	= (temp_material(material_num7,5)*inx(4)+
     &	temp_material(material_num8,5)*(1-inx(4)))*q
      !da(1,4)	= 5.00*q	!Γ谷
	da(2,4)	= (temp_material(material_num7,6)*inx(4)+
     &	temp_material(material_num8,6)*(1-inx(4)))*q
      !da(2,4)	= 5.00*q	!Ｌ谷
	da(3,4)	= (temp_material(material_num7,7)*inx(4)+
     &	temp_material(material_num8,7)*(1-inx(4)))*q
      !da(3,4)	= 5.00*q	!Ｘ谷
c
c---( 光学フォノンの変位ポテンシャル(eV/m) (d(i,j) = Ｄij))---	 
	d(1,1,4) = 0.0			!ΓtoΓ
	d(2,1,4) = (temp_material(material_num7,9)*inx(4)+
     &	temp_material(material_num8,9)*(1-inx(4)))*1.0e10*q
      !d(2,1,4) = 5.06e10*q	!ΓtoＬ
	d(3,1,4) = (temp_material(material_num7,10)*inx(4)+
     &	temp_material(material_num8,10)*(1-inx(4)))*1.0e10*q
      !d(3,1,4) = 4.98e10*q	!ΓtoＸ
	d(2,2,4) = (temp_material(material_num7,11)*inx(4)+
     &	temp_material(material_num8,11)*(1-inx(4)))*1.0e10*q
      !d(2,2,4) = 5.75e10*q	!ＬtoＬ
	d(3,2,4) = (temp_material(material_num7,12)*inx(4)+
     &	temp_material(material_num8,12)*(1-inx(4)))*1.0e10*q
      !d(3,2,4) = 4.68e10*q	!ＬtoＸ
	d(3,3,4) = (temp_material(material_num7,13)*inx(4)+
     &	temp_material(material_num8,13)*(1-inx(4)))*1.0e10*q
      !d(3,3,4) = 2.80e10*q	!ＸtoＸ
c
c---(	有極性光学フォノンのエネルギー )---
	hwo(1:3,4) = temp_material(material_num7,14)*
     &	inx(4)+temp_material(material_num8,14)*(1-inx(4))
      !hwo(1:3,4) = 0.0430
c	hwo(1,4) = 0.0
c	hwo(2,4) = 0.0404
c	hwo(3,4) = 0.0
c
c---(	バンド端散乱のフォノンエネルギー(eV/m) (hw(i,j) = hwij))---
	hw(1,1,4) = 0.0		!ΓtoΓ
	hw(2,1,4) = (temp_material(material_num7,19)*inx(4)+
     &	temp_material(material_num8,19)*(1-inx(4)))
      !hw(2,1,4) = 0.0278	!ΓtoＬ
	hw(3,1,4) = (temp_material(material_num7,20)*inx(4)+
     &	temp_material(material_num8,20)*(1-inx(4)))
      !hw(3,1,4) = 0.0299	!ΓtoＸ
	hw(2,2,4) = (temp_material(material_num7,21)*inx(4)+
     &	temp_material(material_num8,21)*(1-inx(4)))
      !hw(2,2,4) = 0.029	!ＬtoＬ
	hw(3,2,4) = (temp_material(material_num7,22)*inx(4)+
     &	temp_material(material_num8,22)*(1-inx(4)))
      !hw(3,2,4) = 0.0293	!ＬtoＸ
	hw(3,3,4) = (temp_material(material_num7,23)*inx(4)+
     &	temp_material(material_num8,23)*(1-inx(4)))
      !hw(3,3,4) = 0.0299	!ＸtoＸ
c
c

c

c---(合金散乱)---
c	ea = 0.0		!合金散乱非考慮
c	ea(4) = temp_material(material_num7,24) !ea(4) = 0.0
	ea(4) = (temp_material(material_num7,24)*inx(1)+
     &	temp_material(material_num8,24)*(1-inx(1)))
c
c---( X谷、L谷の最大値_I.I.用)---
	hiXL(4,2) = temp_read(para_num4,33)
	hiXL(4,3) = temp_read(para_num4,34)

c---( 衝突電離 )---
c	eth = 1000		!衝突電離非考慮		!ここを変更すればよい(2006/12/09 Hara)!材料4を使わないときは上を変えれば良い
	eth(4) = temp_read(4,28)		!1.69		!JAP94(2003)4096
	a(4)	= temp_material(material_num7,25) !a(4) = 1.5e14
	II_S(4)	= a(4)		
	b(4)	= temp_material(material_num7,26) !b(4) = 2.0
c
c	eth = 0.80		!JAP94(2003)4096
c	a	= 2.0e12
c	b	= 2.0

c	eth = 0.75		!JAP96(2004)5650
c	a	= 2.0e10
c	b	= 2.5
c
c---------------------InSb NotUse( In0.53Ga0.47As�AチャネルTop and Bottom Layer)---------

	inx(5) = temp_read(para_num5,1)	
	cx(5) = inx(5)*(1-inx(5))
c	write(*,*) inx(5)
c
c---( 格子定数,体積 )---
	lc(5) =(temp_read(para_num5,29)*inx(5)+5.65325*(1-inx(5)))*1.0e-10
	 !temp_read(5,29)*1.0e-10		!5.872*1e-10       !格子定数
	va(5) = lc(5)**3/4.0

c---( 禁制帯幅 )---
	egmin(5) = temp_read(para_num5,27)	!0.675
c
c---( フォノン散乱の諸パラメータ )---
	rou(5) = temp_material(material_num9,1)*
     &	inx(5)+temp_material(material_num10,1)*(1-inx(5))
     	!rou(5) = 5469	!半導体の比重(kg/m^3)
	sv(5)  = temp_material(material_num9,2)*
     &	inx(5)+temp_material(material_num10,2)*(1-inx(5))	!修正11/07/25原
     	!sv(5)  = 4742	!半導体中の音速(m/s)
c
c---( 誘電率 )---
	eps(5)	= (temp_material(material_num9,3)*
     &	inx(5)+temp_material(material_num10,3)*(1-inx(5)))*ep0
      !eps(5)	= 13.88*ep0		!半導体の誘電率εs
	epf(5)  = (temp_material(material_num9,4)*
     &	inx(5)+temp_material(material_num10,4)*(1-inx(5)))*ep0
      !epf(5)  = 11.34*ep0		!光学的誘電率ε∞
	ep(5)   = 1.0/(1.0/epf(5)-1.0/eps(5))
c
c---( 電子の有効質量 )---
	am(1,5)	= temp_read(para_num5,16)*am0		!0.0459*am0	!Γ谷
	am(2,5)	= temp_read(para_num5,17)*am0		!0.1717*am0	!Ｌ谷
	am(3,5)	= temp_read(para_num5,18)*am0		!0.3551*am0	!Ｘ谷

	am_aniso(1,5,1) = temp_read(para_num5,5)*am0	!0.0459*am0
	am_aniso(1,5,2) = temp_read(para_num5,5)*am0	!0.0459*am0
	am_aniso(1,5,3) = temp_read(para_num5,6)*am0	!0.0459*am0
	am_aniso(2,5,1) = temp_read(para_num5,7)*am0	!0.1717*am0
	am_aniso(2,5,2) = temp_read(para_num5,7)*am0	!0.1717*am0
	am_aniso(2,5,3) = temp_read(para_num5,7)*am0	!0.1717*am0
	am_aniso(3,5,1) = temp_read(para_num5,10)*am0	!0.3551*am0
	am_aniso(3,5,2) = temp_read(para_num5,10)*am0	!0.3551*am0
	am_aniso(3,5,3) = temp_read(para_num5,10)*am0	!0.3551*am0

c
c---( 非放物線性パラメータα )---
	aff(1,5) = temp_read(para_num5,21)		!1.450	!Γ谷
	aff(2,5) = temp_read(para_num5,22)		!0.466	!Ｌ谷
	aff(3,5) = temp_read(para_num5,24)		!.133	!Ｘ谷

	aff_aniso(1,5,1) =  temp_read(para_num5,21)	!1.450
	aff_aniso(1,5,2) =  temp_read(para_num5,21)	!1.450
	aff_aniso(1,5,3) =  temp_read(para_num5,21)	!1.450
	aff_aniso(2,5,1) =  temp_read(para_num5,22)	!0.466
	aff_aniso(2,5,2) =  temp_read(para_num5,22)	!0.466
	aff_aniso(2,5,3) =  temp_read(para_num5,22)	!0.466
	aff_aniso(3,5,1) =  temp_read(para_num5,24)	!0.133
	aff_aniso(3,5,2) =  temp_read(para_num5,24)	!0.133
	aff_aniso(3,5,3) =  temp_read(para_num5,24)	!0.133
c
	hole_am_aniso(1,5,1) = temp_read(para_num5,30)*am0	!0.57162*am0
	hole_am_aniso(1,5,2) = temp_read(para_num5,30)*am0	!0.57162*am0
	hole_am_aniso(1,5,3) = temp_read(para_num5,31)*am0	!0.57162*am0
c
	hole_aff_aniso(1,5,1) = temp_read(para_num5,32)	!55.7159
	hole_aff_aniso(1,5,2) = temp_read(para_num5,32)	!55.7159
	hole_aff_aniso(1,5,3) = temp_read(para_num5,32)	!55.7159

c---( 伝導帯底のエネルギー )---
	ec(1,5)	= 0.000		!Γ谷
	ec(2,5)	= temp_read(para_num5,2)			!0.738		!Ｌ谷
	ec(3,5)	= temp_read(para_num5,3)			!1.079		!Ｘ谷
c
c---( 音響フォノンの変位ポテンシャル (da(i) = Ξd))---
	da(1,5)	= (temp_material(material_num9,5)*inx(5)+
     &	temp_material(material_num10,5)*(1-inx(5)))*q
      !da(1,5)	= 5.887*q	!Γ谷
	da(2,5)	= (temp_material(material_num9,6)*inx(5)+
     &	temp_material(material_num10,6)*(1-inx(5)))*q
      !da(2,5)	= 10.80*q	!Ｌ谷
	da(3,5)	= (temp_material(material_num9,7)*inx(5)+
     &	temp_material(material_num10,7)*(1-inx(5)))*q
      !da(3,5)	= 9.657*q	!Ｘ谷
c
c---( 光学フォノンの変位ポテンシャル(eV/m) (d(i,j) = Ｄij))---	 
	d(1,1,5) = 0.0			!ΓtoΓ
	d(2,1,5) = (temp_material(material_num9,9)*inx(5)+
     &	temp_material(material_num10,9)*(1-inx(5)))*1.0e10*q
c      d(2,1,5) = 7.83e10*q	!ΓtoＬ
	d(3,1,5) = (temp_material(material_num9,10)*inx(5)+
     &	temp_material(material_num10,10)*(1-inx(5)))*1.0e10*q
c      d(3,1,5) = 11.32e10*q	!ΓtoＸ
	d(2,2,5) = (temp_material(material_num9,11)*inx(5)+
     &	temp_material(material_num10,11)*(1-inx(5)))*1.0e10*q
c      d(2,2,5) = 6.40e10*q	!ＬtoＬ
	d(3,2,5) = (temp_material(material_num9,12)*inx(5)+
     &	temp_material(material_num10,12)*(1-inx(5)))*1.0e10*q
c      d(3,2,5) = 6.80e10*q	!ＬtoＸ
	d(3,3,5) = (temp_material(material_num9,13)*inx(5)+
     &	temp_material(material_num10,13)*(1-inx(5)))*1.0e10*q
c      d(3,3,5) = 8.54e10*q	!ＸtoＸ
c
c---(	有極性光学フォノンのエネルギー )---
	hwo(1:3,5) = temp_material(material_num9,14)*
     &	inx(5)+temp_material(material_num10,14)*(1-inx(5))
c      hwo(1:3,5) = 0.0328
c	hwo(1,5) = 0.0
c	hwo(2,5) = 0.0328
c	hwo(3,5) = 0.0
c
c---(	バンド端散乱のフォノンエネルギー(eV/m) (hw(i,j) = hwij))---
	hw(1,1,5) = 0.0			!ΓtoΓ
	hw(2,1,5) = (temp_material(material_num9,19)*inx(5)+
     &	temp_material(material_num10,19)*(1-inx(5)))
c      hw(2,1,5) = 0.0254		!ΓtoＬ
	hw(3,1,5) = (temp_material(material_num9,20)*inx(5)+
     &	temp_material(material_num10,20)*(1-inx(5)))
c      hw(3,1,5) = 0.0257		!ΓtoＸ
	hw(2,2,5) = (temp_material(material_num9,21)*inx(5)+
     &	temp_material(material_num10,21)*(1-inx(5)))
c      hw(2,2,5) = 0.0248		!ＬtoＬ
	hw(3,2,5) = (temp_material(material_num9,22)*inx(5)+
     &	temp_material(material_num10,22)*(1-inx(5)))
c      hw(3,2,5) = 0.0302		!ＬtoＸ
	hw(3,3,5) = (temp_material(material_num9,23)*inx(5)+
     &	temp_material(material_num10,23)*(1-inx(5)))
c      hw(3,3,5) = 0.0284		!ＸtoＸ
c
c
c---( 格子定数,体積 )---
	lc(5) = temp_read(para_num5,29)*1.0e-10		!5.8687e-10
	va(5) = lc(5)**3/4.0
c
c---(合金散乱)---
c	ea = 0.0		!合金散乱非考慮
c	ea(5) = temp_material(material_num9,24) !ea(5) = 1.5		!08/11/10 竹岸
c	ea(5) = (temp_material(material_num9,24)*inx(5)/0.47+
c     &	temp_material(material_num10,24)*(1-inx(5))/0.47)
	ea(5) = (temp_material(material_num9,24)*inx(5)+	!11/07/25 原
     &	temp_material(material_num10,24)*(1-inx(5)))
c
c---( X谷、L谷の最大値_I.I.用)---
	hiXL(5,2) = temp_read(para_num5,33)
	hiXL(5,3) = temp_read(para_num5,34)

c---( 衝突電離 )---
c	eth = 1000		!衝突電離非考慮		
	eth(5) = temp_read(para_num5,28)		!0.808793		!JAP94(2003)4096
	a(5)	= temp_material(material_num9,25) !a(5)	= 1.0e12
	II_S(5)	= a(5)		
	b(5)	= temp_material(material_num9,26) !b(5)	= 2.0
c
c	eth = 0.80		!JAP94(2003)4096
c	a	= 2.0e12
c	b	= 2.0

c	eth = 0.75		!JAP96(2004)5650
c	a	= 2.0e10
c	b	= 2.5
c
c
c---(conduction band offsets by Nextnano)---	2017/12/1 鈴木貴博
	!dltecは基準材料から各材料へのΔEcです。ΓΓ,LLなど各バレー間の値
	!（多分バッファか何かに設定するのがデフォだと思いますが）

	dltec(1,1)=temp_read(para_num1,35)
	dltec(2,1)=temp_read(para_num1,36)
	dltec(3,1)=temp_read(para_num1,37)

	dltec(1,2)=temp_read(para_num2,35)
	dltec(2,2)=temp_read(para_num2,36)
	dltec(3,2)=temp_read(para_num2,37)

	dltec(1,3)=temp_read(para_num3,35)
	dltec(2,3)=temp_read(para_num3,36)
	dltec(3,3)=temp_read(para_num3,37)

	dltec(1,4)=temp_read(para_num4,35)
	dltec(2,4)=temp_read(para_num4,36)
	dltec(3,4)=temp_read(para_num4,37)

	dltec(1,5)=temp_read(para_num5,35)
	dltec(2,5)=temp_read(para_num5,36)
	dltec(3,5)=temp_read(para_num5,37)

	close(199)
	close(200)
	close(201)
	close(202)
	close(203)
	close(204)
	close(205)
	close(206)
	close(207)
	close(208)
	close(209)

c-----パラメータの読み込みの確認100826------
	open(600,file='test_para.txt')
	do roop_num=1,5
		write(600,*) 'inx(',roop_num,')=',inx(roop_num)
		write(600,*) 'lc(',roop_num,')=',lc(roop_num)
		write(600,*) 'egmin(',roop_num,')=',egmin(roop_num)
		write(600,*) 'am(1,',roop_num,')=',am(1,roop_num)/am0
		write(600,*) 'am(2,',roop_num,')=',am(2,roop_num)/am0
		write(600,*) 'am(3,',roop_num,')=',am(3,roop_num)/am0
		write(600,*) 'aff(1,',roop_num,')=',aff(1,roop_num)
		write(600,*) 'aff(2,',roop_num,')=',aff(2,roop_num)
		write(600,*) 'aff(3,',roop_num,')=',aff(3,roop_num)
		write(600,*) 'ec(1,',roop_num,')=',ec(1,roop_num)
		write(600,*) 'ec(2,',roop_num,')=',ec(2,roop_num)
		write(600,*) 'ec(3,',roop_num,')=',ec(3,roop_num)
		write(600,*) 'ea(,',roop_num,')=',ea(roop_num)
		write(600,*) 'eth(,',roop_num,')=',eth(roop_num)
	
	end do
	write(*,*) '次いきます'
c	read(*,*) roop_num2
	close(600)
c
c
c-----チャネルの不純物散乱のレートを0にするための小細工----c
c	allocate (dn3(npart))	!09/2/19 不純物散乱 川端→竹岸修正
		dn3(1) = dconc(1)
		dn3(2) = dconc(2)
		dn3(3) = dconc(3)
		dn3(4) = dconc(4)
		dn3(5) = 0
		dn3(6) = dconc(6)	
		dn3(7) = dconc(7)	
		dn3(8) = dconc(8)	
		dn3(9) = dconc(9)
		dn3(10) = dconc(10)
		dn3(11) = dconc(11)	
c
c------------------------------------------------------------
c	narea==2の材料がチャネルであるという前提でつけたし
c　　　20070502 よこ
c　　　チャネル内の電子濃度によるあたらな散乱レートをつくる
c		dn1	   = minval(dconc)	!undope時の不純物濃度
c		dn2	   = 0.0e21			!電子濃度
c		dn3(1) = 1.0e21			!以下チャネル内電子濃度のパラメータ9つ
c		dn3(2) = 5.0e21			!濃度は手打ち
c		dn3(3) = 1.0e22 
c		dn3(4) = 5.0e22
c		dn3(5) = 1.0e23
c		dn3(6) = 5.0e23
c		dn3(7) = 1.0e24
c		dn3(8) = 5.0e24
c		dn3(9) = 1.0e25
c
c	do 10 ia = 1,narea	
	do 10 ipart = 1,npart	!07/8/4 不純物散乱
cc---		ia:1->In(0.52)Al(0.48)As, 2->In(0.53)Ga(0.47)As, 3->InxGa(1-x)Asを表す
cc---		ia:4->InPを表す(2006/12/09)
c
c---		ia:1->InSb, 2->Al(0.85)In(0.15)Sb, 3->Al(0.75)In(0.25)Sb, 4->AlSb 
c--		!dn1:不純物散乱レートを計算する時の不純物濃度	 --
c--		!通常チャネル領域の濃度を与えておく				 --
c		電極の散乱レートは不要
c------------ !07/8/4 不純物散乱--------------------
		if((ipart.eq.(npart-1)).or.(ipart.eq.(npart-2)))cycle !ipart=9,10=電極
c		if((ipart.eq.3).or.(ipart.eq.4)
c     &				.or.(ipart.eq.5).or.(ipart.eq.11))then
		if(ipart.eq.5)then		!cap & channel層
			ia=1
c		elseif(ipart.eq.1)then					
		elseif((ipart.eq.1).or.(ipart.eq.2).or.(ipart.eq.3)			!AlInSb
     &		.or.(ipart.eq.4).or.(ipart.eq.6)
     &		.or.(ipart.eq.7).or.(ipart.eq.8)
     &		.or.(ipart.eq.11))then						
			ia=2
c		elseif(ipart.eq.7)then
c			ia=3
c		elseif(ipart.eq.2)then
c			ia=4
c	    elseif((ipart.eq.6).or.(ipart.eq.8))then
c			ia=5
		endif
c-------------------------------------------------------
	do iv=2,nvalley
	do jv=1,iv-1
		d(jv,iv,ia)	= d(iv,jv,ia)
		hw(jv,iv,ia)	= hw(iv,jv,ia)
	enddo
	enddo
c
	do 20 itp = 1,ntenum			!各温度
c
	temp = btmp + dtmp*(float(itp)-0.5)
c
	!qeps   = q/eps*dx*dx	!？
	cl(ia)     = rou(ia)*sv(ia)*sv(ia)
      bktq(itp)   = bk*temp/q		!ボルツマンファクター	
      !qh     = q/h
c
	do iv = 1, nvalley
		eg(iv,ia)	= ec(iv,ia)+egmin(ia)
		af(iv,ia)  = aff(iv,ia)
		af2(iv,ia)	= 2.0*aff(iv,ia)
		af4(iv,ia)	= 4.0*aff(iv,ia)
		smh(iv,ia)	= sqrt(2.0*am(iv,ia))*sqrt(q)/h
		hhm(iv,ia)	= h/am(iv,ia)/q*h/2.
		hm(iv,ia)	= h/am(iv,ia)
c
		dos(iv)	= (sqrt(2.0*am(iv,ia))*sqrt(q)/h)**3/4.0/pi/pi
c	有極性光学フォノン散乱パラメータ
		wo(iv)	= hwo(iv,ia)*q/h
		no(iv)	= 1.0/(exp(hwo(iv,ia)/bktq(itp))-1.0)
		poe(iv)	= q/8.0/pi/ep(ia)*q*wo(iv)*(no(iv)+1.0)	!放出
		poa(iv)	= poe(iv)*no(iv)/(1.0+no(iv))			!吸収
c	音響フォノン散乱パラメータ
		aco(iv)	= 2.0*pi*da(iv,ia)/q*da(iv,ia)*bktq(itp)/h*q/cl(ia)
c	非有極性光学フォノン・バンド間散乱パラメータ
		do jv=1,nvalley
			w(jv,iv)	= hw(jv,iv,ia)*q/h
			if(hw(jv,iv,ia).eq.0) then
				n(jv,iv) = 0.0
				ope(jv,iv)= 0.0
				opa(jv,iv)= ope(jv,iv)*n(jv,iv)/(1.0+n(jv,iv))
			else
				n(jv,iv) = 1.0/(exp(hw(jv,iv,ia)/bktq(itp))-1.0)
				ope(jv,iv)= pi*d(jv,iv,ia)/w(jv,iv)*d(jv,iv,ia)
     &					/rou(ia)/q*(n(jv,iv)+1.0)
				opa(jv,iv)= ope(jv,iv)*n(jv,iv)/(1.0+n(jv,iv))
			endif
		enddo
	enddo
c
c
c---( 不純物散乱のためのパラメータ )---
c      qd21   = q*dn1/bktq(itp)/eps(ia)
c      bimp   = 2.0*pi*dn1*q*q/h*q/eps(ia)/eps(ia)
	qd21   = q*dconc(ipart)/bktq(itp)/eps(ia)		!07/8/1不純物散乱
      bimp   = 2.0*pi*dn3(ipart)*q*q/h*q/eps(ia)/eps(ia) !07/8/1不純物散乱
c---( 合金散乱のためのパラメータ )---
c	alloy(1:nvalley)=pi**3/h*va(ia)*ea(ia)**2*cx(ia)*dos(1:nvalley)*q	!!半導体デバイスシミュレーション
c	alloy(1:nvalley)=(3.0/8.0)*pi**3/h*va(ia)*ea(ia)**2*cx(ia)
c     &												*dos(1:nvalley)*q	!!APL61(1992)1202
c	alloy(1:nvalley)=2.0*pi/h*va(ia)*ea(ia)**2*cx(ia)
c     &												*dos(1:nvalley)*q	!!半導体物理
	alloy(1:nvalley)=((3*(pi**3.0))/(16*h))*va(ia)*ea(ia)**2*cx(ia)
     &												*dos(1:nvalley)*q	!!Gonzalez ref IEEE TED Vol.38 No.3(1991) 08/11/10 竹岸
c




c---( 散乱レート計算アルゴリズム )---
c	===================================================================================================
	do 30 ien=1,nenergy
	do 40 iv=1,nvalley
	do 50 ie=1,nemax
      ei=de(ien)*(dble(ie))	!-0.5)
      sei=sqrt(ei)
c
	call getswk(
     &		  iv,ei,sei,qd21,bimp,
     &		  af(1,ia),smh(1,ia),ec(1,ia),hwo(1,ia),hw(1,1,ia),
     &		  dos,poa,poe,aco,ope,opa,z,
c     &		  swk(1,ie,iv,ien,itp,ia),
     &		  swk(1,ie,iv,ien,itp,ipart),	!07/8/4 不純物散乱
     &		  escat(1,iv,ia),iband(1,iv,ia),iarg(1,ia),alloy,
     &		  eth(ia),a(ia),b(ia))
   50 continue
   40 continue
   30 continue
   20 continue
c
	call param_heat(hwo,hw,hescat(1,1,ia))
c
   10	continue
c

	deallocate (wo,da,no,w,d,n)
	deallocate (z,dos,poe,poa,aco,ope,opa)
	deallocate (alloy,inx,cx,lc,va,ea)
	deallocate (eth,a,b)
c
c---( 散乱レートのファイル出力 )---
      open(unit=15,file='swkpara.txt')
      open(unit=16,file='swk.txt')
	write(15,'(I)') nscat,nemax,nvalley,nenergy
	write(15,'(f)') de(1:nenergy)
	do ie=1,nemax
c		write(16,'(e)') swk(1:nscat,ie,1:nvalley,1:nenergy,1,7)
		write(16,'(e)') swk(1:nscat,ie,1:nvalley,1:nenergy,1,nchannel2)
	enddo
	close(15)
	close(16)
c
c---( 散乱レートの総和の計算 )---
	do iscat=2,nscat
c	do ia = 1,narea
	do ipart = 1,npart			!07/8/4 不純物散乱			
	do itp = 1,ntenum			!各温度
	do ien=1,nenergy
	do iv=1,nvalley
      do ie=1,nemax

c		swk(iscat,ie,iv,ien,itp,ia)=swk(iscat,ie,iv,ien,itp,ia)
c    &							+swk(iscat-1,ie,iv,ien,itp,ia)
		swk(iscat,ie,iv,ien,itp,ipart)=swk(iscat,ie,iv,ien,itp,ipart)
     &							+swk(iscat-1,ie,iv,ien,itp,ipart)	!07/8/4 不純物散乱	
	enddo
	enddo
	enddo
	enddo
	enddo
	enddo
c
c	allocate(gm(nvalley,nenergy,narea))
	allocate(gm(nvalley,nenergy,npart))	!07/8/4 不純物散乱	
	gm=0.0
	pgm = huge(pgm)
	nemin = int(nemax*0.01)+1
c	do ia = 1,narea
	do ipart = 1,npart	!07/8/4 不純物散乱	
	do ien = 1,nenergy
	do iv  = 1,nvalley
c		gm(iv,ien,ia)=
c     &	maxval(swk(nscat,nemin:nemax,iv,ien,1:ntenum,ia))
		gm(iv,ien,ipart)=
     &	maxval(swk(nscat,nemin:nemax,iv,ien,1:ntenum,ipart))
	enddo
	enddo
	enddo
c
	where(gm.ne.0.0)pgm = 1/gm	!高速化のためgmの逆数pgmを定義
c
	deallocate(gm)
c
c	do ia = 1,narea
	do ipart = 1,npart	!07/8/4 不純物散乱	
	do ien = 1,nenergy
	do iv  = 1,nvalley
c		swk(1:nscat,1:nemax,iv,ien,1:ntenum,ia) =
c     &		min(1.0,swk(1:nscat,1:nemax,iv,ien,1:ntenum,ia)
c     &										* pgm(iv,ien,ia))
		swk(1:nscat,1:nemax,iv,ien,1:ntenum,ipart) =
     &		min(1.0,swk(1:nscat,1:nemax,iv,ien,1:ntenum,ipart)
     &										* pgm(iv,ien,ipart))	!07/8/4 不純物散乱	
	enddo
	enddo
	enddo
c
c
	return
	end
c
c
c===================================================================================================
	subroutine	getswk(
     &			  iv,ei,sei,qd21,bimp,
     &			  af,smh,ec,hwo,hw,
     &			  dos,poa,poe,aco,ope,opa,z,
     &			  swk,escat,iband,iarg,alloy,eth,a,b)
	include 'arraysize.fi'
	real(8) pi,q,h,bk,ep0,am0
	parameter(pi  = 3.141592654, q   = 1.60219e-19)
	parameter(h   = 1.05459e-34, bk  = 1.38066e-23)
	parameter(ep0 = 8.85419e-12, am0 = 9.10953e-31)
	real,dimension(nvalley):: smh,af,ec
	real	swk(nscat)
	real,	dimension (nscat)	:: escat
	integer(1),dimension (nscat)	:: iarg
	integer(1),dimension (nscat)	:: iband
	real(8)	ei,sei
	integer	iv
	real,dimension(nvalley):: z,dos,poe,poa,aco
	real,dimension(nvalley,nvalley):: ope,opa
c
	integer	jv
	real,	dimension (nvalley):: hwo
	real,	dimension (nvalley,nvalley):: hw
c
c	----------------------------------------
	real	qd21,bimp
	real(8)	zj,ef,sef,ak,qq,wk
	real(8)	qmax,qmin
	integer dj,iscat
c
c---(合金散乱)---
	real(8),dimension(nvalley) :: alloy
	real eth,a,b
	real alen
c
	iscat=1
c
c---( 有極性光学フォノン散乱 )---
cc	---放出---
	ef= ei-hwo(iv)
      if (ef.gt.0.0) then
		sef=sqrt(ef)
		qmax=sef+sei
		qmin=sei-sef
c		write(*,*) ei,ef,hwo(iv)
		swk(iscat)=poe(iv)*smh(iv)*sei/ei/q*log(qmax/qmin)
		escat(iscat)=-hwo(iv)
		iarg(iscat)=2
		iband(iscat)=iv
      endif
	iscat=iscat+1
cc	---吸収---
	ef=ei+hwo(iv)
	sef=sqrt(ef)
	qmax=sef+sei
	qmin=sef-sei
	if(ei.ne.0.0)then
		swk(iscat)=poa(iv)*smh(iv)*sei/ei/q*log(qmax/qmin)
	else
		swk(iscat)=0.0
	endif
	escat(iscat)=hwo(iv)
	iarg(iscat)=2
	iband(iscat)=iv
	iscat=iscat+1
c
c---( 音響フォノン散乱 )---
      ef=ei
      sef=sqrt(ef*(1.0+af(iv)*ef))
      swk(iscat)=aco(iv)*sef*dos(iv)*(1.0+2.0*af(iv)*ef)
	escat(iscat)=0
	iarg(iscat)=1
	iband(iscat)=iv
	iscat=iscat+1
c---( 不純物散乱 )---
      ef=ei
      sef=sqrt(ef*(1.0+af(iv)*ef))
      ak=smh(iv)*sef
      qq=qd21
	qq=qq*(4.0*ak*ak+qd21)
      wk=bimp/qq*sef*dos(iv)*(1.0+2.0*af(iv)*ef)
      swk(iscat)=wk
	escat(iscat)=0
	iarg(iscat)=3
	iband(iscat)=iv
	iscat=iscat+1
c
c---( バンド間フォノン散乱, from iv to jv)---
c	iv ... 遷移元バンド、	  jv ... 遷移先バンド
	dj = iscat
	do jv=1,nvalley
c		jev = ivbs+jv
c
		zj=z(jv)				!遷移先の谷数
		if(iv.eq.jv)zj=zj-1
cc	---放出 release---
		ef=ei-hw(iv,jv)+ec(iv)-ec(jv)
		if ((zj.gt.0.0).and.(ef.gt.0.0)) then
			sef=sqrt(ef*(1.0+af(jv)*ef))
			swk(iscat)=zj*ope(iv,jv)*sef*dos(jv)*(1.0+2.0*af(jv)*ef)
			continue
		else
			swk(iscat) = 0.0
		endif
		escat(iscat)=-hw(iv,jv)+ec(iv)-ec(jv)
		iarg(iscat)=1
		iband(iscat)=jv
		iscat=iscat+1
cc	---吸収 absorption---
		ef=ei+hw(iv,jv)+ec(iv)-ec(jv)
		if ((zj.gt.0.0).and.(ef.gt.0.0)) then
			sef=sqrt(ef*(1.0+af(jv)*ef))
			swk(iscat)=zj*opa(iv,jv)*sef*dos(jv)*(1.0+2.0*af(jv)*ef)
		else
			swk(iscat) = 0.0
		endif
		escat(iscat)=hw(iv,jv)+ec(iv)-ec(jv)
		iarg(iscat)=1
		iband(iscat)=jv
		iscat=iscat+1
	enddo
c
c---( 合金散乱 )---
      ef=ei
      sef=sqrt(ef*(1.0+af(iv)*ef))
      swk(iscat)=alloy(iv)*sef*(1.0+2.0*af(iv)*ef)
	escat(iscat)=0
	iarg(iscat)=1
	iband(iscat)=iv
	iscat=iscat+1
c
c---(衝突電離)---
	ef=ei
	alen=ei+ec(iv)
	if (alen.gt.eth) then
	  swk(iscat) = a*(((alen-eth)/eth)**b)
c	  swk(iscat) = a*(alen-eth)**b
	else 
	  swk(iscat)=0.0
	end if
	ef=0
	escat(iscat)=0
	iarg(iscat) =4
	iband(iscat)=1
	end
c
c===================================================================================================
	subroutine param_heat(hwo,hw,hescat)
c
c---	input	---
	include 'arraysize.fi'
	real, dimension (nvalley):: hwo
	real, dimension (nvalley,nvalley):: hw
c---	output	---
	real, dimension (nscat,nvalley):: hescat
c
	integer	iv,jv,iscat
c
c--- フォノンエネルギーのテーブルを作成 ---
	do iv = 1, nvalley
		hescat( 1,iv) =  hwo(iv)		!有極性光学フォノン散乱 放出
		hescat( 2,iv) = -hwo(iv)		!有極性光学フォノン散乱 吸収
		hescat( 3,iv) =  0				!音響フォノン散乱
		hescat( 4,iv) =  0				!不純物散乱
		do jv = 1, nvalley
			iscat = 5+(jv-1)*2
			hescat(iscat  ,iv) =  hw(jv,iv)		!バンド間フォノン散乱 放出
			hescat(iscat+1,iv) = -hw(jv,iv)		!バンド間フォノン散乱 吸収
		enddo
		hescat( 5+nvalley*2,iv) =  0	!合金散乱
		hescat( 5+nvalley*2+1,iv) =  0	!衝突電離
	enddo
c
	return
	end