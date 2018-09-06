c-----仮のEfを決定するサブルーチン(Γ谷のみ)-----
	subroutine fermi_calc(efermi,ka3,
     &					  tel1,tel2,tel3,efermi1,efermi2,efermi3,
     &					  avsumconc,avsumtel,avsumtei1,ntab1,etab1,
     &                    avsumteiA,avsumconcA,lhet,am,aff)		!竹岸変更

	implicit none
c	===引数===
c---変数配列用パラメータ---
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c---基本パラメータ---
	integer	ix,iz,iv,ia,i1,ind,flag3 !非放物線性
	integer(1)	ka,ka3(0:nx,0:nz,nvalley)
	real, dimension (nvalley,narea)::am	
	real(8) aff(nvalley,narea)					!120201
c---濃度・Ef・Tel---
	real efermi(0:nx,0:nz,nvalley)
	real,dimension(:),save,allocatable	:: ntabb
c
	real,	dimension(0:nx,0:nz)	:: tel1,efermi1
	real,	dimension(0:nx,0:nz)	:: tel2,efermi2
	real,	dimension(0:nx,0:nz)	:: tel3,efermi3
c
c-----08/8/6 竹岸追加-----
	real avsumconc(0:nx,0:nz,nvalley)
	real avsumconcA(0:nx)          !101221電子濃度チャネル平均2
	real avsumtel(0:nx,0:nz,nvalley)
	real(8) avsumtei1(0:nx,0:nz,nvalley)
	real(8) avsumteiA(0:nx)                 !101221電子エネルギーチャネル平均2
	real ini_fermi,ini_tel		!Fermi_Level,Telの初期値
c
	real ntab1(0:300000,5)		!濃度(番号、材料)
	real etab1(0:300000,5)		!エネルギー(番号、材料)
	integer(2)	lhet(nlayer)
c-------------------------
c
	if(.not. allocated(ntabb))then
		allocate(ntabb(0:5))
		efermi  = -1.0595		!粒子のいないところはka=4のX谷Eg/2
		efermi1 = -1.0595		!粒子のいないところはka=4のX谷Eg/2
		efermi2 = -1.0595		!粒子のいないところはka=4のX谷Eg/2
		efermi3 = -1.0595		!粒子のいないところはka=4のX谷Eg/2
	endif
c
c-----Γ谷のみ考慮しているためEc=0と考え、状態密度関数の(E-Ec)**(1/2)のEc=0
c	全ての谷を考慮する場合、Ec=0以外も考える必要あり
	iv=1				!Γ谷のみ考慮

c-----------101221電子エネルギー，濃度チャネル平均------------
c-----プログラムに含める：シート平均　プログラムから外す：メッシュ平均 !120126homma
c	do iz = lhet(nchannel1)+1,lhet(nchannel2)-1
c	     do ix=0,nx
c	          avsumconc(ix,iz,1) =avsumconcA(ix)
c	          avsumtei1(ix,iz,1) =avsumteiA(ix)
c		 enddo
c      enddo
c-------------------------------------------------------------

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	do iz=38,46			!Channel層のEfとTel算出 InAs
	do iz = lhet(nchannel1)+1,lhet(nchannel2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		do ix=1,nx-1
c
	ka = ka3(ix,iz,iv)		!素材No.→1～4
c
c-----粒子のいないところのフェルミレベルと電子温度-----
      if(ka.eq.0)then
	goto 1201
      endif
c
	if(avsumconc(ix,iz,iv).eq.0.0)then
		efermi(ix,iz,iv) = -1.0595  !ka=4のX谷のEgの1/2の値(=minEg)
		avsumtel(ix,iz,iv)=0.0		!粒子のいないメッシュは電子温度0[K]
		goto 1201
	endif
c------------------------------------------------------
c

c
c	write(*,*) 'ka',ka
	if(ka.eq.1)then			!InAs	11/04/21原
		goto 12
	elseif(ka.eq.2)then		!In(0.52)Al(0.48)As	 11/04/21原
		goto 13
	elseif(ka.eq.3)then		!!InAs	　11/04/21原
		goto 14
	elseif(ka.eq.4)then		!!InAs
		goto 15
	elseif(ka.eq.5)then     !!InAs
		goto 16
	else
	continue
		write(99,*)'適する材料がみつかりません'
		stop
	endif
ccccccccccccccccccccccccccccccccccccccccccccccccc
c	以降設定方法が分からない!!2011/03/22　原
cccccccccccccccccccccccccccccccccccccccccccccccc
c
c-----InAs (old_In(0.52)Al(0.48)As)-----11/04/21原
c   12 continue
c		write(99,*)'InAlAsに入りました'
c		stop
   12 continue    
	if((avsumconc(ix,iz,iv).gt.0.0).and.
     &(avsumconc(ix,iz,iv).le.1e18))then
		avsumconc(ix,iz,iv) = 1e18
	elseif(avsumconc(ix,iz,iv).ge.1e25)then
		avsumconc(ix,iz,iv) = 1e25
	endif
c
	do i1=0,100000,10
		if(avsumconc(ix,iz,iv).le.ntab1(i1,1))then
			avsumconc(ix,iz,iv) = ntab1(i1,1)
			exit
		endif
	enddo
c
	if(avsumtei1(ix,iz,iv).le.etab1(i1,1))then
		avsumtei1(ix,iz,iv) = etab1(i1,1)
	endif
c
c-----Efの初期値の設定-----
	ini_fermi = -0.30
	ini_tel   = 300
	if((avsumconc(ix,iz,iv).ge.5e24)
     &.and.(avsumconc(ix,iz,iv).le.1e25))then
		ini_fermi = 0.20
	elseif((avsumconc(ix,iz,iv).ge.1e24)
     &.and.(avsumconc(ix,iz,iv).lt.5e24))then
		ini_fermi = 0.10
	elseif((avsumconc(ix,iz,iv).ge.5e23)
     &.and.(avsumconc(ix,iz,iv).lt.1e24))then
		ini_fermi = 0.03
	endif
c
	goto 1230
c
c-----Al(0.52)In(0.48)As (old_In(0.53)Ga(0.47)As)----- 11/04/21原
   13 continue
	if((avsumconc(ix,iz,iv).gt.0.0).and.
     &(avsumconc(ix,iz,iv).le.1e18))then
		avsumconc(ix,iz,iv) = 1e18
	elseif(avsumconc(ix,iz,iv).ge.1e25)then
		avsumconc(ix,iz,iv) = 1e25
	endif
c
	do i1=0,100000,10
		if(avsumconc(ix,iz,iv).le.ntab1(i1,2))then
			avsumconc(ix,iz,iv) = ntab1(i1,2)
			exit
		endif
	enddo
c
	if(avsumtei1(ix,iz,iv).le.etab1(i1,2))then
		avsumtei1(ix,iz,iv) = etab1(i1,2)
	endif
c
c-----Efの初期値の設定-----
	ini_fermi = -0.30  ! -0.05
	ini_tel   = 300
	if((avsumconc(ix,iz,iv).ge.5e24)
     &.and.(avsumconc(ix,iz,iv).le.1e25))then
		ini_fermi = 0.20
	elseif((avsumconc(ix,iz,iv).ge.1e24)
     &.and.(avsumconc(ix,iz,iv).lt.5e24))then
		ini_fermi = 0.10
	elseif((avsumconc(ix,iz,iv).ge.5e23)
     &.and.(avsumconc(ix,iz,iv).lt.1e24))then
		ini_fermi = 0.03
	endif
c
	goto 1230
c
c-----InAs (old_InAs①)-----
   14 continue
	if((avsumconc(ix,iz,iv).gt.0.0).and.
     &(avsumconc(ix,iz,iv).le.1e18))then
		avsumconc(ix,iz,iv) = 1e18
	elseif(avsumconc(ix,iz,iv).ge.1e25)then
		avsumconc(ix,iz,iv) = 1e25
	endif
c
c-----平均濃度と平均エネルギーの補正-----
	do i1=0,100000,10
		if(avsumconc(ix,iz,iv).le.ntab1(i1,3))then
			avsumconc(ix,iz,iv) = ntab1(i1,3)
			exit
		endif
	enddo
c
	if(avsumtei1(ix,iz,iv).le.etab1(i1,3))then
		avsumtei1(ix,iz,iv) = etab1(i1,3)
	endif
c
c-----Efの初期値の設定-----
	ini_fermi = -0.30
	ini_tel   = 300
	if((avsumconc(ix,iz,iv).ge.5e24)
     &.and.(avsumconc(ix,iz,iv).le.1e25))then
		ini_fermi = 0.30
	elseif((avsumconc(ix,iz,iv).ge.1e24)
     &.and.(avsumconc(ix,iz,iv).lt.5e24))then
		ini_fermi = 0.20
	elseif((avsumconc(ix,iz,iv).ge.5e23)
     &.and.(avsumconc(ix,iz,iv).lt.1e24))then
		ini_fermi = 0.05
	endif
c
	goto 1230
c
c-----InP-----
   15 continue
		write(99,*)'InPに入りました'
		stop
c
c
c-----channel材料-----
   16 continue
	if((avsumconc(ix,iz,iv).gt.0.0).and.
     &(avsumconc(ix,iz,iv).le.1e18))then
		avsumconc(ix,iz,iv) = 1e18
	elseif(avsumconc(ix,iz,iv).ge.1e25)then
		avsumconc(ix,iz,iv) = 1e25
	endif
c
c-----平均濃度と平均エネルギーの補正-----
	do i1=0,100000,10
		if(avsumconc(ix,iz,iv).le.ntab1(i1,5))then
			avsumconc(ix,iz,iv) = ntab1(i1,5)
			exit
		endif
	enddo
c
	if(avsumtei1(ix,iz,iv).le.etab1(i1,5))then
		avsumtei1(ix,iz,iv) = etab1(i1,5)
	endif
c
c-----Efの初期値の設定-----
	ini_fermi = -0.30
	ini_tel   = 300
	if((avsumconc(ix,iz,iv).ge.5e24)
     &.and.(avsumconc(ix,iz,iv).le.1e25))then
		ini_fermi = 0.20
	elseif((avsumconc(ix,iz,iv).ge.1e24)
     &.and.(avsumconc(ix,iz,iv).lt.5e24))then
		ini_fermi = 0.10
	elseif((avsumconc(ix,iz,iv).ge.5e23)
     &.and.(avsumconc(ix,iz,iv).lt.1e24))then
		ini_fermi = 0.05
	endif
c
	goto 1230



c-----平均化された濃度とエネルギーからEfとTelを求める 竹岸-----
 1230 continue
	call newton_raphson(ka,ini_fermi,ini_tel,avsumconc,avsumtel,
     &                    avsumtei1,efermi,ix,iz,iv,ind,flag3,am,aff) !非放物線性 変更
c
 1201 continue
c
c-----収束されなかったとき-----
	flag3 = 0	!非放物線性のフラグ
	if(ind.ne.0)then
		ini_fermi = ini_fermi + 0.001   !非放物線性 変更
		if((ini_fermi.lt.0.85).and.(ini_tel.lt.2500))then
			goto 1230
		endif
		if((ini_fermi.gt.0.85).and.(ini_tel.lt.2500))then
			ini_fermi = -0.30
			ini_tel   = ini_tel + 100
			goto 1230
		endif
	endif
c
c-----NaN計算されたとき-----
	if(isnan(efermi(ix,iz,iv)))then
		ini_fermi = ini_fermi + 0.001   !非放物線性 変更
		if((ini_fermi.lt.0.85).and.(ini_tel.lt.2500))then
			goto 1230
		endif
		if((ini_fermi.gt.0.85).and.(ini_tel.lt.2500))then
			ini_fermi = -0.30
			ini_tel   = ini_tel + 100
			goto 1230
		endif
	endif
c
c-----エラー出力-----
	if(ind.ne.0)then
		write(3,*)'N=',avsumconc(ix,iz,iv),'E=',avsumtei1(ix,iz,iv)
		write(3,*)'ix=',ix,'iz=',iz,'ka=',ka,'ind=',ind
		write(3,*)'ini_fermi=',ini_fermi
		write(3,*)'Tel=',avsumtel(ix,iz,iv),'Ef=',efermi(ix,iz,iv)
	endif
c
	if(isnan(efermi(ix,iz,iv)))then
		write(3,*)'N=',avsumconc(ix,iz,iv),'E=',avsumtei1(ix,iz,iv)
		write(3,*)'ix=',ix,'iz=',iz,'ka=',ka,'ind=',ind
		write(3,*)'ini_fermi=',ini_fermi
		write(3,*)'Tel=',avsumtel(ix,iz,iv),'Ef=',efermi(ix,iz,iv)
	endif
c--------------------
		enddo
	enddo
c
	return
	end