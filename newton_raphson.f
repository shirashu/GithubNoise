c-----Newton_Raphson法により収束解を求めるサブルーチン-----
	subroutine newton_raphson(ka,ini_fermi,ini_tel,avsumconc,avsumtel,
     &				avsumtei1,efermi,ix,iz,iv,ind,flag3,am,aff)		!非放物線性


	implicit none
c
c---変数配列用パラメータ---
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c
	integer	n,m,i,ind,flag3		!非放物線性
	integer	ix,iz,iv
	integer(1)	ka
	real x(20),s(20),ae
	real ini_fermi,ini_tel
	real efermi(0:nx,0:nz,nvalley)
	real avsumtel(0:nx,0:nz,nvalley)
	real avsumconc(0:nx,0:nz,nvalley)
	real(8) avsumtei1(0:nx,0:nz,nvalley)
	real, dimension (nvalley,narea)::am			!120201
	real(8) aff(nvalley,narea)					!120201
c
	real(8) bk,q,bkq
	parameter(bk = 1.38066e-23,q = 1.60219e-19)
c
	n=2
	x(1)=ini_tel		!Telの初期値
	x(2)=ini_fermi		!Efの初期値
	ae=0.001			!相対許容誤差 非放物線性を入れた後は0.001に変更した
	m=100000			!最大繰り返し回数
c
c-----x(2)をitaに変換-----
	bkq=bk/q
	x(2)=x(2)/(bkq*x(1))
c-------------------------
c
	call newton_raphson_calc(ka,avsumconc,avsumtei1,
     &						 ix,iz,iv,n,x,ae,m,s,ind,flag3,am,aff)		!非放物線性
c
c-----stopしないための処理-----
	if(flag3.eq.1)then	!非放物線性
		return
	endif
c------------------------------
c
c-----itaからEfを求める-----
	bkq=bk/q
	s(2)=bkq*s(1)*s(2)
c---------------------------
c
	avsumtel(ix,iz,iv) = s(1)
	efermi(ix,iz,iv)   = s(2)
c
	return
	end