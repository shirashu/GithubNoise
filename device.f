	subroutine device(np1,dx,dz,
     &				cxpart1,czpart1,cxpart2,czpart2,
     &                lxpart1,lzpart1,lxpart2,lzpart2,
     &				cxrecess1,czrecess1,cxrecess2,czrecess2,
     &                lxrecess1,lzrecess1,lxrecess2,lzrecess2,
     &				dconc,xmax,zmax,
     &				lhet,iarea,twodeg,
     &				cxpole1,cxpole2,lnpole1,lnpole2,vb,
     &				melpos,lg,spnum,ncon)
c
      implicit none
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c
	integer	np1
	real dx,dz
	real,		dimension (npart)	:: cxpart1,czpart1,cxpart2,czpart2
	integer(2),	dimension (npart)	:: lxpart1,lzpart1,lxpart2,lzpart2
c
c	----リセス---
	real,		dimension (nrecess)	:: cxrecess1,czrecess1
	real,		dimension (nrecess)	:: cxrecess2,czrecess2
	integer(2),	dimension (nrecess)	:: lxrecess1,lzrecess1
	integer(2),	dimension (nrecess)	:: lxrecess2,lzrecess2
c
	real,		dimension (npart)	:: dconc
	real xmax,zmax
c
	integer(2),	dimension (nlayer)	:: lhet
	real,		dimension (nlayer)	:: twodeg
	integer(1),	dimension (nlayer)	:: iarea

	integer(1),	dimension (npole)	:: melpos
	real,		dimension (npole)	:: cxpole1,cxpole2
	integer(2),	dimension (npole)	:: lnpole1,lnpole2
	real,		dimension (npole)	:: vb
	real lg
c
	real spnum
	integer ncon
c
c	--- 内部変数 ---
	integer i
	real	mmax,dd,dmax,dmax2
	real,	dimension (npole)	:: dpole,daddr
	real,	dimension (nlayer)	:: dhet
c	character(80) form
c
c---( デバイス構造 )---
c
	read(8,*) dx		!1.0e-8
	read(8,*) dz		!2.5e-9
	read(8,*) (cxpart1(i),	i = 1, npart)
	read(8,*) (cxpart2(i),	i = 1, npart)
	read(8,*) (czpart1(i),	i = 1, npart)
	read(8,*) (czpart2(i),	i = 1, npart)
	read(8,*) (cxrecess1(i),	i = 1, nrecess)
	read(8,*) (cxrecess2(i),	i = 1, nrecess)
	read(8,*) (czrecess1(i),	i = 1, nrecess)
	read(8,*) (czrecess2(i),	i = 1, nrecess)
	read(8,*) (dconc(i),	i = 1, npart)	!m^-3 (*1.0e-6 cm^-3)	不純物濃度
	read(8,*) (dhet(i), 	i = 1, nlayer)
	read(8,*) (iarea(i),	i = 1, nlayer)	!iarea(iar):1->InAlAs, 2->InGaAsを表す
	read(8,*) (twodeg(i),	i = 1, nlayer)	!2DEGシート電荷量(*1.0e-4 cm^-3)
	read(8,*) (dpole(i), 	i = 1, npole)	!各電極の長さ
	read(8,*) (melpos(i),	i = 1, npole)	!各電極の付く位置(1:上、2:左、3:右)
	read(8,*) (daddr(i),	i = 1, npole)	!各電極間の距離	dpole+lair=xmax 
	read(8,*) (vb(i),		i = 1, npole)  	!各電極のビルドイン ポテンシャル
c
c	---(デバイスサイズ)---
	xmax = maxval(cxpart2)
	zmax = maxval(czpart2)
c
c	---(メッシュサイズ)---
	nx=nint(xmax/dx)
	nz=nint(zmax/dz)
c
	do i=1,nlayer
		dhet(i)=max(0.0,min(zmax,dhet(i)))
		lhet(i)=max(0,min(nz,nint(dhet(i)/dz)))
c	write(*,*) i, lhet(i)
	enddo
c
c	---( 領域別メッシュ設定 )---
	do i=1,npart
	cxpart1(i)=max(0.0,min(xmax,cxpart1(i)))
	cxpart2(i)=max(0.0,min(xmax,cxpart2(i)))
	czpart1(i)=max(0.0,min(zmax,czpart1(i)))
	czpart2(i)=max(0.0,min(zmax,czpart2(i)))
	enddo
c
c	---(リセス)---
	do i=1,nrecess
	cxrecess1(i)=max(0.0,min(xmax,cxrecess1(i)))
	cxrecess2(i)=max(0.0,min(xmax,cxrecess2(i)))
	czrecess1(i)=max(0.0,min(zmax,czrecess1(i)))
	czrecess2(i)=max(0.0,min(zmax,czrecess2(i)))
	enddo
c
	lxpart1(1:npart) = nint(cxpart1(1:npart)/dx)
	lzpart1(1:npart) = nint(czpart1(1:npart)/dz)
	lxpart2(1:npart) = nint(cxpart2(1:npart)/dx)
	lzpart2(1:npart) = nint(czpart2(1:npart)/dz)
c
c	---(リセス)---
	lxrecess1(1:nrecess) = nint(cxrecess1(1:nrecess)/dx)
	lzrecess1(1:nrecess) = nint(czrecess1(1:nrecess)/dz)
	lxrecess2(1:nrecess) = nint(cxrecess2(1:nrecess)/dx)
	lzrecess2(1:nrecess) = nint(czrecess2(1:nrecess)/dz)
c	write(*,*) lxrecess1(1),lxrecess2(1),lzrecess1(1),lzrecess2(1)
c	write(*,*) lxrecess1(2),lxrecess2(2),lzrecess1(2),lzrecess2(2)
c
c	---(電極)---
	lg   = dpole(2)
	do i=1,npole
		cxpole1(i) = daddr(i)
		cxpole2(i) = daddr(i)+dpole(i)
		if(melpos(i).eq.1)then
			mmax = xmax
			dd  = dx
		else
			mmax = zmax
			dd  = dz
		endif
		lnpole1(i)	= nint(cxpole1(i)/dd)
c		cxpole1(i)	= max(0.0,min(mmax,dd*(float(lnpole1(i)))))
		lnpole2(i)	= nint(cxpole2(i)/dd)
c		cxpole2(i)	= max(0.0,min(mmax,dd*(float(lnpole2(i)))))
c		write(*,*) lnpole1(i),lnpole2(i)
	enddo

c	dmax = maxval(dconc)
c	dmax = dconc(4)		!dconc(4)?
c	dmax = dconc(3)		!11/03/19 原
	dmax = dconc(ndope)		!ドープ層に合わせる 11/03/19原
	dmax2 = dconc(npart-1)			!06/12/28
	spnum = max(1.0,dmax*dx*dz/float(np1))
c
	open(unit=570,file='spnum.txt')		!! 11/04/20原
	write(*,*) spnum, '<-spnum'
	write(570,*) spnum, 'spnum'	!!　11/04/20原
c	
c	ncon = nint(dmax*dx*dz/spnum/2.0)
	ncon = nint(dmax2*dx*dz/spnum/2.0)	!06/12/28
c
	close(570)		!!　11/04/20原
c
	return
	end