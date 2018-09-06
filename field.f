c
c---- ( �d�E�v�Z�p���v���O���� ) ----
c
	subroutine field(x,z,pdx,pdz,u,fx,fz,dhet,kl,
     &				cxr1,czr1,cxr2,czr2,lxr1,lzr1,lxr2,lzr2)
c	
	implicit none
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx, nz,ntenum
	real,		dimension (nrecess)	:: cxr1,czr1,cxr2,czr2
	integer(2),	dimension (nrecess)	:: lxr1,lzr1,lxr2,lzr2
	real	x, z, pdx, pdz, u(0:nx,0:nz)
	real	fx, fz
	real	dhet(0:nlayer)
	integer(1)	kl
c
	integer ix, iz, ii
	real	xp, zp, xn, zn
	real,	dimension (nrecess)	:: hxd3
c
c
c	if (kv.eq.0) return
c
c---- ( �����i�q�_�̌��� ) ----
c
	xp = x*pdx
	zp = z*pdz
c
c	ix = max(0,min(ifix(xp),nx-1)) !���̃A���S���Y���̕�������(��7.2%)
c	iz = max(0,min(ifix(zp),nz-1)) !���̃A���S���Y���̕�������(��7.2%)
c
c	ix = ifix(xp)
c	if(ix.lt.0)then
c		ix = 0
c	elseif(ix.gt.nx-1)then
c		ix = nx-1
c	endif
c	iz = ifix(zp)
c	if(iz.lt.0)then
c		iz = 0
c	elseif(iz.gt.nz-1)then
c		iz = nx-1
c	elseif(z.eq.dhet(kl-1))then
c		iz = ifix(zp+0.5)	!��Ɋۂ�
c	elseif(z.eq.dhet(kl))then
c		iz = ifix(zp-0.5)	!���Ɋۂ�
c	endif
c
c	----(���Z�X)----
	ix = ifix(xp)
	iz = ifix(zp)
c
	if(ix.lt.0)then						
		ix = 0								!���O�̈�
	elseif(ix.gt.nx-1)then				
		ix = nx-1							!�E�O�̈�
	endif
c	
	if(iz.lt.0)then							!��O�̈�
		iz = 0
	elseif(iz.gt.nz-1)then					!���O�̈�
		iz = nz-1
	endif
c
	if(z.eq.dhet(kl-1))then		!���E��
		iz = ifix(zp+0.5)	!��Ɋۂ�
	elseif(z.eq.dhet(kl))then
		iz = ifix(zp-0.5)	!���Ɋۂ�
	endif
c
	hxd3 = lxr1+(lxr2-lxr1)/2.0
	do ii = 1, nrecess
	  if((iz.ge.lzr1(ii)).and.(iz.lt.lzr2(ii))) then				
	    if(ix.le.hxd3(ii)) then
	      ix = min(ix,lxr1(ii))		!���Z�X���̈�
	    else
	      ix = max(lxr2(ii),ix)		!���Z�X�E�̈�
	    endif
	  endif
	enddo
c
c-----------------------------------------------------
c	ix = max(0,min(ifix(xp),nx-1))
c	iz = max(0,min(ifix(zp),nz-1))
c
	xn = xp-float(ix)
	zn = zp-float(iz)
c
c---- ( �_(x,z)�̓d�E���x(fx,fz) ) ----
	fx = ((u(ix,iz)-u(ix+1,iz))*(1-zn)
     &	      +(u(ix,iz+1)-u(ix+1,iz+1))*zn)*pdx
	fz = ((u(ix,iz)-u(ix,iz+1))*(1-xn)
     &		  +(u(ix+1,iz)-u(ix+1,iz+1))*xn)*pdz

	return
	end