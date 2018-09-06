c-----逆行列を求めるサブルーチン(完全ピボット法)-----
	subroutine gauss_jordan(ka,avsumconc,avsumtei1,ix,iz,iv,g,n,ind)
	implicit none
c
c---変数配列用パラメータ---
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c
	integer	i,j,k,m,n,p,q,nm,kw,nw,ind
	real g(20,22),l(20)
	real eps,d,aw,mw
c
	integer	ix,iz,iv
	integer(1)	ka
	real avsumconc(0:nx,0:nz,nvalley)
	real(8) avsumtei1(0:nx,0:nz,nvalley)
c
	ind=0
	m=1
	if(22.lt.n+m)then
		ind=10
		return
	endif
	nm=n+m
	eps=1e-7
	do j=1,n
		l(j)=j
	enddo
	d=1.0d0
	do 90 k=1,n
		aw=abs(g(k,k))
		p=k
		q=k
		do j=k,n
			do i=k,n
				if(aw.lt.abs(g(i,j)))then
				aw=abs(g(i,j))
					p=i
					q=j
				endif
			enddo
		enddo
		if(aw.lt.eps)then
			ind=11
			return
		elseif(k.ne.p)then
			d=-d
			do j=k,nm
				aw=g(k,j)
				g(k,j)=g(p,j)
				g(p,j)=aw
			enddo
		endif
		if(k.ne.q)then
			d=-d
			do i=1,n
				aw=g(i,k)
				g(i,k)=g(i,q)
				g(i,q)=aw
			enddo
			mw=l(k)
			l(k)=l(q)
			l(q)=mw
		endif
		d=d*g(k,k)
		kw=k+1
		do j=kw,nm
			g(k,j)=g(k,j)/g(k,k)
		enddo
		do 80 i=1,n
		if(i.ne.k)then
			do j=kw,nm
				g(i,j)=g(i,j)-g(i,k)*g(k,j)
			enddo
		endif
   80		continue
   90 continue
	nw=n+1
	do 120 j=nw,nm
		do i=1,n
			p=l(i)
			g(p,n)=g(i,j)
		enddo
		do i=1,n
			g(i,j)=g(i,n)
		enddo
  120 continue
	return
	end