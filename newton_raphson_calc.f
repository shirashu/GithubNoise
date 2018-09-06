	subroutine newton_raphson_calc(ka,avsumconc,avsumtei1,
     &					ix,iz,iv,n,x,ae,m,s,ind,flag3,am,aff)				!非放物線性

	implicit none
c
c---変数配列用パラメータ---
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c
	integer	i,j,k,m,n,ind,flag3		!非放物線性
	integer	ix,iz,iv
	integer(1)	ka
	real avsumconc(0:nx,0:nz,nvalley)
	real(8) avsumtei1(0:nx,0:nz,nvalley)
	real(8) aff(nvalley,narea)		!非放物線性
	real x(20),s(20),f(20),g(20,22),w(20)
	real e,h,s1,ae
	real, dimension (nvalley,narea)::am			!120201
c
	f=0.0
	g=0.0
	e=0.0
	s=0.0
	ind=0
c
	if(n.gt.20)then
		ind=33
		return
	endif
      h=1.0d0/500.0d0
	do 90 k=1,m
		call fermi_formula(ka,avsumconc,avsumtei1,ix,iz,iv,x,f,
     &             	       ind,aff,flag3,am)	!非放物線性
c-----stopしないための処理-----
		if(flag3.eq.1)then	!非放物線性
			return
		endif
c------------------------------
			do j=1,n
				w(j)=f(j)
			enddo
			do 30 i=1,n
				s1=x(i)
				x(i)=x(i)*(1.0d0+h)
			call fermi_formula(ka,avsumconc,avsumtei1,ix,iz,iv,x,f,
     &						   ind,aff,flag3,am)	!非放物線性
c-----stopしないための処理-----
			if(flag3.eq.1)then	!非放物線性
				return
			endif
c------------------------------
				do j=1,n
					g(j,i)=(f(j)-w(j)) / (h*s1)
				enddo
				x(i)=s1
   30	continue
			do j=1,n
				g(j,n+1)=-w(j)
			enddo
c
			call gauss_jordan(ka,avsumconc,avsumtei1,ix,iz,iv,g,n,ind)
			if(ind.ne.0) return
			do j=1,n
				s(j)=x(j)+g(j,n+1)
			enddo
			do j=1,n
				e=abs(g(j,n+1) / s(j))
				if(e.gt.ae) goto 70
			enddo
			m=k
			return
   70	do j=1,n
		x(j)=s(j)
	enddo
   90 continue
	ind=34
c
	return
	end