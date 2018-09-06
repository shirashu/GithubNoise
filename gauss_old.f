	subroutine gauss(
     &				dt,spnum,jpnum,
     &				dx,dz,xmax,lhet,
     &				cxpart1,cxpart2,lxpart1,lxpart2,
     &				cxelrd1,cxelrd2,lnelrd1,lnelrd2,melpos,
     &				hm,hhm,af4,eps,p,kp,u,cur)
 	implicit none
c---変数配列用パラメータ---
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c
c---基本パラメータ---
	real	dt,spnum
	integer	jpnum
c---デバイス構造---
	real	dx,dz,xmax
	integer(2),dimension(narea)::lhet
	real,	dimension(npart):: cxpart1, cxpart2
	integer(2),dimension(npart)::lxpart1,lxpart2
c---電極---
	real	cxelrd1(npole),cxelrd2(npole)
	integer(2)	lnelrd1(npole),lnelrd2(npole)
	integer(1)	melpos(npole)
c---領域別パラメータ---
	real,	dimension (nvalley,narea)::	hhm,hm,af4
	real,	dimension (narea)::	eps
c---粒子状態---
	real	p(6,npmax)
	integer(1)	kp(3,npmax)
c---デバイス内状態---
	real	u(0:nx,0:nz)
	real	cur(npole)
c
	real q
	parameter(q	 = 1.60219e-19)
c
c---ローカル変数 ---
	real,allocatable,save	:: udt(:,:)
	real(8)	ss1,ss2,sv1,sv2
	real(8)	rs1,rs2,rv1,rv2
	real(8)	sk,gk
	real	v1,v2
	real(8)	sc1,sc2
	real	space
	real	right(2),left(2)
	integer(2)	ispace
	integer(2)	iright(2),ileft(2)
	integer	n,iz
	integer ka,kv
c-------------------
c
	if (.not. allocated(udt)) then
		allocate(udt(0:nx,0:nz))
		udt = u
	endif
	ss1 = 0.0
	ss2 = 0.0
	sv1 = 0.0
	sv2 = 0.0
	rs1 = 0.0
	rs2 = 0.0
	rv1 = 0.0
	rv2 = 0.0
c
	ispace	= 0 
	space	= float(ispace)*dx
	if(melpos(1).eq.1) left(1) = max(cxelrd2(1),cxpart1(1))
	if(melpos(1).eq.2) left(1) = max(space ,cxpart1(1))
	right(1) = min(cxelrd1(2),cxpart2(1))
	left(2) = max(cxelrd2(2),cxpart1(1))
	if(melpos(3).eq.1) right(2) = min(cxelrd1(3),cxpart2(1))
	if(melpos(3).eq.3) right(2) = min(xmax-space,cxpart2(1))
	do n=1,jpnum
c
c		ゲート電極の左側
		if((p(5,n).ge.left(1)).and.(p(5,n).le.right(1)))then
			sk = p(1,n)*p(1,n)+p(2,n)*p(2,n)+p(3,n)*p(3,n)
			kv = kp(1,n)
			ka = kp(3,n)
			gk = hhm(kv,ka)*sk
			v1 = p(1,n)*hm(kv,ka)/sqrt(1.0+af4(kv,ka)*gk)
			sv1 = sv1 + v1	!書き換え：意味同等
		endif
c
c		ゲート電極の右側
		if((p(5,n).ge.left(2)).and.(p(5,n).le.right(2)))then
			sk = p(1,n)*p(1,n)+p(2,n)*p(2,n)+p(3,n)*p(3,n)
			kv = kp(1,n)
			ka = kp(3,n)
			gk = hhm(kv,ka)*sk
			v2 = p(1,n)*hm(kv,ka)/sqrt(1.0+af4(kv,ka)*gk)
			sv2 = sv2 + v2	!書き換え：意味同等
		endif
c
	end do
c
	if(melpos(1).eq.1) ileft(1) = max(lnelrd2(1),lxpart1(1))
	if(melpos(1).eq.2) ileft(1) = max(ispace ,lxpart1(1))
	iright(1) = min(lnelrd1(2),lxpart2(1))
c
	ileft(2) = max(lnelrd2(2),lxpart1(1))
	if(melpos(3).eq.1) iright(2) = min(lnelrd1(3),lxpart2(1))
	if(melpos(3).eq.3) iright(2) = min(int2(nx)-ispace,lxpart2(1))
c
	do iz=0,nz
c
		ka = 1
		do while(ka<narea)
			if(iz.lt.lhet(ka))then
				exit
			endif
			ka=ka+1
		enddo
c
c		ゲート電極の左側
		if(ileft(1).lt.iright(1))then
			ss1 = ss1 + eps(ka)*u  (ileft (1),iz)
			ss1 = ss1 - eps(ka)*udt(ileft (1),iz)
			ss1 = ss1 - eps(ka)*u  (iright(1),iz)
			ss1 = ss1 + eps(ka)*udt(iright(1),iz)
		endif
c
c		ゲート電極の右側
		if(ileft(2).lt.iright(2))then
			ss2 = ss2 + eps(ka)*u  (ileft (2),iz)
			ss2 = ss2 - eps(ka)*udt(ileft (2),iz)
			ss2 = ss2 - eps(ka)*u  (iright(2),iz)
			ss2 = ss2 + eps(ka)*udt(iright(2),iz)
		endif
c
	end do
c
	if(right(1).ne.left(1))
     &		sc1 = ((-q*spnum*sv1)+(dz/dt*ss1))/(right(1)-left(1))
	if(right(2).ne.left(2))
     &		sc2 = ((-q*spnum*sv2)+(dz/dt*ss2))/(right(2)-left(2))
c
	cur(1)= sc1			!ソース電流
	cur(3)=-sc2			!ドレイン電流
	cur(2)=-sc1+sc2		!ゲート電流=－(ドレイン電流＋ソース電流)
c
	udt=u		!全体配列
	end subroutine gauss
