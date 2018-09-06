	subroutine	charge(jpnum,dx,dz,spnum,p,kp,cn,icn,jpot,ibord,cloud,
     &					lxr1,lzr1,lxr2,lzr2)
c
c===( �d�ו��z�̌v�Z )====
c
c === �ϐ���� ===
c	--- ���� ---
c	p(1-6,n) ... ���q���(1-3:k���W(kx,ky,kz)[m^-1],4:�U������[s],5-6:�ʒu(x,y))[m]
c	iv(n) ... ���q�̏������Ă���J
c	imax ... �ő�L��������
c	jpnum ... �f�o�C�X���ɂ��钴���q�̐�
c	spnum ... �����q�̗��q��
c	cn(i) ... �i�q�_(ix,iz)�ɂ�����d�q�Z�x[m^-3]
c	dx,dz ... ���b�V���T�C�Y[m]
c	nx,nz ... �i�q�_�̕����̐�
	implicit none
c	
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c
	integer jpnum
	real dx,dz,spnum
	integer(2),	dimension (nrecess)	:: lxr1,lzr1,lxr2,lzr2
	real p(6,npmax)
	integer(1) kp(3,npmax)
	real cn(0:nx,0:nz)
	integer	icn,jpot
	integer	ibord
	real cloud(-ngx:ngx,-ngz:ngz,ngn)
c
	real x,z,x2,z2,xb,zb,dxyz
	real cn2,cn3
	allocatable	::	cn2(:,:),cn3(:,:)
	integer ii
	integer ix,iz,n,idx,idz
	integer iix,iiz,jn
	integer jxmin,jxmax,jzmin,jzmax
	integer	mesh,nngx,nngz,np1
	real	dnp2,sq
c
c---( �Z���d�׉_�@�ŗ��q���z�����b�V���ɐU�蕪�� )----
      do n=1,jpnum
c		if (kp(1,n).eq.0) cycle
		x=p(5,n)/dx
		z=p(6,n)/dz
c		ix=max(0,min(ifix(x),nx-1))
		ix = ifix(x)
		if(ix.lt.0)then
			ix = 0
		elseif(ix.gt.nx-1)then
			ix = nx-1
		endif
c		iz=max(0,min(ifix(z),nz-1))
		iz = ifix(z)
		if(iz.lt.0)then
			iz = 0
		elseif(iz.gt.nz-1)then
			iz = nz-1
		endif

		xb=float(ix)
		zb=float(iz)
		x2=1.0-(x-xb)
		z2=1.0-(z-zb)
		cn(ix  ,iz  )=cn(ix  ,iz  )+(    x2)*(    z2)
		cn(ix  ,iz+1)=cn(ix  ,iz+1)+(    x2)*(1.0-z2)
		cn(ix+1,iz  )=cn(ix+1,iz  )+(1.0-x2)*(    z2)
		cn(ix+1,iz+1)=cn(ix+1,iz+1)+(1.0-x2)*(1.0-z2)
	enddo
c
	icn = icn + 1
c
	if(icn.eq.jpot)then		!�o�͒i�K
c
		goto 10		!�K�E�X���z�d�׉_�@�͎g��Ȃ�
c
		jxmin =	0	!ixd1(1)	!+5
		jxmax = nx	!ixd2(1)	!-5
		jzmin =	ibord;jzmax = nz
c
		allocate(cn2(jxmin:jxmax,jzmin:jzmax))
		cn2 = 0.0
		dnp2 = (1.0e+21*dx*dz)/spnum*icn /2.0	!���q������郁�b�V���L��
		sq = sqrt(max(0.0,dnp2))
		mesh = nint(((1.0/sq)-1.0)/2.0)
		nngz=min(ngz,mesh)
		nngx=min(ngx,mesh)
		mesh=(2*nngx+1)*(2*nngz+1)
		do ix=jxmin,jxmax
		do iz=jzmin,jzmax
			if(cn(ix,iz).eq.0.0)		cycle
			do idz = -nngz,+nngz
				iiz = iz+idz
				if(iiz.lt.jzmin)iiz=jzmin+(jzmin-iiz)
				if(iiz.gt.jzmax)iiz=jzmax-(iiz-jzmax)
			do idx = -nngx,+nngx
				iix = ix+idx
				if(iix.lt.jxmin)iix=jxmin+(jxmin-iix)
				if(iix.gt.jxmax)iix=jxmax-(iix-jxmax)
				cn2(iix,iiz)=cn2(iix,iiz)+cn(ix,iz)/mesh
			enddo
			enddo
		enddo
		enddo
c
c---( �Q�����K�E�X���z�ŕ��ω� )----
c	�Z�x�ɂ���čL���镝�ύX
		allocate(cn3(jxmin:jxmax,jzmin:jzmax))
		cn3 = 0.0
		np1 = nint((1.0e+24*dx*dz)/spnum)
		do ix=jxmin,jxmax
		do iz=jzmin,jzmax
c
			if(cn(ix,iz).eq.0.0)		cycle
c
			do jn = 1,ngn-1
				if(cloud(0,0,jn)*np1.lt.(cn2(ix,iz)))exit
			enddo
c
			mesh = 2**(jn-1)
			nngz=min(mesh,ngz)
			nngx=min(mesh,ngx)
			do idz = -nngz,+nngz
				iiz = iz+idz
				if(iiz.lt.jzmin)iiz=jzmin+(jzmin-iiz)
				if(iiz.gt.jzmax)iiz=jzmax-(iiz-jzmax)
			do idx = -nngx,+nngx
				iix = ix+idx
				if(iix.lt.jxmin)iix=jxmin+(jxmin-iix)
				if(iix.gt.jxmax)iix=jxmax-(iix-jxmax)
				cn3(iix,iiz)=cn3(iix,iiz)+cloud(idx,idz,jn)*cn(ix,iz)
			enddo
			enddo
			cn(ix,iz)=0.0
		enddo
		enddo
		cn(jxmin:jxmax,jzmin:jzmax) = !cn(jxmin:jxmax,jzmin:jzmax) +
     &								  cn3(jxmin:jxmax,jzmin:jzmax)
c---------------------------
c
		deallocate(cn2,cn3)
10	continue
c
c---( �L�����A��->�L�����A�Z�x�̕ϊ� )----
		dxyz	= spnum/dx/dz/float(icn)
		cn	= cn * dxyz		!cn:�S�̔z��
		cn(0 ,0:nz)=cn(0 ,0:nz)*2.0		!���[
		cn(nx,0:nz)=cn(nx,0:nz)*2.0		!�E�[
c		cn(0:nx,0 )=cn(0:nx,0 )*2.0		!��[
		cn(0:nx,nz)=cn(0:nx,nz)*2.0		!���[
c	----(���Z�X)----
		cn(0:minval(lxr1),0) = cn(0:minval(lxr1),0)*2.0		!����[
		cn(maxval(lxr2):nx,0) = cn(maxval(lxr2):nx,0)*2.0	!�E��[

		do ii=1,nrecess
c			cn(lxr1(ii):lxr2(ii),lzr1(ii))=cn()*2
			cn(lxr1(ii):lxr2(ii),lzr2(ii))=cn(lxr1(ii):lxr2(ii),lzr2(ii))*2.0 !����
			cn(lxr1(ii),lzr1(ii):lzr2(ii))=cn(lxr1(ii),lzr1(ii):lzr2(ii))*2.0 !���[
			cn(lxr2(ii),lzr1(ii):lzr2(ii))=cn(lxr2(ii),lzr1(ii):lzr2(ii))*2.0 !�E�[
			cn(lxr1(ii),lzr2(ii)) = cn(lxr1(ii),lzr2(ii))/3.0				  !���p
			cn(lxr2(ii),lzr2(ii)) = cn(lxr2(ii),lzr2(ii))/3.0				  !�E�p
		enddo	
	endif
	end