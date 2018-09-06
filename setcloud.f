	subroutine setcloud(dx,dz,cloud)
c	�d�׉_���K�E�X���z�Őݒ肷��T�u���[�`��
	implicit none
c
	include 'arraysize.fi'
	real dx,dz
	real cloud(-ngx:ngx,-ngz:ngz,ngn)	!�񎟌��d�׉_�i�[�z��
c
c  ---�����ϐ�---
	integer mesh		!���蓖�Ă郁�b�V����
	real	eps			!cloud(mesh,0)�̒��S�Ƃ̔�
	real	oval		!�ȉ~�W��
	real(8)	a,alpha		!a:�z��̍��v,alpha:���z�̍L����
	integer	ix,iz,jn	!ix,iz:���W,jn:�e�[�u��No.
	real dx2,dz2		!dx2,dz2:���z�̍L����w�W
	integer	nngx,nngz	!�K�E�X���z�e�[�u����
c
	oval = dx/dz/5	!/10		!�ȉ~�W��
	dx2 = dx*dx/oval
	dz2 = dz*dz*oval
	eps = 0.001			!cloud(mesh,0)�̒��S�Ƃ̔�
c
	cloud = 0.0			!cloud:�S�̔z��
	do jn=1,ngn
		mesh = 2**(jn-1)
		alpha = (dx*mesh)*(dx*mesh)/log(eps)
		a = 0.0				!=��(cloud)�ŏ�����
		nngz=min(mesh,ngz)
		nngx=min(mesh,ngx)
c
c
		do iz=-nngz,nngz
		do ix=-nngx,nngx
			cloud(ix,iz,jn)=exp(((dx2*ix*ix)+(dz2*iz*iz))/alpha)
			a=a+cloud(ix,iz,jn)
		enddo
		enddo
c
		do iz=-nngz,nngz
		do ix=-nngx,nngx
			cloud(ix,iz,jn)=cloud(ix,iz,jn)/a	!�W�����@��(cloud)=1.0��
		enddo
		enddo
c
	enddo
c
	end