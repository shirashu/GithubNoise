c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c �ʎq�␳���v�Z����T�u���[�`��
c
c �����|�e���V�����̗̈敪��(�Q����)
c 2007/1/23, Hara, Ver.1.0
c ����F�\�ʂɋ߂Â����Ƃ��́A���̒l�����߂邽�߂̏���(EPVer3.0�ɑΉ�)
c
c	atype=0���ʏ�̗̈�Aatype=1���\�ʂ̊i�q�_
c	atype=2���\�ʂ���P���ꂽ�i�q�_
c	atype=3���\�ʂ���2���ꂽ�i�q�_
c
c	Ver1.0	�����[�X�o�[�W����
c
c	Input : lxr1,lzr1,lxr2,lzr2
c	Output: atype
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	subroutine area_type(lxr1,lzr1,lxr2,lzr2,atype)
	implicit none
c
c===����===
c---�ϐ��z��p�p�����[�^---
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c----(���Z�X)---
	integer(2),dimension (nrecess) :: lxr1,lzr1,lxr2,lzr2
c
	integer(1),	dimension (0:nx,0:nz)	:: atype
c	
c----(���[�J���ϐ�)----
	integer ix,iz,ii
c
	atype = 0
c
	do iz=0,nz
	  atype(2,iz)=3
	  atype(nx-2,iz)=3
	enddo
	do ix=0,nx
	  atype(ix,2)=3
	  atype(ix,nz-2)=3
	enddo
c	
	do ii=1,nrecess
	  do iz=lzr1(ii),lzr2(ii)
		atype(lxr1(ii)-2,iz)=3
		atype(lxr2(ii)+2,iz)=3
	  enddo
	  do ix=lxr1(ii),lxr2(ii)
		atype(ix,lzr2(ii)+2)=3
	  enddo
	enddo
c
	do iz=0,nz
	  atype(1,iz)=2
	  atype(nx-1,iz)=2
	enddo
	do ix=0,nx
	  atype(ix,1)=2
	  atype(ix,nz-1)=2
	enddo
c	
	do ii=1,nrecess
	  do iz=lzr1(ii),lzr2(ii)
		atype(lxr1(ii)-1,iz)=2
		atype(lxr2(ii)+1,iz)=2
	  enddo
	  do ix=lxr1(ii),lxr2(ii)
		atype(ix,lzr2(ii)+1)=2
	  enddo
	enddo
c
	do iz=0,nz
	  atype(0,iz)=1
	  atype(nx,iz)=1
	enddo
	do ix=0,nx
	  atype(ix,0)=1
	  atype(ix,nz)=1
	enddo
c
	do ii=1,nrecess
	  do iz=lzr1(ii),lzr2(ii)
		atype(lxr1(ii),iz)=1
		atype(lxr2(ii),iz)=1
	  enddo
	  do ix=lxr1(ii),lxr2(ii)
		atype(ix,lzr2(ii))=1
	  enddo
	enddo
c
	return
	end
