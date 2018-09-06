	subroutine temperature(mxpart1,mxpart2,ntemp,btmp,dtmp)
c---	�f�o�C�X�����x���z���t�@�C������ǂݍ��݃e�[�u���ɂ���T�u���[�`��
c
	implicit none
c === �ϐ���� ===
c ---	input ---
c	ntenum ... ���x���z������
c	mxpart1,mxpart2 ... 	�e�G���A�̍��W[mesh]
c ---	output ---
c	btmp ... �x�[�X���x(�f�o�C�X���̍Œቷ�x)
c	dtmp ... ntemp1������̉��x�㏸����T
c	ntemp(ix,iz) ... �e���b�V����btmp����̉��x�㏸���e�[�u��
c
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c
	integer(2) mxpart1,mxpart2
	integer(2) ntemp(0:nx,0:nz)		!���x���z�e�[�u��
	real	btmp,dtmp
c
c###	���[�J���ϐ�  ###
	integer	ix,iz
	real	tmin,tmax
	real,allocatable ::temp(:,:)
	real	tem
	character(80) form
c
c---( �f�o�C�X�\�� )---
	allocate(temp(0:nx,0:nz))
c	ntemp = 1
	ntenum = ntemax
	btmp = 300.0		!�x�[�X���x[K]
	dtmp = 0.2			!���x���݃�T
c
	tmax = 0.0
	tmin = huge(tmin)
c---	�t�@�C�����牷�x���z�ǂݍ���(n+�̈�͏���)   ---
c---	�G���[��goto 100,�t�@�C��end��goto 200	---
	do iz = 0,nz
	do ix = mxpart1,mxpart2
		read(8,'(f15.7)',ERR=100,END=200)tem
		tmin = min(tmin,tem)
		tmax = max(tmax,tem)
		temp(ix,iz)=tem
	enddo
c
c---	n+�̈�̃f�[�^��[   ---
	do ix = mxpart2+1,nx
		temp(ix,iz)=temp(mxpart2,iz)
	enddo
	do ix = 0,mxpart1-1
		temp(ix,iz)=temp(mxpart1,iz)
	enddo
	enddo
c
	btmp = tmin		!�x�[�X���x[K]
	dtmp = (tmax-tmin)/(ntenum)	!0.2			!���x���݃�T
c
	if(tmax.eq.tmin)then
		form= '(''���x�ω����݂��܂���'')'
		write(*,form)
		write(99,form)
		goto 400
		return
	elseif((tmax.lt.tmin).or.(tmin.lt.0.0))then
		form= '(''���x�ݒ肪���������Ȃ��Ă��܂��B'')'
		write(*,form)
		write(99,form)
		stop
	endif
c
c---	�e�[�u���ϊ�	 ---
	do ix = 0,nx
	do iz = 0,nz
		ntemp(ix,iz) = ifix((temp(ix,iz)-btmp)/dtmp)+1
		ntemp(ix,iz) = min(int2(ntenum),max(int2(1),ntemp(ix,iz)))
	enddo
	enddo
	deallocate(temp)
	return
c
c---	�G���[����	 ---
200	continue
100	continue
	form = "(x,'���x���z�t�@�C�����ǂݍ��߂܂���ł����B')"
	write(*,form)
	write(99,form)
	if(tmax.ge.tmin) btmp = tmin
c
400	continue
	dtmp = 0.0			!���x���݃�T
	ntemp = 1			!ntemp:�S�̔z��
	ntenum = 1
c
	form = "(x,'�S��',F9.4,'K�Ōv�Z���܂��B')"
	write(*,form)btmp
	write(99,form)btmp
c
	return
c	
	end
