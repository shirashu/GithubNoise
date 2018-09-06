	subroutine condition(dt,np1,bias,dvi,nig,nic,de,jpot)

	include 'arraysize.fi'
	common /step/jinit,jstat,jstep,jdisp
	integer	jinit,jstat,jstep,jdisp

	real	dt
	integer	np1
	real,dimension(npole):: bias,dvi
	integer nig,nic
	real	de(nenergy)
	integer	jpot

	real	dem(nenergy)
	real	dinit,dstat,dstep
	integer i
	character(80) form

c	---(�V�~�����[�V��������)---
	read(8,*) dt		!2.0e-15	!�P�X�e�b�v�̎�����
	read(8,*) dinit		!50.0e-12	!�������ɗv���鎞�ԃX�e�b�v
	read(8,*) dstat		!0�Œ�		!�ϑ��X�^�[�g���ԃX�e�b�v
	read(8,*) dstep		!200.0e-12	!��Ԃ�����̃X�e�b�v��
	read(8,*) jpot		!1			!�|�e���V�����v�Z��
	read(8,*) np1		!400		!�P���b�V��������̗��q��
	read(8,*) (bias(i), i = 1, 3)	!�����lsgd
	read(8,*) (dvi(i), i = 1, 3)	!������
	read(8,*) nig		!5
	read(8,*) nic		!5
	read(8,*) (dem(i), i = 1, nenergy)	!demax

	de	= dem/float(nemax)	!�S�̔z��

	jinit = nint(dinit/dt)
	jstat = nint(dstat/dt)
	jstep = nint(dstep/dt)
	jpot  = max(1,jpot)

	return
	end
