	subroutine outinit(dt)
c
      implicit none
	common /step/jinit,jstat,jstep,jdisp
	common /outp/joutc,jouti,jouta,jsave,jheat,jcput,joutpi,sw
c
	integer	jinit,jstat,jstep,jdisp
	integer	joutc,jouti,jouta,jsave,jheat,jcput,joutpi,sw(8)
c
	real	ddisp,doutc,douti,douta,dsave,dheat,dcput,doutpi
	real	dt
c
	integer i
c
	integer jdata
c
c k:�o�̓X�C�b�`	�d��, �G��, �|�e���V����, �G�l���M�[, �d�q���x, �Փ˓d��, ���M��,cpu����
c	output.data�ǂݍ���
	read(8,*) (sw(i),	i = 1, 8)
	read(8,*) ddisp			!��ʏo�̓X�e�b�v
	read(8,*) doutc			!�d���o�̓X�e�b�v
	read(8,*) douti			!�o�ߌv�Z�o�̓X�e�b�v
	read(8,*) douta			!���Ϗ�ԏo�̓X�e�b�v
	read(8,*) dsave			!�o�C�i���[�o�̓X�e�b�v
	read(8,*) dheat			!���M���o�̓X�e�b�v
	read(8,*) dcput			!cpu���ԏo�̓X�e�b�v
	read(8,*) doutpi			!�o�ߌv�Z�o�̓X�e�b�v	

	jdisp	= max(0,nint(ddisp/dt))
c	write(*,*)  jdisp
	joutc	= jdata(doutc ,dt)
	jouti	= jdata(douti ,dt)
	jouta	= jdata(douta ,dt)
	jsave	= jdata(dsave ,dt)
	jheat	= jdata(dheat ,dt)
	jcput	= jdata(dcput ,dt)
	joutpi	= jdata(doutpi ,dt)
	end subroutine outinit


	integer function jdata(ddata,dt)
	real ddata,dt

	if(ddata.gt.0.0)then
		jdata	  = max(nint(ddata /dt),1)
	elseif(ddata.lt.0.0)then
		jdata	  = -1
	endif

	end function jdata