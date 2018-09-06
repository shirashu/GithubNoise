	subroutine exam(dx,dz,dt,d,eps,bktq,hm,jpot,pgm)
      implicit none
c
	include 'arraysize.fi'
c
	real dx,dz,dt,d,eps,bktq,hm
	integer	jpot
c	real,	dimension (nvalley,nenergy,narea)	:: pgm
	real,	dimension (nvalley,nenergy,npart)	:: pgm	!07/8/4 �s�����U��
c
	real*8	omegap,ramdad
	real*8	lmax,vmax
	real	epsdt
c	integer	iv,ien,ia
	integer	iv,ien,ipart	!07/8/4 �s�����U��
	character(80) form
c
	real q,h
	parameter(q   = 1.60219e-19)
	parameter(h   = 1.05459e-34)
c
	vmax	= 1.0e+6	!m/s = *100cm/s(3e5)
c
	omegap	= sqrt(eps*h/hm/q/q/d)
	ramdad	= sqrt(eps*bktq/q/d)
	lmax	= vmax*dt
c
	form="(A2,'���傫�����܂��B',E10.3,'�ȉ��ɂ��Ă�������')"
	if(omegap.lt.dt)then
		write(99,form)'dt',omegap
		write(*,form)'dt',omegap
		stop
	endif
c
	if(ramdad.lt.dx)then
		write(99,form)'dx',ramdad
		write(*,form)'dx',ramdad
		stop
	endif
c
	if(ramdad.lt.dz)then
		write(99,form)'dz',ramdad
		write(*,form)'dz',ramdad
		stop
	endif
c
	if(lmax.gt.dx)then
		write(99,form)'dx',lmax
		write(*,form)'dx',lmax
		stop
	endif
c
	epsdt = epsilon(dt)
	epsdt = (dt/float(jpot))*(epsilon(dt))
	do iv = 1, nvalley
	do ien= 1, nenergy
c	do ia = 1, narea
	do ipart = 1,npart			!07/8/4 �s�����U��
c		�d�ɂł�pgm�͗p�ӂ��Ă��Ȃ��̂ŃT�C�N���I
		if((ipart.eq.6).or.(ipart.eq.7))cycle	!07/8/4 �s�����U��
		if(epsdt.gt.pgm(iv,ien,ipart))then
			form="('pgm(',I2,2(I2,','),')���v�Z�@�C�v�V�����ȉ��ł�')"
c			write(99,form)iv,ien,ia
c			write(* ,form)iv,ien,ia
			write(99,form)iv,ien,ipart
			write(* ,form)iv,ien,ipart
			write(99,*)'�i�v���[�v�ɂȂ�\��������܂�'
			write(* ,*)'�i�v���[�v�ɂȂ�\��������܂�'
			stop
		endif
	enddo
	enddo
	enddo
c
	end
