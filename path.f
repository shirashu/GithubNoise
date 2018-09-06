	subroutine path(dx,z,dhet,cxpole1,cxpole2,
     &			pass,pass_r,x_start_1,x_start_r,xx_goal,scatpoint,			!120817sato
     &			x_mean_free_path_sum,x_mean_free_path_count,
     &			z_start_1,z_start_r,zz_goal,mean_free_path_sum,
     &			mean_free_path_sum2,xcenter)

c---�ϐ��z��p�p�����[�^---
	include 'arraysize.fi'
	common	/arraydata/nx,nz
	integer	nx,nz

c---��{�p�����[�^---
	real	dx
	real	z		!07/03/15
	real	dhet(0:nlayer)						!�w�e���E��(m)
	real	cxpole1(npole),cxpole2(npole)		!�Q�[�g�d��(2)

c-----	�̈�ʉ߂��闱�q���J�E���g120817sato
	real(8),dimension (10) :: pass,pass_r		!�J�E���^�Cr�̓��W�F�N�V�����΍�
	real	x_start_1				!���[�v������ϐ�
	real	x_start_r				!rejection�΍�
	real	xx_goal					!drift���I���ʒu��ۑ�����ϐ�
	integer	scatpoint			!�U���D��U������
	real	x_mean_free_path_sum
	integer x_mean_free_path_count

	real	z_start_1				!���[�v������ϐ�
	real	z_start_r				!rejection�΍�
	real	zz_goal					!drift���I���ʒu��ۑ�����ϐ�
	real	mean_free_path_sum
	real,dimension (0:nx,3) :: 	mean_free_path_sum2		!���ώ��R�s��x���W	 1:�J�E���^,2:x�������ώ��R�s��,3:���ώ��R�s��
	integer xcenter				!x_start��xx_goal�̒��S�_
		
c-----	���[�J���ϐ�  -------
	real	regionA,regionB						!�w�肷��̈�

c==================================================================================

	regionA = cxpole1(2)		!!�Ƃ肠�����Q�[�g�������w��
	regionB = cxpole2(2)

c==================================================================================

c	�h���t�g�J�n�ʒu��������ԂȂ�xx_goal�����ĕ��ώ��R�s�����J�E���g���Ȃ�
	if(x_start_1.eq.0.0)then
		x_start_1 = xx_goal
		x_start_r = x_start_1
		scatpoint = 0			!���ώ��R�s�����J�E���g���郋�[�v���΂�����
	endif
	if(z_start_1.eq.0.0)then
		z_start_1 = zz_goal
		z_start_r = z_start_1
		scatpoint = 0			!���ώ��R�s�����J�E���g���郋�[�v���΂�����
	endif

c-----pass(1)A�O����̈�ʉ�,(2)A�O����̈��,(3)�̈������B�O,(4)�̈������A�O,
c-----pass(5)�̈���A��,(6)�̈�A����B�O,(7)�̈�A����A�O	 (10)�̈���A��(��)
c-----pass(8)B�O����̈�ʉ�,(9)B�O����̈��
					
		!	���ώ��R�s��  !

	if(scatpoint.ge.1)then		!�U�����N�������i��U��(scatpoint=0)�ɂȂ�Ȃ��������j

	if((dhet(nchannel1).lt.z).and.(z.lt.dhet(nchannel2)))then	!channel���Ɍ���

	if((x_start_1.lt.regionA).and.(xx_goal.gt.regionB)) then		!x |AB| xx
		pass(1) = pass(1) +1.0
		pass_r(1) = pass_r(1) +1.0	
c		write(*,*)'pass(1)',pass(1)
c		write(852,*)'pass(1)',pass(1)
c		write(*,*)'pass_r(1)',pass_r(1)
c		write(852,*)'pass_r(1)',pass_r(1)
	endif
	if((x_start_1.lt.regionA).and.(xx_goal.gt.regionA)
     &			.and.(xx_goal.lt.regionB)) then		!x |A xx B|
		pass(2) = pass(2) +1.0
		pass_r(2) = pass_r(2) +1.0	
c		write(*,*)'pass(2)',pass(2)
c		write(852,*)'pass(2)',pass(2)
c		write(*,*)'pass(2)_r',pass_r(2)	
c		write(852,*)'pass_r(2)',pass_r(2)
	
	endif	 
	if((x_start_1.gt.regionA).and.(x_start_1.lt.regionB)
     &			.and.(xx_goal.gt.regionB)) then		!|A x B| xx
			pass(3) = pass(3) +1.0
			pass_r(3) = pass_r(3) +1.0
		if(pass(10).gt.0)then
			pass(6) = pass(6) +pass(10)
			pass_r(6) = pass_r(6) +pass_r(10)
			pass(10) = 0.0
			pass_r(10) = 0.0
		endif
	endif
	if((xx_goal.lt.regionA).and.(x_start_1.gt.regionA)
     &			.and.(x_start_1.lt.regionB)) then		!xx |A x B|
			pass(4) = pass(4) +1.0
			pass_r(4) = pass_r(4) +1.0
		if(pass(10).gt.0)then
			pass(7) = pass(7) +pass(10)
			pass_r(7) = pass_r(7) +pass_r(10)
			pass(10) = 0.0
			pass_r(10) = 0.0
		endif
	endif	
	if((x_start_1.gt.regionA).and.(x_start_1.lt.regionB)
     &		.and.(xx_goal.gt.regionA).and.(xx_goal.lt.regionB)) then		!|A x xx B|
		pass(5) = pass(5) +1.0
		pass_r(5) = pass_r(5) +1.0
		pass(10) = pass(10) +1.0
		pass_r(10) = pass_r(10) +1.0	
	endif	
	
	if((xx_goal.lt.regionA).and.(x_start_1.gt.regionB)) then		!xx |AB| x
		pass(8) = pass(8) +1.0
		pass_r(8) = pass_r(8) +1.0	
	endif
	if((x_start_1.gt.regionB).and.(xx_goal.gt.regionA)
     &			.and.(xx_goal.lt.regionB)) then		!|A xx B| x
		pass(9) = pass(9) +1.0
		pass_r(9) = pass_r(9) +1.0	
	endif
	
		x_mean_free_path_count = x_mean_free_path_count +1			!�J�E���^
		x_mean_free_path_sum = x_mean_free_path_sum+abs(xx_goal-x_start_1)		!���ώ��R�s���̍��v���J�E���^�Ŋ���(output)
		mean_free_path_sum = mean_free_path_sum
     &		+sqrt((xx_goal-x_start_1)**2+(zz_goal-z_start_1)**2)

		xcenter = nint((x_start_1+xx_goal)/dx/2)		!drift���钆�S�_(�����^)
		mean_free_path_sum2(xcenter,1) 
     &			= mean_free_path_sum2(xcenter,1) +1.0		!�J�E���^
		mean_free_path_sum2(xcenter,2) 	
     &			= mean_free_path_sum2(xcenter,2) +abs(xx_goal-x_start_1)
		mean_free_path_sum2(xcenter,3) = mean_free_path_sum2(xcenter,3) 
     &		+sqrt((xx_goal-x_start_1)**2+(zz_goal-z_start_1)**2)

	endif	 !channel���Ɍ���

	x_start_1 = xx_goal		!�U�������ӏ������̃h���t�g�J�n�ʒu�ɂ���
	z_start_1 = zz_goal		!�U�������ӏ������̃h���t�g�J�n�ʒu�ɂ���
		
	endif	!scatpoint

	return
	end