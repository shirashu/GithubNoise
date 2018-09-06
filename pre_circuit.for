	SUBROUTINE pre_circuit(IG1,IDS1,ISS1,ID1,VO,VG1,V1,VD1,VS1,CT,t,
     &			icstp,Vd_add,Vg_add,idcount,istp,V_in,V_out,Wg,jc_on)
c
c###########################2014/11/10 takahashi########################
	implicit none
c
	INTEGER N_2,N1_2,N_1,N1_1,EIR,icstp,i,jc_on,g_or_d
c
	real*8 EPS
	real ZS,ZO,Lch,Cc
	real CT					!��H�V�~�����[�V�����̃X�e�b�v����dt
	real t					!�{��MC�̎�����
	real pi,icstpc
	real Z2					!S�p���v�Z�p
	real tw,t0,gt,pt		!gaussian palse�p�ϐ�
c
	parameter(N_1   = 2)	!gate�̔z��i�C���ړ_�i�[�p�j
	parameter(N1_1  = 2)
	parameter(N_2   = 3)	!Drein	
	parameter(N1_2  = 3)
c
	parameter(EPS = 3.52e-15)
	parameter(pi  = 3.141592)
c	parameter(icstpc = 10000)
	parameter(icstpc = 40000)	!MC�o�͓d���i�[�p�z��i�{�񂵃X�e�b�v�����j
c###########################################################################
c	�C���ړ_������MNA�Ŏg�p����z��̐錾
c	 IP:�s�|�b�h�i�[�p�AWK�F�s��v�Z�i�[�p�z��@�@�@�@���s��v�Z
c	 A,B�FMNA�L�q�p�z��i��H�̋L�q�j
c###########################################################################
	INTEGER IP_2(N_2+1),IP_1(N_1+1)		
	real*8 A2(N1_2,N1_2),B2(N_2+1),WK1_2(N_2),WK2_2(N_2)
	real*8 A1(N1_1,N1_1),B1(N_1+1),WK1_1(N_1),WK2_1(N_1)
c--------nodel current,voltage------------------------------------------
	real VG1(-1:icstpc),IG1(-1:icstpc),V1(-1:icstpc)
	real VS1(-1:icstpc),ISS1(-1:icstpc),ID1(-1:icstpc)
	real IDS1(-1:icstpc),VD1(-1:icstpc),IO(-1:icstpc),VO(-1:icstpc)	
c--------Sparameter--------------------------------------------
	real P(-1:icstpc)						!�������p���X�l�i�[�p�z��
	real p1(-1:icstpc),p2(-1:icstpc),p3(-1:icstpc),p4(-1:icstpc)
	real S11(-1:icstpc),S12(-1:icstpc),S21(-1:icstpc),S22(-1:icstpc)
c----------current source----------------------------------------------------
	integer	IER
	real V_in,V_out,Vd_add,Vg_add,Vdd,V_io,V_oo		!V_*o��o��old
	integer idcount,istp
	real Wg
c----------bias circuit ----------------------------------------------------- 
	open(unit=310,file='ExCirPara.data')
c
c----extrinsic--------------
	read(310,*) V_in
	read(310,*) V_out
	read(310,*) Vdd
	read(310,*) ZS
	read(310,*) ZO
	read(310,*) Lch
	read(310,*) Cc
	read(310,*) g_or_d		!!g_or_d=�@0:gate,1drein�փp���X�����
c--------------------------
	close(310)
c#######################################################
c
c=====(�K�E�V�A���p���X���͕�)======================
c	t=istp*dt
c	CT=dtc*dt	!!!!!!!!!!!!!!!!!!!!!!!!!!!�Œ� ��H�̃X�e�b�v������
	tw=1.00e-12
	t0=3.00e-12
c
	if(icstp.eq.0) then
      write (*, *) '�K�E�V�A���p���X�̐ݒ�'
	write (*, *) 'tw=',tw
	write (*, *) 't0=',t0
	write(503,*) 't,Gausswave'
	endif
5000	format(2(E12.5,','))
c
	gt=CT*(icstp+1)
	pt=(gt-t0)*(gt-t0)
	P(icstp) = 0.1*exp(-pt/(tw*tw))
c	P(icstp) = 0.0
	P(0)=0.0
c	write(503,2000) gt,P(icstp)
c	
c========Drain or Gate�Ƀp���X�����=============================================
	if(g_or_d.eq.0)then
		if(icstp.eq.0)then
			write(*,*)'Gate�Ƀp���X�����'
		endif
	V_io=V_in+P(icstp-1)					!V_io=V_in(icstp-1)
	V_in=V_in+P(icstp)
	endif
c
	if(g_or_d.eq.1)then
		if(icstp.eq.0)then
		write(*,*)'Drein�Ƀp���X�����'
		endif
	V_oo = V_out + P(icstp-1)				!V_oo:V_out(icstp-1)
	V_out=V_out+P(icstp)
	endif
c
c-----(Gate)---------------------------------
	VG1(-1)=VG1(0)
	V1(-1)=VG1(0)
	IG1(icstp) = -IG1(icstp)
c	VG1(icstp)=V_in-(ZS*IG1(icstp))	!MC�ւ̎����͂����܂�
c
	A1(1,1)=1
	A1(1,2)=0
	A1(2,1)=1
	A1(2,2)=-1
c
	B1(1)=V_in-ZS*IG1(icstp)
	B1(2)=CT/Cc*IG1(icstp)+V1(icstp-1)-VG1(icstp-1)
c
	CALL GLU1( A1, N_1, N1_1, EPS, WK1_1, IP_1, IER )
	CALL GSLV1( A1, N_1, N1_1, B1, IP_1 )
c
c---------------------------------------------------
	V1(icstp) =B1(1)								!MC�ւ̎�����					
	VG1(icstp)=B1(2)
c-----(Drain)---------------------------------------
c#############��������#######################	
	ID1(-1)=IDS1(0)
	VD1(-1)=Vdd
c############################################
	A2(1,1)=1
	A2(1,2)=0
	A2(1,3)=ZO
	A2(2,1)=1
	A2(2,2)=Lch/CT
	A2(2,3)=0
	A2(3,1)=0
	A2(3,2)=1
	A2(3,3)=1
c
	B2(1)=V_out
	B2(2)=Vdd + Lch/CT*ID1(icstp-1)
	B2(3)=IDS1(icstp)
c
	CALL GLU1( A2, N_2, N1_2, EPS, WK1_2, IP_2, IER )
	CALL GSLV1( A2, N_2, N1_2, B2, IP_2 )

c---------------------------------------------------
	VD1(icstp)=B2(1)								!MC�ւ̎�����					
	ID1(icstp)=B2(2)
	IO(icstp) =B2(3)
c------------source---------------------------------
	VS1(icstp)=0									!MC�ւ̎�����
	ISS1(icstp)=ISS1(icstp)
	
	if(icstp.eq.0)then					
		write(229,*) 't,Vin,VG1,VD1,VS1'
		write(230,*) 't,IG1,IDS1,ID1,ISS1,IO'
		write(500,*) 'V1,I1,V2,I2'
			if(g_or_d.eq.0)then
				write(501,*) 't,S11,S21'
			endif
			if(g_or_d.eq.1)then
				write(502,*) 't,S12,S22'
			endif
c---------(�i�s�g&Spara�Z�o(�v�ύX))----------------------------------------
c	V��I�͎��Ԏ�������Ă���̂ŃJ�E���^�̒l�����炵�Ă�
c	S11,S21�����߂鎞��Gate�Ƀp���X�����
c	S12,S22�����߂鎞��Drein�Ƀp���X�����
c---------------------------------------------------------------------------
		Z2=1/(2*sqrt(ZO))
c
	i=icstp
c
		if(g_or_d.eq.0)then						!Gate�Ƀp���X���
			p1(i)=Z2*P(i)						!=Z2*V_in
			p2(i)=Z2*(P(i)-ZS*IG1(i))			
c			p3(i)=Z2*((VD1(i)-V_out)+ZO*IO(i))
			p4(i)=Z2*-(ZO*IO(i))
c
			S11(i)=abs( p2(i) / p1(i) )
			S21(i)=abs( p4(i) / p1(i) )
		endif
c
		if(g_or_d.eq.1)then						!Drein�Ƀp���X���
c			p1(i)=Z2*((VG1(i)-V_in)+ZS*IG1(i))
			p2(i)=Z2*(P(i)-ZS*IG1(i))
			p3(i)=Z2*P(i)
			p4(i)=Z2*-(ZO*IO(i))
c
			S12(i)=abs( p2(i) / p3(i) )
			S22(i)=abs( p4(i) / p3(i) )
		endif
c
2000	format(4(E12.5,','))
3000	format(6(E12.5,','))
4000	format(5(E12.5,','))
c
c      write(229,4000) t*1.0e+12,V_in,VG1(icstp)
c     &				,VD1(icstp),VS1(icstp)	!sec -> ps
c      write(230,3000) t*1.0e+12,IG1(icstp),IDS1(icstp),ID1(icstp),
c     &				ISS1(icstp),IO(icstp)
c-------(�i�s�g�A���˔g�����o��)----------------------------------
c		write(500,*) 't,VG1,VD1,IG1,IDS1,a1,b1,a2,b2'
c		write(501,*) 't,S11,S12,S21,S22'
c	(unit=500,file='a1_b2.csv')
c	if(g_or_d.eq.0)then
c	write(500,2000) P(icstp),IG1(icstp),ZO*IO(icstp),IO(icstp)
c	endif
c
c	if(g_or_d.eq.1)then
c	write(500,2000) IG1(icstp)*ZS,IG1(icstp),P(icstp),IO(icstp)
c	endif
c-------(Spara�����o��)-------------------------------------------
c	(unit=501,file='Spara.csv')
c	if(g_or_d.eq.0)then
c	write(501,4000) t*1.0e+12,S11(icstp),S21(icstp)
c	endif
c
c	if(g_or_d.eq.1)then
c	write(502,4000) t*1.0e+12,S12(icstp),S22(icstp)
c	endif
c-----------------------------------------------------------------
c	write(401,*) IDS1(icstp)/Wg,'	',IG1(icstp)/Wg,'	',
c     &				ISS1(icstp)/Wg			!current.txt
c	write(402,*) IDS1(icstp),'	',IG1(icstp),'	',
c     &				ISS1(icstp)			!current.txt
c
	endif
c
	return
	end	
c##################################################################
