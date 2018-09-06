	SUBROUTINE S_EX_CIRCUIT_G(IG1, IDS1, ISS1, VG1, VD1, VS1, CT, t,
     &						icstp,Vd_add,Vg_add,idcount,istp,V_in,
     &						V_out,Wg)
*-------------------回路行列計算------------------------------
*        入力 - -
*             A(N1,N)  R *8  : 2次元配列の係数行列
*             N        I *4  : 行数
*             N1       I *4  : Aの行数
*             EPS      R *8  : 特異性の判定値　通常3.52D-15
*        OUTPUT - -
*             A(N1,N)        : ガウス消去法の結果
*             IP(N)    I *4  : PIVOT番号
*             IER      I *4  : エラーコード
*							 = 0,  FOR NORMAL EXECUTION.
*                              = 1,  FOR SINGULAR MATRIX.
*                              = 3,  FOR INVALID ARGUEMENT.
*-------------------説明------------------------
*入力回路
*-----------------------------------------------
	implicit none
c
	INTEGER N_3,N1_3,N_2,N1_2,N_1,N1_1,EIR,icstp
c
	real*8 EPS
	real R,CG,L1,L2,RG,LG,RD,LD,RS,LS,VDC,VAC,f,CS,V1
	real CT
	real t
	real pi,icstpc
c
	parameter(N_1   = 2)
	parameter(N1_1  = 2)
	parameter(N_2   = 2)
	parameter(N1_2  = 2)
	parameter(N_3   = 4)
	parameter(N1_3  = 4)
c
	parameter(EPS = 3.52e-15)
	parameter(pi  = 3.141592)
c	parameter(icstpc = 1000000)
	parameter(icstpc = 500000)
c
	INTEGER IP_3(N_3+1),IP_2(N_2+1),IP_1(N_1+1)
	real*8 A3(N1_3,N1_3),B3(N_3+1),WK1_3(N_3),WK2_3(N_3)
	real*8 A2(N1_2,N1_2),B2(N_2+1),WK1_2(N_2),WK2_2(N_2)
	real*8 A1(N1_1,N1_1),B1(N_1+1),WK1_1(N_1),WK2_1(N_1)
	real VD2(-1:icstpc)
	real VG2(-1:icstpc)
c--------extinsic------------------------------------------
	real VG1(-1:icstpc),IG1(-1:icstpc),VGL(-1:icstpc)
	real VS1(-1:icstpc),ISS1(-1:icstpc),VSL(-1:icstpc)
	real IDS1(-1:icstpc),VD1(-1:icstpc),VDL(-1:icstpc)
	real ISL(-1:icstpc),ISC(-1:icstpc)
c----------------------------------------------------------
	integer IER
	real V_in,V_out,Vd_add,Vg_add
	integer idcount,istp
	real Wg

c########nari###################################################
	real Vo(-1:icstpc)
	real IC(-1:icstpc)
	real Idd(-1:icstpc)
	real IR(-1:icstpc)
c
	open(unit=310,file='S_ExCirPara.data')
c
c----extrinsic--------------
	read(310,*) RS
	read(310,*) LS
	read(310,*) CS
c--------------------------
c
	close(310)
c#######################################################
c
c-----初期値----------------------------------------------------------------
	IG1(-1)=IG1(0)
c	VG2(-1)= VG1(-1)+RG*IG1(0)
c	VD2(-1)=VD1(-1)+RD*IDS1(0)	
	IDS1(-1)=IDS1(0)
	ISS1(-1)=ISS1(0)
	ISL(-1)=ISS1(0)
	VS1(-1)=ISS1(0)*RS

c####### 04.4.14 nari#####################################################
c	open(unit=311,file='InitiaCir2.data')
cc
c	read(311,*) V2(-1)	!入力回路に使用
c	read(311,*) I2(-1)	!入力回路に使用
c	read(311,*) IR(-1)	!出力回路に使用
c	read(311,*) Idd(-1)	!出力回路に使用
c	read(311,*) Vdd		!出力回路に使用
cc
c	close(311)
c###########nari#################################################


c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c	I2(-1)=Ig(0)
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c-----------------------------------------------------------
c------入力電圧---------------------
c	VD2(icstp)=VD2(icstp-1)
c	VG2(icstp)=VG2(icstp-1)
c	VG2(0)=V_in
c	VD2(0)=V_out
	if((mod(istp,idcount).eq.0).and.(t.ne.0))then
		VD1(icstp)=VD1(icstp-1)+Vd_add
		VG1(icstp)=VG1(icstp-1)+Vg_add
		write(*,*) 'VD2(icstp)=',VD2(icstp)
	endif
c------入力回路---------------------

c############# Gate ######################################################
c	A1(1,1) = 0
c	A1(1,2) = -CT/LG
cc
c	A1(2,1) = -1
c	A1(2,2) = 1
cc
c	B1(1) = -CT/LG*VG2(icstp) +IG1(icstp) -IG1(icstp-1) 	!我流
c	B1(2) = RG*IG1(icstp)
c	
c	CALL GLU1( A1, N_1, N1_1, EPS, WK1_1, IP_1, IER ) 
c	CALL GSLV1( A1, N_1, N1_1, B1, IP_1 )
c---------------------------------------------------
c
c	VG1(icstp) = B1(1)
c	VGL(icstp) = B1(2)
c
c############# Gate END ##################################################

c############# Drain ######################################################
c	A2(1,1) = -CT/LD
c	A2(1,2) = 0
cc
c	A2(2,1) = 1
c	A2(2,2) = -1
cc
c	B2(1) = -CT/LD*VD2(icstp) +IDS1(icstp) -IDS1(icstp-1) 	!我流
c	B2(2) = RD*IDS1(icstp)
cc	
c	CALL GLU1( A2, N_2, N1_2, EPS, WK1_2, IP_2, IER )
c	CALL GSLV1( A2, N_2, N1_2, B2, IP_2 )
cc---------------------------------------------------
c
c	VDL(icstp) = B2(1)
c	VD1(icstp) = B2(2)
c############# Drain END ######################################################

C############# Source ##############################################
	A3(1,1) = 1
	A3(1,2) = 1
	A3(1,3) = 0
	A3(1,4) = 0
c
	A3(2,1) = -RS
	A3(2,2) = 0
	A3(2,3) = 1
	A3(2,4) = -1
c
	A3(3,1) = -1
	A3(3,2) = 0
	A3(3,3) = 0
	A3(3,4) = CT/LS
c
c
	A3(4,1) = 0
	A3(4,2) = CT/CS
	A3(4,3) = -1
	A3(4,4) = 0
c
	B3(1) = ISS1(icstp) 	!我流
	B3(2) = 0
	B3(3) = -ISL(icstp-1) 	!我流
	B3(4) = -VS1(icstp-1)
c	
c	CALL GLU2(A,N,N1,EPS,WK1,WK2,IP,IER)
c	CALL GLU4(A,N,N1,EPS,WK1,WK2,IP,IER)
c	CALL GSLV4(A,N,N1,B,IP)
	CALL GLU1( A3, N_3, N1_3, EPS, WK1_3, IP_3, IER )
	CALL GSLV1( A3, N_3, N1_3, B3, IP_3 )
*---------------------------------------------------
c
	ISL(icstp) = B3(1)
	ISC(icstp) = B3(2)
	VS1(icstp) = B3(3)
	VSL(icstp) = B3(4)
C############# Source END ##############################################
c
	if(icstp.eq.0)then
		write(229,*) 't,VG1,VD1,VS1,VSL'
		write(230,*) 't,IG1,IDS1,ISS1,ISL,ISC'
	endif
c
c	if(mod(icstp,200).eq.0) then
c
2000	format(5(E12.5,','))
3000	format(6(E12.5,','))
      write(229,2000) t*1.0e+12,VG1(icstp),
     &				VD1(icstp),VS1(icstp),VSL(icstp)				!sec -> ps

      write(230,3000) t*1.0e+12,IG1(icstp),IDS1(icstp),ISS1(icstp),
     &				ISL(icstp),ISC(icstp)

	write(401,*) IDS1(icstp)/Wg,'	',IG1(icstp)/Wg,'	',
     &				ISS1(icstp)/Wg			!Ex_current_ypara_NoWg.txt
	write(402,*) IDS1(icstp),'	',IG1(icstp),'	',
     &				ISS1(icstp)			!Ex_current_ypara_kakeruWg.txt


c	endif
c
c
	return
	end
c