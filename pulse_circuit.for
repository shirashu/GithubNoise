	SUBROUTINE pulse_circuit(IG1,IDS1,ISS1,ID1,VO,VG1,V1,VD1,VS1,CT,t,
     &			icstp,Vd_add,Vg_add,idcount,istp,V_in,V_out,Wg,jc_on)
c
c###########################2014/11/10 takahashi########################
	implicit none
c
	INTEGER N_2,N1_2,N_1,N1_1,EIR,icstp,i,jc_on,g_or_d
c
	real*8 EPS
	real ZS,ZO,Lch,Cc
	real CT					!回路シミュレーションのステップ時間dt
	real t					!本回しMCの実時間
	real pi,icstpc
	real Z2					!Sパラ計算用
	real tw,t0,gt,pt		!gaussian palse用変数
c
	parameter(N_1   = 2)	!gateの配列（修正接点格納用）
	parameter(N1_1  = 2)
	parameter(N_2   = 3)	!Drein	
	parameter(N1_2  = 3)
c
	parameter(EPS = 3.52e-15)
	parameter(pi  = 3.141592)
c	parameter(icstpc = 1000000)
	parameter(icstpc = 40000)	!MC出力電流格納用配列（本回しステップ数分）
c###########################################################################
c	修正接点方程式MNAで使用する配列の宣言
c	 IP:ピポッド格納用、WK：行列計算格納用配列　　　　←行列計算
c	 A,B：MNA記述用配列（回路の記述）
c###########################################################################
	INTEGER IP_2(N_2+1),IP_1(N_1+1)		
	real*8 A2(N1_2,N1_2),B2(N_2+1),WK1_2(N_2),WK2_2(N_2)
	real*8 A1(N1_1,N1_1),B1(N_1+1),WK1_1(N_1),WK2_1(N_1)
c--------nodel current,voltage------------------------------------------
	real VG1(-1:icstpc),IG1(-1:icstpc),V1(-1:icstpc)
	real VS1(-1:icstpc),ISS1(-1:icstpc),ID1(-1:icstpc)
	real IDS1(-1:icstpc),VD1(-1:icstpc),IO(-1:icstpc),VO(-1:icstpc)	
c--------Sparameter--------------------------------------------
	real P(-1:icstpc)						!印加するパルス値格納用配列
	real p1(-1:icstpc),p2(-1:icstpc),p3(-1:icstpc),p4(-1:icstpc)
	real S11(-1:icstpc),S12(-1:icstpc),S21(-1:icstpc),S22(-1:icstpc)
c----------current source----------------------------------------------------
	integer	IER
	real V_in,V_out,Vd_add,Vg_add,Vdd,V_io,V_oo		!V_*oのoはold
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
	read(310,*) g_or_d		!!g_or_d=　0:gate,1dreinへパルスを印加
c--------------------------
	close(310)
c#######################################################
c
c=====(ガウシアンパルス入力部)======================
c	t=istp*dt
c	CT=dtc*dt	!!!!!!!!!!!!!!!!!!!!!!!!!!!固定 回路のステップ実時間
	tw=1.00e-12
	t0=3.00e-12
c
	if(icstp.eq.0) then
      write (*, *) 'ガウシアンパルスの設定'
	write (*, *) 'tw=',tw
	write (*, *) 't0=',t0
	write(503,*) 't,Gausswave'
	endif
5000	format(2(E12.5,','))
c
	gt=CT*icstp
	pt=(gt-t0)*(gt-t0)
	P(icstp) = 0.1*exp(-pt/(tw*tw))
c	P(icstp) = 0.0
	P(0)=0.0
c	write(503,2000) gt,P(icstp)
c
c========Drain or Gateにパルスを印加=============================================
	if(g_or_d.eq.0)then
		if(icstp.eq.0)then
			write(*,*)'Gateにパルスを印加'
		endif
	V_io=V_in+P(icstp-1)					!V_io=V_in(icstp-1)
	V_in=V_in+P(icstp)
	endif
c
	if(g_or_d.eq.1)then
		if(icstp.eq.0)then
		write(*,*)'Dreinにパルスを印加'
		endif
	V_oo = V_out + P(icstp-1)				!V_oo:V_out(icstp-1)
	V_out=V_out+P(icstp)					!V_outはVddを含む(直流カットCcの代わり)
	endif
	write(503,2000) gt,P(icstp),V_in,V_out
c-----(Gate)---------------------------------
	VG1(-1)=VG1(0)
	V1(-1)=VG1(0)
	IG1(icstp) = -IG1(icstp)
c	VG1(icstp)=V_in-(ZS*IG1(icstp))	!MCへの次入力が求まる
c
	A1(1,1)=1
	A1(1,2)=0
	A1(2,1)=1
	A1(2,2)=-1
c
	B1(1)=V_in - ZS*IG1(icstp)
	B1(2)=CT/Cc*IG1(icstp) + V1(icstp-1) - VG1(icstp-1)
c
	CALL GLU1( A1, N_1, N1_1, EPS, WK1_1, IP_1, IER )
	CALL GSLV1( A1, N_1, N1_1, B1, IP_1 )
c
c---------------------------------------------------
	V1(icstp) =B1(1)						!MCへの次入力					
	VG1(icstp)=B1(2)
c-----(Drain)---------------------------------------
c#############初期条件#######################	
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
	VD1(icstp)=B2(1)						!MCへの次入力					
	ID1(icstp)=B2(2)
	IO(icstp) =B2(3)
c------------source---------------------------------
	VS1(icstp)=0							!MCへの次入力
c-------------------------------平均化--------------------------------------
c	VG1(icstp)=(VG1(icstp)+VG1(icstp-1))/2
c	VD1(icstp)=(VD1(icstp)+VD1(icstp-1))/2
c	IG1(icstp)=(IG1(icstp)+IG1(icstp-1))/2
c	IDS1(icstp)=(IDS1(icstp)+IDS1(icstp-1))/2
c	ID1(icstp)=(ID1(icstp)+ID1(icstp-1))/2
c	ISS1(icstp)=(ISS1(icstp)+ISS1(icstp-1))/2
c	IO(icstp)=(IO(icstp)+IO(icstp-1))/2
c
c---------(進行波&Spara算出(要変更))----------------------------------------
c	VとIは時間軸がずれているのでカウンタの値をずらしてる
c	S11,S21を求める時はGateにパルスを印加
c	S12,S22を求める時はDreinにパルスを印加
c---------------------------------------------------------------------------
		Z2=1/(2*sqrt(ZO))
c
	i=icstp
c
		if(g_or_d.eq.0)then						!Gateにパルス印加
			p1(i)=Z2*P(i)						!=Z2*V_in
			p2(i)=Z2*(P(i)-ZS*IG1(i))			
c			p3(i)=Z2*((VD1(i)-V_out)+ZO*IO(i))
			p4(i)=Z2*-(ZO*IO(i))
c
			S11(i)=abs( p2(i) / p1(i) )
			S21(i)=abs( p4(i) / p1(i) )
		endif
c
		if(g_or_d.eq.1)then						!Dreinにパルス印加
c			p1(i)=Z2*((VG1(i)-V_in)+ZS*IG1(i))
			p2(i)=Z2*(P(i)-ZS*IG1(i))
			p3(i)=Z2*P(i)
			p4(i)=Z2*-(ZO*IO(i))
c
			S12(i)=abs( p2(i) / p3(i) )
			S22(i)=abs( p4(i) / p3(i) )
		endif
c
1000  format(3(E12.5,','))
2000	format(4(E12.5,','))
3000	format(6(E12.5,','))
4000	format(5(E12.5,','))
c
	if((gt-t0).ge.5e-16)then
c
      write(229,4000) t*1.0e+12,V_in,VG1(icstp)
     &				,VD1(icstp),VS1(icstp)	!sec -> ps
      write(230,3000) t*1.0e+12,IG1(icstp),IDS1(icstp),ID1(icstp),
     &				ISS1(icstp),IO(icstp)
c-------(進行波、反射波書き出し)----------------------------------
c		write(500,*) 't,VG1,VD1,IG1,IDS1,a1,b1,a2,b2'
c		write(501,*) 't,S11,S12,S21,S22'
c	(unit=500,file='a1_b2.csv')
	if(g_or_d.eq.0)then
	write(500,2000) P(icstp),IG1(icstp),ZO*IO(icstp),IO(icstp)
	endif
c
	if(g_or_d.eq.1)then
	write(500,2000) IG1(icstp)*ZS,IG1(icstp),P(icstp),IO(icstp)
	endif
c-------(Spara書き出し)-------------------------------------------
c	(unit=501,file='Spara.csv')
	if(g_or_d.eq.0)then
	write(501,1000) t*1.0e+12,S11(icstp),S21(icstp)
	endif
c
	if(g_or_d.eq.1)then
	write(502,1000) t*1.0e+12,S12(icstp),S22(icstp)
	endif
c-----------------------------------------------------------------
	write(401,*) IDS1(icstp)/Wg,'	',IG1(icstp)/Wg,'	',
     &				ISS1(icstp)/Wg			!current.txt
	write(402,*) IDS1(icstp),'	',IG1(icstp),'	',
     &				ISS1(icstp)			!current.txt
c
	endif
c
	return
	end	
c##################################################################