c-----オーミック電極から発生するエネルギーに関するサブルーチン-----
c-----チャネル層の材料だけ-----
	subroutine ohmic_theory(de2,de3,i2max,i3max,E2,E3,n2,n3,aff,am)	!非放物線性
	implicit none
c===引数===
c---変数配列用パラメータ---
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c
	integer i,i2max,i3max
	real Ze1,Ze2,fe						!非放物線性
	real(8) a1,a2,he,ke,am1,am2
	real(8) aff(nvalley,narea)			!非放物線性
	real(8) aff1,aff2
	real E,Tel,Ef,Ec,de2,de3
	real dx
	real Ne1,Ne2						!非放物線性
c
	real b1,b2,b3,b4,c1,c2,c3,c4,d4		!非放物線性
	real	,dimension(:),allocatable :: nn2
	real	,dimension(:),allocatable :: nn3
	real E2(0:i2max),n2(0:i2max)
	real E3(0:i3max),n3(0:i3max)

	real, dimension (nvalley,narea)::am			!120201
c
	real(8) pi,q,h,bk,am0
	parameter(pi = 3.141592654,q = 1.60219e-19)
	parameter(h  = 6.626e-34,bk  = 1.38066e-23)
	parameter(am0 = 9.10953e-31)
c
c-----入力部分InAs(In(0.53)Ga(0.47)As)-----
	Ec=0.00			!Ec[eV]
	Ef=0.050		!Ef[eV] 120126homma
	Tel=300.0		!Tel[K]
	am1=am(1,1)	!param1を引用（コンポジットチャネル上下）


c--------------------------------------
c
c-----非放物線性α-----
c	aff(1,1) = 4.97462067510536 !param1を引用
	aff1=aff(1,1)		!param1を引用（コンポジットチャネル上下）
c----------------------
c
	he=h/q
	ke=bk/q
	a1=(pi/2.0)*(8.0*am1/h*q/h)**(3.0/2.0)
c
	E=0.0
	Ne1=0.0
	Ne2=0.0
	b1=0.0
	b2=0.0
	b3=0.0
	b4=0.0
	c1=0.0
	c2=0.0
	c3=0.0
	c4=0.0
c
c	write(*,*) 'i2max',i2max,'de2',de2
c
	allocate(nn2(0:i2max))
c
	do i=0,i2max
		E=Ec+de2*i
		b1=Ne1
		Ze1=a1*sqrt(E-Ec)
c
c-----Fermi_Dirac Distribution-----
		fe=1.0+dexp((E-Ef)/(ke*Tel))
		fe=1.0d0/fe
c----------------------------------
		Ne1=Ze1*fe
c-----台形法による積分-----
		b2=Ne1
		b3=(b1+b2)*de2/2.0
c
		b4=b4+b3	!b4=積分値
c
		c1=Ne2
c
		Ze2=a1*((E-Ec)**(3.0/2.0))*(5*aff1/2)		!非放物線性
		Ne2=Ze2*fe
c-----台形法による積分-----
		c2=Ne2
		c3=(c1+c2)*de2/2.0
c
		nn2(i)=b3+c3	!出力用c3値確保
c
		c4=c4+c3
c
		d4=b4+c4	!c4=積分値
	enddo
c
c	write(*,*) 'In(0.53)Ga(0.47)As電極濃度=',d4
c
	do i=0,i2max
		E2(i)=Ec+de2*i
		n2(i)=nn2(i)/d4		!粒子数割合出力
c		write(*,*)i,E2(i),n2(i)
	enddo
c
	deallocate(nn2)
c
c
c-----入力部分(InAs)-----
	Ec=0.00			!Ec[eV]
	Ef=0.050		!Ef[eV]
	Tel=300.0		!Tel[K]
	am2=am(1,1)	!param1を引用（コンポジットチャネルチャネル中）		!120201
c------------------------
c
c-----非放物線性α-----
c	aff(1,1) = 4.97462067510536
	aff2=aff(1,1)		!param1を引用（コンポジットチャネルチャネル中）		!120201	
c----------------------
c
	he=h/q
	ke=bk/q
	a1=(pi/2.0)*(8.0*am2/h*q/h)**(3.0/2.0)
c
	E=0.0
	Ne1=0.0
	Ne2=0.0
	b1=0.0
	b2=0.0
	b3=0.0
	b4=0.0
	c1=0.0
	c2=0.0
	c3=0.0
	c4=0.0
c
c	write(*,*) 'i3max',i3max,'de3',de3
c
	allocate(nn3(0:i3max))
c
	do i=0,i3max
		E=Ec+de3*i
		b1=Ne1
		Ze1=a1*sqrt(E-Ec)
c
c-----Fermi_Dirac Distribution-----
		fe=1.0+dexp((E-Ef)/(ke*Tel))
		fe=1.0d0/fe
c----------------------------------
		Ne1=Ze1*fe
c-----台形法による積分-----
		b2=Ne1
		b3=(b1+b2)*de2/2.0
c
		b4=b4+b3	!b4=積分値
c
		c1=Ne2
c
		Ze2=a1*((E-Ec)**(3.0/2.0))*(5*aff2/2)		!非放物線性
		Ne2=Ze2*fe
c-----台形法による積分-----
		c2=Ne2
		c3=(c1+c2)*de3/2.0
c
		nn3(i)=b3+c3	!出力用c3値確保
c
		c4=c4+c3
c
		d4=b4+c4	!c4=積分値
	enddo
c
c	write(*,*) 'In(0.53)Ga(0.47)As電極濃度=',d4
c
	do i=0,i3max
		E3(i)=Ec+de3*i
		n3(i)=nn3(i)/d4		!粒子数割合出力
c		write(*,*)i,E3(i),n3(i)
	enddo
c
	deallocate(nn3)

	return
	end