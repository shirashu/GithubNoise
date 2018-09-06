	subroutine fermi_formula(ka,avsumconc,avsumtei1,ix,iz,iv,x,f,
     &						 ind,aff,flag3,am)		!��������� ind,aff,flag3�ǉ�

	implicit none
c
c---�ϐ��z��p�p�����[�^---
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c
	integer(1)	ka
	real avsumconc(0:nx,0:nz,nvalley)
	real(8) avsumtei1(0:nx,0:nz,nvalley)
	real(8) aff(nvalley,narea)	!���������
	integer	ix,iz,iv
	real x(20),f(20)
	integer	ind,flag3			!���������
c
	real(8) pi,q,h,bk,am0
	parameter(pi = 3.141592654,q  = 1.60219e-19)
	parameter(h  = 1.05459e-34,bk = 1.38066e-23)
	parameter(am0 = 9.10953e-31)
c
c	real am(2)
	real, dimension (nvalley,narea)::am			!120201
c
c-----(�t�F���~�ϕ�)----
	real bkq,h2
	real const1,const2									!���������
	real a1,b1,c1,a2,b2,c2,a3,b3,c3						!���������
	real j1,j2,j3,ita									!���������
	real f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13		!���������
	real Fj1,Fj2,Fj3,Fjn1,Fjn2,Fjn3						!���������
	real N,Energy
c
	if(ka.eq.1)then

c-----InAs -----
c		am(1) = 3.084113076986626E-002*am0	!���J
c		am(2) = 2.87863695866510*am0	!�k�J		
c		write(*,*)'Channel�ޗ��ł͂���܂���B'
c		am(1) = 1.645605041074556E-002*am0	!���J
c		am(2) = 1.97981284185862*am0	!�k�J		
	elseif(ka.eq.2)then
c-----AlInAs ----!11/06/16
c		am(1) = 0.083*am0	!���J
c		am(2) = 0.304*am0	!�k�J
c		am(1) = 5.174696570551261E-002*am0	!���J
c		am(2) = 2.16675363260798*am0	!�k�J
	elseif(ka.eq.3)then
c-----InAs ----!11/06/16
c		am(1) = 3.084113076986626E-002*am0	!���J
c		am(2) = 2.87863695866510*am0	!�k�J	
c		am(1) = 5.908619328810235E-002*am0	!���J
c		am(2) = 2.20669055934314*am0	!�k�J
	elseif(ka.eq.4)then
c-----AlSb (InP)--------------------11/06/06
c		am(1) = 0.08*am0	!���J
c		am(2) = 0.25*am0	!�k�J
c		am(1) = 0.104524959705449*am0	!���J
c		am(2) = 4.66042392191266*am0	!�k�J
		write(*,*)'Channel�ޗ��ł͂���܂���B'
c	elseif(ka.eq.5)then
c-----InGaAs�A-------------------11/06/06
c		am(1) = 4.067762830590998E-002*am0	!param5�����p
c		am(2) = 0.175634779831713*am0	!�k�J
c----------------------------
	else
		write(*,*)'�K����ޗ����݂���܂���B'
		stop
	endif
c
c-----��������W���� ���݂�In(0.53)Ga(0.47)As����(���J)-----
c-----�@���� InAs(���J)�֕ύX11/03/22 ���@------------------
c	aff(1,1)= 2.31681736990881 	!�R�����g�A�E�g120201
c-----------------------------------------------------------
	iv=1
c
	if(x(2).lt.-100)then
		x(2)=-100d0
c		write(*,*)'x(2)�̒l���s���ł�'
	endif
c
	bkq=bk/q		!kB/q
	h2=2*pi*h		!�v�����N�萔�ɕϊ�
c
c	-----�Z�x�v�Z-----
	j1 = 0.5d0		!�t�F���~�ϕ��p�����[�^
c
	a1 = ( 1.0 + 15.0/4.0*(j1+1) + 1.0/40.0*(j1+1)**2.0 ) ** 0.5
	b1 = 1.8 + 0.61*j1
	c1 = 2.0 + ( 2.0-sqrt(2.0) ) * 2.0**(-j1)
c
	const1=(pi/2.0)*((8.0*am(iv,ka))/h2*q/h2)**(3.0/2.0)
c
      j2 = 1.5d0		!�t�F���~�ϕ��p�����[�^
c
	a2 = ( 1.0 + 15.0/4.0*(j2+1) + 1.0/40.0*(j2+1)**2.0 ) ** 0.5
	b2 = 1.8 + 0.61*j2
	c2 = 2.0 + ( 2.0-sqrt(2.0) ) * 2.0**(-j2)
c
	const2=(pi/2.0)*(((8.0*am(iv,ka))/h2*q/h2)**(3.0/2.0))*
     &                                       ((5*aff(1,1))/2)
c
c	-----�Z�x�ɑ΂���Fermi_Integral�J�n-----
	f1 = (j1+1)*(2.0**(j1+1))
	f2 =(b1+x(2)+((((abs(x(2)-b1))**c1)+(a1**c1))**(1.0/c1)))
     &														**(j1+1)
	f3 = exp(-x(2))
	f4 = sqrt(pi) / 2.0
	Fj1= (f1/f2+f3/f4) ** -1.0
c
	f5 = (j2+1)*(2.0**(j2+1))
	f6 =(b2+x(2)+((((abs(x(2)-b2))**c2)+(a2**c2))**(1.0/c2)))
     &														**(j2+1)
	f7 = exp(-x(2))
	f8 = 3.0*sqrt(pi) / 4.0
	Fj2= (f5/f6+f7/f8) ** -1.0
c
	Fjn1 = (bkq*x(1))**(1.5d0) * Fj1
	Fjn2 = (bkq*x(1))**(2.5d0) * Fj2
c
	f(1)=((const1 * Fjn1 + const2 * Fjn2)/avsumconc(ix,iz,iv))-1
c
c	-----�G�l���M�[�v�Z-----
	j3 = 2.5d0		!�t�F���~�ϕ��p�����[�^
c
	a3 = ( 1.0 + 15.0/4.0*(j3+1) + 1.0/40.0*(j3+1)**2.0 ) ** 0.5
	b3 = 1.8 + 0.61*j3
	c3 = 2.0 + ( 2.0-sqrt(2.0) ) * 2.0**(-j3)
c
c	-----�G�l���M�[�ɑ΂���Fermi_Integral�J�n-----
	f9 = (j3+1)*(2.0**(j3+1))
	f10 =(b3+x(2)+((((abs(x(2)-b3))**c3)+(a3**c3))**(1.0/c3)))
     &													**(j3+1)
	f11 = exp(-x(2))
	f12 = 15.0*sqrt(pi) / 8.0
c
	f13 = f9/f10+f11/f12
c
c-----stop���Ȃ����߂̏���-----
	if(f13.le.0)then
		ind = 11
		flag3 = 1
		return
	endif
c------------------------------
c
	Fj3= 1/f13
c
	Fjn3 = (bkq*x(1))**(3.5d0) * Fj3
      Energy = const1 * Fjn2 + const2 * Fjn3
c
	f(2)=(Energy/avsumconc(ix,iz,iv))-avsumtei1(ix,iz,iv)
c
	return
	end