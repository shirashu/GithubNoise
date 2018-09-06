	subroutine fermi_calc_table(ntab1,etab1,aff,am)	!非放物線性
	implicit none
c=====引数=====
c-----変数配列用パラメータ-----
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c-----ローカル変数-----
	integer i,iv,k
	integer(1) ka
c-----有効質量-----
c	real am(1)
	real, dimension (nvalley,narea)::am			!120201
	real(8) aff(nvalley,narea)	!非放物線性
c-----メッシュ&谷別電子濃度-----
	real ntab1(0:300000,5)		!濃度(番号、材料)
	real fermi1(0:300000,5)		!Ef(番号、材料)
	real etab1(0:300000,5)		!Energy(番号、材料)
c-----フェルミ積分-----
	real bkq,h2
	real const1,const2,tell							!非放物線性
	real a1,b1,c1,a2,b2,c2,a3,b3,c3					!非放物線性
	real j1,j2,j3,ita,ita1							!非放物線性
	integer it
	real f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12		!非放物線性
	real Fj1,Fj2,Fj3,Fjn1,Fjn2,Fjn3					!非放物線性
	real efn(300000),ntab(300000)
	real ee
c-----基本パラメータ-----
	real(8) pi,q,h,bk,am0
	parameter(pi=3.141592654,q=1.60219e-19)
	parameter(h=1.05459e-34, bk  = 1.38066e-23)
	parameter(am0 = 9.10953e-31)
c
	tell=300.0d0
c
c	do ka= 1,2				!材料loop 2011/03/22 原
        ka = 1			!2011/03/22
	 do i=0,100000,10	!濃度loop
c
c-----材料別有効質量-----
c       am(1) = 3.084113076986626E-002*am0	!In(0.53)Ga(0.47)As Γ谷 !param5を引用
c	if(ka.eq.1)then
c		am(1) = 1.645605041074556E-002*am0	!2011/03/22
c	elseif(ka.eq.2)then		!2011/03/22
c	    am(1) = 0.083*am0	!2011/06/30
c	elseif(ka.eq.3)then
c		am(1) = 0.035*am0	! Γ谷
c	elseif(ka.eq.4)then
c		am(1) = 0.08*am0	!InP Γ谷
c	else
c		write(*,*)'適する材料がありません'
c		stop
c	endif					!2011/03/22
c------------------------
c
c-----非放物線性α-----
c	aff(1,1) = 2.31681736990881  !param1を引用
c----------------------

c
	ntab1(i,ka) = 1e22 + 1e21*i
c
	iv = 1
c
	bkq=bk/q		!kB/q
	h2=2*pi*h		!プランク定数に直す
c
	j1 = 0.5d0		!フェルミ積分パラメータ
c
	a1 = ( 1.0 + 15.0/4.0*(j1+1) + 1.0/40.0*(j1+1)**2.0 ) ** 0.5
	b1 = 1.8 + 0.61*j1
	c1 = 2.0 + ( 2.0-sqrt(2.0) ) * 2.0**(-j1)
c
	const1=(pi/2.0)*((8.0*am(1,1))/h2*q/h2)**(3.0/2.0)
c
      j2 = 1.5d0		!フェルミ積分パラメータ
c
	a2 = ( 1.0 + 15.0/4.0*(j2+1) + 1.0/40.0*(j2+1)**2.0 ) ** 0.5
	b2 = 1.8 + 0.61*j2
	c2 = 2.0 + ( 2.0-sqrt(2.0) ) * 2.0**(-j2)
c
	const2=(pi/2.0)*(((8.0*am(1,1))/h2*q/h2)**(3.0/2.0))*
     &								(5*aff(1,1)/2)		!非放物線性
c
c-----ita<0の場合の濃度テーブル算出 dita=0.01-----
		k=1
		do ita= -150.0 , 0.0 , 0.01
			efn(k) = ita * bkq * tell
c-----フェルミ積分Fjを求める-----
			f1 = (j1+1)*(2.0**(j1+1))
			f2 =(b1+ita+((((abs(ita-b1))**c1)+(a1**c1))**(1.0/c1)))
     &														**(j1+1)
			f3 = exp(-ita)
			f4 = sqrt(pi) / 2.0
			Fj1= (f1/f2+f3/f4) ** -1.0
c
			f5 = (j2+1)*(2.0**(j2+1))
			f6 =(b2+ita+((((abs(ita-b2))**c2)+(a2**c2))**(1.0/c2)))
     &														**(j2+1)
			f7 = exp(-ita)
			f8 = 3.0*sqrt(pi) / 4.0
			Fj2= (f5/f6+f7/f8) ** -1.0
c-----Fjをもとに濃度テーブルを求める-----
			Fjn1 = (bkq*tell)**(1.5d0) * Fj1	!Nのとき
			Fjn2 = (bkq*tell)**(2.5d0) * Fj2	!Nのとき
c
			ntab(k) = const1 * Fjn1 + const2 * Fjn2
c
				if(ntab(k).ge.ntab1(i,ka))then
			fermi1(i,ka) = ( (ntab1(i,ka)-ntab(k-1))/
     &				(ntab(k)-ntab(k-1))*(efn(k)-efn(k-1)) )+efn(k-1)
				go to 1111
				endif
			k=k+1
		enddo
c
c-----ita>0の場合の濃度テーブル算出 dita=0.1-----
		do ita= 0.1 , 1000.0 , 0.1
			efn(k) = ita * bkq * tell
c-----フェルミ積分Fjを求める-----
			f1 = (j1+1)*(2.0**(j1+1))
			f2 =(b1+ita+((((abs(ita-b1))**c1)+(a1**c1))**(1.0/c1)))
     &														**(j1+1)
			f3 = exp(-ita)
			f4 = sqrt(pi) / 2.0
			Fj1= (f1/f2+f3/f4) ** -1.0
c
			f5 = (j2+1)*(2.0**(j2+1))
			f6 =(b2+ita+((((abs(ita-b2))**c2)+(a2**c2))**(1.0/c2)))
     &														**(j2+1)
			f7 = exp(-ita)
			f8 = 3.0*sqrt(pi) / 4.0
			Fj2= (f5/f6+f7/f8) ** -1.0
c-----Fjをもとに濃度テーブルを求める-----
			Fjn1 = (bkq*tell)**(1.5d0) * Fj1	!Nのとき
			Fjn2 = (bkq*tell)**(2.5d0) * Fj2	!Nのとき
c
			ntab(k) = const1 * Fjn1 + const2 * Fjn2
c
				if(ntab(k).ge.ntab1(i,ka))then
			fermi1(i,ka) = ( (ntab1(i,ka)-ntab(k-1))/
     &				(ntab(k)-ntab(k-1))*(efn(k)-efn(k-1)) )+efn(k-1)
				go to 1111
				endif
			k=k+1
		enddo
 1111 continue
c
c-----エネルギー計算-----
	j3 = 2.5d0		!フェルミ積分パラメータ
c
	a3 = ( 1.0 + 15.0/4.0*(j3+1) + 1.0/40.0*(j3+1)**2.0 ) ** 0.5
	b3 = 1.8 + 0.61*j3
	c3 = 2.0 + ( 2.0-sqrt(2.0) ) * 2.0**(-j3)
c
	ita1=fermi1(i,ka)/(bkq*tell)
c
c-----Fermi_Integral開始-----
		    f5 = (j2+1)*(2.0**(j2+1))
			f6 =(b2+ita1+((((abs(ita1-b2))**c2)+(a2**c2))**(1.0/c2)))
     &														**(j2+1)
			f7 = exp(-ita1)
			f8 = (3.0d0*sqrt(pi)) / 4.0
			Fj2= (f5/f6+f7/f8) ** -1.0
			f9 = (j3+1)*(2.0d0**(j3+1))
		f10 =(b3+ita1+((((abs(ita1-b3))**c3)+(a3**c3))**(1.0d0/c3)))
     &													**(j3+1)
			f11 = exp(-ita1)
			f12 = (15.0d0*sqrt(pi)) / 8.0d0
			Fj3= (f9/f10+f11/f12) ** -1.0d0
c
			Fjn2 = (bkq*tell)**(2.5d0) * Fj2
			Fjn3 = (bkq*tell)**(3.5d0) * Fj3
			ee = const1 * Fjn2 + const2 * Fjn3
			etab1(i,ka) = ee / ntab1(i,ka)
c-----エネルギー計算終了-----
c-----書き出すのは番号、材料、Γ谷の濃度テーブル、フェルミレベル、電子温度、平均エネルギー
c	write(*,*) i,ka,ntab1(i,ka),fermi1(i,ka),tell,etab1(i,ka)
	enddo
c	enddo	!2011/03/22
c
	return
	end