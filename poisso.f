
	subroutine poisso(dx, dz, lhet, twodeg,iarea,
     &				vb, vi, lnpole1, lnpole2, melpos,
     &				eps, eg, c_ratio,
     &				u, cn, cp, dopem, maceps,
     &				lxr1,lzr1,lxr2,lzr2,dltec)
      implicit none
c
c===物理定数===
	real	qe
	parameter(qe=1.60219e-19)
c
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c
c	nx,nz	... 格子の数
c	npole   	... 電極の数
c	nlayer   	... エリア数
c	nvalley		... 谷の数
c	lxr1,lzr1	...リセス左端・左上端
c	lxr1,lzr1	...リセス右端・右下端
c
c---デバイス構造---
	real	dx,dz		!メッシュサイズ[m]
	integer(2),	dimension (nrecess)	:: lxr1,lzr1,lxr2,lzr2
	integer(2)	lhet(nlayer)	!ヘテロ界面の深さ[m]
	integer(1)	iarea(nlayer)
c---2DEG電荷量---
	real	twodeg(nlayer)		!2DEGシート電荷量(*1.0e-4 cm^-3)
c---電極---
	real,	dimension(npole)	:: vb,vi			!vb:ビルドイン電圧,vi:電圧
	integer(2),dimension(npole)	:: lnpole1,lnpole2	!電極の幅
	integer(1),dimension(npole)	:: melpos		!
c---領域別パラメータ---
	real,	dimension(narea)	:: eps		!誘電率
	real,	dimension(nvalley,narea):: eg
	real	c_ratio
c---デバイス内状態---
	real,	dimension(0:nx,0:nz)	:: u
	real,	dimension((nx+1)*(nz+1))	:: cn,cp,dopem
c	dopem(j) ...j番目の領域の不純物濃度
c	
c===poisso用変数==========
	real	maceps
c
c====ローカル変数====
c====poisso用変数==========
	integer,save	::	nb1,nb2
	integer,save	::	ierr,m1,nxnz
	integer,save,allocatable	::	ib1(:),ib2(:)
c
	real,save,allocatable	::	bv(:),alp(:),bet(:),p(:),f(:),q(:)

	integer	i,ii,ia,il,ix,iz,intyp
	integer	nxnz2
	integer	lh(0:nlayer)
	real,save	:: vidt(3)
	real	altwodeg
	real,save	:: c_ratio_old
c
c---( 表面電荷 )---	
	real,dimension(0:nrecess) :: vsfr 
	real vsf,vsb
	parameter(vsf=0.4, vsb=0.8)
c---( 裏面電荷 )--- 11/08/06原
	real vsbr
	real, dimension	(nvalley,narea):: dltec
c------------------------------------------
c=========================================================================================================
c
	nxnz=(nx+1)*(nz+1)	!全メッシュ数
	m1=nx+1				!nx+1
c
c---( poisson方程式の定数 )----
	if (.not. allocated(p)) then	!最初だけ実行
		allocate (p(nxnz))
		lh(0)=0;lh(1:nlayer)=lhet(1:nlayer)
		do il = 1, nlayer
			ia = iarea(il)
			p(lh(il-1)*m1+1:lh(il)*m1)	= eps(ia)	!誘電率をpに格納
		enddo
		ia = iarea(nlayer)
		p(lh(nlayer)*m1+1:nxnz)	= eps(ia)
		allocate (f(nxnz))	!(0:nx,0:nz))
		allocate (q(nxnz))
		q=0.0
	endif
c
	f=(dopem-cn+cp)*qe		!f,dopem,cn:全体配列
c
c----------------------------------------------------------
c
c=======( 境界条件{横方向拡張} )======================
	if (.not. allocated(bv)) then	!最初だけ実行
		nxnz2 = (nx+nz)*2	!!!!!!!!!!
		allocate (bv(nxnz2),alp(nxnz2),bet(nxnz2))
		allocate (ib1(nxnz2),ib2(nxnz2))
		vidt=huge(vidt)		!vidt:全体配列
		bv=0;alp=0;bet=0	!初期化
	endif
c
c	境界条件の設定・再設定
c	（電圧に変化があったら実行）
c	if(    (vi(1).ne.vidt(1))
c     &   .or.(vi(2).ne.vidt(2))
c     &   .or.(vi(3).ne.vidt(3))
c     &   .or.(c_ratio.ne.c_ratio_old)
c     &   )then
	if(npole.ge.5)then
		vb(4) = vb(1)
		vi(4) = vi(1)
		vb(5) = vb(3)
		vi(5) = vi(3)
	endif
	nb1 = 0;nb2 = 0
c
c---( 表面電荷 )--- 2DEGより算出 2011/3/25原	
c	Al025 Vg0.0V		9.6263360E+11 cm-2
c	InAs Vg0.0V			6.0841971E+11 cm-2
c
	vsfr(0) = 24.0e15*qe
	vsfr(1) = 24.0e15*qe
	vsfr(2) = 24.0e15*qe
c
c---( 裏面電荷 )--- 2DEGより算出 2011/8/6原
c	Al025 Vg0.0V		2.9901177E+11 cm-2
c	InAs Vg0.0V			4.6512364E+11 cm-2
c	InAs AlSb			8.2250168E+22 cm-2	
c
	vsbr	= 1.0e15*qe
c
c---( 境界条件の設定 )---
	do 10 i=1,nxnz			!全メッシュ数まで
		ix = (mod(i-1,m1))	!あまり分、座標x  (m1=nx+1)
		iz = ((i-1)/m1)		!割った整数分、座標z
c
c	---デバイス左上端・右上端---
		if((i.eq.1).or.(i.eq.m1))then
			if(i.eq. 1)	intyp = 2	!電極左端
			if(i.eq.m1)	intyp = 3	!電極右端
			do ii=1,npole
				if(lnpole1(ii).ne.lnpole2(ii))then	!左端右端等しくないとき
				!!-------リセスがあれば入らない　2017/05/22藤沢 ↓
					if((melpos(ii).eq.1).and.		!ゲートのとき
     &				  (ix.ge.lnpole1(ii)).and.(ix.le.lnpole2(ii)))then
						goto 1000				!ディレクレ境界へ
				!!-------　↑
					elseif((melpos(ii).eq.intyp).and.			
     &				  (iz.ge.lnpole1(ii)).and.(iz.le.lnpole2(ii)))then
						goto 1000				!ディレクレ境界へ
					endif	
				endif
			enddo
			goto 3000
c
c		elseif	(iz.eq.0)then						!上端
c			do ii=1,npole
c				if((melpos(ii).eq.1).and.(lnpole1(ii).ne.lnpole2(ii)).and.
c    &				(ix.ge.lnpole1(ii)).and.(ix.le.lnpole2(ii)))then
c					goto 1000
c				endif
c			enddo
c			goto 3000
c
c	---デバイス上端---
		elseif (iz.eq.0)then
			if ((ix.le.minval(lxr1)).or.(ix.ge.maxval(lxr2)))then
				do ii=1,npole
			if((melpos(ii).eq.1).and.(lnpole1(ii).ne.lnpole2(ii)).and.
     &				(ix.ge.lnpole1(ii)).and.(ix.le.lnpole2(ii)))then
						goto 1000
					endif
				enddo
			elseif((ix.gt.minval(lxr1)).and.(ix.lt.maxval(lxr2)))then
				go to 10
			endif 
			goto 3000
c
c	---リセス内左端---
		elseif ((iz.gt.lzr1(1)).and.(iz.lt.lzr2(1))
     &						   .and.(ix.eq.lxr1(1)))then
			go to 2000
		elseif ((iz.gt.lzr1(2)).and.(iz.lt.lzr2(2))
     &						   .and.(ix.eq.lxr1(2)))then !2段目アウト
c		write(*,*) 'err'
c-----06/11/9 川端　ゲートとリセスの接触のため追加 ---------
			if (lnpole1(2).eq.lxr1(2))then
				ii=2	!要注意　ゲート電極を指定
				goto 1000
			endif
c------------------------------------------------------------
			go to 2000
c
c	---リセス内右端---
		elseif ((iz.gt.lzr1(1)).and.(iz.lt.lzr2(1))
     &						   .and.(ix.eq.lxr2(1)))then
			go to 2000
		elseif ((iz.gt.lzr1(2)).and.(iz.lt.lzr2(2))
     &						   .and.(ix.eq.lxr2(2)))then !2段目アウト
c		write(*,*) 'err'
c-----06/11/9 川端　ゲートとリセスの接触のため追加 -----------
			if (lnpole2(2).eq.lxr2(2))then
				ii=2	!要注意　ゲート電極を指定
				goto 1000
			endif
c-------------------------------------------------------------
			go to 2000
c
c	---リセス内上端(一段目)---
		elseif ((iz.eq.lzr2(1)).and.(ix.ge.lxr1(1))
     &						   .and.(ix.le.lxr1(2)))then
c-----06/11/9 川端　ゲートとリセスの接触のため追加 -----------2段目アウト
			if ((lnpole1(2).eq.lxr1(2)).and.(ix.eq.lnpole1(2)))then
				ii=2	!要注意　ゲート電極を指定
				goto 1000
			endif
c-------------------------------------------------------------			
			go to 5000
		elseif ((iz.eq.lzr2(1)).and.(ix.ge.lxr2(2))
     &						   .and.(ix.le.lxr2(1)))then
c-----06/11/9 川端　ゲートとリセスの接触のため追加 -----------2段目アウト
			if ((lnpole2(2).eq.lxr2(2)).and.(ix.eq.lnpole2(2)))then
				ii=2	!要注意　ゲート電極を指定
				goto 1000
			endif
c-------------------------------------------------------------
			go to 5000
c
c	---リセス内上端(二段目)---
ccc		elseif ((iz.eq.lzr2(1)).and.(ix.ge.lxr1(1)).and.(ix.le.lxr2(1)))then !もともとなし
		elseif ((iz.eq.lzr2(2)).and.(ix.ge.lxr1(2))
     &						   .and.(ix.le.lxr2(2)))then !2段目アウト
			do ii=1,npole
				if((melpos(ii).eq.1).and.(lnpole1(ii).ne.lnpole2(ii)).and.	
     &				(ix.ge.lnpole1(ii)).and.(ix.le.lnpole2(ii)))then
					goto 1000
				endif
			enddo
			go to 6000 
c
c	---デバイス左右端---
		elseif	((ix.eq.0).or.(ix.eq.nx))then
			if(ix.eq. 0)	intyp = 2
			if(ix.eq.nx)	intyp = 3
			do ii=1,npole
				if((melpos(ii).eq.intyp).and.
     &				(iz.ge.lnpole1(ii)).and.(iz.le.lnpole2(ii)))then
					goto 1000
				endif
			enddo
			goto 2000
c
c	---デバイス下端---
		elseif	(iz.eq.nz)then
			goto 4000
c
c	---デバイス中身---
		else
			goto 10
		endif
c
c	---境界条件適用---
c	---Dirichletの境界条件---
 1000		continue
		nb1=nb1+1
		ib1(nb1)=i
		bv(nb1)=-vb(ii)+vi(ii)
c
c	---HEMT用境界---
		if((ix.eq.0).or.(ix.eq.nx))then		!HEMTの左端または右端
			if(iz.le.lhet(1))then			!cap層なら
				ia = iarea(1)
c				bv(nb1)=bv(nb1)+0.7*(eg(1,ia)-eg(1,2)) !バンドギャップ？変更必要
				bv(nb1)=bv(nb1)-dltec(1,ia)		!2017/12/1 鈴木
c		write(*,*) iz,ia, bv(nb1)
c	pause
			else
				do il = 2, nlayer
					if((lhet(il-1).lt.iz).and.(iz.le.lhet(il)))then
						ia = iarea(il)
						bv(nb1)=bv(nb1)-dltec(1,ia)
c						if((ia.eq.4))then
c						  bv(nb1)=bv(nb1)+0.42*(eg(1,ia)-eg(1,2))
c						else 
c						  bv(nb1)=bv(nb1)+0.7*(eg(1,ia)-eg(1,2)) !11/06/29原 
c						endif
c		write(*,*) iz,ia, bv(nb1)
c		pause
					exit
					endif
				enddo
			endif
		endif
		goto 10
c
c	---自然境界条件---
 2000		continue
		nb2=nb2+1
		ib2(nb2)=i
		alp(nb2)=0.0
		bet(nb2)=0.0
		goto 10
c
c	---InGaAs表面電荷(上端)---
 3000		continue
		nb2=nb2+1
		ib2(nb2)=i
		alp(nb2)=0.0
		bet(nb2)= vsfr(0)
c		bet(nb2)= 6.20e16*qe
c		bet(nb2)=sqrt(2.0*p(i)*vsf*qe*1.0e+25)
		goto 10
c
c	---InAlAs表面電荷(リセス内一段目)---
5000		continue
		nb2=nb2+1
		ib2(nb2)=i
		alp(nb2)=0.0
		bet(nb2)= vsfr(1)
c		bet(nb2)= 3.50e16*qe
		goto 10
c
c	---InAlAs表面電荷(リセス内二段目)---
6000		continue
		nb2=nb2+1
		ib2(nb2)=i
		alp(nb2)=0.0
		bet(nb2)= vsfr(2)
c		bet(nb2)= 4.00e16*qe
		goto 10
c
c	---InAlAs裏面電荷(下端)---
 4000		continue
		nb2=nb2+1
		ib2(nb2)=i
		alp(nb2)=0.0
c		bet(nb2)=sqrt(2.0*p(i)*vsb*qe*1.0e+21)
		bet(nb2) = vsbr
		goto 10
c
c	---境界条件を適用しない---
   10	continue
	vidt(1:3)=vi(1:3)
	c_ratio_old = c_ratio 
c
c===================================================
c
c---( poisson方程式のsolverをcall )----
	call poison(nxnz, m1, dx, dz, 
     &			p, q, f, nb1, ib1, bv, nb2, ib2, alp, bet,
     &			u, ierr, maceps, lxr1, lzr1, lxr2, lzr2)

c	do ix=0,nx
c      write(160,1111)(-u(ix,iz),iz=0,nz)
c      end do 
c1111   format (E15.7,400(',',E15.7))	
	return
	end