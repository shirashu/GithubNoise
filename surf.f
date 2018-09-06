c
c---- ( 表面計算用内部副プログラム ) ----
c
	subroutine surf(akx,akz,x,z,kv,jspsum,
     &				cxpole1,cxpole2,melpos,xmax,zmax,
     &				cxr1,czr1,cxr2,czr2) !2011/3/25原
c
	implicit none
	include 'arraysize.fi'
	real akx,akz,x,z
	integer(1) kv
	integer(4) jspsum(npole)
	real cxpole1(npole),cxpole2(npole)
	integer(1)	melpos(npole)
	real xmax,zmax
	integer i
c-----------------------(リセス)------------------------
	integer ii
	real, dimension (nrecess) :: cxr1,czr1,cxr2,czr2
c-------------------------------------------------------
c	real rnd
c	real pi
c	parameter(pi  = 3.14159265358979323846)
c
c----- kl不正(drift)の回避用 2011/03/23原 ---------------
c	integer(1) kl
c	integer(1) klc
c	real dhet(0:nlayer)
c--------------------------------------------------------


	if(kv.eq.0) return
c
c----------------------------(リセス)--------------------------------------
c	デバイス上端側(リセス領域除く)にはみ出し
	if(z.le.0)then
c	  if((x.lt.cxr1(1)).or.(x.gt.cxr2(1)))) then
c
	!SMP$	ASSERT(ITERCNT(1))
		do i=1,npole
c			if(i.eq.2)cycle		!ゲート電極は除く
			if((melpos(i).eq.1).and.
     &				(x.ge.cxpole1(i)).and.(x.le.cxpole2(i)))then
				kv=0
				jspsum(i)=jspsum(i)+1
				return
			endif
		enddo
		z=-z
		z=max(0.0,min(zmax,z))
		akz=-akz
c	デバイス上端側(リセス領域)にはみ出し
c	  else	
c		z=2*czr2(1)-z
c		akz=-akz
c	  endif
c	デバイス下端側にはみ出し
	elseif(z.ge.zmax)then	
		z=zmax-(z-zmax)
		z=max(0.0,min(zmax,z))
		akz=-akz
	end if
c--------------------------------------------------------------------------
c	デバイス左端側にはみ出し
	if(x.le.0.0)then
c
	!SMP$	ASSERT(ITERCNT(1))
		do i=1,npole
			if((melpos(i).eq.2).and.
     &				(z.ge.cxpole1(i)).and.(z.le.cxpole2(i)))then
				kv=0
				jspsum(i)=jspsum(i)+1
				return
			endif
		enddo
		x=-x
		x=max(0.0,min(xmax,x))
		akx=-akx

c	デバイス右端側にはみ出し
	elseif(x.ge.xmax)then
c	
	!SMP$	ASSERT(ITERCNT(1))
		do i=1,npole
			if((melpos(i).eq.3).and.
     &				(z.ge.cxpole1(i)).and.(z.le.cxpole2(i)))then
				kv=0
				jspsum(i)=jspsum(i)+1
				return
			endif
		enddo
		x=xmax-(x-xmax)
		akx=-akx
c----------------------------(リセス)------------------------------------
c	リセス領域にはみ出し
c	06/12/11 ゲート側面も電子を吸収にさせるため変更
c　　　　　　また、ゲート下面からの吸収がないため追加
	else
c	write(*,*) x,z
	  do i=1,nrecess
	    if((x.ge.cxr1(i)).and.(x.le.cxr2(i))
     &		.and.(z.ge.czr1(i)).and.(z.lt.czr2(i))) then
c	if(((x.ge.cxr1(i)).and.(x.le.cxr2(i)))
c     &		.and.((z.ge.czr1(i)).and.(z.lt.czr2(i)))) then
c			
		  if((abs(cxr1(i)-x).le.abs(cxr2(i)-x))
     &		  .and.(abs(cxr1(i)-x).le.abs(czr2(i)-z))) then
			if(cxr1(i).eq.cxpole1(2)) then
			  kv=0
			  jspsum(2)=jspsum(2)+1
c	write(*,*) 'kv=0_1', n	!確認用　2011/03/25原
			  return
			else 
			  x=2*cxr1(i)-x
			  x=max(0.0,min(xmax,x))
			  akx=-akx
c	write(*,*) '-akx1', n	!確認用　2011/03/25原
			endif
c
			elseif((abs(cxr2(i)-x).lt.abs(cxr1(i)-x))
     &		  .and.(abs(cxr2(i)-x).lt.abs(czr2(i)-z))) then
			if(cxr2(i).eq.cxpole2(2)) then
			  kv=0
			  jspsum(2)=jspsum(2)+1
c	write(*,*) 'kv=0_2', n	!確認用　2011/03/25原
			  return
			else 
		      x=2*cxr2(i)-x
			  x=max(0.0,min(xmax,x))
			  akx=-akx
c	write(*,*) '-akx2', n	!確認用　2011/03/25原
			endif

		  elseif((abs(czr2(i)-z).lt.abs(cxr1(i)-x))
     &		  .and.(abs(czr2(i)-z).lt.abs(cxr2(i)-x))) then
			if((x.ge.cxpole1(2)).and.(x.le.cxpole2(2))) then
			  kv=0
			  jspsum(2)=jspsum(2)+1
			  return
c	write(*,*) 'kv=0_3', n	!確認用　2011/03/25原
			else 
c	write(*,*) 'zref_i', n, z, kl	!確認用　2011/03/25原
			  z=2*czr2(i)-z
			  z=max(0.0,min(zmax,z))
			  akz=-akz
cccccc kl不正(drift)の回避策 2011/03/23原 ccccccccccccccccc
c			klc = 0
c				do klc=1,nchannel2
c				  if((z.gt.dhet(klc-1)).and.(z.lt.dhet(klc))) then
c					kl = klc
c					exit
c				  endif
c				enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	write(*,*) 'zref_f', n, z, kl	!確認用　2011/03/25原
			endif
		  endif
		endif
	  enddo
c----------------------------(リセス)------------------------------------
c	リセス領域にはみ出し
c	else
c	  do i=1,nrecess
c		if((x.ge.cxr1(i)).and.(x.le.cxr2(i))
c     &		.and.(z.ge.czr1(i)).and.(z.lt.czr2(i))) then
c
c		  if((abs(cxr1(i)-x).le.abs(cxr2(i)-x))
c     &		  .and.(abs(cxr1(i)-x).le.abs(czr2(i)-z))) then
c			x=2*cxr1(i)-x
c			x=max(0.0,min(xmax,x))
c			akx=-akx
c	
c		  elseif((abs(cxr2(i)-x).lt.abs(cxr1(i)-x))
c     &		  .and.(abs(cxr2(i)-x).lt.abs(czr2(i)-z))) then
c			x=2*cxr2(i)-x
c			x=max(0.0,min(xmax,x))
c			akx=-akx
c
c		  elseif((abs(czr2(i)-z).lt.abs(cxr1(i)-x))
c    &		  .and.(abs(czr2(i)-z).lt.abs(cxr2(i)-x))) then
c			z=2*czr2(i)-z
c			z=max(0.0,min(zmax,z))
c			akz=-akz
c
c		  endif
c		endif
c	  enddo
c------------------------------------------------------------
c		  if(abs(cxr1(i)-x).le.abs(cxr2(i)-x)) then
c			x=2*cxr1(i)-x
c			x=max(0.0,min(xmax,x))
c			akx=-akx
c		  else
c			x=2*cxr2(i)-x
c			x=max(0.0,min(xmax,x))
c			akx=-akx
c		  endif
c		endif
c	  enddo
c--------------------------------------------------------------------------
	end if
c
c	if(z.le.0)then	!デバイス上端側にはみ出し
c
	!SMP$	ASSERT(ITERCNT(1))
c		do i=1,npole
c			if(i.eq.2)cycle		!ゲート電極は除く
c			if((melpos(i).eq.1).and.
c    &				(x.ge.cxpole1(i)).and.(x.le.cxpole2(i)))then
c				kv=0
c				jspsum(i)=jspsum(i)+1
c				return
c			endif
c		enddo
c		z=-z
c		akz=-akz
c	elseif(z.ge.zmax)then	!デバイス下端側にはみ出し
c		z=zmax-(z-zmax)
c		akz=-akz
c	end if
	return
	end subroutine surf