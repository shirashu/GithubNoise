c
c-----粒子が境界面で反射または進入するイベント-----
	subroutine border_roughness2(
     &					hhm,af,af2,af4,eg,
     &					akx,aky,akz,kv,kl,kl2,ka,iarea,lhet,iiz,iix, !2006/12/09 Hara
     &                    xxx,vvv,basho_reflection,
     &                      basho_roughness,
     &                    split,delta,lambda,count_reflection,
     &                    count_roughness,average,
     &					epA,u)			!ラフネス散乱用ep障壁	
	implicit none
c
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
	real,	dimension (nvalley,narea)	:: hhm,af,af2,af4,eg
	real	akx,aky,akz
	integer(1)	kv,kl,kl2,ka,iiz	!2006/12/09 Hara
	integer		iix
	integer(1)	iarea(nlayer)
	integer(2)	lhet(nlayer)	!2006/12/09 Hara
c
c	---	ローカル変数	---
c	skx,sky,skz:各方向波数の二乗を格納する変数
c	vb:障壁高さを格納する変数
c	ex,ey,ez:各方向運動エネルギーを格納する変数
	real(8)	skx,sky,skz
	real	vb
	real(8)	ex,ey,ez
	integer(1)	ka2
c----ラフネス散乱----------------------------------------------
      real rA,rB,INC,OUT,Pspe,Pdif,aaa,bbb,ccc,ddd,
     &	 delta,lambda,kakudo_OUT,kakudo_INC,zure,average
	real	kkk,kkk2
	integer nang,split,count_reflection,count_roughness,
     &        vvv,xxx,
     &		count_roughness1 
	integer,dimension (0:nx,2) :: basho_roughness,basho_reflection
	real(8) pi
	parameter(pi = 3.141592654)
	real rnd

c----ep-cp障壁用---(ラフネス散乱)
	real	vbcp	!cp障壁高さ
	real,	dimension (0:nx,0:nz)	:: epA
	real,	dimension (0:nx,0:nz)	:: u
	real	ek_epcp	!CPとEPの差分を運動エネルギーとして足す(ラフネス散乱レート)
	real	sk_epcp	!運動エネルギーから求める波数kの二乗
	real	ak_epcp	!運動エネルギーから求める波数k	 (x,y,zどの方向でも等しい)

c
c=================================================================================================
c
c	--- z方向の運動エネルギー計算 ---
	if(af4(kv,ka).ne.0.0)then
		skz = dble(akz)*dble(akz)
		ez=(sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*skz)-1.0)/af2(kv,ka)		!eV
	else
		ez=hhm(kv,ka)*akz*akz
	endif
c	---------------------------------
c
	ka2 = iarea(kl2)
c	vb =0.7*(eg(kv,ka2)-eg(kv,ka))		!Vb:障壁高さ[eV]
c	if(iz.lt.(nlayer-2)-5) then	!ドリフト前位置がclassical領域のとき(2006/08/08改良)!2006/12/09 Hara
c	if(iz.lt.(nlayer-4)-4) then	!ドリフト前位置がclassical領域のとき(2006/08/08改良)!2006/12/22 Hara
c		vb =0.7*(eg(kv,ka2)-eg(kv,ka))		!Vb:障壁高さ[eV]
c	else						!ドリフト前位置がEP領域のとき
c		vb = 0.0			!なぜ 0 eV ?? 2011/03/23原
c	endif
c
c----(cp-ep)のヘテロ障壁------------
	vbcp =0.7*(eg(kv,ka2)-eg(kv,ka))		!Vb:障壁高さ[eV]	
	open(unit=1212,file='zztriela.txt')	!確認用!	

	if(iiz.eq.lhet(nchannel1))then
		vb = -(u(iix,iiz) - epA(iix,iiz)) !cp-epの障壁
	endif
	if(iiz.eq.lhet(nchannel2))then
		vb = -(u(iix,iiz+1) - epA(iix,iiz)) !cp-epの障壁
	endif
	if(vb.lt.0)then	!電極付近のCP,EPは不安定　エラー回避sato 要らないかも
		vb=0
	endif


c----(ヘテロ障壁のラフネス散乱に足すep-cp差分の運動エネルギー)
	if(iiz.eq.lhet(nchannel1))then
		ek_epcp = -(epA(iix,iiz)-(u(iix,iiz+1)-vbcp))
	endif

	if(iiz.eq.lhet(nchannel2))then
		ek_epcp = -(epA(iix,iiz)-(u(iix,iiz)-vbcp)) 
	endif

	if(ek_epcp.lt.0)then	!電極付近のCP,EPは不安定　エラー回避sato
		ek_epcp=0
		vb=0
	endif
	
	if(af4(kv,ka).ne.0.0)then
		sk_epcp = ((af2(kv,ka)*ek_epcp+1)**2-1)/(af4(kv,ka)*hhm(kv,ka))
	else
		sk_epcp = ek_epcp/hhm(kv,ka)		
	endif
		
		ak_epcp = sqrt(sk_epcp/3)		!ラフネス散乱レートに足す波数


	if(ez.gt.vb)then					!障壁Vbより粒子のエネルギーが 高いか
c	===	障壁を乗り越える場合（進入）	===
c		ez = ez-vb			!z方向のエネルギーから障壁エネルギーをひく
		skx = dble(akx)*dble(akx)
		sky = dble(aky)*dble(aky)
		if(af4(kv,ka).ne.0.0)then
			ex=(sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*skx)-1.0)/af2(kv,ka)		!eV
			ey=(sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*sky)-1.0)/af2(kv,ka)		!eV
		else
			ex=hhm(kv,ka)*skx		!αが０の時は単純計算
			ey=hhm(kv,ka)*sky
		endif
c
		kl = kl2
		ka = ka2
c		---	求めたエネルギーから波数決定（α＝０でも計算可能）
		akx=sign(sngl(sqrt(ex*(1.0+af(kv,ka2)*ex)/hhm(kv,ka2))),akx)
		aky=sign(sngl(sqrt(ey*(1.0+af(kv,ka2)*ey)/hhm(kv,ka2))),aky)
		akz=sign(sngl(sqrt(ez*(1.0+af(kv,ka2)*ez)/hhm(kv,ka2))),akz)
c		---	sign(a,b) = aの絶対値にbの符号を掛けた値
c
	else
c	===	障壁を乗り越えない場合（反射）	===
          delta  = 6.5e-10		!凸凹高さ	格子定数の整数倍
		lambda = 200.0e-10		!凸凹広がり
		split  = 100000		!反射角分割数(akx,akzの精度100万に合わせている)→時間短縮のため10万

		INC  = atan(abs(akx)/abs(akz))
		kkk  = sqrt((akx*akx)+(akz*akz))
		kkk2  = sqrt((akx+ak_epcp)**2+(akz+ak_epcp)**2)		!ラフネス散乱レート計算でcpとepの差分を足している
		Pspe = exp(-4*delta*delta*kkk2*kkk2*cos(INC))
		rA   = rnd()

	    kakudo_INC = INC*180/pi

c	    write(*,*)
c	    write(*,*)
c         write(*,*)
c         write(*,*)
c         write(*,*)
c	    write(*,*) '(iix,iiz)=','(',iix,',',iiz,')' 
c         write(*,*) 'Pspe=',Pspe
c         write(*,*) 'kkk=',kkk
	
c---------反射回数カウント-----------------------------
	
	    vvv = iix

		if(iiz.eq.lhet(nchannel1)) then
			if(akz.le.0)then	!チャネル→ヘテロ障壁の場合のみカウント
		    count_reflection = count_reflection + 1
			basho_reflection(vvv,1) = basho_reflection(vvv,1) + 1
			endif
		endif
		if(iiz.eq.lhet(nchannel2)) then
			if(akz.ge.0)then	!チャネル→ヘテロ障壁の場合のみカウント
			count_reflection = count_reflection + 1
			basho_reflection(vvv,2) = basho_reflection(vvv,2) + 1
			endif
		endif


c---------ラフネス散乱回数カウント-----------------------------
		if(rA.gt.Pspe) then
	        xxx = iix

		if(iiz.eq.lhet(nchannel1)) then
			if(akz.le.0)then	!チャネル→ヘテロ障壁の場合のみカウント
			count_roughness = count_roughness + 1
			basho_roughness(xxx,1) = basho_roughness(xxx,1) + 1
			endif
		endif	
		if(iiz.eq.lhet(nchannel2)) then
			if(akz.ge.0)then	!チャネル→ヘテロ障壁の場合のみカウント	
			count_roughness = count_roughness + 1
			basho_roughness(xxx,2) = basho_roughness(xxx,2) + 1
			endif
		endif	
		
		count_roughness1 = count_roughness1 +1


        	           
c	        write(*,*) count_roughness,'/',count_reflection
c           write(*,*) '==========================='

			aaa = 0
	        bbb = 0

			do nang=0,split 
				OUT  = -(pi/2)+(nang*pi/split)
				Pdif = (cos(INC)+cos(OUT))**2/
     &			       (1+(lambda*lambda*kkk2*kkk2*			!ラフネス散乱レート計算でcpとepの差分を足している
     &				   (sin(INC)+sin(OUT))**2)/2)	
				aaa  = aaa + Pdif
			end do

              rB   = rnd()

			do nang=0,split
				OUT  = -(pi/2)+(nang*pi/split)
				Pdif = (cos(INC)+cos(OUT))**2/
     &				   (1+(lambda*lambda*kkk2*kkk2*			!ラフネス散乱レート計算でcpとepの差分を足している
     &				   (sin(INC)+sin(OUT))**2)/2)
                  bbb  = bbb + Pdif
				ccc  = bbb/aaa
				if(rB.le.ccc) then
					go to 4545
	            end if
			end do

4545			if(akx.ge.0) then			!akxプラスの場合
				if(akz.ge.0) then
c					write(*,*) 'akx入=',akx
c	                write(*,*) 'akz入=',akz
					akx = -kkk*sin(OUT)
					akz = -kkk*cos(OUT)
					kakudo_OUT = OUT*180/pi
c	                write(*,*) 'akx出=',akx
c	                write(*,*) 'akz出=',akz
c                     write(*,*) '散乱前',kakudo_INC,'度'
c					write(*,*) '散乱後',kakudo_OUT,'度'
				else
c                     write(*,*) 'akx入=',akx
c	                write(*,*) 'akz入=',akz
					akx = -kkk*sin(OUT)
					akz =  kkk*cos(OUT)
 					kakudo_OUT = OUT*180/pi
c                     write(*,*) 'akx出=',akx
c	                write(*,*) 'akz出=',akz
c                     write(*,*) '散乱前',kakudo_INC,'度'
c					write(*,*) '散乱後',kakudo_OUT,'度'
				end if
			else
				if(akz.ge.0) then
c                     write(*,*) 'akx入=',akx
c	                write(*,*) 'akz入=',akz
					akx =  kkk*sin(OUT)
					akz = -kkk*cos(OUT)  	 
					kakudo_out = OUT*180/pi
c                     write(*,*) 'akx出=',akx
c	                write(*,*) 'akz出=',akz
c                     write(*,*) '散乱前',kakudo_INC,'度'
c					write(*,*) '散乱後',kakudo_OUT,'度'
				else
c                     write(*,*) 'akx入=',akx
c	                write(*,*) 'akz入=',akz
                      akx =  kkk*sin(OUT)
					akz =  kkk*cos(OUT)  	 
					kakudo_OUT = OUT*180/pi
c                     write(*,*) 'akx出=',akx
c	                write(*,*) 'akz出=',akz
c                     write(*,*) '散乱前',kakudo_INC,'度'
c					write(*,*) '散乱後',kakudo_OUT,'度'
                  end if

			end if 	
     
	        zure = -(OUT+INC)*180/pi
	        ddd  = ddd + zure
c 	        write(*,*) '==========================='
c              write(*,*) '位相差',zure,'度'
	        average = ddd/count_roughness1
c			write(*,*) 'Avereage',average,'度'  
				  
          else
          
			akz = -akz

	    end if

	endif
c
      return
	end