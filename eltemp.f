c-----各セル空間の平均濃度と平均エネルギーを計算するサブルーチン-----
	subroutine eltemp(adkx,adky,adkz,
     &				  jpnum,spnum,p,kp,
     &				  af2,af4,hhm,ec,
     &				  dx,dz,xmax,zmax,iarea,lhet,
     &				  tel,mconc,cn,
c     &				  ka3)
     &				  ka3,tel1,ict,avsumtel,avsumtel1,		!竹岸変更
     &				  avsumconc,avsumtei1,sflag,			!竹岸変更
     &				  efermi,efermi1,avsumconc1,avsumtei11,tel2,
     &                  epA,epB,epC,u,eg,						!120126homma
     &                  avsumtei1_1,avsumconc_1,avsumteiA,avsumconcA)	!竹岸変更

	implicit none
c===引数===
c---変数配列用パラメータ---
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c---基本パラメータ---
	real(8) pi,q,h,bk
	parameter(pi = 3.141592654, q  = 1.60219e-19)
	parameter(h  = 1.05459e-34, bk = 1.38066e-23)
	real	spnum
	integer	jpnum
c---デバイス構造---
	real	dx,dz,xmax,zmax
	integer(1)	iarea(nlayer)
	integer(2)  lhet(nlayer) !追加　11/06/22　原
c---領域別パラメータ---
	real,	dimension (nvalley,narea)	:: hhm,af2,af4
	real,	dimension (nvalley,narea)	:: ec
c----(リセス)---
	integer ii
	integer(2),dimension (nrecess) :: lxr1,lzr1,lxr2,lzr2
c---粒子状態---
	real,	dimension   (6,npmax)	:: p
	integer(1),dimension (3,npmax)	:: kp
c===ローカル変数===
	real pdx,pdz
	integer(1)	kv,ken,kl,ka
	real akx,aky,akz
	integer n,eflag
	integer	ix,iz,il,iv
c---kdに関するパラメータ---	
	real adkx(0:nx,0:nz,0:nvalley)
	real adky(0:nx,0:nz,0:nvalley)
	real adkz(0:nx,0:nz,0:nvalley)
	real tkx,tky,tkz
c----(メッシュ&谷別電子濃度)---
	real,	dimension (0:nx,0:nz)	:: cn
	real cnmesh(0:nx,0:nz)
	real kvmesh(0:nx,0:nz,nvalley)
	real kvmesh2(0:nx)  !101221
	real mconc(0:nx,0:nz,nvalley)
	real,	dimension (0:nx,0:nz)	:: epA, epB, epC ,ka4
c
c-----08/8/6 竹岸-----
	real avsumtel(0:nx,0:nz,nvalley)
	real avsumconc(0:nx,0:nz,nvalley)
	real avsumconc_1(0:nx,0:nz,nvalley)  !101221電子濃度チャネル平均
	real avsumconcA(0:nx)                !101221電子濃度チャネル平均2
	real(8) avsumtei1(0:nx,0:nz,nvalley)
      real(8) avsumtei1_1(0:nx,0:nz,nvalley)  !101221電子エネルギーチャネル平均
	real(8) avsumteiA(0:nx)                 !101221電子エネルギーチャネル平均2
	real efermi(0:nx,0:nz,nvalley)
c
	real,save,allocatable :: sumtel
	dimension :: sumtel(:,:,:)
	real,save,allocatable :: sumconc
	dimension :: sumconc(:,:,:)
	real(8),save,allocatable :: sumtei1
	real(8),save,allocatable :: sumtei2  !101221
	dimension :: sumtei1(:,:,:)
	dimension :: sumtei2(:)
	real,	dimension (0:nx,0:nz) :: efermi1		!Γ谷のフェルミレベル
	real,	dimension(0:nx,0:nz) :: avsumtel1		!Γ谷の電子温度
	real,	dimension(0:nx,0:nz) :: avsumconc1		!Γ谷の電子濃度
	real(8),	dimension(0:nx,0:nz) :: avsumtei11	!Γ谷の電子エネルギー
c
	integer	ict,sflag,ict2
	real,	dimension(0:nx,0:nz) :: tel1
	real,	dimension(0:nx,0:nz) :: tel2
c------------------
	real,	dimension (0:nx,0:nz)	:: u	!120126homma
      real	eg(nvalley,narea)				!120126homma
c
	real tel(0:nx,0:nz,nvalley)
	real(8) tei1(0:nx,0:nz,nvalley)
	real(8) tei2(0:nx)  !101221
c
	real(8) ei1,sk,sq
	real bkq23
	integer(1) ka3(0:nx,0:nz,nvalley)
	real(8),save,allocatable :: ei2(:,:,:)

	if (.not. allocated(ei2))then				!最初だけ実行
		allocate(ei2(0:nx,0:nz,nvalley))
		allocate(sumtel(0:nx,0:nz,nvalley))		!08/8/6 竹岸
		allocate(sumconc(0:nx,0:nz,nvalley))	!08/8/6 竹岸
		allocate(sumtei1(0:nx,0:nz,nvalley))	!08/8/6 竹岸
		allocate(sumtei2(0:nx))	!101221
		ei2 = 0.0
		tel = 0.0
		ka3 = 0.0
		mconc = 0.0
c-----08/8/6 竹岸-----
		avsumtel = 0.0
		avsumconc = 0.0
		sumconc = 0.0
		sumtel = 0.0
		sumtei1= 0.0
	    sumtei2= 0.0
		tel1 = 0.0
		avsumtel1 = 0.0
		avsumconc1 = 0.0
		avsumtei11 = 0.0
c--------------------
	endif

	tel2=0.0					!08/10/10 竹岸
	kvmesh=0.0 ; cnmesh=0.0		!必要
	tei1 = 0.0					!必要

	tkx = 0.0 ; tky = 0.0 ; tkz = 0.0
	sk = 0.0

	bkq23 = 2.0/(3.0*bk/q)
	pdx = 1.0/dx
	pdz = 1.0/dz

	do n=1,jpnum
		if(kp(1,n).eq.0)cycle

		ei1 = 0.0
		kv	= kp(1,n)	!谷No.
		kl  = kp(3,n)	!層No.
		ka	= iarea(kl)	!素材No.

		ix = min(nx,max(0,nint(p(5,n)*pdx)))
		iz = min(nz,max(0,nint(p(6,n)*pdz)))
		ka3(ix,iz,kv) = ka 

		akx	= p(1,n)
		aky	= p(2,n)
		akz	= p(3,n)

		tkx = abs(p(1,n)) - abs(adkx(ix,iz,kv))
		tky = abs(p(2,n)) - abs(adky(ix,iz,kv))
		tkz = abs(p(3,n)) - abs(adkz(ix,iz,kv))

		sk = tkx*tkx + tky*tky + tkz*tkz

		if(af4(kv,ka).ne.0.0)then
			sq = sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*sk)
c			ei1=(sq-1.0)/af2(kv,ka)+ec(kv,ka)
			ei1=(sq-1.0)/af2(kv,ka)					!08/11/4 竹岸
c     &			-epA(ix,iz)+u(ix,iz)-0.7*(eg(1,ka)-eg(1,2))		!120126homma
		else
c			ei1=hhm(kv,ka)*sk+ec(kv,ka)
			ei1=hhm(kv,ka)*sk						!08/11/4 竹岸
c     &			-epA(ix,iz)+u(ix,iz)-0.7*(eg(1,ka)-eg(1,2))		!120126homma
		endif

c	------------谷占有率集計-------------
		kvmesh(ix,iz,kv) = kvmesh(ix,iz,kv)+1.0		! メッシュ内谷別粒子数
		cnmesh(ix,iz) = cnmesh(ix,iz)+1.0			! メッシュ内粒子総和

c	-----エネルギーEv(k-kd(r))の総和-----
		tei1(ix,iz,kv)=tei1(ix,iz,kv)+ei1
		ei2(ix,iz,kv)= ei1

c     ------エネルギーチャネル平均-------
      if(kv.eq.1)then
	  if((lhet(nchannel1)+1.le.iz).and.
     &	(iz.le.lhet(nchannel2))) then !channel 11/06/22原
              tei2(ix)=tei2(ix)+ei1
			kvmesh2(ix)=kvmesh2(ix)+1.0
	  endif
	endif

	enddo

	do iv=1,nvalley
		do iz=0,nz
			do ix=0,nx
c	リセス内への処理
				eflag=0
				do ii = 1,nrecess
					if((ix.gt.lxr1(ii)).and.(ix.lt.lxr2(ii))
     &					.and.(iz.ge.lzr1(ii)).and.(iz.lt.lzr2(ii)))then
					tei1(ix,iz,iv)=0.0
					mconc(ix,iz,kv) = 0.0
					eflag=1
					exit
					endif
				enddo

			if(eflag.eq.1)cycle
c	リセス外での処理
			if(kvmesh(ix,iz,iv).eq.0.0)then
				tei1(ix,iz,iv)=0.0
				cycle
			endif

			if(cnmesh(ix,iz).ne.0.0) then

				tei1(ix,iz,iv)=tei1(ix,iz,iv)/kvmesh(ix,iz,iv) !エネルギー粒子数平均

				kvmesh(ix,iz,iv) = kvmesh(ix,iz,iv) / cnmesh(ix,iz)		!占有率

				tel(ix,iz,iv) = bkq23 * tei1(ix,iz,iv)			!谷別電子温度算出(Boltzmann)
				mconc(ix,iz,iv) = cn(ix,iz) * kvmesh(ix,iz,iv)	!谷別電子濃度算出

			endif

c	-----08/8/6 竹岸 濃度,エネルギーのステップ平均化-----
		sumconc(ix,iz,iv) = sumconc(ix,iz,iv) + mconc(ix,iz,iv)
		sumtei1(ix,iz,iv) = sumtei1(ix,iz,iv) + tei1(ix,iz,iv)
c	----------------------------------------------

			enddo
		enddo
	enddo

c----------電子エネルギー平均加算-------
      do ix=0,nx
        tei2(ix)=tei2(ix)/kvmesh2(ix)
	  sumtei2(ix)=sumtei2(ix)+tei2(ix)
      enddo
c
c
c-----08/6/3 竹岸 濃度,エネルギーのステップ平均化-----
c	sflag=0の時にはそのままの値,sflag=1の時にはステップ平均を行い
c	各セル空間における電子濃度と電子温度を更新する
	sflag=0
	ict2=0
	if((mod(abs(ict),sumcnt).eq.0))then
		sflag=1
		do iv=1,nvalley
			do iz=0,nz
				do ix=0,nx
				avsumconc(ix,iz,iv) = sumconc(ix,iz,iv) / sumcnt
				avsumtei1(ix,iz,iv) = sumtei1(ix,iz,iv) / sumcnt
				enddo
			enddo
		enddo

	    avsumconc_1 = avsumconc
          avsumtei1_1 = avsumtei1

		sumconc = 0.0
		sumtei1 = 0.0


          avsumconcA = 0.0   !avsumconcA初期化101215

	   do iz = lhet(nchannel1)+1,lhet(nchannel2) !修正11/06/22原 
	     do ix=0,nx
	          avsumconcA(ix)=avsumconcA(ix)+avsumconc(ix,iz,1)
		 enddo
		     ict2=ict2+1
         enddo  

	    do ix=0,nx
	          avsumconcA(ix)=avsumconcA(ix)/ict2
                avsumteiA(ix) =sumtei2(ix)/sumcnt
		enddo
	          sumtei2 = 0.0
	endif

c-----------------------------------------------------
c
c-----08/8/6 竹岸 出力関係-----
	do iz=0,nz
		do ix=0,nx
			avsumtel1(ix,iz)  = avsumtel(ix,iz,1)	!outputで利用 Γ谷のみ
			efermi1(ix,iz)    = efermi(ix,iz,1)		!outputで利用 Γ谷のみ
			avsumconc1(ix,iz) = avsumconc_1(ix,iz,1)	!outputで利用 Γ谷のみ
			avsumtei11(ix,iz) = avsumtei1_1(ix,iz,1)	!outputで利用 Γ谷のみ
			tel2(ix,iz)		  = tel(ix,iz,2)		!outputで利用 L谷の電子温度
		enddo
	enddo

	return
	end