	subroutine setcloud(dx,dz,cloud)
c	電荷雲をガウス分布で設定するサブルーチン
	implicit none
c
	include 'arraysize.fi'
	real dx,dz
	real cloud(-ngx:ngx,-ngz:ngz,ngn)	!二次元電荷雲格納配列
c
c  ---内部変数---
	integer mesh		!割り当てるメッシュ数
	real	eps			!cloud(mesh,0)の中心との比
	real	oval		!楕円係数
	real(8)	a,alpha		!a:配列の合計,alpha:分布の広がり
	integer	ix,iz,jn	!ix,iz:座標,jn:テーブルNo.
	real dx2,dz2		!dx2,dz2:分布の広がり指標
	integer	nngx,nngz	!ガウス分布テーブル幅
c
	oval = dx/dz/5	!/10		!楕円係数
	dx2 = dx*dx/oval
	dz2 = dz*dz*oval
	eps = 0.001			!cloud(mesh,0)の中心との比
c
	cloud = 0.0			!cloud:全体配列
	do jn=1,ngn
		mesh = 2**(jn-1)
		alpha = (dx*mesh)*(dx*mesh)/log(eps)
		a = 0.0				!=Σ(cloud)で初期化
		nngz=min(mesh,ngz)
		nngx=min(mesh,ngx)
c
c
		do iz=-nngz,nngz
		do ix=-nngx,nngx
			cloud(ix,iz,jn)=exp(((dx2*ix*ix)+(dz2*iz*iz))/alpha)
			a=a+cloud(ix,iz,jn)
		enddo
		enddo
c
		do iz=-nngz,nngz
		do ix=-nngx,nngx
			cloud(ix,iz,jn)=cloud(ix,iz,jn)/a	!標準化　Σ(cloud)=1.0に
		enddo
		enddo
c
	enddo
c
	end