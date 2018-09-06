c					call charge_heat(
c     &						 x,z,pdx,pdz,hef,hef_mesh,hef_scat,
c     &						 jp,iscat)
	subroutine charge_heat(x, z, pdx, pdz, hef, hef_mesh, hef_scat,
     &				 jp,iscat)
c
c---( 引数x,zの位置にhefの熱が発生したとして、これをセル電荷雲法を用いて、
c	周囲のメッシュに割り振るサブルーチン )----
c	
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c
c
	real x,z,pdx,pdz
	real hef
	real,	dimension (0:nx,0:nz) :: hef_mesh
	real,	dimension (0:nx,nscat,nvalley) :: hef_scat
	integer jp,iscat
c
	integer ix
	integer	iz
	real xx,zz
	real x2,z2
c
	xx=x*pdx		!粒子の位置をメッシュで標準化
	zz=z*pdz		!粒子の位置をメッシュで標準化
	ix=max(0,min(ifix(xx),nx-1))	!粒子のいる左側のメッシュ中心
	iz=max(0,min(ifix(zz),nz-1))	!粒子のいる上側のメッシュ中心
	x2=1.0-(xx-float(ix))       	!メッシュ中心からのズレ
	z2=1.0-(zz-float(iz))       	!メッシュ中心からのズレ
	hef_mesh(ix  ,iz  )=
     &	hef_mesh(ix  ,iz  )+(     x2 *     z2 )*hef
     	hef_mesh(ix  ,iz+1)=
     &	hef_mesh(ix  ,iz+1)+(     x2 *(1.0-z2))*hef
	hef_mesh(ix+1,iz  )=
     &	hef_mesh(ix+1,iz  )+((1.0-x2)*     z2 )*hef
	hef_mesh(ix+1,iz+1)=
     &	hef_mesh(ix+1,iz+1)+((1.0-x2)*(1.0-z2))*hef

	hef_scat(ix  ,iscat,jp)=
     &	hef_scat(ix  ,iscat,jp)+(    z2)*hef
     	hef_scat(ix+1,iscat,jp)=
     &	hef_scat(ix+1,iscat,jp)+(1.0-z2)*hef

	return
	end
