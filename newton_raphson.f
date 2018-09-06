c-----Newton_Raphson�@�ɂ������������߂�T�u���[�`��-----
	subroutine newton_raphson(ka,ini_fermi,ini_tel,avsumconc,avsumtel,
     &				avsumtei1,efermi,ix,iz,iv,ind,flag3,am,aff)		!���������


	implicit none
c
c---�ϐ��z��p�p�����[�^---
	include 'arraysize.fi'
	common	/arraydata/nx,nz,ntenum
	integer	nx,nz,ntenum
c
	integer	n,m,i,ind,flag3		!���������
	integer	ix,iz,iv
	integer(1)	ka
	real x(20),s(20),ae
	real ini_fermi,ini_tel
	real efermi(0:nx,0:nz,nvalley)
	real avsumtel(0:nx,0:nz,nvalley)
	real avsumconc(0:nx,0:nz,nvalley)
	real(8) avsumtei1(0:nx,0:nz,nvalley)
	real, dimension (nvalley,narea)::am			!120201
	real(8) aff(nvalley,narea)					!120201
c
	real(8) bk,q,bkq
	parameter(bk = 1.38066e-23,q = 1.60219e-19)
c
	n=2
	x(1)=ini_tel		!Tel�̏����l
	x(2)=ini_fermi		!Ef�̏����l
	ae=0.001			!���΋��e�덷 �������������ꂽ���0.001�ɕύX����
	m=100000			!�ő�J��Ԃ���
c
c-----x(2)��ita�ɕϊ�-----
	bkq=bk/q
	x(2)=x(2)/(bkq*x(1))
c-------------------------
c
	call newton_raphson_calc(ka,avsumconc,avsumtei1,
     &						 ix,iz,iv,n,x,ae,m,s,ind,flag3,am,aff)		!���������
c
c-----stop���Ȃ����߂̏���-----
	if(flag3.eq.1)then	!���������
		return
	endif
c------------------------------
c
c-----ita����Ef�����߂�-----
	bkq=bk/q
	s(2)=bkq*s(1)*s(2)
c---------------------------
c
	avsumtel(ix,iz,iv) = s(1)
	efermi(ix,iz,iv)   = s(2)
c
	return
	end