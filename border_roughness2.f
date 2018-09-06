c
c-----���q�����E�ʂŔ��˂܂��͐i������C�x���g-----
	subroutine border_roughness2(
     &					hhm,af,af2,af4,eg,
     &					akx,aky,akz,kv,kl,kl2,ka,iarea,lhet,iiz,iix, !2006/12/09 Hara
     &                    xxx,vvv,basho_reflection,
     &                      basho_roughness,
     &                    split,delta,lambda,count_reflection,
     &                    count_roughness,average,
     &					epA,u)			!���t�l�X�U���pep���	
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
c	---	���[�J���ϐ�	---
c	skx,sky,skz:�e�����g���̓����i�[����ϐ�
c	vb:��Ǎ������i�[����ϐ�
c	ex,ey,ez:�e�����^���G�l���M�[���i�[����ϐ�
	real(8)	skx,sky,skz
	real	vb
	real(8)	ex,ey,ez
	integer(1)	ka2
c----���t�l�X�U��----------------------------------------------
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

c----ep-cp��Ǘp---(���t�l�X�U��)
	real	vbcp	!cp��Ǎ���
	real,	dimension (0:nx,0:nz)	:: epA
	real,	dimension (0:nx,0:nz)	:: u
	real	ek_epcp	!CP��EP�̍������^���G�l���M�[�Ƃ��đ���(���t�l�X�U�����[�g)
	real	sk_epcp	!�^���G�l���M�[���狁�߂�g��k�̓��
	real	ak_epcp	!�^���G�l���M�[���狁�߂�g��k	 (x,y,z�ǂ̕����ł�������)

c
c=================================================================================================
c
c	--- z�����̉^���G�l���M�[�v�Z ---
	if(af4(kv,ka).ne.0.0)then
		skz = dble(akz)*dble(akz)
		ez=(sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*skz)-1.0)/af2(kv,ka)		!eV
	else
		ez=hhm(kv,ka)*akz*akz
	endif
c	---------------------------------
c
	ka2 = iarea(kl2)
c	vb =0.7*(eg(kv,ka2)-eg(kv,ka))		!Vb:��Ǎ���[eV]
c	if(iz.lt.(nlayer-2)-5) then	!�h���t�g�O�ʒu��classical�̈�̂Ƃ�(2006/08/08����)!2006/12/09 Hara
c	if(iz.lt.(nlayer-4)-4) then	!�h���t�g�O�ʒu��classical�̈�̂Ƃ�(2006/08/08����)!2006/12/22 Hara
c		vb =0.7*(eg(kv,ka2)-eg(kv,ka))		!Vb:��Ǎ���[eV]
c	else						!�h���t�g�O�ʒu��EP�̈�̂Ƃ�
c		vb = 0.0			!�Ȃ� 0 eV ?? 2011/03/23��
c	endif
c
c----(cp-ep)�̃w�e�����------------
	vbcp =0.7*(eg(kv,ka2)-eg(kv,ka))		!Vb:��Ǎ���[eV]	
	open(unit=1212,file='zztriela.txt')	!�m�F�p!	

	if(iiz.eq.lhet(nchannel1))then
		vb = -(u(iix,iiz) - epA(iix,iiz)) !cp-ep�̏��
	endif
	if(iiz.eq.lhet(nchannel2))then
		vb = -(u(iix,iiz+1) - epA(iix,iiz)) !cp-ep�̏��
	endif
	if(vb.lt.0)then	!�d�ɕt�߂�CP,EP�͕s����@�G���[���sato �v��Ȃ�����
		vb=0
	endif


c----(�w�e����ǂ̃��t�l�X�U���ɑ���ep-cp�����̉^���G�l���M�[)
	if(iiz.eq.lhet(nchannel1))then
		ek_epcp = -(epA(iix,iiz)-(u(iix,iiz+1)-vbcp))
	endif

	if(iiz.eq.lhet(nchannel2))then
		ek_epcp = -(epA(iix,iiz)-(u(iix,iiz)-vbcp)) 
	endif

	if(ek_epcp.lt.0)then	!�d�ɕt�߂�CP,EP�͕s����@�G���[���sato
		ek_epcp=0
		vb=0
	endif
	
	if(af4(kv,ka).ne.0.0)then
		sk_epcp = ((af2(kv,ka)*ek_epcp+1)**2-1)/(af4(kv,ka)*hhm(kv,ka))
	else
		sk_epcp = ek_epcp/hhm(kv,ka)		
	endif
		
		ak_epcp = sqrt(sk_epcp/3)		!���t�l�X�U�����[�g�ɑ����g��


	if(ez.gt.vb)then					!���Vb��藱�q�̃G�l���M�[�� ������
c	===	��ǂ����z����ꍇ�i�i���j	===
c		ez = ez-vb			!z�����̃G�l���M�[�����ǃG�l���M�[���Ђ�
		skx = dble(akx)*dble(akx)
		sky = dble(aky)*dble(aky)
		if(af4(kv,ka).ne.0.0)then
			ex=(sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*skx)-1.0)/af2(kv,ka)		!eV
			ey=(sqrt(1.0+af4(kv,ka)*hhm(kv,ka)*sky)-1.0)/af2(kv,ka)		!eV
		else
			ex=hhm(kv,ka)*skx		!�����O�̎��͒P���v�Z
			ey=hhm(kv,ka)*sky
		endif
c
		kl = kl2
		ka = ka2
c		---	���߂��G�l���M�[����g������i�����O�ł��v�Z�\�j
		akx=sign(sngl(sqrt(ex*(1.0+af(kv,ka2)*ex)/hhm(kv,ka2))),akx)
		aky=sign(sngl(sqrt(ey*(1.0+af(kv,ka2)*ey)/hhm(kv,ka2))),aky)
		akz=sign(sngl(sqrt(ez*(1.0+af(kv,ka2)*ez)/hhm(kv,ka2))),akz)
c		---	sign(a,b) = a�̐�Βl��b�̕������|�����l
c
	else
c	===	��ǂ����z���Ȃ��ꍇ�i���ˁj	===
          delta  = 6.5e-10		!�ʉ�����	�i�q�萔�̐����{
		lambda = 200.0e-10		!�ʉ��L����
		split  = 100000		!���ˊp������(akx,akz�̐��x100���ɍ��킹�Ă���)�����ԒZ�k�̂���10��

		INC  = atan(abs(akx)/abs(akz))
		kkk  = sqrt((akx*akx)+(akz*akz))
		kkk2  = sqrt((akx+ak_epcp)**2+(akz+ak_epcp)**2)		!���t�l�X�U�����[�g�v�Z��cp��ep�̍����𑫂��Ă���
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
	
c---------���ˉ񐔃J�E���g-----------------------------
	
	    vvv = iix

		if(iiz.eq.lhet(nchannel1)) then
			if(akz.le.0)then	!�`���l�����w�e����ǂ̏ꍇ�̂݃J�E���g
		    count_reflection = count_reflection + 1
			basho_reflection(vvv,1) = basho_reflection(vvv,1) + 1
			endif
		endif
		if(iiz.eq.lhet(nchannel2)) then
			if(akz.ge.0)then	!�`���l�����w�e����ǂ̏ꍇ�̂݃J�E���g
			count_reflection = count_reflection + 1
			basho_reflection(vvv,2) = basho_reflection(vvv,2) + 1
			endif
		endif


c---------���t�l�X�U���񐔃J�E���g-----------------------------
		if(rA.gt.Pspe) then
	        xxx = iix

		if(iiz.eq.lhet(nchannel1)) then
			if(akz.le.0)then	!�`���l�����w�e����ǂ̏ꍇ�̂݃J�E���g
			count_roughness = count_roughness + 1
			basho_roughness(xxx,1) = basho_roughness(xxx,1) + 1
			endif
		endif	
		if(iiz.eq.lhet(nchannel2)) then
			if(akz.ge.0)then	!�`���l�����w�e����ǂ̏ꍇ�̂݃J�E���g	
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
     &			       (1+(lambda*lambda*kkk2*kkk2*			!���t�l�X�U�����[�g�v�Z��cp��ep�̍����𑫂��Ă���
     &				   (sin(INC)+sin(OUT))**2)/2)	
				aaa  = aaa + Pdif
			end do

              rB   = rnd()

			do nang=0,split
				OUT  = -(pi/2)+(nang*pi/split)
				Pdif = (cos(INC)+cos(OUT))**2/
     &				   (1+(lambda*lambda*kkk2*kkk2*			!���t�l�X�U�����[�g�v�Z��cp��ep�̍����𑫂��Ă���
     &				   (sin(INC)+sin(OUT))**2)/2)
                  bbb  = bbb + Pdif
				ccc  = bbb/aaa
				if(rB.le.ccc) then
					go to 4545
	            end if
			end do

4545			if(akx.ge.0) then			!akx�v���X�̏ꍇ
				if(akz.ge.0) then
c					write(*,*) 'akx��=',akx
c	                write(*,*) 'akz��=',akz
					akx = -kkk*sin(OUT)
					akz = -kkk*cos(OUT)
					kakudo_OUT = OUT*180/pi
c	                write(*,*) 'akx�o=',akx
c	                write(*,*) 'akz�o=',akz
c                     write(*,*) '�U���O',kakudo_INC,'�x'
c					write(*,*) '�U����',kakudo_OUT,'�x'
				else
c                     write(*,*) 'akx��=',akx
c	                write(*,*) 'akz��=',akz
					akx = -kkk*sin(OUT)
					akz =  kkk*cos(OUT)
 					kakudo_OUT = OUT*180/pi
c                     write(*,*) 'akx�o=',akx
c	                write(*,*) 'akz�o=',akz
c                     write(*,*) '�U���O',kakudo_INC,'�x'
c					write(*,*) '�U����',kakudo_OUT,'�x'
				end if
			else
				if(akz.ge.0) then
c                     write(*,*) 'akx��=',akx
c	                write(*,*) 'akz��=',akz
					akx =  kkk*sin(OUT)
					akz = -kkk*cos(OUT)  	 
					kakudo_out = OUT*180/pi
c                     write(*,*) 'akx�o=',akx
c	                write(*,*) 'akz�o=',akz
c                     write(*,*) '�U���O',kakudo_INC,'�x'
c					write(*,*) '�U����',kakudo_OUT,'�x'
				else
c                     write(*,*) 'akx��=',akx
c	                write(*,*) 'akz��=',akz
                      akx =  kkk*sin(OUT)
					akz =  kkk*cos(OUT)  	 
					kakudo_OUT = OUT*180/pi
c                     write(*,*) 'akx�o=',akx
c	                write(*,*) 'akz�o=',akz
c                     write(*,*) '�U���O',kakudo_INC,'�x'
c					write(*,*) '�U����',kakudo_OUT,'�x'
                  end if

			end if 	
     
	        zure = -(OUT+INC)*180/pi
	        ddd  = ddd + zure
c 	        write(*,*) '==========================='
c              write(*,*) '�ʑ���',zure,'�x'
	        average = ddd/count_roughness1
c			write(*,*) 'Avereage',average,'�x'  
				  
          else
          
			akz = -akz

	    end if

	endif
c
      return
	end