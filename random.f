c
c---( ���������֐��̒�` )----
c---	������PC�ƁAVPP5000/3�ŗ��p
c---	rand_spc.f�Ɣr�����p ---
c
	subroutine init_rand()
	integer*4	iseed
	interface
	   subroutine init_genrand(seed)
c
	    !DEC$ ATTRIBUTES C, ALIAS :'_init_genrand' :: init_genrand
		integer*4 seed
	   END SUBROUTINE init_genrand
	end interface
c
c
c---( ���������֐��̏����l )----
	iseed  = 38467
c	k = 0
	call init_genrand(iseed)
c
	return
	end
c
c===( ���������֐� )===
c
      real function rnd()
c
      interface
         double precision function genrand_real3()
c
		!DEC$ ATTRIBUTES C, ALIAS :'_genrand_real3' :: genrand_real3
         END function genrand_real3
      end interface
c
   10 rnd = sngl(genrand_real3())
      if((rnd.ge.1.0).or.(rnd.le.0.0)) goto 10
c
      return
      end

