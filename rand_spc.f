c
c---( IBM pSeries 690 ��p���������֐��̒�` )----
c--- random.f�Ɣr�����p ---
c
	subroutine init_rand()
c
	call random_seed(generator = 2)
c
      return
      end
c
c===( �X�p�R���p���������֐� )===
c
      real function rnd()
c
c
   10 continue
	call random_number(rnd)
      if((rnd.ge.1.0).or.(rnd.le.0.0)) goto 10
c
      return
      end
