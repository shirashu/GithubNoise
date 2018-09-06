c
c---( IBM pSeries 690 専用乱数発生関数の定義 )----
c--- random.fと排他利用 ---
c
	subroutine init_rand()
c
	call random_seed(generator = 2)
c
      return
      end
c
c===( スパコン用乱数発生関数 )===
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
