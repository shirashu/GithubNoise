	subroutine enercalc(sk,ei,hhm,af2,af4)
	implicit none
	real(8)	sk,ei
	real	hhm,af2,af4
c
	real(8)	gk,sq
c
	gk=hhm*sk				!hhm=hh/(2m*q)
	if(af4.ne.0.0)then
		sq=1.0+af4*gk
		if((sq.ge.0.0).and.(sq.le.huge(sq)))then
			sq=sqrt(sq)
		else
			write(99,*)'sqの値が不正です(enercalc)',sq,gk,sk
			write( *,*)'sqの値が不正です(enercalc)',sq,gk,sk
			stop
		endif
		ei=(sq-1.0)/af2		!eV
	else
		ei=gk
	endif
      return
	end