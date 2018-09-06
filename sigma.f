	subroutine sigma(x,sum,reg)
c
	real x,sum,reg
c
	reg  = reg + x
	temp = sum
	sum  = sum + reg
	temp = sum - temp
	reg  = reg - temp
	return
	end subroutine
c