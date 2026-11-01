program a1

real :: a, b, c 

open (unit=15, file='in_1.txt')
open (unit=16, file='out_1.txt')

read (15,*) a, b, c
b = b*2
print *, a, b, c

write (16,*) a, b, c

stop
end

