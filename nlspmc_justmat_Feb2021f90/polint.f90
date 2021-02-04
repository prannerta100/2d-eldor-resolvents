    SUBROUTINE polint(xa,ya,n,x,y,dy)
    include 'stdio.inc'
    integer :: n,NMAX
    double precision :: dy,x,y,xa(n),ya(n)
    PARAMETER (NMAX=10)
    integer :: i,m,ns
    double precision :: den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
    ns=1
    dif=abs(x-xa(1))
!      write(luttyo,*) "n,x=",n,x
!      write(luttyo,*) "line 1 xa,ya"
!      write(luttyo,*) "line 2 c,d"
            do 11 i=1,n
    !        write(luttyo,*) xa(i),ya(i)
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
            ns=i
            dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
    !        write(luttyo,*) c(i),d(i)
        11 continue
                        y=ya(ns)
    !      write(luttyo,*) "initial approximation ",y 
        ns=ns-1
        do 13 m=1,n-1
            do 12 i=1,n-m
                ho=xa(i)-x
                hp=xa(i+m)-x
                w=c(i+1)-d(i)
                den=ho-hp
                if(den.eq.0.)pause 'failure in polint'
                den=w/den
            !          write(luttyo,*) "den=",den
                d(i)=hp*den
                c(i)=ho*den
            !          write(luttyo,*) "c,d:",c(i),d(i)
                12 continue
                if (2*ns.lt.n-m)then
                    dy=c(ns+1)
                else
                    dy=d(ns)
                    ns=ns-1
                endif
                y=y+dy
            !        write(luttyo,*) "dy=",dy
                13 continue
                return
                END
