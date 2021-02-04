! NLSPMC Version 1.0 2/5/99
!*********************************************************************

!                         SUBROUTINE SWITCH
!                         =================

!        This subroutine shifts the second half of the FT,
!        making it the first half.

!*********************************************************************

    subroutine switch(temp,npt)

    integer :: i,j,k

    complex*16 temp(npt),tmp

!#####################################################################

    k=npt/2
    do 10 i=1,k
        tmp=temp(i+k)
        temp(i+k)=temp(i)
        temp(i)=tmp
    10 END DO

    return
    end subroutine switch
