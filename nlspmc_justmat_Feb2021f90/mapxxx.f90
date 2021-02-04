    subroutine mapxxx(xxx,nvar,idir)

! subroutine to map infinite range variables handled by levenburg
! routine into finite range physical parameters.
! This routine is only called from lmnls.f.

!    xxx: 	variables passed by LM routine
!    bvar: 	Number of parameters in xxx
!    idir:	direction to transform, 1=to physical values, else = LM vars.

    integer :: i,nvar,idir
    double precision :: xxx(*)

    include 'limits.inc'
    include 'parms.inc'

    do 10 i=1,nvar
    ! check if this variable has min and max set:
        if(abs(prmax(i))+abs(prmin(i)) > 1D-7) then
            if(idir == 1) then
            ! map to physical values for spectral calculation
                xxx(i)=prmin(i)+(prmax(i)-prmin(i))* &
                (1.D0-1.D0/(1.D0+exp(4.*(xxx(i)- &
                (prmin(i)+prmax(i))/2.D0)/(prmax(i)-prmin(i)))))
            else
            ! map to infinite range for LM routine
                if(xxx(i) >= prmax(i) .OR. xxx(i) <= prmin(i)) then
                ! values should be within limits, force them.
                    write(*,*)'error parameter limits exceeded', &
                    i,xxx(i),prmax(i),prmin(i)
                    if(xxx(i) >= prmax(i)) xxx(i)=prmax(i)* 0.99999D0
                    if(xxx(i) <= prmax(i)) xxx(i)=prmin(i)* 1.00001D0
                end if
            ! ok to do the mapping:
                xxx(i)=(prmin(i)+prmax(i))/2.+(prmax(i)-prmin(i))/4.* &
                log(1./(1.-(xxx(i)-prmin(i))/(prmax(i)-prmin(i)))-1.)
            end if
        end if
    10 END DO
    return
    end subroutine mapxxx
