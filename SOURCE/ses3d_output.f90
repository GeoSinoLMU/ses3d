!==============================================================================
! store field values at lower order collocation points
!==============================================================================

subroutine ses3d_store(field,name,n,deg,tag)
use parameters
use variables
implicit none
include 'mpif.h'

    ! field - field variable
    ! name  - string containing the name of the field variable
    ! n     - iteration number
    ! deg   - polynomial degree of the stored field

    !======================================================================
    ! local variables
    !======================================================================

    integer, intent(in) :: n, deg, tag
    integer :: i, j, k, ip, iq, ir

    real :: cknots(0:7)
    real :: lgll, lgll_dummy
    real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd), intent(in) :: field
    real, dimension(0:nx,0:ny,0:nz,0:deg,0:deg,0:deg) :: cfield

    character(len=10) :: dummy, ns
    character(len=1000) :: fn
    character(len=*), intent(in) :: name

    !======================================================================
    ! compress only if deg<lpd
    !======================================================================

    if (deg<lpd) then

       !======================================================================
       ! determine collocation points (knots) and integration weights
       !======================================================================

       if (deg==0) then
               
           cknots(0)=0.0

       elseif (deg==1) then
               
           cknots(0)=-1.0
           cknots(1)=1.0               

       elseif (deg==2) then

           cknots(0)=-1.0
           cknots(1)=0.0
           cknots(2)=1.0

       elseif (deg==3) then

           cknots(0)=-1.0
           cknots(1)=-0.4472135954999579
           cknots(2)=0.4472135954999579
           cknots(3)=1.0

       elseif (deg==4) then

           cknots(0)=-1.0
           cknots(1)=-0.6546536707079772
           cknots(2)=0.0
           cknots(3)=0.6546536707079772
           cknots(4)=1.0

       elseif (deg==5) then

           cknots(0)=-1.0
           cknots(1)=-0.7650553239294647
           cknots(2)=-0.2852315164806451
           cknots(3)=0.2852315164806451
           cknots(4)=0.7650553239294647
           cknots(5)=1.0

       elseif (deg==6) then

           cknots(0)=-1.0
           cknots(1)=-0.8302238962785670
           cknots(2)=-0.4688487934707142
           cknots(3)=0.0
           cknots(4)=0.4688487934707142
           cknots(5)=0.8302238962785670
           cknots(6)=1.0

       elseif (deg==7) then

           cknots(0)=-1.0
           cknots(1)=-0.8717401485096066
           cknots(2)=-0.5917001814331423
           cknots(3)=-0.2092992179024789
           cknots(4)=0.2092992179024789
           cknots(5)=0.5917001814331423
           cknots(6)=0.8717401485096066
           cknots(7)=1.0

       endif

       !======================================================================
       ! generate lower degree interpolant
       !======================================================================

       cfield(:,:,:,:,:,:)=0.0

       do i=0,deg
       do j=0,deg
       do k=0,deg

           do ip=0,lpd
           do iq=0,lpd
           do ir=0,lpd

               lgll_dummy=lgll(lpd,ip,cknots(i))*lgll(lpd,iq,cknots(j))*lgll(lpd,ir,cknots(k))
               cfield(0:nx,0:ny,0:nz,i,j,k)=cfield(0:nx,0:ny,0:nz,i,j,k)+field(0:nx,0:ny,0:nz,ip,iq,ir)*lgll_dummy

           enddo
           enddo
           enddo

       enddo
       enddo
       enddo

    endif

    !======================================================================
    ! store field variable
    !======================================================================

    call int2str(my_rank,dummy)

    if (n==1) then
            
        fn=ffd(1:len_trim(ffd))//name(1:len_trim(name))//'_'//dummy(1:len_trim(dummy))

        write(99,*) fn

        open(unit=tag,file=fn,action='write',form='unformatted')

    endif

    if (deg<lpd) then

        write(tag) cfield

    else

        write(tag) field

    endif

    if (n==nt) then

        close(unit=tag)

    endif

end subroutine ses3d_store


!==============================================================================
! read field values at lower order collocation points and interpolate
! to higher order collocation points
!==============================================================================

subroutine ses3d_restore(field,name,idx,deg,tag)
use parameters
use variables
implicit none
include 'mpif.h'

    !======================================================================
    ! local variables
    !======================================================================

    integer, intent(in) :: idx, deg, tag
    integer :: i, j, k, il, im, in

    real :: lgll, lgll_dummy
    real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd), intent(inout) :: field
    real, dimension(0:nx,0:ny,0:nz,0:deg,0:deg,0:deg) :: cfield

    character(len=*), intent(in) :: name
    character(len=10) :: dummy
    character(len=1000) :: fn

    !======================================================================
    ! read lower degree polynomial collocation point values
    !======================================================================

    call int2str(my_rank,dummy)

    !- open file in first iteration

    if (idx==1) then

        fn=ffd(1:len_trim(ffd))//name(1:len_trim(name))//'_'//dummy(1:len_trim(dummy))

        open(unit=tag,file=fn,action='read',form='unformatted',position='append')
        backspace(unit=tag)

    endif

    !- read and check if interpolation is necessary

    if (deg<lpd) then
        
        read(tag) cfield
        backspace(unit=tag)
        if (idx<nt) then
            backspace(unit=tag)
        endif

    else

        read(tag) field
        backspace(unit=tag)
        if (idx<nt) then
            backspace(unit=tag)
        endif

    endif

    !- close file in last iteration

    if (idx==nt) then

        close(unit=tag)

    endif

    !======================================================================
    ! interpolate to higher polynomial collocation point values if necessary
    !======================================================================

    if (deg<lpd) then

        field=0.0

        do i=0,lpd
        do j=0,lpd
        do k=0,lpd

            do il=0,deg
            do im=0,deg
            do in=0,deg

                lgll_dummy=lgll(deg,il,knots(i))*lgll(deg,im,knots(j))*lgll(deg,in,knots(k))
                field(0:nx,0:ny,0:nz,i,j,k)=field(0:nx,0:ny,0:nz,i,j,k)+cfield(0:nx,0:ny,0:nz,il,im,in)*lgll_dummy

            enddo
            enddo
            enddo

        enddo
        enddo
        enddo

    endif

end subroutine ses3d_restore


subroutine ses3d_output(n)
use parameters
use variables
implicit none
include 'mpif.h'

    !======================================================================
    ! local variables
    !======================================================================

    integer, intent(in) :: n
    integer :: i, j, k
    character(len=10) :: dummy, ns
    character(len=1000) :: fn
    character(len=12) :: sname

    !======================================================================
    ! write displacement fields and preconditioner to files
    !======================================================================

    call int2str(my_rank,dummy)
    call int2str(n,ns)

    if (output_displacement==1) then

        fn=ofd(1:len_trim(ofd))//'vx_'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns))
        open(unit=10,file=fn,action='write',form='unformatted')

        fn=ofd(1:len_trim(ofd))//'vy_'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns))
        open(unit=11,file=fn,action='write',form='unformatted')

        fn=ofd(1:len_trim(ofd))//'vz_'//dummy(1:len_trim(dummy))//'_'//ns(1:len_trim(ns))
        open(unit=12,file=fn,action='write',form='unformatted')

        write(10) vx
        write(11) vy
        write(12) vz

        close(unit=10)
        close(unit=11)
        close(unit=12)

    endif

    !======================================================================
    ! write Frechet derivatives
    !======================================================================

    if (adjoint_flag==2) then

        !- Kernels for elastic parameters.=============================

if (USE_DIAG_GRAD == 1) then
        call diag_grad
endif

        fn=ofd(1:len_trim(ofd))//'grad_rho_'//dummy(1:len_trim(dummy))
        open(unit=10,file=fn,status='replace',action='write',form='unformatted')

        fn=ofd(1:len_trim(ofd))//'grad_cp_'//dummy(1:len_trim(dummy))
        open(unit=12,file=fn,status='replace',action='write',form='unformatted')

        fn=ofd(1:len_trim(ofd))//'grad_csh_'//dummy(1:len_trim(dummy))
        open(unit=14,file=fn,status='replace',action='write',form='unformatted')

        fn=ofd(1:len_trim(ofd))//'grad_csv_'//dummy(1:len_trim(dummy))
        open(unit=15,file=fn,status='replace',action='write',form='unformatted')

        write(10) grad_rho
        write(12) grad_cp
        write(14) grad_csh
        write(15) grad_csv

        close(unit=10)
        close(unit=12)
        close(unit=14)
        close(unit=15)

        !- Kernels for visco-elastic parameters. ======================

!       if (is_diss==1) then

!           fn=ofd(1:len_trim(ofd))//'grad_Q_mu_'//dummy(1:len_trim(dummy))
!           open(unit=10,file=fn,status='replace',action='write',form='unformatted')

!           fn=ofd(1:len_trim(ofd))//'grad_alpha_mu_'//dummy(1:len_trim(dummy))
!           open(unit=11,file=fn,status='replace',action='write',form='unformatted')

!           fn=ofd(1:len_trim(ofd))//'grad_Q_kappa_'//dummy(1:len_trim(dummy))
!           open(unit=12,file=fn,status='replace',action='write',form='unformatted')

!           fn=ofd(1:len_trim(ofd))//'grad_alpha_kappa_'//dummy(1:len_trim(dummy))
!           open(unit=13,file=fn,status='replace',action='write',form='unformatted')

!           write(10) grad_Q_mu
!           write(11) grad_alpha_mu
!           write(12) grad_Q_kappa
!           write(13) grad_alpha_kappa

!           close(unit=10)
!           close(unit=11)
!           close(unit=12)
!           close(unit=13)

!       endif


    endif

    !======================================================================
    ! write saving vector to files if forward adjoint simulation
    !======================================================================

    if (adjoint_flag==1) then

        fn=ffd(1:len_trim(ffd))//'saving_vector_'//dummy(1:len_trim(dummy))
        open(unit=10,file=fn,action='write')

        do k=1,nt

            write(10,*) saving_vector(k)

        enddo

        close(unit=10)

    endif

    !======================================================================
    ! write seismograms to files
    !======================================================================

    if ((nr>0) .and. (adjoint_flag<2)) then

        do i=1,nr

            sname=station_name_local(i)

            open(unit=51,file=ofd(1:len_trim(ofd))//sname//'.x',action='write')
            open(unit=52,file=ofd(1:len_trim(ofd))//sname//'.y',action='write')
            open(unit=53,file=ofd(1:len_trim(ofd))//sname//'.z',action='write')

            write(51,*) 'theta component seismograms'
            write(52,*) 'phi component seismograms'
            write(53,*) 'r component seismograms'

            write(51,*) 'nt=', nt
            write(51,*) 'dt=', dt

            write(52,*) 'nt=', nt
            write(52,*) 'dt=', dt

            write(53,*) 'nt=', nt
            write(53,*) 'dt=', dt

            write(51,*) 'receiver location (colat [deg],lon [deg],depth [m])'
            write(52,*) 'receiver location (colat [deg],lon [deg],depth [m])'
            write(53,*) 'receiver location (colat [deg],lon [deg],depth [m])'

            write(51,*) 'x=', recloc(1,i)*180/pi, 'y=', recloc(2,i)*180/pi, 'z=', recloc(3,i)
            write(52,*) 'x=', recloc(1,i)*180/pi, 'y=', recloc(2,i)*180/pi, 'z=', recloc(3,i)
            write(53,*) 'x=', recloc(1,i)*180/pi, 'y=', recloc(2,i)*180/pi, 'z=', recloc(3,i)

            write(51,*) 'source location (colat [deg],lon [deg],depth [m])'
            write(52,*) 'source location (colat [deg],lon [deg],depth [m])'
            write(53,*) 'source location (colat [deg],lon [deg],depth [m])'

            write(51,*) 'x=', xxs*180/pi, 'y=', yys*180/pi, 'z=', zzs
            write(52,*) 'x=', xxs*180/pi, 'y=', yys*180/pi, 'z=', zzs
            write(53,*) 'x=', xxs*180/pi, 'y=', yys*180/pi, 'z=', zzs

            do j=1,nt

                write(51,*) seismogram_x(i,j)
                write(52,*) seismogram_y(i,j)
                write(53,*) seismogram_z(i,j)

            enddo

            close(unit=51)
            close(unit=52)
            close(unit=53)

        enddo

    endif

end subroutine ses3d_output


subroutine record_seismograms
use parameters
use variables
implicit none
include 'mpif.h'

    !======================================================================
    ! local variables
    !======================================================================

    integer :: idx
    integer :: i, j, k
    real :: dummy
    real :: lgll

    !======================================================================
    ! interpolate displacement field
    !======================================================================

    if (nr>0) then

        do idx=1,nr

            do i=0,lpd
            do j=0,lpd
            do k=0,lpd

                dummy=lgll(lpd,i,recloc_std(1,idx))*lgll(lpd,j,recloc_std(2,idx))*lgll(lpd,k,recloc_std(3,idx))

                seismogram_x(idx,it)=seismogram_x(idx,it)+vx(rx(idx),ry(idx),rz(idx),i,j,k)*dummy
                seismogram_y(idx,it)=seismogram_y(idx,it)+vy(rx(idx),ry(idx),rz(idx),i,j,k)*dummy
                seismogram_z(idx,it)=seismogram_z(idx,it)+vz(rx(idx),ry(idx),rz(idx),i,j,k)*dummy

            enddo
            enddo
            enddo

        enddo

    endif

end subroutine record_seismograms

!
!    ---- Diagnostics (optional)
!

subroutine diag_minmax(data, dmin, dmax)
use parameters
use variables
implicit none
include 'mpif.h'

    real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd), intent(in) :: data
    real, intent(out) :: dmin, dmax

    call mpi_reduce(minval(data),dmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
    call mpi_reduce(maxval(data),dmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)

end subroutine diag_minmax

subroutine diag_model
use parameters
use variables
implicit none
include 'mpif.h'

    real :: rho_min, rho_max
    real :: mu_min, mu_max
    real :: lambda_min, lambda_max
    real :: A_min, A_max
    real :: B_min, B_max
    real :: C_min, C_max
    real :: tau_min, tau_max

    call diag_minmax(1.0/rhoinv, rho_min, rho_max)
    call diag_minmax(mu, mu_min, mu_max)
    call diag_minmax(lambda, lambda_min, lambda_max)
    call diag_minmax(A, A_min, A_max)
    call diag_minmax(B, B_min, B_max)
    call diag_minmax(C, C_min, C_max)
    if (is_diss == 1) then
        call diag_minmax(tau, tau_min, tau_max)
    endif

    if (my_rank == 0) then
        write (*, *) 'DIAG: Model'
        write (*, *) 'DIAG:     rho:    min', rho_min, 'max', rho_max
        write (*, *) 'DIAG:     mu:     min', mu_min, 'max', mu_max
        write (*, *) 'DIAG:     lambda: min', lambda_min, 'max', lambda_max
        write (*, *) 'DIAG:     A:      min', A_min, 'max', A_max
        write (*, *) 'DIAG:     B:      min', B_min, 'max', B_max
        write (*, *) 'DIAG:     C:      min', C_min, 'max', C_max
        if (is_diss == 1) then
            write (*, *) 'DIAG:     tau:    min', tau_min, 'max', tau_max
        endif
    endif

end subroutine diag_model

subroutine diag_fw(tag, p_vx, p_vy, p_vz, p_exx, p_eyy, p_ezz, p_exy, p_exz, p_eyz)
use parameters
use variables
implicit none
include 'mpif.h'

    character(len=*), intent(in) :: tag
    real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd), intent(in) :: p_vx, p_vy, p_vz
    real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd), intent(in) :: p_exx, p_eyy, p_ezz
    real, dimension(0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd), intent(in) :: p_exy, p_exz, p_eyz

    real :: vx_min, vx_max
    real :: vy_min, vy_max
    real :: vz_min, vz_max
    real :: exx_min, exx_max
    real :: eyy_min, eyy_max
    real :: ezz_min, ezz_max
    real :: exy_min, exy_max
    real :: exz_min, exz_max
    real :: eyz_min, eyz_max

    call diag_minmax(p_vx, vx_min, vx_max)
    call diag_minmax(p_vy, vy_min, vy_max)
    call diag_minmax(p_vz, vz_min, vz_max)
    call diag_minmax(p_exx, exx_min, exx_max)
    call diag_minmax(p_eyy, eyy_min, eyy_max)
    call diag_minmax(p_ezz, ezz_min, ezz_max)
    call diag_minmax(p_exy, exy_min, exy_max)
    call diag_minmax(p_exz, exz_min, exz_max)
    call diag_minmax(p_eyz, eyz_min, eyz_max)

    if (my_rank == 0) then
        write(*, *) 'DIAG: Iteration', it, ': ', tag
        write(*, *) 'DIAG:     vx_fw:  min', vx_min, 'max', vx_max
        write(*, *) 'DIAG:     vy_fw:  min', vy_min, 'max', vy_max
        write(*, *) 'DIAG:     vz_fw:  min', vz_min, 'max', vz_max
        write(*, *) 'DIAG:     exx_fw: min', exx_min, 'max', exx_max
        write(*, *) 'DIAG:     eyy_fw: min', eyy_min, 'max', eyy_max
        write(*, *) 'DIAG:     ezz_fw: min', ezz_min, 'max', ezz_max
        write(*, *) 'DIAG:     exy_fw: min', exy_min, 'max', exy_max
        write(*, *) 'DIAG:     exz_fw: min', exz_min, 'max', exz_max
        write(*, *) 'DIAG:     eyz_fw: min', eyz_min, 'max', eyz_max
    endif

end subroutine diag_fw

subroutine diag_store_fw
use parameters
use variables
implicit none
include 'mpif.h'

    call diag_fw('store FW', vx, vy, vz, dxux, dyuy, dzuz, (dxuy+dyux)/2, (dxuz+dzux)/2, (dyuz+dzuy)/2)

end subroutine diag_store_fw
      
subroutine diag_restore_fw
use parameters
use variables
implicit none
include 'mpif.h'

    call diag_fw('restore FW', vx_fw, vy_fw, vz_fw, exx_fw, eyy_fw, ezz_fw, exy_fw, exz_fw, eyz_fw)

end subroutine diag_restore_fw

subroutine diag_grad
use parameters
use variables
implicit none
include 'mpif.h'

    real :: grad_rho_min, grad_rho_max
    real :: grad_cp_min, grad_cp_max
    real :: grad_csh_min, grad_csh_max
    real :: grad_csv_min, grad_csv_max

    call diag_minmax(grad_rho, grad_rho_min, grad_rho_max)
    call diag_minmax(grad_cp, grad_cp_min, grad_cp_max)
    call diag_minmax(grad_csh, grad_csh_min, grad_csh_max)
    call diag_minmax(grad_csv, grad_csv_min, grad_csv_max)

    if (my_rank == 0) then
        write (*, *) 'DIAG: Iteration', it, ': gradients'
        write (*, *) 'DIAG:     grad_rho:  min', grad_rho_min, 'max', grad_rho_max
        write (*, *) 'DIAG:     grad_cp:   min', grad_cp_min, 'max', grad_cp_max
        write (*, *) 'DIAG:     grad_csh:  min', grad_csh_min, 'max', grad_csh_max
        write (*, *) 'DIAG:     grad_csv:  min', grad_csv_min, 'max', grad_csv_max
    endif

end subroutine diag_grad
