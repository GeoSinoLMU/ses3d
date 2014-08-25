! THIS IS A VERSION OF SES3D_GRAD.F90 FOR THE COMPUTATION OF KERNELS FOR
! FOR VISCO-ELASTIC PARAMETERS. THIS CODE IS STILL A RESEARCH VERSION.

!*****************************************************************************
! computation of Frechet kernels *********************************************
!*****************************************************************************

subroutine ses3d_grad
use parameters
use variables
implicit none
include 'mpif.h'

    !======================================================================
    !- Local variables
    !======================================================================

    integer :: k
    real, dimension(1:nrdiss,0:nx_max,0:ny_max,0:nz_max,0:lpd,0:lpd,0:lpd) :: M_dummy
    
    !- These are the derivatives of the weights D with respect to alpha. They are hardcoded here, and may need to be changed.

    real, dimension(1:nrdiss) :: dPdalpha
    real :: alpha, Q_kappa

    dPdalpha(1)=-3.06
    dPdalpha(2)=-1.54
    dPdalpha(3)=2.56

    !- Also, alpha and Q_kappa are hardcoded right now. 
    
    alpha=0.3

    Q_kappa=57823.0

    !======================================================================
    !- Compute kernels by integrating on the fly.
    !======================================================================

    if (saving_vector(nt+1-it)>0) then

        !==============================================================
        ! read dynamic fields
        !==============================================================

        !- read velocity fields ---------------------------------------

        call ses3d_restore(vx_fw,'vx',it,4,1021)
        call ses3d_restore(vy_fw,'vy',it,4,1022)
        call ses3d_restore(vz_fw,'vz',it,4,1023)

        !- read strain rates ------------------------------------------

        call ses3d_restore(exx_fw,'exx',it,4,2021)
        call ses3d_restore(eyy_fw,'eyy',it,4,2022)
        call ses3d_restore(ezz_fw,'ezz',it,4,2023)

        call ses3d_restore(exy_fw,'exy',it,4,2024)
        call ses3d_restore(exz_fw,'exz',it,4,2025)
        call ses3d_restore(eyz_fw,'eyz',it,4,2026)

        !==============================================================
        !- Preliminaries.==============================================
        !==============================================================

        !- Adjoint strain tensor.

        exx_ad=dxux
        eyy_ad=dyuy
        ezz_ad=dzuz

        exy_ad=0.5*(dxuy+dyux)
        exz_ad=0.5*(dxuz+dzux)
        eyz_ad=0.5*(dyuz+dzuy)

        !- Traces of forward and adjoint strain tensor.
        
        tr_e_ad=exx_ad+eyy_ad+ezz_ad
        tr_e_fw=exx_fw+eyy_fw+ezz_fw

        !- Scalar product of forward and adjoint strain tensor.
        
        e_ddot_e=exx_ad*exx_fw+eyy_ad*eyy_fw+ezz_ad*ezz_fw+2.0*exy_ad*exy_fw+2.0*exz_ad*exz_fw+2.0*eyz_ad*eyz_fw


        !==============================================================
        !- Frechet kernels for elastic parameters. ====================
        !============================================================== 
        !- These are kernels for *absolute* perturbations. ============
        !============================================================== 

        grad_cp=grad_cp+samp_ad*dt*(2.0*rho*cp*tr_e_fw*tr_e_ad)/Jac

        grad_rho=grad_rho+samp_ad*dt*(vx_fw*vx+vy_fw*vy+vz_fw*vz)/Jac
        grad_rho=grad_rho+samp_ad*dt*((cp*cp-2.0*cs*cs)*tr_e_fw*tr_e_ad+2.0*cs*cs*e_ddot_e)/Jac

        grad_csh=grad_csh+samp_ad*dt*(4.0*rho*cs*(e_ddot_e-tr_e_ad*tr_e_fw-2.0*exz_fw*exz_ad-2.0*eyz_fw*eyz_ad))/Jac
        grad_csv=grad_csv+samp_ad*dt*(8.0*rho*cs*(exz_fw*exz_ad+eyz_fw*eyz_ad))/Jac


        !==============================================================
        !- Frechet kernels for visco-elastic parameters. ==============
        !==============================================================
        !- These are kernels for *relative* perturbations. ============
        !==============================================================

        if (is_diss==1) then

            !- Shear Q ============================================

            do k=1,nrdiss

                M_dummy(k,:,:,:,:,:,:)=Mxx(k,:,:,:,:,:,:)*exx_fw(:,:,:,:,:,:)+Myy(k,:,:,:,:,:,:)*eyy_fw(:,:,:,:,:,:)+Mzz(k,:,:,:,:,:,:)*ezz_fw(:,:,:,:,:,:) &
                        +Mxy(k,:,:,:,:,:,:)*exy_fw(:,:,:,:,:,:)+Mxz(k,:,:,:,:,:,:)*exz_fw(:,:,:,:,:,:)+Myz(k,:,:,:,:,:,:)*eyz_fw(:,:,:,:,:,:) &
                        -tr_e_fw(:,:,:,:,:,:)*(Mxx(k,:,:,:,:,:,:)+Myy(k,:,:,:,:,:,:)+Mzz(k,:,:,:,:,:,:))/3.0

            enddo

            do k=1,nrdiss

                grad_Q_mu(:,:,:,:,:,:)=grad_Q_mu(:,:,:,:,:,:)+samp_ad*dt*M_dummy(k,:,:,:,:,:,:)*tau_p(k)/Jac
                grad_alpha_mu(:,:,:,:,:,:)=grad_alpha_mu(:,:,:,:,:,:)+samp_ad*dt*M_dummy(k,:,:,:,:,:,:)*dPdalpha(k)*tau_p(k)/(D_p(k)*Jac)

            enddo

            !- Bulk Q =============================================

            do k=1,nrdiss

                M_dummy(k,:,:,:,:,:,:)=(Mxx(k,:,:,:,:,:,:)+Myy(k,:,:,:,:,:,:)+Mzz(k,:,:,:,:,:,:))*tr_e_fw

            enddo

            do k=1,nrdiss

                grad_Q_kappa(:,:,:,:,:,:)=grad_Q_kappa(:,:,:,:,:,:)+samp_ad*dt*M_dummy(k,:,:,:,:,:,:)*tau_p(k)/Jac
                grad_alpha_kappa(:,:,:,:,:,:)=grad_alpha_kappa(:,:,:,:,:,:)+samp_ad*dt*M_dummy(k,:,:,:,:,:,:)*dPdalpha(k)*tau_p(k)/(D_p(k)*Jac)

            enddo

            !- Normalisations =====================================================

            if (it==nt) then

                grad_Q_mu=2.0*mu*grad_Q_mu/QQ
                grad_alpha_mu=-2.0*alpha*mu*grad_alpha_mu/QQ
                grad_Q_kappa=kappa*grad_Q_kappa/Q_kappa
                grad_alpha_kappa=-alpha*kappa*grad_alpha_kappa/Q_kappa

            endif

        endif           

    endif

end subroutine ses3d_grad
