!-------------------------------------------------------------------------------
!Oscar Barragan, Jun 2019
!Last modified July 2024
!-------------------------------------------------------------------------------
!This routine computes the gamma terms for the Quasi-Periodic kernel
!To be used to compute the covariance matrix of a multidimensional GP regression
subroutine get_QP_gammas(x1,x2,pars,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1,nx2,npars)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2, npars
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  real(kind=mireal), intent(in) ::  pars(0:npars-1)
  real(kind=mireal), intent(out), dimension(0:nx1-1,0:nx2-1) ::  gamma_g_g, gamma_dg_dg, gamma_dg_g, gamma_g_dg
  !
  real(kind=mireal), dimension(0:nx1-1,0:nx2-1) :: titj, tau, sintau, alpha, beta
  real(kind=mireal) :: le, lp, P, le2, lp2, P2

  le = pars(0)
  lp = pars(1)
  P  = pars(2)
!  le = le*P

  le2 = le*le
  lp2 = lp*lp
  P2  = P*P

  call fcdist(x1,x2,titj,nx1,nx2)

  tau = 2.*pi*titj/P

  sintau = sin(tau)

  alpha = pi*sintau/(2.*P*lp2)

  beta = titj/le2

  gamma_g_g = - (sin(tau/2.))**2/2./lp2 - titj*titj/2./le2
  gamma_g_g = exp(gamma_g_g)

  gamma_g_dg = - gamma_g_g * (alpha + beta)

  gamma_dg_g = - gamma_g_dg

  gamma_dg_dg = - alpha*alpha &
                + pi*pi*cos(tau)/P2/lp2         &
                - tau*sintau/(2.*lp2*le2)     &
                - beta*beta                     &
                + 1./le2
  gamma_dg_dg = gamma_g_g * gamma_dg_dg

end subroutine get_QP_gammas

!This routine computes the gamma terms for the Periodic kernel
!To be used to compute the covariance matrix of a multidimensional GP regression
subroutine get_PK_gammas(x1,x2,pars,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1,nx2,npars)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2, npars
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  real(kind=mireal), intent(in) ::  pars(0:npars-1)
  real(kind=mireal), intent(out), dimension(0:nx1-1,0:nx2-1) ::  gamma_g_g, gamma_dg_dg, gamma_dg_g, gamma_g_dg
  !
  real(kind=mireal), dimension(0:nx1-1,0:nx2-1) :: titj, tau, sintau, alpha, beta
  real(kind=mireal) :: lp, P, lp2, P2

  lp = pars(0)
  P  = pars(1)

  lp2 = lp*lp
  P2  = P*P

  call fcdist(x1,x2,titj,nx1,nx2)

  tau = 2.*pi*titj/P

  sintau = sin(tau)

  alpha = pi*sintau/(2.*P*lp2)

  gamma_g_g = - (sin(tau/2.))**2/2./lp2
  gamma_g_g = exp(gamma_g_g)

  gamma_g_dg = - gamma_g_g * alpha

  gamma_dg_g = - gamma_g_dg

  gamma_dg_dg = - alpha*alpha &
                + pi*pi*cos(tau)/P2/lp2

  gamma_dg_dg = gamma_g_g * gamma_dg_dg

end subroutine get_PK_gammas

!Super QP kernel with two sinosoidal terms
subroutine get_SQP_gammas(x1,x2,pars,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1,nx2,npars)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2, npars
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  real(kind=mireal), intent(in) ::  pars(0:npars-1)
  real(kind=mireal), intent(out), dimension(0:nx1-1,0:nx2-1) ::  gamma_g_g, gamma_dg_dg, gamma_dg_g, gamma_g_dg
  !
  real(kind=mireal), dimension(0:nx1-1,0:nx2-1) :: titj, tau1, tau2, sintau1, sintau2
  real(kind=mireal) :: le, lp1, lp2, P1, P2

  le  = pars(0)
  lp1 = pars(1)
  P1  = pars(2)
  lp2 = pars(3)
  P2  = pars(4)

  call fcdist(x1,x2,titj,nx1,nx2)

  tau1 = 2.*pi*titj/P1
  tau2 = 2.*pi*titj/P2

  sintau1 = sin(tau1)
  sintau2 = sin(tau2)


  gamma_g_g = - (sin(tau1/2.))**2/2./lp1/lp1 - (sin(tau2/2.))**2/2./lp2/lp2 - titj*titj/2./le/le
  gamma_g_g = exp(gamma_g_g)

  gamma_g_dg = - gamma_g_g * (pi*sintau1/(2.*P1*lp1*lp1) + pi*sintau2/(2.*P2*lp2*lp2) + titj/le/le)

  gamma_dg_g = - gamma_g_dg

  gamma_dg_dg = 1./le/le + pi*pi*cos(tau1)/(P1*P1*lp1*lp1) + pi*pi*cos(tau2)/(P2*P2*lp2*lp2)
  gamma_dg_dg = gamma_dg_dg - (titj/le/le + pi*sintau1/(2.*lp1*lp1*P1) + pi*sintau2/(2.*lp2*lp2*P2))**2
  gamma_dg_dg = gamma_g_g * gamma_dg_dg

end subroutine get_SQP_gammas

subroutine get_EXP_gammas(x1,x2,pars,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1,nx2,npars)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2, npars
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  real(kind=mireal), intent(in) ::  pars(0:npars-1)
  real(kind=mireal), intent(out), dimension(0:nx1-1,0:nx2-1) ::  gamma_g_g, gamma_dg_dg, gamma_dg_g, gamma_g_dg
  !
  real(kind=mireal), dimension(0:nx1-1,0:nx2-1) :: titj
  real(kind=mireal) :: l !Lambda parameter for the square exponential kernel

  l = pars(0)

  call fcdist(x1,x2,titj,nx1,nx2)

  gamma_g_g = - 0.5*titj*titj/l/l
  gamma_g_g = exp(gamma_g_g)

  gamma_g_dg = titj/l/l
  gamma_g_dg = gamma_g_g * gamma_g_dg

  gamma_dg_g = - gamma_g_dg

  gamma_dg_dg = ( 1./l/l - titj*titj/l/l/l/l)
  gamma_dg_dg = gamma_g_g * gamma_dg_dg

end subroutine get_EXP_gammas

subroutine get_M52_gammas(x1,x2,pars,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1,nx2,npars)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2, npars
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  real(kind=mireal), intent(in) ::  pars(0:npars-1)
  real(kind=mireal), intent(out), dimension(0:nx1-1,0:nx2-1) ::  gamma_g_g, gamma_dg_dg, gamma_dg_g, gamma_g_dg
  !
  real(kind=mireal), dimension(0:nx1-1,0:nx2-1) :: titj, expt, dt, sgn
  real(kind=mireal) :: sq5
  real(kind=mireal) :: l !Lambda parameter for the square exponential kernel


  l = pars(0)

  call fcdist(x1,x2,titj,nx1,nx2)

  sq5  = sqrt(5.)
  dt   = sq5*abs(titj)/l
  expt = exp(-dt)

  !Find the signs required for the absolute value derivative
  sgn = titj/abs(titj)
  !Fill with zero all the divisions by zero
  where (sgn .ne. sgn)
    sgn = 0.0
  end where

  gamma_g_g = expt * ( 1. + dt + dt*dt/3. )

  gamma_g_dg = 1. / 3. * sgn * expt * dt * (dt + 1.)
  gamma_g_dg = sq5 / l * gamma_g_dg

  gamma_dg_g = - gamma_g_dg

  gamma_dg_dg = 1. / 3. * (dt*dt - dt - 1.) * expt
  gamma_dg_dg =  - 5. / l / l * gamma_dg_dg

end subroutine get_M52_gammas

!This subroutine allow us to select what gammas do we want
subroutine get_gammas(x1,x2,pars,kernel,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1,nx2,npars)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2, npars
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  real(kind=mireal), intent(in) ::  pars(0:npars-1)
  character(len=2), intent(in)  :: kernel
  real(kind=mireal), intent(out), dimension(0:nx1-1,0:nx2-1) ::  gamma_g_g, gamma_dg_dg, gamma_dg_g, gamma_g_dg
  !

  if (kernel == 'MQ') then !Multi-Quasi-periodic Kernel
       call get_QP_gammas(x1,x2,pars,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1,nx2,3)
  else if (kernel == 'SQ') then !Super QP Kernel
       call get_SQP_gammas(x1,x2,pars,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1,nx2,5)
  else if (kernel == 'ME') then !Multi-Exponential Kernel
      call get_EXP_gammas(x1,x2,pars,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1,nx2,1)
  else if (kernel == 'MM') then !Multi-Matern 5/2 Kernel
      call get_M52_gammas(x1,x2,pars,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1,nx2,1)
  end if

end subroutine

subroutine MultidimCov(pars,x1,x2,kernel,ndim,cov,nx1,nx2)
  use constants
  implicit none

  !
  integer, intent(in) :: nx1, nx2, ndim
  real(kind=mireal), intent(in) :: pars(0:ndim*2+5-1)
  real(kind=mireal), intent(in) :: x1(0:nx1-1)
  real(kind=mireal), intent(in) :: x2(0:nx2-1)
  character(len=2), intent(in)  :: kernel
  real(kind=mireal), intent(out) :: cov(0:nx1-1,0:nx2-1)
  !

  call Multin(pars,x1,x2,kernel(1:2),ndim,cov,nx1,nx2)

end subroutine

subroutine MultidimCov_rectangular(pars,x1,x2,dl1,dl2,kernel,ndim,cov,nx1,nx2)
  use constants
  implicit none

  !
  integer, intent(in) :: nx1, nx2, ndim
  real(kind=mireal), intent(in) :: pars(0:ndim*2+5-1)
  real(kind=mireal), intent(in) :: x1(0:nx1-1)
  real(kind=mireal), intent(in) :: x2(0:nx2-1)
  integer, intent(in) :: dl1(0:nx1-1),dl2(0:nx2-1)
  character(len=2), intent(in)  :: kernel
  real(kind=mireal), intent(out) :: cov(0:nx1-1,0:nx2-1)
  !

  call Multin_rectangular(pars,x1,x2,dl1,dl2,kernel(1:2),ndim,cov,nx1,nx2)

end subroutine

subroutine Multin(pars,x1,x2,kernel,m,cov,nx1,nx2)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2, m
  real(kind=mireal), intent(in) :: pars(0:m*2+5-1) !m*2 amplitudes, 3 parameters for the QP kernel
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  character(len=2), intent(in)  :: kernel
  real(kind=mireal), intent(out) :: cov(0:nx1-1,0:nx2-1)
 !
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: gamma_g_g, gamma_dg_dg, gamma_g_dg, gamma_dg_g
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1,0:m-1,0:m-1) :: kas
  real(kind=mireal) :: Amplitudes(0:m*2-1)
  real(kind=mireal) :: hyperpars(0:4)
  integer :: i, j, npars

  !Read the parameters
  Amplitudes(:) = pars(0:m*2-1)
  npars = size(pars) - m*2
  hyperpars = pars(m*2:)

  call get_gammas(x1(0:nx1/m-1),x2(0:nx2/m-1),hyperpars,kernel,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1/m,nx2/m,npars)

  if ( nx1 .ne. nx2  ) then !Compute the K's for not squared matrices

  do i = 0, m - 1
    do j = 0, m - 1
      kas(:,:,i,j) = amplitudes(i*2)*amplitudes(j*2)*gamma_g_g + &
                     amplitudes(i*2+1)*amplitudes(j*2+1)*gamma_dg_dg + &
                     ( amplitudes(i*2)*amplitudes(j*2+1) - amplitudes(i*2+1)*amplitudes(j*2) )* gamma_g_dg
      cov(i*nx1/m:(i+1)*nx1/m-1,j*nx2/m:(j+1)*nx2/m-1) = kas(:,:,i,j)
    end do
  end do

  else !Compute the K's for square matrices

  do i = 0, m - 1
    do j = i, m - 1
      kas(:,:,i,j) = amplitudes(i*2)*amplitudes(j*2)*gamma_g_g + &
                     amplitudes(i*2+1)*amplitudes(j*2+1)*gamma_dg_dg + &
                     ( amplitudes(i*2)*amplitudes(j*2+1) - amplitudes(i*2+1)*amplitudes(j*2) )* gamma_g_dg
      cov(i*nx1/m:(i+1)*nx1/m-1,j*nx2/m:(j+1)*nx2/m-1) = kas(:,:,i,j)
    end do
  end do

  do i = 1, m - 1
    do j = 0, i - 1
      kas(:,:,i,j) = transpose(kas(:,:,j,i))
      cov(i*nx1/m:(i+1)*nx1/m-1,j*nx2/m:(j+1)*nx2/m-1) = kas(:,:,i,j)
    end do
  end do

  end if

end subroutine


subroutine Multin_rectangular_backup(pars,x1,x2,dl1,dl2,kernel,m,cov,nx1,nx2)
  ! Parameters:
  ! pars - kernel hyper-parameters
  ! x - vector including all the time stamps of the big matrix (dim nx)
  ! dl - dim_label - labels for each element of the x vector indicating to which dimension belongs (dim nx)
  ! kernel - kernel label
  ! m - number of dimensions to deal with
  ! cov - Big matrix that will be used for the multi-GP regression
  ! nx - total number of observations

  use constants
  implicit none

  integer, intent(in) :: nx1, nx2, m
  real(kind=mireal), intent(in) :: pars(0:m*2+5-1) ! m*2 amplitudes, 3 parameters for the QP kernel
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  integer, intent(in) :: dl1(0:nx1-1),dl2(0:nx2-1)
  character(len=2), intent(in)  :: kernel
  real(kind=mireal), intent(out) :: cov(0:nx1-1,0:nx2-1)

  real(kind=mireal), dimension(0:nx1-1,0:nx2-1) :: gamma_g_g, gamma_dg_dg, gamma_g_dg, gamma_dg_g
  real(kind=mireal) :: Amplitudes(0:m*2-1)
  real(kind=mireal) :: hyperpars(0:4)
  integer :: i, j, npars

  ! Read the parameters
  Amplitudes(:) = pars(0:m*2-1)
  npars = size(pars) - m*2
  hyperpars = pars(m*2:)

  ! Compute all the gammas for the whole K big matrix
  call get_gammas(x1, x2, hyperpars, kernel, gamma_g_g, gamma_g_dg, gamma_dg_g, gamma_dg_dg, nx1, nx2, npars)


 ! Compute the lower triangle and diagonal of the covariance matrix
  do i = 0, nx1 - 1
    do j = 0, i
      cov(i, j) = Amplitudes(dl1(i) * 2) * Amplitudes(dl2(j) * 2) * gamma_g_g(i, j) + &
      Amplitudes(dl1(i) * 2 + 1) * Amplitudes(dl2(j) * 2 + 1) * gamma_dg_dg(i, j) + &
      (Amplitudes(dl1(i) * 2) * Amplitudes(dl2(j) * 2 + 1) - Amplitudes(dl1(i) * 2 + 1) * Amplitudes(dl2(j) * 2)) * gamma_g_dg(i, j)
      if (i /= j) then
        cov(j, i) = cov(i, j)  ! Mirror the lower triangle to the upper triangle
      end if
    end do
  end do


end subroutine

subroutine Multin_rectangular(pars, x1, x2, dl1, dl2, kernel, m, cov, nx1, nx2)
  ! Parameters:
  ! pars - kernel hyper-parameters
  ! x - vector including all the time stamps of the big matrix (dim nx)
  ! dl - dim_label - labels for each element of the x vector indicating to which dimension belongs (dim nx)
  ! kernel - kernel label
  ! m - number of dimensions to deal with
  ! cov - Big matrix that will be used for the multi-GP regression
  ! nx - total number of observations

  use constants
  implicit none

  integer, intent(in) :: nx1, nx2, m
  real(kind=mireal), intent(in) :: pars(0:m*2+5-1) ! m*2 amplitudes, 3 parameters for the QP kernel
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  integer, intent(in) :: dl1(0:nx1-1), dl2(0:nx2-1)
  character(len=2), intent(in) :: kernel
  real(kind=mireal), intent(out) :: cov(0:nx1-1, 0:nx2-1)

  real(kind=mireal), dimension(0:nx1-1, 0:nx2-1) :: gamma_g_g, gamma_dg_dg, gamma_g_dg, gamma_dg_g
  real(kind=mireal) :: Amplitudes(0:m*2-1)
  real(kind=mireal) :: hyperpars(0:4)
  integer :: i, j, npars
  integer, dimension(:), allocatable :: is, js
  integer :: idx, jdx

  ! Read the parameters
  Amplitudes(:) = pars(0:m*2-1)
  npars = size(pars) - m*2
  hyperpars = pars(m*2:)

  ! Compute all the gammas for the whole K big matrix
  call get_gammas(x1, x2, hyperpars, kernel, gamma_g_g, gamma_g_dg, gamma_dg_g, gamma_dg_dg, nx1, nx2, npars)

  ! Initialize the covariance matrix to zero
  cov = 0.0

  if ( nx1 == nx2) then

  ! Loop over each unique dimension pair (i, j) and apply where statements
  do i = 0, m - 1
    do j = 0, i
      ! Integer indices for current dimensions i and j
      is = pack([(idx, idx=0, nx1-1)], dl1 == i)
      js = pack([(jdx, jdx=0, nx2-1)], dl2 == j)

      ! Apply calculations for the corresponding rows and columns
      do idx = 1, size(is)
        do jdx = 1, size(js)
          cov(is(idx), js(jdx)) = Amplitudes(i*2) * Amplitudes(j*2) * gamma_g_g(is(idx), js(jdx)) + &
                                  Amplitudes(i*2+1) * Amplitudes(j*2+1) * gamma_dg_dg(is(idx), js(jdx)) + &
             (Amplitudes(i*2) * Amplitudes(j*2+1) - Amplitudes(i*2+1) * Amplitudes(j*2)) * gamma_g_dg(is(idx), js(jdx))
          if (i /= j) then
            cov(js(jdx), is(idx)) = cov(is(idx), js(jdx))  ! Mirror the lower triangle to the upper triangle
          end if
        end do
      end do
    end do
  end do

  else

   do i = 0, m - 1
    do j = 0, m - 1
      ! Integer indices for current dimensions i and j
      is = pack([(idx, idx=0, nx1-1)], dl1 == i)
      js = pack([(jdx, jdx=0, nx2-1)], dl2 == j)

      ! Apply calculations for the corresponding rows and columns
      do idx = 1, size(is)
        do jdx = 1, size(js)
          cov(is(idx), js(jdx)) = Amplitudes(i*2) * Amplitudes(j*2) * gamma_g_g(is(idx), js(jdx)) + &
                                  Amplitudes(i*2+1) * Amplitudes(j*2+1) * gamma_dg_dg(is(idx), js(jdx)) + &
             (Amplitudes(i*2) * Amplitudes(j*2+1) - Amplitudes(i*2+1) * Amplitudes(j*2)) * gamma_g_dg(is(idx), js(jdx))
        end do
      end do
    end do
  end do

  end if


  !ma1 = matmul(transpose(reshape(Amplitudes(dl1(:)*2),[nx1,1])),reshape(Amplitudes(dl2(:2)*2),[nx2,1]))
  !print *,reshape(Amplitudes(dl1(:)*2),[1,nx1])
  !print *,shape(reshape(Amplitudes(dl1(:)*2),[nx1,1]) )
  !ma1 = matmul(reshape(Amplitudes(dl1(:)*2),[nx1,1]),reshape(Amplitudes(dl2(:)*2),[1,nx2]) )
  !ma2 = matmul(reshape(Amplitudes(dl1(:)*2+1),[nx1,1]),reshape(Amplitudes(dl2(:)*2+1),[1,nx2]) )
  !ma3 = matmul(reshape(Amplitudes(dl1(:)*2),[nx1,1]),reshape(Amplitudes(dl2(:)*2+1),[1,nx2]) )
  !ma4 = matmul(reshape(Amplitudes(dl1(:)*2+1),[nx1,1]),reshape(Amplitudes(dl2(:)*2),[1,nx2]) )

  !cov = ma1 * gamma_g_g + ma2 * gamma_dg_dg + ( ma3 - ma4) * gamma_g_dg
  !print *, ma1
  !stop

end subroutine
