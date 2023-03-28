! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, mc.f90, is copyright (c) 2019-2023, Miguel A. Caro
! HND X
! HND X   This file was authored by Tigany Zarrouk
! HND X
! HND X   TurboGAP is distributed in the hope that it will be useful for non-commercial
! HND X   academic research, but WITHOUT ANY WARRANTY; without even the implied
! HND X   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! HND X   ASL for more details.
! HND X
! HND X   You should have received a copy of the ASL along with this program
! HND X   (e.g. in a LICENSE.md file); if not, you can write to the original
! HND X   licensor, Miguel Caro (mcaroba@gmail.com). The ASL is also published at
! HND X   http://github.com/gabor1/ASL
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   Miguel A. Caro. Phys. Rev. B 100, 024112 (2019)
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module mc


  use neighbors


  contains


  subroutine monte_carlo_insert(e_new, e_prev, T, mu, m, volume, volume_bias, N_exch)
    implicit none

    real*8, intent(in) :: e_new, e_prev, T, mu, m, volume_bias
    integer intent(in) :: N_exch
    real*8 :: lam
    real*8 :: kB = 8.617333262e-5, hbar = 6.582119569e-1, pi=3.1415926535
    real*8 intent(out) :: p_accept
    ! mass has units eV*fs^2/A^2
    ! hbar in eV.fs
    ! lam is thermal debroglie wavelength
    lam = sqrt( ( 2.0 * pi * hbar*hbar ) / ( m * kB * T ) )

    p_accept =  (V*volume_bias) / (lam**3.0 * (N_exch + 1)) * &
                exp( -( e_new - e_prev - mu  ) / (kB * T) )

  end subroutine monte_carlo_insert


  subroutine monte_carlo_remove(e_new, e_prev, T, mu, m, volume, volume_bias, N_exch)
    implicit none

    real*8, intent(in) :: e_new, e_prev, T, mu, m, volume_bias
    integer intent(in) :: N_exch
    real*8 :: lam
    real*8 :: kB = 8.617333262e-5, hbar = 6.582119569e-1, pi=3.1415926535
    real*8 intent(out) :: p_accept
    ! mass has units eV*fs^2/A^2
    ! hbar in eV.fs
    ! lam is thermal debroglie wavelength
    lam = sqrt( ( 2.0 * pi * hbar*hbar ) / ( m * kB * T ) )

    p_accept = (lam**3.0 * N_exch ) / (V*volume_bias)  * &
                exp( -( e_new - e_prev + mu  ) / (kB * T) )

  end subroutine monte_carlo_remove

  subroutine monte_carlo_move(e_new, e_prev, T)
    implicit none

    real*8, intent(in) :: e_new, e_prev, T,
    real*8 :: kB = 8.617333262e-5
    real*8 intent(out) :: p_accept
    p_accept = exp( -( e_new - e_prev ) / (kB * T) )

  end subroutine monte_carlo_move

  subroutine monte_carlo_volume(e_new, e_prev, T, V_new, V_prev, P, N_exch):
    real*8, intent(in) :: e_new, e_prev, T, V_new, V_prev, P
    real*8 :: kB = 8.617333262e-5, beta
    integer intent(in) :: N_exch
    real*8 intent(out) :: p_accept

    beta = (1./(kB * T))
    p_accept = exp( - beta * ( (e_new - e_prev) + &
                      P * ( V_new-V_prev ) + &
                     -(N_exch+1) * log( Vf/V )/ beta ) )
  end subroutine monte_carlo_volume


end module mc
