module qs3_shape_mod
    !---------------------------
    ! Three-quadratic-surface (3QS) shapes
    !   J. R. Nix, Nucl. Phys. A 130 (1969) 241
    !   P. Moller et al., Phys. Rev. C 79 (2009) 064304, Sec. II.A.3, Eqs. (3)-(15)
    !   P. Moller et al., ADNDT 109-110 (2016) 1 (FRDM2012), Sec. 2.2.2 and Eq. (75)
    !
    ! The surface is made of three smoothly joined quadratic surfaces of revolution about z:
    !   rho^2 = a1^2 - a1^2/c1^2 (z - l1)^2,   l1 - c1 <= z <= z1   (left end body, 1)
    !   rho^2 = a3^2 - a3^2/c3^2 (z - l3)^2,   z1 <= z <= z2        (middle body, 3)
    !   rho^2 = a2^2 - a2^2/c2^2 (z - l2)^2,   z2 <= z <= l2 + c2   (right end body, 2)
    ! The middle body is a spheroid (sigma2 > 0) or a hyperboloid of one sheet (sigma2 < 0).
    ! rho and d rho/dz are continuous at z1 and z2, the volume is 4/3 pi R^3 and the centre
    ! of mass is at z = 0.
    ! Shape coordinates (Eqs. 4, 5):
    !   u = [ (a1^2 + a2^2)/2 ]^(1/2)
    !   sigma1 = (l2 - l1)/u,  sigma2 = a3^2/c3^2,  sigma3 = (a1^2/c1^2 + a2^2/c2^2)/2
    !   alpha1 = (l1 + l2)/(2u) (fixed by the centre of mass),
    !   alpha2 = (a1^2 - a2^2)/u^2,  alpha3 = a1^2/c1^2 - a2^2/c2^2
    !
    ! Construction. With P_i(z) = rho^2 of surface i (a quadratic polynomial in z),
    ! smooth joining at z1 means that P1 - P3 has a double root at z1:
    !   P1 - P3 = k1 (z - z1)^2,  P2 - P3 = k2 (z - z2)^2,  k_i = sigma2 - a_i^2/c_i^2.
    ! Subtracting, P1 - P2 = k1 (z - z1)^2 - k2 (z - z2)^2; the z^2 terms agree identically,
    ! the z^1 and z^0 terms give two equations for z1 and z2 (a quadratic equation). The root
    ! with l1 - c1 <= z1 < z2 <= l2 + c2 and rho^2 > 0 on the middle body is taken.
    ! Each piece is stored as rho^2 = p0 + p1 z + p2 z^2 on [zb(i-1), zb(i)], so that the
    ! limiting case sigma2 = 0 (c3 -> infinity) needs no special treatment.
    !---------------------------
    ! 【日本語での説明】
    ! 3 つの 2 次曲面 (回転面) をなめらかにつないだ形 (3QS) を表すモジュール.
    ! z 軸まわりに軸対称. 左の端の物体 (1), 首の部分 (3), 右の端の物体 (2) からなる.
    ! 【形の式】 (PRC 79 式 3)
    !   rho^2 = a1^2 - (a1/c1)^2 (z - l1)^2   (l1 - c1 <= z <= z1)
    !   rho^2 = a3^2 - (a3/c3)^2 (z - l3)^2   (z1 <= z <= z2)    sigma2 < 0 なら一葉双曲面 (首)
    !   rho^2 = a2^2 - (a2/c2)^2 (z - l2)^2   (z2 <= z <= l2 + c2)
    ! z1, z2 で rho と d rho/dz が連続, 体積 4/3 pi R^3, 重心 z = 0.
    ! 【座標】 (PRC 79 式 4, 5)
    !   u = sqrt((a1^2 + a2^2)/2),  sigma1 = (l2 - l1)/u,  sigma2 = a3^2/c3^2,
    !   sigma3 = (a1^2/c1^2 + a2^2/c2^2)/2,  alpha1 = (l1 + l2)/(2u) (重心で決まる),
    !   alpha2 = (a1^2 - a2^2)/u^2,  alpha3 = a1^2/c1^2 - a2^2/c2^2
    ! 【作り方】P_i(z) を曲面 i の rho^2 (z の 2 次式) とすると, z1 でなめらかにつながる条件は
    ! 「P1 - P3 が z1 に 2 重根をもつ」こと: P1 - P3 = k1 (z - z1)^2, k1 = sigma2 - a1^2/c1^2.
    ! 同様に P2 - P3 = k2 (z - z2)^2. 2 式の差 P1 - P2 = k1 (z - z1)^2 - k2 (z - z2)^2 の
    ! 1 次と 0 次の係数から z1, z2 の 2 次方程式が出る. l1 - c1 <= z1 < z2 <= l2 + c2 で,
    ! 首の部分で rho^2 > 0 になる根をとる. 各部分は rho^2 = p0 + p1 z + p2 z^2 で持つ.
    !
    ! Surface parametrisation for surface_shape%eval: u = 0 is the right pole (z = l2 + c2),
    ! u = pi the left pole. [0, pi] is split at ub(1), ub(2) (multiples of pi/8, so that the
    ! panels of potential_gh_mod do not straddle the joints); in each end body
    ! z - z_pole ~ (1 - cos(pi t / 2)), so rho ~ t is smooth at the pole, the middle body is
    ! linear in t. N = dr/du x dr/dv = (-dz/du) (rho cos v, rho sin v, -P'(z)/2).
    !---------------------------
    use iso_fortran_env, only: real64
    use constant_mod, only: pi
    use shape_mod, only: surface_shape, gauss_legendre
    implicit none
    private

    integer, parameter :: dp = real64

    ! reference nucleus for the renormalised quadrupole moment (Fig. 2 of PRC 79)
    real(dp), parameter, public :: qs3_Z_ref = 94.0_dp, qs3_A_ref = 240.0_dp, qs3_r0_ref = 1.16_dp

    type, extends(surface_shape), public :: qs3_shape
        logical :: valid = .false.
        real(dp) :: R = 0.0_dp               ! fm, radius of the sphere of the same volume
        ! shape coordinates, Eq. (5)
        real(dp) :: sigma1 = 0.0_dp, sigma2 = 0.0_dp, sigma3 = 0.0_dp
        real(dp) :: alpha1 = 0.0_dp, alpha2 = 0.0_dp, alpha3 = 0.0_dp
        ! end bodies (fm) and the auxiliary unit u (fm)
        real(dp) :: a1 = 0.0_dp, c1 = 0.0_dp, l1 = 0.0_dp
        real(dp) :: a2 = 0.0_dp, c2 = 0.0_dp, l2 = 0.0_dp
        real(dp) :: u = 0.0_dp
        ! zb(0) = l1 - c1, zb(1) = z1, zb(2) = z2, zb(3) = l2 + c2  (fm)
        real(dp) :: zb(0:3) = 0.0_dp
        ! rho^2 = p(0,i) + p(1,i) z + p(2,i) z^2 on [zb(i-1), zb(i)]; i = 1 left, 2 middle, 3 right
        real(dp) :: p(0:2, 3) = 0.0_dp
        ! eval: u-break points (u = 0 right pole, u = pi left pole)
        real(dp) :: ub(2) = 0.0_dp
        contains
            procedure :: eval => qs3_eval
            procedure :: inside => qs3_inside
            procedure :: volume => qs3_volume
            procedure :: relative_area => qs3_relative_area
            procedure :: rho2 => qs3_rho2
            procedure :: q2_hat => qs3_q2_hat
            procedure :: q2_barn => qs3_q2_barn
            procedure :: neck_diameter => qs3_neck_diameter
            procedure :: b_w => qs3_b_w
            procedure :: eps_f => qs3_eps_f
            procedure :: alpha_g => qs3_alpha_g
            procedure :: axial_nodes => qs3_axial_nodes
            procedure :: profile_length => qs3_profile_length
            procedure :: scaled => qs3_scaled
            procedure :: reflection_symmetric => qs3_reflection_symmetric
    end type qs3_shape

    public :: qs3_from_sigma
    public :: sigma3_alpha3_from_eps_f
    public :: alpha2_from_alpha_g

    contains

        !==========================================================================
        ! construction
        !==========================================================================

        function qs3_from_sigma(R, sigma1, sigma2, sigma3, alpha2, alpha3) result(s)
            !---------------------------
            ! 3QS shape of volume 4/3 pi R^3 with its centre of mass at z = 0.
            ! s%valid = .false. if no shape exists for these coordinates.
            !---------------------------
            ! 【目的】(sigma1, sigma2, sigma3, alpha2, alpha3) から体積 4/3 pi R^3, 重心 z = 0 の
            ! 3QS 形状を作る. 形がないとき (端の物体の外で接続する, 首が切れる など) は valid = .false.
            ! 【手順】u = 1 の単位で l1 = -sigma1/2, l2 = sigma1/2 として z1, z2 を求め,
            ! 重心を原点に移し (これが alpha1 を決める), 体積が合うように全体を拡大する.
            implicit none
            real(dp), intent(in) :: R, sigma1, sigma2, sigma3, alpha2, alpha3
            type(qs3_shape) :: s
            real(dp), parameter :: tiny = 1.0e-12_dp
            real(dp) :: AA1, AA2, a1sq, a2sq, c1, c2, l1, l2, k1, k2, B, Cc
            real(dp) :: qa, qb, qc, disc, roots(2), z1, z2, P1(0:2), P2(0:2), P3(0:2)
            real(dp) :: best(0:2, 3), best_zb(0:3), vol, zc, lam
            integer :: nr, ir, n_good

            s%valid = .false.
            s%R = R
            s%sigma1 = sigma1
            s%sigma2 = sigma2
            s%sigma3 = sigma3
            s%alpha2 = alpha2
            s%alpha3 = alpha3
            AA1 = sigma3 + 0.5_dp * alpha3          ! a1^2/c1^2
            AA2 = sigma3 - 0.5_dp * alpha3          ! a2^2/c2^2
            a1sq = 1.0_dp + 0.5_dp * alpha2         ! u = 1
            a2sq = 1.0_dp - 0.5_dp * alpha2
            if (AA1 <= 0.0_dp .or. AA2 <= 0.0_dp .or. a1sq <= 0.0_dp .or. a2sq <= 0.0_dp) return
            if (sigma1 < 0.0_dp) return
            c1 = sqrt(a1sq / AA1)
            c2 = sqrt(a2sq / AA2)
            l1 = -0.5_dp * sigma1
            l2 = 0.5_dp * sigma1
            k1 = sigma2 - AA1
            k2 = sigma2 - AA2
            P1 = [a1sq - AA1 * l1 * l1, 2.0_dp * AA1 * l1, -AA1]
            P2 = [a2sq - AA2 * l2 * l2, 2.0_dp * AA2 * l2, -AA2]
            ! P1 - P2 = k1 (z - z1)^2 - k2 (z - z2)^2:
            !   z^1: k2 z2 - k1 z1 = B,   z^0: k1 z1^2 - k2 z2^2 = Cc
            B = AA1 * l1 - AA2 * l2
            Cc = P1(0) - P2(0)
            if (abs(k1) < tiny .and. abs(k2) < tiny) return

            n_good = 0
            if (abs(k2) >= abs(k1)) then
                ! z2 = (B + k1 z1)/k2:  k1 (k2 - k1) z1^2 - 2 B k1 z1 - (B^2 + Cc k2) = 0
                qa = k1 * (k2 - k1)
                qb = -2.0_dp * B * k1
                qc = -(B * B + Cc * k2)
            else
                ! z1 = (k2 z2 - B)/k1:  k2 (k2 - k1) z2^2 - 2 B k2 z2 + (B^2 - Cc k1) = 0
                qa = k2 * (k2 - k1)
                qb = -2.0_dp * B * k2
                qc = B * B - Cc * k1
            end if
            call solve_quadratic(qa, qb, qc, roots, nr)
            do ir = 1, nr
                if (abs(k2) >= abs(k1)) then
                    z1 = roots(ir)
                    z2 = (B + k1 * z1) / k2
                else
                    z2 = roots(ir)
                    z1 = (k2 * z2 - B) / k1
                end if
                if (.not. (z1 >= l1 - c1 .and. z1 <= l1 + c1)) cycle
                if (.not. (z2 >= l2 - c2 .and. z2 <= l2 + c2)) cycle
                if (.not. (z1 < z2)) cycle
                P3 = P1 - k1 * [z1 * z1, -2.0_dp * z1, 1.0_dp]
                if (min_on(P3, z1, z2) <= tiny) cycle
                n_good = n_good + 1
                best(:, 1) = P1
                best(:, 2) = P3
                best(:, 3) = P2
                best_zb = [l1 - c1, z1, z2, l2 + c2]
                if (n_good == 1) exit
            end do
            if (n_good == 0) return

            s%p = best
            s%zb = best_zb
            s%a1 = sqrt(a1sq)
            s%c1 = c1
            s%l1 = l1
            s%a2 = sqrt(a2sq)
            s%c2 = c2
            s%l2 = l2
            s%u = 1.0_dp

            ! centre of mass to z = 0, then scale to the volume 4/3 pi R^3
            vol = pi * poly_moment(s, 0)
            zc = pi * poly_moment(s, 1) / vol
            call shift(s, zc)
            s%alpha1 = 0.5_dp * (s%l1 + s%l2) / s%u
            lam = (4.0_dp / 3.0_dp * pi * R**3 / vol)**(1.0_dp / 3.0_dp)
            call scale(s, lam)
            call set_breaks(s)
            s%valid = .true.
        end function qs3_from_sigma

        pure subroutine solve_quadratic(a, b, c, x, n)
            ! real roots of a x^2 + b x + c = 0 (linear if a = 0), n = number of roots
            real(dp), intent(in) :: a, b, c
            real(dp), intent(out) :: x(2)
            integer, intent(out) :: n
            real(dp) :: disc, q, scale_ab
            x = 0.0_dp
            n = 0
            scale_ab = max(abs(b), abs(c), 1.0e-300_dp)
            if (abs(a) <= 1.0e-13_dp * scale_ab) then
                if (b /= 0.0_dp) then
                    n = 1
                    x(1) = -c / b
                end if
                return
            end if
            disc = b * b - 4.0_dp * a * c
            if (disc < 0.0_dp) return
            q = -0.5_dp * (b + sign(sqrt(disc), b))
            n = 2
            x(1) = q / a
            if (q /= 0.0_dp) then
                x(2) = c / q
            else
                x(2) = x(1)
            end if
        end subroutine solve_quadratic

        pure function min_on(P, za, zb) result(m)
            ! minimum of the quadratic P on [za, zb]
            real(dp), intent(in) :: P(0:2), za, zb
            real(dp) :: m, zv
            m = min(peval(P, za), peval(P, zb))
            if (P(2) > 0.0_dp) then
                zv = -P(1) / (2.0_dp * P(2))
                if (zv > za .and. zv < zb) m = min(m, peval(P, zv))
            end if
        end function min_on

        pure function peval(P, z) result(v)
            real(dp), intent(in) :: P(0:2), z
            real(dp) :: v
            v = P(0) + z * (P(1) + z * P(2))
        end function peval

        pure function poly_moment(s, k) result(m)
            ! int z^k rho^2 dz over the three pieces (Gauss-Legendre, exact for these degrees)
            ! k = 0, 1, 2; k = -1 gives int rho^4 dz
            type(qs3_shape), intent(in) :: s
            integer, intent(in) :: k
            real(dp) :: m, z, h, c, P
            real(dp), parameter :: x5(5) = [-0.9061798459386640_dp, -0.5384693101056831_dp, 0.0_dp, &
                                            0.5384693101056831_dp, 0.9061798459386640_dp]
            real(dp), parameter :: w5(5) = [0.2369268850561891_dp, 0.4786286704993665_dp, 0.5688888888888889_dp, &
                                            0.4786286704993665_dp, 0.2369268850561891_dp]
            integer :: i, j
            m = 0.0_dp
            do i = 1, 3
                h = 0.5_dp * (s%zb(i) - s%zb(i - 1))
                c = 0.5_dp * (s%zb(i) + s%zb(i - 1))
                do j = 1, 5
                    z = c + h * x5(j)
                    P = peval(s%p(:, i), z)
                    if (k >= 0) then
                        m = m + h * w5(j) * z**k * P
                    else
                        m = m + h * w5(j) * P * P
                    end if
                end do
            end do
        end function poly_moment

        pure subroutine shift(s, zc)
            ! z -> z - zc
            type(qs3_shape), intent(inout) :: s
            real(dp), intent(in) :: zc
            integer :: i
            real(dp) :: q(0:2)
            do i = 1, 3
                q = s%p(:, i)
                s%p(0, i) = q(0) + q(1) * zc + q(2) * zc * zc
                s%p(1, i) = q(1) + 2.0_dp * q(2) * zc
            end do
            s%zb = s%zb - zc
            s%l1 = s%l1 - zc
            s%l2 = s%l2 - zc
        end subroutine shift

        pure subroutine scale(s, lam)
            ! all lengths times lam:  rho^2(z) -> lam^2 rho^2(z/lam)
            type(qs3_shape), intent(inout) :: s
            real(dp), intent(in) :: lam
            s%p(0, :) = lam * lam * s%p(0, :)
            s%p(1, :) = lam * s%p(1, :)
            s%zb = lam * s%zb
            s%a1 = lam * s%a1
            s%c1 = lam * s%c1
            s%l1 = lam * s%l1
            s%a2 = lam * s%a2
            s%c2 = lam * s%c2
            s%l2 = lam * s%l2
            s%u = lam * s%u
        end subroutine scale

        pure subroutine set_breaks(s)
            ! u-break points: each piece gets a multiple of pi/8 roughly proportional to its length
            type(qs3_shape), intent(inout) :: s
            real(dp) :: L
            integer :: n_right, n_mid
            L = s%zb(3) - s%zb(0)
            n_right = max(1, min(6, nint(8.0_dp * (s%zb(3) - s%zb(2)) / L)))
            n_mid = max(1, min(7 - n_right, nint(8.0_dp * (s%zb(2) - s%zb(1)) / L)))
            s%ub(1) = pi * n_right / 8.0_dp
            s%ub(2) = pi * (n_right + n_mid) / 8.0_dp
        end subroutine set_breaks

        function qs3_scaled(this, R) result(s)
            ! the same shape with the volume 4/3 pi R^3
            ! 【目的】同じ形で体積を 4/3 pi R^3 にしたもの (一粒子ポテンシャルの半径 R_pot 用).
            class(qs3_shape), intent(in) :: this
            real(dp), intent(in) :: R
            type(qs3_shape) :: s
            s = this
            if (.not. s%valid) return
            call scale(s, R / this%R)
            s%R = R
        end function qs3_scaled

        !==========================================================================
        ! Moller's coordinates (eps_f1, eps_f2, alpha_g) -> (sigma3, alpha3, alpha2)
        !==========================================================================

        pure subroutine sigma3_alpha3_from_eps_f(eps_f1, eps_f2, sigma3, alpha3)
            ! Eqs. (6), (7): a/c of a spheroid with Nilsson's eps is (3 - 2 eps)/(3 + eps)
            ! 【目的】端の物体の Nilsson 変形 eps_f1, eps_f2 から sigma3, alpha3 (PRC 79 式 6, 7).
            real(dp), intent(in) :: eps_f1, eps_f2
            real(dp), intent(out) :: sigma3, alpha3
            real(dp) :: t1, t2
            t1 = ((3.0_dp - 2.0_dp * eps_f1) / (3.0_dp + eps_f1))**2
            t2 = ((3.0_dp - 2.0_dp * eps_f2) / (3.0_dp + eps_f2))**2
            sigma3 = 0.5_dp * (t1 + t2)
            alpha3 = t1 - t2
        end subroutine sigma3_alpha3_from_eps_f

        pure function alpha2_from_alpha_g(alpha_g, sigma3, alpha3) result(alpha2)
            ! Eq. (11), from alpha_g = (a1^2 c1 - a2^2 c2)/(a1^2 c1 + a2^2 c2)  (Eq. 10)
            ! 【目的】質量非対称 alpha_g から alpha2 (PRC 79 式 11).
            ! a_i^3 / sqrt(A_i) の比が (1 + alpha_g)/(1 - alpha_g) なので
            ! t = (a1/a2)^2 = [ (1 + alpha_g)^2 (2 sigma3 + alpha3) / ((1 - alpha_g)^2 (2 sigma3 - alpha3)) ]^(1/3),
            ! alpha2 = 2 (t - 1)/(t + 1).
            real(dp), intent(in) :: alpha_g, sigma3, alpha3
            real(dp) :: alpha2, t
            t = ((alpha_g + 1.0_dp)**2 * (2.0_dp * sigma3 + alpha3) &
               / ((alpha_g - 1.0_dp)**2 * (2.0_dp * sigma3 - alpha3)))**(1.0_dp / 3.0_dp)
            alpha2 = 2.0_dp * (t - 1.0_dp) / (t + 1.0_dp)
        end function alpha2_from_alpha_g

        !==========================================================================
        ! geometry
        !==========================================================================

        pure function qs3_rho2(this, z) result(r2)
            ! rho^2(z); negative outside [zb(0), zb(3)]
            class(qs3_shape), intent(in) :: this
            real(dp), intent(in) :: z
            real(dp) :: r2
            integer :: i
            if (z < this%zb(0) .or. z > this%zb(3)) then
                r2 = -1.0_dp
                return
            end if
            i = 2
            if (z < this%zb(1)) i = 1
            if (z > this%zb(2)) i = 3
            r2 = peval(this%p(:, i), z)
        end function qs3_rho2

        pure subroutine piece_map(this, i, t, z, dzdt)
            ! z(t), t in [0, 1], on piece i, measured from the side nearer to u = 0
            ! (i = 3: t = 0 at the right pole; i = 2: t = 0 at z2; i = 1: t = 0 at z1, t = 1 left pole)
            class(qs3_shape), intent(in) :: this
            integer, intent(in) :: i
            real(dp), intent(in) :: t
            real(dp), intent(out) :: z, dzdt
            real(dp) :: Lp
            Lp = this%zb(i) - this%zb(i - 1)
            select case (i)
            case (3)
                z = this%zb(3) - Lp * (1.0_dp - cos(0.5_dp * pi * t))
                dzdt = -Lp * 0.5_dp * pi * sin(0.5_dp * pi * t)
            case (2)
                z = this%zb(2) - Lp * t
                dzdt = -Lp
            case default
                z = this%zb(0) + Lp * (1.0_dp - cos(0.5_dp * pi * (1.0_dp - t)))
                dzdt = -Lp * 0.5_dp * pi * sin(0.5_dp * pi * (1.0_dp - t))
            end select
        end subroutine piece_map

        pure subroutine qs3_eval(this, u, v, r, N)
            ! surface point and N = dr/du x dr/dv (see the module header)
            class(qs3_shape), intent(in) :: this
            real(dp), intent(in) :: u, v
            real(dp), intent(out) :: r(3), N(3)
            real(dp) :: t, z, dzdt, dtdu, P, dPdz, rho, cv, sv
            integer :: i
            if (u <= this%ub(1)) then
                i = 3
                t = u / this%ub(1)
                dtdu = 1.0_dp / this%ub(1)
            else if (u <= this%ub(2)) then
                i = 2
                t = (u - this%ub(1)) / (this%ub(2) - this%ub(1))
                dtdu = 1.0_dp / (this%ub(2) - this%ub(1))
            else
                i = 1
                t = (u - this%ub(2)) / (pi - this%ub(2))
                dtdu = 1.0_dp / (pi - this%ub(2))
            end if
            t = max(0.0_dp, min(1.0_dp, t))
            call piece_map(this, i, t, z, dzdt)
            P = max(0.0_dp, peval(this%p(:, i), z))
            dPdz = this%p(1, i) + 2.0_dp * this%p(2, i) * z
            rho = sqrt(P)
            cv = cos(v)
            sv = sin(v)
            r = [rho * cv, rho * sv, z]
            N = (-dzdt * dtdu) * [rho * cv, rho * sv, -0.5_dp * dPdz]
        end subroutine qs3_eval

        pure function qs3_inside(this, x) result(yes)
            class(qs3_shape), intent(in) :: this
            real(dp), intent(in) :: x(3)
            logical :: yes
            yes = x(1) * x(1) + x(2) * x(2) <= this%rho2(x(3))
        end function qs3_inside

        function qs3_volume(this) result(vol)
            ! exact: pi int rho^2 dz
            class(qs3_shape), intent(in) :: this
            real(dp) :: vol
            vol = pi * poly_moment(this, 0)
        end function qs3_volume

        function qs3_relative_area(this) result(B_s)
            ! B_s = S / (4 pi R^2),  S = 2 pi int sqrt(P + P'^2/4) dz  (smooth on each piece)
            ! 【目的】相対表面積 B_s (FRDM2012 式 70). 各部分で被積分関数はなめらか.
            class(qs3_shape), intent(in) :: this
            real(dp) :: B_s, S, h, c, z, P, dPdz, Rv
            integer, parameter :: nq = 48
            real(dp) :: x(nq), w(nq)
            integer :: i, j
            call gauss_legendre(nq, x, w)
            S = 0.0_dp
            do i = 1, 3
                h = 0.5_dp * (this%zb(i) - this%zb(i - 1))
                c = 0.5_dp * (this%zb(i) + this%zb(i - 1))
                do j = 1, nq
                    z = c + h * x(j)
                    P = max(0.0_dp, peval(this%p(:, i), z))
                    dPdz = this%p(1, i) + 2.0_dp * this%p(2, i) * z
                    S = S + h * w(j) * sqrt(P + 0.25_dp * dPdz * dPdz)
                end do
            end do
            S = 2.0_dp * pi * S
            Rv = (3.0_dp * this%volume() / (4.0_dp * pi))**(1.0_dp / 3.0_dp)
            B_s = S / (4.0_dp * pi * Rv * Rv)
        end function qs3_relative_area

        function qs3_q2_hat(this) result(q)
            ! q = int (2 z^2 - rho^2) dV / (V R^2), z from the centre of mass (at 0)
            !   int z^2 dV = pi int z^2 P dz,  int rho^2 dV = pi/2 int P^2 dz
            ! 【目的】無次元の四重極モーメント q = int (2 z^2 - rho^2) dV / (V R^2).
            ! Q2 = Z R^2 q [e fm^2] (PRC 79 式 14 で lambda = 2, 体積 = 4/3 pi r0^3 A).
            class(qs3_shape), intent(in) :: this
            real(dp) :: q, vol
            vol = this%volume()
            q = pi * (2.0_dp * poly_moment(this, 2) - 0.5_dp * poly_moment(this, -1)) / (vol * this%R**2)
        end function qs3_q2_hat

        function qs3_q2_barn(this) result(Q2)
            ! charge quadrupole moment in b of 240Pu with the same shape (Fig. 2 of PRC 79):
            !   Q2 = Z_ref (r0 A_ref^(1/3))^2 q / 100
            ! 【目的】同じ形の 240Pu の電気四重極モーメント [b] (PRC 79 図 2 の規格化).
            class(qs3_shape), intent(in) :: this
            real(dp) :: Q2
            Q2 = qs3_Z_ref * (qs3_r0_ref * qs3_A_ref**(1.0_dp / 3.0_dp))**2 * this%q2_hat() / 100.0_dp
        end function qs3_q2_barn

        function qs3_neck_diameter(this) result(d)
            ! d = 2 min rho on the middle body [fm] (the neck diameter when the middle body has a waist)
            ! 【目的】首の直径 d = 2 x (首の部分 z1 <= z <= z2 での rho の最小値) [fm].
            class(qs3_shape), intent(in) :: this
            real(dp) :: d
            d = 2.0_dp * sqrt(max(0.0_dp, min_on(this%p(:, 2), this%zb(1), this%zb(2))))
        end function qs3_neck_diameter

        function qs3_b_w(this, a_d) result(B_w)
            !---------------------------
            ! shape dependence of the A^0 and Wigner terms, FRDM2012 Eq. (75):
            !   B_W = (1 - S3/S1)^2 a_d + 1   (sigma2 <= 0),   B_W = 1   (sigma2 >= 0)
            ! S3: area of the neck cross section, S1: maximum cross section pi a^2 of the
            ! smaller end body (the one with the smaller completed volume a^2 c).
            !---------------------------
            ! 【目的】A^0 項と Wigner 項の形状因子 B_W (FRDM2012 式 75). 切断直前 (S3 -> 0) で 1 + a_d.
            class(qs3_shape), intent(in) :: this
            real(dp), intent(in) :: a_d
            real(dp) :: B_w, S3, S1
            if (this%sigma2 >= 0.0_dp) then
                B_w = 1.0_dp
                return
            end if
            S3 = 0.25_dp * this%neck_diameter()**2
            if (this%a1**2 * this%c1 <= this%a2**2 * this%c2) then
                S1 = this%a1**2
            else
                S1 = this%a2**2
            end if
            B_w = (1.0_dp - min(1.0_dp, S3 / S1))**2 * a_d + 1.0_dp
        end function qs3_b_w

        function qs3_eps_f(this, i) result(e)
            ! Nilsson eps of end body i: a/c = (3 - 2 eps)/(3 + eps)  ->  eps = 3 (c - a)/(2 c + a)
            class(qs3_shape), intent(in) :: this
            integer, intent(in) :: i
            real(dp) :: e
            if (i == 1) then
                e = 3.0_dp * (this%c1 - this%a1) / (2.0_dp * this%c1 + this%a1)
            else
                e = 3.0_dp * (this%c2 - this%a2) / (2.0_dp * this%c2 + this%a2)
            end if
        end function qs3_eps_f

        function qs3_alpha_g(this) result(ag)
            ! Eq. (10)
            class(qs3_shape), intent(in) :: this
            real(dp) :: ag
            ag = (this%a1**2 * this%c1 - this%a2**2 * this%c2) / (this%a1**2 * this%c1 + this%a2**2 * this%c2)
        end function qs3_alpha_g

        pure function qs3_reflection_symmetric(this) result(yes)
            ! z -> -z symmetry (alpha2 = alpha3 = 0)
            class(qs3_shape), intent(in) :: this
            logical :: yes
            yes = abs(this%alpha2) < 1.0e-12_dp .and. abs(this%alpha3) < 1.0e-12_dp
        end function qs3_reflection_symmetric

        function qs3_profile_length(this) result(L)
            ! length of the meridian from pole to pole, int sqrt(dz^2 + d rho^2) [fm]
            class(qs3_shape), intent(in) :: this
            real(dp) :: L
            real(dp) :: Lp(3)
            call piece_lengths(this, Lp)
            L = sum(Lp)
        end function qs3_profile_length

        subroutine piece_lengths(this, Lp)
            ! meridian length of each piece (Gauss-Legendre in t; smooth with the pole mapping)
            class(qs3_shape), intent(in) :: this
            real(dp), intent(out) :: Lp(3)
            integer, parameter :: nq = 32
            real(dp) :: x(nq), wx(nq), t, zz, dzdt, P, dPdz
            integer :: i, j
            call gauss_legendre(nq, x, wx)
            Lp = 0.0_dp
            do i = 1, 3
                do j = 1, nq
                    t = 0.5_dp * (x(j) + 1.0_dp)
                    call piece_map(this, i, t, zz, dzdt)
                    P = max(0.0_dp, peval(this%p(:, i), zz))
                    dPdz = this%p(1, i) + 2.0_dp * this%p(2, i) * zz
                    ! ds/dt = |dz/dt| sqrt(1 + (P'/2)^2 / P) = |dz/dt| sqrt(P + P'^2/4) / sqrt(P)
                    if (P > 0.0_dp) then
                        Lp(i) = Lp(i) + 0.5_dp * wx(j) * abs(dzdt) * sqrt(1.0_dp + 0.25_dp * dPdz**2 / P)
                    end if
                end do
            end do
        end subroutine piece_lengths

        subroutine qs3_axial_nodes(this, n_total, z, rho, pz, w)
            !---------------------------
            ! rings for surface integrals of the axially symmetric surface:
            ! about n_total Gauss-Legendre nodes in t, shared among the three pieces in
            ! proportion to their meridian lengths (at least 8 per piece),
            ! oint F dS = sum_k w_k int dphi F(rho_k cos phi, rho_k sin phi, z_k) with the vector
            ! area element (rho_k cos phi, rho_k sin phi, pz_k) w_k dphi, pz = -P'(z)/2.
            ! The integrands are smooth on each piece (the curvature jumps at z1, z2).
            !---------------------------
            ! 【目的】軸対称な表面の面積分のためのリング. 合計約 n_total 本を, 各部分の
            ! 子午線の長さに比例して分ける (各部分 8 本以上).
            ! dS = (rho cos phi, rho sin phi, -P'(z)/2) w dphi. frldm_mod の B_1, B_3 で使う.
            class(qs3_shape), intent(in) :: this
            integer, intent(in) :: n_total
            real(dp), allocatable, intent(out) :: z(:), rho(:), pz(:), w(:)
            real(dp), allocatable :: x(:), wx(:)
            real(dp) :: t, zz, dzdt, Lp(3)
            integer :: i, j, m, np(3)
            call piece_lengths(this, Lp)
            do i = 1, 3
                np(i) = max(8, nint(n_total * Lp(i) / sum(Lp)))
            end do
            allocate(z(sum(np)), rho(sum(np)), pz(sum(np)), w(sum(np)))
            m = 0
            do i = 3, 1, -1
                if (allocated(x)) deallocate(x, wx)
                allocate(x(np(i)), wx(np(i)))
                call gauss_legendre(np(i), x, wx)
                do j = 1, np(i)
                    m = m + 1
                    t = 0.5_dp * (x(j) + 1.0_dp)
                    call piece_map(this, i, t, zz, dzdt)
                    z(m) = zz
                    rho(m) = sqrt(max(0.0_dp, peval(this%p(:, i), zz)))
                    pz(m) = -0.5_dp * (this%p(1, i) + 2.0_dp * this%p(2, i) * zz)
                    w(m) = abs(dzdt) * 0.5_dp * wx(j)
                end do
            end do
        end subroutine qs3_axial_nodes

end module qs3_shape_mod
