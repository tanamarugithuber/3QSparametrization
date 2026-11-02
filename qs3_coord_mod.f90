module qs3_coord_mod
    !---------------------------
    ! Moller's five fission shape coordinates for 3QS shapes, PRC 79 (2009) 064304, Sec. II.A.3:
    !   Q2     : charge quadrupole moment of 240Pu with the same shape [b]  (elongation)
    !   d      : neck diameter
    !   eps_f1, eps_f2 : Nilsson eps of the left and right end bodies (nascent fragments)
    !   alpha_g: (M1 - M2)/(M1 + M2), M_i = volume of end body i completed to a spheroid
    ! (eps_f1, eps_f2, alpha_g) give sigma3, alpha3, alpha2 in closed form (Eqs. 6, 7, 11).
    ! For given (Q2, eps_f1, eps_f2, alpha_g):
    !   * sigma1(sigma2) is found numerically so that Q2 is obtained ("for each combination d,
    !     alpha_g, eps_f1, eps_f2 we determine by a numerical procedure 45 values of sigma1 so
    !     that 45 preselected values of Q2 are obtained");
    !   * the neck diameter d(sigma2) = 2 min rho of the middle body grows with sigma2 from the
    !     smallest neck the parametrization allows (zero at scission, or a finite value for short
    !     shapes) up to a maximum d_max near sigma2 ~ 0; beyond it (bulging middle bodies) d
    !     decreases again because the "neck" is then the joint with an end body.
    !     We take the increasing branch [d_min, d_max] and choose n_d values of sigma2 so that d
    !     is equally spaced on it ("we select 15 values of sigma2 so that for each combination of
    !     Q2, eps_f1, eps_f2 and alpha_g we obtain a suitable spacing for d in the range between
    !     the smallest and largest value possible").
    ! The paper does not give the precise neck definition and spacing (Ref. [28] of FRDM2012 for
    ! "details"); the choice above is ours and is documented in README_3QS.md.
    !---------------------------
    ! 【日本語での説明】
    ! PRC 79 の 5 つの形の座標 (Q2, d, eps_f1, eps_f2, alpha_g) から 3QS の sigma, alpha を求める.
    !   eps_f1, eps_f2, alpha_g -> sigma3, alpha3, alpha2 は式 6, 7, 11 で直接決まる.
    !   Q2 は sigma1 を数値的に合わせる (sigma2 ごとに).
    !   d (首の直径 = 首の部分での 2 rho の最小値) は sigma2 とともに増え, sigma2 ~ 0 付近で最大
    !   (d_max) になったあと, 中央がふくらむ形では端の物体とのつなぎ目が最小になるので減る.
    !   増える側 [d_min, d_max] を n_d 等分するように sigma2 を選ぶ (論文の「15 個の sigma2 を
    !   d が最小値と最大値の間でほどよい間隔になるように選ぶ」の我々の解釈).
    !---------------------------
    use iso_fortran_env, only: real64
    use qs3_shape_mod, only: qs3_shape, qs3_from_sigma, sigma3_alpha3_from_eps_f, alpha2_from_alpha_g
    implicit none
    private

    integer, parameter :: dp = real64

    real(dp), parameter, public :: sigma1_max = 12.0_dp     ! scan range of sigma1
    real(dp), parameter, public :: sigma1_step = 0.1_dp
    real(dp), parameter, public :: sigma2_lowest = -1000.0_dp ! lower end of the sigma2 scan
    integer, parameter :: n_sig2_scan = 64

    ! the increasing branch d(sigma2) at fixed (Q2, eps_f1, eps_f2, alpha_g)
    type, public :: neck_branch
        logical :: ok = .false.
        real(dp) :: R = 0.0_dp, Q2 = 0.0_dp
        real(dp) :: sigma3 = 0.0_dp, alpha2 = 0.0_dp, alpha3 = 0.0_dp
        real(dp) :: sig2_min = 0.0_dp, sig2_max = 0.0_dp   ! sigma2 at d_min and d_max
        real(dp) :: d_min = 0.0_dp, d_max = 0.0_dp         ! fm
    end type neck_branch

    public :: moller_to_sigma
    public :: shape_at_q2
    public :: find_neck_branch
    public :: shape_at_neck

    contains

        pure subroutine moller_to_sigma(eps_f1, eps_f2, alpha_g, sigma3, alpha3, alpha2)
            ! Eqs. (6), (7), (11)
            real(dp), intent(in) :: eps_f1, eps_f2, alpha_g
            real(dp), intent(out) :: sigma3, alpha3, alpha2
            call sigma3_alpha3_from_eps_f(eps_f1, eps_f2, sigma3, alpha3)
            alpha2 = alpha2_from_alpha_g(alpha_g, sigma3, alpha3)
        end subroutine moller_to_sigma

        function shape_at_q2(R, Q2, sigma2, sigma3, alpha2, alpha3) result(s)
            !---------------------------
            ! 3QS shape with the quadrupole moment Q2 [b, 240Pu normalisation] at fixed sigma2:
            ! scan sigma1 = 0, 0.1, ... for the first interval of valid shapes that brackets Q2,
            ! then bisection. s%valid = .false. if Q2 is not reached.
            !---------------------------
            ! 【目的】sigma2 を固定して Q2 になる sigma1 を求める (走査 + 二分法).
            implicit none
            real(dp), intent(in) :: R, Q2, sigma2, sigma3, alpha2, alpha3
            type(qs3_shape) :: s, s_prev, s_mid
            real(dp) :: x, x_prev, lo, hi, f_lo, f_mid
            integer :: i, it, n_scan

            n_scan = nint(sigma1_max / sigma1_step)
            s_prev%valid = .false.
            x_prev = 0.0_dp
            do i = 0, n_scan
                x = sigma1_step * i
                s = qs3_from_sigma(R, x, sigma2, sigma3, alpha2, alpha3)
                if (s%valid .and. s_prev%valid) then
                    if ((s_prev%q2_barn() - Q2) * (s%q2_barn() - Q2) <= 0.0_dp) then
                        lo = x_prev
                        hi = x
                        f_lo = s_prev%q2_barn() - Q2
                        do it = 1, 60
                            s_mid = qs3_from_sigma(R, 0.5_dp * (lo + hi), sigma2, sigma3, alpha2, alpha3)
                            if (.not. s_mid%valid) exit
                            f_mid = s_mid%q2_barn() - Q2
                            if (f_mid * f_lo <= 0.0_dp) then
                                hi = 0.5_dp * (lo + hi)
                            else
                                lo = 0.5_dp * (lo + hi)
                                f_lo = f_mid
                            end if
                            if (hi - lo < 1.0e-12_dp) exit
                        end do
                        s = qs3_from_sigma(R, 0.5_dp * (lo + hi), sigma2, sigma3, alpha2, alpha3)
                        return
                    end if
                end if
                s_prev = s
                x_prev = x
            end do
            s%valid = .false.
        end function shape_at_q2

        function find_neck_branch(R, Q2, eps_f1, eps_f2, alpha_g) result(br)
            !---------------------------
            ! increasing branch of d(sigma2) at fixed Q2, eps_f1, eps_f2, alpha_g:
            !   1. sigma2 = sinh(y), y equally spaced from asinh(sigma2_lowest) to
            !      asinh(min(a1^2/c1^2, a2^2/c2^2)) (dense around sigma2 = 0);
            !   2. the largest d over the scan, refined by golden-section search -> d_max;
            !   3. going down in sigma2 from there while the shape exists and d decreases;
            !      if the scan ends at a non-existing shape the boundary is refined by bisection
            !      -> d_min.
            !---------------------------
            ! 【目的】Q2, eps_f1, eps_f2, alpha_g を固定したときの d(sigma2) の増加する枝
            ! [d_min, d_max] と, その両端の sigma2 を求める.
            implicit none
            real(dp), intent(in) :: R, Q2, eps_f1, eps_f2, alpha_g
            type(neck_branch) :: br
            real(dp) :: y_lo, y_hi, sg(n_sig2_scan), dd(n_sig2_scan), sig_hi
            real(dp) :: a, b, c1g, c2g, f1, f2, lo, hi
            logical :: ok(n_sig2_scan)
            type(qs3_shape) :: s
            real(dp), parameter :: gr = 0.6180339887498949_dp
            integer :: i, ipk, ilo, it

            br%ok = .false.
            br%R = R
            br%Q2 = Q2
            call moller_to_sigma(eps_f1, eps_f2, alpha_g, br%sigma3, br%alpha3, br%alpha2)
            sig_hi = min(br%sigma3 + 0.5_dp * br%alpha3, br%sigma3 - 0.5_dp * br%alpha3)
            y_lo = asinh(sigma2_lowest)
            y_hi = asinh(sig_hi)
            do i = 1, n_sig2_scan
                sg(i) = sinh(y_lo + (y_hi - y_lo) * (i - 1) / real(n_sig2_scan, dp))
                ok(i) = .false.
                dd(i) = -1.0_dp
                s = shape_at_q2(R, Q2, sg(i), br%sigma3, br%alpha2, br%alpha3)
                if (s%valid) then
                    ok(i) = .true.
                    dd(i) = s%neck_diameter()
                end if
            end do
            if (.not. any(ok)) return
            ipk = maxloc(dd, dim=1)

            ! golden-section refinement of the maximum between the neighbouring samples
            a = sg(max(1, ipk - 1))
            b = sg(min(n_sig2_scan, ipk + 1))
            c1g = b - gr * (b - a)
            c2g = a + gr * (b - a)
            f1 = neck_at(c1g)
            f2 = neck_at(c2g)
            do it = 1, 40
                if (f1 > f2) then
                    b = c2g
                    c2g = c1g
                    f2 = f1
                    c1g = b - gr * (b - a)
                    f1 = neck_at(c1g)
                else
                    a = c1g
                    c1g = c2g
                    f1 = f2
                    c2g = a + gr * (b - a)
                    f2 = neck_at(c2g)
                end if
                if (b - a < 1.0e-8_dp * max(1.0_dp, abs(a))) exit
            end do
            if (max(f1, f2) > dd(ipk)) then
                if (f1 > f2) then
                    br%sig2_max = c1g
                    br%d_max = f1
                else
                    br%sig2_max = c2g
                    br%d_max = f2
                end if
            else
                br%sig2_max = sg(ipk)
                br%d_max = dd(ipk)
            end if

            ! go down along the increasing branch
            ilo = ipk
            do while (ilo > 1)
                if (.not. ok(ilo - 1)) exit
                if (dd(ilo - 1) > dd(ilo)) exit
                ilo = ilo - 1
            end do
            br%sig2_min = sg(ilo)
            br%d_min = dd(ilo)
            if (ilo > 1) then
                if (.not. ok(ilo - 1)) then
                    ! boundary of existence (scission or Q2 out of reach) by bisection
                    lo = sg(ilo - 1)
                    hi = sg(ilo)
                    do it = 1, 50
                        s = shape_at_q2(R, Q2, 0.5_dp * (lo + hi), br%sigma3, br%alpha2, br%alpha3)
                        if (s%valid) then
                            hi = 0.5_dp * (lo + hi)
                            br%sig2_min = hi
                            br%d_min = s%neck_diameter()
                        else
                            lo = 0.5_dp * (lo + hi)
                        end if
                        if (hi - lo < 1.0e-10_dp * max(1.0_dp, abs(hi))) exit
                    end do
                end if
            end if
            br%ok = .true.

            contains

                function neck_at(x) result(dn)
                    real(dp), intent(in) :: x
                    real(dp) :: dn
                    type(qs3_shape) :: ss
                    ss = shape_at_q2(R, Q2, x, br%sigma3, br%alpha2, br%alpha3)
                    dn = -1.0_dp
                    if (ss%valid) dn = ss%neck_diameter()
                end function neck_at
        end function find_neck_branch

        function shape_at_neck(br, d) result(s)
            !---------------------------
            ! shape on the branch br with the neck diameter d [fm], d_min <= d <= d_max
            ! (bisection in sigma2; d increases with sigma2 on the branch)
            !---------------------------
            ! 【目的】枝 br の上で首の直径が d [fm] になる形 (sigma2 の二分法).
            implicit none
            type(neck_branch), intent(in) :: br
            real(dp), intent(in) :: d
            type(qs3_shape) :: s, sm
            real(dp) :: lo, hi, mid
            integer :: it

            s%valid = .false.
            if (.not. br%ok) return
            if (d <= br%d_min) then
                s = shape_at_q2(br%R, br%Q2, br%sig2_min, br%sigma3, br%alpha2, br%alpha3)
                return
            end if
            if (d >= br%d_max) then
                s = shape_at_q2(br%R, br%Q2, br%sig2_max, br%sigma3, br%alpha2, br%alpha3)
                return
            end if
            lo = br%sig2_min
            hi = br%sig2_max
            do it = 1, 60
                mid = 0.5_dp * (lo + hi)
                sm = shape_at_q2(br%R, br%Q2, mid, br%sigma3, br%alpha2, br%alpha3)
                if (.not. sm%valid) then
                    lo = mid
                    cycle
                end if
                if (sm%neck_diameter() < d) then
                    lo = mid
                else
                    hi = mid
                end if
                if (hi - lo < 1.0e-11_dp * max(1.0_dp, abs(mid))) exit
            end do
            s = shape_at_q2(br%R, br%Q2, 0.5_dp * (lo + hi), br%sigma3, br%alpha2, br%alpha3)
        end function shape_at_neck

end module qs3_coord_mod
