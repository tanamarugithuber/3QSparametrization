program test_qs3
    !---------------------------
    ! checks of the 3QS shape (qs3_shape_mod) and of Moller's coordinates (qs3_coord_mod)
    !   1. volume, centre of mass, continuity of rho and d rho/dz at z1, z2
    !   2. eps_f1, eps_f2, alpha_g of the constructed shape equal the input (Eqs. 6, 7, 10, 11)
    !   3. Q2 and d of shape_at_q2 / shape_at_neck equal the targets
    !   4. volume and area by the generic surface integrals (eval) agree with the exact ones
    !   5. mirror: (eps_f1, eps_f2, alpha_g) and (eps_f2, eps_f1, -alpha_g) give the same Q2, d, B_s
    !---------------------------
    use iso_fortran_env, only: real64
    use constant_mod, only: pi
    use shape_mod, only: surface_shape
    use qs3_shape_mod
    use qs3_coord_mod
    implicit none
    integer, parameter :: dp = real64
    real(dp), parameter :: R = 1.16_dp * 240.0_dp**(1.0_dp / 3.0_dp)
    type(qs3_shape) :: s, sm
    type(neck_branch) :: br, brm
    real(dp) :: s3, a3, a2, zc, h, e1, e2, ag, d, vol_gen, area_gen
    integer :: n_fail, j
    real(dp) :: tests(5, 3)

    n_fail = 0
    tests(1, :) = [0.0_dp, 0.0_dp, 0.0_dp]
    tests(2, :) = [0.3_dp, 0.0_dp, 0.2_dp]
    tests(3, :) = [-0.2_dp, 0.5_dp, 0.1_dp]
    tests(4, :) = [0.175_dp, 0.4_dp, -0.3_dp]
    tests(5, :) = [0.1_dp, 0.1_dp, 0.66_dp]

    print '(a)', " --- 1, 2: geometry at Q2 = 50 b, sigma2 = -0.5"
    do j = 1, 5
        call moller_to_sigma(tests(j, 1), tests(j, 2), tests(j, 3), s3, a3, a2)
        s = shape_at_q2(R, 50.0_dp, -0.5_dp, s3, a2, a3)
        if (.not. s%valid) then
            print '(a, 3f8.3)', "   no shape for ", tests(j, :)
            n_fail = n_fail + 1
            cycle
        end if
        zc = pi * centre(s) / s%volume()
        h = 1.0e-7_dp * R
        call check("volume / (4/3 pi R^3) - 1", s%volume() / (4.0_dp / 3.0_dp * pi * R**3) - 1.0_dp, 1.0e-12_dp)
        call check("centre of mass / R", zc / R, 1.0e-8_dp)
        call check("rho^2 jump at z1 / R^2", (s%rho2(s%zb(1) - h) - s%rho2(s%zb(1) + h)) / R**2, 1.0e-6_dp)
        call check("rho^2 jump at z2 / R^2", (s%rho2(s%zb(2) - h) - s%rho2(s%zb(2) + h)) / R**2, 1.0e-6_dp)
        call check("slope jump at z1", slope(s, s%zb(1) - h) - slope(s, s%zb(1) + h), 1.0e-5_dp)
        call check("slope jump at z2", slope(s, s%zb(2) - h) - slope(s, s%zb(2) + h), 1.0e-5_dp)
        e1 = s%eps_f(1)
        e2 = s%eps_f(2)
        ag = s%alpha_g()
        call check("eps_f1 - input", e1 - tests(j, 1), 1.0e-12_dp)
        call check("eps_f2 - input", e2 - tests(j, 2), 1.0e-12_dp)
        call check("alpha_g - input", ag - tests(j, 3), 1.0e-12_dp)
        call check("Q2 - 50 b", s%q2_barn() - 50.0_dp, 1.0e-8_dp)
        vol_gen = generic_volume(s)
        call check("volume (eval, Gauss) / exact - 1", vol_gen / s%volume() - 1.0_dp, 1.0e-6_dp)
        area_gen = generic_area(s) / (4.0_dp * pi * R**2)
        call check("B_s (eval, Gauss) - B_s", area_gen - s%relative_area(), 1.0e-5_dp)
        print '(a, 3f7.3, a, 2f8.4, a, f8.4, a, f8.3, a, f8.5, a, f7.4)', "   eps_f1, eps_f2, alpha_g =", tests(j, :), &
              "  sigma1, alpha1 =", s%sigma1, s%alpha1, "  d/R =", s%neck_diameter() / R, &
              "  Q2 =", s%q2_barn(), "  B_s =", s%relative_area(), "  B_W =", s%b_w(0.9_dp)
    end do

    print '(a)', " --- 3, 5: neck branch and mirror symmetry"
    do j = 2, 4
        br = find_neck_branch(R, 60.0_dp, tests(j, 1), tests(j, 2), tests(j, 3))
        brm = find_neck_branch(R, 60.0_dp, tests(j, 2), tests(j, 1), -tests(j, 3))
        if (.not. br%ok .or. .not. brm%ok) then
            print *, "  no neck branch"
            n_fail = n_fail + 1
            cycle
        end if
        print '(a, 3f7.3, a, 2f9.4, a, 2f11.4)', "   eps_f1, eps_f2, alpha_g =", tests(j, :), &
              "  d_min/R, d_max/R =", br%d_min / R, br%d_max / R, "  sigma2 range", br%sig2_min, br%sig2_max
        call check("d_min mirror", (br%d_min - brm%d_min) / R, 1.0e-6_dp)
        call check("d_max mirror", (br%d_max - brm%d_max) / R, 1.0e-6_dp)
        d = br%d_min + 0.37_dp * (br%d_max - br%d_min)
        s = shape_at_neck(br, d)
        sm = shape_at_neck(brm, d)
        call check("d - target", (s%neck_diameter() - d) / R, 1.0e-8_dp)
        call check("Q2 - 60 b", s%q2_barn() - 60.0_dp, 1.0e-6_dp)
        call check("B_s mirror", s%relative_area() - sm%relative_area(), 1.0e-9_dp)
        call check("sigma1 mirror", s%sigma1 - sm%sigma1, 1.0e-7_dp)
    end do

    if (n_fail == 0) then
        print '(a)', " all tests passed"
    else
        print '(a, i0, a)', " ", n_fail, " test(s) FAILED"
        error stop 1
    end if

contains

    subroutine check(name, val, tol)
        character(*), intent(in) :: name
        real(dp), intent(in) :: val, tol
        if (abs(val) > tol .or. val /= val) then
            print '(a, a, es12.3, a, es9.1)', "   FAIL ", name, val, "  tol", tol
            n_fail = n_fail + 1
        end if
    end subroutine check

    function centre(s) result(m)
        ! int z rho^2 dz by fine Gauss-Legendre per piece (independent of the module routine)
        type(qs3_shape), intent(in) :: s
        real(dp) :: m, z
        integer :: i, k
        integer, parameter :: n = 20000
        m = 0.0_dp
        do i = 1, 3
            do k = 1, n
                z = s%zb(i - 1) + (s%zb(i) - s%zb(i - 1)) * (k - 0.5_dp) / n
                m = m + (s%zb(i) - s%zb(i - 1)) / n * z * max(0.0_dp, s%p(0, i) + s%p(1, i) * z + s%p(2, i) * z * z)
            end do
        end do
    end function centre

    function slope(s, z) result(sl)
        ! d rho^2/dz from the piece containing z
        type(qs3_shape), intent(in) :: s
        real(dp), intent(in) :: z
        real(dp) :: sl
        integer :: i
        i = 2
        if (z < s%zb(1)) i = 1
        if (z > s%zb(2)) i = 3
        sl = (s%p(1, i) + 2.0_dp * s%p(2, i) * z) / R
    end function slope

    function generic_volume(sh) result(v)
        class(surface_shape), intent(in) :: sh
        real(dp) :: v
        v = gen_int(sh, 1)
    end function generic_volume

    function generic_area(sh) result(a)
        class(surface_shape), intent(in) :: sh
        real(dp) :: a
        a = gen_int(sh, 2)
    end function generic_area

    function gen_int(sh, what) result(val)
        ! midpoint rule in u (fine), exact in v by axial symmetry
        class(surface_shape), intent(in) :: sh
        integer, intent(in) :: what
        real(dp) :: val, r(3), N(3), u
        integer :: i
        integer, parameter :: nu = 40000
        val = 0.0_dp
        do i = 1, nu
            u = pi * (i - 0.5_dp) / nu
            call sh%eval(u, 0.0_dp, r, N)
            if (what == 1) then
                val = val + dot_product(r, N) / 3.0_dp
            else
                val = val + norm2(N)
            end if
        end do
        val = val * (pi / nu) * 2.0_dp * pi
    end function gen_int

end program test_qs3
