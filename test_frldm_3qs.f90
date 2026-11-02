module wrapped_shape_mod
    ! a 3QS shape seen only through the generic surface_shape interface
    ! (forces the general double surface integral of frldm_mod, for comparison)
    use iso_fortran_env, only: real64
    use shape_mod, only: surface_shape
    use qs3_shape_mod, only: qs3_shape
    implicit none
    integer, parameter, private :: dp = real64
    type, extends(surface_shape) :: wrapped_shape
        type(qs3_shape) :: s
        contains
            procedure :: eval => w_eval
            procedure :: inside => w_inside
    end type wrapped_shape
contains
    pure subroutine w_eval(this, u, v, r, N)
        class(wrapped_shape), intent(in) :: this
        real(dp), intent(in) :: u, v
        real(dp), intent(out) :: r(3), N(3)
        call this%s%eval(u, v, r, N)
    end subroutine w_eval
    pure function w_inside(this, x) result(yes)
        class(wrapped_shape), intent(in) :: this
        real(dp), intent(in) :: x(3)
        logical :: yes
        yes = this%s%inside(x)
    end function w_inside
end module wrapped_shape_mod

program test_frldm_3qs
    !---------------------------
    ! B_1, B_3 of axially symmetric shapes by ring sums (frldm_mod::axial_sums)
    !   1. sphere: rings vs the analytic B_1^(0), B_3^(0) of FRDM2012 Eq. (74)
    !   2. 3QS shapes: axial sums (default for qs3_shape) vs the general double surface integral
    !      (the same shape wrapped as a plain surface_shape), and convergence in n
    !   3. E_mac(3QS) - E_mac(sphere) and B_W along a symmetric path towards scission (240Pu)
    !---------------------------
    use iso_fortran_env, only: real64
    use constant_mod, only: pi
    use nucleus_mod, only: nucleus_property
    use shape_mod, only: epsilon_shape_from, epsilon_shape, gauss_legendre
    use qs3_shape_mod
    use qs3_coord_mod
    use frldm_mod
    use wrapped_shape_mod
    implicit none
    integer, parameter :: dp = real64
    type(nucleus_property) :: nuc
    type(frldm_variables) :: fr
    type(qs3_shape) :: s
    type(wrapped_shape) :: ws
    type(epsilon_shape) :: sph
    type(neck_branch) :: br
    real(dp) :: R, A, x0, y0, B1_0, B3_0, B1_ax, B3_ax, B1_gen, B3_gen, E_sph, E, t0, t1
    real(dp) :: s3, a3, a2
    integer :: n_fail, j, n

    n_fail = 0
    nuc%Z = 94
    nuc%N = 146
    nuc%A = 240
    A = 240.0_dp
    nuc%R_a = nuc%R0 * A**(1.0_dp / 3.0_dp)
    R = nuc%R_a
    fr%verbose = .false.

    !--------------------------- 1. sphere
    x0 = R / a_Yukawa
    y0 = R / a_den
    B1_0 = 1.0_dp - 3.0_dp / x0**2 + (1.0_dp + x0) * (2.0_dp + 3.0_dp / x0 + 3.0_dp / x0**2) * exp(-2.0_dp * x0)
    B3_0 = 1.0_dp - 5.0_dp / y0**2 * (1.0_dp - 15.0_dp / (8.0_dp * y0) + 21.0_dp / (8.0_dp * y0**3) &
         - 0.75_dp * (1.0_dp + 4.5_dp / y0 + 7.0_dp / y0**2 + 3.5_dp / y0**3) * exp(-2.0_dp * y0))
    call sphere_rings(64, 96, B1_ax, B3_ax)
    print '(a, 2f12.8, a, 2es10.2)', " sphere, rings     B_1, B_3 =", B1_ax, B3_ax, "   - Eq. (74):", &
          B1_ax - B1_0, B3_ax - B3_0
    call check("sphere B_1", B1_ax - B1_0, 1.0e-4_dp)
    call check("sphere B_3", B3_ax - B3_0, 1.0e-8_dp)

    !--------------------------- 2. 3QS: axial vs general
    sph = epsilon_shape_from(R, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp)
    call fr%calculate_frldm_energy(nuc, sph, E_sph)
    call moller_to_sigma(0.2_dp, 0.0_dp, 0.15_dp, s3, a3, a2)
    s = shape_at_q2(R, 60.0_dp, -0.6_dp, s3, a2, a3)
    if (.not. s%valid) error stop "no shape"
    ws%s = s
    print '(a)', " 3QS (eps_f1, eps_f2, alpha_g) = (0.2, 0, 0.15), Q2 = 60 b, sigma2 = -0.6"
    do n = 32, 128, 32
        call cpu_time(t0)
        call fr%calculate_b_1_b_3_surface(nuc, s, n)
        call cpu_time(t1)
        if (n == 64) then
            B1_gen = fr%B_1
            B3_gen = fr%B_3
        end if
        B1_ax = fr%B_1
        B3_ax = fr%B_3
        print '(a, i4, a, 2f13.8, a, f8.3, a)', "   n =", n, "  rings  : B_1, B_3 =", B1_ax, B3_ax, &
              "   (", t1 - t0, " s cpu)"
    end do
    ! n = 64 (the default) against n = 128
    call check("B_1(n = 64) - B_1(n = 128)", B1_gen - B1_ax, 5.0e-5_dp)
    call check("B_3(n = 64) - B_3(n = 128)", B3_gen - B3_ax, 1.0e-7_dp)
    ! the general double surface integral converges slowly here (rho'' jumps at z1, z2)
    call cpu_time(t0)
    call fr%calculate_b_1_b_3_surface(nuc, ws, 64)
    call cpu_time(t1)
    print '(a, i4, a, 2f13.8, a, f8.3, a)', "   n =", 64, "  general: B_1, B_3 =", fr%B_1, fr%B_3, &
          "   (", t1 - t0, " s cpu)"
    call check("B_1 general(n = 64) - rings", fr%B_1 - B1_ax, 2.0e-2_dp)
    call check("B_3 general(n = 64) - rings", fr%B_3 - B3_ax, 2.0e-3_dp)

    !--------------------------- 3. symmetric path towards scission, 240Pu
    print '(a)', " 240Pu, eps_f1 = eps_f2 = 0, alpha_g = 0, Q2 = 80 b: E_mac - E_mac(sphere) on the neck branch"
    br = find_neck_branch(R, 80.0_dp, 0.0_dp, 0.0_dp, 0.0_dp)
    print '(a, 2f9.4)', "   d_min/R, d_max/R:", br%d_min / R, br%d_max / R
    do j = 0, 4
        s = shape_at_neck(br, br%d_min + (br%d_max - br%d_min) * j / 4.0_dp)
        call fr%calculate_frldm_energy(nuc, s, E)
        print '(a, f8.4, a, f9.4, a, f8.4, a, f9.4, a, f7.4, a, f10.4)', "   d/R =", s%neck_diameter() / R, &
              "  sigma2 =", s%sigma2, "  sigma1 =", s%sigma1, "  B_s =", fr%B_s, "  B_W =", fr%B_w, &
              "  E - E_sph =", E - E_sph
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

    subroutine sphere_rings(n1, n2, B1, B3)
        ! rings on a sphere (Gauss-Legendre in theta), Richardson as in frldm_mod
        integer, intent(in) :: n1, n2
        real(dp), intent(out) :: B1, B3
        real(dp) :: s1a, s3a, s1b, s3b
        call sums(n1, s1a, s3a)
        call sums(n2, s1b, s3b)
        s1a = (real(n2, dp)**3 * s1b - real(n1, dp)**3 * s1a) / (real(n2, dp)**3 - real(n1, dp)**3)
        B1 = a_Yukawa**2 * s1a / (8.0_dp * pi**2 * nuc%R0**2 * a_Yukawa**4) * A**(-2.0_dp / 3.0_dp)
        B3 = -a_den * s3b * 15.0_dp / (32.0_dp * pi**2 * nuc%R0**5) * A**(-5.0_dp / 3.0_dp)
    end subroutine sphere_rings

    subroutine sums(n, S1, S3)
        integer, intent(in) :: n
        real(dp), intent(out) :: S1, S3
        real(dp) :: x(n), w(n), th(n), z(n), rho(n), pz(n), wr(n)
        call gauss_legendre(n, x, w)
        th = 0.5_dp * pi * (x + 1.0_dp)
        z = -R * cos(th)
        rho = R * sin(th)
        pz = z
        wr = R * sin(th) * 0.5_dp * pi * w
        call axial_sums(z, rho, pz, wr, 2 * n, S1, S3)
    end subroutine sums

end program test_frldm_3qs
