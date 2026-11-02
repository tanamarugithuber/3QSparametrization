module shell_bcs_mod
    !---------------------------
    ! Microscopic corrections: Strutinsky shell correction + pairing correction
    ! (Lipkin-Nogami as in FRDM2012, or BCS)
    !   P. Moller et al., ADNDT 109-110 (2016) 1  (FRDM2012), Secs. 2.7-2.10
    !   P. Moller and J. R. Nix, Nucl. Phys. A 536 (1992) 20  (thesis/Nuclear_Pairing_Model.pdf)
    !
    ! 1. single-particle levels
    !    V = V1 + V_so + V_C (Eqs. 80-92) for the shape (eps2, gamma, eps4, eps6) of radius
    !    R_pot (shape_mod), V1 and V_C on the Gauss-Hermite points by surface integrals
    !    (potential_gh_mod), spin-orbit term from V1 by integration by parts, diagonalised
    !    in the deformed HO basis with hw0 = C_cur / A^(1/3) (ho_basis_mod).
    ! 2. Strutinsky shell correction (Eq. 78)
    !    E_shell = sum_{i <= N} e_i - E_smooth,
    !    smoothing function: Gaussian times the correction polynomial of order p (Hermite form)
    !      f(x) = exp(-x^2)/sqrt(pi) * sum_{m=0}^{p/2} c_m H_2m(x),  c_m = (-1)^m / (4^m m!)
    !    range gamma_S = C_S hw0 B_s (Eq. 111), smooth Fermi level from N_smooth(lambda) = N.
    ! 3. pairing correction (Eq. 79)
    !    E_pair = E_pc - E_pc_average,
    !    E_pc from the Lipkin-Nogami (Eqs. 94-101) or BCS (Moller-Nix Eq. 5) equations
    !    between levels N1 and N2,
    !    G from the effective-interaction gap Delta_G = r_mic B_s / N^(1/3) (Eqs. 102-107),
    !    E_pc_average (Eq. 108, or Moller-Nix Eq. 27 for BCS) for the same G, Delta_G and
    !    level density; for odd particle numbers it contains the average gap
    !    r_mac B_s / N^(1/3) of the macroscopic energy (Eqs. 63, 64).
    !    Pairing window (Moller-Nix Sec. 2.5): N2 = last bound level or last level within
    !    5 MeV above the Fermi level, whichever is higher; the same number of levels below.
    !    Odd particle numbers: no blocking, quasi-particle energy of the odd level added.
    !
    ! Each eigenvalue of the 2 n_states x 2 n_states Hamiltonian holds one nucleon;
    ! Kramers pairs (time-reversal partners) form the doubly degenerate levels of BCS.
    ! Constants are taken from micro_constant_mod.
    !---------------------------
    ! 【日本語での説明】
    ! 微視的補正 E_s+p = E_shell + E_pair (FRDM2012 2.7-2.10 節) を計算する.
    !   E_shell : Strutinsky 殻補正 = (占有準位のエネルギー和) - (平滑化したエネルギー)   (式 78)
    !   E_pair  : 対補正 = E_pc - <E_pc>  (Lipkin-Nogami 式 94-110, または BCS)          (式 79)
    ! 【物理量】
    !   e(:)         : 一粒子準位 [MeV] (昇順, 1 準位に核子 1 個; 時間反転の対が縮退)
    !   gamma_S      : Strutinsky の平滑化幅 C_S hbar omega0 B_s [MeV]              (式 111)
    !   p            : 補正多項式の次数 (8)
    !   lambda_smooth: 平滑化した準位分布のフェルミ準位 [MeV]
    !   g_smooth     : そこでの平滑化準位密度 (スピン両方) [1/MeV]
    !   rho          : 二重縮退準位の平均密度 g_smooth / 2 [1/MeV]
    !   Delta_G      : 有効相互作用の対ギャップ r_mic B_s / N^(1/3) [MeV]           (式 102)
    !   G            : 対相互作用の強さ [MeV]                                       (式 107)
    !   Delta, lambda: 対ギャップとフェルミ準位 (BCS/LN の解) [MeV]
    !   lambda2      : Lipkin-Nogami の粒子数ゆらぎ定数 [MeV]                       (式 98)
    !   N1, N2       : 対相関を考える準位の範囲 (二重縮退準位の番号)
    !   E_pc, <E_pc> : 対相関エネルギー + 準粒子エネルギーとその平均 [MeV]   (式 101, 108)
    use iso_fortran_env, only: real64
    use constant_mod, only: pi, e2
    use nucleus_mod, only: nucleus_property
    use ho_basis_mod, only: ho_basis_type, frequencies_from_eps_gamma, so_strength, eigenvalues_hamiltonian
    use potential_gh_mod, only: ellipsoid_shape, epsilon_shape, epsilon_shape_from, surface_quadrature, &
                                gauss_legendre
    use micro_constant_mod, only: microscopic_variables, a_pot, N_bas, mi_p, C_s, r_mic, mi_r_mac
    use shape_mod, only: surface_shape
    use qs3_shape_mod, only: qs3_shape
    implicit none
    private

    integer, parameter :: dp = real64

    !---------------------------
    ! results for one kind of nucleon
    !---------------------------
    ! 中性子または陽子 1 種類についての結果 (各成分の意味は上の物理量の一覧を参照)
    type, public :: nucleon_correction
        real(dp) :: e_sum = 0.0_dp          ! MeV, sum of the occupied single-particle energies
        real(dp) :: e_smooth = 0.0_dp       ! MeV, Strutinsky smooth energy
        real(dp) :: lambda_smooth = 0.0_dp  ! MeV, Fermi level of the smooth distribution
        real(dp) :: g_smooth = 0.0_dp       ! 1/MeV, smooth level density at lambda_smooth (both spins)
        real(dp) :: e_shell = 0.0_dp        ! MeV, shell correction
        real(dp) :: delta_G = 0.0_dp        ! MeV, effective-interaction pairing gap
        real(dp) :: G = 0.0_dp              ! MeV, pairing strength
        real(dp) :: delta = 0.0_dp          ! MeV, pairing gap (0: no non-trivial solution)
        real(dp) :: lambda = 0.0_dp         ! MeV, Fermi level of the pairing solution
        real(dp) :: lambda2 = 0.0_dp        ! MeV, Lipkin-Nogami number-fluctuation constant
        integer :: N1 = 0, N2 = 0           ! pairing window (doubly degenerate levels)
        real(dp) :: e_pc = 0.0_dp           ! MeV, pairing correlation + quasi-particle energy
        real(dp) :: e_pc_ave = 0.0_dp       ! MeV, its average
        real(dp) :: e_pair = 0.0_dp         ! MeV, pairing correction
        real(dp), allocatable :: e(:)       ! MeV, single-particle energies (ascending, one nucleon each)
    end type nucleon_correction

    ! 原子核全体の結果: e_micro = e_shell + e_pair が main で E_mac に足す微視的補正
    type, public :: shell_bcs_result
        type(nucleon_correction) :: neutron, proton
        real(dp) :: B_s = 1.0_dp            ! relative surface area of the shape
        real(dp) :: gamma_S = 0.0_dp        ! MeV, Strutinsky range
        real(dp) :: e_shell = 0.0_dp        ! MeV, neutron + proton shell correction
        real(dp) :: e_pair = 0.0_dp         ! MeV, neutron + proton pairing correction
        real(dp) :: e_micro = 0.0_dp        ! MeV, e_shell + e_pair
    end type shell_bcs_result

    public :: shell_bcs_calculation
    public :: strutinsky_correction
    public :: bcs_correction
    public :: ln_correction
    public :: ellipsoid_relative_area

    contains

        subroutine shell_bcs_calculation(nucleus, mv, res, n_gh, pairing, qs3, n_basis)
            !---------------------------
            ! shell + pairing correction of the nucleus at the deformation
            ! (nucleus%epsilon2, nucleus%gamma_deg, nucleus%epsilon3, nucleus%epsilon4, nucleus%epsilon6).
            ! mv receives the potential quantities (micro_constant_mod).
            ! n_gh   : Gauss-Hermite points per axis (default 44)
            ! pairing: 'LN' (Lipkin-Nogami, default, as FRDM2012) or 'BCS'
            ! qs3    : if present, the 3QS shape is used instead of the epsilon shape (it is
            !          rescaled to the radius R_pot); the HO basis is that of the spheroid with
            !          the same volume and quadrupole moment q = int (2z^2 - rho^2) dV/(V R^2)
            !          = 2/5 (c^2 - a^2)/R^2, a^2 c = R^3.
            ! n_basis: number of oscillator shells (default N_bas = 12 of micro_constant_mod)
            !---------------------------
            ! 【3QS】qs3 を渡すと, その形 (半径 R_pot に拡大縮小) のポテンシャルを使う.
            ! 振動子基底は同じ体積と四重極モーメントの回転楕円体 (半軸 a, c) の振動数
            ! hbar omega_z = hbar omega0 R/c, hbar omega_perp = hbar omega0 R/a にする.
            ! 細長い形では n_basis (殻の数) を増やす必要がある.
            ! 【日本語】
            ! 【目的】変形 (eps2, gamma, eps3, eps4, eps6) での殻補正と対補正を
            ! 中性子・陽子それぞれ計算する.
            ! 【手順】
            !   1. calculate_microscopic_variables で R_pot, V_p, V_n, lambda_p, lambda_n などを求める.
            !   2. 振動子基底を作り, 半径 R_pot の epsilon 形状で Y, C を Gauss-Hermite 点上に求める.
            !   3. 中性子 h = T - V_n Y + V_so,  陽子 h = T - V_p Y + e rho_c C + V_so を対角化する.
            !   4. Strutinsky 補正 (strutinsky_correction) と対補正 (ln_correction / bcs_correction).
            implicit none
            type(nucleus_property), intent(in) :: nucleus
            type(microscopic_variables), intent(inout) :: mv
            type(shell_bcs_result), intent(out) :: res
            integer, intent(in), optional :: n_gh
            character(*), intent(in), optional :: pairing
            type(qs3_shape), intent(in), optional :: qs3
            integer, intent(in), optional :: n_basis

            real(dp) :: eps, gam, hw(3), cr
            integer :: ngh, ns, nb
            logical :: use_ln, d2h
            type(ho_basis_type) :: basis
            class(surface_shape), allocatable :: shape
            type(surface_quadrature) :: q
            real(dp), allocatable :: xg(:), yg(:), zg(:)
            real(dp), allocatable :: Y(:,:,:), dYx(:,:,:), dYy(:,:,:), dYz(:,:,:), C(:,:,:)
            real(dp), allocatable :: T(:,:), Vm(:,:), Axy(:,:), Ayz(:,:), Azx(:,:)
            complex(dp), allocatable :: H(:,:)

            ngh = 44
            if (present(n_gh)) ngh = n_gh
            use_ln = .true.
            if (present(pairing)) use_ln = (pairing /= 'BCS' .and. pairing /= 'bcs')

            eps = nucleus%epsilon2
            gam = nucleus%gamma_deg * pi / 180.0_dp

            ! potential quantities, Eqs. (82)-(87), (89), (90), (92), (93)
            call mv%calculate_microscopic_variables(nucleus)

            !---------------------------
            ! basis, shape and potentials on the Gauss-Hermite points
            !---------------------------
            nb = nint(N_bas)
            if (present(n_basis)) nb = n_basis
            if (present(qs3)) then
                cr = spheroid_c_over_R(qs3%q2_hat())
                hw = mv%hbar_omega0 * [sqrt(cr), sqrt(cr), 1.0_dp / cr]   ! R/a = (c/R)^(1/2), R/c
                allocate(shape, source=qs3%scaled(mv%R_pot))
                d2h = qs3%reflection_symmetric()
            else
                hw = frequencies_from_eps_gamma(mv%hbar_omega0, eps, gam)
                allocate(shape, source=epsilon_shape_from(mv%R_pot, eps, gam, nucleus%epsilon3, &
                                                          nucleus%epsilon4, nucleus%epsilon6))
                d2h = nucleus%epsilon3 == 0.0_dp
            end if
            call basis%initialize(hw, nb, ngh)
            ns = basis%n_states
            call q%initialize(shape, a_pot, 8, 16, 8)
            xg = basis%b(1) * basis%xi
            yg = basis%b(2) * basis%xi
            zg = basis%b(3) * basis%xi
            allocate(Y(ngh, ngh, ngh), dYx(ngh, ngh, ngh), dYy(ngh, ngh, ngh), dYz(ngh, ngh, ngh), C(ngh, ngh, ngh))
            ! the octant filling needs the z reflection, which eps3 /= 0 (or a 3QS mass asymmetry) breaks
            call q%on_grid(xg, yg, zg, Y, dYx, dYy, dYz, C, d2h)

            !---------------------------
            ! single-particle spectra
            !---------------------------
            allocate(T(ns, ns), Vm(ns, ns), Axy(ns, ns), Ayz(ns, ns), Azx(ns, ns))
            allocate(H(2 * ns, 2 * ns))
            allocate(res%neutron%e(2 * ns), res%proton%e(2 * ns))
            call basis%kinetic_matrix(T)

            call basis%potential_matrix(-mv%V_n * Y, Vm)
            call basis%spin_orbit_matrices_by_parts(-mv%V_n * Y, Axy, Ayz, Azx)
            call basis%build_hamiltonian(T, Vm, Axy, Ayz, Azx, so_strength(mv%lambda_n), H)
            call eigenvalues_hamiltonian(H, res%neutron%e)

            call basis%potential_matrix(-mv%V_p * Y + mv%e_rho_c * C, Vm)
            call basis%spin_orbit_matrices_by_parts(-mv%V_p * Y, Axy, Ayz, Azx)
            call basis%build_hamiltonian(T, Vm, Axy, Ayz, Azx, so_strength(mv%lambda_p), H)
            call eigenvalues_hamiltonian(H, res%proton%e)

            !---------------------------
            ! Strutinsky and pairing
            !---------------------------
            res%B_s = shape%relative_area()
            res%gamma_S = C_s * mv%hbar_omega0 * res%B_s
            call correct_one(res%neutron, nucleus%N)
            call correct_one(res%proton, nucleus%Z)
            res%e_shell = res%neutron%e_shell + res%proton%e_shell
            res%e_pair = res%neutron%e_pair + res%proton%e_pair
            res%e_micro = res%e_shell + res%e_pair

        contains

            ! 【目的】1 種類の核子 (n_part 個) について E_shell と E_pair を求める.
            ! delta_odd は奇数核子数のとき <E_pc> に足す平均ギャップ r_mac B_s / N^(1/3)
            ! (巨視的エネルギー式 63, 64 と同じ項なので, 両者で打ち消し合う).
            subroutine correct_one(nc, n_part)
                type(nucleon_correction), intent(inout) :: nc
                integer, intent(in) :: n_part
                real(dp) :: delta_odd
                call strutinsky_correction(nc%e, n_part, res%gamma_S, nint(mi_p), &
                                           nc%e_sum, nc%e_smooth, nc%lambda_smooth, nc%g_smooth)
                nc%e_shell = nc%e_sum - nc%e_smooth
                nc%delta_G = r_mic * res%B_s / real(n_part, dp)**(1.0_dp / 3.0_dp)
                delta_odd = mi_r_mac * res%B_s / real(n_part, dp)**(1.0_dp / 3.0_dp)
                if (use_ln) then
                    call ln_correction(nc%e, n_part, 0.5_dp * nc%g_smooth, nc%delta_G, nc%G, nc%delta, &
                                       nc%lambda, nc%lambda2, nc%N1, nc%N2, nc%e_pc, nc%e_pc_ave, delta_odd)
                else
                    call bcs_correction(nc%e, n_part, 0.5_dp * nc%g_smooth, nc%delta_G, nc%G, nc%delta, &
                                        nc%lambda, nc%N1, nc%N2, nc%e_pc, nc%e_pc_ave, delta_odd)
                end if
                nc%e_pair = nc%e_pc - nc%e_pc_ave
            end subroutine correct_one

        end subroutine shell_bcs_calculation

        pure function spheroid_c_over_R(q) result(cr)
            !---------------------------
            ! c/R of the spheroid (semi-axes a, a, c; a^2 c = R^3) with
            ! q = int (2 z^2 - rho^2) dV/(V R^2) = 2/5 (c^2 - a^2)/R^2:
            ! with x = c/R, 2/5 (x^2 - 1/x) = q, solved by Newton's method (unique root x > 0).
            !---------------------------
            ! 【目的】体積と四重極モーメント q が同じ回転楕円体の c/R (3QS の振動子基底用).
            real(dp), intent(in) :: q
            real(dp) :: cr, f, df
            integer :: it
            cr = max(1.0_dp, 1.0_dp + q)
            do it = 1, 100
                f = 0.4_dp * (cr * cr - 1.0_dp / cr) - q
                df = 0.4_dp * (2.0_dp * cr + 1.0_dp / (cr * cr))
                cr = cr - f / df
                if (abs(f) < 1.0e-14_dp) exit
            end do
        end function spheroid_c_over_R

        !==========================================================================
        ! Strutinsky shell correction
        !==========================================================================

        subroutine strutinsky_correction(e, n_part, gam, p, e_sum, e_smooth, lambda, g_lambda)
            !---------------------------
            ! e(:)     : single-particle energies, ascending, one nucleon each
            ! n_part   : number of nucleons
            ! gam      : smoothing range gamma_S (MeV),  p : order of the correction polynomial (even)
            ! e_sum    : sum of the lowest n_part energies
            ! e_smooth : smooth energy  int^lambda e g(e) de,  with  int^lambda g(e) de = n_part
            ! g_lambda : smooth level density at lambda
            !---------------------------
            ! 【日本語】
            ! 【目的】Strutinsky 殻補正の 2 つの項 e_sum = sum_{i<=N} e_i と e_smooth を求める.
            ! 【式】平滑化準位密度 g(e) = (1/gamma) sum_i f((e - e_i)/gamma),
            !   f(x) = e^(-x^2)/sqrt(pi) sum_{m=0}^{p/2} c_m H_2m(x)  (ガウス関数 x 補正多項式)
            !   フェルミ準位 lambda は int^lambda g(e) de = N を二分法で解いて求め,
            !   e_smooth = int^lambda e g(e) de とする. E_shell = e_sum - e_smooth.
            implicit none
            real(dp), intent(in) :: e(:), gam
            integer, intent(in) :: n_part, p
            real(dp), intent(out) :: e_sum, e_smooth, lambda, g_lambda
            real(dp) :: lo, hi, n_lo
            integer :: it

            if (n_part > size(e)) stop "strutinsky_correction: more nucleons than levels"
            e_sum = sum(e(1:n_part))

            ! bracket and bisect N_smooth(lambda) = n_part
            lo = e(1) - 10.0_dp * gam
            hi = e(n_part) + 10.0_dp * gam
            n_lo = smooth_number(e, lo, gam, p) - n_part
            if (n_lo > 0.0_dp .or. smooth_number(e, hi, gam, p) - n_part < 0.0_dp) &
                stop "strutinsky_correction: smooth Fermi level not bracketed"
            do it = 1, 200
                lambda = 0.5_dp * (lo + hi)
                if ((smooth_number(e, lambda, gam, p) - n_part) * n_lo > 0.0_dp) then
                    lo = lambda
                else
                    hi = lambda
                end if
                if (hi - lo < 1.0e-12_dp) exit
            end do
            lambda = 0.5_dp * (lo + hi)
            e_smooth = smooth_energy(e, lambda, gam, p)
            g_lambda = smooth_density(e, lambda, gam, p)
        end subroutine strutinsky_correction

        ! 【目的】エルミート多項式 H_0(x) ... H_nmax(x).  漸化式 H_{n+1} = 2x H_n - 2n H_{n-1}
        pure subroutine hermite(x, nmax, Hn)
            ! physicists' Hermite polynomials H_0 .. H_nmax at x
            implicit none
            real(dp), intent(in) :: x
            integer, intent(in) :: nmax
            real(dp), intent(out) :: Hn(0:nmax)
            integer :: n
            Hn(0) = 1.0_dp
            if (nmax >= 1) Hn(1) = 2.0_dp * x
            do n = 1, nmax - 1
                Hn(n + 1) = 2.0_dp * x * Hn(n) - 2.0_dp * n * Hn(n - 1)
            end do
        end subroutine hermite

        ! 【目的】Strutinsky の補正多項式の係数 c_m = (-1)^m / (4^m m!)
        pure function corr_coef(m) result(c)
            ! c_m = (-1)^m / (4^m m!)
            implicit none
            integer, intent(in) :: m
            real(dp) :: c
            integer :: k
            c = 1.0_dp
            do k = 1, m
                c = -c / (4.0_dp * k)
            end do
        end function corr_coef

        ! 【目的】平滑化準位密度 g(lambda) [1/MeV] (スピン両方, 1 準位 1 核子として数える).
        pure function smooth_density(e, lambda, gam, p) result(g)
            ! g(lambda) = (1/gam) sum_i f((lambda - e_i)/gam)
            implicit none
            real(dp), intent(in) :: e(:), lambda, gam
            integer, intent(in) :: p
            real(dp) :: g, x, s, Hn(0:p)
            integer :: i, m
            g = 0.0_dp
            do i = 1, size(e)
                x = (lambda - e(i)) / gam
                if (abs(x) > 12.0_dp) cycle
                call hermite(x, p, Hn)
                s = 0.0_dp
                do m = 0, p / 2
                    s = s + corr_coef(m) * Hn(2 * m)
                end do
                g = g + exp(-x * x) * s
            end do
            g = g / (sqrt(pi) * gam)
        end function smooth_density

        ! 【目的】平滑化した粒子数 N(lambda) = int_{-inf}^lambda g(e) de (解析的に積分した式).
        pure function smooth_number(e, lambda, gam, p) result(nn)
            ! N(lambda) = sum_i [ erfc(-x)/2 - exp(-x^2)/sqrt(pi) sum_{m>=1} c_m H_{2m-1}(x) ]
            ! (int_{-inf}^x exp(-t^2) H_n(t) dt = -exp(-x^2) H_{n-1}(x), n >= 1)
            implicit none
            real(dp), intent(in) :: e(:), lambda, gam
            integer, intent(in) :: p
            real(dp) :: nn, x, s, Hn(0:p)
            integer :: i, m
            nn = 0.0_dp
            do i = 1, size(e)
                x = (lambda - e(i)) / gam
                if (x > 12.0_dp) then
                    nn = nn + 1.0_dp
                    cycle
                else if (x < -12.0_dp) then
                    cycle
                end if
                call hermite(x, p, Hn)
                s = 0.0_dp
                do m = 1, p / 2
                    s = s + corr_coef(m) * Hn(2 * m - 1)
                end do
                nn = nn + 0.5_dp * erfc(-x) - exp(-x * x) * s / sqrt(pi)
            end do
        end function smooth_number

        ! 【目的】平滑化したエネルギー E(lambda) = int_{-inf}^lambda e g(e) de [MeV].
        pure function smooth_energy(e, lambda, gam, p) result(es)
            ! E(lambda) = sum_i [ e_i n_i(x) + gam int_{-inf}^x t f(t) dt ],
            ! int_{-inf}^x t f(t) dt = -exp(-x^2)/sqrt(pi) sum_m c_m [H_2m(x)/2 + 2m H_{2m-2}(x)]
            implicit none
            real(dp), intent(in) :: e(:), lambda, gam
            integer, intent(in) :: p
            real(dp) :: es, x, s1, s2, Hn(0:p)
            integer :: i, m
            es = 0.0_dp
            do i = 1, size(e)
                x = (lambda - e(i)) / gam
                if (x > 12.0_dp) then
                    es = es + e(i)
                    cycle
                else if (x < -12.0_dp) then
                    cycle
                end if
                call hermite(x, p, Hn)
                s1 = 0.0_dp
                s2 = 0.5_dp * Hn(0)
                do m = 1, p / 2
                    s1 = s1 + corr_coef(m) * Hn(2 * m - 1)
                    s2 = s2 + corr_coef(m) * (0.5_dp * Hn(2 * m) + 2.0_dp * m * Hn(2 * m - 2))
                end do
                es = es + e(i) * (0.5_dp * erfc(-x) - exp(-x * x) * s1 / sqrt(pi)) &
                        - gam * exp(-x * x) * s2 / sqrt(pi)
            end do
        end function smooth_energy

        !==========================================================================
        ! Pairing corrections: BCS and Lipkin-Nogami
        !==========================================================================

        subroutine pairing_window(e, n_part, ek, N1, N2, nk, n_int)
            !---------------------------
            ! doubly degenerate levels ek(k) = mean of the Kramers pair (2k-1, 2k)
            ! window (Moller-Nix Sec. 2.5): N2 = last bound level or last level within 5 MeV
            ! above the Fermi level, whichever is higher; the same number of levels below.
            ! nk: sharp occupation (2, 1, 0), n_int: nucleons in the window (Eq. 94)
            !---------------------------
            ! 【日本語】
            ! 【目的】二重縮退準位 ek と対相関の窓 N1..N2, 各準位の鋭い占有数 nk を決める.
            ! ek(k) は時間反転の対 (2k-1, 2k) の平均. n_int は窓内の核子数 (式 94).
            implicit none
            real(dp), intent(in) :: e(:)
            integer, intent(in) :: n_part
            real(dp), allocatable, intent(out) :: ek(:), nk(:)
            integer, intent(out) :: N1, N2
            real(dp), intent(out) :: n_int
            integer :: nlev, kF, k

            nlev = size(e) / 2
            allocate(ek(nlev))
            do k = 1, nlev
                ek(k) = 0.5_dp * (e(2 * k - 1) + e(2 * k))
            end do
            kF = (n_part + 1) / 2          ! last (partly) occupied level
            N2 = kF
            do k = kF, nlev
                if (ek(k) < 0.0_dp .or. ek(k) <= ek(kF) + 5.0_dp) N2 = k
            end do
            N2 = max(N2, kF + 1)
            N1 = max(1, kF - (N2 - kF) + 1)
            if (N2 > nlev) stop "pairing_window: pairing window exceeds the basis"
            allocate(nk(N1:N2))
            do k = N1, N2
                if (2 * k <= n_part) then
                    nk(k) = 2.0_dp
                else if (2 * k - 1 == n_part) then
                    nk(k) = 1.0_dp
                else
                    nk(k) = 0.0_dp
                end if
            end do
            n_int = real(n_part - 2 * (N1 - 1), dp)
        end subroutine pairing_window

        ! 【目的】平均的な原子核で対ギャップが Delta_G になるように G を決める (式 106, 107).
        ! y1, y2 は窓の下端・上端の, 平均準位密度 rho で測ったフェルミ準位からのエネルギー.
        subroutine pairing_strength(n_part, N1, N2, rho, delta_G, G, y1, y2)
            ! G from the average gap equation (Eq. 107):
            !   1/G = rho/2 [asinh(y2/Delta_G) - asinh(y1/Delta_G)],  y1, y2 from Eq. (106)
            implicit none
            integer, intent(in) :: n_part, N1, N2
            real(dp), intent(in) :: rho, delta_G
            real(dp), intent(out) :: G, y1, y2
            y1 = (-0.5_dp * n_part + N1 - 1) / rho
            y2 = (-0.5_dp * n_part + N2) / rho
            G = 2.0_dp / (rho * (asinh(y2 / delta_G) - asinh(y1 / delta_G)))
        end subroutine pairing_strength

        subroutine solve_gap(eps, n_int, G, delta, lambda)
            !---------------------------
            ! BCS gap and number equations (Eqs. 94-96) for the levels eps(:):
            !   2/G = sum 1/sqrt((eps_k - lambda)^2 + delta^2),  2 sum v_k^2 = n_int
            ! gap equation by bisection in delta, lambda(delta) from the number equation.
            ! delta = 0 if there is no non-trivial solution.
            !---------------------------
            ! 【日本語】
            ! 【目的】準位 eps について BCS のギャップ方程式と粒子数方程式を解き,
            ! Delta と lambda を求める.
            ! Delta について二分法, 各 Delta で lambda を粒子数方程式から二分法で決める.
            implicit none
            real(dp), intent(in) :: eps(:), n_int, G
            real(dp), intent(out) :: delta, lambda
            real(dp) :: lo, hi
            integer :: it

            lo = 1.0e-6_dp
            hi = 30.0_dp
            if (gap_function(lo) >= 0.0_dp) then
                delta = 0.0_dp
                lambda = fermi_level(1.0e-6_dp)
                return
            end if
            if (gap_function(hi) < 0.0_dp) stop "solve_gap: gap not bracketed"
            do it = 1, 200
                delta = 0.5_dp * (lo + hi)
                if (gap_function(delta) < 0.0_dp) then
                    lo = delta
                else
                    hi = delta
                end if
                if (hi - lo < 1.0e-12_dp) exit
            end do
            delta = 0.5_dp * (lo + hi)
            lambda = fermi_level(delta)

        contains

            function fermi_level(d) result(lam)
                ! 2 sum v_k^2 = n_int for fixed gap d (monotone in lambda)
                real(dp), intent(in) :: d
                real(dp) :: lam, a, b
                integer :: j
                a = minval(eps) - 50.0_dp
                b = maxval(eps) + 50.0_dp
                do j = 1, 200
                    lam = 0.5_dp * (a + b)
                    if (sum(1.0_dp - (eps - lam) / sqrt((eps - lam)**2 + d * d)) < n_int) then
                        a = lam
                    else
                        b = lam
                    end if
                    if (b - a < 1.0e-13_dp) exit
                end do
                lam = 0.5_dp * (a + b)
            end function fermi_level

            function gap_function(d) result(f)
                ! 2/G - sum 1/sqrt((eps_k - lambda)^2 + d^2); increases with d, zero at the solution
                real(dp), intent(in) :: d
                real(dp) :: f, lam
                lam = fermi_level(d)
                f = 2.0_dp / G - sum(1.0_dp / sqrt((eps - lam)**2 + d * d))
            end function gap_function

        end subroutine solve_gap

        subroutine bcs_correction(e, n_part, rho, delta_G, G, delta, lambda, N1, N2, e_pc, e_pc_ave, delta_odd)
            !---------------------------
            ! e(:)    : single-particle energies, ascending, one nucleon each (Kramers pairs)
            ! rho     : average density of doubly degenerate levels = g_smooth(lambda_smooth)/2
            ! delta_G : effective-interaction pairing gap
            ! delta_odd (optional, default delta_G): average gap added to e_pc_ave for odd
            !           particle numbers; must be the same term as in the macroscopic energy
            ! out: G, BCS gap delta and Fermi level lambda, window N1..N2 (doubly degenerate
            !      levels numbered from the bottom), e_pc (Moller-Nix Eq. 5), e_pc_ave (Eq. 27)
            ! Odd particle numbers: no blocking, quasi-particle energy of the odd level added.
            !---------------------------
            ! 【日本語】
            ! 【目的】BCS による対補正 E_pair = e_pc - e_pc_ave を求める (Moller-Nix 1992).
            !   v_k^2 = [1 - (e_k - lambda) / sqrt((e_k - lambda)^2 + Delta^2)] / 2  (占有確率)
            !   e_pc  = sum (2 v_k^2 - n_k) e_k - Delta^2/G - G/2 sum (2 v_k^4 - n_k)  (Moller-Nix 式 5)
            !   閉殻付近では Delta = 0 (解が潰れる) になり, e_pc = 0 となる.
            implicit none
            real(dp), intent(in) :: e(:), rho, delta_G
            integer, intent(in) :: n_part
            real(dp), intent(out) :: G, delta, lambda, e_pc, e_pc_ave
            integer, intent(out) :: N1, N2
            real(dp), intent(in), optional :: delta_odd
            real(dp), allocatable :: ek(:), v2(:), nk(:)
            real(dp) :: y1, y2, n_int
            integer :: kF
            logical :: odd

            odd = mod(n_part, 2) == 1
            kF = (n_part + 1) / 2
            call pairing_window(e, n_part, ek, N1, N2, nk, n_int)
            call pairing_strength(n_part, N1, N2, rho, delta_G, G, y1, y2)

            ! average pairing correlation + quasi-particle energy (Moller-Nix Eq. 27)
            e_pc_ave = 0.5_dp * rho * ((y2 - G) * (y2 - sqrt(y2**2 + delta_G**2)) &
                                     + (y1 - G) * (y1 + sqrt(y1**2 + delta_G**2))) &
                     + 0.25_dp * G * rho * delta_G * (atan(y2 / delta_G) - atan(y1 / delta_G))
            if (odd) then
                if (present(delta_odd)) then
                    e_pc_ave = e_pc_ave + delta_odd
                else
                    e_pc_ave = e_pc_ave + delta_G
                end if
            end if

            call solve_gap(ek(N1:N2), n_int, G, delta, lambda)
            if (delta == 0.0_dp) then
                ! no non-trivial solution: sharp distribution, no correlation energy
                lambda = 0.5_dp * (ek(kF) + ek(min(kF + 1, size(ek))))
                e_pc = 0.0_dp
                return
            end if
            allocate(v2(N1:N2))
            v2 = 0.5_dp * (1.0_dp - (ek(N1:N2) - lambda) / sqrt((ek(N1:N2) - lambda)**2 + delta**2))

            ! pairing correlation + quasi-particle energy (Moller-Nix Eq. 5)
            e_pc = sum((2.0_dp * v2 - nk) * ek(N1:N2)) - delta**2 / G &
                 - 0.5_dp * G * sum(2.0_dp * v2**2 - nk)
            if (odd) e_pc = e_pc + sqrt((ek(kF) - lambda)**2 + delta**2)
        end subroutine bcs_correction

        subroutine ln_correction(e, n_part, rho, delta_G, G, delta, lambda, lambda2, N1, N2, &
                                 e_pc, e_pc_ave, delta_odd)
            !---------------------------
            ! Lipkin-Nogami pairing correction (FRDM2012 Eqs. 94-101, 107-110), the pairing
            ! model of FRDM2012. Arguments as in bcs_correction, plus
            !   lambda2 : number-fluctuation constant (Eq. 98)
            ! The coupled equations are solved as one equation for lambda2 (see below); the
            ! BCS equations for the shifted levels eps_k = e_k + (4 lambda2 - G) v_k^2 (Eq. 97)
            ! are solved inside. Unlike BCS, the solution never collapses (magic numbers).
            ! E_pc  : Eq. (101); odd particle numbers: no blocking, E_i of Eq. (100) added.
            ! <E_pc>: Eq. (108) with lambda2 from Eqs. (109), (110):
            !   lambda2_ave = G/4 (A - C)/(B - C)
            !   A = (rho Delta_G/4)^2 [(2/(G rho))^2 - (ln((y2^2+Delta_G^2)/(y1^2+Delta_G^2))/2)^2]
            !     = (sum u^3 v)(sum u v^3)
            !   B = (rho Delta_G/4)^2 [atan(y2/Delta_G) - atan(y1/Delta_G)]^2 = (sum u^2 v^2)^2
            !   C = rho Delta_G/32 [Delta_G (y2/(y2^2+Delta_G^2) - y1/(y1^2+Delta_G^2))
            !       + atan(y2/Delta_G) - atan(y1/Delta_G)] = sum u^4 v^4
            !---------------------------
            ! 【日本語】
            ! 【目的】Lipkin-Nogami (LN) による対補正 E_pair = e_pc - e_pc_ave (FRDM2012 の方法).
            !   LN は粒子数の射影を近似的に取り入れた BCS で, 閉殻でも Delta > 0 の解をもつ.
            !   準位をずらした eps_k = e_k + (4 lambda2 - G) v_k^2 (式 97) について BCS を解き,
            !   lambda2 (式 98) が自己無撞着になるまで二分法で探す.
            !   e_pc = BCS と同じ項 - 4 lambda2 sum u_k^2 v_k^2  (式 101)
            !   e_pc_ave は平均的な原子核での値 (式 108-110).
            implicit none
            real(dp), intent(in) :: e(:), rho, delta_G
            integer, intent(in) :: n_part
            real(dp), intent(out) :: G, delta, lambda, lambda2, e_pc, e_pc_ave
            integer, intent(out) :: N1, N2
            real(dp), intent(in), optional :: delta_odd
            real(dp), allocatable :: ek(:), nk(:), eps(:), v2(:), u2(:)
            real(dp) :: y1, y2, n_int, d, at, aa, bb, cc, l2_ave, lo, hi, f_ln, f_prev, l2_max
            integer :: kF, it
            integer, parameter :: n_scan = 200
            logical :: odd, bracketed, restart

            odd = mod(n_part, 2) == 1
            kF = (n_part + 1) / 2
            call pairing_window(e, n_part, ek, N1, N2, nk, n_int)
            call pairing_strength(n_part, N1, N2, rho, delta_G, G, y1, y2)

            ! average quantities, Eqs. (108)-(110)
            d = delta_G
            at = atan(y2 / d) - atan(y1 / d)
            aa = (rho * d / 4.0_dp)**2 * ((2.0_dp / (G * rho))**2 &
                 - (0.5_dp * log((y2**2 + d**2) / (y1**2 + d**2)))**2)
            bb = (rho * d / 4.0_dp)**2 * at**2
            cc = rho * d / 32.0_dp * (d * (y2 / (y2**2 + d**2) - y1 / (y1**2 + d**2)) + at)
            l2_ave = 0.25_dp * G * (aa - cc) / (bb - cc)
            e_pc_ave = 0.5_dp * rho * ((y2 - G) * (y2 - sqrt(y2**2 + d**2)) &
                                     + (y1 - G) * (y1 + sqrt(y1**2 + d**2))) &
                     + 0.25_dp * (G - 4.0_dp * l2_ave) * rho * d * at
            if (odd) then
                if (present(delta_odd)) then
                    e_pc_ave = e_pc_ave + delta_odd
                else
                    e_pc_ave = e_pc_ave + delta_G
                end if
            end if

            ! microscopic Lipkin-Nogami equations: F(lambda2) = lambda2(Eq. 98) - lambda2 = 0
            ! by bisection. For given lambda2 the BCS equations with the shifted levels are
            ! solved self-consistently in v_k^2 (inner_solution). As delta -> 0 the lambda2 of
            ! Eq. (98) diverges, so F > 0 for small lambda2 (collapsed: F = +inf) and F < 0 for
            ! large lambda2; the solution always has delta > 0.
            ! Very light nuclei: with a window of a few (near-)degenerate levels the iteration
            ! started from the previous v_k^2 can fall into the collapsed solution for every
            ! lambda2 tried. Then lambda2 is scanned on [0, 2 max(<lambda2>, G)], each point
            ! iterated from the sharp occupations, and the bisection does the same (restart).
            ! If there is still no sign change, BCS (lambda2 = 0) is used with a warning.
            allocate(eps(N1:N2), v2(N1:N2), u2(N1:N2))
            v2 = 0.5_dp * nk
            restart = .false.
            lo = 0.0_dp
            hi = max(l2_ave, G)
            bracketed = .false.
            do it = 1, 60
                if (ln_function(hi) < 0.0_dp) then
                    bracketed = .true.
                    exit
                end if
                lo = hi
                hi = 2.0_dp * hi
            end do
            if (.not. bracketed) then
                restart = .true.
                l2_max = 2.0_dp * max(l2_ave, G)
                f_prev = ln_function(0.0_dp)
                do it = 1, n_scan
                    hi = l2_max * it / n_scan
                    f_ln = ln_function(hi)
                    if (f_prev > 0.0_dp .and. f_ln <= 0.0_dp) then
                        lo = l2_max * (it - 1) / n_scan
                        bracketed = .true.
                        exit
                    end if
                    f_prev = f_ln
                end do
            end if
            if (.not. bracketed) then
                print '(a, i4, a)', " ln_correction: no Lipkin-Nogami solution for ", n_part, &
                                    " nucleons, BCS used (lambda2 = 0)"
                call bcs_correction(e, n_part, rho, delta_G, G, delta, lambda, N1, N2, e_pc, e_pc_ave, delta_odd)
                lambda2 = 0.0_dp
                return
            end if
            do it = 1, 100
                lambda2 = 0.5_dp * (lo + hi)
                if (ln_function(lambda2) > 0.0_dp) then
                    lo = lambda2
                else
                    hi = lambda2
                end if
                if (hi - lo < 1.0e-10_dp) exit
            end do
            lambda2 = 0.5_dp * (lo + hi)
            f_ln = ln_function(lambda2)
            eps = ek(N1:N2) + (4.0_dp * lambda2 - G) * v2

            if (delta == 0.0_dp) then
                lambda2 = 0.0_dp
                e_pc = 0.0_dp
                return
            end if
            u2 = 1.0_dp - v2
            ! Eq. (101)
            e_pc = sum((2.0_dp * v2 - nk) * ek(N1:N2)) - delta**2 / G &
                 - 0.5_dp * G * sum(2.0_dp * v2**2 - nk) - 4.0_dp * lambda2 * sum(u2 * v2)
            if (odd) e_pc = e_pc + sqrt((eps(kF) - lambda)**2 + delta**2) + lambda2

        contains

            ! 【目的】l2 で解いた BCS 解から式 (98) の lambda2 を計算し, l2 との差を返す (0 が解).
            function ln_function(l2) result(f)
                ! lambda2 of Eq. (98) minus l2, for the self-consistent solution at l2
                real(dp), intent(in) :: l2
                real(dp) :: f, s22, s44
                call inner_solution(l2)
                if (delta == 0.0_dp) then
                    f = huge(1.0_dp)
                    return
                end if
                u2 = 1.0_dp - v2
                s22 = sum(u2 * v2)
                s44 = sum(u2**2 * v2**2)
                ! (sum u^3 v)(sum u v^3) = (sum u v u^2)(sum u v v^2)
                f = 0.25_dp * G * (sum(sqrt(u2 * v2) * u2) * sum(sqrt(u2 * v2) * v2) - s44) &
                    / (s22**2 - s44) - l2
            end function ln_function

            ! 【目的】lambda2 = l2 を固定して, ずらした準位の BCS 方程式を v_k^2 について反復で解く.
            subroutine inner_solution(l2)
                ! BCS equations for eps_k = e_k + (4 l2 - G) v_k^2, iterated in v_k^2
                ! (starts from the v_k^2 of the previous call, or from the sharp occupations
                ! when restart is set)
                real(dp), intent(in) :: l2
                real(dp) :: v2_new(size(v2))
                integer :: j
                if (restart) v2 = 0.5_dp * nk
                do j = 1, 500
                    eps = ek(N1:N2) + (4.0_dp * l2 - G) * v2
                    call solve_gap(eps, n_int, G, delta, lambda)
                    if (delta == 0.0_dp) then
                        v2_new = 0.5_dp * nk
                    else
                        v2_new = 0.5_dp * (1.0_dp - (eps - lambda) / sqrt((eps - lambda)**2 + delta**2))
                    end if
                    if (maxval(abs(v2_new - v2)) < 1.0e-11_dp) exit
                    v2 = 0.5_dp * (v2 + v2_new)
                end do
                v2 = v2_new
            end subroutine inner_solution

        end subroutine ln_correction

        !==========================================================================
        ! relative surface area of an ellipsoid
        !==========================================================================

        ! 【目的】楕円体の相対表面積 B_s = S / (4 pi R^2), R^3 = R_x R_y R_z (検算・テスト用).
        function ellipsoid_relative_area(shape) result(B_s)
            ! surface area / (4 pi R^2), R^3 = R_x R_y R_z
            implicit none
            type(ellipsoid_shape), intent(in) :: shape
            real(dp) :: B_s
            integer, parameter :: nq = 64
            real(dp) :: tq(nq), wq(nq), th, ph, a, b, c, st, s
            integer :: i, j

            a = shape%semi(1)
            b = shape%semi(2)
            c = shape%semi(3)
            call gauss_legendre(nq, tq, wq)
            s = 0.0_dp
            ! theta in [0, pi/2], phi in [0, pi/2]: one octant, times 8
            do j = 1, nq
                ph = 0.25_dp * pi * (tq(j) + 1.0_dp)
                do i = 1, nq
                    th = 0.25_dp * pi * (tq(i) + 1.0_dp)
                    st = sin(th)
                    s = s + wq(i) * wq(j) * st * sqrt((b * c * st * cos(ph))**2 + (a * c * st * sin(ph))**2 &
                                                      + (a * b * cos(th))**2)
                end do
            end do
            s = 8.0_dp * s * (0.25_dp * pi)**2
            B_s = s / (4.0_dp * pi * (a * b * c)**(2.0_dp / 3.0_dp))
        end function ellipsoid_relative_area

end module shell_bcs_mod
