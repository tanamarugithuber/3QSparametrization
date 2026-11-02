program pes3qs
    !---------------------------
    ! Potential-energy surface of one nucleus on Moller's five-dimensional 3QS grid
    !   P. Moller et al., Phys. Rev. C 79 (2009) 064304, Sec. II.A.3 and Fig. 2
    !   coordinates: Q2 (elongation), d (neck), eps_f1, eps_f2 (end-body deformations),
    !                alpha_g (mass asymmetry)
    !   E_pot = E_mac(shape) - E_mac(sphere) [+ E_shell + E_pair if micro = .true.]
    !   (all energies relative to the macroscopic energy of the spherical shape)
    ! E_mac: FRLDM of frldm_mod with the 3QS shape (B_1, B_3 by ring sums, B_W of Eq. (75)).
    ! E_s+p: shell_bcs_mod with the 3QS shape (optional; about 10 - 30 s per shape).
    !
    ! Usage:  ./pes3qs.exe [input.nml]      (namelist &pes3qs_input, see README_3QS.md)
    ! Output: one line per grid point (columns in the header of the output file), read by
    !         saddle3qs.exe. With restart = .true. the (Q2, eps_f1, eps_f2, alpha_g) columns
    !         already in the file are skipped.
    !---------------------------
    ! 【日本語での説明】
    ! 1 つの原子核について, PRC 79 の 5 次元 3QS 格子 (Q2, d, eps_f1, eps_f2, alpha_g) の
    ! 各点でポテンシャルエネルギーを計算し, ファイルに書く. saddle3qs.exe がこれを読んで鞍点を探す.
    !   E_pot = E_mac(形) - E_mac(球) [+ E_s+p]   (球の巨視的エネルギーを基準にする, 論文と同じ)
    ! Q2 の値は (Q2/b)^(1/2) = sqrtq2_min + k sqrtq2_step (k = 0, ..., n_q2 - 1).
    ! d は (Q2, eps_f1, eps_f2, alpha_g) ごとに [d_min, d_max] を n_d 等分した値 (qs3_coord_mod).
    ! (Q2, eps_f1, eps_f2, alpha_g) の組ごとに OpenMP で並列に計算する.
    !---------------------------
    !$ use omp_lib
    use iso_fortran_env, only: real64, output_unit
    use constant_mod, only: pi
    use nucleus_mod, only: nucleus_property
    use shape_mod, only: epsilon_shape, epsilon_shape_from
    use qs3_shape_mod, only: qs3_shape
    use qs3_coord_mod, only: neck_branch, find_neck_branch, shape_at_neck
    use frldm_mod, only: frldm_variables
    use micro_constant_mod, only: microscopic_variables
    use shell_bcs_mod, only: shell_bcs_calculation, shell_bcs_result
    implicit none
    integer, parameter :: dp = real64
    integer, parameter :: max_list = 64

    ! namelist
    integer :: Z = 94, N = 146
    real(dp) :: sqrtq2_min = 2.0_dp, sqrtq2_step = 0.25_dp
    integer :: n_q2 = 45
    integer :: n_d = 15
    integer :: n_eps_f = 15
    real(dp) :: eps_f(max_list) = 0.0_dp
    real(dp) :: alpha_g_min = -0.02_dp, alpha_g_step = 0.02_dp
    integer :: n_alpha_g = 35
    logical :: same_eps_f = .false.     ! eps_f2 = eps_f1 only (4D)
    logical :: micro = .false.
    integer :: n_gh = 44, n_basis = 12
    character(len=8) :: pairing = 'LN'
    integer :: n_quad = 64              ! FRLDM B_1, B_3 resolution (frldm_mod)
    character(len=256) :: out = 'pes3qs.dat'
    logical :: restart = .false.
    namelist /pes3qs_input/ Z, N, sqrtq2_min, sqrtq2_step, n_q2, n_d, n_eps_f, eps_f, alpha_g_min, &
        alpha_g_step, n_alpha_g, same_eps_f, micro, n_gh, n_basis, pairing, n_quad, out, restart

    type(nucleus_property) :: nuc
    type(frldm_variables) :: fr
    type(epsilon_shape) :: sph
    real(dp) :: E_sph, t_start, t_end
    character(len=256) :: fname
    integer :: u_in, ios, u_out, n_combo, ic, iq, ie1, ie2, ia, n_e2, n_done_combo, n_valid
    integer, allocatable :: combo(:, :)
    logical, allocatable :: done(:, :, :, :)
    logical :: exists

    eps_f(1:15) = [-0.2_dp, -0.15_dp, -0.1_dp, 0.0_dp, 0.1_dp, 0.15_dp, 0.175_dp, 0.2_dp, &
                   0.225_dp, 0.25_dp, 0.275_dp, 0.3_dp, 0.35_dp, 0.4_dp, 0.5_dp]       ! Eq. (8)
    fname = 'pes3qs.nml'
    if (command_argument_count() >= 1) call get_command_argument(1, fname)
    inquire(file=trim(fname), exist=exists)
    if (exists) then
        open(newunit=u_in, file=trim(fname), status='old', action='read')
        read(u_in, nml=pes3qs_input, iostat=ios)
        if (ios /= 0) error stop "pes3qs: error in the namelist &pes3qs_input"
        close(u_in)
    else
        print '(3a)', " ", trim(fname), " not found: default grid for 240Pu"
    end if
    if (n_eps_f > max_list) error stop "pes3qs: n_eps_f too large"

    nuc%Z = Z
    nuc%N = N
    nuc%A = Z + N
    nuc%R_a = nuc%R0 * real(nuc%A, dp)**(1.0_dp / 3.0_dp)
    nuc%mass = 0.0_dp
    fr%verbose = .false.
    sph = epsilon_shape_from(nuc%R_a, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp)
    call fr%calculate_frldm_energy(nuc, sph, E_sph)

    n_e2 = n_eps_f
    if (same_eps_f) n_e2 = 1
    allocate(done(n_q2, n_eps_f, n_e2, n_alpha_g))
    done = .false.
    if (restart) call read_done()

    n_combo = 0
    allocate(combo(4, n_q2 * n_eps_f * n_e2 * n_alpha_g))
    do ia = 1, n_alpha_g
        do ie2 = 1, n_e2
            do ie1 = 1, n_eps_f
                do iq = 1, n_q2
                    if (done(iq, ie1, ie2, ia)) cycle
                    n_combo = n_combo + 1
                    combo(:, n_combo) = [iq, ie1, ie2, ia]
                end do
            end do
        end do
    end do

    print '(a, i0, a, i0, a, i0)', " Z = ", Z, ", N = ", N, ", A = ", nuc%A
    print '(a, f12.4, a)', " E_mac(sphere) = ", E_sph, " MeV (reference)"
    print '(a, 5i5, a, i0, a)', " grid (Q2, d, eps_f1, eps_f2, alpha_g):", n_q2, n_d, n_eps_f, n_e2, n_alpha_g, &
          "   (", n_combo * n_d, " points to compute)"
    if (micro) print '(a, i0, a, i0, 2a)', " with shell + pairing corrections: n_gh = ", n_gh, ", n_basis = ", &
          n_basis, ", pairing = ", trim(pairing)

    if (restart .and. any(done)) then
        open(newunit=u_out, file=trim(out), status='old', position='append', action='write')
    else
        open(newunit=u_out, file=trim(out), status='replace', action='write')
        call write_header()
    end if

    call cpu_time(t_start)
    t_start = wall()
    n_done_combo = 0
    n_valid = 0
    !$omp parallel do default(shared) private(ic) schedule(dynamic)
    do ic = 1, n_combo
        call compute_combo(combo(:, ic))
    end do
    !$omp end parallel do
    t_end = wall()
    close(u_out)
    print '(a, i0, a, i0, a, f10.1, a)', " done: ", n_valid, " shapes of ", n_combo * n_d, " exist; ", &
          t_end - t_start, " s"
    print '(2a)', " output: ", trim(out)

contains

    function wall() result(t)
        real(dp) :: t
        integer(8) :: c, r
        call system_clock(c, r)
        t = real(c, dp) / real(r, dp)
    end function wall

    pure real(dp) function q2_of(iq)
        integer, intent(in) :: iq
        q2_of = (sqrtq2_min + sqrtq2_step * (iq - 1))**2
    end function q2_of

    pure real(dp) function alpha_g_of(ia)
        integer, intent(in) :: ia
        alpha_g_of = alpha_g_min + alpha_g_step * (ia - 1)
    end function alpha_g_of

    subroutine write_header()
        write(u_out, '(a, i0, a, i0)') "# 3QS potential-energy surface (PRC 79, 064304), Z = ", Z, ", N = ", N
        write(u_out, '(a, 5i6)') "# dims (Q2, d, eps_f1, eps_f2, alpha_g):", n_q2, n_d, n_eps_f, n_e2, n_alpha_g
        write(u_out, '(a, f14.6, a)') "# E_mac(sphere) = ", E_sph, " MeV; energies below are relative to it"
        write(u_out, '(a, l2)') "# micro (shell + pairing):", micro
        write(u_out, '(a)') "# iq id ie1 ie2 ia ok  Q2[b] d[fm] d/R eps_f1 eps_f2 alpha_g" // &
            "  sigma1 sigma2 sigma3 alpha1 alpha2 alpha3  B_s B_1 B_3 B_W  E_mac-E_sph E_s+p E_pot [MeV]"
    end subroutine write_header

    subroutine read_done()
        ! (iq, ie1, ie2, ia) with all n_d lines in the existing output
        integer :: u, k(6), cnt(n_q2, n_eps_f, n_e2, n_alpha_g)
        character(len=1024) :: line
        inquire(file=trim(out), exist=exists)
        if (.not. exists) return
        cnt = 0
        open(newunit=u, file=trim(out), status='old', action='read')
        do
            read(u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
            read(line, *, iostat=ios) k
            if (ios /= 0) cycle
            if (k(1) < 1 .or. k(1) > n_q2 .or. k(3) < 1 .or. k(3) > n_eps_f .or. k(4) < 1 .or. k(4) > n_e2 &
                .or. k(5) < 1 .or. k(5) > n_alpha_g) cycle
            cnt(k(1), k(3), k(4), k(5)) = cnt(k(1), k(3), k(4), k(5)) + 1
        end do
        close(u)
        done = cnt >= n_d
        print '(a, i0, a)', " restart: ", count(done), " (Q2, eps_f1, eps_f2, alpha_g) combinations already done"
    end subroutine read_done

    subroutine compute_combo(c)
        ! all n_d necks for one (Q2, eps_f1, eps_f2, alpha_g)
        integer, intent(in) :: c(4)
        type(neck_branch) :: br
        type(qs3_shape) :: s
        type(frldm_variables) :: frl
        type(microscopic_variables) :: mv
        type(shell_bcs_result) :: res
        real(dp) :: e1, e2, ag, Q2, d, E_mac, E_sp
        character(len=512) :: lines(n_d)
        integer :: id, nv

        Q2 = q2_of(c(1))
        e1 = eps_f(c(2))
        e2 = eps_f(c(3))
        if (same_eps_f) e2 = e1
        ag = alpha_g_of(c(4))
        frl%verbose = .false.
        br = find_neck_branch(nuc%R_a, Q2, e1, e2, ag)
        nv = 0
        do id = 1, n_d
            s%valid = .false.
            if (br%ok) then
                if (n_d > 1) then
                    d = br%d_min + (br%d_max - br%d_min) * (id - 1) / real(n_d - 1, dp)
                else
                    d = br%d_max
                end if
                s = shape_at_neck(br, d)
            end if
            if (.not. s%valid) then
                write(lines(id), '(5i4, i3, f9.3, 2f9.4, 3f8.4, 22x)') c(1), id, c(2), c(3), c(4), 0, Q2, &
                    0.0_dp, 0.0_dp, e1, e2, ag
                cycle
            end if
            nv = nv + 1
            call frl%calculate_frldm_energy(nuc, s, E_mac)
            E_sp = 0.0_dp
            if (micro) then
                call shell_bcs_calculation(nuc, mv, res, n_gh, trim(pairing), qs3=s, n_basis=n_basis)
                E_sp = res%e_micro
            end if
            write(lines(id), '(5i4, i3, f9.3, 2f9.4, 3f8.4, 6es15.7, 4f10.6, 3f12.4)') c(1), id, c(2), c(3), c(4), 1, &
                s%q2_barn(), s%neck_diameter(), s%neck_diameter() / s%R, e1, e2, ag, &
                s%sigma1, s%sigma2, s%sigma3, s%alpha1, s%alpha2, s%alpha3, &
                frl%B_s, frl%B_1, frl%B_3, frl%B_w, E_mac - E_sph, E_sp, E_mac - E_sph + E_sp
        end do
        !$omp critical (pes_output)
        do id = 1, n_d
            write(u_out, '(a)') trim(lines(id))
        end do
        flush(u_out)
        n_done_combo = n_done_combo + 1
        n_valid = n_valid + nv
        if (mod(n_done_combo, max(1, n_combo / 20)) == 0 .or. n_done_combo == n_combo) then
            print '(a, i0, a, i0, a, f9.1, a)', "   ", n_done_combo, " / ", n_combo, " combinations, ", &
                wall() - t_start, " s"
            flush(output_unit)
        end if
        !$omp end critical (pes_output)
    end subroutine compute_combo

end program pes3qs
