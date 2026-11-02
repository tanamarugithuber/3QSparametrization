program saddle3qs
    !---------------------------
    ! Fission saddle point on a 3QS potential-energy surface written by pes3qs.exe,
    ! by the immersion method of P. Moller et al., Phys. Rev. C 79 (2009) 064304, Sec. III
    ! (immersion_mod).
    !
    ! Usage:  ./saddle3qs.exe [saddle.nml]     (namelist &saddle3qs_input, see README_3QS.md)
    !   file    : output of pes3qs.exe
    !   energy  : 'pot' (E_pot, default) or 'mac' (E_mac - E_mac(sphere) only)
    !   entry_point: grid indices (iq, id, ie1, ie2, ia) of the entry point; 0 = automatic:
    !             the lowest point with iq <= entry_iq_max
    !   exit_point : grid indices of the exit point; 0 = automatic: the lowest point with
    !             iq >= exit_iq_min (default: the largest Q2 of the grid)
    !   noise   : amplitude of the random numbers added to E (MeV, default 0.001), seed
    !   n_saddles: number of saddles (the lowest one and n_saddles - 1 more by the dam method)
    !   path_file: the path of the water from the entry to the exit (the fission path through
    !             the lowest saddle) is written there
    !---------------------------
    ! 【日本語での説明】
    ! pes3qs.exe の出力 (5 次元格子上のポテンシャル) を読み, 入口 (基底状態側) から出口
    ! (大きい Q2 の谷) への一番低い鞍点を immersion 法で求める (immersion_mod).
    ! 入口・出口を指定しなければ, 入口は iq <= entry_iq_max の最低点, 出口は iq >= exit_iq_min の最低点.
    ! n_saddles > 1 ならダム法で次の鞍点も探す. 鞍点を通る水の道を path_file に書く.
    !---------------------------
    use iso_fortran_env, only: real64, real32
    use immersion_mod
    implicit none
    integer, parameter :: dp = real64, sp = real32

    character(len=256) :: file = 'pes3qs.dat', path_file = 'saddle3qs_path.dat'
    character(len=8) :: energy = 'pot'
    integer :: entry_point(5) = 0, exit_point(5) = 0
    integer :: entry_iq_max = 1, exit_iq_min = 0
    real(dp) :: noise = 0.001_dp
    integer :: seed = 20090601
    integer :: n_saddles = 1
    namelist /saddle3qs_input/ file, energy, entry_point, exit_point, entry_iq_max, exit_iq_min, noise, seed, &
        n_saddles, path_file

    character(len=256) :: fname
    character(len=2048) :: line
    integer :: dims(5), u, ios, k(6), np, i, i_entry, i_exit, i_sad, n_path, nf, j, ZZ, NN, A
    real(dp) :: v(19), E_sph, Es
    real(dp), allocatable :: E(:), E_raw(:), Es_list(:)
    real(sp), allocatable :: Q2(:), dR(:), d(:), e1(:), e2(:), ag(:)
    logical, allocatable :: valid(:), mask(:)
    integer, allocatable :: path(:), is_list(:)
    logical :: exists, have_dims
    integer :: sub(5)

    fname = 'saddle3qs.nml'
    if (command_argument_count() >= 1) call get_command_argument(1, fname)
    inquire(file=trim(fname), exist=exists)
    if (exists) then
        open(newunit=u, file=trim(fname), status='old', action='read')
        read(u, nml=saddle3qs_input, iostat=ios)
        if (ios /= 0) error stop "saddle3qs: error in the namelist &saddle3qs_input"
        close(u)
    end if

    !--------------------------- read the surface
    open(newunit=u, file=trim(file), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop "saddle3qs: cannot open the surface file"
    have_dims = .false.
    ZZ = 0
    NN = 0
    E_sph = 0.0_dp
    do
        read(u, '(a)', iostat=ios) line
        if (ios /= 0) exit
        if (line(1:1) /= '#') cycle
        if (index(line, '# dims') == 1) then
            read(line(index(line, ':') + 1:), *) dims
            have_dims = .true.
        else if (index(line, 'Z =') > 0 .and. index(line, 'N =') > 0) then
            read(line(index(line, 'Z =') + 3:), *) ZZ
            read(line(index(line, 'N =') + 3:), *) NN
        else if (index(line, 'E_mac(sphere)') > 0) then
            read(line(index(line, '=') + 1:), *) E_sph
        end if
    end do
    if (.not. have_dims) error stop "saddle3qs: no '# dims' line in the surface file"
    A = ZZ + NN
    np = product(dims)
    allocate(E(np), valid(np), Q2(np), dR(np), d(np), e1(np), e2(np), ag(np))
    valid = .false.
    E = huge(1.0_dp)
    Q2 = 0.0
    dR = 0.0
    d = 0.0
    e1 = 0.0
    e2 = 0.0
    ag = 0.0
    rewind(u)
    do
        read(u, '(a)', iostat=ios) line
        if (ios /= 0) exit
        if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
        read(line, *, iostat=ios) k
        if (ios /= 0) cycle
        i = sub_to_index(k(1:5), dims)
        if (k(6) == 1) then
            read(line, *, iostat=ios) k, v
            if (ios /= 0) cycle
            valid(i) = .true.
            if (trim(energy) == 'mac') then
                E(i) = v(17)
            else
                E(i) = v(19)
            end if
            Q2(i) = real(v(1), sp)
            d(i) = real(v(2), sp)
            dR(i) = real(v(3), sp)
            e1(i) = real(v(4), sp)
            e2(i) = real(v(5), sp)
            ag(i) = real(v(6), sp)
        end if
    end do
    close(u)
    print '(a, i0, a, i0, a, 5i4)', " surface: Z = ", ZZ, ", N = ", NN, ", dims", dims
    print '(a, i0, a, i0)', " shapes that exist: ", count(valid), " of ", np
    if (count(valid) == 0) error stop "saddle3qs: empty surface"

    allocate(E_raw(np))
    E_raw = E
    if (noise > 0.0_dp) call add_noise(E, valid, noise, seed)

    !--------------------------- entry and exit
    allocate(mask(np))
    if (all(entry_point > 0)) then
        i_entry = sub_to_index(entry_point, dims)
        if (.not. valid(i_entry)) error stop "saddle3qs: the entry point does not exist"
    else
        do i = 1, np
            call index_to_sub(i, dims, sub)
            mask(i) = valid(i) .and. sub(1) <= max(1, entry_iq_max)
        end do
        if (.not. any(mask)) error stop "saddle3qs: no entry point"
        i_entry = minloc(E, dim=1, mask=mask)
    end if
    if (all(exit_point > 0)) then
        i_exit = sub_to_index(exit_point, dims)
        if (.not. valid(i_exit)) error stop "saddle3qs: the exit point does not exist"
    else
        if (exit_iq_min <= 0) exit_iq_min = dims(1)
        do i = 1, np
            call index_to_sub(i, dims, sub)
            mask(i) = valid(i) .and. sub(1) >= exit_iq_min
        end do
        if (.not. any(mask)) error stop "saddle3qs: no exit point"
        i_exit = minloc(E, dim=1, mask=mask)
    end if
    call show(" entry", i_entry, E_raw(i_entry))
    call show(" exit ", i_exit, E_raw(i_exit))

    !--------------------------- saddle(s)
    call find_saddle(E, valid, dims, i_entry, i_exit, Es, i_sad, path, n_path)
    if (i_sad == 0) error stop "saddle3qs: the exit cannot be reached from the entry"
    print '(a)', " --- lowest saddle between entry and exit (immersion)"
    call show(" saddle", i_sad, E_raw(i_sad))
    print '(a, f10.4, a)', "   barrier E_saddle - E_entry = ", E_raw(i_sad) - E_raw(i_entry), " MeV"
    if (A > 0) print '(a, f7.2, a, f7.2)', "   M_H / M_L (Eq. 13) = ", 0.5 * A * (1.0 + abs(ag(i_sad))), &
        " / ", 0.5 * A * (1.0 - abs(ag(i_sad)))

    if (n_saddles > 1) then
        allocate(Es_list(n_saddles), is_list(n_saddles))
        call find_saddles_with_dams(E, valid, dims, i_entry, i_exit, n_saddles, Es_list, is_list, nf)
        print '(a, i0, a)', " --- dam method: ", nf, " saddle(s)"
        do j = 1, nf
            write(line, '(a, i0)') " saddle ", j
            call show(trim(line), is_list(j), E_raw(is_list(j)))
        end do
    end if

    open(newunit=u, file=trim(path_file), status='replace', action='write')
    write(u, '(a)') "# path of the water from the entry to the exit through the lowest saddle"
    write(u, '(a)') "# step  iq id ie1 ie2 ia   Q2[b]   d[fm]    d/R  eps_f1 eps_f2 alpha_g   E[MeV]"
    do j = 1, n_path
        call index_to_sub(path(j), dims, sub)
        write(u, '(i6, 5i4, f9.3, 2f8.3, 3f7.3, f10.4)') j, sub, Q2(path(j)), d(path(j)), dR(path(j)), &
            e1(path(j)), e2(path(j)), ag(path(j)), E_raw(path(j))
    end do
    close(u)
    print '(2a)', " path: ", trim(path_file)

contains

    subroutine show(label, i, Ei)
        character(*), intent(in) :: label
        integer, intent(in) :: i
        real(dp), intent(in) :: Ei
        integer :: s(5)
        call index_to_sub(i, dims, s)
        print '(2a, 5i4, a, f8.3, a, f7.3, a, 2f7.3, a, f6.3, a, f10.4, a)', label, "  (iq id ie1 ie2 ia) =", s, &
            "  Q2 =", Q2(i), " b  d/R =", dR(i), "  eps_f =", e1(i), e2(i), "  alpha_g =", ag(i), &
            "  E =", Ei, " MeV"
    end subroutine show

end program saddle3qs
