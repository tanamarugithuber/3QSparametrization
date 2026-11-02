program test_immersion
    !---------------------------
    ! immersion_mod on a 2D landscape with two minima joined by two passes of different height
    !   E(x, y) = (x^2 - 1)^2 + 0.5 y^2 + 2 exp(-(x^2 + y^2)/0.3) + 0.3 y
    ! (entry at x = -1, exit at x = +1; the passes are at y < 0 (lower) and y > 0).
    !   1. the priority flood gives the same saddle as Moller's literal water filling
    !      (level raised in steps, scans repeated until nothing changes; PRC 79 Sec. III.C)
    !   2. the saddle is the highest point of the returned path and lies at y < 0
    !   3. the dam method (Sec. III.E) then finds the second pass at y > 0
    !   4. 5D: a separable landscape where the saddle is known
    !---------------------------
    use iso_fortran_env, only: real64
    use immersion_mod
    implicit none
    integer, parameter :: dp = real64
    integer, parameter :: nx = 81, ny = 61
    real(dp) :: E(nx * ny), x, y, Es, Ew, Es2(3)
    logical :: valid(nx * ny)
    integer :: i, j, i_entry, i_exit, i_sad, sub(2), n_fail, n_path, is2(3), nf, sub5(5)
    integer, allocatable :: path(:)
    integer :: dims5(5)
    real(dp), allocatable :: E5(:)
    logical, allocatable :: v5(:)

    n_fail = 0
    valid = .true.
    do j = 1, ny
        do i = 1, nx
            x = -2.0_dp + 4.0_dp * (i - 1) / (nx - 1)
            y = -1.5_dp + 3.0_dp * (j - 1) / (ny - 1)
            E(i + (j - 1) * nx) = (x * x - 1.0_dp)**2 + 0.5_dp * y * y + 2.0_dp * exp(-(x * x + y * y) / 0.3_dp) &
                                  + 0.3_dp * y
        end do
    end do
    call add_noise(E, valid, 1.0e-6_dp, 12345)
    i_entry = minloc(E(1:nx * ny), dim=1, mask=[(index_x(i) < nx / 2, i = 1, nx * ny)])
    i_exit = minloc(E(1:nx * ny), dim=1, mask=[(index_x(i) > nx / 2 + 1, i = 1, nx * ny)])

    call find_saddle(E, valid, [nx, ny], i_entry, i_exit, Es, i_sad, path, n_path)
    call index_to_sub(i_sad, [nx, ny], sub)
    print '(a, f10.6, a, 2f8.4)', " flood : saddle E =", Es, "  at (x, y) =", -2.0_dp + 4.0_dp * (sub(1) - 1) / (nx - 1), &
          -1.5_dp + 3.0_dp * (sub(2) - 1) / (ny - 1)
    Ew = water_filling(1.0e-7_dp)
    print '(a, f10.6)', " water : saddle E =", Ew
    call check("flood - water filling", Es - Ew, 2.0e-7_dp)
    call check("E(saddle point) - saddle energy", E(i_sad) - Es, 0.0_dp)
    call check("max E on the path - saddle", maxval(E(path)) - Es, 0.0_dp)
    if (.not. (-1.5_dp + 3.0_dp * (sub(2) - 1) / (ny - 1) < 0.0_dp)) then
        print *, "  FAIL: the lowest pass should be at y < 0"
        n_fail = n_fail + 1
    end if
    if (path(1) /= i_entry .or. path(n_path) /= i_exit) then
        print *, "  FAIL: path end points"
        n_fail = n_fail + 1
    end if

    call find_saddles_with_dams(E, valid, [nx, ny], i_entry, i_exit, 2, Es2, is2, nf)
    call index_to_sub(is2(2), [nx, ny], sub)
    print '(a, 2f10.6, a, 2f8.4)', " dams  : saddles E =", Es2(1:2), "  second at (x, y) =", &
          -2.0_dp + 4.0_dp * (sub(1) - 1) / (nx - 1), -1.5_dp + 3.0_dp * (sub(2) - 1) / (ny - 1)
    if (nf /= 2) then
        print *, "  FAIL: second saddle not found"
        n_fail = n_fail + 1
    else if (.not. (-1.5_dp + 3.0_dp * (sub(2) - 1) / (ny - 1) > 0.0_dp .and. Es2(2) > Es2(1))) then
        print *, "  FAIL: second saddle should be the higher pass at y > 0"
        n_fail = n_fail + 1
    end if

    ! 5D: E = (x1^2 - 1)^2 + sum_k c_k (x_k - s_k)^2; saddle at x1 = 0, x_k = nearest grid point to s_k
    dims5 = [21, 5, 6, 7, 4]
    allocate(E5(product(dims5)), v5(product(dims5)))
    v5 = .true.
    do i = 1, size(E5)
        call index_to_sub(i, dims5, sub5)
        x = -2.0_dp + 4.0_dp * (sub5(1) - 1) / (dims5(1) - 1)
        E5(i) = (x * x - 1.0_dp)**2 + 0.3_dp * (sub5(2) - 2)**2 + 0.2_dp * (sub5(3) - 4)**2 &
              + 0.1_dp * (sub5(4) - 7)**2 + 0.4_dp * (sub5(5) - 1)**2 + 0.01_dp * sub5(2) * x
    end do
    v5(sub_to_index([11, 2, 4, 7, 2], dims5)) = .false.   ! an unphysical point
    call add_noise(E5, v5, 1.0e-6_dp, 7)
    i_entry = sub_to_index([6, 2, 4, 7, 1], dims5)
    i_exit = sub_to_index([16, 2, 4, 7, 1], dims5)
    call find_saddle(E5, v5, dims5, i_entry, i_exit, Es, i_sad)
    call index_to_sub(i_sad, dims5, sub5)
    print '(a, f10.6, a, 5i3)', " 5D    : saddle E =", Es, "  at", sub5
    if (any(sub5 /= [11, 2, 4, 7, 1])) then
        print *, "  FAIL: 5D saddle position"
        n_fail = n_fail + 1
    end if

    if (n_fail == 0) then
        print '(a)', " all tests passed"
    else
        print '(a, i0, a)', " ", n_fail, " test(s) FAILED"
        error stop 1
    end if

contains

    pure integer function index_x(k)
        integer, intent(in) :: k
        index_x = mod(k - 1, nx) + 1
    end function index_x

    subroutine check(name, val, tol)
        character(*), intent(in) :: name
        real(dp), intent(in) :: val, tol
        if (abs(val) > tol .or. val /= val) then
            print '(a, a, es12.3, a, es9.1)', "   FAIL ", name, val, "  tol", tol
            n_fail = n_fail + 1
        end if
    end subroutine check

    function water_filling(tol) result(level)
        ! Moller's procedure: raise the level in steps (1, 0.1, 0.01, ...), at each level repeat
        ! the scan "dry point below the level with a wet neighbour -> wet" until nothing changes
        real(dp), intent(in) :: tol
        real(dp) :: level, step
        logical :: wet(nx * ny), changed
        integer :: k, nb(4), n_nb, m
        step = 1.0_dp
        level = E(i_entry)
        do while (step > tol)
            do
                wet = .false.
                wet(i_entry) = .true.
                do
                    changed = .false.
                    do k = 1, nx * ny
                        if (wet(k) .or. E(k) > level + step) cycle
                        call neighbours(k, [nx, ny], nb, n_nb)
                        do m = 1, n_nb
                            if (wet(nb(m))) then
                                wet(k) = .true.
                                changed = .true.
                                exit
                            end if
                        end do
                    end do
                    if (.not. changed) exit
                end do
                if (wet(i_exit)) exit
                level = level + step
            end do
            step = 0.1_dp * step
        end do
        level = level + 10.0_dp * step
    end function water_filling

end program test_immersion
