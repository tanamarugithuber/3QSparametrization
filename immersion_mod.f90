module immersion_mod
    !---------------------------
    ! Saddle points on a potential-energy grid by immersion
    !   P. Moller et al., Phys. Rev. C 79 (2009) 064304, Secs. III.B-III.E
    !
    ! The energy is given on a mesh of nd coordinates (stored column-major in E(:)); points
    ! with valid = .false. ("unphysical" shapes) are never wet. Nearest neighbours are the
    ! 2 nd points that differ by one step in one coordinate.
    !
    ! Moller's procedure fills the entry point (ground state, isomer, ...) with water in
    ! steps of the water level and marks a point wet when it is below the level and has a wet
    ! neighbour; the saddle height is the level at which the exit point becomes wet, and
    ! with ever smaller steps the last point to become wet before the exit is the saddle.
    ! Taking the level steps to zero, the points become wet in the order of
    !   L(j) = min over paths entry -> j of (max E along the path),
    ! which is what is computed here, exactly, by a priority flood (a Dijkstra search with
    ! "max" instead of "+"): the exit gets wet at L(exit) = the saddle energy, and the saddle
    ! point is the highest point on the path through which the water reached the exit.
    ! As in the paper, a small random number (default 0.001 R MeV, R uniform in [-1, 1])
    ! can be added to every energy so that the saddle point is unique.
    !
    ! Additional saddles (Sec. III.E, "dam"): the energy at the saddle point is raised to
    ! 1000 MeV and the search repeated. If the new saddle is next to the dam, the water
    ! only flows around it: that point is added to the dam, until the saddle found is not
    ! adjacent to the dam; it is the next saddle.
    !---------------------------
    ! 【日本語での説明】
    ! 格子上のポテンシャル曲面の鞍点を immersion 法 (水を満たす方法, PRC 79 III 章) で求める.
    ! 入口 (基底状態など) に水を入れて水位を少しずつ上げ, 水位より低く隣が濡れている点を
    ! 濡らしていく. 出口が濡れた水位が鞍点の高さ, 最後に濡れた (出口へ水が流れ込んだ)
    ! 点が鞍点. 水位の刻みを 0 にした極限では, 点 j が濡れる水位は
    !   L(j) = min_{入口から j への道} (道の上の E の最大値)
    ! になる. これを優先度つきキュー (ダイクストラ法の + を max にしたもの) で正確に求める.
    ! 鞍点は, 水が出口へ届いた道の上で E が最大の点. 論文と同じく, 鞍点を一意にするため
    ! 各点に小さな乱数 (既定 0.001 R MeV) を足せる.
    ! 【ダム】見つけた鞍点のエネルギーを 1000 MeV にして探し直す. 新しい鞍点がダムの隣なら
    ! 水がダムを回り込んでいるだけなので, その点もダムにする. ダムの隣でない鞍点が出たら,
    ! それが次の鞍点 (III.E 節).
    !---------------------------
    use iso_fortran_env, only: real64
    implicit none
    private

    integer, parameter :: dp = real64
    real(dp), parameter, public :: dam_height = 1000.0_dp

    public :: flood_levels
    public :: find_saddle
    public :: find_saddles_with_dams
    public :: index_to_sub
    public :: sub_to_index
    public :: neighbours
    public :: add_noise

    contains

        pure subroutine index_to_sub(i, dims, sub)
            ! flat index (1-based, column-major) -> subscripts
            integer, intent(in) :: i, dims(:)
            integer, intent(out) :: sub(size(dims))
            integer :: k, r
            r = i - 1
            do k = 1, size(dims)
                sub(k) = mod(r, dims(k)) + 1
                r = r / dims(k)
            end do
        end subroutine index_to_sub

        pure function sub_to_index(sub, dims) result(i)
            integer, intent(in) :: sub(:), dims(:)
            integer :: i, k, stride
            i = 1
            stride = 1
            do k = 1, size(dims)
                i = i + (sub(k) - 1) * stride
                stride = stride * dims(k)
            end do
        end function sub_to_index

        pure subroutine neighbours(i, dims, nb, n_nb)
            ! the up to 2 nd nearest neighbours of point i
            integer, intent(in) :: i, dims(:)
            integer, intent(out) :: nb(2 * size(dims)), n_nb
            integer :: k, stride, sub(size(dims))
            call index_to_sub(i, dims, sub)
            n_nb = 0
            stride = 1
            do k = 1, size(dims)
                if (sub(k) > 1) then
                    n_nb = n_nb + 1
                    nb(n_nb) = i - stride
                end if
                if (sub(k) < dims(k)) then
                    n_nb = n_nb + 1
                    nb(n_nb) = i + stride
                end if
                stride = stride * dims(k)
            end do
        end subroutine neighbours

        subroutine add_noise(E, valid, amplitude, seed)
            ! E = E + amplitude * R, R uniform in [-1, 1] (PRC 79 Sec. III.C)
            real(dp), intent(inout) :: E(:)
            logical, intent(in) :: valid(:)
            real(dp), intent(in) :: amplitude
            integer, intent(in) :: seed
            integer :: n, k
            integer, allocatable :: s(:)
            real(dp), allocatable :: r(:)
            call random_seed(size=n)
            allocate(s(n))
            s = [(seed + 37 * k, k = 1, n)]
            call random_seed(put=s)
            allocate(r(size(E)))
            call random_number(r)
            where (valid) E = E + amplitude * (2.0_dp * r - 1.0_dp)
        end subroutine add_noise

        subroutine flood_levels(E, valid, dims, i_entry, L, pred, i_stop)
            !---------------------------
            ! L(j): water level at which j becomes wet when water is poured in at i_entry,
            ! pred(j): the neighbour from which j got wet (0 for the entry and dry points).
            ! If i_stop is given, the flood stops when i_stop gets wet (L and pred are then
            ! final for every point with L(j) < L(i_stop)); points not reached have L = huge.
            !---------------------------
            ! 【目的】入口 i_entry から水を入れたとき点 j が濡れる水位 L(j) と, 水が来た隣の点 pred(j).
            real(dp), intent(in) :: E(:)
            logical, intent(in) :: valid(:)
            integer, intent(in) :: dims(:), i_entry
            real(dp), intent(out) :: L(size(E))
            integer, intent(out) :: pred(size(E))
            integer, intent(in), optional :: i_stop
            real(dp), allocatable :: hk(:)
            integer, allocatable :: hi(:)
            logical, allocatable :: done(:)
            integer :: nh, cur, nb(2 * size(dims)), n_nb, k, j
            real(dp) :: lev, new

            L = huge(1.0_dp)
            pred = 0
            allocate(done(size(E)), hk(16), hi(16))
            done = .false.
            nh = 0
            L(i_entry) = E(i_entry)
            call push(E(i_entry), i_entry)
            do while (nh > 0)
                call pop(lev, cur)
                if (done(cur)) cycle
                if (lev > L(cur)) cycle
                done(cur) = .true.
                if (present(i_stop)) then
                    if (cur == i_stop) exit
                end if
                call neighbours(cur, dims, nb, n_nb)
                do k = 1, n_nb
                    j = nb(k)
                    if (done(j) .or. .not. valid(j)) cycle
                    new = max(lev, E(j))
                    if (new < L(j)) then
                        L(j) = new
                        pred(j) = cur
                        call push(new, j)
                    end if
                end do
            end do

            contains

                subroutine push(key, idx)
                    ! binary min-heap on hk, hi
                    real(dp), intent(in) :: key
                    integer, intent(in) :: idx
                    integer :: c, p
                    real(dp), allocatable :: tk(:)
                    integer, allocatable :: ti(:)
                    if (nh == size(hk)) then
                        allocate(tk(2 * nh), ti(2 * nh))
                        tk(1:nh) = hk
                        ti(1:nh) = hi
                        call move_alloc(tk, hk)
                        call move_alloc(ti, hi)
                    end if
                    nh = nh + 1
                    c = nh
                    do while (c > 1)
                        p = c / 2
                        if (hk(p) <= key) exit
                        hk(c) = hk(p)
                        hi(c) = hi(p)
                        c = p
                    end do
                    hk(c) = key
                    hi(c) = idx
                end subroutine push

                subroutine pop(key, idx)
                    real(dp), intent(out) :: key
                    integer, intent(out) :: idx
                    real(dp) :: lk
                    integer :: li, c, ch
                    key = hk(1)
                    idx = hi(1)
                    lk = hk(nh)
                    li = hi(nh)
                    nh = nh - 1
                    c = 1
                    do
                        ch = 2 * c
                        if (ch > nh) exit
                        if (ch < nh) then
                            if (hk(ch + 1) < hk(ch)) ch = ch + 1
                        end if
                        if (lk <= hk(ch)) exit
                        hk(c) = hk(ch)
                        hi(c) = hi(ch)
                        c = ch
                    end do
                    if (nh > 0) then
                        hk(c) = lk
                        hi(c) = li
                    end if
                end subroutine pop
        end subroutine flood_levels

        subroutine find_saddle(E, valid, dims, i_entry, i_exit, E_saddle, i_saddle, path, n_path)
            !---------------------------
            ! saddle between i_entry and i_exit: E_saddle = L(i_exit), i_saddle = the highest
            ! point on the path along which the water reached i_exit; path(1:n_path) is that
            ! path from the entry to the exit. i_saddle = 0 if the exit cannot be reached.
            !---------------------------
            ! 【目的】入口と出口の間の鞍点の高さ E_saddle と位置 i_saddle, 水の通り道 path.
            real(dp), intent(in) :: E(:)
            logical, intent(in) :: valid(:)
            integer, intent(in) :: dims(:), i_entry, i_exit
            real(dp), intent(out) :: E_saddle
            integer, intent(out) :: i_saddle
            integer, allocatable, intent(out), optional :: path(:)
            integer, intent(out), optional :: n_path
            real(dp), allocatable :: L(:)
            integer, allocatable :: pred(:), tmp(:)
            integer :: j, n

            allocate(L(size(E)), pred(size(E)))
            call flood_levels(E, valid, dims, i_entry, L, pred, i_exit)
            i_saddle = 0
            E_saddle = huge(1.0_dp)
            if (present(n_path)) n_path = 0
            if (L(i_exit) == huge(1.0_dp)) return
            E_saddle = L(i_exit)
            ! walk back from the exit
            allocate(tmp(size(E)))
            n = 0
            j = i_exit
            do while (j /= 0)
                n = n + 1
                tmp(n) = j
                if (j == i_entry) exit
                j = pred(j)
            end do
            i_saddle = tmp(1)
            do j = 1, n
                if (E(tmp(j)) > E(i_saddle)) i_saddle = tmp(j)
            end do
            if (present(path)) then
                allocate(path(n))
                path = tmp(n:1:-1)
            end if
            if (present(n_path)) n_path = n
        end subroutine find_saddle

        subroutine find_saddles_with_dams(E, valid, dims, i_entry, i_exit, n_want, E_s, i_s, n_found)
            !---------------------------
            ! the lowest saddle and up to n_want - 1 further saddles by the dam method
            ! (PRC 79 Sec. III.E). E is not changed.
            !---------------------------
            ! 【目的】一番低い鞍点と, ダム法で見つかる次の鞍点を合わせて最大 n_want 個.
            real(dp), intent(in) :: E(:)
            logical, intent(in) :: valid(:)
            integer, intent(in) :: dims(:), i_entry, i_exit, n_want
            real(dp), intent(out) :: E_s(n_want)
            integer, intent(out) :: i_s(n_want), n_found
            real(dp), allocatable :: Ed(:)
            logical, allocatable :: dam(:)
            integer :: nb(2 * size(dims)), n_nb, k, i_sad, it
            real(dp) :: Es
            logical :: next_to_dam

            allocate(Ed(size(E)), dam(size(E)))
            Ed = E
            dam = .false.
            n_found = 0
            E_s = 0.0_dp
            i_s = 0
            call find_saddle(Ed, valid, dims, i_entry, i_exit, Es, i_sad)
            if (i_sad == 0) return
            n_found = 1
            E_s(1) = Es
            i_s(1) = i_sad
            do while (n_found < n_want)
                ! dam at the last saddle, extended while the water only flows around it
                dam(i_s(n_found)) = .true.
                Ed(i_s(n_found)) = dam_height
                do it = 1, size(E)
                    call find_saddle(Ed, valid, dims, i_entry, i_exit, Es, i_sad)
                    if (i_sad == 0) return
                    if (Es >= dam_height) return
                    call neighbours(i_sad, dims, nb, n_nb)
                    next_to_dam = .false.
                    do k = 1, n_nb
                        if (dam(nb(k))) next_to_dam = .true.
                    end do
                    if (.not. next_to_dam) exit
                    dam(i_sad) = .true.
                    Ed(i_sad) = dam_height
                end do
                n_found = n_found + 1
                E_s(n_found) = Es
                i_s(n_found) = i_sad
            end do
        end subroutine find_saddles_with_dams

end module immersion_mod
