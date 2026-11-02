  # MACROSCOPIC-MICROSCOPIC法

## 概要

このプロジェクトは、原子核物理学におけるFRLDM+Shell+BCSを解くためのプログラム群です．

## 3QS 形状による核分裂の鞍点

P. Möller et al., Phys. Rev. C 79, 064304 (2009) の三つの 2 次曲面（3QS）の形で FRLDM（と殻補正・対補正）のポテンシャル曲面を 5 次元の格子 (Q2, d, εf1, εf2, αg) で計算し，immersion 法で核分裂の鞍点を求めるプログラム（pes3qs.exe，saddle3qs.exe）があります．使い方と中身は [README_3QS.md](README_3QS.md) を見てください．

```bash
bash build_3qs.sh
./pes3qs.exe pes3qs.nml && ./saddle3qs.exe saddle3qs.nml
```

## オブジェクト指向プログラミングに慣れていない人向け

このコードは Fortran のオブジェクト指向の書き方（派生型，型束縛手続き `call frldm%calculate_frldm_energy(...)` など）を使っています．慣れていない人は，先に解説書を読んでください．

- [Fortran で学ぶオブジェクト指向プログラミング（PDF）](docs/fortran_oop_guide.pdf)

モジュールとカプセル化，派生型と型束縛手続き，継承，多態性（`class`，`select type`），抽象型と `deferred`，演算子の多重定義を，すぐに動かせる短い例で説明しています．液滴模型項・殻補正・対相関を抽象型で部品化する例（簡略化したトイモデル）や，複数ファイルのコンパイル順，よくある落とし穴もまとめています．

## 使い方 (main.f90)

main は 1 つの原子核について，FRLDM の巨視的エネルギーと，Strutinsky 法のシェル補正 + Lipkin–Nogami（または BCS）の対補正を計算し，質量超過と結合エネルギーを出力します．

- M_th = E_mac + E_shell + E_pair + E_zp（質量超過，MeV）
- B = Z M_H + N M_n − M_th（結合エネルギー）
- 形は FRDM の ε パラメータ (ε2, γ, ε3, ε4, ε6) です．ε1 は重心が原点に来るように ε3 から自動で決まります．
- `search = .true.` にすると E_mac + E_s+p が最小になる形を探し，その形で結果を出します．
  - `search_method = 'pattern'`（既定）：(ε2, ε3, ε4, ε6) の刻み `d_eps`（既定 0.01）の格子上で，複数の初期点からパターン探索をします（FRDM2012 3 章の手順 3 と同じ離散最小化）．今の点から各軸 ±s の 8 近傍を計算して一番低い近傍へ移り，下がらなくなったら s を `d_eps_start`（既定 0.04）→ 0.02 → 0.01 と半分にします．計算した点は覚えておき，計算し直しません．初期点は入力の形，(ε2, ε4) = (±0.25, ±0.04) と球形で，それぞれ ε3 = 0 と 0.10 です（ε6 = 0）．E(−ε3) = E(ε3) なので ε3 は 0 以上だけ調べます．γ は `gamma_deg` に固定です．
  - 探索範囲は `eps2_min, eps2_max, eps3_max, eps4_min, eps4_max, eps6_min, eps6_max` です．min = max にするとその変数は固定になります（例：`eps3_max = 0` で ε3 = 0 のまま）．
  - `search_method = 'grid'`：これまでどおり (ε2, ε4) の格子（刻み `d_eps2, d_eps4`）を全部計算します．ε3, ε6 は入力の値に固定です．
  - `n_gh_search` を正にすると，探索中だけ Gauss–Hermite 点数をその値にします（最後の形は `n_gh` で計算し直します）．
- `zpe = .true.`（既定）なら，最後の形で零点エネルギー E_zp（FRDM2012 Sec. 2.11，zero_point_mod）を足します．ε2 を ±0.05, ±0.10, ±0.15 ずらした 6 点（ε2 ≥ 0.17 なら γ を ±5°, ±10°, ±15° ずらした 6 点も）で E_pot を計算し，7 点の最小二乗放物線からかたさ C を求めます．ε2 < 0.17 では E_0,γ = E_0,ε2 です．形を 6〜12 点余分に計算するので，時間は 7〜13 倍かかります．`zpe = .false.` で零点エネルギーなしになります．

```bash
bash build.sh        # FFTW3 と LAPACK が必要
./main.exe input.nml # input.nml がなければ 208Pb（球形）で計算します
./main.exe input_240Pu.nml
```

入力は namelist `&frldm_input` です（例：²⁴⁰Pu を FRDM2012 の基底状態の形で計算）．

```
&frldm_input
  Z = 94, N = 146,
  eps2 = 0.22, gamma_deg = 0.0, eps3 = 0.0, eps4 = -0.05, eps6 = 0.04,
  pairing = 'LN',          ! 'LN'（既定，FRDM2012）または 'BCS'
  n_gh = 44,               ! 一粒子問題の Gauss-Hermite 点数（各軸）
  use_grid = .false.,      ! .true. で B_1, B_3 を格子 + Hockney 法で計算（検算用）
  h = 0.5, n_times = 3,    ! そのときの格子間隔 [fm] と箱の大きさ
  search = .false.,        ! 基底状態の探索
  search_method = 'pattern', ! 'pattern'（(ε2, ε3, ε4, ε6) の離散最小化）または 'grid'
  eps2_min = -0.30, eps2_max = 0.45, eps3_max = 0.20,
  eps4_min = -0.12, eps4_max = 0.12, eps6_min = -0.08, eps6_max = 0.08,
  d_eps = 0.01, d_eps_start = 0.04,   ! パターン探索の最後と最初の刻み
  d_eps2 = 0.05, d_eps4 = 0.04,       ! 'grid' の刻み
  zpe = .true.             ! 零点エネルギーを足す
/
```

計算時間は 1 つの形あたり約 4.5 秒です（4 スレッド，ε3 = 0 のとき）．一粒子問題は固有値だけを求めます（固有ベクトルは使わないので計算しません．以前は固有ベクトルも求めていて約 11 秒でした）．

### FRDM2012 の表の形でまとめて計算する

FRDM2012 の Table 1（9318 核）の ε2, ε3, ε4, ε6，E_s+p，M_th，M_exp（AME2003），M_th^FL を data/frdm2012_table.dat にしてあります（tools/parse_frdm2012_table.py で thesis/FRLDM_2012.pdf から作りました）．tools/run_frdm2012_shapes.py は，この表の形（γ = 0）で main.exe を核ごとに走らせ，E_mac，E_s+p，E_zp，M_th を表の値と並べて書き出します．

```bash
bash build.sh
python3 tools/run_frdm2012_shapes.py --zmax 20 --measured --jobs 4   # Z ≤ 20 で実験値のある核
python3 tools/run_frdm2012_shapes.py --zmin 82 --zmax 82 --no-zpe   # Pb 同位体，零点エネルギーなし
```

- 結果は data/frldm_at_frdm2012_shapes.dat（`--out` で変更）に 1 行 1 核で追記されます．列は Z, N, A, ε2, ε3, ε4, ε6, E_mac, E_s+p, E_zp, M_th（本コード），M_th^FL（表），M_th（FRDM2012），M_exp，M_th − M_th^FL，M_exp − M_th，flags です．
- flags は，LN の解がなく BCS を使った（`LN->BCS`），零点エネルギーのかたさが負（`C_eps2<=0`，`C_gamma<=0`），main が異常終了した（`failed`）のしるしです．
- 入力ファイルと main の出力は `--workdir`（既定 runs_frdm2012/）に残ります．結果ファイルにすでにある核は飛ばすので，途中で止めても同じコマンドで続きから計算できます．
- `--jobs` 個の核を同時に計算し，各核は `--threads`（既定 1）スレッドです．零点エネルギーありで 1 核あたり 7〜13 形を計算します．

### FRDM2012 との比較

FRDM2012 の表の値と比べました（E_mac の表の値は M_th^FL − E_mic^FL）．

| 核 | 形 (ε2, ε4, ε6) | E_mac 本コード | E_mac 表 | E_s+p 本コード | E_s+p 表 |
|---|---|---|---|---|---|
| ²⁰⁸Pb | (0, 0, 0) | −8.70 | −8.72 | −13.94 | −13.93 |
| ¹³²Sn | (0, 0, 0) | −64.97 | −64.98 | −12.83 | −12.82 |
| ²⁴⁰Pu | (0.22, −0.05, 0.04) | — | — | −4.93 | −4.95 |
| ⁴⁰Ca | (0, 0, 0) | −37.51 | −37.52 | −0.52 | −0.52 |
| ⁴⁸Ca | (0, 0, 0) | −43.34 | −43.34 | −1.95 | −1.94 |

²⁴⁰Pu は表の E_mac が球形の値なので比べていません．零点エネルギーを入れると，²⁰⁸Pb の M_th は −21.26 MeV（表の M_th^FL は −21.31），²⁴⁰Pu は 50.90 MeV（表は 50.69）になります．

## 参考文献

[Macro-Micro Frame Work](https://www.sciencedirect.com/science/article/pii/S0092640X1600005X)  
[Axial Asymmetry](https://iopscience.iop.org/article/10.1088/0031-8949/8/1-2/003)  
[Shell+Pairing](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.5.1050)  
[Pairing](https://www.sciencedirect.com/science/article/pii/037594749290244E)  
[Shell1](https://www.sciencedirect.com/science/article/pii/0375947468906994)  
[Shell2](https://www.sciencedirect.com/science/article/pii/0375947467905106)  
[Nilsson Oscillator](https://www.sciencedirect.com/science/article/pii/0375947469908094)




## 各ファイルの依存関係
A --> B は「B が A を use している」ことを表します（各ファイルの `use` 文から作成）．間接的に参照しているものは省略しています（例：constant_mod.f90 はほぼすべてのファイルが use しています）．

### FRLDM+Shell+BCS で使うファイル

```mermaid
graph TD;
    constant_mod.f90-->nucleus_mod.f90;
    constant_mod.f90-->shape_mod.f90;
    nucleus_mod.f90-->grid_mod.f90;
    shape_mod.f90-->grid_mod.f90;
    grid_mod.f90-->frldm_mod.f90;
    shape_mod.f90-->frldm_mod.f90;
    calculate_potential_hockney.f90-->frldm_mod.f90;
    constant_mod.f90-->calculate_potential_hockney.f90;
    constant_mod.f90-->ho_basis_mod.f90;
    shape_mod.f90-->potential_gh_mod.f90;
    nucleus_mod.f90-->micro_constant_mod.f90;
    micro_constant_mod.f90-->shell_bcs_mod.f90;
    ho_basis_mod.f90-->shell_bcs_mod.f90;
    potential_gh_mod.f90-->shell_bcs_mod.f90;
    nucleus_mod.f90-->zero_point_mod.f90;
    micro_constant_mod.f90-->zero_point_mod.f90;
    frldm_mod.f90-->main.f90;
    shell_bcs_mod.f90-->main.f90;
    zero_point_mod.f90-->main.f90;
```

### 3QS の鞍点で使うファイル

```mermaid
graph TD;
    shape_mod.f90-->qs3_shape_mod.f90;
    qs3_shape_mod.f90-->qs3_coord_mod.f90;
    qs3_shape_mod.f90-->frldm_mod.f90;
    qs3_shape_mod.f90-->shell_bcs_mod.f90;
    qs3_coord_mod.f90-->pes3qs.f90;
    frldm_mod.f90-->pes3qs.f90;
    shell_bcs_mod.f90-->pes3qs.f90;
    immersion_mod.f90-->saddle3qs.f90;
```

- qs3_shape_mod：3QS 形状（surface_shape を拡張）．qs3_coord_mod：PRC 79 の座標 (Q2, d, εf1, εf2, αg) から 3QS の σ, α．
- frldm_mod は 3QS の形で B_1, B_3 を軸対称のリング和で計算し，B_W に FRDM2012 式 (75) を使います．shell_bcs_mod は `qs3=` を渡すと 3QS の形のポテンシャルを使います．
- immersion_mod：immersion 法（ダム法つき）．

- 形状：shape_mod（FRDM の ε パラメータの形）を FRLDM と一粒子ポテンシャルの両方で使います．
- FRLDM：frldm_mod が B_1, B_3 を形の表面積分で計算します（`use_grid` のときは grid_mod + calculate_potential_hockney）．
- Shell+pairing：ho_basis_mod（三軸変形 HO 基底）と potential_gh_mod（求積点上の V1, VC）で一粒子準位を求め，shell_bcs_mod でシェル補正と対補正を計算します．
- 零点エネルギー：zero_point_mod が慣性とかたさから E_zp を計算します．形をずらした E_pot は main の `add_zero_point` が集めます．
- main がすべてを呼んで合計します．

### テストプログラム

```mermaid
graph TD;
    shape_mod.f90-->test_shape.f90;
    frldm_mod.f90-->test_frldm.f90;
    grid_mod.f90-->test_grid.f90;
    ho_basis_mod.f90-->test_ho_basis.f90;
    ho_basis_mod.f90-->test_potential_gh.f90;
    potential_gh_mod.f90-->test_potential_gh.f90;
    calculate_potential_hockney.f90-->test_hockney.f90;
    shell_bcs_mod.f90-->test_shell_bcs.f90;
    zero_point_mod.f90-->test_zero_point.f90;
    qs3_coord_mod.f90-->test_qs3.f90;
    frldm_mod.f90-->test_frldm_3qs.f90;
    immersion_mod.f90-->test_immersion.f90;
```

### 使わないファイル (not_use/)

FRLDM+Shell+BCS の計算に使わないファイルは not_use/ にまとめています（ビルドには含めません）．

- 旧方式（格子 + Lanczos / 虚時間）：Micro_potential_mod.f90，EV_Lanczos.f90，Imarginary_time_mod.f90，calculate_potential_mod.f90，fundamental/（CG_method_mod.f90，3Dderivative_mod.f90，3Dderivative_cp_mod.f90，calculate_vector_cp.f90，common_func_mod.f90，FFT_integral_mod.f90）
- FALR 法のクーロン計算：calculate_coulomb_FFT.f90，test_coulomb_fft.f90（Hockney 法に置き換えました）
- その他：io_mod.f90（空），data.text（以前の出力）

## 原子核の形 (shape_mod.f90)

FRDM2012 の ε パラメータ (ε2, γ, ε3, ε4, ε6)（Eqs. (17)–(25)）の形を，伸縮座標での Nilsson 等ポテンシャル面として作ります．

- `epsilon_shape_from(R, eps2, gamma, eps3, eps4, eps6)`（γ はラジアン）で，体積 4πR³/3 の形を返します．
- ε3 ≠ 0 のときは，重心の条件 ∫cosθ_t F⁻² dΩ = 0 を二分法で解いて ε1 を決めます（FRDM2012 p. 35 の Notes added in proof）．結果は `shape%eps1` に入ります．詳しくは「[八重極変形 ε3 と重心の ε1](#八重極変形-ε3-と重心の-ε1)」を見てください．
- どの形も x, y の鏡映について対称です．ε3 は z の鏡映対称性を破るので，体積・表面積分は 1/4 領域（x, y ≥ 0）で行います．
- `eval(u, v, r, N)` で表面の点と法線 N = ∂r/∂u × ∂r/∂v，`inside(x)` で内外判定，`relative_area()` で B_s（Eq. 70），`volume()` で体積を返します．
- ε4 = ε6 = 0 では `ellipsoid_from_eps_gamma` の楕円体と一致します．²⁴⁰Pu（ε2 = 0.22, ε4 = −0.05, ε6 = 0.04）の β2, β4, β6 = 0.237, 0.086, −0.024 は FRDM2012 の表と一致します．

```bash
gfortran -O2 constant_mod.f90 shape_mod.f90 test_shape.f90 -o test_shape.exe
./test_shape.exe
```

## FRLDM (frldm_mod.f90)

`call frldm%calculate_frldm_energy(nucleus, shape, E_mac)` で FRDM2012 の Eq. (62) を計算します．E_mac は質量超過としての値で，M_H Z + M_n N を含みます．各項は `frldm%E_v`，`frldm%E_s` などに入ります．

- B_1, B_3 は二重体積積分を Gauss の定理で二重表面積分 −∮∮(dS·dS′) g(|r−r′|) に直して計算します（∇²g が元の核）．B_1 は g = −a² e^{−σ/a}，B_3 は g = a[s/2 − e^{−s}/2 + 2(1 − e^{−s})/s]（s = σ/a）です．
- B_3 の被積分関数は滑らかで，球で 1e−10 の精度です．B_1 は σ = 0 の折れ目のため誤差が n^{−3} で減るので，n と 3n/2 の 2 点から Richardson 外挿します（球で 3e−5）．1 回の計算は約 2 秒です．
- 格子 + Hockney 法（`g_mod` を渡したとき）も残してあります．表面のセルを体積比で扱うので B_3 は格子間隔 0.4 fm で 1e−3 まで合いますが，B_1 は表面がぼやける分だけ O(h) 小さく出ます（¹⁶O，h = 0.5 fm で 5%）．
- B_W = 1（基底状態の一体形状）としています．

```bash
gfortran -O2 -fopenmp -I/usr/include constant_mod.f90 nucleus_mod.f90 shape_mod.f90 qs3_shape_mod.f90 grid_mod.f90 calculate_potential_hockney.f90 frldm_mod.f90 test_frldm.f90 -lfftw3 -o test_frldm.exe
./test_frldm.exe
```

## 三軸変形調和振動子基底 (ho_basis_mod.f90)

三軸変形核の一粒子準位を変形調和振動子基底の行列対角化で求めるモジュールです．ポテンシャル・∇V₁はGauss–Hermite求積点 `basis%point(axis,i)` 上で与えます．LAPACK (dstev, zheev) が必要です．

```bash
gfortran -O2 -fopenmp constant_mod.f90 ho_basis_mod.f90 test_ho_basis.f90 -llapack -o test_ho_basis.exe
./test_ho_basis.exe
```

## 求積点上のポテンシャル (potential_gh_mod.f90)

折りたたみ湯川ポテンシャル V₁，その勾配 ∇V₁，クーロンポテンシャル V_C を，生成形状の表面積分で任意の点（Gauss–Hermite求積点など）上に計算するモジュールです（FRDM2012 式(81),(91)）．ポアソン方程式を解かないので，無限遠の境界条件は厳密に満たされます．

- `Y = V₁/(-V₀)`，`∇Y`，`C = V_C/(e²ρ_c)` を返します．`V₁ = -V₀ Y`，`V_C = e²ρ_c C`（ρ_c = Z/(4πR_pot³/3)）です．
- 形状は shape_mod の `surface_shape` を使います（`epsilon_shape`，`ellipsoid_shape`）．
- 表面近くの点は，表面パネルを適応的に分割して計算します．∇Y の 1/s 特異性は恒等式で取り除いてあります．
- `on_grid(..., d2h=.true.)` で鏡映対称性を使い，1/8 の点だけを計算します．shell_bcs_mod は ε3 = 0 のときだけ d2h を使います．

```bash
gfortran -O2 -fopenmp constant_mod.f90 shape_mod.f90 ho_basis_mod.f90 potential_gh_mod.f90 test_potential_gh.f90 -llapack -o test_potential_gh.exe
./test_potential_gh.exe
```

スピン軌道項：∇V₁ は表面で折れ曲がるため，⟨(∂_jV₁)∂_k⟩ の形のままでは Gauss–Hermite 求積が収束せず，球形核でも (2j+1) 重縮退が最大 0.04 MeV 分裂していました．部分積分した A_jk = −∫V₁(∂_jφ_a ∂_kφ_b − ∂_kφ_a ∂_jφ_b) を使う `spin_orbit_matrices_by_parts` に変えたので，分裂は求積点の数とともに小さくなります（²⁰⁸Pb で n_gh = 44 のとき 0.01 MeV 以下，56 のとき 0.003 MeV 以下）．

## 空間格子 (grid_mod.f90)

格子の範囲は変形パラメータから決めます．`nucleus%epsilon2`，`nucleus%gamma_deg`（度），`nucleus%epsilon4`，`nucleus%epsilon6` を設定してから `initialize_grid` を呼びます．格子は main の `use_grid = .true.` のとき（B_1, B_3 の検算）だけ使います．

- 形は shape_mod の `epsilon_shape` で，`grid%shape` に入ります．`grid%semi` はその各軸方向の広がりです．
- γ = 0 は z 軸まわりの軸対称（ε > 0 で長楕円），γ = 60° は x 軸まわりの扁平形です．
- 内外判定は形の `inside` で行います．表面にかかるセルは n_sub³（既定 6³）点のうち内側の割合を値にするので，体積の誤差は 1e−4 程度です．
- `grid%cubic = .true.`（既定）では，一番長い半軸に合わせた立方体の箱にします．以前の CG_method_mod が n_x = n_y = n_z を仮定していたためです（Hockney 法は直方体の箱でも計算できます）．`.false.` にすると軸ごとに箱の大きさを決めます．

```bash
gfortran -O2 -fopenmp constant_mod.f90 nucleus_mod.f90 shape_mod.f90 grid_mod.f90 test_grid.f90 -o test_grid.exe
./test_grid.exe
```

## クーロンポテンシャル (not_use/calculate_coulomb_FFT.f90，現在は未使用)

クーロンポテンシャル V(r) = ∫ρ(r′)/|r − r′| d³r′ を，Lauritsch–Reinhard の FALR 法（thesis/FFTcoulomb.pdf）で計算します．箱の境界でディリクレ条件を課す CG_method_poisson3D の代わりに，無限遠で V → 0 という正しい境界条件になります．

1. 多重極モーメント q_lm（l ≤ 4）を求め，同じモーメントをもつ ρ_long = Σ a_lm ρ_l(r) S_lm を作ります．
2. V_long は l ごとの動径方程式をグリーン関数で解きます．
3. 残りの ρ_short = ρ − ρ_long は多重極をもたないので，FFT（V = 4πρ̃/k²）で解いても周期像の影響がほとんど残りません．k = 0 成分は −(2π/3)∫r²ρ_short で決めます．

論文は l_cut = 2 ですが，変形核は十六重極モーメントが大きいので l_cut = 4 にしています．呼び方は `call calculate_coulomb(rho, R, h, V)` です（R は参照密度の半径で，`nucleus%R_a` を渡します）．FFTW3 が必要です（Ubuntu は `sudo apt install libfftw3-dev`，Windows は MSYS2 の `pacman -S mingw-w64-x86_64-fftw`）．

```bash
gfortran -O2 -fopenmp -I/usr/include constant_mod.f90 not_use/calculate_coulomb_FFT.f90 not_use/test_coulomb_fft.f90 -lfftw3 -o test_coulomb_fft.exe
./test_coulomb_fft.exe
```

## 自由空間の畳み込み (calculate_potential_hockney.f90)

Hockney 法で，クーロン ∫ρ/|r−r′|，湯川 ∫ρ e^{−|r−r′|/a}/|r−r′|，指数 ∫ρ e^{−|r−r′|/a} の 3 種類のポテンシャルを計算します．frldm_mod の格子版の B_1, B_3（`use_grid`）はこれを使います．

- 格子を各軸 2 倍の箱に入れて残りを 0 で埋め，FFT で巡回畳み込みをとります．元の格子の 2 点間の距離は 2 倍の箱の半分を超えないので，自由空間の和と厳密に一致します．箱の境界条件は不要です．
- 核は各セルでの平均値を使います（自己セルは 1/r の立方体積分 2.3800773…/h）．このため，セルごとに一定の密度（grid_mod の内外判定 `density_index3D`）に対しては厳密です．滑らかな密度を点の値で渡すときは `smooth = .true.` にします．
- B_1 の核 (2a/r − 1)e^{−r/a} と B_3 の核 1/r − (1/r + 1/(2a))e^{−r/a} を同じ方法で畳み込むので，以前の O(N²) の二重和が不要になりました．
- メモリは格子の 8 倍です．FFTW3 が必要です．

```bash
gfortran -O2 -fopenmp -I/usr/include constant_mod.f90 calculate_potential_hockney.f90 test_hockney.f90 -lfftw3 -o test_hockney.exe
./test_hockney.exe
```

## シェル補正と対補正 (shell_bcs_mod.f90)

`call shell_bcs_calculation(nucleus, mv, res, n_gh, pairing)` で，(ε2, γ, ε4, ε6) の形について中性子・陽子のシェル補正と対補正を計算します（ε3 も使えます）（FRDM2012 Eqs. (78)–(111)，Möller–Nix 1992）．定数は micro_constant_mod のものを使います．

1. 一粒子準位：ho_basis_mod の変形 HO 基底（ħω0 = C_cur/A^{1/3}，N_bas = 12 殻）と potential_gh_mod の V1, VC（半径 R_pot の ε 形状）からハミルトニアンを作り，対角化します．スピン軌道項は部分積分形です．δ̄, ε̄, R_pot, V_p, V_n などは `mv%calculate_microscopic_variables`（Eqs. (82)–(93)）で `mv` に入ります．
2. Strutinsky 法：ガウス関数 × p 次の補正多項式（p = 8）で平滑化します．幅は γ_S = C_S ħω0 B_s です．E_shell = Σe_i − Ẽ．
3. 対補正（既定は Lipkin–Nogami，`pairing = 'BCS'` で BCS）：有効相互作用ギャップ Δ_G = r_mic B_s / N^{1/3} と平滑化した準位密度から G を決めます（Eq. 107）．
   - Lipkin–Nogami（Eqs. (94)–(101)）：λ2 についての 1 つの方程式にまとめ，二分法で解きます（内側でずらした準位 ε_k = e_k + (4λ2 − G)v_k² の BCS 方程式を解きます）．Δ → 0 で λ2 が発散するので，魔法数でも解がつぶれません．平均値は Eqs. (108)–(110) です．
   - 軽い核（⁷Li, ⁸Be など）では対の窓に縮退した準位が 2 つしかなく，前の v_k² から反復すると，ためしたどの λ2 でも解がつぶれることがあります．そのときは λ2 を [0, 2 max(⟨λ2⟩, G)] で 200 点しらべ（各点は鋭い占有数から反復），符号が変わる区間で二分法をします．それでも解がなければ警告を出して BCS（λ2 = 0）を使います．重い核では今までと同じ結果です．
   - BCS：Möller–Nix Eqs. (5), (27)．魔法数など解がないときは Δ = 0 とします．
   - 対の窓は「最後の束縛準位」か「フェルミ準位から 5 MeV 上」の高い方まで，下にも同じ数の準位をとります．奇数粒子は blocking せず，準粒子エネルギーを加えます．平均の奇数項は FRLDM と同じ Δ̄ = r_mac B_s / N^{1/3}（Eqs. (63), (64)）です．
   - r_mic = 3.2 MeV は Lipkin–Nogami 用の値です．BCS で使う場合は合わせ直す必要があります．

```bash
gfortran -O2 -fopenmp constant_mod.f90 nucleus_mod.f90 shape_mod.f90 qs3_shape_mod.f90 micro_constant_mod.f90 ho_basis_mod.f90 potential_gh_mod.f90 shell_bcs_mod.f90 test_shell_bcs.f90 -llapack -o test_shell_bcs.exe
./test_shell_bcs.exe
```

## 八重極変形 ε3 と重心の ε1

FRDM2012 の形（Eqs. (17)–(25)）に，八重極変形 ε3 と，重心を原点に置くための ε1 を入れています．

### 形の式

伸縮座標 ξ_k = √f_k x_k の極座標 (ρ_t, θ_t, φ_t) で，表面は ρ_t² F(θ_t, φ_t) = C² です．

```
F = f_x n_x² + f_y n_y² + f_z n_z² + 2ε1 P1(cosθ_t) + 2ε3 P3(cosθ_t) + 2ε4 V4 + 2ε6 P6(cosθ_t)
```

- P3 は z → −z で符号を変えるので，ε3 > 0 では +z 側がふくらみ，−z 側がとがった梨の形になります．ε3 と −ε3 は z を反転した同じ形なので，エネルギーも同じです．
- ε3 だけを入れると重心が原点からずれます．P1(cosθ_t) = cosθ_t は形をほぼ平行移動させる項なので，ε1 で重心を原点に戻します．
- C は体積が 4πR³/3 になるように決めます（Eq. 25）．

### ε1 の決め方（FRDM2012 p. 35 の Notes added in proof）

重心は ∫z dV = det(D) f_z^{−1/2} C⁴/4 · ∫cosθ_t F⁻² dΩ なので，次の式を解きます．

```
M(ε1) = ∫ cosθ_t F(θ_t, φ_t)⁻² dΩ = 0
```

- M は C と R によらないので，体積を合わせる前に解けます．
- dM/dε1 = −4∫cos²θ_t F⁻³ dΩ < 0 なので，M は単調減少で，解は 1 つです．
- `epsilon_shape_from` は [−0.02, 0.02] から区間を倍々に広げて符号が変わる区間を探し，二分法で解きます（幅 1e−14 まで）．
- ε3 = 0 なら ε1 = 0 のままです．入力では ε1 を与えません．main は求めた ε1 を「eps1 (centre of mass at the origin)」として表示します．

### 対称性と積分

- ε3 ≠ 0 の形も x, y の各面については鏡映対称ですが，z の鏡映対称性は破れます．
- そのため，体積・表面積・FRLDM の B_1, B_3（`surface_sums`）の積分は，x, y ≥ 0 の 1/4 領域（θ ∈ [0, π]，φ ∈ [0, π/2]）で行って 4 倍します．以前の 1/8 領域（z ≥ 0 も仮定）は使っていません．
- 一粒子ポテンシャル（`q%on_grid`）は，ε3 = 0 のときだけ 1/8 の点を計算して鏡映で埋めます（d2h）．ε3 ≠ 0 では全点を計算するので，その部分の時間は約 8 倍になります．
- 振動子基底は (ε2, γ) だけで決めます（ε3 はポテンシャルの側だけに入ります）．パリティの混ざりはハミルトニアンの対角化で自然に入ります．

### 確認したこと

- ²⁰⁸Pb（ε3 = 0）の M_th は ε3 を入れる前と同じでした．
- ²²⁴Ra（ε2 = 0.12，ε4 = −0.04）では，ε3 = ±0.08 で M_th が 20.28 → 18.83 MeV に下がり，+ε3 と −ε3 で同じ値になります．
- test_shape で，ε1 を入れると重心が 1e−10 fm 以内で原点に来ること（入れないと原点からずれること），ε3 → −ε3 で ε1 の符号が反転すること，1/4 領域の積分が全表面の積分と一致することを確かめています．

### 使い方

```
&frldm_input
  Z = 88, N = 136,
  eps2 = 0.12, eps3 = 0.08, eps4 = -0.04,
/
```

基底状態の探索（`search = .true.`，既定の `search_method = 'pattern'`）は ε3 も含めて (ε2, ε3, ε4, ε6) で最小化します．

## 零点エネルギー (zero_point_mod.f90，main の add_zero_point)

基底状態の形 ε_gs で求めた E_pot = E_mac + E_s+p はポテンシャル曲面の谷底の値です．形は量子力学的にゆらぐので，実際の基底状態はそこから零点振動のエネルギーだけ上にあります．FRDM2012（Sec. 2.11，Eqs. (112)–(118)）はこれを

```
M_th = E_mac(ε_gs) + E_s+p(ε_gs) + E_zp,    E_zp = E_0,ε2 + E_0,γ
```

として足します．ε2 方向（核分裂方向）と γ 方向（三軸非対称方向）の 2 つの振動を考え，ε3, ε4, ε6 方向は入れません．

### 式

基底状態のまわりで E(λ) ≈ E(λ_gs) + ½C_λ(λ − λ_gs)²（λ = ε2 または γ）とすると，1 次元の調和振動子になります．

```
E_0,λ = ½ ħω_λ,    ħω_λ = 𝒦 ħc (C_λ / B_λ^irr)^{1/2},    𝒦 = 0.2475（frK）
B_ε2^irr = (2/15) (1 + 2/9 ε2²) / (1 − 2/3 ε2)² · (1 − 1/3 ε2² − 2/27 ε2³)^{−4/3} · M0 R0²     (116)
B_γ^irr  = (2/15) [(1 − 2/3 ε2)/(1 + 1/3 ε2)]^{2/3} · [ln((1 + 1/3 ε2)/(1 − 2/3 ε2))]² · M0 R0²  (117)
```

- B^irr は非圧縮・渦なし流れの慣性です．本当の慣性はこれより大きく，その分を 1 つの定数 𝒦 で表します（B = B^irr/𝒦²）．
- M0 = A m_nuc（`nucleus%mass`），R0 = r0 A^{1/3}（`nucleus%R_a`）です．B は MeV fm² で，ħc を掛けると ħω が MeV で出ます．
- 球では B_ε2^irr = (2/15) M0 R0²，B_γ^irr → 0 になります（γ は角度なので，動く量が ε2 γ に比例するため）．

### かたさ C の求め方（Eq. (118) とその下の文）

基底状態から等間隔に 7 点 λ_gs + k h（k = −3, …, 3）をとり，E_pot に最小二乗で放物線をフィットして C = 2 ×（2 次の係数）とします．等間隔なので重み付き和で書けます．

```
C = [5(E₃ + E₋₃) + 0(E₂ + E₋₂) − 3(E₁ + E₋₁) − 4E₀] / (42h²)
```

- 殻構造による細かい凹凸をならした曲率で，局所的な 2 階微分（7 点中心差分）とは違います．
- ε2 方向は h = 0.05，γ 方向は h = 5°（微分はラジアンで行います）です．ずらすのはその変数だけで，ほかの変形は基底状態の値に固定します．
- C ≤ 0（その方向で極小になっていない．軽くやわらかい核で起こります）のときは振動数が決まらないので，止めずに警告を出し，ε2 方向は E_0,ε2 = 0，γ 方向は E_0,γ = E_0,ε2 とします（例：³²S を表の形で計算すると γ 方向が C ≤ 0）．`zero_point_mod` を直接呼んだときは今までどおり停止します．

### 手順（main の add_zero_point）

1. 最後の形（入力の形，または `search` で見つけた形）で E_pot を計算します．このときの E_mac と E_s+p を控えておきます．
2. ε2 を ±0.05, ±0.10, ±0.15 ずらした 6 点で E_pot を計算し，7 点から `zero_eps2%cal_energy` で E_0,ε2 を求めます．
3. ε2 ≥ 0.17 なら，γ を ±5°, ±10°, ±15° ずらした 6 点で同じように `zero_gamma%cal_energy` で E_0,γ を求めます．ε2 < 0.17 では，球に近く γ 方向の曲面が求めにくいので，E_0,γ = E_0,ε2 とします．
4. 次のように表示します（²⁰⁸Pb の例）．

```
 --- zero-point energy (Sec. 2.11) [MeV]
   E_pot(eps2 - 0.15 .. eps2 + 0.15)     :   -14.0765  -17.2285  -21.1505  -22.6370  -21.0528  -16.0063   -8.9127
   E_0,eps2              :       0.6870
   E_0,gamma             :       0.6870   (= E_0,eps2, eps2 < 0.17)
   E_zp                  :       1.3739
 --- result with the zero-point energy [MeV]
   E_mac (ground state)  :      -8.7002
   E_s+p (ground state)  :     -13.9368
   E_zp                  :       1.3739
   mass excess M_th      :     -21.2631
```

形を 6 点（γ も計算するときは 12 点）余分に計算するので，時間は 7〜13 倍になります．`zpe = .false.` で零点エネルギーを足さずに計算します．

### プログラムの構成

- `zero_point` は抽象型で，`cal_energy(nucleus, eps2, h, E_gs_lambda(7), E_0_lambda)` を deferred で持ちます．`zero_eps2` と `zero_gamma` がそれを拡張し，慣性の式 (116), (117) だけが違います．
- 公開している関数は `B_eps2_irr(nucleus, eps2)`，`B_gamma_irr(nucleus, eps2)`，`stiffness_fit(E, h)` です．

### 確認したこと

- `stiffness_fit` は任意の放物線で C = 2a を厳密に返し，x⁴ を加えると最小二乗の値になります（test_zero_point）．
- B_ε2^irr(0) = (2/15) M0 R0²，小さい ε2 で B_γ^irr ≈ (2/15) ε2² M0 R0² になります．
- ²⁰⁸Pb の E_zp = 1.37 MeV で，FRDM2012 の表の E_mic^FL − E_s+p = 1.34 MeV（球なので E_zp にあたる）と 0.03 MeV で合います．M_th は −21.26 MeV（表の M_th^FL は −21.31）です．
- ²⁴⁰Pu の M_th は 50.90 MeV で，表の 50.69 MeV より 0.2 MeV 高く出ます．ずらすときに ε4, ε6 を固定するか最小化し直すか，M0 の定義などが原因の候補ですが，まだ分かっていません．

```bash
gfortran -O2 constant_mod.f90 nucleus_mod.f90 micro_constant_mod.f90 zero_point_mod.f90 test_zero_point.f90 -o test_zero_point.exe
./test_zero_point.exe
```

## 残っていること

- 零点エネルギー：基底状態が三軸非対称な核で慣性を補間する処理（FRDM2012 Sec. 2.11）は入れていません．
- 零点エネルギー：²⁴⁰Pu で表より 0.2 MeV 高い原因は分かっていません．
- 基底状態の探索は軸対称 (ε2, ε3, ε4, ε6) だけです．γ を含む探索，最小値の補間，重い核で分裂アイソマーの極小を基底状態としない規則（FRDM2012 3 章の手順 5）は入れていません．
- D2h 対称性によるハミルトニアンのブロック対角化（ho_basis_mod の TODO）．速くなりますが結果は変わりません．

## 運用

1. [プロジェクトのフォークを作成]
2. [変更を加える]
3. [プルリクエストを送信]



## 更新履歴

* 2026/4/21: [README.md を作成]
* 2026/5/14: [依存関係を明記]
* 2026/10/12: [main を FRLDM+Shell+LN に対応，ε4, ε6 の形，B_1, B_3 の表面積分，Lipkin–Nogami，スピン軌道項の部分積分]
* 2026/10/15: [Fortran のオブジェクト指向プログラミングの解説書 (docs/fortran_oop_guide.pdf) を追加]
* 2026/10/18: [八重極変形 ε3 と重心の ε1，振動子基底の打ち切りの修正，零点エネルギー (zero_point_mod)]
* 2026/10/24: [README に ε3 と零点エネルギーの説明を追加]
