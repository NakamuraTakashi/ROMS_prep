# prepOCN2ROMS リファクタリング実装計画 — BRY / HIS 出力部の統合

対象ブランチ: `refactor/prepOCN2ROMS`
対象ファイル: `src/prepOCN2ROMS.F90`（必要に応じ `src/mod_roms_netcdf.F90`）

> **状態（2026-06-25）: 段階 1〜6 実装完了、master マージ可否のユーザー判断待ち。**
> 進捗の詳細は [`prepOCN2ROMS_refactor_checklist.md`](prepOCN2ROMS_refactor_checklist.md)。
> 出力部の二重実装を解消し単一 region パラメータ化エンジンへ統合（`src/prepOCN2ROMS.F90`
> は master 比 約524行純減）。非 WET_DRY は TokyoBay2/MOVE-JPN で master と byte-identical。
> ROMS+WET_DRY のみ §3.1 の重みトリガー修正で挙動が変わる（実走データ未入手で論理担保）。

## 1. 目的

- `BRY_MODE`（現 1511 行〜）と `HIS_MODE`/`INI_MODE`（現 2435 行〜）の出力部に**重複**している
  「ドナー読み込み＋内挿＋書き込み」を、**領域 (region) パラメータ化**で 1 本化する。
- BRY を **time 外側・region 内側**に再構成し、**単一 bry ファイル**へ `bry_time` を
  UNLIMITED 追記する方式にする（4 分割は不要）。
- I/O は **境界別の小領域リクエストを維持**（閉境界はスキップ）。データ量・転送途中切断の観点で
  全領域取得より有利なため。OPeNDAP のリトライ／再オープン機構はそのまま有効。
- 出力ファイル変数 `HIS_FILE` / `BRY_FILE` を共通名 `OUT_FILE` に統一する。

## 2. 現状の構造（精読結果）

共有エンジン（〜1510 行：ROMS 格子読込・ドナー格子読込・時刻リスト構築）は既に共通。
**出力部のみ**が 2 つに分かれている。

### BRY_MODE（1511–2433）
```
createNetCDFbry/bry2(BRY_FILE, ..., SNWE)        ! 1 ファイルに全有効境界, bry_time は既に UNLIMITED
do ibry=1,4                                       ! ← 境界が外側
  if (SNWE(ibry)==0) cycle
  境界範囲を設定 (LBri/UBri/LBrj/UBrj, LBui.., LBvi..)
  作業配列・スライス配列・重み配列を allocate
  do itime=1,Nt                                   ! ← 時間が内側
    if (ドナーファイル切替 iNCm<iNC):
        donor 格子 load → donor IJ box 探索 (Irdg_min.. 等) → 重み計算 (weight2D/3D_grid3_2)
    read donor box の itime → interp(境界範囲) → 境界スライス抽出 → 書込(index itime)
  enddo
enddo
```

### HIS_MODE / INI_MODE（2435–末尾）
```
createNetCDFini(HIS_FILE, ...)                    ! ocean_time UNLIMITED
全領域範囲を設定 (LBri=0,UBri=L, ...)
作業配列・重み配列を allocate
#if INI_MODE: Nt = 1
do itime=1,Nt                                     ! ← 時間が外側
  if (ドナーファイル切替): donor load → IJ box 探索 → 重み計算
  read donor → interp(全領域) → 全フィールド書込(index itime)
enddo
```

### 重要な事実
- **HIS の read+interp ブロックは、BRY の対応ブロックと bounds の値以外バイト一致**。
  差異は実質 2 点のみ:
  1. **領域範囲の値**（境界部分集合 vs 全領域）。
  2. **書き込み**（BRY: スライス `zeta_<dir>` を `start=(/1,itime/)` ／ HIS: 全 field `zeta` を
     `start=(/1,1,itime/)`）。
- 重みは **(領域, ドナーファイル)** の関数。`if(iNCm<iNC)` でファイル切替時に再計算。
- `createNetCDFbry(..., snwe_flag(4))` は **境界フラグで変数を出し分け済み**
  （`snwe=(1,0,0,0)` なら south の変数のみ生成）。`bry_time` は既に `NF90_UNLIMITED`。

## 3. 設計：領域 (region) 抽象

`region(ibry)`, `ibry = 1..n` を導入する。

- **HIS / INI**: `n = 1`、`region(1) = 全領域`。
- **BRY**: `n = 4`、`region(ibry) = 各境界`（`SNWE(ibry)==0` は `cycle`）。

各 `region` が保持する情報（Fortran 派生型 `T_region` の allocatable 成分にする）。
**現状 `do ibry` ループ内で alloc/dealloc を繰り返している配列を `T_region` 成分に移し、
`region(1:n)` を確保して使い回す**（反復 alloc/dealloc を廃止）。確保タイミングで 2 種に分かれる:

| 成分 | サイズ依存 | 確保タイミング |
|------|------|------|
| ターゲット範囲 `LBri/UBri/…`, `Nxr/Nyr/…`（スカラ） | region のみ | **(A)** region 設定時に 1 回 |
| 内挿の重み `ID_cnt2Dr/w_cnt2Dr/ID_cnt3Dr/w_cnt3Dr`（+ u, v） | region のみ | **(A)** region 設定時に 1 回（**配列は再確保せず**、値のみファイル切替時に再計算） |
| 作業配列 `zeta, ubar, vbar, u, v, t`（+ `ull, uu, vu, vll, uv, vv`） | region のみ | **(A)** region 設定時に 1 回 |
| (BRY) スライス `zeta_bry, ubar_bry, vbar_bry, u_bry, v_bry, t_bry` | region のみ | **(A)** region 設定時に 1 回 |
| ドナー IJ box 添字 `Irdg_min/max, …`, `Iudg.., Ivdg..` | region × ファイル | **(B)** ファイル切替時 |
| ドナーデータ `zeta_dg, u_dg, v_dg, t_dg`（+ `ull_dg` 等） | ドナー box（region × ファイル） | **(B)** ファイル切替時に realloc |

**region に入れない（共有のまま）**: ドナー側の**格子全体**の座標
`latr_dg / lonr_dg / z_r_dg / rmask_dg …` は、その時刻のドナーファイルの格子そのもので
全境界で共通。各 region の `seek_IJrange` はこの共有座標から自分の box を探すだけなので、
`T_region` には入れず**ファイル load 時に 1 回確保**（現状どおり）。

> (A) はサイズが領域だけで決まるので region 設定時に確保したまま保持。
> (B) はドナー box がファイルで変わり得るのでファイル切替時に確保し直す。
> いずれも周囲シェル分で小さく、4 領域同時保持でもメモリは軽微。

### 3.1 重みの再計算トリガー（重要・要修正）

重み（`ID_cnt/w_cnt`）の **値**を再計算する条件は次のとおり:

```
need_weights = ( itime==1 .or. ドナーファイル切替 )  .or.  WET_DRY
```

- **非 WET_DRY（HYCOM/JCOPE/MRICOM/FORA、および WET_DRY 無効の ROMS）**: 海域・陸域の配置は
  固定なので **初回 1 回**（複数ファイルで格子が変わり得る場合のみファイル切替時に再計算）。
- **ROMS_MODEL + WET_DRY**: 親の wet/dry 配置が**毎ステップ変化**するため、**毎ステップ再計算**。
  あわせて WET_DRY 時はドナーの wetdry マスクを**毎ステップ `idt(itime)` で読み**、そのマスクで
  重みを計算する。
- 配列（A）自体は再確保しない。再計算するのは**値だけ**。

> **現状の問題（要修正）**: 現コードは重みを `if(iNCm < iNC)`（ファイル切替時のみ）で計算しており、
> WET_DRY でも**固定時刻のマスク**で 1 回計算するだけ。一方データ側は毎ステップ wetdry で
> マスクしている（`u_dg *= umask_wet`）ため、**重みが wet/dry の時間変化に追従しない不整合**がある。
> 本リファクタで上記トリガーに改める（ファイル切替検出と重み再計算を分離し、WET_DRY 時は毎ステップ
> 再計算）。ROMS+WET_DRY のマスク供給経路は分岐が深いので実装時に丁寧に追跡する。

## 4. 目標構造

```fortran
! ---- 出力ファイル生成（モード固有） ----
  OUT_FILE = trim( OUT_prefix )//OUT_suffix
#if defined BRY_MODE
  call createNetCDFbry( OUT_FILE, ..., romsvar, SNWE )   ! 単一ファイル・全有効境界
#else
  call createNetCDFini( OUT_FILE, ... )
#endif

! ---- region 設定（ibry=1..n の範囲確定・配列確保） ----
#if defined BRY_MODE
  n_region = 4
#else
  n_region = 1
#endif
  do ibry = 1, n_region
#if defined BRY_MODE
    if (SNWE(ibry)==0) cycle
    region(ibry) の境界範囲を設定（現 1572–1636 の case）
#else
    region(1) の全領域範囲を設定（現 HIS の範囲）
#endif
    region(ibry) の作業・スライス・重み配列を allocate
  enddo

! ---- 時間ループ（共有） ----
#if defined INI_MODE
  Nt = 1
#endif
  do itime = 1, Nt
    iNC 更新
    if (ドナーファイル切替 iNCm<iNC) then
      donor 格子 load                                     ! ← 共有
      do ibry = 1, n_region
        region(ibry) の donor IJ box 探索 + 重み計算       ! ← 共有（region でパラメータ化）
      enddo
    endif
    ocean_time(1) = bry_time(itime)
    時刻座標を index itime に追記（bry_time / ocean_time）
    do ibry = 1, n_region
      read donor box(ibry) の itime + interp(region ibry)  ! ← 共有
#if defined BRY_MODE
      境界スライス抽出 → zeta_<dir> 等を index itime に書込
#else
      全フィールド zeta/u/v/temp/salt を index itime に書込
# if defined NAOTIDE || defined NAOTIDEJ
      （INI: zeta に潮汐加算）
# endif
#endif
    enddo
  enddo
```

## 5. 共有 / モード固有の切り分け

- **完全共有**：時間ループ骨格、ファイル切替検出、donor 格子 load、donor IJ box 探索、
  重み計算、donor read、interp。← ここの重複が消えるのが本リファクタの主目的。
- **region でパラメータ化**：境界範囲、donor box、重み、配列サイズ。
- **CPP モード固有（最小限に残す）**：
  - 出力生成（`createNetCDFbry` vs `createNetCDFini`）
  - `n_region`（4 vs 1）と境界範囲設定（case vs 全領域）
  - 書き込み（スライス `zeta_<dir>` vs 全 field `zeta`）
  - 時刻変数名（`bry_time` vs `ocean_time`）
  - `INI_MODE` の `Nt=1` と NAOTIDE 潮汐加算

## 6. 段階的実装手順（各段階でコンパイル＋TokyoBay 実走の回帰確認）

1. **`OUT_FILE` 共通化（挙動不変）**
   `HIS_FILE`/`BRY_FILE` → `OUT_FILE`、`HIS_prefix`/`BRY_prefix` の扱いを整理。
   生成呼び出しのみ `#if BRY_MODE: createNetCDFbry #else: createNetCDFini`。
   → 出力が master と一致することを確認（最小・低リスクの先行コミット）。

2. **`T_region` 派生型を導入し、HIS を `region(1)`（n=1, 全領域）で書き直す**
   HIS は元々 1 領域なので、派生型の土台を**挙動不変**で固められる。
   → HIS/INI 出力が master と一致することを確認。

3. **BRY を region 化（time 外側へ再構成）**
   現 `do ibry`（外）→ region 設定（ループ前）に移し、`do itime`（外）/ `do ibry`（内）へ。
   スライス抽出＋書込を `#if BRY_MODE` 内に。単一 `OUT_FILE` へ追記。
   → BRY 出力が master と一致することを確認（下記 7 の方法）。

4. **read+interp+(donor read) の物理マージ**
   手順 2・3 で構造が揃った BRY/HIS の重複ブロックを 1 本化（region パラメータ化、
   書込のみ `#if` 分岐）。`#elif defined HIS_MODE || INI_MODE` と `#if BRY_MODE` の
   出力セクションを単一セクションへ統合。

5. **網羅コンパイル＋回帰**
   `BRY/HIS/INI × HYCOM(+LOCAL)/JCOPE/ROMS/MOVEJPN/FORA`、`WET_DRY`、`NAOTIDEJ` の
   組合せをコンパイル。代表ケースを実走し master と一致確認。

## 7. 検証方法

- **挙動保存の確認**：同一入力に対し、master でビルドした実行ファイルと本ブランチの
  実行ファイルの出力 netCDF を比較。
  - `cdo diffn old.nc new.nc`（全変数の差分統計）または
  - `ncks`/`ncdump` ＋ 数値比較。差分ゼロ（丸め誤差以内）を確認。
- **代表データ**：`COAWST_DATA/TokyoBay/TokyoBay2`（MOVEJPN ローカル、HYCOM）、
  必要に応じ FORA（OPeNDAP）。
- ステージ 1・2 は厳密一致、ステージ 3・4 は BRY の単一ファイル化に伴い**ファイル構成は
  同じ（単一 bry ファイル）**なので、変数・値の一致を確認できる。

## 8. 注意点・リスク

- **重みの再計算粒度（§3.1 参照）**：非 WET_DRY は初回（複数ファイルはファイル切替時）に
  **全 region** を計算して `T_region` に保持し、以後は再計算しない。**ROMS+WET_DRY のみ毎ステップ
  再計算**（毎ステップ wetdry マスクを読む）。現状はこの区別がなくファイル切替時のみ計算＝
  WET_DRY が時間追従しない不整合があるため、本リファクタで修正する。
- **I/O 方針の維持**：donor read は region ごとの**小領域**（閉境界スキップ）を維持する。
  全領域取得（HIS 流の 1 リクエスト化）は採用しない（データ量・切断リスクで不利な場合が多い）。
- **分岐の多さ**：read ブロックは `ROMS_MODEL`(refine)・`WET_DRY`・各ドナー（HYCOM/JCOPE/
  MRICOM）で分岐が多い。**先に HIS(n=1) で土台を固めてから BRY を載せる**ことで、
  既存の動作を壊しにくくする。
- **`INI_MODE`**：`n=1`（全領域）、`Nt=1`、NAOTIDE 潮汐加算は zeta 書込直前に挿入（現状維持）。
- **時刻変数名**：`bry_time`（bry ファイル）/ `ocean_time`（his ファイル）はファイル形式に
  紐づくため、`OUT_FILE` を共通化しても変数名は CPP で出し分ける。
- **3 つのオリジナル `*OCN2ROMS.F90` と既存 `run_*.sh` は触らない**（CLAUDE.md の方針どおり）。
- **共有 `mod_*` の変更は最小に**。今回は基本 `prepOCN2ROMS.F90` 内で完結する見込み。

## 9. 想定される成果

- 出力部の重複（read+interp の二重実装）が解消し、モード固有は「生成・書込・少数の定数」だけになる。
- BRY が単一ファイル＋UNLIMITED 追記となり、HIS と同一の時間ループ構造に揃う。
- I/O 効率（境界別小領域・閉境界スキップ）は維持。
- 以降の機能追加（新ドナー・新モード）が 1 箇所の共通エンジンへの追加で済むようになる。
