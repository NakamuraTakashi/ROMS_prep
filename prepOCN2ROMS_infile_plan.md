# prepOCN2ROMS: ドナー入力ファイル情報の構造体を mod_infile に統一（実装計画）

対象ブランチ: `refactor/prepOCN2ROMS`
対象ファイル: `src/mod_infile.F90`, `src/prepOCN2ROMS.F90`
進捗: [`prepOCN2ROMS_infile_checklist.md`](prepOCN2ROMS_infile_checklist.md)

## 1. 目的・背景

`src/prepOCN2ROMS.F90` は HYCOM / MRICOM / ROMS の各ドナーで「入力ファイル情報＋時刻リスト」を
**別々の `T_NC` 派生型**（CPP 分岐で3定義）と、ドナーごとに重複した時刻リスト構築ロジックで
扱っている。共通コア（`time_all, Nt, ItS, ItE`）は同じなのに型もコードも分かれており見通しが悪い。
`frcATM2ROMS` 系でも同種の構造体を使うため、**1つの共通型＋共通ルーチンを別モジュール
`mod_infile.F90` に切り出して再利用**する。`mod_infile.F90`（`T_INFILE` ＋ `infile_check_time`）は
作りかけであり、本リファクタで完成させ `prepOCN2ROMS` を載せ替える。

## 2. 確定済みの設計判断（ユーザー合意）

- **A**: `T_INFILE` は **コアのみ**（`NAME(:), time_all(:), Nt, ItS, ItE`）。格子情報
  (`lat/lon/depth, Nxr_dg…`) は型に含めず、ドナー側で別管理。
- **B**: `time_all` は **ジュリアン日に正規化**して格納。`infile_check_time` を拡張し、統合時刻列に
  加えて **`iNCt`(全体step→ファイル番号) / `idt`(ファイル内index)** も返す。未完成バグも修正。
- **C**: 対象は `prepOCN2ROMS.F90` ＋ `mod_infile.F90` のみ。元 `bry/ini/hisOCN2ROMS.F90` などは
  据え置き（将来削除予定）。`frcATM2ROMSv2.F90` は未完成のため今回は対象外（触らない・壊れても可）。
- **D**: 出力時刻のジュリアン日往復変換による **最下位ビット差は許容**。検証は `cmp` 完全一致ではなく
  **`cdo diffn` 差0** を基準にする（既存 regress harness は cdo diffn フォールバック済み）。
- **JCOPE/MRICOM**（1時刻=1ファイル / 1日=1ファイル）も解析後は `INFILE(j)%NAME`＋`time_all(1)=その時刻`
  ＋`Nt=1` として格納し、ROMS/HYCOM と同じ `infile_check_time` 経路に統一。

## 3. 現状の要点（調査済み）

- `T_NC` 定義3つ: `prepOCN2ROMS.F90` ROMS版 / HYCOM版(+lat/lon/depth/Nxr_dg…) / MRICOM版(+`MRICOM_FILE(5)`)。
  共通コアは `time_all,Nt,ItS,ItE`。
- 時刻リスト構築は各ドナーブロック（`!= Read grid coordinate =` 以下の
  `#if ROMS / #elif HYCOM / #elif JCOPE / #elif MRICOM`）内に散在。代表ロジック = ROMS の
  重複トリミング→ItE/ItS探索→Nt集計(iNCs/iNCe)→`bry_time/idt/iNCt` 構築。HYCOM/MRICOM/JCOPE も同型。
- HYCOM 格子は別ループで `NC(iNC)%lat/lon/depth/Nxr_dg…` に格納し、出力フェーズの GRIDLOAD で
  `latr_dg/z_r_dg` 構築に使用。
- `mod_infile.F90`: `T_INFILE{NAME,time_all,Nt,ItS,ItE}` ＋ `infile_check_time(Nfile,t_Start,t_End,Nt,time)`。
  **末尾ループにバグ**（`do j=js,je` 内で `j` を上書き＝未定義動作）。

## 4. 実装方針

### 4.1 `mod_infile.F90` を完成（共通基盤）
- `T_INFILE` はコアのまま（`NAME(:),time_all(:),Nt,ItS,ItE`）。`time_all` はジュリアン日。
- `infile_check_time(Nfile, t_Start, t_End, Nt, time, iNCt, idt)`:
  - `time(:)`(jdate), `optional` の `iNCt(:)/idt(:)` を `allocatable, intent(out)` で返す。
  - トリミング→ItE/ItS→Nt集計→`time/iNCt/idt` 構築（ROMS ロジックを移植）。
  - **末尾ループのバグ修正**（ループ変数を上書きしない）。

### 4.2 `prepOCN2ROMS.F90` を `mod_infile` に載せ替え
- 先頭で `use mod_infile`。3つの `T_NC` 定義と `NC(:)` を撤去（JCOPE は元々 T_NC なし）。
- 各ドナーの時刻フェーズを「`INFILE` 構築 → `infile_check_time(...,iNCt,idt)`」に統一:
  - ファイル名 → `INFILE(j)%NAME(:)`（ROMS=`ROMS_HISFILE`, HYCOM=`HYCOM_FILE`, MRICOM=5本,
    JCOPE=日付生成名）。namelist 変数は従来通り読み、`NAME` にコピー。
  - `time_all` をジュリアン日格納（例 ROMS: `d_jdate_Ref + ocean_time/86400`）。
  - 出力 `bry_time`(秒) は統合 `time`(jdate) から換算: `bry_time = (time - d_jdate_Ref)*86400`。
  - JCOPE/MRICOM は `Nt=1, time_all(1)=その時刻`。
- **HYCOM 格子の外出し**: `prepOCN2ROMS.F90` 内に軽量型 `T_HYCOM_GRID{lat,lon,depth,Nxr_dg,Nyr_dg,Nzr_dg}`
  ＋ `HYgrid(:)` を新設（`mod_infile` には入れない）。格子読込先と GRIDLOAD 参照を
  `NC(iNC)%lat…` → `HYgrid(iNC)%lat…` に置換。

## 5. 段階分け（各段階でコンパイル＋該当回帰）

- **S0**: 本計画と進捗チェックリストを md に書き出す（このファイル群）。
- **S1**: `mod_infile` 完成（`infile_check_time` 拡張[新引数 optional]＋バグ修正）。単体コンパイル。
- **S2**: `prepOCN2ROMS` に `use mod_infile`、**MRICOM パス**を INFILE 化（MRICOM `T_NC` 削除）。
  → MOVE-JPN で `cdo diffn` 回帰（BRY/HIS/INI）。
- **S3**: ROMS パスを INFILE 化（ROMS `T_NC` 削除）→ コンパイル。
- **S4**: HYCOM パスを INFILE 化＋`HYgrid` 外出し（HYCOM `T_NC` 削除）→ コンパイル。
- **S5**: JCOPE パスを INFILE 化 → コンパイル。
- **S6**: 網羅コンパイル＋最終 `cdo diffn` 回帰。

## 6. 検証

- 回帰: `Projects/TokyoBay/regress_prepOCN2ROMS.sh all`（MOVE-JPN/MRICOM, 2日=ファイル切替あり）で
  BRY/HIS/INI とも **`cdo diffn` 差0** を確認（`cmp` 完全一致は時刻 ULP 差で崩れ得るため不問。
  time 座標が cdo diffn で差検出される場合は harness を微調整）。
- 網羅コンパイル: BRY/HIS/INI × JCOPE/ROMS/ROMS+WET_DRY/FORA/MOVEJPN/HYCOM_LOCAL(+FAST_READ)、
  HIS/INI×NAOTIDE/NAOTIDEJ。
- HYCOM/FORA/ROMS はローカルデータが無いためコンパイルのみ（MRICOM のみ数値回帰）。

## 7. リスク・注意

- `time_all` ジュリアン日化で出力時刻が ULP オーダーで変わり得る（D で許容、cdo diffn 差0 基準）。
- `T_NC` 撤去は3ドナーに波及するが CPP 分岐なので S2〜S5 で1ドナーずつ移行・各段階コンパイル可能。
- HYCOM 格子は「ファイル毎に異なり得る」ため `HYgrid(iNCs:iNCe)` で全ファイル分保持。
- `infile_check_time` の新引数は optional にして frcATMv2（未完成・対象外）の既存呼出を壊さない。
- 元 `bry/ini/hisOCN2ROMS.F90` 等は touch しない。
