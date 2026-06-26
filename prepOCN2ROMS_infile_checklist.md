# prepOCN2ROMS × mod_infile 統一 進捗チェックリスト

計画: [`prepOCN2ROMS_infile_plan.md`](prepOCN2ROMS_infile_plan.md)
ブランチ: `refactor/prepOCN2ROMS`

回帰: `Projects/TokyoBay/` で `./regress_prepOCN2ROMS.sh [his|bry|ini|all]`
（working tree と master をビルドし同一入力で実行。今回は時刻のジュリアン日正規化に伴い
**`cdo diffn` 差0** を合格基準とする＝`cmp` 完全一致は問わない。harness は cdo diffn フォールバック済み）

各段階は「コンパイル成功」＋（MRICOM は）「master と cdo diffn 差0」を満たしたらチェック。

## S0: 計画・チェックリストを md 化 ✅
- [x] `prepOCN2ROMS_infile_plan.md` / `prepOCN2ROMS_infile_checklist.md` 作成
- [x] コミット

## S1: mod_infile 完成 ✅
- [x] `T_INFILE` をコア確定（`NAME(:), time_all(:), Nt, ItS, ItE`、`time_all` はジュリアン日）※既存のまま
- [x] `infile_check_time(Nfile, t_Start, t_End, Nt, time, iNCt, idt)` へ拡張
      （`iNCt/idt` は `optional, allocatable, intent(out)`。両方 present のときのみ確保・構築）
- [x] 統合 `time/iNCt/idt` 構築ロジック（ROMS 流）を実装、**末尾ループのバグ修正**
      （旧 `j = i + ...` のループ変数上書きを廃止、`ItS==-1` ファイルは cycle）
- [x] `mod_infile` 単体コンパイル成功
- [x] コミット

## S2: prepOCN2ROMS MRICOM パスを INFILE 化 ✅
- [x] `use mod_infile` 追加
- [x] MRICOM の時刻フェーズを INFILE 化:
      - FORP: `INFILE`構築＋`infile_check_time(...,iNCt,idt)` に集約、`time_all` を hour→jdate 正規化
      - MOVEJPN/FORA: 日付スキャンで `INFILE%NAME(5)/time_all(1)=jdate/Nt=1`、`iNCt/idt` 直接構築
- [x] `bry_time` を jdate から秒換算（FORP=`(time-ref)*86400`、MOVEJPN=スキャン後 in-place 換算）
- [x] MRICOM の `T_NC`/`NC(:)` 定義削除、grid 読込・出力の `NC%MRICOM_FILE/%ItS` を `INFILE%NAME/%ItS` に置換
- [x] harness の COMMON に `mod_infile.F90` を追加
- [x] コンパイル成功（MOVEJPN/FORA/FORP × BRY/HIS/INI）
- [x] MOVE-JPN で master と **byte-identical**（cmp 完全一致。cdo diffn 差0 基準を上回る）
- [x] コミット

## S3: ROMS パスを INFILE 化 ✅
- [x] ROMS の時刻フェーズを INFILE＋`infile_check_time` に置換（`ocean_time` 秒→jdate 正規化）、
      `bry_time = (time-ref)*86400` に。`ROMS_HISFILE(:)` は維持し `INFILE%NAME(1)` にコピー
- [x] ROMS の `T_NC`/`NC(:)` 定義削除、重複 `allocate(NC(NCnum))` も撤去（INFILE は時刻フェーズで確保）
- [x] コンパイル成功（ROMS × BRY/HIS/INI、+WET_DRY）。MOVEJPN も無回帰
- [x] コミット

## S4: HYCOM パスを INFILE 化＋HYgrid 外出し ✅
- [x] `T_HYCOM_GRID{lat,lon,depth,Nxr_dg,Nyr_dg,Nzr_dg}` ＋ `HYgrid(:)` を新設
- [x] HYCOM 格子読込先・GRIDLOAD 参照を `NC%lat/lon/depth/N*_dg` → `HYgrid%…`、`NC%ItS` → `INFILE%ItS` に置換
- [x] HYCOM の時刻フェーズを INFILE＋`infile_check_time` に置換、`time_all` を hour(since 2000)→jdate 正規化、
      `bry_time=(time-ref)*86400`。`HYCOM_FILE(:)` は維持し `INFILE%NAME(1)` にコピー
- [x] HYCOM の `T_NC`/`NC(:)` 定義削除、重複 `allocate(NC(NCnum))` 撤去
- [x] コンパイル成功（HYCOM_LOCAL(+FAST_READ) × BRY/HIS/INI）。MOVEJPN/ROMS も無回帰
- [x] コミット
  - 注意: 非 HYCOM_LOCAL（GOFS）の time キャッシュ `.dat` は jdate で書く形に変えたため、
    旧 hour ベースの `.dat` とは非互換（`SKIP_CHECK_TIME` で旧キャッシュを読む場合は再生成が必要）。
    HYCOM_LOCAL は `.dat` 不使用で影響なし。

## S5: JCOPE パスを INFILE 化
- [ ] JCOPE の時刻フェーズ（日付生成名）を INFILE＋infile_check_time に統一（`Nt=1`）
- [ ] コンパイル成功（JCOPE × BRY/HIS/INI）
- [ ] コミット

## S6: 仕上げ
- [ ] 網羅コンパイル（BRY/HIS/INI × JCOPE/ROMS/ROMS+WET_DRY/FORA/MOVEJPN/HYCOM_LOCAL(+FAST_READ)、
      HIS/INI×NAOTIDE/NAOTIDEJ）
- [ ] 最終 cdo diffn 回帰（MOVE-JPN）
- [ ] CLAUDE.md 追記、master マージ可否をユーザー確認
- [ ] コミット

---
進捗メモ:
- （随時更新）
