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

## S1: mod_infile 完成
- [ ] `T_INFILE` をコア確定（`NAME(:), time_all(:), Nt, ItS, ItE`、`time_all` はジュリアン日）
- [ ] `infile_check_time(Nfile, t_Start, t_End, Nt, time, iNCt, idt)` へ拡張
      （`iNCt/idt` は `optional, allocatable, intent(out)`）
- [ ] 統合 `time/iNCt/idt` 構築ロジック（ROMS 流）を実装、**末尾ループのバグ修正**
- [ ] `mod_infile` 単体コンパイル成功
- [ ] コミット

## S2: prepOCN2ROMS MRICOM パスを INFILE 化
- [ ] `use mod_infile` 追加
- [ ] MRICOM の時刻フェーズを `INFILE` 構築＋`infile_check_time` に置換（`Nt=1` 形式）
- [ ] `bry_time` を統合 `time`(jdate) から秒換算
- [ ] MRICOM の `T_NC` 定義削除（NC 参照を INFILE/格子別管理に置換）
- [ ] コンパイル成功（MOVEJPN × BRY/HIS/INI）
- [ ] MOVE-JPN で master と **cdo diffn 差0**（BRY/HIS/INI）
- [ ] コミット

## S3: ROMS パスを INFILE 化
- [ ] ROMS の時刻フェーズを INFILE＋infile_check_time に置換、`bry_time` 秒換算
- [ ] ROMS の `T_NC` 定義削除
- [ ] コンパイル成功（ROMS × BRY/HIS/INI、+WET_DRY）
- [ ] コミット

## S4: HYCOM パスを INFILE 化＋HYgrid 外出し
- [ ] `T_HYCOM_GRID{lat,lon,depth,Nxr_dg,Nyr_dg,Nzr_dg}` ＋ `HYgrid(:)` を新設
- [ ] HYCOM 格子読込先・GRIDLOAD 参照を `NC%…` → `HYgrid%…` に置換
- [ ] HYCOM の時刻フェーズを INFILE＋infile_check_time に置換、`bry_time` 秒換算
- [ ] HYCOM の `T_NC` 定義削除
- [ ] コンパイル成功（HYCOM_LOCAL(+FAST_READ) × BRY/HIS/INI）
- [ ] コミット

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
