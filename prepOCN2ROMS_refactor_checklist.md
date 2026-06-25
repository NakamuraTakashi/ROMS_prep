# prepOCN2ROMS リファクタリング 進捗チェックリスト

計画: [`prepOCN2ROMS_refactor_plan.md`](prepOCN2ROMS_refactor_plan.md)
ブランチ: `refactor/prepOCN2ROMS`

回帰テスト: `Projects/TokyoBay/` で `./regress_prepOCN2ROMS.sh [his|bry|ini|all]`
（working tree と master をビルドして同一入力で実行し byte-identical を確認。入力は `regress_TokyoBay2_movejpn.in`、要 MOVE-JPN ローカルデータ）

各段階は「コンパイル成功」＋「master と出力一致（回帰）」を満たしたらチェック。

## 段階 1: `OUT_FILE` 共通化（挙動不変）✅
- [x] `HIS_FILE` / `BRY_FILE` → `OUT_FILE` に統一（宣言・全参照）
- [x] 生成呼び出しのみ CPP 分岐（`#if BRY_MODE: createNetCDFbry #else: createNetCDFini`）※元から CPP 分岐済み
- [x] 全モード組合せでコンパイル成功（BRY/HIS/INI × HYCOM/JCOPE/ROMS/MOVEJPN/FORA, WET_DRY）
- [x] master と実行コードが一致（=挙動不変）を確認（`git diff` が純粋なリネームのみ）
- [x] コミット

## 段階 2: `T_region` 派生型を導入し HIS を region(1) で書き直す ✅
- [x] `T_region` 型定義（(A)領域サイズ成分／(B)ドナー box 成分）
- [x] (A) 作業/重み/スライス配列の宣言を `pointer` 化（`region` は `target`）
- [x] HIS/INI を `region(1)`（n=1, 全領域）で書き直し（A配列を region(1) に確保しポインタ別名）
      ※(B)ドナー配列は当面 allocatable のまま（段階3で region 化）
- [x] コンパイル成功（HIS/INI/BRY × MOVEJPN/HYCOM/ROMS）
- [x] HIS 出力が master と **byte-identical**（`cdo diffn` 差分0／`cmp` 一致）。INI は HIS と同一コード経路＋コンパイル確認。
- [x] （追加確認）BRY 出力も master と **byte-identical**（pointer 化の影響なし）
- [x] コミット

## 段階 3: BRY を region 化（time 外側へ再構成・単一ファイル追記）
- [x] 3a: (B) ドナー配列を `pointer` 化（HIS は直接確保のまま＝byte-identical 確認済み）
- [x] 3b-1: BRY の A 配列を region(ibry) 化（ibry ループ維持）。`allocate(region(1:4))` 追加。BRY byte-identical 確認済み
- [x] 3b-2: BRY のループを `do itime`(外)/`do ibry`(内) へ swap、(B)ドナーを region(ibry) 化
- [x] `do itime`（外）/ `do ibry`（内）へ再構成、単一 `OUT_FILE` へ追記（元から単一ファイル）
- [x] スライス抽出＋書込は `#if BRY_MODE` 内に存続（HIS と統合は段階4）
- [x] コンパイル成功（BRY × JCOPE/ROMS/ROMS+WET_DRY/FORA。HYCOM は master と同じ既知の
      別件エラー＝HYCOM_TIME 生成 include 未指定で本リファクタとは無関係）
- [x] BRY 出力が master と **byte-identical**（TokyoBay2/MOVEJPN）。HIS/INI も byte-identical を再確認
- [x] コミット

## 段階 4: read+interp+(donor read) の物理マージ ✅
- [x] BRY/HIS の重複ブロックを 1 本化（region パラメータ化、書込のみ `#if` 分岐）
- [x] 出力セクションを単一セクションへ統合。HIS/INI 専用の HIS セクション（旧 `#elif HIS||INI`）を削除し、
      BRY セクションを `n_region`(BRY=4/HIS=1) でパラメータ化した単一エンジンに。差分は CPP で局所化:
      ①ファイル生成(createNetCDFbry/bry2 vs createNetCDFini)＋接尾辞/接頭辞、②時刻変数(bry_time 一括 vs
      ocean_time を itime 毎)、③region 数・境界範囲(case vs 全域)、④書込(スライス zeta_<dir> 等 vs 全 field)、
      ⑤INI の Nt=1、⑥NAOTIDE 潮汐加算(zeta_tide は HIS/INI のみ確保)
- [x] コンパイル成功（BRY/HIS/INI × JCOPE/ROMS/ROMS+WET_DRY/FORA/MOVEJPN、HIS/INI×NAOTIDEJ。
      HYCOM は master と同じ既知の別件 include エラーのみ）
- [x] BRY/HIS/INI とも master と **byte-identical**（TokyoBay2/MOVEJPN、`cmp`）
- [x] コミット

## 段階 5: 重み再計算トリガーの修正（§3.1）
- [ ] ファイル切替検出と重み再計算を分離
- [ ] 非 WET_DRY: 初回（複数ファイルは切替時）に全 region 計算・保持
- [ ] ROMS+WET_DRY: 毎ステップ再計算（毎ステップ wetdry マスクを `idt(itime)` で読む）
- [ ] コンパイル成功（特に ROMS_MODEL + WET_DRY）
- [ ] WET_DRY ケースの妥当性確認（時間追従するか）／非 WET_DRY は master と一致
- [ ] コミット

## 段階 6: 仕上げ
- [ ] 全 CPP 組合せの網羅コンパイル
- [ ] 代表データ（TokyoBay2/3, MOVEJPN/HYCOM, 必要なら FORA）で実走回帰
- [ ] CLAUDE.md / 計画 md の追記・更新
- [ ] master へのマージ可否を確認（ユーザー判断）

---
進捗メモ:
- 段階1 完了: `BRY_FILE`/`HIS_FILE` を `OUT_FILE` に統一（純粋リネーム、挙動不変）。全モードコンパイル確認済み。
- 段階2 完了: `T_region` 型導入、(A)配列を pointer 化、HIS/INI を region(1) 経由に。TokyoBay2/MOVEJPN(2日)で master と HIS・BRY とも byte-identical を確認。(B)ドナー配列は段階3で region 化予定。
- 段階3a 完了: (B)ドナー配列を pointer 化（HIS は直接確保で byte-identical）。
- 段階3b-1 完了: BRY の A 配列を region(ibry) 化（`allocate(region(1:4))` 追加、ibry ループ維持）。BRY byte-identical。
- 段階3b-2 完了: BRY を `do itime`(外)/`do ibry`(内) へ swap。setup（bounds+A確保）を前出し、GRIDLOAD
  (HYCOM 共有座標 load / ファイル open) を ibry 外へ hoist、(B)ドナーを region(ibry) 化＋bare ポインタへ alias、
  box スカラ(Irdg_min..,Nxr_dg..)を region(ibry) に保存し非ファイル切替 step で復元。ファイル切替検出は
  原コードどおり `if(itime<Nt){ if(iNC==iNCt(itime+1)) cycle }` を使用（`.or.` は短絡しないので
  `iNCt(Nt+1)` 範囲外参照を避ける）。TokyoBay2/MOVEJPN で BRY/HIS/INI とも master と byte-identical。
- 段階4 完了: 出力部の二重実装を解消。HIS 専用セクションを削除し、旧 BRY セクションを `n_region` で
  パラメータ化した単一エンジンへ統合（read+interp+donor-read は完全共有、差分は createNetCDF/時刻変数/
  region 範囲/書込/Nt=1/NAOTIDE のみ CPP 局所化）。BRY-active 経路は不変なので BRY は構成上不変、
  HIS/INI の `#else` 分岐を新規追加。TokyoBay2/MOVEJPN で BRY/HIS/INI とも master と byte-identical、
  全 mode×donor 行列でコンパイル確認。次は段階5（§3.1 重み再計算トリガー: 非WET_DRYは初回のみ、
  ROMS+WET_DRYは毎ステップ）。

### 3b-2 手術手順（次回実施・要点メモ）
現 BRY: `do ibry { setup; iNC=0; do itime { if(iNCm<iNC){GRIDLOAD(1741-1810)+BOX(1811-1958)}; RW(1960-2405); cycle/donor dealloc(2406-2434) } enddo }`
目標: setup を前出し、loop を swap、GRIDLOAD は ibry 外へ hoist:
```
do ibry { if(SNWE)cycle; bounds→region(ibry); allocate region(ibry)%A } enddo   ! setup 前出し（retarget は外す）
iNC=0
do itime
  iNCm=iNC; iNC=iNCt(itime)
  if(iNCm<iNC): GRIDLOAD（共有 latr_dg 等。ibry 外で1回）   ! ←hoist 必須（ibry内だと共有座標が二重 allocate でエラー）
  do ibry { if(SNWE)cycle
     A/donor ポインタを region(ibry) に retarget; box スカラ Irdg_min 等を region(ibry) から復元
     if(iNCm<iNC): BOX = seek→region(ibry)%Irdg_min 等; allocate region(ibry)%donor + alias; weights→region(ibry)
     RW（素の名前のまま）
  } enddo
  if(itime==Nt .or. iNC/=iNCt(itime+1)): 共有座標 deallocate（1回）＋ region(ibry)%donor を全 ibry 解放
enddo
```
要点: ①GRIDLOAD/BOX 分割点 = 1810(#endif)/1811(Seek)。②box スカラ(Irdg_min..,Nxr_dg..)は region(ibry) に保存し毎 ibry 復元（非ファイル切替step で stale を防ぐ）。③ドナー配列は region(ibry)%＋alias。④`cycle` ベースの解放を `if(itime==Nt .or. iNC/=iNCt(itime+1))` に。⑤行長 ≤132 厳守（production は -ffree-line-length-none 無し）。検証: TokyoBay2/MOVEJPN で BRY byte-identical。
