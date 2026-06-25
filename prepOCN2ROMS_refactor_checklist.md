# prepOCN2ROMS リファクタリング 進捗チェックリスト

計画: [`prepOCN2ROMS_refactor_plan.md`](prepOCN2ROMS_refactor_plan.md)
ブランチ: `refactor/prepOCN2ROMS`

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
- [ ] 3b-2: BRY のループを `do itime`(外)/`do ibry`(内) へ swap、(B)ドナーを region(ibry) 化
- [ ] `do itime`（外）/ `do ibry`（内）へ再構成、単一 `OUT_FILE` へ UNLIMITED 追記
- [ ] スライス抽出＋書込を `#if BRY_MODE` 内に
- [ ] コンパイル成功（BRY × 各ドナー）
- [ ] BRY 出力が master と一致（`cdo diffn`）
- [ ] コミット

## 段階 4: read+interp+(donor read) の物理マージ
- [ ] BRY/HIS の重複ブロックを 1 本化（region パラメータ化、書込のみ `#if` 分岐）
- [ ] 出力セクションを単一セクションへ統合（`#if BRY_MODE` … `#else` … `#endif`）
- [ ] コンパイル成功（全組合せ）
- [ ] BRY/HIS/INI とも master と一致（`cdo diffn`）
- [ ] コミット

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
