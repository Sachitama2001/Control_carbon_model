# VISIT native実行による2000年7月15日状態の取得

## 結果

固定ソース `Sachitama2001/VISIT-matrix` commit
`3285bd8e131a932e338b59892751648fd9edcc7b` を、提供されたTKY気象・
窒素沈着・GHGデータでspinupから再実行した。基準日は
2000年7月15日（0始まりdoy=196）である。

成果xlsxは
`artifacts/visit_matrix_workbook/VISIT_TKY_matrix_equations_native_20000715.xlsx`。
従来未設定だった37状態と15履歴をすべてnative実行値で埋めた。
気象値は当日の `f_loct_proc` 後を採用した。未設定で残るのは、
原典で値自体が未定義の `z_ngas2_aa` と `z_curry_undefined` の2個である。

Excel式に同じ開始状態・履歴・当日気象を与え、1日後の37状態を
native Cと照合した。最大規格化誤差は `5.44e-16`、最大絶対誤差は
`1.14e-13`で、浮動小数の丸めの範囲で一致した。

## 1日内の時点

VISITは当日開始に一括更新するのではない。診断版 `ansis.c` は次の3時点を
`TKY_native_snapshot_20000715.tsv` へ出力する。

1. `start`: 当日 `f_loct_proc` の前。Excelの状態 `X` と履歴 `h` に採用。
2. `after_loct`: 当日気象の設定、放射、水収支、生理特性の更新後。
   Excelの気象入力 `z` に採用。
3. `end`: `daily_scheme` 後。Excelの `m_end_*` との照合対象。

水貯留は `f_loct_proc` 内ですでに更新される。そのため `after_loct` の水量を
開始状態として用いると二重更新になる。xlsx内の `09_native実行` に
3時点を並べ、`98A_native照合` に全37状態の比較を収録した。

## ソースと実行条件

- 式と実行順序: 上記の固定2013年ソース。
- 気象: `ext_28-73.dat`。
- 窒素沈着: `ndepo_1980-2003_tky.txt`。
- GHG: `AtmGHG_timeseries.txt`。
- TKYサイト・パラメータ: 固定ソース側 `Config_TKY.xlsx` の
  `SAVE as ...` sheetからテキストを復元。
- 16バイオームパラメータ: `parameter_VISITc_16.xlsx` の
  `parameter_VISITc_16.txt` sheetから復元。コードが要求する
  `parameter_S1b.txt` はリポジトリに欠けるため、この復元値を同名で使用した。
- 選択経路: `FLUX_SCHEME=0`。ユーザー要件とExcelブックの対象経路に合わせた。
- spinup: `SPUPT=2000`, 気象期間1948–2010年。

実行用複製は `artifacts/visit_native_20000715/run_pinned2013/` に保存した。
`diagnostic_changes.patch` に固定ソースからの全変更、`native_manifest.sha256` に
実行物・入力・診断出力のSHA-256を記録した。同じバイナリを2回実行し、
診断TSVがバイト単位で一致することも確認した。

最終xlsxを生成した後のリポジトリ全体テストは163件成功、7件スキップである。
スキップは任意のnative・依存関係経路であり、ブック検査の失敗ではない。

## 実行で判明した注意点

### 保水容量の上書き

Configの上層保水容量64.21 mmは実行時に使われない。
`init_site.c::f_init_site` は `soil_physics.c::f_soil_saxton` を呼び、
砂・粘土率から上層を104.810272811 mmに上書きする。深層はConfigの全層
808.87 mmから上層を引いた704.059727189 mmである。従来ブックの値では
水の1日更新が一致せず、実行値に修正した後に37状態が一致した。

### TKYの水蒸気圧分岐

`location_proc.c` のCEAMIPサイト判定に、
`strcmp(grid->site_id, "CEAMIP_TMK")` だけ `==0` がない。TKYではこの式が真になり、
通常の気圧と比湿による式ではなく `vp = vps * specific_humidity` が実行される。
基準日のnative値は0.3532280953 hPaである。xlsxは実行同値性を優先してこの値を採用し、
意図されたと思われる24.57 hPaへは無断修正していない。

### 提供2011年版との区別

提供ディレクトリのCソースは2011年版で、同梱バイナリはmacOS Mach-Oのため
Linux上では直接実行できなかった。再コンパイル時は `NFILE=3` に対して
16バイオームを書く配列外アクセスがあり、`NFILE=17` としないと発散した。
修正後は概ね提供出力を再現したが完全一致ではなく、NSC状態もない。
そのため最終xlsxの状態値には用いず、固定2013年版の実行値を採用した。

## 残る未定義量

- `z_ngas2_aa`: `n2o_emit.c::f_n2o_emit_ngas_2` の初期化されない局所変数。
- `z_curry_undefined`: `ch4_oxy.c::f_ch4oxy_curry` の条件の間で値が定義されない分岐。

両者は基準日の37状態更新照合には影響しない診断経路である。
値を推測せず、引き続き未設定とした。
