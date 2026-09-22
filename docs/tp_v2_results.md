# v2実験結果：spinupと関数形比較

2026-09-22。結果は `artifacts/temperature_precipitation/v2/spinup_comparison/`。
設定は [tp_v2_spinup.json](../configs/tp_v2_spinup.json)、
式と出典、設定方法は [実装説明](tp_v2_spinup_and_functions.md) を参照する。
[日本語PDF](../artifacts/temperature_precipitation/v2/spinup_comparison/report/tp_v2_results.pdf) と
[TeX原稿](tp_v2_results.tex) も保存した。

## 初期化の疑問への回答

前回は一定気象の平衡を連立方程式で直接求めて初期化した。
今回追加した時間積分によるspinupでも、同じ平衡と同じ実験軌跡を再現した。
負になっていた量は `I_C` であり、貯留炭素ではない。

今回の全設定で、全8状態が非負で容量内にあることを確認した。
ゼロへの事後補正による質量追加は0である。
spinupは375～400年で収束。直接平衡解との規格化距離は最大2.66e-7、
前回と同じ式の4設定について、20年間の軌跡差は最大2.13e-7だった。
spinupの収束閾値は有限なので、直接平衡解とのわずかな差は残る。

| 元と同じ関数の設定 | spinup年数 | 前回軌跡との最大規格化差 |
|---|---:|---:|
| 通常ターンオーバーのみ | 400 | 1.19e-7 |
| 温度枯死 | 375 | 2.13e-7 |
| 水分枯死 | 400 | 9.96e-8 |
| 加算枯死 | 375 | 1.82e-7 |

各差は `abs(v2-v1)/max(1,abs(v1))` で全時刻・全8プールを比較したもの。
気象は基準27°C・6 mm/day、実験は31°C・2 mm/dayを用いた0/T/P/TP。
すべて同じspinup終端から分岐する。

![全8プールの非負性](../artifacts/temperature_precipitation/v2/spinup_comparison/baseline_water/nonnegative_pools.png)

負の時刻はspinup終盤、時刻0が実験開始。全軸の下端を0にした。
全375～400年の初期化軌跡は各設定の `spinup_all_pools.png` と `spinup.csv` にある。
気象と符号つき差分も含む図は `states_forcing_and_contrast.png`。

## 関数形の比較

以下は水分依存枯死を共通に使い、それぞれの関数形でspinupした結果。
単位はMg C/ha。すべて仮定パラメータによる理論実験である。

| 変更箇所 | spinup後の全炭素 | 20年後 I_C |
|---|---:|---:|
| 現行式 | 596.430 | −9.901 |
| VISIT型の温度応答 | 596.430 | −28.751 |
| VISIT型の群落光合成 | 620.455 | −33.937 |
| VISITの土壌非気孔制限を追加 | 468.143 | +5.360 |
| VISIT型の蒸発散供給制限 | 608.197 | +7.382 |
| VISIT型のbaseflowを追加 | 596.430 | +1.182 |
| 温度＋群落光合成＋蒸発散供給制限 | 627.608 | −1.233 |

![関数形比較](../artifacts/temperature_precipitation/v2/spinup_comparison/function_comparison.png)

前回の負のI_Cは初期化の失敗ではなかった。一方、その値と符号は過程関数の選択に
依存した。従って現時点では「この相互作用が森林一般で頑健に生じる」とは言えない。
関数の形に加え、各設定の基準炭素・水状態も変わっている。今後は基準点の利得を
揃えた比較や、パラメータ範囲の制約で、この二つの寄与を分ける必要がある。
速度族・吸引域解析は今回行っていないため、R-tippingの有無はこの表から判定しない。

## 検査結果

- 10設定、各4実験、計40本。各本実験を異なるsolverと半分の最大刻みで再計算した。
- 全体テスト190 passed、7 skipped。v2の専用テスト14件。
- 原典Cを実際にコンパイルし、温度応答6条件、GPP12条件、供給制限4条件を照合。
- 受理ステップ・保存点・途中の補間点について、集計した監査点数は1,095,865。
- 負の状態、容量超過、質量のゼロ補正：いずれも0。
- 本実験の刻み・solver比較による最大規格化差：2.41e-8以下。
- 水分枯死の基準設定ではspinupもBDF・最大5日で再計算し、終端差は7.27e-11以下。
- 最大瞬間収支残差：2.00e-15以下（C・Wそれぞれのフラックス単位）。

全監査点での最小値は以下の通り。

| 状態 | C_L | C_S | C_R | C_O | W_L | W_S | W_R | W_O |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 最小値 | 2.29755 | 29.32815 | 7.20841 | 80.00000 | 0.29401 | 1.15068 | 0.68205 | 104.32888 |

詳細は `nonnegative_minima.csv` と各 `spinup_audit.json`、`manifest.json`。
途中点の検査は有限標本であり、連続時間全域の区間証明ではない。
方程式の非負境界・水容量境界の構造テストを併用している。

## 再現とカスタマイズ

```bash
MPLCONFIGDIR=/tmp/control-carbon-mpl PYTHONPATH=src \
  python examples/run_tp_v2.py --config configs/tp_v2_spinup.json
MPLCONFIGDIR=/tmp/control-carbon-mpl PYTHONPATH=src \
  python examples/audit_tp_v2.py --run artifacts/temperature_precipitation/v2/生成された実験名
```

設定ファイルの `baseline`、`target`、`initial_state`、`years`、`spinup`、
`parameters`、`models[].processes` を編集する。
気温は°C、降水はmm/day、期間はyear（365.25日）を用いる。
複数の関数形は同じ `models` 配列で一括比較できる。

今回使った `spinup_comparison` は保持し、次回は新しい出力ディレクトリを使う。
旧v1のコード・結果は上書きしていない。v1のソースハッシュと保存版計画書のハッシュも
照合した。現行計画書には初期化・非負性に関する§4.4だけを追加している。
