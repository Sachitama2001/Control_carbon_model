# 最小炭素–水モデルのR-tipping研究ノート

[本文PDF](water_tipping_minimal.pdf) / [TeX](water_tipping_minimal.tex)

今回の位置づけは、ユーザーの依頼に基づく独立した理論検討です。
既存のVISITc・8状態モデルの較正や不具合修正を完了したものではありません。
既存コードの係数は変更していません。

## 結果

- 水利用を `W*C^2` とするKlausmeier型の最小モデルを採用。森林の自己維持に正の密度依存があるという仮説で、森林一般の事実ではありません。
- 準定常水の1状態モデルでは、静的森林平衡が残る間の単調な降水減少だけではR-tippingを起こせないことを証明。
- 降水 `P=P0*lambda` と水損失 `e=e0*lambda^2` の同時移動では、有限終点までの臨界速度を変数分離から導出。これは降水減少シナリオではありません。
- 例示値 `a=2.05, m=0.1/year, Lambda=2` で `r_crit=0.00507617318728/year`。森林の予測値ではありません。
- 動的水を残した2状態では降水のみの減少でも転換例を確認。ただし静的限界から0.1%上、`m/e0=1.8`という限定的条件です。
- 月・年の区間流れ写像、原ODE、解析閾値を照合。専用11件、全体129件のテストが成功。
- 移動の始終端を滑らかにしても転換例は残ります。ただし2状態例では同じ移動期間でも結果が変わり、変化の形への感度が大きいことを確認しています。

## 再現

リポジトリのルートで実行します。

```bash
python -m pytest tests/test_minimal_water_tipping.py -q
python examples/analyze_minimal_water_tipping.py
```

`numpy`、`scipy`、図の生成には`matplotlib`が必要です。
`artifacts/minimal_water_tipping/summary.json`に数値と条件、同じ場所に図PDFを生成します。
TeXの図はこれらを参照するため、TeXだけを別の場所にコピーした場合は図パスの調整が必要です。

`docs/`でLuaLaTeXを2回実行します。

```bash
lualatex -interaction=nonstopmode -halt-on-error water_tipping_minimal.tex
lualatex -interaction=nonstopmode -halt-on-error water_tipping_minimal.tex
```

キャッシュへの書込み制限がある場合は、[既存のビルド手順](carbon_water_matrix_report_build.md)と同様に
`TEXMFVAR`と`TEXMFCACHE`を`/tmp/carbon_water_report_build/`配下へ指定してください。

## 原典境界

VISITa参照コミット：`5c513196f21c1b1efb9ead540e3d40865b15e07e`。
公開ソースを`/tmp/control_carbon_visita_source`へ取得して確認しましたが、全球モデルは実行していません。
この理論モデルにはVISITaの係数を移植していません。
月次ループ、水収支、枯死の実験フラグと符号の監査はTeX本文第8節にあります。

## 次の判定事項

森林の更新・生残が低炭素量で自己維持できなくなる根拠と、水・損傷の緩和時間を調べること。
水分依存枯死の候補式は提示していますが、その式を使ったR-tippingは未検証です。
季節軌道からのR-tipping、全吸引域の証明、実林分の臨界速度は未実施です。
