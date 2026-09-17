# 日本語レポートのビルド

本文は [carbon_water_matrix_report.tex](carbon_water_matrix_report.tex)、
閲覧用は [carbon_water_matrix_report.pdf](carbon_water_matrix_report.pdf) です。
LuaLaTeX、LuaTeX-ja、本文で指定する標準LaTeXパッケージを使用します。
文献一覧はTeX内にあるためBibTeX/Biberは不要です。

`docs/` で以下を2回実行します（目次・式番号・文献参照の解決）。

```bash
lualatex -interaction=nonstopmode -halt-on-error carbon_water_matrix_report.tex
```

キャッシュの既定ディレクトリに書き込めない環境では、書込み可能な場所を指定します。
今回の検証では以下を使用しました。中間ファイルはリポジトリに含めていません。

```bash
mkdir -p /tmp/carbon_water_report_build
env TEXMFVAR=/tmp/carbon_water_report_build/texmf-var TEXMFCACHE=/tmp/carbon_water_report_build/texmf-cache lualatex -interaction=nonstopmode -halt-on-error -output-directory=/tmp/carbon_water_report_build carbon_water_matrix_report.tex
```

初回はフォントキャッシュの生成に時間がかかる場合があります。
2回目以降も参照の再実行警告が残る場合は、警告が消えるまで再実行します。

本稿は共通の数理枠組み、原典の監査、月・年単位の解析方針を記述したものです。
既存の試作コードの境界処理や物理単位を修正したリリースではありません。
今回の原典照合で見つかった未解決事項は本文第8節に記載しています。
