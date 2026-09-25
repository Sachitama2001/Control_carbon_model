# 入れ子で採用式を読むHTML

閲覧するファイルは **`artifacts/model_explorer/index.html`**。
Windowsでは `D:\ct\Control_carbon_model\artifacts\model_explorer\index.html`。
VS CodeのHTMLビューアでも、この生成済みファイルを開く。
`docs/model_explorer/template.html` は生成用の雛形であり、閲覧用ではない。

## 読み方

1. `x ｜ 何が蓄えられているか` の三角をクリックする。
2. 炭素 → 植物 → 樹木 → 葉、とその場で入れ子を展開する。
3. 状態の更新式の下で、各流量、依存する量、採用式をさらに開く。
4. `B ｜ どの箱へ配るか` → 植物への炭素配分 → 動的配分 → 目標LAI → 具体式、と進む。

式の記号リンクは該当見出しへ移動する。展開には各見出しの三角を使う。
親の説明は残り、ページの差し替えは行わない。検索機能は削除した。
展開にJavaScript、ネット接続、サーバーは不要。HTML一つをコピーして閲覧できる。
VS Code拡張によってはHTML標準機能にも制限があるため、その場合は生成済みHTMLを通常のブラウザで開く。

## 式と範囲

旧VISIT TKY日次経路の37状態を既存ワークブックから取得する。
各状態の更新式と流量式は既存Graphの式を埋め込み、再計算・改変しない。
分類の式、モデル固有の式、固定配分の説明用仮定を区別する。
原典抜粋、単位、出典をHTML内に保存する。

深い依存式は「式中の記号と定義」付録にすべて収録する。
式の識別子リンクで移動できる。ビューアが閉じた付録を自動展開しない場合は、付録を開いて参照する。
IFなどの関数名は既存式木の表記。未設定の入力を0にはしない。
数値シミュレーションや別モデルへの計算切替は含まない。

## 再生成と検証

リポジトリルートから実行する。

```bash
python examples/build_model_explorer.py
python -m pytest -q tests/test_model_explorer.py
```

既存workbookの依存関係と、固定commitの `../VISIT-matrix/visit_local` が必要。
原典・出力先は `--source` と `--output` で指定できる。
出力HTMLとmanifestを上書きするが、Excelや研究モデルは変更しない。

- 説明・出典：`examples/build_model_explorer.py`
- 入れ子構造と採用式の抽出：`examples/render_nested_explorer.py`
- HTMLの外枠と見た目：`docs/model_explorer/template.html`

任意の実ブラウザ検証はPlaywrightを使用する（閲覧者には不要）。
`PLAYWRIGHT_MODULE` と `PLAYWRIGHT_BROWSERS_PATH` で一時環境を指定可能。

```bash
node tests/browser_model_explorer.cjs artifacts/model_explorer/index.html
```

JavaScriptとネット接続を無効にし、状態・流量・配分・原典の入れ子展開、
キーボード操作、390px幅、深い定義へのリンクを検査する。
検証JSON・スクリーンショットは生成HTMLと同じフォルダへ保存する。
利用者のVS Code拡張そのもののテストではない。
