# 炭素＋窒素の最小R-tipping検討

2026-09-15。独立したモデル非依存の研究例。VISITの転記や森林の較正ではない。

- [日本語レポート（TeX）](nitrogen_tipping_minimal.tex)
- [日本語レポート（PDF）](nitrogen_tipping_minimal.pdf)
- 実装：`src/control_carbon/minimal_nitrogen_tipping.py`
- 再現：`examples/analyze_minimal_nitrogen_tipping.py`
- テスト：`tests/test_minimal_nitrogen_tipping.py`

固定C:N比を使い、炭素生産に対応する窒素吸収を明示する。
単調な窒素制限だけの2状態モデルでは、最終環境に正の平衡があれば永続的R-tippingはない。
窒素過剰による生産低下を仮定すると、安定な植生枝が残ったまま、急な窒素流入増加で裸地側へ移れる。
急変時の吸引域侵入は不等式で示し、緩変時の追随と合わせて可能性を導出した。
窒素過剰の応答式はAndrews (1968), 式(1)の関数型を借りた仮説であり、森林の証拠ではない。

線形ランプ3→10では、試験係数で移行時間の境界が約27.20年。
時刻の中心をそろえれば、総窒素投入量も同じ条件で速い移行と遅い移行を比較できる。
枯死物を加えた3状態では再循環の遅れが境界を変える。
月・年は連続系の出力間隔であり、年刻みEuler法や年平均への厳密な集約ではない。

## 再現

リポジトリ直下で実行する。

```bash
python examples/analyze_minimal_nitrogen_tipping.py
python -m pytest tests/test_minimal_nitrogen_tipping.py -q
```

`artifacts/minimal_nitrogen_tipping/`に図PDFと`summary.json`を生成する。
LuaLaTeXと日本語フォントを使用し、`docs/`から次を2回実行する。

```bash
lualatex -interaction=nonstopmode -halt-on-error nitrogen_tipping_minimal.tex
```

今後は、窒素不足と過剰のどちらを対象にするかを定め、森林データから関数型・尺度・時間を制約する。
微生物競合や可変C:N比を導入する前に、この否定例と条件付き肯定例を対照として残す。
