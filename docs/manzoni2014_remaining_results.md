# Manzoni (2014)：図3・図4・Appendix Aの再現結果

2026-10-07。[第一段階の図2再現](manzoni2014_reproduction.md)の続き。
対象は DOI [10.1002/2013WR015236](https://doi.org/10.1002/2013WR015236)。
原著の図を転載せず、モデル式から計算した独立の図を保存した。
既存のnative VISIT、炭素・森林草原被覆モデルを変更・結合していない。

## 1. 完了範囲と出典の境界

| 対象 | 今回できたこと | 原著との数値一致を主張しない部分 |
|---|---|---|
| Fig. 3a | 原著に記載されたconductance比の中央値78 MPaによる根/VPD閾値。39、156 MPaは明示した感度ケース | 元観測データの百分位5曲線と頻度分布 |
| Fig. 3b | 比78 MPa、根−0.5 MPa、D=0.01〜0.06での形質閾値 | 原著の観測点・図のピクセル一致 |
| Fig. A1 | 原著captionの根0 MPa、半減値(−3,−2) MPa、D=0.02で線形/sigmoidを比較 | A〜Dの形状指数の厳密な指定がcaptionにないため、明示した(4,12),(12,12),(4,4),(12,4)を使う。観測点を補わない |
| Fig. A2 | sigmoidの有限平衡が0、1、2、3個となる相図と全2状態系の安定性 | low/high Dやunsafe形質の具体値がcaptionだけでは確定しないため、下記の明示したケースを使う |
| Fig. 4 | 土壌→根の連続条件、損傷係数0.6、活動枝・fold、土壌係数感度 | **ユーザー承認の仮定値による感度実験**であり、原著Fig. 4の数値再現ではない |

実装：[manzoni_extensions.py](../src/control_carbon/manzoni_extensions.py)。
検査：[test_manzoni_extensions.py](../tests/test_manzoni_extensions.py)。
実行：[reproduce_manzoni2014_remaining.py](../examples/reproduce_manzoni2014_remaining.py)。
元の線形モデルでは流束・貯留収支・Jacobianの共通部分だけを抽出し、既存検査を維持した。

## 2. 図3：活動平衡を失う閾値

図2で検証した活動foldの式を、根境界・形質の組合せへ展開した。
Fig. 3aの中央値以外の二曲線は、原著の百分位ではなく中央値の半分/倍の感度ケース。
conductance比を増すと同じ根境界での臨界Dが小さくなり、破綻域が拡大する。

Fig. 3bでは各Dと気孔半減値に対して木部半減値を解き、
foldのDと設定Dの差を検査した。最大誤差は
\(2.11\times10^{-12}\) mol/mol。内部fold境界の直前は活動平衡が二つ、
直後はゼロであることも確認した。
気孔のゼロconductance境界や形質順序が限界となる部分を、
内部foldとは区別してプロットした。

根−0.5 MPaの補助線は根境界を示す線であり、
気孔の完全閉鎖条件そのものではない。
線形式の完全閉鎖は \(\psi_l=2\psi_{50,s}\) である。
「活動平衡ゼロ」は非通水の全平衡集合が存在しないという意味ではない。

## 3. Appendix A：sigmoidと有限平衡数

論文(A1),(A2)を実装した：
\[
g_x(\psi_l)=\frac{g_{x,\max}}{1+(\psi_l/\psi_{50,x})^{\alpha_x}},\qquad
g_s(\psi_l)=\frac{g_{s,\max}}{1+(\psi_l/\psi_{50,s})^{\alpha_s}}.
\]
負のポテンシャル同士の比は非負なので、非整数の正の形状指数も扱う。
コンダクタンスはlog-domain/expitで評価し、冪のoverflowを避ける。
葉ポテンシャル0で指数<1の場合のJacobianは特異であるため、明示的にエラーとする。

### 探索範囲を切らない根分離

\(u=-\psi_l,\ R=-\psi_r,\ x=|\psi_{50,x}|,\ s=|\psi_{50,s}|,\)
\(t=g_{s,\max}D/g_{x,\max}\) とすると、
\[
F(u)=(u-R)\{1+(u/s)^{\alpha_s}\}
      -t\{1+(u/x)^{\alpha_x}\}=0,\qquad u>R.
\]
これは高々5項の一般化多項式であり、\(z=\log u\) に変換すると指数和となる。
同じ指数の係数を合算し、最小の指数をくくり出す。
その導関数は一項少ない指数和になるため、再帰的に停留点を分離する。
停留点間の単調区間で根を括り、停留点上の接する根も扱う。

最小/最大指数項が他項の絶対値の和より大きくなる位置を計算し、
両端の無限遠の裾を処理する。固定の「−100 MPaまで」や有限グリッドの符号変化で
平衡数を決めない。これは浮動小数点の根分離であり、区間演算による厳密数値証明ではない。
表現可能範囲を超える有限根は明示的なoverflowとして扱う。

平衡を乾燥側まで数えると、今回の細分格子では \(|\psi_l|\) が
\(1.32\times10^{13}\) MPaに達する例もある。
専用の境界近傍検査ではさらに約 \(1.81\times10^{127}\) MPaの形式的な有限根を確認した。
これらは**一定VPD・経験的sigmoidを外挿した数式上の根**であり、
実在の植物水分状態として解釈しない。
Appendix Aが説明する、極端な乾燥で葉内の水蒸気駆動力を一定VPDと見なせない問題に対応する。

### 裾と安定性

給水/蒸散の比は乾燥極限で \(u^{1+\alpha_s-\alpha_x}\) に比例する。
- \(\alpha_s<\alpha_x-1\)：比がゼロへ向かい、0または2個の有限平衡が生じ得る。
- \(\alpha_s>\alpha_x-1\)：比が大きくなり、1または3個が生じ得る。
- 等号では主係数も重要で、指数の比較だけでは決めない。

元の2×2 Jacobianを実装・有限差分で確認した。
極端な乾燥根ではコンダクタンスと固有値がunderflowするので、
平衡でのJacobianを正の \(g_x\) で割った行列でも判定する。
これは固有値の符号を保存するが、絶対的な時間尺度を与えない。
同時にlog(給水/蒸散)を確認し、両流束が数値的にゼロになることを
平衡残差が小さい証拠として利用しない。

### A1比較とA2相図

A1の線形モデルの活動安定平衡は−0.98560376 MPa。
sigmoidのA〜Dの安定平衡は順に
−1.11273183、−1.09226027、−1.03444368、−1.02300010 MPa。
Dの形状では鞍点もあり、−4.38898035 MPa。

線形モデルの−6 MPaは木部のゼロconductance境界であり、
その先は気孔もゼロになる。元の2状態系では第一段階で説明した平衡退化があるため、
これを孤立した2次元鞍点とは呼ばない。
A1の第二比較図は「線形ゼロconductance境界とsigmoid鞍点との差」として表示した。
白い領域は有限のsigmoid鞍点がない領域、色は大きな裾の差も見えるsymlog表示。

A2では根−0.5 MPa、D=0.01/0.04、
safe半減値(−3,−2) MPaとunsafe半減値(−1.5,−2) MPaを使った。
\(\alpha_x,\alpha_s\in[1,15]\) を刻み0.25、各57点で検査した。
4条件合計12,996組の有限平衡数は0〜3個。
活動安定枝・鞍点・乾燥側安定枝がある場合、全Jacobianの安定性が順に
stable/saddle/stableとなる。
半減値がsafe順序でも、sigmoidの裾まで含めた大域的一意性は自動的に成立しない。

## 4. 図4：土壌連続条件と固定の水理損傷

本文§2の尺度関係を、明示的な係数を必要とする形で実装した：
\[
g_{sr}(\psi_{\rm soil})=g_{sr,\rm sat}
(\psi_{\rm sat}/\psi_{\rm soil})^2,\qquad \psi_{\rm soil}\le\psi_{\rm sat}<0.
\]
飽和より湿潤側を勝手にclampしない。
根節点の瞬時連続条件から
\[
\psi_r=\frac{g_{sr}\psi_{\rm soil}+2g_x\psi_x}{g_{sr}+2g_x}.
\]
この根ポテンシャルを2貯留ODEへ代入する。
状態依存の根境界なので、Jacobianには
\(\partial\psi_r/\partial\psi_x,\partial\psi_r/\partial\psi_l\) の両方を含める。
平衡での給水は直列抵抗の形：
\[
E_{\rm supply}=\frac{g_xg_{sr}}{g_x+g_{sr}}
(\psi_{\rm soil}-\psi_l).
\]
根に新しい貯留状態は追加していない。土壌は処方した境界であり、土壌水量ODEもない。

損傷係数0.6は本文§3.3に従い、
\(g_{x,\max}\mapsto0.6g_{x,\max}\)、
\(\psi_{50,x}\mapsto0.6\psi_{50,x}\) とする。
気孔形質、容量は変えず、損傷・回復の時間発展は追加しない。
基準の木部/気孔半減値(−3,−2) MPaは、損傷で(−1.8,−2) MPaへ移る。

ユーザー承認の感度値は
\(\psi_{\rm sat}=-0.1\) MPa、
\(g_{sr,\rm sat}=0.03,0.1,0.3\) m/day/MPa。
これらを原著の較正値とは呼ばない。
土壌−0.5 MPaでの損傷後foldは：

| \(g_{sr,\rm sat}\) [m/day/MPa] | 臨界D [mol/mol] | 葉fold [MPa] |
|---:|---:|---:|
| 0.03 | 0.01115087603 | −3.18191569 |
| 0.10 | 0.01506316443 | −3.01483483 |
| 0.30 | 0.01719337400 | −2.90534515 |

土壌抵抗を増すと臨界Dが低下する。
土壌conductanceを十分大きくすると処方根境界モデルへ収束することも検査した。
3D図は活動枝とfoldを表示し、非通水平衡の連続集合を孤立した安定点として描いていない。

## 5. 検証結果

- 新規拡張28件、第一段階27件、既存の水理・被覆炭素11件：**66 passed**。
- 整数指数ケースで、指数和の根分離と独立な通常多項式の全根を比較した。
  非整数、二重根、無限遠に近い裾、同指数の係数合算も検査した。
- sigmoidと土壌結合の全Jacobian、内部水移動の相殺、根連続条件を検査した。
- sigmoidの二つの安定先と、土壌結合のfold前後の活動/乾燥限界への時間積分を検査した。
  時間積分の許容誤差細分化でも到達先・終状態が保たれた。
- A2刻み0.5→0.25の共通3,364点で平衡数が**すべて一致**。
  A1の共通点の平衡差も数値一致した。これは未サンプルの相境界を厳密に決めた証明ではない。
- Appendixの最大 \(|\log(E_{\rm supply}/E)|\)：\(2.55\times10^{-12}\)。
- 土壌結合の最大ODE平衡残差：\(1.98\times10^{-14}\) MPa/day。
- 土壌foldのゼロ固有値誤差：最大 \(4.24\times10^{-12}\) day\(^{-1}\)。
- 生成図を目視確認した。全リポジトリのテストは今回実行していない。

## 6. 再実行と成果物

```bash
python examples/reproduce_manzoni2014_remaining.py --shape-step .25 --output artifacts/manzoni2014/new_remaining_run
python -m pytest -q tests/test_manzoni_extensions.py tests/test_manzoni_hydraulics.py tests/test_hydraulic_relations.py tests/test_cover_carbon_baseline.py
```

出力先が既存ならエラーとし、上書きしない。前回と同じ依存だけを使い、追加インストールなし。
最終結果：[保存ディレクトリ](../artifacts/manzoni2014/remaining_final_20261007/)。

- [図3 PDF](../artifacts/manzoni2014/remaining_final_20261007/figure3_thresholds.pdf)
- [Appendix A1 PDF](../artifacts/manzoni2014/remaining_final_20261007/appendix_A1.pdf)
- [Appendix A2 PDF](../artifacts/manzoni2014/remaining_final_20261007/appendix_A2_counts.pdf)
- [図4感度実験 PDF](../artifacts/manzoni2014/remaining_final_20261007/figure4_soil_sensitivity.pdf)
- [格子・曲線・平衡差・残差・形質・土壌係数・ソースSHA-256・環境](../artifacts/manzoni2014/remaining_final_20261007/results.json)

刻み0.5の初回結果は `artifacts/manzoni2014/remaining_20261007/`、
細分化の中間結果は `artifacts/manzoni2014/remaining_refined_20261007/` に保持した。
保存先はgit-ignoredであり、ローカルに永続保存している。

## 7. 残る制限と研究への接続

原著の観測形質データ、Fig. 3百分位、Fig. 4の土壌係数、
絶対時間尺度を決める容量は未同定のまま。
原著の全図・全観測点を厳密に数値再現したとは主張しない。
形質・根/土壌境界・容量の実測に照らして使用範囲を制限することが、
炭素モデルとの接続に先立つ課題になる。

今回検証したのは凍結条件でのfold・吸引域・有限平衡数・土壌抵抗による閾値変化。
外力経路と速さを変えるR-tipping、個体死亡、回復、気温/RHからNPPや被覆転換への
写像は実装・実証していない。再現できた水理破綻を、根拠なしに森林→草原率へ換算しない。
