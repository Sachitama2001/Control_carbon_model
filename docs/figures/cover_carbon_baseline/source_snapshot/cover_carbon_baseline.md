# 被覆率・炭素11状態基準モデル

2026-10-02。分類：**model-agnostic derived analysis**。本書は
[森林・草原被覆率と気候の研究計画](forest_grass_cover_climate_research_plan.md)
の段階Aを式として固定した最初の実装・証明記録である。気象応答、較正、
気候駆動R-tippingは含まない。コードは
[`cover_carbon_baseline.py`](../src/control_carbon/cover_carbon_baseline.py)、
専用検査は [`test_cover_carbon_baseline.py`](../tests/test_cover_carbon_baseline.py)。

## 相図と関数の図集

[3ページの図集PDF](figures/cover_carbon_baseline/cover_carbon_baseline_atlas.pdf)。

1. [被覆率パラメータ相図・位相線](figures/cover_carbon_baseline/01_cover_phase.png)：固定した
   \(\phi_0,\theta,w\) の下で \((\alpha,\phi_1)\) を変えた安定平衡数と、
   基準係数での不安定平衡を境とする二つの吸引域。
2. [採用関数の概形](figures/cover_carbon_baseline/02_function_shapes.png)：
   ロジスティック森林→草原率、草原→森林の対抗項、面積率付き炭素生産、
   生体炭素プールに作用する群別転換損失。
3. [炭素capacity](figures/cover_carbon_baseline/03_carbon_capacity.png)：森林・草原由来の
   5プール別capacityと総量を被覆率に対して表示。安定な被覆平衡では、
   固定被覆capacityが全系の炭素平衡に一致する。

炭素係数は形を説明するための未較正試験値であり、生態学的な推定値ではない。
各図の計算値とパラメータ、ソースハッシュは同じディレクトリの
`plot_data.json`・`manifest.json`に保存する。再生成するには:

```bash
python -m pip install -e '.[analysis]'
python examples/plot_cover_carbon_baseline.py
```

## 今回分かったこと

既存の草原被覆率モデルに、面積転換による炭素の損失・枯死物への移動と、
群別の炭素収支を接続した。固定気象・固定係数下で内部被覆平衡ごとに炭素平衡が
一意に定まり、炭素部分系が安定なら被覆率の安定平衡は全11状態系の安定平衡へ
持ち上がる。既存の三つの内部被覆平衡を使う検算例では、安定・不安定・安定の
平衡がそれぞれ11状態系にも得られた。

## 仮定と状態

状態は \(y=(g,x)\)、\(g\in[0,1]\) と
\[
x=(L_T,S_T,R_T,D_T,H_T,L_G,S_G,R_G,D_G,H_G)^\top
\]
である。\(g\) は草原面積率、\(1-g\) は森林面積率。森林区画と草地区画は
排他的で、プール量は全て地表面積当たり kg C m\(^{-2}\)、計画上の時間規約は年。
ただし既存の被覆係数を較正された年率へ変換したわけではなく、現実の時間尺度は未確定。
T/Gは炭素の由来を示す。\(D\) は枯死物、\(H\) は腐植である。

被覆式は Kumar & Dutta (2026), DOI
[10.1098/rspa.2025.0803](https://doi.org/10.1098/rspa.2025.0803) の
(2.3), (4.1) に対応する既存解析をそのまま再利用する。
\[
\dot g=(1-g)\{\phi(g)-\alpha g\},\quad
\phi(g)=\phi_0+\frac{\phi_1-\phi_0}{1+\exp[-(g-\theta)/w]}.
\]
\(\phi\) は森林から草原への実効面積転換率、\(\alpha(1-g)\) は草原から森林への率。
これらを個体死亡率・燃焼率とは解釈しない。

炭素側は本計画の解析用仮定であり、上記論文・VISIT・SPITFIREの式を転記したものではない。
成熟区画の面積当たりNPPを固定 \(P_T,P_G>0\) とし、
\[
\mu_T=(1-g)P_T,\qquad \mu_G=gP_G,\qquad
r_T=\phi(g),\quad r_G=\alpha(1-g).
\]
各群 \(i\in\{T,G\}\)、生体プール \(j\in\{L,S,R\}\) に対して
\[
\dot C_{ij}=p_{ij}\mu_i-(d_{ij}+r_i)C_{ij},
\]
\[
\dot D_i=\sum_j[d_{ij}+(1-b_{ij})r_i]C_{ij}-k_{Di}D_i,\qquad
\dot H_i=\rho_i k_{Di}D_i-k_{Hi}H_i.
\]
ここで \(p_{ij}>0,\sum_jp_{ij}=1\)、\(d_{ij},k_{Di},k_{Hi}>0\)、
\(b_{ij}\in[0,1]\)、\(\rho_i\in(0,1]\)。\(b_{ij}\) は転換時に大気へ出る割合で、
残りは由来別枯死物へ移る。草原から森林への基準転換では
\(b_{Gj}=0\) とする。転換した面積の炭素密度が区画内で一様という仮定から、
面積流束と炭素密度の積は森林→草原で \(\phi C_{Tj}\)、草原→森林で
\(\alpha(1-g)C_{Gj}\) となる。転換後に元の生体炭素を新しい植生へ移植しない。
topkill、萌芽、葉量ゼロからの再生も表現しない。

すべての炭素係数は明示的な `CarbonParameters` として渡す。実装に生態学的な
既定値は置いていない。テスト中の係数は式・数値恒等式を検査するための
**未較正の試験値**である。

## 炭素予算と不変性

全炭素 \(X=\sum_{ij}C_{ij}+\sum_i(D_i+H_i)\) について、
内部移動は相殺し、
\[
\dot X=\mu_T+\mu_G-R_h-E_{tr},
\]
\[
R_h=\sum_i[(1-\rho_i)k_{Di}D_i+k_{Hi}H_i],\qquad
E_{tr}=\sum_{ij}b_{ij}r_iC_{ij}.
\]
従ってNEP \(=\mu_T+\mu_G-R_h\) と炭素貯留の変化率を区別する。
被覆率は \(g=0\) で \(\dot g=\phi_0>0\)、\(g=1\) で \(\dot g=0\) なので
\([0,1]\) が前方不変。炭素の各成分がゼロの面では生産・流入が非負であり、
炭素の非負直交象限も前方不変である。

有限な係数の下で \(\mu_i\leq P_i\)、\(r_i\) は有界。各生体プールは
\(\dot C_{ij}\le p_{ij}P_i-d_{ij}C_{ij}\) で有界となる。生体プールが有界なら
枯死物への入力も有界で \(k_{Di}>0\) により \(D_i\) は有界、続いて
\(k_{Hi}>0\) により \(H_i\) も有界となる。

## 被覆平衡の炭素系への持ち上げ

固定 \(g\) における炭素系は \(\dot x=q(g)+M(g)x\)。群ごとの行列は
生体3プール→枯死物→腐植の下三角形で、対角要素は
\[
-(d_{ij}+r_i),\quad -k_{Di},\quad -k_{Hi}.
\]
全て負なので \(M(g)\) はHurwitzで可逆。各 \(g\in[0,1]\) に唯一の凍結平衡
\(x^*(g)=-M(g)^{-1}q(g)\) がある。成分表示は
\[
C_{ij}^*(g)=\frac{p_{ij}\mu_i(g)}{d_{ij}+r_i(g)},\quad
D_i^*(g)=\frac{\sum_j[d_{ij}+(1-b_{ij})r_i(g)]C_{ij}^*(g)}{k_{Di}},\quad
H_i^*(g)=\frac{\rho_i k_{Di}}{k_{Hi}}D_i^*(g).
\]

内部被覆平衡 \(g^*\) では \(F(g^*)=0\)、\(F(g)=(1-g)(\phi(g)-\alpha g)\)。
全系ヤコビアンは
\[
J(g^*,x^*)=
\begin{pmatrix}
F'(g^*)&0\\
\partial_g(q+Mx)&M(g^*)
\end{pmatrix},\qquad
F'(g^*)=(1-g^*)[\phi'(g^*)-\alpha].
\]
よってその固有値は \(F'(g^*)\) と \(M(g^*)\) の固有値の和集合である。
\(F'(g^*)<0\) の被覆平衡は全11状態で局所漸近安定、
\(F'(g^*)>0\) の平衡は炭素系が安定でも不安定方向を一つ持つ。
これは三角構造による厳密な持ち上げであり、炭素から被覆へのフィードバックを
加えた場合には成立しない。

さらに \(P_G>0\) の下で
\[
C_{Gj}^*(g)=\frac{p_{Gj}gP_G}{d_{Gj}+\alpha(1-g)}
\]
は \(0<g<1\) で厳密に増加する。従って異なる内部被覆平衡の炭素ベクトルは
異なる。一方、総炭素量の大小は一般には定まらない。森林由来プールと草原由来
プールの相対生産・損失率に依存する。各群の総量は
\[
X_i^*(g)=\sum_j C_{ij}^*(g)+
\left(1+\frac{\rho_i k_{Di}}{k_{Hi}}\right)
\frac{\sum_j[d_{ij}+(1-b_{ij})r_i(g)]C_{ij}^*(g)}{k_{Di}},
\qquad X^*(g)=X_T^*(g)+X_G^*(g).
\]
従って二枝の総量差の符号は \(X^*(g_H)-X^*(g_L)\) の係数条件で決まり、
この一般形だけでは符号を主張できない。

## capacity / potential

\[
c_{cap}(g)=-M(g)^{-1}q(g),\qquad c_{pot}=c_{cap}(g)-x.
\]
これは被覆率を現在値に固定した炭素部分系のcapacityで、被覆率まで平衡化した
全系平衡とは区別する。全系平衡では \(g=g^*\) も満たし \(x^*=c_{cap}(g^*)\)。
ゆえに平衡上のpotentialはどちらの安定枝でもゼロであり、枝の判別には使えない。

## 検証・未解決点・次の作業

`pytest -q tests/test_cover_carbon_baseline.py` は被覆安定平衡の持ち上げ、
炭素行列のHurwitz性、capacity残差、炭素予算、解析ヤコビアンの有限差分一致、
非負性の境界符号を検査する。コードは気象 \(T,R\) をまだ入力せず、
\(P_i,\alpha,\phi_0,\phi_1,\theta,w,d,k,\rho,b\) の気象依存も未定である。
相対湿度、VPD、土壌水分、燃料水分の対応づけやR-tippingを主張しない。

次の成果物は、温度・相対湿度から一つの過程へ至る応答を一次文献で制約する
応答台帳である。その応答が定まる前に気象ランプや臨界速度の探索は行わない。
