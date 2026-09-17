# 炭素循環と水循環の完全な数式仕様

> 2026-09-15追記：この文書は前回の実装説明の記録です。
> 日更新の数値積分としての位置づけ、供給制限式の期間水量への換算、
> PM・土壌ポテンシャルの単位上の留保、クリッピング範囲での水理容量の解釈は、
> [日本語レポート（TeX）](carbon_water_matrix_report.tex)／
> [PDF](carbon_water_matrix_report.pdf)で補足・訂正しました。
> 特に、平坦化したポテンシャル関数の逆微分を容量0と解釈することはできません。
> 今回は文書化であり、対応する試作コードの変更は行っていません。

## 1. この文書が区別する二つの数理対象

本リポジトリには、時間表現の異なる二つの対象がある。

1. **native VISITc水文写像**：C原典の逐次一日更新
   \(z_{k+1}=F_{\rm VISITc}(z_k,u_k)\)。これは微分方程式ではない。
2. **4成分炭素–水結合モデル**：日を時間単位とする8状態の連続微分方程式
   \(\dot x=F(x,u)\)。これはVISITcそのものではなく、VISITcと植物水理文献に
   根拠を持つ縮約統合モデルである。

原典の逐次代数を無理にODEと呼ぶと、同日内の計算順、クリッピング、baseflowの
二重控除が消える。そのため本書では、最初にnative写像を完全に示し、その後に
研究用ODEを示す。

## 2. 記号、単位、向き

炭素状態は

\[
\boldsymbol C=(C_l,C_t,C_r,C_s)^\top
\quad[\mathrm{Mg\ C\ ha^{-1}}],
\]

水状態は

\[
\boldsymbol W=(W_l,W_t,W_r,W_s)^\top
\quad[\mathrm{mm}],
\]

時間は日である。添字 \(l,t,r,s\) はleaf、stem、root、soilを表す。
水ポテンシャル \(\psi\) はMPa、内部水輸送 \(q\) はmm day\(^{-1}\)、
水理コンダクタンス \(G\) はmm day\(^{-1}\) MPa\(^{-1}\)である。

正の内部フラックスは

\[
s\to r,\qquad r\to t,\qquad t\to l
\]

の向きとする。計算値が負なら逆流であり、ゼロへ切らない。

## 3. native VISITcの一日水文写像

### 3.1 状態と強制

native状態を

\[
z_k=(S_k,U_k,L_k)
\]

とする。\(S\) は `snwa`、\(U\) は `sw30`、\(L\) は `sww` である。
\(U,L\)を厳密に非重複な上下層と解釈することは原典からは保証されない。

外部入力は降水 \(P\)、気温 \(T_a\)、深層地温 \(T_d\)、および7個の
Penman–Monteith potential fluxである。

### 3.2 雪と雨

降雪割合、雪、雨は

\[
f_{snow}=\frac{1}{1+\exp[0.75(T_a-2)]},\quad
P_s=f_{snow}P,\quad P_r=(1-f_{snow})P.
\]

融雪は

\[
M=\begin{cases}
\dfrac{1/11}{1+\exp[-0.5(T_a-4)]}
\left(1+\dfrac{10}{0.05S_k+1}\right)S_k,&S_k>0.1,\\
S_k,&S_k\le0.1,
\end{cases}
\]

\[
S_{k+1}=S_k+P_s-M.
\]

### 3.3 供給制限演算子

原典で遮断蒸発、土壌蒸発、蒸散に繰り返し使われる演算子を

\[
\mathcal A(Q,E_p)=
\frac{Q+E_p-\sqrt{(Q+E_p)^2-4(0.85)QE_p}}{2(0.85)}
\]

と置く。\(Q\)は利用可能水量、\(E_p\)はpotential fluxである。

### 3.4 遮断蒸発

tree、C3、C4のLAIを \(L_T,L_3,L_4\)、下層占有率を
\(f_3,f_4\) とする。遮断可能量は

\[
Q_T=\min(P_r,0.25L_T),
\]

\[
Q_3=\min[f_3(P_r-I_T),0.25f_3L_3],
\]

\[
Q_4=\min[f_4(P_r-I_T),0.25f_4L_4].
\]

注意すべき原典仕様は、C4でもC3遮断量を引かず、tree遮断量だけを引くことである。

\[
I_j=\mathcal A(Q_j,I_{p,j}),\quad I=I_T+I_3+I_4,
\]

\[
G=P_r-I+M.
\]

### 3.5 上側bucket、蒸発、草本蒸散

VISITcのcubic runoff演算子を

\[
\mathcal R(g,d)=\max\left(
\max(g^3+d^3,0)^{0.33333}-d,0
\right)
\]

とする。上側の乾燥指数とrunoffは

\[
d_U=FC_U-U_k,\qquad R_1=\mathcal R(G,d_U),
\]

\[
U^{(1)}=U_k+G-R_1.
\]

土壌蒸発とその後の状態は

\[
E_s=\max[\mathcal A(U^{(1)},E_{p,s}),0],\qquad
U^{(2)}=U^{(1)}-E_s.
\]

C3とC4は同じ \(U^{(2)}\) を供給量として個別に評価する。

\[
T_3=\mathcal A(U^{(2)},T_{p,3}),\qquad
T_4=\mathcal A(U^{(2)},T_{p,4}),
\]

\[
U^{(3)}=U^{(2)}-T_3-T_4.
\]

このため二つの合計が水量を超え、後段クリッピングが作動する場合がある。

### 3.6 tree蒸散、baseflow、下側runoff

tree蒸散は

\[
T_T=\mathcal A(L_k,T_{p,T}).
\]

baseflowは

\[
B=\begin{cases}
0,&T_d\le0,\\
0.003L_k,&T_d>0\ \text{and site=QHB},\\
0.001L_k,&\text{otherwise}.
\end{cases}
\]

まず

\[
L^{(1)}=L_k-B
\]

とし、

\[
d_L=FC_L-L^{(1)},\quad R_{2,b}=\mathcal R(R_1,d_L),\quad
R_2=R_{2,b}+B.
\]

### 3.7 store間再配分

\[
r_F=\frac{FC_U}{FC_L},\qquad
X=\frac{L^{(1)}r_F-U^{(3)}}{1+r_F},
\]

最大交換量は

\[
X_{max}=K_{sat}\,1000\,3600\,24.
\]

原典の非対称規則は

\[
X_r=\begin{cases}
\min(0.5X,X_{max}),&X>0,\\
\max(X,-X_{max}),&X\le0.
\end{cases}
\]

である。クリップ前状態を

\[
\widetilde U=U^{(3)}+X_r,\qquad
\widetilde L=L^{(1)}-X_r
\]

とし、

\[
\delta_U=\max(-\widetilde U,0),\qquad
\delta_L=\max(-\widetilde L,0),
\]

\[
U_{k+1}=\max(\widetilde U,0),
\]

\[
L_{k+1}=\max(\widetilde L,0)+R_1-R_2-T_T.
\]

最後の式は、クリップ後に再び負値を生じ得ることにも注意が必要である。

### 3.8 native水収支とbaseflow二重控除

\[
AET=I+E_s+T_T+T_3+T_4.
\]

通常期待する残差を

\[
\varepsilon_W=
(S_{k+1}+U_{k+1}+L_{k+1}-S_k-U_k-L_k)
-(P-AET-R_2)
\]

とすると、原典どおりの代数では

\[
\boxed{\varepsilon_W=-B+\delta_U+\delta_L}
\]

となる。baseflowは一度 \(L\) から直接引かれ、さらに \(R_2\) に含まれて
最後にもう一度引かれる。実装はこれを修正せず再現・診断する。

## 4. VISITc Penman–Monteith依存式

### 4.1 飽和水蒸気圧

\[
e_s(T)=\begin{cases}
6.1078\,10^{7.5T/(237.3+T)},&T>0,\\
6.1078\,10^{9.5T/(265.3+T)},&T\le0.
\end{cases}
\]

単位はhPaである。

### 4.2 飽和水蒸気圧曲線の傾き

\[
\Delta(T)=\begin{cases}
\dfrac{6.1078(2500-2.4T)}{0.4615(273.15+T)^2}
10^{7.5T/(237.3+T)},&T>0,\\
\dfrac{6.1078(2834)}{0.4615(273.15+T)^2}
10^{9.5T/(265.3+T)},&T\le0.
\end{cases}
\]

### 4.3 空気密度と空気力学抵抗

気圧 \(p_a\)、水蒸気圧 \(e_a\) をhPaとすると

\[
\rho_a=1.293\frac{273.15}{273.15+T}
\frac{p_a}{1013.25}\left(1-0.378\frac{e_a}{p_a}\right).
\]

風速はまず \(u'=\max(u,0.1)\) とし、

\[
r_a=\operatorname{clip}\left[
\frac{\log^2(10)}{0.41^2u'},0.1,59.5
\right]\quad[\mathrm{s\ m^{-1}}].
\]

### 4.4 葉炭素からLAI、太陽高度から消散係数

葉炭素 \(C_{fol}\) [Mg C ha\(^{-1}\)]、比葉面積 \(SLA\)
[cm\(^2\) gDM\(^{-1}\)] から、native `lai_mass` は

\[
LAI=\max\left(\frac{SLA\,C_{fol}\,(2.2)}{100\times2},0\right).
\]

ここで2.2はcarbonからdry matterへの係数、100は面積・質量基準変換、2は
片面葉面積への変換である。太陽高度 \(\beta\) [degree] から

\[
k=\frac{k_0}{\operatorname{clip}[\sin(0.0174533\beta),0.3,1]}
\]

を `irr_attn` として用いる。

### 4.5 LAIから被覆率

初期消散係数を \(k^0_j\) とすると

\[
f_T=1-e^{-k_T^0 LAI_T},
\]

\[
f_{C3}=(1-f_T)u_3(1-e^{-k_3^0 LAI_3}),
\]

\[
f_{C4}=(1-f_T)u_4(1-e^{-k_4^0 LAI_4}),qquad
f_g=1-f_T-f_{C3}-f_{C4}.
\]

### 4.6 層別純放射

地表面温度 \(T_s\)、気温 \(T_a\)、雲量 \(c\)、水蒸気圧 \(e_a\)、
Stefan–Boltzmann定数 \(\sigma=5.6703\times10^{-8}\) とすると

\[
L_\uparrow=0.95\sigma(T_s+273.15)^4,
\]

\[
\epsilon_a=0.53+0.06\sqrt{e_a},
\]

\[
L_\downarrow=(1-c)\epsilon_a\sigma(T_a+273.15)^4
+c[\sigma(T_a+273.15)^4-9],
\]

\[
L_{net}=L_\downarrow-L_\uparrow,qquad L_j=f_jL_{net}.
\]

透過率を0.1、放射消散係数を \(k_j\) とすると

\[
a_T=1-e^{-0.9k_TLAI_T},
\]

\[
a_3=(1-a_T)u_3[1-0.9e^{-k_3LAI_3}],
\]

\[
a_4=(1-a_T)u_4[1-0.9e^{-k_4LAI_4}],qquad
a_g=1-a_T-a_3-a_4.
\]

日平均短波 \(R_{SW}\) とalbedo \(\alpha_j\) から

\[
S_j=R_{SW}a_j(1-\alpha_j).
\]

原典の演算を文字どおり写すと

\[
\boxed{R_{n,j}=S_j-L_j}
\]

である。ここでは物理的な符号修正を挟んでいない。

### 4.7 LAI・GPP・CO2・VPDから群落コンダクタンス

`f_canopy_cond` の内部GPP proxyは、\(p_{sat}>0\)なら

\[
b=\frac{k\,LUE\,PPFD_{top}}{p_{sat}},
\]

\[
G=\frac{2p_{sat}}{k}
\log\frac{1+\sqrt{1+b}}
{1+\sqrt{1+be^{-kLAI}}},
\]

それ以外では0である。CO2とVPD係数は

\[
f_{CO2}=\frac{1}{CO2-\Gamma},\qquad
f_{VPD}=\frac{1}{1+VPD/b_2}.
\]

群落コンダクタンスは

\[
g_c=g_{s,b0}LAI+g_{s,b1}f_{CO2}f_{VPD}G.
\]

native `FIX_GSCO2=1` の場合だけ \(CO2\) を350 ppmへ固定する。

### 4.8 PM共通式

\(r_c\)を表面抵抗、\(\gamma=0.667\) hPa K\(^{-1}\)、
\(\lambda=695\) W h kg\(^{-1}\)、日長を \(D\) hourとすると

\[
\mathcal{PM}(R_n,c_p,r_c)=
\max\left[
\frac{D}{\lambda}
\frac{\Delta R_n+c_p\rho_aVPD/r_a}
{\Delta+\gamma(1+r_c/r_a)},0
\right].
\]

遮断蒸発では \(r_c=0\)、\(c_p=0.2813\)。蒸散では

\[
r_c=\frac{1}{g_c\eta},\qquad \eta=0.0224/1000,
\]

かつ \(g_c>0,R_n>0\) のときだけ計算する。

土壌蒸発では

\[
g_g=1000\frac{U}{FC_U}+100,qquad r_s=\frac{1}{g_g\eta}.
\]

さらに原典 `pm_evap` は \(c_p=0.2813\) を代入した直後に
\(c_p=1014.0\) で上書きするため、有効値1014.0を再現する。

## 5. 貯水量と水ポテンシャルの構成則

### 5.1 VISITc土壌保水式

土壌水量を \(W_s\)、field capacityを \(W_s^{sat}\) とし、

\[
\widehat W_s=\max(W_s,0.2)
\]

とする。土性ごとの \((a,b)\) は

\[
(0.121,4.05)\ \text{sandy},\quad
(0.478,5.39)\ \text{medium},\quad
(0.405,11.4)\ \text{fine}.
\]

\[
\psi_{m,s}=-a\left(\frac{\widehat W_s}{W_s^{sat}}\right)^{-b},
\qquad
\psi_s=\psi_{g,s}+\psi_{m,s}.
\]

native上層では \(\psi_g=-0.05\) MPa、whole/deep側では
\(\psi_g=-1.00\) MPaである。8状態縮約モデルは後者を既定値とする。

### 5.2 植物pressure–volume曲線

植物器官 \(i\in\{l,t,r\}\) に対して

\[
R_i=\operatorname{clip}\left(
\frac{W_i}{W_i^{sat}},R_{min},1
\right)
\]

を相対含水量とする。\(\pi_{0,i}<0\) は飽和時浸透ポテンシャル、
\(\epsilon_i>0\) は体積弾性率である。turgor loss pointは

\[
R_{tlp,i}=1+\frac{\pi_{0,i}}{\epsilon_i},
\qquad
\psi_{tlp,i}=left(\frac1{\epsilon_i}+\frac1{\pi_{0,i}}\right)^{-1}.
\]

SurEau-Ecosのpressure–volume式を使い、

\[
\psi_i(R_i)=
\begin{cases}
-\pi_{0,i}-\epsilon_i(1-R_i)+\dfrac{\pi_{0,i}}{R_i},
&R_i\ge R_{tlp,i},\\
\dfrac{\pi_{0,i}}{R_i},&R_i<R_{tlp,i}.
\end{cases}
\]

器官capacitanceは

\[
C_i=\frac{dW_i}{d\psi_i}=W_i^{sat}\frac{dR_i}{d\psi_i},
\]

\[
\frac{dR_i}{d\psi_i}=
\begin{cases}
\dfrac{R_i}{-\pi_{0,i}-\psi_i-\epsilon_i+2\epsilon_iR_i},
&\psi_i\ge\psi_{tlp,i},\\
-\dfrac{\pi_{0,i}}{\psi_i^2},&\psi_i<\psi_{tlp,i}.
\end{cases}
\]

これはSurEau-Ecosのsymplasm式をleaf/stem/rootの各集約storeへ適用したもので、
native VISITc式ではない。\(R_{min}\) は \(R\to0\) で
\(\psi\to-\infty\) となるための明示的な数値正則化である。clipが有効な
\(R\le R_{min}\) と \(R\ge1\) では、実装上の局所capacitanceを0と報告する。

## 6. 8状態炭素–水結合ODE

### 6.1 補助関数

\[
h(x;a)=\frac{\max(x,0)}{\max(x,0)+a}.
\]

\[
\theta_i=\frac{W_i}{W_i^{sat}}.
\]

### 6.2 炭素入力と分解

\[
\mu=\mu_{max}f_{photo}(t)
\sqrt{h(\theta_l;a_p)h(\theta_s;a_p)},
\]

\[
f_d=f_{decomp}(t)h(\theta_s;a_d),qquad
R_h=k_sf_dC_s.
\]

配分率は \(a_l+a_t+a_r=1\) とする。

### 6.3 炭素の4本の微分方程式

\[
\boxed{\dot C_l=a_l\mu-k_lC_l},
\]

\[
\boxed{\dot C_t=a_t\mu-k_tC_t},
\]

\[
\boxed{\dot C_r=a_r\mu-k_rC_r},
\]

\[
\boxed{\dot C_s=k_lC_l+k_tC_t+k_rC_r-k_sf_dC_s}.
\]

従って

\[
\boxed{\frac{d}{dt}(C_l+C_t+C_r+C_s)=\mu-R_h}.
\]

### 6.4 炭素量による水理コンダクタンス制御

\[
g_{sr}=G_{sr}^{max}h(C_r;K_r),
\]

\[
g_{rt}=G_{rt}^{max}\sqrt{h(C_r;K_r)h(C_t;K_t)},
\]

\[
g_{tl}=G_{tl}^{max}\sqrt{h(C_t;K_t)h(C_l;K_l)}.
\]

内部フラックスは

\[
q_{sr}=g_{sr}(\psi_s-\psi_r),
\]

\[
q_{rt}=g_{rt}(\psi_r-\psi_t),
\]

\[
q_{tl}=g_{tl}(\psi_t-\psi_l).
\]

ここで各 \(\psi\) は第5節の構成則から得る。

### 6.5 外部水フラックス

葉からの蒸散は

\[
E_l=T_p(t)h(\theta_l;a_T)h(C_l;K_{C,l}),
\]

土壌蒸発は

\[
E_s=E_p(t)h(\theta_s;a_E).
\]

\[
W_{fc}=f_{fc}W_s^{sat},qquad
D=k_D\max(W_s-W_{fc},0),
\]

\[
R=f_RP(t)\,[\operatorname{clip}(\theta_s,0,1)]^2.
\]

### 6.6 水の4本の微分方程式

\[
\boxed{\dot W_l=q_{tl}-E_l},
\]

\[
\boxed{\dot W_t=q_{rt}-q_{tl}},
\]

\[
\boxed{\dot W_r=q_{sr}-q_{rt}},
\]

\[
\boxed{\dot W_s=P-q_{sr}-E_s-D-R}.
\]

内部輸送は総和から厳密に消えるため、

\[
\boxed{
\frac{d}{dt}(W_l+W_t+W_r+W_s)=P-E_l-E_s-D-R
}.
\]

### 6.7 行列表現

\[
\dot{\boldsymbol C}=\boldsymbol u_C(\boldsymbol W,t)
+M_C(\boldsymbol W,t)\boldsymbol C,
\]

\[
\boldsymbol u_C=(a_l\mu,a_t\mu,a_r\mu,0)^\top,
\]

\[
M_C=\begin{pmatrix}
-k_l&0&0&0\\
0&-k_t&0&0\\
0&0&-k_r&0\\
k_l&k_t&k_r&-k_sf_d
\end{pmatrix}.
\]

水側は

\[
\dot{\boldsymbol W}=S_W\boldsymbol q+\boldsymbol p-\boldsymbol e,
\]

\[
S_W=\begin{pmatrix}
0&0&1\\
0&1&-1\\
1&-1&0\\
-1&0&0
\end{pmatrix},\quad
\boldsymbol q=(q_{sr},q_{rt},q_{tl})^\top,
\]

\[
\boldsymbol p=(0,0,0,P)^\top,qquad
\boldsymbol e=(E_l,0,0,E_s+D+R)^\top.
\]

\(\boldsymbol1^\top S_W=0\) が水収支を保証する。

## 7. 平衡、貯留容量、局所縮約

瞬間的に水を固定した炭素容量は

\[
\boldsymbol C_{cap}=-M_C^{-1}\boldsymbol u_C,
\quad X_c=\boldsymbol1^\top\boldsymbol C_{cap}.
\]

現在量 \(X=\boldsymbol1^\top\boldsymbol C\)、符号付き貯留potentialは

\[
X_p=X_c-X.
\]

これは結合平衡とは異なる。結合平衡は

\[
F_C(C^*,W^*;u)=0,\qquad F_W(C^*,W^*;u)=0
\]

を同時に満たす。

Jacobianを

\[
J=\begin{pmatrix}J_{CC}&J_{CW}\\J_{WC}&J_{WW}\end{pmatrix}
\]

とする。水が安定かつ炭素・外部強制より十分速い場合、

\[
\delta W\simeq-J_{WW}^{-1}J_{WC}\delta C,
\]

\[
\boxed{J_{eff}=J_{CC}-J_{CW}J_{WW}^{-1}J_{WC}}.
\]

この近似の成立範囲を、今後のdynamic-water対quasi-steady-water実験で検証する。

## 8. 実装と検証対応

| 数式群 | 実装 | 検証 |
|---|---|---|
| native一日水収支 | `visitc_hydrology.py` | `visitc_hydrology_bridge.c`との全出力比較 |
| PM・大気・放射 | `visitc_hydrology.py` | `hydro_flows.c`, `radiation.c`直接比較 |
| 群落conductance | `visitc_hydrology.py` | `ecophysiology.c::f_canopy_cond`直接比較 |
| 土壌 \(W\leftrightarrow\psi\) | `hydraulic_relations.py` | 3土性の手計算、dry floor、単調性 |
| 植物P–V・capacitance | `hydraulic_relations.py` | turgor-loss連続性、微分の有限差分比較 |
| 8状態ODE | `coupled_carbon_water.py` | C/W収支、非負軌道、平衡、両方向Jacobian、数値精密化 |

native比較用スクリプトは

```bash
PYTHONPATH=src python examples/validate_visitc_hydrology.py /path/to/VISITc/point
```

で実行する。checkoutのHEADと対象C/headerの未変更性はビルド時に検査される。

## 9. 現段階でなお近似である部分

- 8状態モデルの炭素配分、turnover、最大GPPは例示値である。
- 植物3器官をそれぞれ単一symplasmic reservoirとしている。apoplasm、cavitation、
  vulnerability curve、重力高さ、複数土壌層は未導入である。
- plant saturated waterを現段階では固定parameterとしており、炭素量から器官容積へ
  導くallometryは未較正である。
- 8状態ODEの蒸散・蒸発は理想化式であり、native PMを直接forcing生成器として
  接続するAPIはあるが、縮約ODEの全気象駆動実験は次段階である。
- native VISITc水文は逐次日次写像であり、このODEとの直接差を「水の効果」とは
  解釈できない。過程内容と時間離散化を揃えた比較が必要である。
