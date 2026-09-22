# Temperature–precipitation tipping research plan

Last updated: 2026-09-22

## 1. One-sentence objective

葉・幹・根・土壌に対応する四つの炭素プールと四つの水貯留から始め、
気温と降水を単独または同時に変化させたとき、単独強制では現れない
平衡構造、安定性変化、履歴依存性、または rate-induced tipping が
両気象変数の相互作用によって生じるかを明らかにする。

## 2. Scientific question and contribution

中心的な問いは次である。

> 温度と降水を独立または同時に変化させたとき、単独強制では現れない
> 平衡構造、安定性変化、吸引域移動、経路依存性、または R-tipping が、
> 二つの気象変数の相互作用によって生じるか。

「複合ストレスで炭素量がより減る」という量的増幅だけでは不十分である。
本研究で特に区別するのは次の三段階である。

1. **量的非加算性**: 複合応答が単独応答の和と一致しない。
2. **質的構造変化**: 平衡点数、安定性、fold、basin boundary、回復可能性が変わる。
3. **非自励 tipping**: 同じ始点と終点を持つ気象経路でも、速度または経路順序により
   異なる吸引状態へ到達する。

この順序を守り、一時的な炭素減少、source 化、長い回復時間だけを
`tipping` と呼ばない。

## 3. Position within this repository

この計画は、既存の二つの研究レイヤーを置き換えない。

| レイヤー | 目的 | この計画との関係 |
|---|---|---|
| ネイティブ VISIT/VISITc 監査 | ソース順序、状態、収支、係数を忠実に再現する | 式・パラメータ範囲・過程候補を得る検証系 |
| 既存8状態炭素–水ODE | 植物内部水貯留とSPAC輸送を明示する | 将来の高解像度水理拡張 |
| 本計画の8状態ODE | 二気象ドライバーの低次元理論を解析する | 当面の主たる理論研究系 |

本計画のモデルは `VISIT`、`VISITc`、または特定サイトの校正モデルとは呼ばない。
**source-informed, low-dimensional theoretical model** と記述する。

## 4. Scope fixed for the first model

### 4.1 State variables

炭素状態は四プールとする。

\[
\boldsymbol C=
\begin{pmatrix}
C_L & C_S & C_R & C_O
\end{pmatrix}^{\mathsf T},
\]

ここで、\(L\) は葉、\(S\) は幹、\(R\) は根、\(O\) は集約した土壌有機炭素である。
炭素と同じ四区画について、水貯留を

\[
\boldsymbol W=
\begin{pmatrix}
W_L & W_S & W_R & W_O
\end{pmatrix}^{\mathsf T}.
\]

とする。\(W_L,W_S,W_R,W_O\) はそれぞれ葉、幹、根、土壌の水貯留である。
初期モデル全体は、

\[
\boldsymbol x=(\boldsymbol C^{\mathsf T},\boldsymbol W^{\mathsf T})^{\mathsf T}
\in\mathbb R_{\ge0}^{8}
\]

の8状態である。水ポテンシャル、窒素、リン、個体数、サイズ分布は第一段階には
動的状態として入れない。水ポテンシャルを必要とする場合も、まず水貯留と炭素量
から計算する代数診断とし、状態次元を増やさない。

### 4.2 External inputs

時間変化させる外部気象入力は、

\[
\boldsymbol u(t)=\begin{pmatrix}T(t)&P(t)\end{pmatrix}^{\mathsf T}
\]

だけとする。光、\(\mathrm{CO_2}\)、土壌特性、栄養条件、風、火災、伐採、
土地利用は固定する。VPDは独立の第三入力にせず、必要な温度効果を
蒸発散関数内に集約する。この集約は実データ適用時の限界として明記する。

### 4.3 Time semantics

- 方程式は連続時間ODEとする。
- 外部入力は月別の区分一定値、または月別値を滑らかに補間したものとする。
- 数値積分器の内部刻みは1か月より短くする。
- 年次値は診断用の集計値であり、年次更新式を基本モデルにしない。
- \(P\) の月降水量は、月長を明示してODEの水フラックス単位へ変換する。

月次温度で表す枯死は、数時間の葉組織壊死ではなく、慢性的高温ストレスの
月間ハザードと解釈する。Thermal Death Time のような急性熱損傷を扱う場合は、
日次またはそれより細かい強制へ分岐する。

## 5. Initial ecological domain

対象は、東部から南部アマゾンの季節乾燥を経験する熱帯常緑林を模した
**理想化サイト**とする。第一候補の地理的文脈は Tapajós 周辺である。
ただし第一段階ではサイト観測への最適化を行わず、
`Amazon seasonal-forest-type idealized vegetation` と呼ぶ。

この選択の理由は次である。

- 明瞭な乾季があり、降水操作と長期乾燥の研究蓄積がある。
- 高温と乾燥を別々にも複合的にも考えられる。
- 単一の熱帯常緑広葉樹型として低次元化しやすい。
- 後にフラックス、森林動態、throughfall-exclusion 実験と比較できる。

Amazon の降水除外実験は、乾燥による生産・樹冠・死亡応答と長期調整の双方を
示している。ただし、観測された死亡機構を一つの土壌水分閾値、または一つの
植物器官水分閾値へ直結させない。
根拠と適用範囲は [Nepstad et al. (2002)](https://doi.org/10.1029/2001JD000360)、
[Rowland et al. (2015)](https://doi.org/10.1038/nature15539)、
[Sanchez-Martinez et al. (2025)](https://doi.org/10.1038/s41559-025-02702-x)
を起点に整理する。

## 6. Baseline carbon equations

### 6.1 Photosynthesis

固定光条件下の総炭素入力を、

\[
G(T,\boldsymbol W,C_L)
=G_{\max} f_T^G(T)f_W^G(\boldsymbol W,\boldsymbol C)
\left(1-\exp[-k_L C_L]\right)
\]

とする。最後の項は葉量に対する受光・群落光合成の飽和を表す。
第一候補の温度関数は解析しやすい山型関数、

\[
f_T^G(T)=
\exp\left[-\frac{(T-T_{\mathrm{opt}})^2}{2\sigma_T^2}\right],
\]

水分制限は、

\[
f_W^G(\boldsymbol W,\boldsymbol C)
=\frac{r_G^{n_G}}{r_G^{n_G}+r_{50,G}^{n_G}},
\qquad
r_G=\sqrt{r_Lr_R}
\]

とする。ここで \(r_L,r_R\) は後述する葉・根の相対水貯留であり、幾何平均は
どちらか一方の強い不足を隠しにくい最小候補である。後に最小値、律速型、
peaked-Arrhenius 型やVISIT由来関数へ交換しても、実験設計は
変更しない。

### 6.2 Allocation, turnover, and respiration

固定配分ベクトルを

\[
\boldsymbol a=(a_L,a_S,a_R)^{\mathsf T},
\qquad a_i\ge0,\qquad \sum_i a_i=1
\]

とする。植物維持呼吸は、

\[
R_i(T,C_i)=r_{i,0}q_R(T)C_i,
\qquad
q_R(T)=Q_{10,R}^{(T-T_0)/10}
\]

を基準とする。一定 \(Q_{10}\) の感度試験として、全球観測に基づく
exponential-of-quadratic 応答も比較する
([Huntingford et al. 2017](https://doi.org/10.1038/s41467-017-01774-z))。

通常ターンオーバー率を \(\tau_i\)、枯死率を
\(m_i(T,\boldsymbol W,\boldsymbol C)\) として、

\[
\begin{aligned}
\dot C_L &=a_LG-R_L-(\tau_L+m_L)C_L,\\
\dot C_S &=a_SG-R_S-(\tau_S+m_S)C_S,\\
\dot C_R &=a_RG-R_R-(\tau_R+m_R)C_R,\\
\dot C_O &=
\sum_{i\in\{L,S,R\}}(\tau_i+\eta_i m_i)C_i
-R_h(T,W_O,C_O).
\end{aligned}
\]

\(\eta_i\) は枯死炭素のうちモデル化した土壌プールへ入る割合である。
\(1-\eta_i\) は明示的な系外流出として収支表に残す。

土壌呼吸は、

\[
R_h(T,W_O,C_O)=r_{O,0}q_h(T)f_W^h(W_O)C_O
\]

とする。乾燥が光合成と分解を同時に抑えるため、炭素入力効果と滞留時間効果が
競合し得る。

### 6.3 Carbon budget

全炭素 \(C_\Sigma=C_L+C_S+C_R+C_O\) について、

\[
\dot C_\Sigma
=G-\sum_i R_i-R_h
-\sum_i(1-\eta_i)m_iC_i
\]

が数値誤差内で成立しなければならない。内部配分、ターンオーバー、土壌への
死亡炭素移動は全系収支では相殺される。

## 7. Four-component water balance

### 7.1 Storage and relative water content

植物器官の水保持容量を炭素量に比例させる最小形として、

\[
W_i^{\mathrm{cap}}(C_i)=\omega_i C_i+\varepsilon_i,
\qquad
r_i=\frac{W_i}{W_i^{\mathrm{cap}}(C_i)},
\qquad i\in\{L,S,R\}
\]

とする。\(r_i\) は無次元の相対水貯留である。土壌については、

\[
r_O=\frac{W_O-W_{\mathrm{dry}}}
{W_{\mathrm{fc}}-W_{\mathrm{dry}}}
\]

を第一候補とし、使用時は \([0,1]\) の範囲と超過水をフラックス側で処理する。
炭素量が水容量を変えるため、\(C\to W\) の結合が最小形でも残る。

### 7.2 Internal transport and external fluxes

土壌→根、根→幹、幹→葉の流量を \(q_{OR},q_{RS},q_{SL}\) として、

\[
\begin{aligned}
\dot W_L &=q_{SL}-E_L-L_L^W,\\
\dot W_S &=q_{RS}-q_{SL}-L_S^W,\\
\dot W_R &=q_{OR}-q_{RS}-L_R^W,\\
\dot W_O &=P(t)-q_{OR}-E_O-D-Q
+\sum_{i\in\{L,S,R\}}\zeta_iL_i^W.
\end{aligned}
\]

\(E_L\) は蒸散、\(E_O\) は土壌蒸発、\(D\) は深部排水、\(Q\) は表面流出である。
\(L_i^W\) は器官ターンオーバーまたは枯死に伴う植物水の放出、\(\zeta_i\) は
そのうち土壌へ戻る割合である。第一候補は
\(L_i^W=(\tau_i+m_i)W_i\) とするが、炭素ターンオーバーと組織水放出を同じ時定数に
置く仮定を感度試験する。残り \((1-\zeta_i)L_i^W\) は系外流出である。
内部流量は例えば、

\[
q_{ij}=g_{ij}\,\phi(r_i-r_j)
\]

とし、\(\phi\) には滑らかな非負の流量制限を用いる。これは水ポテンシャル差に
基づくDarcy型輸送の縮約であり、同一視はしない。代替として、一定容量近似から
\(\psi_i=\psi_i(W_i,C_i)\) を代数的に計算し、
\(q_{ij}=K_{ij}(\psi_i-\psi_j)\) とする感度系を用意する。

蒸散と土壌蒸発の最小形は、

\[
E_L=E_{L,0}q_E(T)f_E(r_L)A_L(C_L),
\qquad
E_O=E_{O,0}q_E(T)f_O(r_O),
\]

\[
A_L(C_L)=1-\exp(-k_EC_L)
\]

とする。温度は直接的な熱ストレスだけでなく、蒸散・蒸発を通じて間接的に
植物と土壌の水貯留を低下させる。

内部輸送は全水収支で相殺され、

\[
\frac{d}{dt}(W_L+W_S+W_R+W_O)
=P-E_L-E_O-D-Q
-\sum_{i\in\{L,S,R\}}(1-\zeta_i)L_i^W
\]

が成立しなければならない。各貯留の物理的上限、利用可能水、降水超過分は
フラックス制限へ入れ、solver後の無名のclippingで収支を隠さない。

## 8. Mortality hierarchy

枯死過程は一度に追加せず、次の順序で導入する。

1. 通常ターンオーバーのみ。
2. 温度依存枯死 \(m_T(T)\)。
3. 水分依存枯死 \(m_W(r_P)\)。
4. 両者を加算した基準モデル。
5. 明示的な温度–水分相互作用。
6. 必要な場合だけ損傷記憶または非構造性炭素を追加する。

### 8.1 Chronic temperature mortality

月次モデルの基準候補を、

\[
m_T(T)=
m_{T,\max}
\frac{1}{1+\exp[-k_T(T-T_{50,m})]}
\]

とする。\(T_{50,m}\) は月代表温度に対する経験的パラメータであり、葉の短時間
加熱実験で得た \(T_{50}\) ではない。葉代謝の高温限界はバイオーム間で変わり、
曝露時間にも依存する
([O'Sullivan et al. 2017](https://doi.org/10.1111/gcb.13477);
[Faber et al. 2024](https://doi.org/10.1093/jxb/erae096))。

急性熱損傷を後に扱う場合は、

\[
\log_{10}\tau_f(T)=a-bT,
\qquad
\mathcal D_T(t)=\int_0^t\frac{ds}{\tau_f(T(s))}
\]

という Thermal Death Time 型へ切り替える。月次平均値からこの積分を
復元したとは主張しない。

### 8.2 Water-stress mortality

全植物の相対含水量を、

\[
r_P=
\frac{W_L+W_S+W_R}
{W_L^{\mathrm{cap}}+W_S^{\mathrm{cap}}+W_R^{\mathrm{cap}}}
\]

とする。枯死関数の最小候補は、

\[
m_W(r_P)=
m_{W,\max}
\frac{1}{1+\exp[k_W(r_P-r_{50,m})]}
\]

である。器官のうち最も乾燥したものを重視する
\(r_P=\operatorname{softmin}(r_L,r_S,r_R)\) も感度試験に含める。
これは水ポテンシャルやPLCを解かない経験的縮約である。Amazon の干ばつ死亡には
水理破壊の強い証拠があるため、\(\boldsymbol W\) を機構そのものではなく、
死亡ハザードを予測する集約状態として扱う
([Rowland et al. 2015](https://doi.org/10.1038/nature15539))。

### 8.3 Interaction structure

最初の複合モデルは、

\[
m_i(T,\boldsymbol W)=m_{i,0}+\alpha_{i,T}m_T(T)+\alpha_{i,W}m_W(r_P)
\]

とする。この段階で複合応答が非加算的なら、それは光合成、呼吸、蒸発散、
炭素状態を介した**創発的相互作用**である。

次にのみ、

\[
m_i(T,\boldsymbol W)=m_{i,0}+\alpha_{i,T}m_T+alpha_{i,W}m_W
+\gamma_i m_Tm_W
\]

を試す。\(\gamma_i\) は観測または感度範囲なしに調整して tipping を作らない。
Amazon の長期個体記録では、干ばつと熱の双方が死亡増加および遅延効果と
関連したが、これは上式の特定形を一意に決めない
([Aleixo et al. 2019](https://doi.org/10.1038/s41558-019-0458-0))。

## 9. Factorial forcing experiments

相互作用を定義するため、三実験ではなく基準を含む四実験を一組にする。

| ID | Temperature | Precipitation | Purpose |
|---|---|---|---|
| 0 | \(T_0\) 固定 | \(P_0\) 固定 | 基準平衡・基準周期軌道 |
| T | \(T(t)\) 変動 | \(P_0\) 固定 | 温度ドライバーの総効果 |
| P | \(T_0\) 固定 | \(P(t)\) 変動 | 降水ドライバーの総効果 |
| TP | \(T(t)\) 変動 | \(P(t)\) 変動 | 複合効果 |

「温度のみ」は外部降水を固定する意味であり、温度依存蒸散・蒸発を介した
\(\boldsymbol W\) の
変化を禁止しない。この総効果に加え、機構帰属のため次の pathway-clamped
実験を行う。

- \(\boldsymbol W\) を固定し、温度の直接的な光合成・呼吸・死亡効果だけを測る。
- \(m_T\) を固定し、温度の蒸発散経由効果だけを測る。
- \(m_W\) を固定し、乾燥の光合成・分解効果だけを測る。

診断量 \(Y\) の量的相互作用を、

\[
I_Y(t)=Y_{TP}(t)-Y_T(t)-Y_P(t)+Y_0(t)
\]

と定義する。炭素量、NEP、死亡フラックス、回復時間、最小固有値余裕について
計算する。ただし \(I_Y\ne0\) だけでは tipping と判定しない。

## 10. Frozen-system analysis in the \((T,P)\) plane

各固定 \((T,P)\) に対して、

\[
F(\boldsymbol x;T,P)=0
\]

を解き、次を記録する。

- 正の生存平衡、低炭素平衡、境界平衡の存在。
- ヤコビアン \(J=D_xF\) の固有値と安定性。
- fold、transcritical、Hopf の候補。
- 安定平衡間の不安定状態または basin boundary。
- 全炭素、NEP、呼吸、死亡率、土壌水分。
- frozen carbon capacity と完全な結合平衡の差。

二変数を同時に扱う価値は、座標軸上の一次元分岐図だけでなく、
\((T,P)\) 平面におけるfold曲線、安定性境界、可能ならcuspなどのcodimension-2
構造を調べられる点にある。

## 11. Hypotheses

以下は結論ではなく、反証可能な仮説である。

- **H0 — no qualitative interaction:** 基本的な加算死亡モデルでは複合応答は
  量的に非加算でも、凍結系は一意な安定平衡を保つ。
- **H1 — emergent compound threshold:** 温度による光合成低下・呼吸増加・蒸発散増加と、
  水分による光合成低下・死亡増加が共有状態を介して結合し、単独強制軸にはない
  平衡消失または吸引域境界を生む。
- **H2 — order dependence:** 同じ最終 \((T_1,P_1)\) でも、昇温後に乾燥する経路と
  乾燥後に昇温する経路で、最終状態または最大損傷が異なる。
- **H3 — rate dependence:** 凍結系の生存枝が経路上で安定なままでも、十分速い
  二変数変化が basin boundary crossing を生む領域がある。
- **H4 — memory requirement:** H2またはH3が8状態モデルに現れない場合、器官別水貯留
  だけでは必要な記憶が不足し、NSC、損傷炭素、順化、または第二土壌層が必要である。

H0が支持されても重要な結果である。どの最小過程まで tipping が不可能かを示す
ことは、複雑モデルで見つかった挙動の機構帰属に使える。

## 12. Nonautonomous forcing and R-tipping tests

### 12.1 Constant and seasonal references

まず定常入力の平衡を調べる。月別季節性を入れた後は、年平均平衡ではなく
安定な1年周期軌道を参照軌道とする。

### 12.2 Path families

少なくとも次を比較する。

1. \(T_0\to T_1\) のみ。
2. \(P_0\to P_1\) のみ。
3. 温度と降水の同時直線経路。
4. 昇温後に乾燥する折れ線経路。
5. 乾燥後に昇温する折れ線経路。
6. 同一経路・同一端点で速度だけを変える族。
7. 気象を基準へ戻す有限パルス。

### 12.3 Operational R-tipping criteria

`R-tipping` 候補と呼ぶため、次をすべて確認する。

1. 凍結系の参照吸引子または周期軌道族を定義した。
2. 経路上で追跡対象の吸引子が消滅するB-tippingではない。
3. 始点と終点を固定し、速度の違いで長期的結果が変わる。
4. basin boundary、edge state、または接続軌道に相当する証拠がある。
5. 最終保持時間を延ばしても単なる遅い回復ではない。
6. solver、時間刻み、初期条件、閾値定義に対して頑健である。

変化速度を変えると累積曝露時間も変わる。したがって、

\[
\mathcal B_C(t)=
\int_0^t\left[R_a(s)+R_h(s)-G(s)\right]_+ds
\]

のような累積炭素赤字と、乾燥・高温ハザードの時間積分を併記する。
速度効果と単純なdose効果を同一視しない。

## 13. State-expansion gates

8状態を維持することを既定とし、次のいずれかを満たすときだけ状態を増やす。

| Candidate state | Add only if |
|---|---|
| 非構造性炭素 \(C_{NSC}\) | 呼吸赤字から死亡までの遅れ、または低備蓄のbasin boundaryが中心仮説に必要 |
| 機能葉／損傷葉 | 熱損傷と回復の履歴を炭素状態だけで保持する必要がある |
| 第二土壌水層 | 根域深度または長期乾燥記憶が結論を変える |
| 順化状態 | 慢性高温への順化速度がR-tipping仮説の中心になる |

状態追加前に、追加しないモデルの失敗指標を保存する。

## 14. Work packages and acceptance criteria

### WP1 — eight-state baseline

- 方程式、単位、符号、外部流入出を固定する。
- 非負領域と炭素・水収支を検査する。
- 正の基準平衡と安定性を確認する。
- 温度・水分関数を無効化した極限が手計算と一致する。

**Pass:** 収支残差、平衡残差、刻み細分化誤差が事前閾値を満たす。

### WP2 — one-driver controls

- 温度のみ、降水のみの平衡枝と時間応答を求める。
- 各ドライバーの作用経路をclamp実験で分ける。
- 一意安定系なら、その範囲で「tippingなし」を結果として保存する。

**Pass:** 単独強制の安定性と応答が複合実験前に再現可能である。

### WP3 — compound frozen-system map

- \((T,P)\) 平面で平衡を継続する。
- 固有値、fold候補、吸引域を地図化する。
- 相互作用死亡項なし／ありを分ける。

**Pass:** 質的変化が数値初期値探索だけでなく、継続計算と安定性で確認される。

### WP4 — path and rate experiments

- 同時・逐次経路を比較する。
- 同一端点の速度族を計算する。
- dose、最大ストレス、経路長、保持時間を記録する。

**Pass:** R-tipping基準を満たすか、満たさない理由を明示できる。

### WP5 — ecological grounding

- Amazon型パラメータ範囲を文献から設定する。
- 一点最適化より範囲・感度・無次元群を優先する。
- 選択した関数形に対する代替式で結論を確認する。

**Pass:** 各パラメータに単位、範囲、出典、校正／仮定区分がある。

### WP6 — only justified extensions

- 事前に記録した失敗を解決する最小状態だけを追加する。
- 8状態モデルとの入れ子比較を行う。
- 複雑化でのみ現れる挙動の機構を特定する。

## 15. Required outputs

各実験は少なくとも次を保存する。

- forcing path \(T(t),P(t)\) と単位。
- 全状態、全外部フラックス、全枯死フラックス。
- 炭素・水収支残差。
- NEPと総炭素 \(C_\Sigma\)。
- 平衡残差、固有値、最大実部、条件数。
- \(I_Y\) 相互作用項。
- cumulative carbon deficit と stress dose。
- solver、許容誤差、最大刻み、refinement比較。
- モデルレベル、パラメータ由来、git commitを含むmanifest。

生成物は既存規約に従い `artifacts/temperature_precipitation/` 以下へ置き、
意図的に選んだ小さな表・図以外はgitへ追加しない。

## 16. Non-claims

第一段階では次を主張しない。

- 熱帯林全体またはTapajósの将来死亡率を予測した。
- 月平均気温から急性葉熱損傷を再現した。
- 四つの水貯留が水理破壊を機構的に再現する。
- 複合応答の非加算性だけでR-tippingを示した。
- 理論モデルがVISITまたはVISITcの縮約として一意に導出された。
- 特定のパラメータ集合で tipping が出ることを現実の証拠とした。

## 17. Immediate implementation order

1. この計画と既存コードの対応表を作る。
2. 8状態方程式と純粋なフラックス関数を実装する。
3. 単位、非負性、炭素・水収支のテストを先に作る。
4. 定常 \((T_0,P_0)\) の平衡とヤコビアンを実装する。
5. 枯死なしで温度のみ・降水のみ・複合の4実験を実行する。
6. 温度枯死、水分枯死の順に追加する。
7. \((T,P)\) 平衡面と分岐候補を計算する。
8. その後にのみ速度・経路依存実験を行う。
9. R-tippingが成立しなければ、no-tipping範囲を記録して状態追加gateを評価する。

## 18. Core references for this track

| Reference | Use in this plan |
|---|---|
| [Huntingford et al. (2017)](https://doi.org/10.1038/s41467-017-01774-z) | 植物呼吸の瞬間温度応答と順化の代替式 |
| [O'Sullivan et al. (2017)](https://doi.org/10.1111/gcb.13477) | バイオーム横断の葉代謝高温限界とthermal safety margin |
| [Faber et al. (2024)](https://doi.org/10.1093/jxb/erae096) | 温度×曝露時間、加算的熱損傷、TDT拡張 |
| [Aleixo et al. (2019)](https://doi.org/10.1038/s41558-019-0458-0) | Amazon樹木死亡と干ばつ・熱・遅延応答の観測的背景 |
| [Nepstad et al. (2002)](https://doi.org/10.1029/2001JD000360) | 東部Amazonの降水除外と生産・生物地球化学応答 |
| [Rowland et al. (2015)](https://doi.org/10.1038/nature15539) | 熱帯林干ばつ死亡における水理破壊の機構的制約 |
| [Sanchez-Martinez et al. (2025)](https://doi.org/10.1038/s41559-025-02702-x) | 長期降水除外に対するAmazon森林の調整とモデル時間尺度 |
| [Ashwin et al. (2012)](https://doi.org/10.1098/rsta.2011.0306) | B/N/R-tippingの区別 |
| [Feudel (2023)](https://doi.org/10.5194/npg-30-481-2023) | basin boundaryと不安定状態によるR-tipping診断 |

文献の式を実装する場合は、文脈引用ではなく、式、単位、対象組織、時間尺度、
校正範囲、変更点をコード側のprovenanceへ記録する。
