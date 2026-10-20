# 有理関数E4のcapacity trackingとISS：Phase 1–3中間報告

2026年10月16日。方針は [connectISS.md](connectISS.md) に従う。
本報告ではcapacityの定義と誤差方程式から出発し、外生入力に対するISSを証明した後、
炭素フィードバックを戻して何が変わるかを調べる。
全系平衡の数値評価は、この区別を検査する補助材料として用いる。

**今回の区切りはPhase 1–3である。** 大規模なVPD走査、既存二軌道の全面的な再解析、
新しい論文草稿の本格執筆はまだ行っていない。
理論式の数値確認に限り、81個の凍結平衡と二つの外生入力の検証問題を用いた。
後者をE4の干ばつ実験と混同しない。

## 0. 根拠の水準と一次資料

- **[A]** 仮定と証明を示した解析的結果。
- **[B]** 外向き丸めなどで誤差を保証した計算機支援結果。今回新規の[B]はない。
- **[C]** 通常の浮動小数点計算・有限時間積分による検証。
- **[D]** 未検証の仮説・今後の課題。

作業開始時のブランチは `main`、コミットは
`cab30ea3a639f5568ee4d32006f1f1b750557caa`（`videos and proofs`）。
一次資料は [有理関数モデル](../src/control_carbon/spac_carbon_rational.py)、
[炭素行列](../src/control_carbon/spac_carbon.py)、
[モデル係数を読む関数](../examples/prepare_manuscript_phase_portraits.py) である。
指数関数版の [既存capacity診断](../src/control_carbon/spac_carbon_capacity.py) と
[全系平衡への追従診断](../src/control_carbon/spac_carbon_tracking.py) は参考にしたが、
前者の指数関数の微分を有理関数版に流用していない。

確認した範囲で、モデル式・係数と最新の有理関数版の記述に新たな不一致はなかった。
既報の図の誤接続および数値的な完全気孔閉鎖の注意は
[既存の訂正・減衰解析](e4_rational_equilibria_and_decay_20261016.md) の扱いを維持した。
既存モデル、旧論文、スライド、保存済み結果は変更していない。

## 1. 瞬間capacityを参照として定義する [A]

全状態と炭素状態を

\[
y=(\psi_x,\psi_l,C_L,C_W,C_R,C_S)^{\mathsf T},\qquad
C=(C_L,C_W,C_R,C_S)^{\mathsf T}
\]

とする。炭素収支とcapacityは

\[
\dot C=MC+bA(y),\qquad
b=(b_L,b_W,b_R,0)^{\mathsf T},\qquad
\boxed{K(y)=-M^{-1}bA(y)}.
\]

\[
M=\begin{pmatrix}
-k_L&0&0&0\\
0&-k_W&0&0\\
0&0&-k_R&0\\
\eta_Lk_L&\eta_Wk_W&\eta_Rk_R&-k_S
\end{pmatrix}.
\]

したがって、全系平衡を解く前に

\[
\boxed{\dot C=M(C-K)}
\]

が厳密に得られる。capacityは、現在の入力を固定した場合の炭素収支の参照である。
現在の状態から計算する \(K(y)\)、水理も停止する全系平衡 \(y^*\)、
実軌道の到達先は別の対象である。
\(C=K(y)\) は炭素の時間微分がゼロであることだけを意味する。
水理が変化すれば、そこからcapacityが動き、直ちに炭素とずれることもある。

有理関数光合成とE4の水理は、時間を日として

\[
A=p\,h(\psi_l)\frac{qC_L}{1+qC_L},\qquad
h(\psi_l)=\max\left(0,1+\frac{\psi_l}{4\ {\rm MPa}}\right),
\]
\[
\begin{aligned}
g_x&=g_{x,\max}\max(0,1+\psi_l/(3\ {\rm MPa})),\qquad g_s=g_{s,\max}h,\\
\frac{d\psi_x}{dt}
&=\frac{2g_x(C_W/C_{W,\mathrm{ref}})(-\psi_x)
 -2g_x\frac{C_L+C_{L,f}}{C_{L,\mathrm{ref}}+C_{L,f}}(\psi_x-\psi_l)}{c_x},\\
\frac{d\psi_l}{dt}
&=\frac{2g_x\frac{C_L+C_{L,f}}{C_{L,\mathrm{ref}}+C_{L,f}}(\psi_x-\psi_l)
 -g_s(C_L/C_{L,\mathrm{ref}})D(t)}{c_l},\\
D(t)&=\frac{\mathrm{VPD}(t)}{0.101325\ {\rm MPa}}
\end{aligned}
\]

である。根水ポテンシャルは0 MPaに固定する。
光合成に \(C_L/C_{L,\mathrm{ref}}\) をもう一度掛けることはしない。

基準値は \(p=1.2\ {\rm kg\,C\,m^{-2}\,yr^{-1}}\)、\(q=5\ ({\rm kg\,C\,m^{-2}})^{-1}\)、
\((b_L,b_W,b_R)=(0.3,0.4,0.3)\)、
\((k_L,k_W,k_R,k_S)=(0.8,0.05,0.2,0.04)\ {\rm yr^{-1}}\)、
\((\eta_L,\eta_W,\eta_R)=(0.7,0.9,0.8)\)、
\(g_{x,\max}=0.0086\ {\rm m\,day^{-1}\,MPa^{-1}}\)、\(g_{s,\max}=0.47\ {\rm m\,day^{-1}}\)、
\(c_x=c_l=0.001\ {\rm m\,MPa^{-1}}\)。
\(C_{L,\mathrm{ref}}=0.2981897047957961\)、\(C_{W,\mathrm{ref}}=6.361380368976984\)、
\(C_{L,f}=0.01\) の単位は \({\rm kg\,C\,m^{-2}}\)。
実装では1年365.25日で換算する。

## 2. 外生入力に対する追従誤差のISS [A]

### 定理1：変分定数法による評価

誤差を \(e=C-K\) とする。\(M\) は定数Hurwitz行列とし、
\(K(t)\) は有限時間区間ごとに絶対連続とする。
ここで初めて \(w=-\dot K\) を独立な外生入力として扱うと、
ほとんど至る所で

\[
\boxed{\dot e=Me+w=Me-\dot K}.
\]

誘導行列ノルムについて \(\|\exp(Mt)\|\leq c\exp(-\alpha t)\)、
\(c\geq1,\alpha>0\) とする。
\(\dot K\) が対象区間で本質的に有界なら、\(\tau=t-t_0\) に対して

\[
\begin{aligned}
e(t)&=e^{M\tau}e(t_0)-\int_{t_0}^t e^{M(t-s)}\dot K(s)\,ds,\\
\boxed{\|e(t)\|
&\leq ce^{-\alpha\tau}\|e(t_0)\|
+\frac{c}{\alpha}(1-e^{-\alpha\tau})
 \|\dot K\|_{L^\infty[t_0,t]}} .
\end{aligned}
\]

**証明。** 第一式は線形誤差方程式の変分定数公式である。
被積分項に半群の上界を適用し、
\(\int_{t_0}^t e^{-\alpha(t-s)}ds=(1-e^{-\alpha\tau})/\alpha\) とすれば得られる。
後半の係数を \(c/\alpha\) で上から押さえれば、通常のISS形式になる。
この定理は誤差の状態 \(e\) と入力 \(-\dot K\) の関係であり、
水理・炭素六次元系の安定性定理ではない。□

帰結は次の通りである。

1. \(K\) が一定なら \(e\to0\)。
2. \(\dot K\) が有界なら誤差は有界であり、入力の大きさに比例する上界を持つ。
3. 本質的上限の意味で \(\dot K(t)\to0\) なら \(e(t)\to0\)。
   任意の後刻から同じ評価を適用し、入力の初期部分の影響を指数的に減衰させればよい。
   これは \(K\) 自身の収束や、特定の全系平衡への収束を意味しない。
4. \(\dot K\in L^1[t_0,\infty)\) でも \(e\to0\)。
   積分を有限の初期部分と小さい \(L^1\) 尾部に分ければ示せる。
5. 一時的な大きな入力の寄与は、その入力が止まった後には指数的に減衰する。

各有限区間で積分可能であることを仮定した上での主張である。
単に \(K\) が有界であることから \(\dot K\) が有界だとは言えない。

### プール別に得られる評価

今回の \(M\) はMetzler行列、すなわち非対角成分が非負であるため、
\(\exp(Mt)\) は成分ごとに非負である。よって絶対値を成分ごとに取ると

\[
|e(t)|\leq e^{M\tau}|e(t_0)|
+\int_{t_0}^t e^{M(t-s)}|\dot K(s)|\,ds.
\]

区間内で \(|\dot K(s)|\leq m\) を成分ごとに満たすなら、

\[
\boxed{|e(t)|\leq e^{M\tau}|e(t_0)|
+(-M)^{-1}(I-e^{M\tau})m}.
\]

この評価はプール間移送を残すため、全プールを一つのノルムで押さえる評価より鋭くなり得る。
ただし時間的な符号の打ち消しや入力と誤差の位相関係は捨てている。
成立する上界であることと、精度のよい予測であることは別である。

### Lyapunov関数による別証明

\(Q=Q^{\mathsf T}\succ0\) に対し

\[
M^{\mathsf T}P+PM=-Q,\qquad P=P^{\mathsf T}\succ0
\]

を解く。\(V=e^{\mathsf T}Pe\) とすれば、

\[
\dot V=-e^{\mathsf T}Qe+2e^{\mathsf T}Pw.
\]

\(z=\sqrt V\) に対する上右微分は

\[
D^+z\leq-\frac{\lambda_{\min}(Q)}{2\lambda_{\max}(P)}z
+\sqrt{\lambda_{\max}(P)}\|w\|_2.
\]

これは \(|e^{\mathsf T}Pw|\leq\sqrt V\sqrt{w^{\mathsf T}Pw}\) と
\(e^{\mathsf T}Qe\geq \lambda_{\min}(Q)V/\lambda_{\max}(P)\) による。
比較定理で積分し、\(\sqrt{\lambda_{\min}(P)}\|e\|_2\leq z\) を用いると定理1の定数を

\[
c=\sqrt{\frac{\lambda_{\max}(P)}{\lambda_{\min}(P)}},\qquad
\alpha=\frac{\lambda_{\min}(Q)}{2\lambda_{\max}(P)}
\]

として選べる。変分定数法と同じISSを別の計量で証明している。
Lyapunov方程式の通常の浮動小数点解は、区間保証された証明[B]とは数えない。

### 基準係数での定数 [A：式、C：小数値]

非重み付き1ノルムでは列和から

\[
\alpha_1=\min\{(1-\eta_L)k_L,(1-\eta_W)k_W,
 (1-\eta_R)k_R,k_S\}=0.005\ {\rm yr^{-1}},\qquad c_1=1
\]

を得る。これは理論的に定めた定数で、図に合わせた値ではない。
\(Q=I\) とした2ノルム評価は数値的に
\(c_2\simeq7.127592603\)、\(\alpha_2\simeq0.01618487409\ {\rm yr^{-1}}\) だった。
\(M\) の最も遅い固有値 \(-0.04\ {\rm yr^{-1}}\) を、
そのまま \(c=1\) のノルム減衰率に置き換えてはいけない。

## 3. 炭素フィードバックを戻す [A]

### capacityの時間微分を展開する

入力の方向は一定なので、定数ベクトル \(v_K=-M^{-1}b\) を一度だけ定義すれば

\[
K=v_KA,\qquad \dot K=v_K\dot A,
\qquad
\frac{v_K}{1\ {\rm yr}}=(0.375,8,1.5,20.25)^{\mathsf T}
\]

である。滑らかな領域では

\[
\begin{aligned}
\frac{\partial A}{\partial C_L}
 &=\frac{pq\,h(\psi_l)}{(1+qC_L)^2},\\
\frac{\partial A}{\partial\psi_l}
 &=\begin{cases}
 \displaystyle\frac{p}{4\ {\rm MPa}}\frac{qC_L}{1+qC_L},&\psi_l>-4\ {\rm MPa},\\
 0,&\psi_l<-4\ {\rm MPa},
 \end{cases}\\
\dot K&=v_K\left(
 \frac{\partial A}{\partial C_L}\dot C_L+
 \frac{\partial A}{\partial\psi_l}\dot\psi_l\right).
\end{aligned}
\]

capacityの定義にはVPDの陽な依存がないため、ここに
\((\partial K/\partial\mathrm{VPD})\,\dot{\mathrm{VPD}}\) を加えない。
VPDは水理の時間発展を通じてcapacityを変える。

\(\dot C_L=-k_Le_L\) を代入すると、

\[
\boxed{\dot e=
\left(M+v_Kk_L\frac{\partial A}{\partial C_L}(1,0,0,0)\right)e
-v_K\frac{\partial A}{\partial\psi_l}\dot\psi_l}.
\]

つまり \(-\dot K\) の中に誤差自身が含まれる。
外生入力として扱ったときの行列 \(M\) と、炭素フィードバックを戻した行列は同じではない。
後者を表示すると

\[
\begin{pmatrix}
-k_L(1-v_{K,L}\frac{\partial A}{\partial C_L})&0&0&0\\
v_{K,W}k_L\frac{\partial A}{\partial C_L}&-k_W&0&0\\
v_{K,R}k_L\frac{\partial A}{\partial C_L}&0&-k_R&0\\
\eta_Lk_L+v_{K,S}k_L\frac{\partial A}{\partial C_L}
 &\eta_Wk_W&\eta_Rk_R&-k_S
\end{pmatrix}.
\]

この下三角構造によって、炭素側の条件を全六次元の固有値計算なしに調べられる。

### 定理2：水ポテンシャル変化率を入力とした炭素誤差の評価

考える軌道の区間全体で

\[
0\leq v_{K,L}\frac{\partial A}{\partial C_L}\leq\rho<1
\]

が成立するとする。\(\rho\) は区間に対する定数の上限である。
上の行列の \(\partial A/\partial C_L\) を \(\rho/v_{K,L}\) で置き換えた
定数行列を \(\overline B_\rho\) とする。
この行列はMetzlerかつHurwitzであり、対角成分は
\(-k_L(1-\rho),-k_W,-k_R,-k_S\) である。

正ベクトル \(s_0=(-\overline B_\rho)^{-1}(1,1,1,1)^{\mathsf T}\)、
その無次元化 \(s=s_0/s_{0,L}\) を用い、
\(\|e\|_{\infty,s}=\max_i |e_i|/s_i\)、
\(\alpha_\rho=\min_i1/s_{0,i}>0\) とする。
また、非正の葉水ポテンシャルでは
\(0\leq\partial A/\partial\psi_l\leq p/(4\ {\rm MPa})\) だから、

\[
\begin{aligned}
\|e(t)\|_{\infty,s}
&\leq e^{-\alpha_\rho\tau}\|e(t_0)\|_{\infty,s}\\
&\quad+\frac{1-e^{-\alpha_\rho\tau}}{\alpha_\rho}
 \max_i\frac{v_{K,i}p}{4\ {\rm MPa}\,s_i}
 \|\dot\psi_l\|_{L^\infty[t_0,t]}.
\end{aligned}
\]

**証明。** 誤差の各成分の絶対値について比較不等式を取り、
状態依存行列を成分ごとに \(\overline B_\rho\) で上から押さえる。
\(\overline B_\rho s_0=-(1,1,1,1)^{\mathsf T}\) なので、
重み付き最大ノルムの対数ノルムは \(-\alpha_\rho\) 以下である。
水理からの入力項のノルムを上式の係数で押さえ、比較不等式を積分すれば得られる。□

この結果は「軌道が条件領域にいる間」の一様な推定である。
その領域が全系の軌道に対して不変であることは、別に示す必要がある。
各時刻の凍結行列がHurwitzであることだけに依存した証明ではない。
数値実装の重みは上式で決め、観測された誤差に合うよう最適化していない。
得られたゲインは保守的で、重みや減衰率を生理学的な時定数とは解釈しない。

葉成分だけなら、より直接に

\[
|e_L(t)|\leq e^{-k_L(1-\rho)\tau}|e_L(t_0)|
+\frac{v_{K,L}p}{4\ {\rm MPa}\,k_L(1-\rho)}
(1-e^{-k_L(1-\rho)\tau})\|\dot\psi_l\|_\infty
\]

と評価できる。

### 条件が成立する場所、成立しない場所

基準係数では

\[
v_{K,L}\frac{\partial A}{\partial C_L}
=\frac{2.25\,h(\psi_l)}{(1+5C_L)^2}
\]

である。ただし数値式の \(C_L\) は \({\rm kg\,C\,m^{-2}}\) 単位で表す。
赤い境界線の下側が1未満の領域である：
[炭素フィードバックの条件領域](../artifacts/spac_carbon/rational_e4_iss_20261016/phase123_v2/e4_carbon_feedback_region.png)。

木部通水停止域 \(\psi_l\leq-3\ {\rm MPa}\) では常にこの比が
\(9/16=0.5625\) 以下なので、\(\rho=9/16\) を選べる。
一方、湿潤・低葉量では1を超える。
したがって定理2を状態空間全体に適用することはできない。
比が1を超えると、葉誤差の凍結係数は正になるが、
それだけで全系の最終的な帰結を判断することもできない。

## 4. 全結合系の安定性との接続

### 正平衡では炭素側の条件は満たされる [A]

正平衡の葉収支 \(b_Lpq h/(1+qC_L^*)=k_L\) より、

\[
\boxed{v_{K,L}\left.\frac{\partial A}{\partial C_L}\right|_{y^*}
=\frac{1}{1+qC_L^*}<1}.
\]

したがって、どの正平衡でも、炭素フィードバックを戻した誤差行列はHurwitzである。
正平衡の近傍でも、十分小さい領域には \(\rho<1\) の余裕を持たせられる。
ただし全系の軌道がその近傍に留まるとは限らない。

### なぜ全系はそれでも不安定になれるか [A：構造、C：E4数値]

水理右辺を二成分の \(G(y,\mathrm{VPD})\) と書く。
誤差の変換の微分と、炭素列・水理列をそれぞれ表示すると、

\[
\delta e=(I-v_K(1,0,0,0)\tfrac{\partial A}{\partial C_L})\delta C
-v_K(0,\tfrac{\partial A}{\partial\psi_l})\delta(\psi_x,\psi_l).
\]

ここで括弧内の4×4行列の行列式は
\(1-v_{K,L}\partial A/\partial C_L\) であり、正平衡では正である。
したがって局所的には \((\psi_x,\psi_l,C)\) と \((\psi_x,\psi_l,e)\) は正則に変換できる。
比が1の場所を跨いで大域的な座標とみなしてはいけない。

この座標で、誤差方程式には依然として

\[
-v_K\frac{\partial A}{\partial\psi_l}
 G_l(\psi_x,\psi_l,C,\mathrm{VPD})
\]

が残る。これを独立入力とすることをやめると、水理・炭素の閉ループが復活する。
線形化でも \(\delta\dot e=B_*\delta e-
v_K(\partial A/\partial\psi_l)_*\,\delta G_l\) となり、
\(\delta G_l\) は炭素にも依存する。
炭素行列 \(B_*\) だけの固有値から全系の固有値は決められない。
実装では6×6の座標変換 \(T\) を作り、\(TJT^{-1}\) と元の \(J\) が相似であることを確認した。

| VPD | \(v_{K,L}\partial A/\partial C_L\) | 葉誤差の凍結減衰率（年\(^{-1}\)） | 全系の最大固有値実部（日\(^{-1}\)） |
|---|---:|---:|---:|
| 2.44 kPa | 0.7458843 | 0.2032925 | \(-4.43623\times10^{-5}\) |
| 3.00 kPa | 0.9239079 | 0.0608737 | \(+0.0108486\) |

どちらの平衡でも炭素誤差行列の固有値は負だが、3.00 kPaの全系は不安定である。
この表は[C]であり、連続区間全体の区間保証やHopf分岐の非線形型の証明ではない。
[部分系と全系のスペクトル比較](../artifacts/spac_carbon/rational_e4_iss_20261016/phase123_v2/e4_carbon_iss_vs_full_feedback.png)も参照。
水理を炭素固定で切り出した行列と、誤差座標で切り出す水理行列は異なるので、
両者のゲインを混ぜない。

### 全系への十分条件は、まだ別に必要である

滑らかな領域の全系平衡において、六次元Jacobianが厳密にHurwitzであることを確認できれば、
通常の局所Lyapunov理論から、小さいVPD摂動に対する局所ISSが得られる。
しかし、今回の浮動小数点スペクトルはその厳密な余裕・近傍半径・ゲインを保証していない。
特に3.00 kPaの不安定平衡には、この十分条件は適用できない。

別の道として、小ゲインによる接続を考えられる。
定理2の誤差ゲインを \(\gamma_C\) とし、同じ座標と同じ不変領域で水理側にも

\[
\|\dot\psi_l\|_{[t_0,t],\infty}
\leq\text{初期値による過渡項}
+\gamma_H\|e\|_{[t_0,t],\infty,s}
+\gamma_D\|\delta\mathrm{VPD}\|_{[t_0,t],\infty}
\]

という適切なISS入出力評価があれば、\(\gamma_C\gamma_H<1\) によって
相互依存する上界を閉じられる。
これは**追加仮定付きの接続方法**であり、E4で水理側評価や小ゲイン条件を証明したわけではない。
元の水理・炭素分割では3.00 kPaで炭素固定の水理行列にも正の固有値があり、
「両部分系が安定」という仮定から始める単純な接続は、そのままでは使えない。

## 5. 状態、capacity、全系平衡の三つの速度 [A]

ここまでの解析を、ユーザーの「状態と平衡点の変化速度」という視点に対応させると、
次を別々に測る必要がある。

1. **実ストックの速度**：\(\dot C=M e\)。
2. **状態依存capacityの速度**：
   \(\dot K=v_K[(\partial A/\partial C_L)\dot C_L+
   (\partial A/\partial\psi_l)\dot\psi_l]\)。
3. **凍結全系平衡の移動速度**：
   正則な平衡枝 \(F(y^*(v),v)=0\) に対し
   \(\dot y^*=-J_*^{-1}(\partial F/\partial v)\dot v\)。
   ここで \(v\) は圧力としてのVPDである。

1と2の関係が今回のISSの対象である。
3を追う誤差 \(y-y^*\) と、capacity誤差 \(C-K(y)\) は違う。
VPDが一定でも実状態が変われば \(\dot K\neq0\) となり得る一方、
固定VPDの平衡枝は動かない。
また、水理の変化は方向を持つので、単純に速度ノルムの大小だけを比べても追従は保証されない。
誤差の履歴を表す畳み込みが必要である。

今回、3の導関数も有理関数の消去式による微分と、Jacobianを解く方法の両方で実装し、
81点で照合した。これは次段階で既存10日・10年軌道へ接続するための準備である。

## 6. 非滑らかな境界でも何が残るか [A]

- \(\psi_l=-3\) MPaでは全水理の古典的Jacobianは一般に定義できないが、
  光合成とcapacityは \(\psi_l\) に関して滑らかなので、capacityの微分は計算できる。
- \(\psi_l=-4\) MPaでは、正の \(C_L\) に対するcapacityの古典的偏微分は一意でない。
  局所Lipschitz性から軌道に沿うcapacityは絶対連続で、
  ほとんど至る所の微分および積分形のISS評価は使える。
  境界そのものの正方向時間微分は
  \(\max(0,\dot\psi_l/(4\ {\rm MPa}))\) を使った方向微分として扱う。
  現モデルでは両伝導度がゼロで \(\dot\psi_l=0\) なので、その値もゼロとなる。
- \(C_L=0\) では、有理関数は非負側で有限の微分を持つ。
  方程式はゼロ炭素境界へ拡張できるが、正炭素の対数座標による数値積分とは区別する。
  \(C_L=0\) 自体はmax関数の折れ点ではない。

したがって、古典的Jacobianが使えないことを理由に定理1を捨てる必要はない。
一方、未定義のJacobianを数値的なゼロで埋めて安定性を判定してはいけない。

## 7. 閉鎖後の減衰とISSの接続 [A：既存結果からの帰結]

基準係数で \(\psi_l\leq-3\) MPaへ到達した後は、その領域が不変で、

\[
0<C_L(t)\leq C_L(t_c)e^{-0.35(t-t_c)},\qquad
0\leq A(t)\leq(pq/4)C_L(t)
\]

である。ここだけ時間は年単位とする。
したがって全成分が同じ入力を参照する \(K=v_KA\) は少なくともこの指数率でゼロへ向かう。
VPDが有界なら、閉鎖域では \(|\dot\psi_l|=O(C_L)\)、
\(|\dot C_L|\leq k_LC_L\) なので

\[
|\dot A|=O(C_L)+O(C_L^2),\qquad
\|\dot K\|=O(e^{-0.35(t-t_c)}).
\]

よって \(\dot K\in L^1\)、\(\dot K\to0\) であり、定理1から \(e\to0\) となる。
これは \(K\to0\) と合わせれば \(C\to0\) を意味する。
**誤差が消えることは、健全な状態へ回復することと同義ではない。**

プール別には、例えば

\[
C_W(t)=e^{-k_W(t-t_c)}C_W(t_c)
+b_W\int_{t_c}^t e^{-k_W(t-s)}A(s)\,ds
\]

である。葉・根・土壌にも対応する式を使うと、この基準係数では
木部は \(O(e^{-0.05(t-t_c)})\)、根は \(O(e^{-0.2(t-t_c)})\)、
土壌は \(O(e^{-0.04(t-t_c)})\) となる。
木部初期量が正ならその正の減衰項が残り、土壌も正の初期量または移送があれば
遅い減衰モードの係数が正になる。
これは「入力が少なくとも0.35年\(^{-1}\)で減衰する」という結果と、
正の移送を含む畳み込みを合わせた結論であり、単独の損失率だけから決めたのではない。
したがって木部・土壌では、capacityが小さくなった後もストックと誤差が長く残り得る。

ゼロ炭素平衡は水理座標に自由度のある連続集合である。
ISSの \(e\to0\) は水理座標を特定せず、特定のゼロ平衡点への収束定理ではない。
既存10日ランプは約11.095日で閉鎖域へ到達した[C]が、
「任意の速いランプがここへ入る」という定理[A]はまだない。

## 8. 限定した数値検証 [C]

再現スクリプトは [analyze_e4_capacity_iss.py](../examples/analyze_e4_capacity_iss.py)、
数学的診断は [spac_carbon_iss.py](../src/control_carbon/spac_carbon_iss.py)、
自動テストは [test_spac_carbon_iss.py](../tests/test_spac_carbon_iss.py) に分離した。
結果は [phase123_v2](../artifacts/spac_carbon/rational_e4_iss_20261016/phase123_v2/) に保存した。

検証内容：

1. 六つの状態で \(\dot C=M(C-K)\) と、
   \(\dot e=Me-\dot K=B(y)e-v_K(\partial A/\partial\psi_l)\dot\psi_l\) を独立経路で照合。
   \(\dot K\) は解析微分を主とし、方向差分でも確認した。
   水理の折れ点−3 MPa、気孔の折れ点−4 MPaと閉鎖側も含めた。
2. 滑らかな点では、直接微分した
   \(\partial K/\partial y=v_K(\partial A/\partial y)\) と
   \(\partial K/\partial y=E_C-M^{-1}J_{C,y}\) を照合。
3. 2.44–3.00 kPaの81正平衡で、誤差座標の正則性、相似変換、
   全系と炭素誤差系のスペクトル、平衡枝の導関数を確認。
4. 同化を外部から
   \(A(t)=[0.2+0.04\sin(2\pi t/T)]\ {\rm kg\,C\,m^{-2}\,yr^{-1}}\)
   と与えた検証問題を \(T=2,20\) 年について100年間積分。
   初期値は \(C(0)=K(0)\)。炭素収支と誤差方程式を別々に積分し、
   \(C-K\) と独立積分した誤差の差、および許容誤差を厳しくした結果との差は
   最大 \(2.4\times10^{-11}\ {\rm kg\,C\,m^{-2}}\) 未満だった。
   成分別・1ノルム・Lyapunovの各上界への超過は \(10^{-9}\) の許容値内でゼロだった。

今回の最終実行では、解析capacity変化率と方向差分の最大差は
\(2.73\times10^{-10}\ {\rm kg\,C\,m^{-2}\,day^{-1}}\)、
独立なcapacity勾配の最大差は \(7.11\times10^{-15}\)、
閉ループ誤差恒等式の最大残差は \(3.47\times10^{-18}\) だった。
関連する自動テストは**123件通過**し、既存capacity・全系平衡追従・有理関数E4の回帰も含む。
これらは丸め誤差を包含した証明[B]ではなく、数値検証[C]として記録する。

[外生入力の追従とISS上界](../artifacts/spac_carbon/rational_e4_iss_20261016/phase123_v2/exogenous_capacity_tracking_validation.png)では、
上段に葉ストックとcapacity、下段に誤差と上界を示す。
上段は比較のため両方とも最初の20年、下段は100年全体である。
2ノルム上界はかなり保守的なため、下段はゼロ付近を線形表示する対称対数軸を用いる。
この周期入力の例は定理1の検証であり、VPDを指定した非線形E4の軌道ではない。
その誤差の大きさからE4の回復速度や限界速度を推定しない。

再現コマンド：

```bash
python examples/analyze_e4_capacity_iss.py \
  --output artifacts/spac_carbon/rational_e4_iss_20261016/NEW_RUN
pytest -q tests/test_spac_carbon_iss.py tests/test_spac_carbon_rational.py \
  tests/test_spac_carbon_capacity.py tests/test_spac_carbon_tracking.py \
  tests/test_e4_rational_diagnostics.py
```

各実行では出力先を新規作成し、既存成果を上書きしない。
全ての検証問題の四プール時系列、capacity、誤差、capacity変化率、ISS上界はNPZに、
モデル係数・平衡六状態・解析微分・残差・ライブラリとコミット情報はJSONに保存する。
ソースと出力のSHA-256を記録している。

## 9. 一般化に向けて、今言える範囲

**定数Hurwitz行列の場合 [A]。**
炭素プール数や非線形入力の形によらず、
\(\dot C=MC+u(y,t)\)、\(K=-M^{-1}u\)、絶対連続性があれば定理1は成立する。
平衡点の陽な表示や全列挙は必要ない。成分別評価にはMetzler性が追加で必要である。
定理2の三角構造は今回のE4入力と移送に依存し、そのまま一般化はできない。

**時間・状態依存行列の場合 [A：恒等式、D：モデルごとの安定性確認]。**
\(M(y,t)\) が可逆なら \(MK=-u\) の微分から

\[
\dot K=-M^{-1}(\dot u+\dot M K),\qquad \dot e=M(y,t)e-\dot K.
\]

恒等式自体は残るが、各時点の \(M\) がHurwitzというだけでは、
時変誤差系に定理1と同じ一様減衰があるとは限らない。
状態遷移行列の一様指数安定性、または共通の対数ノルム評価などが必要である。
時間依存計量を使うなら \(M^{\mathsf T}P+PM+\dot P\) を調べる必要がある。
capacityの微分・有界性には、逆行列と係数・入力の微分の制御も必要になる。

PINNやサロゲートによる近似は、この追加条件や全平衡点の個数を自動的には証明しない。

## 10. 結果と未解決課題の対応表

| 主張・課題 | 水準 | 仮定・範囲 | 実装での検査 | 今回の記載位置／論文での配置候補 |
|---|---|---|---|---|
| 厳密なcapacity・誤差方程式 | A | 定数 \(M\)、絶対連続な入力 | 直接収支・解析微分・閉ループ式 | 第1–2節／Model・Exact tracking |
| 外生入力ISS | A | \(M\) Hurwitz、\(\dot K\) の可積分性・有界性 | 独立積分、三種の上界 | 第2節／ISS properties |
| Lyapunov定数の小数値 | C | \(Q=I\)、通常の数値線形代数 | 正定値・方程式残差 | 第2・8節／ISS properties |
| 水理変化率を入力とした炭素誤差ISS | A | 区間全体で \(v_{K,L}\partial A/\partial C_L\leq\rho<1\) | 比較行列・重み・対数ノルム | 第3節／Nonlinear feedback |
| 正平衡で炭素誤差行列がHurwitz | A | 正の炭素平衡、現入力構造 | 81平衡・微分式 | 第4節／Nonlinear feedback |
| 3.00 kPaで全系は不安定 | C | 基準係数の凍結正平衡 | 六次元・変換後のスペクトル | 第4節／補助的な平衡診断 |
| 閉鎖域の炭素・capacity・誤差減衰 | A | 閉鎖域への到達、非負有界VPD | 既存閉鎖テスト、今回の微分検査 | 第7節／Hydraulic closure |
| 10日・10年軌道の包括的ISS比較 | D | 同じ初期値・振幅・保持期間 | 次段階で既存軌道を再解析 | Phase 4／Drought transients |
| 水理側のISSゲインと小ゲイン条件 | D | 同じ座標・不変領域が必要 | 今回未証明 | Phase 3以降／Feedback・Limitations |
| 大域的回復域、限界速度 | D | 現状では仮定・証明不足 | 二軌道から結論しない | Discussion・Limitations |
| 非線形ISS中心の独立論文 | 未着手 | この中間報告の後に構成 | 既存TeXは変更しない | Phase 6 |

## 11. 次の区切り

まずPhase 4として、既存二軌道について実状態、capacity、平衡枝の**三つの速度**を分けて測る。
capacityの変化率を解析微分で求め、各プール誤差とISS上界を同じ時間軸で比較する。
閉鎖域への到達と「追従誤差が小さいこと」を別の指標として扱う。
その後、閉鎖後のモード分解と水理側の局所評価を詰め、
新規の [論文草稿予定](connectISS.md) の構成に進む。
大域的な安定性・限界速度の証明が先に済んでいるとは扱わない。
