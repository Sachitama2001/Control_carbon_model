# E4 有理関数版：相図・nullcline・capacity 診断

## 対象と実験条件

この診断は、光合成を有理関数

\[
A=p\,\frac{g_s(\psi_l)}{g_{s,\max}}\frac{qC_L}{1+qC_L}
\]

とした**E4だけ**を対象にする。元の指数関数版および既存係数は変更していない。

動的実験は、2.44 kPaの安定な正平衡から始め、VPDを滑らかなstep
\(s^2(3-2s)\)で3.00 kPaまで上げ、1年間保持し、同じ形で2.44 kPaへ戻す。
上昇・下降時間をそれぞれ10日と10年にした二例を比べる。
25°Cの飽和VPD約3.167 kPaより低い範囲である。
帰還後は400年間積分する。ただし葉炭素が \(10^{-6}\) に達した場合はそこで数値積分を止める。
この下限到達は死亡の定義ではなく、停止後の運動も計算していない。

## 「原義」のcarbon storage capacity

炭素収支を \(\dot C=MC+u\) としたとき、ここで描くcapacityは各時刻・各状態に対して

\[
K(t)=-M^{-1}u(t),\qquad
u(t)=(b_L,b_W,b_R,0)^{\mathsf T}A(t)
\]

と定義する従来の瞬間的な凍結係数capacityである。\(M\) は定数で、\(A(t)\) は実際の葉水ポテンシャルと葉炭素量を使う。
したがって \(K(t)\) は**VPDだけの関数ではない**。VPDが実状態を動かし、その状態変化を通じてcapacityが変わる。
また \(K(t)\) と実炭素量 \(C(t)\)、capacityの時間微分、全系の平衡枝は別の量である。
全系平衡では \(K=C\) となる一方、過渡状態では一般に一致しない。

## 図と読み方

図はそれぞれ PNG・PDF・SVG を同じ出力フォルダに保存した。

1. [E4平衡枝と安定性](../artifacts/spac_carbon/rational_e4_diagnostics_20261016/figures_v3/e4_equilibrium_branch_and_stability.png)  
   横軸がVPD、縦軸が正の全系平衡の \(C_L^*\)。平衡曲線は基準係数で単調減少する。
   色分けと黒菱形は有限点スペクトル走査での局所安定性・複素固有値対の横断を示す。
   2.44–3.00 kPaの灰色帯には最初の数値的横断（約2.8453 kPa）が含まれる。
   これは有限格子から絞り込んだ局所スペクトル診断であって、連続区間全体の安定性証明や周期軌道の証明ではない。
   二つ目の横断は約4.792 kPaで、25°Cの飽和VPDより高い。

2. [葉水ポテンシャル–葉炭素 nullcline の断面](../artifacts/spac_carbon/rational_e4_diagnostics_20261016/figures_v3/e4_leaf_nullclines_slices.png)  
   青線は \(\dot\psi_l=0\)、橙破線は \(\dot C_L=0\)。薄い背景色は葉炭素の増減符号、星は全六次元系の平衡である。
   表示される二平面は、各星の \(\psi_x,C_W,C_R,C_S\) を固定して切り出したもの。
   したがって、2次元閉鎖系の相図でも、大域的な回復境界でもない。
   3.00 kPaの星は全系として局所不安定な平衡に対応する。

3. [水理 nullcline の断面](../artifacts/spac_carbon/rational_e4_diagnostics_20261016/figures_v3/e4_water_nullclines_slices.png)  
   青線が \(\dot\psi_x=0\)、橙破線が \(\dot\psi_l=0\)。各パネルでは炭素4プールをそのVPDの全系平衡値に固定する。
   星は二つの水理nullclineの交点であり、対応する六状態全系平衡でもある。
   ここでも、炭素が動く完全な六次元軌道を表す図ではない。

4. [VPDパルスと水ポテンシャル軌道](../artifacts/spac_carbon/rational_e4_diagnostics_20261016/figures_v3/e4_vpd_pulses_and_xylem_trajectories.png)  
   上段は外から与えたVPD、下段は実際の木部・葉水ポテンシャルと、その時刻の凍結全系平衡が持つ \(\psi_x^*\)。
   10日パネルは上昇・下降が時間軸上で急に見えるが、実験上のランプ幅は明示した10日である。
   10年パネルでは実軌道が移動平衡から遅れる様子を見比べられる。
   \(\psi_x^*\) は移動する平衡の座標であり、capacityではない。

5. [葉水ポテンシャル–葉炭素の相軌道](../artifacts/spac_carbon/rational_e4_diagnostics_20261016/figures_v3/e4_leaf_carbon_phase_trajectories.png)・
   [3次元相軌道 \((\psi_l,C_L,C_W)\)](../artifacts/spac_carbon/rational_e4_diagnostics_20261016/figures_v3/e4_three_dimensional_phase_trajectories.png)  
   色はVPD、星は出発点の安定平衡。10年ランプでは水分・葉炭素・木部炭素の変化と履歴が見える。
   これらは実際の六次元軌道の射影で、射影面内のベクトル場や吸引域境界ではない。
   10日ランプの終点にある下向き三角は、積分を止めた \(C_L=10^{-6}\) の数値境界である。

6. [4プールのcapacity \(K_i\) と実ストック \(C_i\) の時間変化](../artifacts/spac_carbon/rational_e4_diagnostics_20261016/figures_v3/e4_classical_capacity_and_carbon_pools.png)・
   [同じ軌道の400年全期間](../artifacts/spac_carbon/rational_e4_diagnostics_20261016/figures_v3/e4_capacity_and_stocks_full_horizon.png)  
   実線は \(C_i\)、破線はその時刻のcapacity \(K_i\)。早期図はゼロ値も表示できる対称対数軸、
   全期間図は対数軸を使う。10日ランプで気孔閉鎖後にcapacityがゼロになる点もここで読める。
   10年ランプの帰還後は \(K_i\) とストックが徐々に再接近する。
   全期間図では、長い回復時間を含む400年の応答を確認できる。
   これらの曲線は capacity が常に状態を引き寄せることの証明ではない。

7. [capacityとストックの差](../artifacts/spac_carbon/rational_e4_diagnostics_20261016/figures_v3/e4_capacity_minus_stock.png)  
   縦軸は固定基準値で割った \((K_i-C_i)/s_i\)。初期平衡では0、VPD変化中はプールごとに異なる符号・遅れを示す。
   これは状態依存capacityと実ストックのずれを直接表示する。和が0かどうかだけで四プールの一致を判定してはいけない。

8. [各プールの \(K_i\) 対 \(C_i\)](../artifacts/spac_carbon/rational_e4_diagnostics_20261016/figures_v3/e4_capacity_vs_current_carbon.png)  
   軸を固定基準値で規格化し、破線 \(K_i=C_i\) を示す。時間色が動的経路の進行を表す。
   平衡では対角線上にある一方、過渡中に離れることを確かめられる。

## この二つのランプ実験で観察されたこと

- 3.00 kPaの凍結全系平衡では、最大固有値実部は約 \(+0.01085\ {\rm day^{-1}}\)。
  これは局所不安定な平衡を指すが、有限時間の強制軌道がその平衡に乗ることは意味しない。
- **10日ランプ**では、VPDを元に戻した後も葉量が減少し、約13.94年で数値葉量下限に達した。
  この時点の \(\psi_l\) は約−4 MPaで気孔が閉じ、瞬間capacityは4プールともゼロ。
  しかし残っているストックは \(C_W\approx0.725\)、\(C_S\approx2.764\ {\rm kg\,C\,m^{-2}}\) など正である。
  よって「\(K=0\)」は「既存ストックが即時ゼロ」ではない。
  ここで積分を終了しており、その後の完全な炭素軌道は主張しない。
- **10年ランプ**では、帰還後400年の時点（実験開始から421年）で積分は完了し、
  炭素ストックは出発平衡の約0.03%以内に戻った。
  これはこの初期値・外力・有限期間の数値結果であり、大域的回復や全初期値での安定性の証明ではない。
- 二例の違いは、選んだパラメータではVPD変化速度に対する応答が大きく異なることを示す。
  二つのランプだけで、一般的なrate-induced tippingや回復／退出境界を確立したとは解釈しない。

## 再現と検証

実行コマンド：

```bash
python examples/analyze_spac_carbon_rational_e4_diagnostics.py \
  --output artifacts/spac_carbon/rational_e4_diagnostics_20261016/NEW_RUN
```

全状態・VPD・capacity・\(K-C\)・凍結平衡枝は各 `*.npz` と
[results.json](../artifacts/spac_carbon/rational_e4_diagnostics_20261016/figures_v3/results.json) に保存した。
両軌道で
\[
\dot C=M(C-K)
\]
の実装残差は \(1.5\times10^{-19}\) 未満であり、凍結全系平衡では \(K=C\) の最大成分差が \(1.8\times10^{-15}\) 未満だった。
計算入力のハッシュも結果JSONに記録してある。
