# 非線形炭素–水理結合系における Carbon Storage Capacity の拡張
## 研究計画書・実装指針（Stage C：Capacity-first）

作成日：2026-10-13  
対象リポジトリ：`Sachitama2001/Control_carbon_model`  
確認した `main`：`2ef9438c080d8183edc1efeff6e3a3b9fad62eb1`  
位置づけ：**Stage A/B1/B2/B3 E4 の直後に接続する研究計画**。既存結果は保全し、まだ新規数値実験を行ったわけではない。

---

## 0. 要旨：この研究の主題は「心材率」でも「Hopfを見つけること」でもない

**中心的な研究質問**：Luoらの carbon storage capacity（CSC）と carbon storage potential（CSP）の考え方を、入力・移動係数が状態依存し、さらに植物の水理と炭素が双方向に結合する非線形系へ拡張するとき、capacityを何と定義し、どの条件で平衡・吸引先・追従目標として解釈できるか。どの条件でその解釈が破綻するか。

この問いに対する貢献は、複雑な木部生理モデルの発明ではない。既存の水理–炭素低次元モデルを**診断可能な反例・ベンチマーク**として使い、定義・安定性・時間変化・吸引域に関する必要条件と十分条件、または境界例を明確にすることである。

- **主目的**：非線形系のCSCの定義・自己整合性・一意性・局所／大域安定性・tracking・非平衡吸引子を体系化する。
- **応用例**：Manzoniベースの水理2状態と炭素4プールからなるE0–E4モデル。
- **補助目的**：木部feedbackの関数形に研究結論が依存しすぎないかを小規模に感度確認する。
- **今回扱わない主題**：心材形成の詳細メカニズムの開発、樹種較正、森林死亡率の推定、Leuning/FvCBの導入、VISIT全系への直結。

**成果の成功基準**：たとえ真のR-tippingや安定周期軌道を発見できなくても、非線形CSCに対して「どの量をcapacityと呼べるか」「どの命題が成立しないか」「何を確認しなければ追従先とは呼べないか」が論理的・計算的に整理されていれば研究成果とする。

## 1. リポジトリ現在地：何を既に得ていて、何が欠けているか

1. **既存の理論核（2026-09-17）**：`docs/capacity_qse_attractor_exploration.tex`。Luo型の線形恒等式、状態依存capacity写像、固定点、局所Jacobian、CSPがLyapunov関数と限らないこと、複数QSEと速度誘起枝移動の数学的反例を整理済み。新しい計画はこの理論を再利用・修正拡張する。
2. **Stage A（2026-10-09）**：水→炭素の一方向結合。条件付き炭素正平衡の一意性・閾値、凍結入力capacity写像の意味、有限温度パルス後の回復条件。
3. **Stage B1（2026-10-09）**：葉・木部の構造炭素→通水能力／蒸散需要feedback（E0–E4）。全Jacobian、fold・分枝、速度依存の葉量領域退出。
4. **Stage B2（2026-10-12）**：E2/E3の凍結正平衡枝の連続区間安定性、subcritical Hopf、不安定周期軌道を検証。ただし高速変化後の別吸引先は未定義。
5. **単位点検・Stage B3 E4（2026-10-13）**：外向けVPDをMPaへ統一し、25℃の物理的範囲を点検。E2/E3のHopfはその範囲外、E4は約0.002664752905 MPaに数値的subcritical Hopf。小振幅周期解は**不安定**。基本E4の木部指数・損失仮定に感度あり。関連回帰163 passedと報告。

**未完了の中心課題**：既に得られた分岐・周期解・速度実験を、Luo型CSC/CSPに関する定理・診断図・反例の言葉で再表現すること。E4の周期解が不安定な事実から「別吸引先」やR-tippingを導いてはならない。

既存コード・レポートは変更せず、新規Stage Cを独立シリーズとして実装する。

## 2. 数学的枠組みと用語の厳密な区別

### 2.1 Luo型の線形基準

炭素状態 `C∈R_+^n`、外生入力 `u`、移送・損失行列 `M` を用いて

\[
\dot C=u+MC,\quad C_{\mathrm{cap}}=-M^{-1}u,\quad
C_{\mathrm{pot}}=-M^{-1}\dot C=C_{\mathrm{cap}}-C.
\]

ここで `M` は可逆とする。固定外力・一定係数かつ `M` がHurwitzなら、`C_cap` は一意な大域吸引平衡である。時間依存外力なら `C_cap(t)` は瞬間的な基準点であり、一般には実際の軌道ではない。

Luoらの原著を、状態依存を一律に禁じているものとして扱わない。**本研究の新規部分は、非線形系でのcapacityの複数の意味と成立条件を明示的に分類・検証する点**である。

### 2.2 非線形炭素方程式での正確な恒等式

\[
\dot C=u(C,z,v)+M(C,z,v)C,\qquad \dot z=H(z,C,v),
\]

ここで `z=(ψ_x,ψ_l)`、`v` は外向けVPD [MPa]。`M` が可逆な領域で

\[
\boxed{K(C,z,v)=-M(C,z,v)^{-1}u(C,z,v)}
\]

を**点ごとの凍結入力capacity写像**と呼ぶ。すると

\[
\dot C=M(C,z,v)[C-K(C,z,v)],\qquad
P(C,z,v)=K(C,z,v)-C=-M^{-1}\dot C.
\]

この表示は恒等式であり、`K` 自体が平衡・長期吸引先・軌道であることを含意しない。`M` を転置したり、`P` の符号を取り違えたりしない。Luoとの呼称の対応を文書冒頭で固定する。

### 2.3 **少なくとも5種類の対象を区別する**

| 対象 | 定義 | 意味・注意 |
|---|---|---|
| (I) 瞬間capacity写像 `K(C,z,v)` | 現在の状態で入力・係数を凍結した `-M^{-1}u` | 非線形系では一般に平衡ではない |
| (II) 水理条件付き炭素平衡 `C*(z,v)` | `C=K(C,z,v)` を解く | 水理状態 `z` を固定する条件付き根。複数あり得る |
| (III) 全結合系の凍結平衡 `y*(v)` | `C=K` かつ `H=0` | 安定性は全Jacobianで判定。複数枝を許す |
| (IV) 凍結系の吸引集合 `A(v)` | 平衡、安定周期軌道など | 単一のcapacity点に還元できるとは限らない |
| (V) 非自律系の実軌道 `y(t)` | `v=v(t)` による解 | 凍結系の枝を追跡する保証は別問題 |

用語方針： (I) を**capacity写像**、(III) を**自己整合的な全系平衡**、安定枝を**吸引平衡枝**と呼び分ける。安定周期軌道が得られた場合は「非平衡吸引集合」であり「新しい点capacity」とは呼ばない。対象が一意でない場合は**集合値のcapacity候補**を必要に応じて導入するが、Luoの元の定義を無言で置き換えない。

### 2.4 局所安定性：capacity写像の傾きが何を変えるか

炭素単独の自己整合平衡 `C*` のJacobianは、十分な滑らかさの下で

\[
J_{CC,*}=M_*[I-D_C K|_*].
\]

`M_*` がHurwitzであるだけでは `J_CC` は安定と限らない。さらに全結合系のJacobianは

\[
J_* = \begin{pmatrix}J_{HH}&J_{HC}\\J_{CH}&J_{CC}\end{pmatrix},
\]

であり、`J_CC` の安定性だけで `J_*` の安定性を判断してはならない。E4の木部feedbackはまさにこの違いを検査する材料となる。`M` が状態依存でも、平衡で残差 `C-K=0` のため上の `J_CC` 表示が成立する。スカラーの `K'<1` 条件を一般行列に無条件で転用しない。

### 2.5 carbon storage potentialは「距離」か？

`P=K-C=-M^{-1}\dot C` はLuo型のdisequilibriumの代数表示だが、非線形の局所安定性・大域収束を常に測るLyapunov関数ではない。以下を別々に調べる：

- `||P||` または総和 `1ᵀP` は単調減少するか。
- `1ᵀP=0` でも `P≠0` となる反例はあるか（プール間相殺）。
- `P=0` でも `H≠0` となり、水理系が動いている例はあるか。
- 全結合系の平衡や周期軌道に対して、より適切な追従誤差は何か。

**推奨する明確な診断**：`||P||`、`||C-C*(z,v)||`、`||y-y*(v)||`、`||F(y,v)||`、全Jacobianのスペクトル、吸引域ラベルを同時記録する。単位の異なる量は尺度行列を明示して正規化する。総炭素量だけで平衡判定をしない。

## 3. 研究仮説と検証可能な命題

**H1（定義）**：非線形系でも `K=-M^{-1}u` という点ごとの恒等式は有効だが、Kを一意の貯留目標とは一般に解釈できない。固定点の数と全系安定性を別に調べる必要がある。既存の一般式・反例から出発して、E0–E4へ具体化する。

**H2（局所性）**：安定な平衡枝の存在とその局所安定性は、環境変化に対する追従を保証しない。安定性余裕、枝の移動速度、非正規増幅、吸引域からの距離等を区別して評価する。

**H3（非平衡集合）**：Hopfが存在しても、安定周期軌道が見つかるとは限らない。基本E4の既検証の局所周期解は不安定であり、「capacity点ではない新たな吸引先」の実例ではない。安定周期軌道が得られない場合でも、それを否定結果として記録する。

**H4（モデル頑健性）**：H1/H2の一般的な定義・命題と、E4のHopf数値位置・局所型・速度依存性を区別する。木部feedback関数形・炭素速度・水理容量を限定して変更し、どれが構造的／パラメータ依存かを検査する。

**H5（研究の適用限界）**：`L_min` 退出は別吸引状態ではない。モデルに回復・死亡レジームがない現状では真のR-tippingを結論にしない。研究主目的を達成するために、死亡機構を無理に追加する必要はない。

## 4. 次段階の作業計画：capacityの問いから順序を決める

### WP-C0：現状の固定と、capacity関連コード・図の棚卸し【必須・最初】

- `main` SHA `2ef9438...` を基準に固定し、Stage A/B1/B2/B3 E4の再現コマンドと試験結果を再確認する。過去artifact（git-ignored）をGitHub上に存在するものと決めつけない。
- `docs/capacity_qse_attractor_exploration.tex` の命題・符号規約と `docs/spac_carbon_research_logic_report.md` のcapacity区分を照合。
- `src/control_carbon/spac_carbon.py`、`spac_carbon_two_way.py`、`spac_carbon_mpa.py`、`spac_carbon_e4.py` で `u,M,K,P` の実装位置と単位を台帳化。
- **成果物**：`docs/spac_carbon_capacity_definition_ledger.md`。式、定義域、必要仮定、反例、証明済み／数値確認／未検証を明示。

### WP-C1：capacity診断器を既存E0–E4に付ける【第一優先】

- 既存微分方程式を変更せず、`K(C,z,v)` と `P=K-C`、条件付き炭素根、全結合平衡、全系残差を計算する純関数を新設。
- `P=0` と `F(y,v)=0` の違い、`P` の総和とベクトルの違いを、再現可能な例で表示。
- 各モデルで凍結 `v` の安定枝・鞍点・形式的ゼロ炭素根を適切に区別。`L=0` の生理的意味を仮定で補わない。
- **受入条件**：固定点恒等式、水・炭素収支、E0→既存Stage Aへの一致、解析Jacobianとの差分一致、正領域維持。
- **成果物**：`src/control_carbon/spac_carbon_capacity.py`、`tests/test_spac_carbon_capacity.py`、capacity診断CSV。

### WP-C2：追従問題を「capacity写像 vs 移動平衡枝」で再定式化【第二優先】

- 凍結平衡枝 `y*(v)` を選び、実軌道との追従誤差 `e(t)=y(t)-y*(v(t))` を定義する。
- 正確な誤差式：

\[
\dot e = F(y^*(v)+e,v)-F(y^*(v),v)-\frac{dy^*}{dv}\dot v.
\]

  線形化すれば `edot ≈ J_*(v)e - y*'(v) vdot`。右辺第2項は**動く平衡点を追うための強制項**。固有値の実部が負でも `v'(t)` が大きければ誤差が大きくなり得る。
- `dy*/dv` は平衡条件 `F(y*,v)=0` から、`J_*` 可逆なら `dy*/dv=-J_*^{-1}F_v` と導出する。fold近傍では特異性に注意。
- E2/E3の既存「10日 vs 10年」はVPD 0.00303975 MPaまでの実験として正確に再掲。E4で同じ結論を仮定しない。E4は新規の限定的速度比較を追加する場合に、主解析範囲 `v≤0.003 MPa` と、旧経路の上限超過を区別する。
- `||P||`、全状態誤差、`K` の変動速度、固有値・非正規性、応答時間尺度を時系列で比較。誤差が小さい条件／大きくなる条件の**十分条件候補**を分析する。
- **成果物**：追従誤差の導出ノート、平衡枝・capacity・実軌道を分離した図、速度比較の再現結果。R-tippingを目的化せずtrackingの範囲を明確にする。

### WP-C3：非線形性によりcapacityの解釈が変わる事例を整理【第三優先】

- **基準1（線形）**：固定外力で一意capacity・大域吸引。時間変化で位相遅れ。解析的に説明できる小モデルを使用。
- **基準2（炭素非線形）**：既存 `capacity_feedback.py` の正炭素反例で複数自己整合capacityと枝間の速度遷移を示す。これは森林の実証ではなく数学的反例。
- **基準3（E0–E3）**：水理feedback導入によって炭素方向の`K`と全系安定性が分離すること、foldや領域退出がどのcapacity概念に対応するかを検査。
- **基準4（E4）**：Hopf前後で平衡・`K`・周期解の関係を計測。**局所不安定周期軌道は吸引先ではない**ことを明示。安定周期吸引子・他の平衡の存在を探索する場合も、見つからなければ「未発見／指定範囲で否定」と記録。
- 既存の一意性定理、縮小写像条件、Lyapunov／比較原理、入力–状態安定性に基づく追従評価を適用できる範囲から探す。厳密に導出できない部分は数値診断として区別。
- **成果物**：「capacity解釈の成立条件」マトリクス、成立例・反例・未解決例を並べた図表。

### WP-C4：構造feedbackの頑健性――生理の細部は従属課題【低コストで実施】

- E4の `F_W=(W/W_ref)^β_W` は**合成の構造仮説**。既存の感度結果（β=0.75/1/1.25、木部損失、容量）をまず利用する。
- 追加検査は少数に制限し、`F_W` の弾性 `ε_W=dln F_W/dln W` を共通指標として、べき乗型・飽和型など**単調性と正則性を保つ1–2種類**を比較する。すべてのモデルで `F_W(W_ref)=1` を揃え、基準平衡と局所微分が一致する／しない場合を区別。
- 必要なら文献の `G~K_s A_s/ℓ` を**係数・次数の妥当性を点検する物理的制約**として使う。ただし心材率を新しい動的状態にすること、詳細な辺材成長モデルの開発は行わない。
- 評価対象は「Hopfを必ず維持するか」ではなく、**capacityの定義・trackingに関する結論がどれほど普遍的か**。Hopfが消えても研究失敗としない。
- **成果物**：少数シナリオ表（平衡の数、安定性、capacity診断、追従誤差、感度・適用限界）。

### WP-C5：最終的な論文命題と停止点【主成果】

最低限、以下を厳密に結論分類する。

1. Luo型の点ごとのcapacity恒等式は、どの仮定の下で非線形系にも成立するか。
2. capacity写像の固定点と、条件付き炭素平衡・結合系全平衡の違い。
3. `M` の安定性・炭素サブ系安定性・全結合系安定性・非自律追従性の非同値性。
4. `P=K-C` が単一の「吸引先までの距離」や大域Lyapunov関数ではない具体例。
5. 複数平衡、周期軌道、動く吸引域があるときのcapacity報告単位・図示・診断手順。
6. どの結果が定理、計算機援用保証、数値実験、未検証仮説なのか。

**停止条件**：全ての証拠が「新しいtippingや死亡を無理に作らず」に論文として整理可能になったら第一稿へ進む。生理的死亡状態や高度な水理allometryの追加は必須条件にしない。

## 5. 優先して導く定理／反例の候補

### 命題1（恒等式と固定点）
`M` 可逆の非線形炭素方程式でcapacity写像の恒等式、および `C=K` と `dot C=0` の同値性を明記する。ただし `dot z=0` は自動ではない。証明は代数的。

### 命題2（平衡でのJacobian分解）
条件付き炭素平衡では `J_CC=M(I-D_C K)`。水理との双方向結合では全JacobianのHurwitz性が必要。Stage Aの三角性を外した際の違いを示す。

### 命題3（追従誤差の強制項）
安定平衡枝の近傍で `e'=J_*(v)e-y*'(v)v'+O(||e||²)`。十分な一様安定性・平衡枝の滑らかさ・局所領域条件が得られる範囲で `|v'|` に依存する誤差評価を試みる。点ごとの固有値が左半平面にあるだけでは、時変線形系の一様安定性や大域追従を自動的に保証しないことに注意。

### 命題4（potentialの限界）
既存の非線形スカラー反例でpotentialの絶対値が増加する条件を示し、E4でも `||P||` と平衡距離が乖離する具体的な軌道を探索する。見つからなければ既存反例とE4の検証範囲を分けて記録する。

### 命題5（非平衡吸引集合への拡張：条件付き）
周期軌道が**安定な吸引集合**として確認できたときに限り、凍結外力下の長期炭素貯留を時平均・振幅・位相など集合値／統計値で報告する定義候補を比較する。E4で確認済みの不安定周期解をこの実例に使わない。

## 6. 事前に決める評価指標と図

| 診断 | 指標・比較 | 誤認を避けるポイント |
|---|---|---|
| capacity identity | `M(C-K)-dot C` の残差 | 恒等式の成立は吸引性を意味しない |
| 炭素disequilibrium | `P=K-C`、`1ᵀP`、`||P||` | 総和ゼロとベクトルゼロを区別 |
| 全系平衡 | `||F(y,v)||`、`||y-y*(v)||` | `P=0`でも水理が非平衡の可能性 |
| 全系安定性 | 全Jacobian固有値・Hurwitz余裕 | `M` の固有値だけでは不十分 |
| 追従性能 | 平衡枝の速度 `|y*'(v) v'|`、軌道誤差 | 追従失敗 ≠ 別吸引先への遷移 |
| 吸引域 | 多初期値・限定断面・領域退出 | 有限断面から大域境界を断言しない |
| 長期対象 | 正平衡・退出・（もし確認できれば）安定周期解 | 不安定周期解は吸引子でない |
| 物理的感度 | `F_W`弾性・水理容量・炭素損失 | 合成パラメータの閾値を観測値と呼ばない |

基本の図：①線形capacityと非線形K写像、②`K`・`C*(z,v)`・全平衡・実Cの同時図、③potentialと全状態平衡距離、④凍結分枝／吸引域／速度ramp、⑤E4の局所Hopfとcapacity診断、⑥構造関数感度。図の軸と単位、安定・不安定の線種、数値的根拠の種別を明記。

## 7. 実装上の制約・受入基準

- **既存Stage A/B1/B2/B3の方程式とコードを破壊的変更しない**。新規Pythonモジュール・テスト・runner・文書で構成。
- 外向けVPDは `v [MPa]`、内部旧コードの `D [mol/mol]` との変換 `D=v/p_atm` を一箇所で一貫して管理。`p_atm=0.101325 MPa` を合成設定として明示。
- 基本解析範囲 `0≤v≤0.003 MPa`。旧B1速度経路の終点 `0.00303975 MPa` は上限を少し超えるため、過去結果の再現と新規主解析を混同しない。温度25℃は湿度整合性診断に限り、追加の温度生理モデルを含めない。
- E0–E3で能動3状態、E4で能動4状態。受動プールの分離はモデルの依存関係に基づき利用する。
- 非負性、水・炭素収支、`M` 可逆性と条件数、Jacobian差分、解の残差、積分精度、平衡枝の根の漏れ、時間刻み・許容誤差の収束検査を必須。
- 条件付き吸引根が消える／foldで`J`が特異になる領域では、滑らかな枝を無理に継続しない。結果と失敗条件を明記。
- `L_min` は領域退出イベントであり死亡閾値ではない。人工的なclip、残留conductance、回復入力の追加によって欲しい結果を作らない。
- 不安定周期軌道、Hopf、単なる大振幅過渡、別吸引先へのR-tippingを厳密に区別。
- 証拠のラベルを統一：【解析的証明】【計算機援用区間保証】【精度確認済み数値】【仮説・未実施】。
- 自動生成artifactは新規パスへ保存し、Git無視下の旧画像等をGitHubに存在するとみなさない。JSONに基準SHA、入力・パラメータ・ソースhash・solver設定を記録。

想定ファイル：

```
docs/spac_carbon_capacity_definition_ledger.md
src/control_carbon/spac_carbon_capacity.py
tests/test_spac_carbon_capacity.py
examples/analyze_spac_carbon_capacity.py
docs/spac_carbon_stage_c_results.md
artifacts/spac_carbon/stage_c_capacity_<timestamp>/
```

## 8. コーディングAIにそのまま渡す最初の指示

> **研究目的はLuo型carbon storage capacityの非線形系への拡張であり、新しい植物生理機構やtippingを作ることではありません。** `main` の `2ef9438c080d8183edc1efeff6e3a3b9fad62eb1` を基準に、まず `docs/capacity_qse_attractor_exploration.tex`、`docs/spac_carbon_research_logic_report.md`、`docs/spac_carbon_units_and_scales.md`、`docs/spac_carbon_stage_b3_e4_results.md` と対応実装を読んでください。Stage CのWP-C0とWP-C1を先に実施し、Luoの `capacity`／`potential`、非線形の点ごとのcapacity写像、条件付き炭素平衡、全結合系平衡、吸引集合、非自律追従軌道を区別する定義台帳とcapacity診断器を新設してください。既存実装・方程式・設定を破壊的に変更しないこと。`P=K-C=-M^{-1}dot C`の恒等式、`P=0`と全系平衡の非同値性、符号・単位・保存則・Jacobianをテストしてください。結果は証明／区間保証／数値／未検証に分類し、`docs/spac_carbon_stage_c_results.md` に日本語で報告してください。WP-C2以降の大型実験は、C0/C1の結果と研究目的への寄与を明示してから行ってください。心材率や詳細allometry、Leuning、死亡レジームは導入しないでください。

## 9. 論文の章立て案と期待される結論

1. **Introduction**：Luo型CSCの有用性と、状態依存・双方向feedbackで生じる概念上の課題。
2. **Methods / Definitions**：線形基準、非線形K写像、条件付き固定点、全系平衡・吸引集合・potential、trackingの定義。
3. **Analytical results**：恒等式、Jacobian、potentialの非Lyapunov性、追従誤差式と成立条件。
4. **Case study**：Manzoni–炭素E0–E4。fold/Hopfはcapacity概念の限界を検査する手段として提示。
5. **Robustness and limits**：木部feedbackの限定感度、VPD整合性、葉ゼロ領域の未定義性。
6. **Discussion**：「非線形モデルでcapacityを報告する前に確認すべき事項」のチェックリスト。Luo理論への追加価値と、生理・実測上の限界。

**想定できる中心結論（現時点では研究仮説を含む）**：

> 非線形炭素循環系では、Luo型capacityは状態依存の代数写像として定義できる一方、その値を一意の安定吸引先とみなすには、自己整合性、全結合系の安定性、吸引域、および時間変化する外力への追従性を別々に検証する必要がある。非線形feedbackと時間スケール差はこれらを非同値にし得る。

**重要**：上記の一般的な考え方の一部はリポジトリで既に示されており、Luoら自身も非線形／動的植生の問題を全く否定していない。論文の独自性は「非線形capacityという言葉を初めて使う」ことではなく、**定義・成立条件・反例・双方向結合ベンチマーク・再現診断を統合する**ことに置く。

## 10. 出典・位置づけ

### 実装・結果に直接依拠する内部資料（2026-10-13時点）
- [`docs/capacity_qse_attractor_exploration.tex`](https://github.com/Sachitama2001/Control_carbon_model/blob/2ef9438c080d8183edc1efeff6e3a3b9fad62eb1/docs/capacity_qse_attractor_exploration.tex) — 非線形capacityの数学的出発点と速度誘起枝遷移の反例。
- [`docs/spac_carbon_stage_a.md`](https://github.com/Sachitama2001/Control_carbon_model/blob/2ef9438c080d8183edc1efeff6e3a3b9fad62eb1/docs/spac_carbon_stage_a.md) — 一方向結合・炭素平衡・capacity。
- [`docs/spac_carbon_stage_b_results.md`](https://github.com/Sachitama2001/Control_carbon_model/blob/2ef9438c080d8183edc1efeff6e3a3b9fad62eb1/docs/spac_carbon_stage_b_results.md) — 双方向feedbackと速度実験。
- [`docs/spac_carbon_stage_b2_results.md`](https://github.com/Sachitama2001/Control_carbon_model/blob/2ef9438c080d8183edc1efeff6e3a3b9fad62eb1/docs/spac_carbon_stage_b2_results.md) — 安定枝・Hopf・周期軌道。
- [`docs/spac_carbon_units_and_scales.md`](https://github.com/Sachitama2001/Control_carbon_model/blob/2ef9438c080d8183edc1efeff6e3a3b9fad62eb1/docs/spac_carbon_units_and_scales.md) — VPD単位と気象範囲。
- [`docs/spac_carbon_stage_b3_e4_results.md`](https://github.com/Sachitama2001/Control_carbon_model/blob/2ef9438c080d8183edc1efeff6e3a3b9fad62eb1/docs/spac_carbon_stage_b3_e4_results.md) — E4の局所非線形解析と感度。
- [`docs/spac_carbon_research_logic_report.md`](https://github.com/Sachitama2001/Control_carbon_model/blob/2ef9438c080d8183edc1efeff6e3a3b9fad62eb1/docs/spac_carbon_research_logic_report.md) — 証明・数値・仮説の境界。

### 文献の理論的出発点
- Luo et al. (2022), *Matrix Approach to Land Carbon Cycle Modeling*, JAMES. DOI: https://doi.org/10.1029/2022MS003008 。
- Luo et al. (2017), *Transient dynamics of terrestrial carbon storage: mathematical foundation and its applications*, Biogeosciences. DOI: https://doi.org/10.5194/bg-14-145-2017 。
- Manzoni, Katul & Porporato (2014), *Water Resources Research*. DOI: https://doi.org/10.1002/2013WR015236 。
- Ashwin et al. (2012), *Tipping points in open systems: bifurcation, noise-induced and rate-dependent examples in the climate system*. DOI: https://doi.org/10.1098/rsta.2011.0306 。

補助の木部文献は、conductanceの正値性・単位・可能な関数形を制約する用途に限定し、研究の主題に昇格させない。
