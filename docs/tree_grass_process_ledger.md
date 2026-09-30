# 木本–草本過程台帳：第一回探索

2026-09-30。[研究計画](tree_grass_carbon_climate_research_plan.md)の段階1。
一次資料の確認、対象選定、未検証の閉鎖仮定を記録する。数値較正は未実施。

## 1. 対象と選定した相互作用

対象候補を**ブラジル中央部Cerradoの自然なサバンナ–森林境界**に絞る。
Tは森林を構成する木本群、Gは開放条件の草本群を代表する。
サバンナ木本を独立に持たないため、二群だけでCerrado全群集を表したとはいえない。
森林伐採でできた湿潤熱帯の境界とは区別する。

選定する作用は、木本葉量増加 → 下層光の減少 → 草本生産・草本由来燃料の減少
→ 火災による木本地上部損失の減少 → 木本葉量の維持、という閉じた相互作用。
これは文献に支えられた**作用の候補**であり、その強さや10プール閉鎖の正しさは未検証。
対照は火災を除いた一方向遮光系とする。水競争・加入抑制を同時には加えない。

## 2. 一次文献と確認範囲

| ID | 文献・確認範囲 | 採用する情報 | 採用しない推論 |
|---|---|---|---|
| E1 | Hoffmann et al. (2005), *Seasonal leaf dynamics across a tree density gradient in a Brazilian savanna*, Oecologia 145, 307–316, [DOI](https://doi.org/10.1007/s00442-005-0129-x), [公的機関PDF](https://www.fs.usda.gov/pnw/pubs/journals/pnw_2005_hoffmann001.pdf)。本文閲覧 | 木本・草本LAIと水分の観測。調査地では木本LAI約3.3でサバンナ草本が排除される | 3.3を普遍的閾値にしない。LAIから炭素への変換には炭素基準SLAが必要。季節相関から温度単独の因果応答を決めない |
| E2 | Hoffmann et al. (2009), *Tree topkill, not mortality, governs the dynamics of savanna–forest boundaries under frequent fire in central Brazil*, Ecology 90, 1326–1337, [DOI](https://doi.org/10.1890/08-0741.1)。出版社抄録確認、本文の係数未抽出 | 地上部枯死、個体死亡、萌芽回復を分離。森林樹種の薄い樹皮とtopkillの関係 | 木本炭素全体を一律に燃焼させない。topkill率を根の死亡率や大気放出率へ置換しない |
| E3 | Hoffmann et al. (2012), *Ecological thresholds at the savanna–forest boundary: how plant traits, resources and fire govern the distribution of tropical biomes*, Ecology Letters 15, 759–768, [DOI](https://doi.org/10.1111/j.1461-0248.2012.01789.x)。出版社抄録確認、統合・総説として使用 | 個体の耐火性と林冠による火災抑制は別の閾値 | 木本/草本の普遍的耐熱順位や年当たり死亡率は供給しない |
| E4 | Thonicke et al. (2010), *The influence of vegetation, fire spread and fire behaviour on biomass burning and trace gas emissions: results from a process-based model*, Biogeosciences 7, 1991–2011, [本文](https://doi.org/10.5194/bg-7-1991-2010)、[訂正](https://doi.org/10.5194/bg-7-2191-2010)。本文§2.2・Appendix Bと訂正を確認 | 燃料量・燃料含水率・着火・燃焼消費・植生損傷を区別する過程構造 | 本研究の単一hazardをSPITFIREの正確な縮約とは呼ばない |
| E5 | [SPITFIRE1.9のモデル検証論文 (2025)](https://doi.org/10.5194/gmd-18-2021-2025)。燃料水分と火災伝播の問題点を確認 | 旧版火災式には実装上・構造上の不確実性がある | 2010年の係数を確認なしに移植しない |

公開CC BYのE4本文・訂正の取得状況とハッシュは
[literature/tree_grass_sources.json](literature/tree_grass_sources.json)に記録。
その他は上記恒久リンクを使用し、閲覧できた範囲を明記した。

## 3. 式・測定・行列への対応

量は地表面積当たり kg C m⁻²、時間は年。
以下の「仮定」は本研究で提案する近似であり、文献の原式ではない。

| 過程 | 式・単位と出典区分 | 必要な測定量・範囲 | 入る場所／数学的役割 | 未検証点 |
|---|---|---|---|---|
| 炭素→LAI | LAI_T = σ_T L_T、σ_T: m²葉/(kg C)。仮定 | 炭素基準SLA、葉C率。数値範囲未決定 | 入力uの状態依存性 | 乾物基準SLAなら炭素質量分率で割る必要がある |
| 上層遮光 | I_G/I_0 = exp(−κ L_T)。κ: (kg C m⁻²)⁻¹。V1の構造を参考にした近似 | 林冠下PPFD、LAI、群落消散係数 | μ_Gに作用。一方向結合 | NPPへ同じ指数を掛けるのはnative GPPの転記ではない |
| 自己遮蔽を含む生産 | μ_T=P_T(1−exp(−a_T L_T))、μ_G=P_G exp(−κ L_T)(1−exp(−a_G L_G))。P: kg C m⁻² yr⁻¹、a: 炭素面密度の逆数。解析用仮定 | 群別NPP、LAI、生育期。P_i(T,θ)は未同定 | Bμ。単位葉量当たり入力が減少する対照 | 呼吸を差し引いた非負NPPの範囲に限定。炭素欠乏を再現しない |
| 草本→火災 | f=λ(T,θ)h(χD_G)、0≤h≤1、h(0)=0。λ: yr⁻¹、χ: 無次元。E1–E4から動機づけた仮定 | 地表細燃料、燃料含水率、着火・焼失面積。hの形・係数未決定 | Mの損失・移動率へ戻りの作用 | D_Gには根由来リターも入る。χを一定とする根拠、燃料連続性、空間平均化が未検証 |
| 木本地上部損傷 | η_ij f C_ij、η∈[0,1]。η f: yr⁻¹ | 樹皮厚、サイズ別topkill、萌芽。E2 | 非燃焼分はMのD行、燃焼分は列の不足 | 原則η_TR=0を地表火災の対照とするが、根死亡ゼロの普遍性は主張しない |
| 燃焼と枯死物 | 損傷炭素のb_ij∈[0,1]を放出、残りをD_iへ。D燃焼=c_i f D_i | 燃焼完全度、残渣。E4は構造参照 | E_fireと内部移動を分離 | 炭化物・立枯木を独立に持たない。Dは統合枯死物という近似 |
| 分解 | k_D D_i、k_H H_i、D→Hはρ k_D D_i、0≤ρ≤1 | 分解・腐植化率、温度/含水率。数値未決定 | M_DD、M_HD、M_HH。Rhは残分 | 由来別の違いの根拠がないため共通率を用いる |

## 4. native VISITの出典境界

V1は `Sachitama2001/VISIT-matrix` の
`3285bd8e131a932e338b59892751648fd9edcc7b`、`visit_local/`。
今回ローカルcheckoutのHEADと以下の遮光代数を再確認した。

- `radiation.c::f_net_rad`（247、285–288、309–313行）:
  I_G = I_top exp(−0.9 k_W LAI_W)。Iは μmol photon m⁻² s⁻¹、
  LAIとk_Wは無次元。入射光式のソース代数として確認。
- `photosynthesis.c::f_gpp`（22–42行）は光・LAIを用いた日GPP。
  native炭素の規約は Mg C ha⁻¹、日フラックスは Mg C ha⁻¹ day⁻¹。
  このNPP仮定と同一ではない。
- 更新順は `experiment.c` のlocation処理→daily処理、
  `location_proc.c::f_loct_proc`で放射・水文・生理を更新し、
  `daily_scheme.c::daily_scheme`（55–57行）で木本→C3→C4、
  `plant_proc.c::plant_process`（134行）からGPPを呼ぶ既存監査を参照。
  同時更新ODEへの厳密変換は行っていない。
- `disturbance.c::disturbance_regime`（31–40行）の指定日伐採と
  `logging_event`は、草本燃料依存hazardの根拠ではない。

元の炭素量を変換する場合は 1 Mg C ha⁻¹ = 0.1 kg C m⁻²。
一定の日率を年率にする場合に限り、1年=365日規約なら365倍。
今回はnative係数を移植していない。VISITcの別スナップショットは使用しない。

## 5. 気象応答の判定と次の調査

外力候補は T（℃）と土壌体積含水率θ（m³ m⁻³）。
P_T(T,θ)、P_G(T,θ)、λ(T,θ)、分解率の**定量関数は未選定**。
気象係数を逆算してR-tippingを作る段階には進まない。
背景の日射・CO₂・土壌・着火・風・湿度は値を含め今後の観測条件に合わせて固定する。

E4は枯死燃料の乾燥履歴と生草の水分を分けている。したがってθを燃料含水率に
そのまま代入することはできない。Tとθだけでλを閉じるには、湿度等を固定した
条件下の経験応答を検証するか、燃料水分の記憶を持つ補助状態が必要である。
SPITFIREのAppendix B2の相対土壌水分も、θそのものとは違う。

**判定**：段階1の対象・相互作用・炭素の行き先は選定済み。
気象応答の較正と燃料代理変数の検証は未完。次の[方程式](tree_grass_carbon_equations.md)と
[解析](tree_grass_bistability_analysis.md)は、その不確実性を明示したモデル非依存の
理論候補・対照であり、source-groundedな完成モデルではない。

問い：戻りの作用を炭素収支に残せるか。仮定：平均火災率・統合枯死物。
得た判断：可能だがD_Gの動態を保持する必要がある。
未検証：燃料含水率、サイズ別損傷、萌芽。
次の一手：E2本文のサイズ別topkill/萌芽と、対象地域の細燃料・含水率資料を照合する。
