# 進捗総括と次期作業計画

> 2026-09-14更新: 現在の優先計画は
> [炭素-水結合の研究計画](carbon_water_research_plan.md)である。本書のP0-P5と
> [決定論的tipping計画](deterministic_tipping_plan.md)は、完了済み原典検証・
> 数学ベンチマーク・後続解析として残す。

## 0. 研究方向の更新と新規実装

主題を「水貯留・SPAC輸送を明示した行列型炭素モデルで、炭素貯留容量、
非平衡性、応答時間を診断する」へ更新した。R-tippingは研究成立の前提ではなく、
通常応答と追従対象を定義した後に条件が整えば検討する。

今回追加した成果は次のとおりである。

- `visit-manager/VISITc`の現行参照commitを旧VISIT-matrixと別に固定した。
- `f_hydrology`の状態更新を、Penman-Monteith potentialを有効入力として
  原典順に部分転記した。
- baseflowが二重に控除される形の水収支残差を、修正せず診断量として露出した。
- leaf/stem/root/soilそれぞれにCとWを持つ8状態ODEを実装した。
- 水の接続行列、炭素貯留容量、Wei型の容量変化分解、結合平衡、Jacobian block、
  Schur補完を実装した。
- 乾燥・再湿潤の理想実験と、収支・数値精密化を含む再現manifestを追加した。
- PDFと文献URLを`literature_index.md`へ整理した。

## 1. プロジェクトの目標

このプロジェクトの最終目標は、陸域生態系炭素循環モデルを制御理論の
共通言語で比較することである。VISITは最初の検証対象であり、最終的には
モデルごとの内部構造が異なっていても、次の量を比較できる枠組みを作る。

- 炭素プールの状態とフラックス構造
- 固定点および準静的平衡（QSE）
- インパルス応答（IRF）
- 極、モード、時定数、留数
- 気象強制からRh、NEP、総炭素への伝達特性
- 強制速度と炭素吸収源・放出源の境界
- 年周期を考慮した周期LTV応答

原則は、VISIT C原典への忠実性を先に確立し、その後に行列化、線形化、
IRF解析を行うことである。

## 2. 現在までに完了した作業

### P0: 実装成果と出典の固定

出典と検証条件の記録は完了している。作業ツリーの変更をGit commitとして
固定する作業は未完了であり、以下はcommit済みという意味ではない。

- VISIT原典をGit commit
  `3285bd8e131a932e338b59892751648fd9edcc7b`へ固定した。
- 各Python実装が参照するCファイル、関数、変数、近似レベルを記録した。
- 追加されたland-cover別パラメータブック
  `parameter_VISITc_16.xlsx`をSHA-256
  `f8748b9dca3a7e7e38a2aa7ce93fd54c22653948f2467fd2ed7c1f6181edbad3`
  で固定した。
- ブックのMOD12 class 0-15と、Tree、Herb3、Herb4、Soil、Canopyの
  構造を確認した。
- TKYにはclass 4「Deciduous Broadleaf Forest」の土壌値を採用した。
- テスト、再現コマンド、既知の原典不具合を文書化した。

注意点として、Excelはpinned Git commit外の補助資料である。またCが参照する
`parameter_S1b.txt`はスナップショットに存在しない。そのため、Excel値は
land-cover別デフォルトの参照値であり、完全なVISIT実行で実際に読まれた値とは
まだ主張しない。

### P1: ERA5からVISIT気象入力まで

完了している。取得経路は公式Copernicus Climate Data Store（CDS）だけを使う。

- 任意地点を最寄りのERA5 0.25度格子へ丸める。
- CDS要求を月別に分割する。
- 瞬時値と積算値に分かれたNetCDFをZIPから読み込む。
- 温度、降水、放射、気圧、雲量、風、土壌温湿度を共通時間別形式へ変換する。
- 欠落時刻を許さず、UTCの完全な24時間を一日へ集約する。
- VISIT `WMODE=1`の意図されたNCEP/NCAR 19列形式を出力する。
- CDS ZIPを期間固有名でキャッシュし、破損検査と原子的置換を行う。
- 入力ZIPとVISIT出力のSHA-256をmanifestへ保存する。
- 複数月結合と365/366日の完全暦年をテストする。

高山付近の要求座標`(36.146, 137.423)`はERA5格子
`(36.25, 137.5)`へ解決され、公式CDSから1日24時間を実取得して、
VISIT 1行19列まで変換できた。

未修正VISIT全体へ入力する場合は、`f_loct_init`が1980-2009年を直接参照する
ため、最低でもその完全期間が必要である。

### P2: VISIT原典CとPython土壌モデルの比較

固定環境・ベースライン`EX_DECTMP=0`の最小土壌炭素ブリッジは完了している。

- 原典`soil_proc.c::f_cycle_soil`を直接コンパイルする。
- 原典`decomposition.c::frl/frh`を直接コンパイルする。
- checkoutのHEADと、コンパイル対象ファイルがpinned commitから未改変である
  ことをビルド時に検査する。
- 初期9プール、日別6 litter入力、albedoを除く土壌16パラメータ、分解環境をCへ渡す。
- 更新後9プール、Rh、`f_tm_l`、`f_tm_h`をCから取得する。
- C原典とPython直接式を一日および40日軌道で比較する。
- C原典のperturbation-minus-baseline応答を100日IRFと比較する。

現在のテストでは、TKYに対応するMOD12 class 4の土壌デフォルトを使用している。

100日litter pulse実験の比較対象は次である。

```text
native変化量 = C原典の摂動実行 - C原典のbaseline実行
IRF変化量    = Python離散状態空間系による摂動予測
誤差         = native変化量 - IRF変化量
```

得られた最大誤差は次のとおりである。

```text
9プール状態差の最大誤差: 8.03e-11 Mg C ha-1
Rh差の最大誤差:          2.94e-15 Mg C ha-1 day-1
```

これはモデル誤差ではなく、約`1e3 Mg C ha-1`のbaselineと摂動軌道を別々に
浮動小数点計算して差し引く際の丸め誤差の水準である。固定係数・負値
クリッピングなしの範囲では、C原典の炭素部分とPython行列系は同じ力学を
表している。

litter pulse当日のRh変化は0で、翌日は
`3.7026e-4 Mg C ha-1 day-1`となった。これはVISITが当日の分解とRhを計算した
後に、その日のlitter入力を追加する更新順序を確認する結果である。

## 3. 現在のテスト基準

2026-09-14の全suiteは101件を収集し、`88 passed, 13 skipped`である。skipは
別checkoutを必要とする旧VISIT-matrixのnative soil/plant比較だけであり、新規の
炭素-水・current VISITcテストは全件成功している。

P0-P5時点の基準は62件である。今回の決定論的非線形基盤とB/Rベンチマークで
20件を追加した。追加テストは解析解、区画収支、平衡枝、安定性、分岐点、
最終流域、保持時間、精密化による判定の再現性を対象とする。

```bash
pytest -q
```

テストは次の範囲を覆う。

- 汎用連続・離散状態空間
- provenanceとsource commit
- `frl/frh`の温度・水分応答と`stype`未定義分岐
- 9プール土壌の一日収支、行列等価性、固定点、IRF
- CDS取得要求、ZIP復号、単位変換、月結合、完全暦年、manifest
- native C対Python直接計算
- native C摂動差分対IRF
- 日変動環境、温度pulse接線予測、非滑らか閾値拒否
- grouped modal residueによるIRF/DC gain再構成
- implicit QSE環境感度と固定点再求解の有限差分比較
- plant turnover、respiration、allocationの手計算とnative C照合

ネットワークを使うCDS実取得は通常のpytestには含めない。通常テストでは決定的な
合成NetCDFを使い、CDSサービスの待ち時間や障害と変換ロジックの退行を分離する。

## 4. まだ完了していないこと

### 土壌ブリッジの範囲

- 日ごとの環境値をCへ渡す機能と温度pulse検証は完了した。
- 環境値を外部指定した土壌単独実験であり、水文系との動的連成は含まない。
- `EX_DECTMP=2/3`の`frh(stype=1/2)`は原典で未定義なので実行しない。
- C-N連成は比較対象外で、炭素更新後の窒素関数をno-opとしてリンクしている。
- 安定同位体はpinned設定の`SCI_SCHEME=0`で無効である。

### VISIT全体実行

- 気象形式判定のCEAMIP `strcmp`に既知の欠落がある。
- `initialize.c`はregion modeでもsite initializerを無条件に呼ぶ。
- 完全な実行時入力一式と`parameter_S1b.txt`がない。
- spin-upと通常実験でうるう日の扱いが一致しない。
- 植物、光合成、水文、フェノロジーを含むnative実行との比較は未完了である。

### 制御理論解析

- 土壌モードの留数・participation factorと凍結固定点の環境感度は実装済み。
- 環境Jacobianと日変動軌道に沿う接線予測は実装済み。
- reduced 18-pool植物・土壌系、NEP、周期LTVは未実装である。

## 5. P3-P5の実装経過と残課題

### P3: 日変動環境を使うnative土壌比較

P3-1とP3-2は完了した。固定係数だけでなく、係数が日々変わる場合と小さな温度摂動を
検証した。

#### P3-1. ブリッジ入力の拡張

各日について次をCへ渡せるようにする。

- `tmp10_soil`
- `tmp200_soil`
- `soilwtr_l`
- `soilwtr_h`
- `soilappr_l`
- `soilappr_w`
- `fieldcap30`
- `fieldcap`

完了条件:

- 固定環境を毎日繰り返した入力が、現在のP2結果と一致する。
- 日変動環境に対するCとPythonの9プール・Rhが丸め誤差で一致する。

結果:

```text
120日の日変動環境における最大状態差: 1.42e-13 Mg C ha-1
最大Rh差:                              5.55e-17 Mg C ha-1 day-1
```

#### P3-2. 温度・水分摂動

baseline環境に小さな温度または水分pulse/stepを与える。

比較するもの:

```text
C原典の摂動差分
Python非線形直接計算の摂動差分
局所Jacobian/IRFによる一次予測
```

完了条件:

- CとPython非線形直接計算が一致する。
- 摂動幅を半分にしたとき、一次近似誤差が滑らかな領域で概ね二次に縮小する。
- `-20 C`閾値や水分・apertureの`min`切替近傍を別扱いにする。

結果: 2 K、1 K、0.5 Kの上層土壌温度pulseに対する最大状態誤差は
`4.09e-4`、`1.02e-4`、`2.54e-5 Mg C ha-1`で、振幅半減ごとにほぼ1/4となった。
`-20 C`閾値を中央差分がまたぐ場合は線形化を拒否する。

#### P3-3. ERA5との接続範囲を決める

ERA5土壌温度はP3の温度入力候補にできる。一方、ERA5体積土壌水分とVISITの
`soilwtr_l/h`は単位も状態定義も同一ではない。変換式を決める前に、VISITの
field capacity、層厚、水収支との対応を原典から整理する。

完了条件:

- 温度はERA5からどのVISIT変数へ入れるか明記する。
- 水分は「外部強制」「初期値」「native水文との比較」のどれに使うか決める。
- 未検証の直接代入を行わない。

結論: ERA5の層平均土壌温度は`tmp10_soil/tmp200_soil`へ接続する。ERA5の
体積土壌水分は、VISIT内部bucketの`sw30/sww`（mm水深）と状態定義・層境界が
異なるため、現段階では診断値として保持し、`soilwtr_l/h`へ直接代入しない。

### P4: 土壌系のモード・QSE解析

基礎実装と数値検証は完了した。

- 離散極から日単位のモード時定数を計算する。
- 各プールのparticipation factorを求める。
- litter入力からRhへのmodal residueを求める。
- 分解スカラーおよび温度への固定点感度を計算する。
- 解析値を有限差分と比較する。
- forcing、Rh、9プール、total litter-total humus射影、誤差を含む図を生成する。

完了条件:

- 固有値、時定数、留数の定義と単位が明示される。
- QSE感度が有限差分の再計算と一致する。
- native、Python非線形、IRFを同じ図と同じ時間規約で比較できる。

現在のTKY fixtureでは、約3.70、5.42、13.56、29.45、66.27、212.07年の
モード群を得た。重複するlitter極は個別固有ベクトルではなくスペクトル射影で
まとめ、IRFとDC gainを留数から再構成できることを確認した。

固定点環境感度は、implicit differentiationと環境を変えて再求解した有限差分が
一致した。温度上昇は平衡炭素量を減らすが、一定litter入力・他の炭素輸出なしの
平衡Rhは入力総量に等しいため、平衡Rhの温度感度は数値ゼロとなる。詳細は
`visit_soil_modal_qse.md`を参照する。

### P5: 植物構造炭素とreduced 18-pool系

turnover、respiration、allocationの個別sliceを実装し、native Cと照合した。

実装済み:

- TKY class 4の`lf0/lc0/lr0/dcd`を使うtree turnover;
- 通常turnover;
- season 3の落葉とLAI低下後の最終落葉;
- crop stage 5の葉・幹・根70% turnover;
- 深層水ポテンシャル閾値による追加葉落下;
- litterfall直後の`fol/stm/rot`炭素減算とsource provenance。
- `f_q10_ar`による日別Q10順化;
- `f_spcfc_resp`によるsapwood/heartwood・fine/coarse root比呼吸;
- 葉・幹・根のmaintenance respirationとgrowth respiration;
- 正負EPP、LAI、休眠season、crop grainを含む`f_allocation`分岐;
- turnover 4分岐、呼吸、allocation 3分岐のnative C直接照合。

これらは`plant_process`内の個別sliceであり、まだ一つの植物一日更新としては
結合していない。GPP、NSCへの蓄積、leaf emergence、survival reallocationも未統合である。

1. GPP/EPPを既知の有効入力としてtreeの`fol/stm/rot`を実装する。
2. `f_lf/f_lc/f_lr`を原典どおり実装する。
3. maintenance/growth respirationを追加する。
4. allocationのpiecewise分岐を実装する。
5. C3/C4へ一般化する。
6. 同じ18構造プールでNSCだけ異なるnative状態を比較し、18プールがMarkov状態か
   判定する。
7. GPP、NPP、Rh、NEP、Ctotを統合する。

完了条件:

- 各非ゼロフラックスがC原典へ追跡できる。
- 一日炭素収支が成立する。
- native小摂動とreduced 18-poolの誤差を振幅別に示せる。

## 6. 推奨する直近の実装順

直近は炭素-水結合とVISITc水文過程の原典検証を優先する。詳細な式、比較条件、
完了基準は[現在の計画](carbon_water_research_plan.md)に記録した。

1. 完了: current VISITcの最小hydrology/PM C bridgeと一日全出力比較。
2. 完了: `pm_incep/pm_evap/pm_transp`、抵抗・放射・LAI・conductance依存の転記。
3. 完了: native soil retentionとreduced plant pressure-volume/capacitanceの導入。
4. 次: dynamic waterとquasi-steady waterを同一炭素式・同一強制で比較する。
5. 次: PM鎖を制御した気象強制へ接続し、パラメータを較正・範囲化する。
6. 周期基準軌道を作り、その後にのみR-tippingの可否を調べる。

以下のtree一日更新は既存VISIT-matrix側の後続原典接続課題として保持する。

1. prescribed GPPと環境からturnover、maintenance、allocation、growth respirationを
  原典順序で結ぶ。
2. `tpc/tpr`が正の場合のNSC優先蓄積を実装する。
3. leaf emergenceとsurvival reallocationを追加する。
4. NSCを含むか否かで翌日状態が変わる実験を行う。
5. treeから土壌`li_tf/li_tc/li_tr`への結合を検証する。

P3により`frl/frh`の時変係数はnative Cまで検証でき、P4で時間スケールとQSE
感度を取り出し、比較図へ固定した。P5では同じ「純粋関数、手計算、native照合」
の順序を植物過程にも適用する。
