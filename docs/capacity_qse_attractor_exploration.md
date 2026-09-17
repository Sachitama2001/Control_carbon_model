# Carbon storage capacity と QSE・吸引軌道の探索

詳細な導出、反例、図、次の検証計画は [TeX](capacity_qse_attractor_exploration.tex) / [PDF](capacity_qse_attractor_exploration.pdf) にまとめた。

## 中心結果

Luo らの一般式自体は状態依存性を排除していない。本稿の狙いは、Luo らが一意な大域吸引平衡を主張したと批判することではなく、状態依存の場合の capacity を固定点写像として明示し、一意性・安定性・吸引域を比較可能にすることである。

状態依存の炭素系

\[
\dot X=u(X,\lambda)+M(X,\lambda)X
\]

に対し、点ごとの Luo 型 capacity を

\[
\widehat X_c(X,\lambda)=-M(X,\lambda)^{-1}u(X,\lambda)
\]

と置くと、平衡は `現在状態 = capacity` という固定点である。

\[
X^*=\widehat X_c(X^*,\lambda).
\]

平衡での全ヤコビアンは

\[
J_*=M_*\left(I-D_X\widehat X_c|_*\right)
\]

となる。このため、凍結した炭素移動行列 `M` が安定でも、capacity の状態フィードバックが全系の安定性を変え得る。capacity 写像が恒等写像と複数回交差すれば、同じ環境でも複数の自己無撞着 capacity が存在する。

一次元では storage potential `p = capacity - stock` について

\[
\frac{d}{dt}\frac{p^2}{2}=k(\widehat c'(x)-1)p^2
\]

である。したがって、state-dependent capacity の下では storage potential は常に減る距離や大域 Lyapunov 関数ではない。

## 用語上の結論

- 別の安定平衡へ移る場合、その平衡も QSE である。正確な表現は「追跡していた QSE 枝とは別の自己無撞着 capacity への遷移」である。
- 季節外力下の周期軌道は瞬間 QSE と異なるが、通常の位相遅れでも生じる。
- 定常環境下で文字どおり QSE でない吸引対象を主張するには、安定周期軌道などを示す必要がある。VISIT ではまだ示していない。

## 再現可能な反例

`capacity_feedback.py` の正の一次元例では、全ての環境値で二つの非ゼロ安定枝が存続する。環境を 20 年で動かすと低位枝、120 年で動かすと高位枝へ達し、臨界継続時間は約 75.003 年である。これは数理的な可能性の証明であり、VISIT の式でも森林の較正値でもない。

```bash
PYTHONPATH=src python examples/analyze_capacity_feedback.py
python -m pytest tests/test_capacity_feedback.py -q
```

## 次の方向

VISIT の日次植物更新をソース順に完成し、一定環境の固定点と季節環境の一年写像を、多初期値・継続・安定性・吸引域の順に調べる。複数吸引状態が確認された場合に限り、気象変数または環境パラメータを一つずつランプして R-tipping を検定する。
