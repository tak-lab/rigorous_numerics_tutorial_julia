# Julia言語を使った精度保証付き数値計算のチュートリアル

精度保証付き数値計算の簡単な方法をJulia言語を用いて紹介する。

# ゼミ参加者（2020-07-09）

- 舩越
- 井藤
- 大谷
- 高安

# gitの基本操作

1. `git pull`
2. （各自ファイルを変更）
3. `git add -A`
4. `git commit -m "コミットメッセージ"`
5. `git push`

# 目次

1. 整数・浮動小数点数（[floating-point.ipynb](https://www.risk.tsukuba.ac.jp/~takitoshi/tutorial/floating-point.html)）
1. 丸め誤差・その他の誤差（[rounding-error.ipynb](https://www.risk.tsukuba.ac.jp/~takitoshi/tutorial/rounding-error.html)）
1. 区間演算（Interval-arithmetic.ipynb）
1. ベクトルの内積・行列ベクトル積・行列積（Blasを使う区間演算->できない！、丸めの変更をしない方法）（nearest-IntervalArithmetic.ipynb）
1. 線型方程式（みんなで）
1. 固有値問題（井藤さん）
1. 高速フーリエ変換の精度保証（井藤さん）
1. 非線形方程式（大谷くん）
1. 常微分方程式の周期解（大谷くん）
