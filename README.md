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
3. `git add (編集したファイル名)`
4. `git commit -m "コミットメッセージ"`
5. `git push`

# 目次

1. 間違える数値計算（[nonrigorous_numerics.ipynb](https://www.risk.tsukuba.ac.jp/~takitoshi/tutorial/nonrigorous_numerics.html)）
1. 整数・浮動小数点数（[floating-point.ipynb](https://www.risk.tsukuba.ac.jp/~takitoshi/tutorial/floating-point.html)）
1. 丸め誤差・その他の誤差（[rounding-error.ipynb](https://www.risk.tsukuba.ac.jp/~takitoshi/tutorial/rounding-error.html)）
1. 区間演算（[interval-arithmetic.ipynb](https://www.risk.tsukuba.ac.jp/~takitoshi/tutorial/interval-arithmetic.html)）
1. ベクトルの内積・行列ベクトル積・行列積（丸めの変更をしないでBlasを使う区間演算）（[interval_dot-mul.ipynb](https://www.risk.tsukuba.ac.jp/~takitoshi/tutorial/interval_dot-mul.html)）
1. 線型方程式（丸めの変更をしないでBlasを使う精度保証）（[verifylss.ipynb](https://www.risk.tsukuba.ac.jp/~takitoshi/tutorial/verifylss.html)）
1. 固有値問題（[verifyalleig.ipynb](https://www.risk.tsukuba.ac.jp/~takitoshi/tutorial/verifyalleig.html)）
1. 高速フーリエ変換の精度保証（井藤さん）
1. 非線形方程式（有限次元：Krawczyk, Interval Newton,大谷くん, verifynlss,ipynb）
1. Newton-Kantorovich type argument (radii-polynomial approach)（大谷くん、有限次元・無限次元両方いける。radii-ppolynomial.ipynb）
1. 常微分方程式の周期解（大谷くんの卒論、ロトカ・ヴォルテラ方程式の周期解をradii-polynomialで検証）

# 構成ファイルの説明（.jlファイルだけ。途中でincludeして使用できる）

- `IntervalLinearAlgebra.jl`: 丸め向きを変更しない区間行列積（`int_mul`）が実装されている。入力のタイプによって区間行列積の場合分けをしている。
