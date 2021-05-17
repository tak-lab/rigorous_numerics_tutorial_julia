# Julia言語を使った精度保証付き数値計算のチュートリアル

精度保証付き数値計算の標準的な方法をJulia言語を用いて紹介する。

# 目次

1. はじめに -間違える数値計算-（[nonrigorous_numerics.ipynb](https://www.risk.tsukuba.ac.jp/~takitoshi/tutorial/nonrigorous_numerics.html)）
1. 整数・浮動小数点数（[floating-point.ipynb](https://www.risk.tsukuba.ac.jp/~takitoshi/tutorial/floating-point.html)）
1. 丸め誤差・その他の誤差（[rounding-error.ipynb](https://www.risk.tsukuba.ac.jp/~takitoshi/tutorial/rounding-error.html)）
1. 区間演算 -精度保証付き数値計算の入り口-（[interval-arithmetic.ipynb](https://www.risk.tsukuba.ac.jp/~takitoshi/tutorial/interval-arithmetic.html)）
1. ベクトルの内積・行列ベクトル積・行列積の区間演算（[interval_dot-mul.ipynb](https://www.risk.tsukuba.ac.jp/~takitoshi/tutorial/interval_dot-mul.html)）
1. 線型方程式の解の精度保証付き数値計算（[verifylss.ipynb](https://www.risk.tsukuba.ac.jp/~takitoshi/tutorial/verifylss.html)）
1. 標準固有値問題の精度保証付き数値解法（[verifyalleig.ipynb](https://www.risk.tsukuba.ac.jp/~takitoshi/tutorial/verifyalleig.html)）
1. 高速フーリエ変換の精度保証（[verifyfft.ipynb](https://www.risk.tsukuba.ac.jp/~takitoshi/tutorial/verifyfft.html)）
1. 非線形方程式（[verifynlss.ipynb](https://www.risk.tsukuba.ac.jp/~takitoshi/tutorial/verifynlss.html)）
1. Newton-Kantorovich型定理 (radii-polynomial approach)（[Newton-Kantorovich.ipynb](https://www.risk.tsukuba.ac.jp/~takitoshi/tutorial/Newton-Kantorovich.html)）
1. フーリエ級数、チェビシェフ級数、離散畳み込み etc.

# 構成ファイルの説明（.jlファイルだけ。途中でincludeして使用している）

- `IntervalLinearAlgebra.jl`: 丸め向きを変更しない区間行列積（`int_mul`）が実装されている。入力のタイプによって区間行列積の場合分けをしている。
- 今後、区間連立一次方程式、非線形方程式、高速フーリエ変換などを適宜追加予定。

## 筆者

本資料は、[高安研究室](http://www.taklab.org/)のゼミ資料として作成しています。これまで貢献した人々は以下の通りです。

- [舩越康太](https://github.com/2754github)
- 井藤佳奈子
- 大谷俊輔
- [高安亮紀](https://www.risk.tsukuba.ac.jp/~takitoshi/)

### gitの基本操作(内部向け)

1. `git pull`
2. （各自ファイルを変更）
3. `git add (編集したファイル名)`
4. `git commit -m "コミットメッセージ"`
5. `git push`
