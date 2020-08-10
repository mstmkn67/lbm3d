# lbm3d

３次元格子ボルツマン法のプログラムです。
[OCTA](http://octa.jp/)のインターフェースであるGOURMETで動きます。
GOURMETのメニューから、
Tool>Environments Setup
から、ダウンロードしたディレクトリを指定して使ってください。

実行ファイルは、Ubuntu18.04.5でコンパイルしています。
lbm3d -I inputfile -O outputfile
で実行してください。
たとえば、
lbm3d -I shear.udf -O shear_s.udf
のようにです。

------------
![paraview screenshot](screenshot.png)
