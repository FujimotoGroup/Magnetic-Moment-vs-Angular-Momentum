## 依存関係
zheev in LAPACK
-qmkl が必要．

# GTS Libraryのインストールのしかた
`sudo apt-get install darcs` でバージョン管理 darcs をインストール

`darcs get http://gerris.dalembert.upmc.fr/darcs/gts-stable` で darcs を使って GTS stable を取得

`sudo apt-get install libglib2.0-dev` でglibをインストール

あとは

`cd gts-stable`

`sh autogen.sh`

`make`

`sudo make install`

