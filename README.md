## mainを実行する場合
プロジェクトのルートディレクトリで`./bin/main`とする．

## 依存関係
LAPACKにある`zheev`を使っている．
コンパイル時にオプション`-qmkl`が必要．

またCGALも使用している．

## CGAL の準備
macだと`brew install cgal`で一発．
ubuntuでも`sudo apt install libcgal-dev`でよいらしい．
