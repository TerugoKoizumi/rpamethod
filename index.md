# RPA計算の手順

基本的に私の修士論文（ tkoizumi@smith2：/home/tkoizumi/Thesis/Master/修士論文.pdf ）
のAppendixに詳細を記載していますが、ここではより端的に分かり易く記述します。

<br />

## 1. ASE・GPAWのインストール

これは以下のURLより手順を踏んで行います。

https://qiita.com/ikutaro_hamada/items/3fcab08ea7e163967c3e

<br />

## 2. GPAWを用いた格子定数の計算

GPAWの使い方をCu(111)表面のエネルギー計算を例に挙げて説明します。
まずCuバルクの格子定数を求めるために、カットオフエネルギーおよびk点メッシュに対する格子定数の収束を図るスクリプトを以下に示します。


### A. カットオオフエネルギーを変化させてそれぞれの格子定数の値を計算する

以下のスクリプトを作成します。ここでは `Cu_fcc.py` とします。
```
import numpy as np
from ase.build import bulk
from gpaw import GPAW, PW, MethfesselPaxton, Davidson

a0 = 3.63
cu = bulk('Cu', 'fcc', a=a0)
cell0 = cu.cell

for ecut in range(200, 801, 40):
    cu.calc = GPAW(mode=PW(ecut),
                   xc='PBE',
                   kpts=(12, 12, 12),
                   occupations=MethfesselPaxton(width=0.1),
                   eigensolver='dav',
                   convergence={'energy': 2.0e-9},
                   txt=f'Cu-{ecut}.txt')
    for eps in np.linspace(-0.02, 0.02, 10):
        cu.cell = (1 + eps) * cell0
        cu.get_potential_energy()
```
一番上のimport,fromの部分でこの計算に必要となるモジュールを読み込みます。

```
a0 = 3.63
cu = bulk('Cu', 'fcc', a=a0)
cell0 = cu.cell
```
↑この3行で初期構造を決定します。

```
    cu.calc = GPAW(mode=PW(ecut),
                   xc='PBE',
                   kpts=(12, 12, 12),
                   occupations=MethfesselPaxton(width=0.1),
                   eigensolver='dav',
                   convergence={'energy': 2.0e-9},
                   txt=f'Cu-{ecut}.txt')
```
↑この部分がGPAWによる計算コードです。どのパラメータもQEと対応させられますが、
単位がQEでは（bohr）（Ry）であるのに対し、GPAWでは（Å）（eV）がデフォルトであることに注意してください。
forループでカットオフエネルギーを変えながら以下のジョブスクリプト(ここでは`run.sh`とします)より計算を回すと、以下のようにテキストファイルが出力されます。

```
#!/bin/sh
#SBATCH -J Cu_bulk
#SBATCH -p i8cpu
#SBATCH -N 1
#SBATCH -n 64

module load intel_mpi/2020.4.304 intel_compiler/2020.4.304 intel_mkl/2020.4.304

source ${HOME}/ASE/ase/bin/activate

srun python3 ./Cu_fcc.py
```
```
sbatch run.sh
```

```
Cu-200.txt  Cu-280.txt  Cu-360.txt  Cu-440.txt  Cu-520.txt  Cu-600.txt  Cu-680.txt  Cu-760.txt
Cu-240.txt  Cu-320.txt  Cu-400.txt  Cu-480.txt  Cu-560.txt  Cu-640.txt  Cu-720.txt  Cu-800.txt
```

続いてこれらの結果をプロットするため次のジョブを実行します。ここでは`lattice.py`とします。
```
import matplotlib.pyplot as plt
import numpy as np
from ase.build import bulk
from ase.eos import EquationOfState
from ase.io import read

def fit(filename):
    configs = read(filename + '@:')
    volumes = [a.get_volume() for a in configs]
    energies = [a.get_potential_energy() for a in configs]
    eos = EquationOfState(volumes, energies)
    v0, e0, B = eos.fit()
    return (4 * v0)**(1 / 3.0)

cutoffs = range(200, 801, 40)
a = [fit(f'Cu-{ecut}.txt') for ecut in cutoffs]

print(*a)
```
```
python3 lattice.py
```
それぞれのカットオフエネルギーに対応する格子定数の値が出力されるので、エクセルでグラフを作ってください。

<br />

### B. k点メッシュを変化させてそれぞれの格子定数の値を計算する

ほとんどAと同じです。`Cu_fcc.py`、`lattice.py`を以下のように変えて同じように実行すればOKです。
```
import numpy as np
from ase.build import bulk
from gpaw import GPAW, PW, MethfesselPaxton, Davidson

a0 = 3.63
cu = bulk('Cu', 'fcc', a=a0)
cell0 = cu.cell

for k in range(4, 21):
    cu.calc = GPAW(mode=PW(640),
                   xc='PBE',
                   kpts=(k, k, k),
                   occupations=MethfesselPaxton(width=0.1),
                   eigensolver='dav',
                   convergence={'energy': 2.0e-9},
                   txt=f'Cu-{k:02}.txt')
    for eps in np.linspace(-0.02, 0.02, 10):
        cu.cell = (1 + eps) * cell0
        cu.get_potential_energy()
```
```
import matplotlib.pyplot as plt
import numpy as np
from ase.eos import EquationOfState
from ase.io import read

def fit(filename):
    configs = read(filename + '@:')
    volumes = [a.get_volume() for a in configs]
    energies = [a.get_potential_energy() for a in configs]
    eos = EquationOfState(volumes, energies)
    v0, e0, B = eos.fit()
    return (4 * v0)**(1 / 3.0)

kpoints = range(4, 21)
a = [fit(f'Cu-{k:02}.txt') for k in kpoints]

print(*a)
```


