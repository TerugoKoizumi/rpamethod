# RPA計算の手順

基本的に私の修士論文（ tkoizumi@smith2：/home/tkoizumi/Thesis/Master/修士論文.pdf ）
のAppendixに詳細を記載していますが、ここではより端的に分かり易く記述します。

<br />

## 1. ASE・GPAWのインストール

これは以下のURLより手順を踏んで行います。

[インストール手順](https://qiita.com/ikutaro_hamada/items/3fcab08ea7e163967c3e)

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

## 3. GPAWを用いたスラブの計算

Cu(111)スラブモデルを作ってその初期構造のエネルギーを計算します。
以下のジョブを実行します。ここでは、`SCF.py`とします。

```
from ase import Atoms
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from gpaw import GPAW, PW, MethfesselPaxton, Davidson
from ase.build import fcc111

a = 3.6434
N = 6
k = 6

slab = fcc111('Cu', (2, 2, N), a=a, vacuum=23)
slab.center(axis=2)
constraint = FixAtoms(mask=slab.positions[:, 2] < 28.259)
slab.set_constraint(constraint)
calc = GPAW(mode=PW(640),
            nbands=280,
            kpts={'size': (k, k, 1), 'gamma': True},
            xc='PBE',
            occupations=MethfesselPaxton(width=0.1),
            eigensolver='dav',
            convergence={'energy': 2.0e-9},
            txt='scf.txt')
slab.calc = calc
e = slab.get_potential_energy()
calc.write('scf.gpw')

print(e)
```
```
sbatch SCF.py
```

ここで
```
a = 3.6434
N = 6
k = 6

slab = fcc111('Cu', (2, 2, N), a=a, vacuum=23)
slab.center(axis=2)
constraint = FixAtoms(mask=slab.positions[:, 2] < 28.259)
slab.set_constraint(constraint)
```
↑ここでスラブモデルを作成しています。2x2のユニットセルで6層、真空層は23Åに設定しています。
また下2行
```
constraint = FixAtoms(mask=slab.positions[:, 2] < 28.259)
slab.set_constraint(constraint)
```
↑では6層のうち下3層を固定しています。計算の前に初期構造を可視化すると下から数えて3層目のCu原子の高さが28.259Å以下であるためこのような指定になります。ちなみに可視化するスクリプトは以下になります。ここでは`view.py`とします。

```
from ase.visualize import view
from ase import Atoms
from ase.build import fcc111

a = 3.6434
N = 6

slab = fcc111('Cu', (2, 2, N), a=a, vacuum=23)
slab.center(axis=2)

view(slab, repeat=(1, 1, 1))
```
```
python3 view.py
```
このスクリプトに固定する2行を足すとそれも可視化されるので一度確認して見てください。

<br />

`SCF.py`より`scf.gpw`が出力されます。これを用いて構造緩和を進めていきます。
以下のスクリプトを実行します。ここでは`RELAX.py`とします。

```
from ase import Atoms
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from gpaw import GPAW, PW, Mixer, FermiDirac, MethfesselPaxton, Davidson
from gpaw import restart
from ase.build import fcc111
from ase.io import write

slab, calc = restart('scf.gpw', txt='1relax.txt')

e = slab.get_potential_energy()

calc.attach(calc.write, 1, '1relax.gpw')

relax = BFGS(slab, trajectory='1relax.traj')
relax.run(fmax=0.0025)
write('Cu_sf.xyz', slab)

print(e)
```
ここで`1relax.gpw`が出力されます。fmaxで収束閾値を決定します。単位はeVなので注意してください。一度で収束しきらない場合は
```
from ase import Atoms
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from gpaw import GPAW, PW, Mixer, FermiDirac, MethfesselPaxton, Davidson
from gpaw import restart
from ase.build import fcc111
from ase.io import write

slab, calc = restart('1relax.gpw', txt='2relax.txt')

e = slab.get_potential_energy()

calc.attach(calc.write, 1, '2relax.gpw')　

relax = BFGS(slab, trajectory='2relax.traj')
relax.run(fmax=0.0025)
write('Cu_sf.xyz', slab)

print(e)
```
とするとrestartできます。収束した場合のみ`Cu_sf.xyz`が出力されます。このxyzファイルは吸着系の計算に使います。
また出力される`1relax.txt`に結果が書いてあるので確認してください。

<br />

## 4. GPAWを用いた水分子および吸着系の計算
計算の方法はほとんど同じです。初期構造の作り方などは修士論文を見ていただけると大体分かるかと思います。

## 5. Cu(111)スラブのRPA計算

3.で行なった計算をrpaによって計算する手法について詳しく説明します。

まず以下のURLよりSiの凝集エネルギーの計算を実行して見てください。

[Siを用いたRPA計算例](https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/energetics/rpa_ex/rpa.html)

この流れでCu(111)スラブの計算も進めていきます。まず以下のスクリプトを作成します。

```
from ase import Atoms
   from ase.build import fcc111
   from gpaw import GPAW, PW, FermiDirac, MixerSum, Davidson
   from gpaw.xc.rpa import RPACorrelation
   from ase.build import molecule
   from ase.parallel import paropen
   from gpaw.hybrids.energy import non_self_consistent_energy as nsc_energy

d = 2.56
a = 2**0.5 * d

slab = fcc111('Cu', a=a, size=(1, 1, 4), vacuum=10.0)
   slab.pbc = True
   slab.calc = GPAW(xc='PBE',
                    mode=PW(600),
                    basis='dzp',
                    eigensolver=Davidson(niter=4),
                    nbands='200%',
                    kpts={'size': (12, 12, 1), 'gamma': True},
                    occupations=FermiDirac(width=0.05),
                    convergence={'density': 1e-5},
                    parallel={'domain': 1},
                    mixer=MixerSum(0.05, 5, 50),
                    txt=f'pbe_output.txt')

   E_pbe = slab.get_potential_energy()

   slab.calc.write('Cu.gpw', mode='all')

   E_hf = nsc_energy('Cu.gpw', 'EXX').sum()
   slab.calc.diagonalize_full_hamiltonian()

   slab.calc.write('Cu.all.gpw', mode='all')

   rpa = RPACorrelation('Cu.all.gpw',
                        txt='RPA.txt',
                        skip_gamma=True,
                        frequency_scale=2.0)

   E_rpa = rpa.calculate(ecut=[200])[0]

   e = E_pbe + E_rpa

   print(e, E_pbe, E_rpa)
```

ここでは1x1のユニットセルで4層、真空層10Åのモデルを利用しています。
```
   E_hf = nsc_energy('Cu.gpw', 'EXX').sum()
   slab.calc.diagonalize_full_hamiltonian()

   slab.calc.write('Cu.all.gpw', mode='all')

   rpa = RPACorrelation('Cu.all.gpw',
                        txt='RPA.txt',
                        skip_gamma=True,
                        frequency_scale=2.0)

   E_rpa = rpa.calculate(ecut=[200])[0]

   e = E_pbe + E_rpa

   print(e, E_pbe, E_rpa)
```
この部分が通常の計算に対する追加項目になります。

```
   E_hf = nsc_energy('Cu.gpw', 'EXX').sum()
   slab.calc.diagonalize_full_hamiltonian()
```
この2行よりハートリーフォックによる交換エネルギーを計算します。この後
```
  slab.calc.write('Cu.all.gpw', mode='all')
```
↑より新たにgpwファイルを出力します。
このgpwファイルを用いてRPA相関エネルギーを計算します。
```
  rpa = RPACorrelation('Cu.all.gpw',
                        txt='RPA.txt',
                        skip_gamma=True,
                        frequency_scale=2.0)

   E_rpa = rpa.calculate(ecut=[200])[0]
```
↑この項目がRPA計算の部分になります。

`ecut=[200]`を`RPACorrelation`の（）内に入れると上手く回らないので`rpa.calculate`のパラメータとして扱ってください。

## 6. GPAWとRPAのマニュアル

[GPAW](https://wiki.fysik.dtu.dk/gpaw/documentation/basic.html#manual-eigensolver)

[RPA](https://wiki.fysik.dtu.dk/gpaw/documentation/xc/rpa.html#rpa)

これらのURL先にパラメータの詳細が記載されているので確認してください。
基本Quantum Espressoと扱いはほとんど同じなので、単位にだけ気を付けてください。

できると思います。






