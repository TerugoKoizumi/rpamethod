# RPA計算の手順

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


