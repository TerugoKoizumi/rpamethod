# RPA計算の手順

## 5. Cu(111)スラブのRPA計算

3.で行なった計算をrpaによって計算する手法について詳しく説明します。

まず以下のURLよりSiの凝集エネルギーの計算を実行して見てください。

https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/energetics/rpa_ex/rpa.html

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

・GPAW

https://wiki.fysik.dtu.dk/gpaw/documentation/basic.html#manual-eigensolver

・RPA

https://wiki.fysik.dtu.dk/gpaw/documentation/xc/rpa.html#rpa

これらのURL先にパラメータの詳細が記載されているので確認してください。
基本Quantum Espressoと扱いはほとんど同じなので、単位にだけ気を付けてください。

できると思います。
