# Ejemplo trabajado: doublecortin (DCX) en Cypriniformes

- `family.faa`: las 79 proteínas de la familia DCX/DCLK (la misma base `fam.dmnd` del censo).
- `proteins.tsv`: cada proteína asignada a DCX, dclk1a, dclk2a, dclk3 u other, con las reglas de nombre de `tests/datos/parse_dmnd_original.py`.
- `genes.tsv`: modo DIAMOND (controles dclk1a y dclk2a; dclk3 excluido porque DIAMOND le asigna lecturas de más).
- `genes_sensible.tsv` y `queries_sensible.faa`: modo sensible (controles dclk1a, dclk2a y dclk3), como en la Tabla S2.
- `runs.tsv`: las cuatro especies de estado conocido.
- `hits/`: tablas DIAMOND de esas cuatro especies, para probar todo sin descargar (15 M spots para *D. cerebrum* y *P. micromegethes*, 40 M para *I. furcatus* y 35 M para *P. dabryanus*).
  - *D. cerebrum* ERR3284980 y *P. micromegethes* SRR11971489, el par de calibración.
  - *I. furcatus* SRR18966194 y *P. dabryanus* SRR23104325, los controles independientes.
  - La de *P. micromegethes* se generó el 2026-11-16 con `readloss test --stream` (no se había guardado la del censo).
