# 1. Adatok szimulálása a MeSS segítségével

## 1.1. A lépés célja

A munkafolyamat első lépése egy ismert összetételű mikrobiális közösségből származó, Oxford Nanopore Technologies (ONT) long-read szekvenciaadatok szimulálása volt. A szimuláció célja olyan adathalmaz előállítása volt, amelyen a későbbi assembly-, binning- és MAG-minőségértékelési lépések vizsgálhatók.

## 1.2. Kiindulási adatok

A szimulációhoz a Human Microbiome Project (HMP) eredetű `gut/SRS104912.tsv` template-et használtuk, amely 720 referencia-genomot tartalmazott.

A közösség komplexitásának csökkentése érdekében a 10 legnagyobb abundanciájú genomot választottuk ki. A kiválasztott genomok relatív abundanciáit a részhalmazon belül normalizáltuk. 

A létrehozott közösség azonosítója `top10` volt.

## 1.3. A szimuláció paraméterei

- Technológia: ONT long-read
- Hibamodell: `QSHMM-ONT-HQ`
- Célzott átlagos readhossz: 15 000 bp
- Beállított szekvenciapontosság: 0,95
- Arányparaméter: `39:24:36`
- Véletlenszám-generátor seed értéke: 42
- Párhuzamos szálak száma: 4
- BAM-kimenet mentése: engedélyezve

## 1.4. A szimuláció futtatása

A MeSS futtatásához az alábbi parancsot használtuk:

```bash
mess run \
  -i communities/top10_community.tsv \
  --accession \
  --tech ont \
  --model QSHMM-ONT-HQ \
  --mean-len 15000 \
  --accuracy 0.95 \
  --ratio 39:24:36 \
  --seed 42 \
  --bam \
  --threads 4 \
  -o results/one_sample
```

A `--seed 42` beállítás a szimuláció reprodukálhatóságát segíti. 

## 1.5. A szimuláció eredményei

A MeSS futása sikeresen befejeződött. A kapott FASTQ-adatok főbb jellemzői:

- Readek száma: 67 341
- Összes szekvenciahossz: 1 000 000 962 bp, azaz megközelítőleg 1 Gb
- Átlagos readhossz: 14 849,8 bp
- Minimális readhossz: 92 bp
- Maximális readhossz: 66 872 bp

A szimuláció kimenete tartalmazta a FASTQ-adatokat és a további vizsgálatokhoz felhasználható BAM-fájlt is.

## 1.6. Reprodukálhatóság ellenőrzése

A szimulációt azonos paraméterekkel és ugyanazzal a seed értékkel kétszer futtattuk le. Azonos eredményeket kaptunk, 'samtools flagstat' statisztikák megegyeztek.

## 1.7. A lépés kimenete és kapcsolata a következő fázissal

A szimuláció eredménye egy ismert összetételű közösségből származó ONT long-read adathalmaz lett. Ezeket a readeket a következő lépésben a Flye programmal, de novo assembly segítségével állítottuk össze hosszabb kontigokká.
