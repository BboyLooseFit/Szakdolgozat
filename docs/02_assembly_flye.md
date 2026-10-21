# 2. De novo assembly a Flye segítségével

## 2.1. A lépés célja

A szimulált ONT long-read szekvenciákból a következő lépésben de novo assembly-t készítettünk a Flye programmal.

A cél az volt, hogy a szimulált mikrobiális közösség szekvenciáiból összeállítást készítsünk, amely a későbbi read mapping és binning alapjául szolgál.

## 2.2. Bemeneti adatok

A bemenet a MeSS által generált FASTQ-fájl volt:

`results/one_sample/fastq/top10.fq.gz`

A fájl 67 341 szimulált ONT readet tartalmazott, megközelítőleg 1 Gb összesített szekvenciahosszal.

## 2.3. A használt eszköz és beállítások

- Eszköz: Flye
- Verzió: 2.9.6-b1802
- Assembly-mód: metagenomikai assembly (`--meta`)
- Bemeneti readtípus: `--nano-hq`
- Becsült readhiba: 0,05
- Párhuzamos szálak száma: 4

A MeSS szimulációjában `QSHMM-ONT-HQ` hibamodellt és 0,95-ös beállított szekvenciapontosságot használtunk. A Flye futtatásánál az 0,05-ös `--read-error` paraméter ennek megfelelően 5%-os becsült readhibát jelentett.

## 2.4. A futtatási parancs

A Flye assembly-t a következő paranccsal futtattuk:

```bash
flye \
  --nano-hq results/one_sample/fastq/top10.fq.gz \
  --meta \
  --read-error 0.05 \
  --threads 4 \
  --out-dir results/assembly/flye
```

A fő kimeneti FASTA-fájl:

'results/assembly/flye/assembly.fasta'

A Flye naplófájlja a futás részleteit és az assembly folyamata során keletkezett információkat tartalmazza.

## 2.5. Az assembly főbb eredményei

Az assembly főbb statisztikai mutatói:

- Kontigok száma: 168
- Teljes assembly-hossz: 37 936 199 bp, azaz körülbelül 37,94 Mb
- Legrövidebb kontig: 5 084 bp
- Átlagos kontighossz: 225 810,7 bp
- Leghosszabb kontig: 5 008 923 bp
- N50: 667 985 bp
- GC-tartalom: 44,65%

A kontigstatisztikák összegzésére a 'seqkit stats -a' parancsot használtuk.

## 2.6. Az eredmények értelmezése

A közel 37,94 Mb-os assembly 168 kontigból állt. A 667 985 bp-os N50 azt jelenti, hogy a kontigok hossz szerint rendezett összegzett hosszának legalább a fele ekkora vagy ennél hosszabb kontigokból áll.

Az assembly hossza, az N50, a kontigszám és a GC-tartalom a de novo assembly szerkezeti és alapvető összetételi jellemzését szolgálja.

Ezek a mutatók önmagukban nem bizonyítják, hogy az összes eredeti mikrobiális genom teljesen és hibamentesen rekonstruálódott. A visszanyerés értékeléséhez ezért további lépésekre volt szükség: read mappingre, binningre, MAG-minőségellenőrzésre és referencia-alapú összehasonlításokra.

## 2.7. A következő munkafolyamathoz való kapcsolat

A Flye assembly szolgált a következő lépések bemeneteként. Először a szimulált readeket visszamappeltük a kontigokra, majd a kontigokat és a BAM-fájlt a SemiBin2 binning során használtuk fel.

A további vizsgálatok célja az volt, hogy felmérjük, a 10 eredeti referencia-genomból milyen minőségű, különálló MAG-jelölteket sikerült visszanyerni.
