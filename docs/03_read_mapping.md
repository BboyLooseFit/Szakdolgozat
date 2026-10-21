# 3. Read mapping a minimap2 segítségével

## 3.1. A lépés célja

A de novo assembly után a szimulált readeket visszamappeltük a Flye által előállított kontigokra.

A cél az volt, hogy megvizsgáljuk, a szimulált readek milyen arányban illeszkednek az elkészült assembly-re, valamint előállítsuk azt a BAM-fájlt, amelyet a későbbi binning során a SemiBin2 is felhasználhatott.

## 3.2. Bemeneti adatok

A folyamat két fő bemenete:

- Szimulált ONT readek: `results/one_sample/fastq/top10.fq.gz`
- Flye assembly: `results/assembly/flye/assembly.fasta`

A FASTQ-fájl 67 341 readet tartalmazott, míg a Flye által létrehozott assembly 168 kontigból állt.

## 3.3. A használt eszközök és beállítások

- Illesztőprogram: minimap2
- Verzió: 2.31-r1302
- Illesztési mód: `map-ont`
- BAM-fájl rendezése és indexelése: samtools
- Rendezett BAM-kimenet: `results/assembly/mapping/top10_vs_assembly.bam`
- Indexfájl: `results/assembly/mapping/top10_vs_assembly.bam.bai`

A minimap2 'map-ont' módját az Oxford Nanopore long-read adatokhoz választottuk.

A minimap2 által előállított SAM-formátumú illesztési adatokat a samtools segítségével rendeztük BAM-formátumba. Ezután indexet készítettünk a BAM-fájlhoz, amely lehetővé teszi a hatékony lekérdezést és a további elemzéseket.

## 3.4. A futtatási parancsok

Először létrehoztuk a kimeneti könyvtárat, majd elvégeztük a read mappinget:

```bash
mkdir -p results/assembly/mapping

minimap2 -ax map-ont \
  results/assembly/flye/assembly.fasta \
  results/one_sample/fastq/top10.fq.gz |
samtools sort -@ 4 \
  -o results/assembly/mapping/top10_vs_assembly.bam

samtools index \
  results/assembly/mapping/top10_vs_assembly.bam
```

A `minimap2 -ax map-ont` parancs az ONT readek és az assembly közötti illesztéseket állította elő. A kimenetet a samtools rendezte BAM-formátumba, majd létrehoztuk a BAM indexét.

## 3.5. Az illesztések statisztikai ellenőrzése

Az illesztések összesített statisztikáit a következő paranccsal vizsgáltuk:

```bash
samtools flagstat \
  results/assembly/mapping/top10_vs_assembly.bam
```

A kontigonkénti lefedettségi mutatók előállításához a következő parancsot használtuk:

```bash
samtools coverage \
  results/assembly/mapping/top10_vs_assembly.bam \
  > results/assembly/mapping/coverage.tsv
```

A 'samtools flagstat' az illesztések általános statisztikáit adja meg, míg a 'samtools coverage' kontigonként szolgáltat információt többek között a lefedettség mértékéről és a leképezett readek számáról.

## 3.6. Az eredmények

A read mapping eredményei alapján:

- Bemeneti readek száma: 67 341
- Elsődleges leképezett readek száma: 67 047
- Elsődleges leképezési arány: 99,56%
- Kimenet: rendezett és indexelt BAM-fájl

Az elsődleges leképezési arány azt mutatja, hogy a bemeneti readek mekkora hányadához tartozott elsődleges illesztés.

A magas elsődleges leképezési arány azt jelzi, hogy a szimulált readek nagy része illeszkedett az összeállított kontigokhoz. Önmagában azonban nem bizonyítja, hogy az assembly teljes vagy hibamentes.

## 3.7. A BAM-fájl szerepe a binningben

A read mapping során létrehozott `top10_vs_assembly.bam` fájlt a SemiBin2 binning bemeneteként is felhasználtuk.

A binning során az assembly kontigjai képezték a csoportosítás alapját, a BAM-fájlból származó read-lefedettségi információ pedig kiegészítő információt biztosított a kontigok csoportosításához.

Fontos megkülönböztetni ezt a BAM-fájlt a MeSS szimuláció során mentett `results/one_sample/bam/top10.bam` fájltól. Az utóbbi a szimulációs kimenet része, míg a `top10_vs_assembly.bam` kifejezetten az eredeti szimulált readek és a Flye-assembly közötti mappingből készült.

## 3.8. Kapcsolat a következő munkafolyamattal

A read mapping eredményeit a SemiBin2 a binning során használta fel. A következő lépésben az assembly kontigjait mikrobiális eredetük szerint csoportosítottuk.
