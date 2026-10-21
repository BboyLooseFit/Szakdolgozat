# 4. Binning a SemiBin2 segítségével

## 4.1. A lépés célja

A Flye assembly után a következő lépésben a kontigokat mikrobiális eredetük szerint csoportosítottuk.

Az így létrejövő bin-ek a további MAG-minőségellenőrzés és referencia-alapú összehasonlítás kiindulási anyagát jelentik.

## 4.2. Bemeneti adatok

A SemiBin2 két fő bemenetet kapott:

- Flye assembly: `results/assembly/flye/assembly.fasta`
- Az assembly-re visszamappelt readeket tartalmazó BAM-fájl: `results/assembly/mapping/top10_vs_assembly.bam`

Az assembly 168 kontigból állt, amelyek összesített hossza 37 936 199 bp volt.

A BAM-fájl a szimulált readek és az assembly kontigjai közötti illesztéseket tartalmazta. Az ebből származó lefedettségi információt a SemiBin2 a kontigok csoportosításához is felhasználta.

## 4.3. A használt eszköz és a beállítások

- Eszköz: SemiBin2
- Verzió: 2.5.0
- Futtatási mód: `single_easy_bin`
- Környezet: `human_gut`
- Szekvenálási típus: `long_read`
- Assembly bemenet: Flye által előállított FASTA
- Mapping bemenet: az assembly-re visszamappelt BAM
- Kimeneti könyvtár: `results/binning/semibin2`

A `single_easy_bin` módot egyetlen metagenomikai minta binningjéhez használtuk. A `human_gut` környezeti beállítást a vizsgált közösség bélmikrobiom-eredete miatt választottuk.

## 4.4. A futtatási parancs

A SemiBin2 programot a következő paranccsal futtattuk:

```bash id="5wdfmk"
mkdir -p results/binning/semibin2

SemiBin2 single_easy_bin \
  --environment human_gut \
  --sequencing-type long_read \
  -i results/assembly/flye/assembly.fasta \
  -b results/assembly/mapping/top10_vs_assembly.bam \
  -o results/binning/semibin2
```

## 4.5. A binning eredményei

A SemiBin2 sikeresen befejezte a futást, és összesen 8 bin-t állított elő.

A binning összesített eredményei:

- Assembly kontigjainak száma: 168
- Bin-ekbe rendelt kontigok száma: 150
- Bin-ekbe nem rendelt kontigok száma: 18
- Bin-ekbe rendelt összes szekvenciahossz: 36 956 013 bp
- Bin-ekbe nem rendelt szekvenciák összesített hossza: 980 186 bp

A bin-ekbe sorolt kontigok az eredeti assembly teljes bázishosszának megközelítőleg 97,4%-át képviselték. A kontigok számát tekintve a 168 kontigból 150 került bin-be, ami körülbelül 89,3%-os aránynak felel meg.

### A nyolc bin főbb jellemzői

| Bin | Kontigok száma | Teljes hossz (bp) |
|---|---:|---:|
| SemiBin_0 | 1 | 5 008 923 |
| SemiBin_1 | 4 | 4 046 964 |
| SemiBin_2 | 3 | 4 857 365 |
| SemiBin_3 | 2 | 2 952 929 |
| SemiBin_4 | 4 | 2 803 628 |
| SemiBin_5 | 17 | 6 063 346 |
| SemiBin_6 | 49 | 4 627 992 |
| SemiBin_7 | 70 | 6 594 866 |

A bin-ek FASTA-formátumú kimenetei a következő könyvtárba kerültek:

`results/binning/semibin2/output_bins/`

A bin-ek különálló, tömörített FASTA-fájlokként, például `SemiBin_0.fa.gz` néven álltak rendelkezésre.

## 4.6. Az eredmények értelmezése

A SemiBin2 nyolc bin-t állított elő, de ebből önmagában nem következik, hogy mind a nyolc teljes, jó minőségű MAG-ot reprezentál, vagy hogy az eredeti tíz referencia-genom mindegyikét sikerült különálló bin-ként visszanyerni.

A bin-ek kontigszáma és összhossza jelentősen különbözött. A `SemiBin_0` egyetlen, körülbelül 5 Mb-os kontigból állt, míg a `SemiBin_7` 70 kontigot tartalmazott. Ezek a különbségek az assembly töredezettségére és a bin-ek szerkezetére utalnak, de önmagukban nem elegendők a biológiai minőség megítéléséhez.

A binning eredményét ezért a következő lépésben külön minőségellenőrzésnek vetettük alá.

## 4.7. Kapcsolat a következő munkafolyamattal

A nyolc bin-t a CheckM2 programmal vizsgáltuk tovább. A CheckM2 eredményei alapján megbecsültük a bin-ek teljességét és szennyezettségét.

Ezután a bin-eket a szimuláció során használt eredeti referencia-genomokhoz hasonlítottuk FastANI segítségével, majd a kiválasztott referencia-genomok mellett QUAST-értékelést is végeztünk.

A cél az volt, hogy ne csak a bin-ek számát határozzuk meg, hanem azt is értékeljük, milyen minőségű genomjelölteket sikerült visszanyerni, és ezek mennyire felelnek meg az eredeti szimulációs közösség összetételének.
