# 6. Referencia-alapú értékelés FastANI és QUAST segítségével

## 6.1. A vizsgálat célja

A SemiBin2 által létrehozott bin-ek és a CheckM2 minőségellenőrzése után a következő lépésben a szimulációhoz használt eredeti referencia-genomokhoz hasonlítottuk a visszanyert szekvenciákat.

Mivel a szimulációs közösség összetétele ismert volt, lehetőség nyílt annak vizsgálatára, hogy a 10 eredeti referencia-genom közül melyeket sikerült különálló bin-ként visszanyerni, valamint milyen mértékű a hasonlóság a visszanyert szekvenciák és az eredeti referencia-genomok között.

Ehhez két eszközt használtunk:

- FastANI: a nukleotidszintű hasonlóság, vagyis az Average Nucleotide Identity (ANI) becslésére.
- QUAST: az assembly-k és bin-ek referencia-alapú szerkezeti és lefedettségi értékelésére.

A két eszköz eltérő tulajdonságokat mér, ezért az eredményeket egymással és a CheckM2 mutatóival együtt értékeltük.

## 6.2. Bemeneti adatok

A FastANI- és QUAST-vizsgálatokhoz a következő adatokat használtuk:

- A SemiBin2 által létrehozott nyolc bin FASTA-fájljai.
- A MeSS által használt tíz eredeti referencia-genom.
- A teljes Flye-assembly, amely 168 kontigból állt.
- A referencia-genomokhoz tartozó taxonómiai metaadatok.

A bin-eket a `results/binning/semibin2/output_bins/` könyvtárban található FASTA-fájlok reprezentálták.

Az eredeti referencia-genomok a `results/one_sample/assembly_finder/download/` könyvtárban voltak elérhetők.

## 6.3. FastANI: a genomok közötti hasonlóság vizsgálata

### A módszer célja

A FastANI az összehasonlított genomok közötti átlagos nukleotidazonosságot becsüli. Az ANI magas értéke azt jelzi, hogy az illeszkedő szekvenciarégiók nagyon hasonlóak egymáshoz.

A FastANI kimenete az ANI mellett az összehasonlításban illeszkedő fragmensek számát és a vizsgált fragmensek teljes számát is tartalmazza. Ezek segítségével nemcsak a hasonlóság mértékét, hanem az egyezés kiterjedését is értékeltük.

### A futtatás

A bin-ek és referencia-genomok közötti összehasonlítást a következő paranccsal végeztük:

```bash
micromamba run -n fastani fastANI \
  --ql results/qc/fastani/query_bins.txt \
  --rl results/qc/fastani/reference_genomes.txt \
  --threads 4 \
  -o results/qc/fastani/bin_vs_reference.tsv
```

Ezután a referencia-genomok és a bin-ek szerepét megcserélve fordított összehasonlítást is végeztünk:

```bash
micromamba run -n fastani fastANI \
  --ql results/qc/fastani/reference_genomes.txt \
  --rl results/qc/fastani/query_bins.txt \
  --threads 4 \
  -o results/qc/fastani/reference_vs_bins.tsv
```

A kétirányú összehasonlítást azért használtuk, mert a FastANI fragmensalapú illeszkedési mutatója a query és a referencia szerepétől függően eltérhet.

### A referencia-genomok hozzárendelése

A vizsgálat alapján nyolc referencia-genomhoz találtunk olyan bin-t, amely megfelelt az alkalmazott, projekt-specifikus hozzárendelési feltételnek.

A feltétel három elemből állt:

1. A bin és a referencia kölcsönösen a legjobb FastANI-találatot adta.
2. Az ANI legalább 95% volt.
3. A referencia fragmenseinek legalább 50%-a illeszkedett a vizsgált bin-hez a FastANI mérőszáma alapján.

Ez a szabály a projektben használt operatív döntési kritérium, nem univerzális MAG-minőségi szabvány.

| Referencia-genom | Legjobb bin | ANI (%) | Referencia-fragmensek illeszkedése (%) | Azonosítható |
|---|---|---:|---:|---:|
| Phocaeicola vulgatus | SemiBin_0 | 100,00 | 99,52 | 1 |
| Segatella copri | SemiBin_1 | 99,98 | 95,18 | 1 |
| Parabacteroides merdae | SemiBin_2 | 99,97 | 96,87 | 1 |
| Bacteroides xylanisolvens | SemiBin_5 | 99,53 | 89,38 | 1 |
| Ruminococcus bicirculans | SemiBin_3 | 99,97 | 99,19 | 1 |
| Bacteroides fragilis | SemiBin_7 | 98,75 | 77,62 | 1 |
| Faecalibacterium prausnitzii | SemiBin_4 | 99,85 | 98,52 | 1 |
| Bacteroides stercoris | SemiBin_6 | 99,39 | 84,31 | 1 |
| Phocaeicola dorei | SemiBin_0 | 95,07 | 70,03 | 0 |
| Bacteroides thetaiotaomicron | SemiBin_7 | 94,66 | 42,95 | 0 |

## 6.4. QUAST: a bin-ek referencia-alapú értékelése

### A módszer célja

A QUAST-tal a bin-ek assembly-jellemzőit és a kiválasztott referencia-genom lefedettségét vizsgáltuk.

A legfontosabb értékelt mutatók:

- Genome fraction (%): a referencia-genomnak az assembly által lefedett hányada.
- Total length (bp): a bin összesített szekvenciahossza.
- Kontigszám: az assembly-t alkotó különálló szekvenciák száma.
- N50 és NG50: az assembly kontighossz-eloszlását jellemző mutatók.
- Misassemblies: a referencia mentén észlelt szerkezeti eltérések száma.

A genome fraction a QUAST által meghatározott referencia-lefedettség, és nem azonos a FastANI fragmenseinek illeszkedési arányával.

### A futtatás

A QUAST-ot külön-külön futtattuk minden bin-re, a FastANI alapján kiválasztott legjobb referencia-genommal. A futtatás általános formája:

```bash
micromamba run -n quast311 quast.py \
  -r REFERENCIA_GENOM.fna \
  BIN.fa \
  --threads 2 \
  -o results/qc/quast/reports/BIN_NEVE
```

A 'REFERENCIA_GENOM.fna' és a 'BIN.fa' helyére a megfelelő referencia- és bin-fájl került.

### A nyolc bin főbb QUAST-eredményei

| Bin | Genome fraction (%) | Teljes hossz (bp) | Kontigok száma | N50 (bp) | NG50 (bp) | Misassembly-k száma |
|---|---:|---:|---:|---:|---:|---:|
| SemiBin_0 | 99,87 | 5 008 923 | 1 | 5 008 923 | 5 008 923 | 1 |
| SemiBin_1 | 98,02 | 4 046 964 | 4 | 1 476 343 | 1 476 343 | 1 |
| SemiBin_2 | 97,21 | 4 857 365 | 3 | 3 181 643 | 3 181 643 | 2 |
| SemiBin_3 | 99,46 | 2 952 929 | 2 | 2 231 097 | 2 231 097 | 0 |
| SemiBin_4 | 98,87 | 2 803 628 | 4 | 718 495 | 718 495 | 0 |
| SemiBin_5 | 89,72 | 6 063 346 | 17 | 560 681 | 431 947 | 2 |
| SemiBin_6 | 85,16 | 4 627 992 | 49 | 146 055 | 180 330 | 1 |
| SemiBin_7 | 76,52 | 6 594 866 | 70 | 137 626 | 178 696 | 1 |

A táblázatban szereplő Genome fraction értékeket minden esetben a kiválasztott referencia-genomhoz viszonyítva határoztuk meg.

A SemiBin_0–SemiBin_4 bin-ek magas referencia-lefedettséget mutattak. A SemiBin_5–SemiBin_7 esetében alacsonyabb genomfrakciót, több kontigot és kisebb N50-értékeket figyeltünk meg.

A misassembly-k száma önmagában nem elegendő a bin minőségének megítéléséhez. A nulla misassembly-t mutató bin is lehet töredékes vagy hiányos, ha a referencia kis hányadát fedi le. Hasonlóképpen, a magasabb misassembly-számot befolyásolhatja a referencia és a vizsgált szekvenciák közötti valódi szerkezeti különbség is.

## 6.5. Kiegészítő vizsgálat a két bizonytalan referenciára

A FastANI alapján a 'Phocaeicola dorei' és a 'Bacteroides thetaiotaomicron' referencia önálló visszanyerése nem volt egyértelmű. Ezért kiegészítő QUAST-vizsgálatban nemcsak a legjobb jelölt bin-t, hanem a teljes Flye-assembly-t is összehasonlítottuk ezekkel a referencia-genomokkal.

| Referencia | Vizsgált szekvenciák | Genome fraction (%) |
|---|---|---:|
| Phocaeicola dorei | SemiBin_0 | 53,87 |
| Phocaeicola dorei | Teljes Flye-assembly | 74,68 |
| Bacteroides thetaiotaomicron | SemiBin_7 | 34,04 |
| Bacteroides thetaiotaomicron | Teljes Flye-assembly | 42,20 |

Mindkét esetben a teljes assembly nagyobb referencia-lefedettséget adott, mint a megfelelő jelölt bin önmagában. Ez arra utal, hogy a bin-en kívüli kontigok között is találhatók a referenciához illeszkedő régiók.

Ezek az eredmények önmagukban nem bizonyítják, hogy a két genom megfelelően rekonstruálódott, vagy kizárólag a binning okozta a hiányt. A részleges egyezés eredhet az assembly töredezettségéből, a genomok közötti hasonlóságból vagy a binning során kialakuló szekvenciaeloszlásból is.

## 6.6. A FastANI- és QUAST-eredmények közös értelmezése

A FastANI- és QUAST-eredmények alapján több bin magas nukleotidszintű hasonlóságot és magas referencia-lefedettséget mutatott az eredeti referencia-genomjához képest.

A SemiBin_0, SemiBin_1, SemiBin_2, SemiBin_3 és SemiBin_4 bin-ek referencia-lefedettsége különösen magas volt. A SemiBin_5 és SemiBin_6 esetében szintén jelentős referencia-egyezést találtunk, de a lefedettség alacsonyabb volt, illetve a kontigok száma nagyobb. A SemiBin_7 töredezettebb bin volt, és a két bizonytalan referencia egyikéhez sem adott egyértelműen kölcsönös, megfelelően kiterjedt egyezést.

Az alkalmazott azonosíthatósági kritérium szerint a 10 eredeti referencia-genom közül 8-hoz tudtunk egyértelmű bin-jelöltet rendelni. Ez a projektben meghatározott szabály szerinti eredmény; nem jelenti automatikusan nyolc teljes, hibamentes vagy minden minőségi szabványnak megfelelő MAG igazolását.

## 6.8. Kapcsolat a következő munkafolyamattal

A referencia-alapú értékelés eredményeit a szimulációs metaadatokkal és a CheckM2-mutatókkal együtt egy tízsoros összesítő táblázatba rendeztük. Ebben minden sor egy eredeti referencia-genomnak felel meg, és külön oszlopokban szerepelnek a szimulációs paraméterek, a legjobb bin-azonosítás, a CheckM2- és QUAST-mutatók, valamint a referencia-genomok alapvető tulajdonságai.

Az összesítő táblázat generálását külön Python-script végzi, így az eredmények egységes formátumban állíthatók elő és később további binningeszközök eredményeivel is kiegészíthetők.
