# 7. A tíz referencia-genomot összefoglaló táblázat előállítása

## 7.1. A lépés célja

A korábbi munkafolyamatok eredményeit egy közös, tízsoros táblázatban foglaltuk össze. A táblázat minden sora egy, a MeSS-szimulációhoz kiválasztott eredeti referencia-genomot reprezentál.

A cél az volt, hogy együttesen áttekinthetők legyenek a szimuláció bemeneti adatai, a referencia-genomok tulajdonságai, a szimulált readek száma, a legjobb bin-hozzárendelések, valamint a visszanyert bin-ek minőségi mutatói.

Az összesítés a korábban előállított MeSS-, FastANI-, CheckM2- és QUAST-eredményeket kapcsolja össze.

## 7.2. A táblázat felépítése

A táblázat tíz sora a tíz eredeti referencia-genomnak felel meg. Az oszlopok több információcsoportot foglalnak össze.

### Azonosító és taxonómiai adatok

- Mintaazonosító.
- Baktériumnév.
- Referencia-genom accession, például GCF-azonosító.

### Szimulációs adatok

- A szimulációhoz használt relatív abundancia.
- A MeSS kimeneti metaadataiban szereplő becsült readszám.
- A BAM-fájlból összesített, az adott referencia-genomhoz rendelt readek száma.
- A teljes mintára vonatkozó readszám.
- A MeSS által rögzített szimulációs paraméterek.

A becsült readszámot és a BAM-ból számított readszámot külön oszlopban tartjuk. Ezek eltérő eredetű mutatók, ezért nem helyettesítik egymást.

### Binning és azonosíthatóság

- A használt binningeszköz.
- A legjobb jelölt bin neve.
- FastANI-érték.
- A referencia-fragmensek FastANI-alapú illeszkedési aránya.
- Az `identifiable_0_1` változó.

Az `identifiable_0_1` értéke a projektben alkalmazott operatív szabály alapján készül. Az azonosíthatóság feltétele a kölcsönös legjobb bin–referencia találat, legalább 95%-os ANI és legalább 50%-os referenciafragmentum-illeszkedés.

Ez projekt-specifikus döntési szabály, nem általánosan elfogadott MAG-minőségi küszöb.

### MAG-minőségellenőrzési adatok

A táblázat a CheckM2 által becsült teljességet (*completeness*) és szennyezettséget (*contamination*) is tartalmazza.

A QUAST eredményeiből szerepel benne a referencia-lefedettség (*genome fraction*), a kontigok száma, a teljes assembly-hossz, az N50 és az NG50.

A CheckM2 'strain heterogeneity' értékét a jelenlegi kimenet nem tartalmazta. Ezért ez a mező 'NA' értéket kapott, amely azt jelenti, hogy a paramétert ebben a munkafolyamatban nem határoztuk meg.

### A referencia-genom tulajdonságai

A táblázat rögzíti az eredeti referencia-genomok méretét és GC-tartalmát is. Ezek a jellemzők az assembly metadata-fájlokból származnak, nem a visszanyert bin-ekből.

## 7.3. Adatforrások

Az összesítő táblázat több korábban előállított adatforrásból készül.

| Adatforrás | Felhasználás |
|---|---|
| `communities/top10_community.tsv` | A szimuláció közösségének referencia-azonosítói és relatív abundanciái |
| `results/one_sample/coverages.tsv` | A MeSS által megadott becsült readszámok és szimulációs metaadatok |
| `results/one_sample/assembly_finder/taxonomy.tsv` | Taxonómiai adatok |
| `results/one_sample/assembly_finder/assembly_summary.tsv` | Referencia-genomok mérete és GC-tartalma |
| `results/qc/fastani/bin_vs_reference.tsv` | Bin–referencia összehasonlítás |
| `results/qc/fastani/reference_vs_bins.tsv` | Fordított referencia–bin összehasonlítás |
| `results/qc/fastani/simulated_reads_by_reference.tsv` | A BAM-ból összesített, referencia-genomhoz rendelt readszámok |
| `results/qc/checkm2/quality_report.tsv` | CheckM2 teljességi és szennyezettségi értékek |
| `results/qc/quast/reports/` | QUAST referencia-alapú assembly-mutatók |

A táblázat létrehozásakor az egyes forrásokból származó adatokat a referencia-genom accession azonosítója, illetve a megfelelő bin-hozzárendelés alapján kapcsoltuk össze.

## 7.4. A táblázat generálása

Az összesítő táblázatot egy Python-script állítja elő:

`scripts/build_10_reference_summary.py`

A script a meglévő metaadatokat és minőségellenőrzési jelentéseket beolvassa, majd elkészíti a tíz referencia-genomot összefoglaló TSV-fájlt.

A generált fájl helye:

`results/qc/10_reference_summary.tsv`

A script jelenlegi változata a pilot elemzésben használt könyvtárstruktúrát és fájlokat várja. A reprodukálható használathoz a szükséges bemeneti fájlokat és szoftverkörnyezetet is biztosítani kell.

## 7.5. Az eredmények értelmezése

Az összesítő táblázat lehetővé teszi, hogy a szimulációs közösség mind a tíz eredeti referencia-genomját külön-külön értékeljük.

Az `identifiable_0_1` mező a választott FastANI-alapú szabály szerinti azonosíthatóságot jelzi. A `1` érték nem jelenti automatikusan, hogy a megfelelő bin teljes, hibamentes vagy minden minőségi követelménynek megfelelő MAG.

A CheckM2-mutatók a marker-gének alapján becsülik a bin teljességét és szennyezettségét, míg a QUAST a kiválasztott referencia-genomhoz viszonyítja az assembly szerkezeti és lefedettségi tulajdonságait.
