# 5. MAG-minőségellenőrzés a CheckM2 segítségével

## 5.1. A lépés célja

A SemiBin2 által létrehozott bin-ek minőségének értékeléséhez a CheckM2 programot használtuk.

A metagenomikai binning során létrejövő bin-ek nem feltétlenül reprezentálnak teljes, tiszta mikrobiális genomokat. Előfordulhat, hogy egyes genomrészek hiányoznak belőlük, vagy különböző mikroorganizmusokból származó kontigokat is tartalmaznak.

A CheckM2 marker-gén alapú megközelítéssel becsüli meg a bin-ek teljességét (*completeness*) és szennyezettségét (*contamination*). Ezek a mutatók segítik a visszanyert MAG-jelöltek minőségének értékelését.

## 5.2. Bemeneti adatok

A CheckM2 bemenetét a SemiBin2 által előállított nyolc bin jelentette.

A bin-ek FASTA-formátumú szekvenciaállományokat tartalmaztak.

## 5.3. A használt eszköz és a beállítások

- Eszköz: CheckM2
- Verzió: 1.1.0
- Bemenet: a SemiBin2 által előállított bin-ek
- Kimeneti könyvtár: `results/qc/checkm2/`
- Eredménytábla: `results/qc/checkm2/quality_report.tsv`
- Adatbázis: CheckM2 `uniref100.KO.1.dmnd`

## 5.4. A futtatás

A CheckM2 a `predict` művelettel értékeli a bemeneti genomállományokat. Az alkalmazott futtatás logikája a következő volt:

```bash
checkm2 predict \
  --input results/binning/semibin2/output_bins \
  --output-directory results/qc/checkm2 \
  --threads 1 \
  --database_path "$HOME/szakdolgozat/checkm2-db/CheckM2_database/uniref100.KO.1.dmnd"
```

A '--threads 1' beállítást a korábbi memóriahiány miatt alkalmaztuk.

## 5.5. A CheckM2 által szolgáltatott mutatók

### Completeness – teljesség

A teljesség azt becsüli meg, hogy a vizsgált bin-ben a várt marker-gének alapján mennyire teljes a genom reprezentációja.

A magasabb érték általában azt jelzi, hogy a bin több, a vizsgált taxonómiai szinten várt markert tartalmaz.

### Contamination – szennyezettség

A szennyezettségi érték a marker-gének mintázata alapján becsüli a bin-ben jelen lévő, várhatóan nem megfelelően reprezentált vagy többszörösen előforduló genetikai anyag mértékét.

A magasabb érték kevert eredetű kontigok jelenlétére vagy más, a bin tisztaságával kapcsolatos problémára utalhat.

## 5.6. Eredmények

A CheckM2 sikeresen lefutott mind a nyolc bin-re. Az eredmények a következő táblázatban foglalhatók össze.

| Bin | Completeness (%) | Contamination (%) |
|---|---:|---:|
| SemiBin_0 | 100,00 | 1,98 |
| SemiBin_1 | 99,96 | 0,50 |
| SemiBin_2 | 96,34 | 5,76 |
| SemiBin_3 | 94,20 | 2,26 |
| SemiBin_4 | 82,83 | 1,29 |
| SemiBin_5 | 78,80 | 7,22 |
| SemiBin_6 | 71,77 | 4,58 |
| SemiBin_7 | 71,74 | 10,01 |

A completeness értéke 71,74–100,00% között, a contamination értéke 0,50–10,01% között változott.

A 'SemiBin_0' és a 'SemiBin_1' mutatta a legmagasabb teljességet, alacsony szennyezettség mellett.

A 'SemiBin_5', 'SemiBin_6' és 'SemiBin_7' teljessége alacsonyabb volt. A 'SemiBin_7' mutatta a legmagasabb szennyezettségi értéket a nyolc bin közül.

## 5.7. Az eredmények értelmezése és korlátai

A CheckM2 eredményei alapján a nyolc bin eltérő teljességű és szennyezettségű volt. A teljesség és a szennyezettség együtt hasznos a MAG-jelöltek minőségének értékeléséhez, de egyik mutató sem bizonyítja önmagában, hogy egy bin egyetlen eredeti referencia-genomnak felel meg.

A CheckM2 nem a szimulációban használt referencia-genommal való közvetlen azonosságot méri. Ezért a kapott értékeket a későbbi FastANI- és QUAST-eredményekkel együtt értékeltük.

## 5.8. Kapcsolat a következő munkafolyamattal

A CheckM2-vizsgálat után a bin-eket a szimuláció során használt, ismert referencia-genomokkal is összehasonlítottuk.

A FastANI segítségével a nukleotidszintű hasonlóságot és a referencia-fragmensek illeszkedését vizsgáltuk. A QUAST-tal pedig a referencia mellett értékeltük az assembly hosszát, a genomfrakciót, a kontigszámot, az N50-et és a szerkezeti eltéréseket.
