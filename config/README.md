# Konfigurációs fájlok

## 1. A könyvtár célja

A `config/` könyvtár a metagenomikai munkafolyamatok futtatásához
szükséges konfigurációs és paraméterfájlokat tartalmazza.

A cél az, hogy a szimulációs és később az elemzési beállításokat
elkülönítsük a futtató szkriptektől. Így a különböző minták és
paraméterkombinációk egységes szerkezetben kezelhetők.

## 2. Szimulációs paraméterek

### `simulation_parameters.tsv`

Ez a fájl a MeSS segítségével végzett szekvenciaadat-szimuláció
paramétereit tartalmazza. Egy sor egy szimulációs mintának felel meg.

A fájl jelenleg az első pilot szimuláció beállításait tartalmazza.

### Oszlopok jelentése

| Oszlop | Jelentés |
|---|---|
| sample_id | A szimulált minta egyedi azonosítója |
| community_file | A szimulációhoz használt közösség bemeneti fájlja |
| technology | A szimulált szekvenálási technológia |
| model | A MeSS által használt szekvenálási hibamodell |
| mean_length_bp | A célzott átlagos readhossz bázispárban |
| accuracy | A MeSS számára megadott szekvenciapontosság |
| ratio | A MeSS futtatásakor megadott arányparaméter |
| seed | A véletlenszám-generátor rögzített kezdőértéke |
| save_bam | A BAM-kimenet előállításának engedélyezése |
| threads | A párhuzamos feldolgozáshoz használt szálak száma |
| output_dir | A szimuláció eredményeinek célkönyvtára |

### Jelenlegi pilot beállítások

| Paraméter | Érték |
|---|---|
| Sample ID | top10 |
| Technológia | ONT |
| Hibamodell | QSHMM-ONT-HQ |
| Célzott átlagos readhossz | 15 000 bp |
| Szekvenciapontosság | 0,95 |
| Ratio | 39:24:36 |
| Seed | 42 |
| BAM-kimenet | Engedélyezve |
| Szálak száma | 4 |

## 3. Fájlformátum és útvonalak

A `simulation_parameters.tsv` tabulátorral elválasztott szöveges
állomány (TSV). A fejléc tartalmazza az oszlopneveket, minden további
sor pedig egy-egy szimulációs mintát ír le.

A `community_file` és az `output_dir` mezők útvonalai az eredeti
MeSS projektkönyvtárához viszonyított relatív útvonalak.

A konfigurációs fájl a saját szakdolgozati repositoryban található,
az eredeti MeSS program és a szimulációs bemeneti adatok külön
projektkönyvtárban maradnak.

## 4. További fejlesztések

A későbbiekben további szimulációs minták és paraméterkombinációk
adhatók hozzá ugyanahhoz a táblázathoz.

A konfigurációs fájlt a tervezett futtató szkript fogja beolvasni,
és az egyes sorok alapján indítja majd el a MeSS szimulációkat.

A további elemzési lépésekhez külön konfigurációs beállítások is
létrehozhatók, például a binningeszközök és a minőségellenőrzési
folyamatok kezelésére.
