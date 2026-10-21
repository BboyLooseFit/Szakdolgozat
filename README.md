###Metagenomikai workflow és MAG-visszanyerés

**1. A projekt háttere**

A szakdolgozat célja metagenomikai bioinformatikai munkafolyamatok vizsgálata és összehasonlítása szimulált, ismert összetételű mikrobiális közösségeken. A projekt során a szekvenciaadatok szimulálásától kezdve az assembly-n és binningen át a visszanyert metagenom-assembled genome-ok (MAG-ok) minőségértékeléséig követjük végig az elemzési folyamatot.

A szimuláció lehetővé teszi, hogy a visszanyert bin-eket a kiindulási közösség ismert referencia-genomjaihoz hasonlítsuk, és értékeljük a különböző bioinformatikai eszközök teljesítményét.

**2. Célkitűzések**

- ONT long-read szekvenciaadatok szimulálása ismert összetételű mikrobiális közösségből.
- De novo assembly készítése a szimulált readekből.
- A readek visszamappelése az assembly-re.
- A kontigok binningje és MAG-jelöltek előállítása.
- A visszanyert bin-ek teljességének és szennyezettségének értékelése.
- A visszanyert genomok összehasonlítása az eredeti referencia-genomokkal.
- További binningeszközök és paraméterkombinációk vizsgálata, majd eredményeik összehasonlítása.

**3. A pilot szimuláció**

Kiindulási adatok

A szimulációhoz a Human Microbiome Project (HMP) eredetű gut/SRS104912.tsv template-et használtuk, amely 720 referencia-genomot tartalmazott. A közösség komplexitásának csökkentésére a 10 legnagyobb abundanciájú genomot választottuk ki, majd a relatív abundanciaértékeket a kiválasztott részhalmazon belül normalizáltuk.

A szimulált közösség azonosítója top10 volt.

Főbb lépések

1. MeSS – adatszimulálás: egy ONT long-read minta generálása rögzített, 42-es seed értékkel.
2. Flye – de novo assembly: az eredeti readekből kontigok összeállítása.
3. minimap2 és samtools – read mapping: a readek visszaillesztése az assembly-re, majd rendezett és indexelt BAM-fájl előállítása.
4. SemiBin2 – binning: a kontigok csoportosítása a BAM-ból származó lefedettségi információ felhasználásával.
5. CheckM2 – minőségellenőrzés: a bin-ek teljességének és szennyezettségének becslése.
6. FastANI és QUAST – referencia-alapú értékelés: a visszanyert bin-ek összehasonlítása az eredeti referencia-genomokkal.
7. Összesítés: a szimulációs metaadatok és a minőségi mutatók egyesítése egy tízsoros táblázatban.

**4. A pilot workflow főbb eredményei**

- Szimulált adatok: 67 341 read, összesen 1 000 000 962 bp szekvencia, átlagosan körülbelül 14,85 kb-os readhossz.
- Flye assembly: 168 kontig, 37,94 Mb összhossz, 667 985 bp N50 és 44,65% GC-tartalom.
- Read mapping: 67 047 elsődlegesen leképezett read, ami 99,56%-os aránynak felel meg.
- SemiBin2: 8 bin; a 168 kontigból 150 került bin-be. A bin-ekbe rendelt kontigok az assembly bázisainak körülbelül 97,4%-át képviselték.
- CheckM2: a teljességi értékek 71,74–100,00%, a szennyezettségi értékek 0,50–10,01% között alakultak.
- FastANI és QUAST: a projekt-specifikus azonosíthatósági szabály alapján 8 referencia-genomhoz sikerült egyértelmű bin-jelöltet rendelni, míg 2 referencia különálló MAG-ként való visszanyerése nem volt igazolható.

Az azonosíthatósági szabály kölcsönös legjobb FastANI-találatot, legalább 95%-os ANI-t és legalább 50%-os referenciafragmentum-illeszkedési arányt követelt meg. Ez a projektben használt döntési szabály, nem univerzális MAG-minőségi szabvány.

A 8 azonosítható referencia nem jelenti automatikusan 8 teljes és hibamentes MAG visszanyerését. A teljességet, szennyezettséget és referencia-lefedettséget külön mutatókkal értékeltük.

Az elemzési lépések részletes leírása a docs/ könyvtárban található.

