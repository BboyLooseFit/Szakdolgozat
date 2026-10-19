# Szakdolgozat
Szakdolgozattal kapcsolatos információk, lépések leírása és a hozzájuk tartozó kódok.

ELSŐ CÉLKITŰZÉS
Egyetlen, ismert összetételű közösségből szimulált ONT-minta előállítása, majd assembly, binning és a visszanyert MAG-jelöltek minőségértékelése.

A szimulált közösségek egy Human Microbiome Project eredetű bélmikrobiom-template alapján készültek. Az azonos abundanciaeloszlás megőrzése mellett a közösség komplexitását a 10 legabundánsabb genom kiválasztásával módosítottuk.

A pilot szimulációhoz a gut/SRS104912 template került kiválasztásra, amely 720 referencia genomot tartalmaz.

ELSŐ WORKFLOW
1. MeSS: 10 referencia-genom relatív abundanciák alapján, egy szimulált minta, rögzített seed=42.
2. Szimulált adatok: kb. 1,0 Gb, 67 341 read, átlagos readhossz kb. 14,85 kb.
3. Flye assembly: 168 kontig, összesen 37,94 Mb, N50 = 667 985 bp, GC = 44,65%.
4. Read mapping: a readek kb. 99,56%-a leképeződött az assembly-re az elsődleges illesztések alapján.
5. SemiBin2: 8 bin; 168 kontigból 150 került bin-be. A bin-ekhez rendelt kontigok az assembly bázisainak kb. 97,4%-át fedik.
6. CheckM2: teljesség 71,74–100%; szennyezettség 0,50–10,01%.
7. FastANI és QUAST: a bin-eket a szimulációban használt referencia-genomokhoz hasonlítottuk.


ELSŐ FŐ EREDMÉNYEK
- A SemiBin_0–SemiBin_4 fő referencia-összehasonlításai magas, 97,2–99,9%-os genomfrakciót mutattak.
- A SemiBin_5–SemiBin_7 töredezettebb bin-ek, alacsonyabb referencia-lefedettséggel és N50-értékekkel.
- A 10 eredeti referencia mellett 8 bin keletkezett. Ezért nem állítható, hogy mind a 10 genom különálló MAG-ként sikeresen visszanyerésre került.
- A GCF_013009555.1 és GCF_024760065.1 referenciák különálló MAG-ként történő visszanyerése nem igazolt. A teljes assembly QUAST-genomfrakciói rendre kb. 74,7% és 42,2%.
