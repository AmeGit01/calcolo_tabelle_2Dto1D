# calcolo_tabelle_2Dto1D

Script in Julia che trasforma i risultati di una simulazione **BASEMENT 2D** in **tabelle
geometrico-idrauliche per sezione**, da usare come input di un modello 1D (moto vario, sezioni
tabellate).

L'idea è semplice. Nel modello 2D si fa crescere la portata a **scalini**, raggiungendo il
moto permanente a ogni scalino, e lungo l'asta si registrano le grandezze integrali su una
serie di **nodestring** (le sezioni). Ogni scalino dà, per ogni sezione, un punto
"livello → area, larghezza, …". Mettendo in fila gli scalini si ottiene, per ogni sezione, una
tabella in funzione del livello idrometrico (WSE) che il modello 1D può interpolare.

## Gli script

Gli script sono in [src/](src/) e sono indipendenti: ognuno si lancia da terminale con
`julia --project=. src/<script>.jl [<cartella>]` e scrive solo nella cartella `output/<cartella>`.

| Script | Legge | Scrive | Scopo |
|---|---|---|---|
| [main_sezioni.jl](src/main_sezioni.jl) | `input/<caso>` | `tabella_sezioni.csv`, `localslope_xcoord.csv`, `topologia_sezioni.csv` | calcola le tabelle delle sezioni per il modello 1D |
| [main_diversione.jl](src/main_diversione.jl) | `input/<caso>` | `scala_diversione.csv` | scala di deflusso portata-portata della diversione Brenta nuovo / Brenta vecchio |
| [grafici.jl](src/grafici.jl) | `output/<caso>` | `grafici/` | grafici di controllo delle tabelle già calcolate |
| [W_check.jl](src/W_check.jl) | `output/<caso>` | `W_check/` | cerca i picchi della larghezza W che creano problemi nel modello 1D |

I primi due lavorano sui risultati 2D, gli altri due **solo sugli output già salvati**, quindi si
possono rilanciare in ogni momento senza rifare il calcolo. Comandi pronti da copiare, dalla
cartella del progetto, nell'ordine in cui vanno lanciati per il caso `Brenta_caldonazzo_borgo`:

```bash
julia --project=. src/main_sezioni.jl Brenta_caldonazzo_borgo
julia --project=. src/main_diversione.jl Brenta_caldonazzo_borgo
julia --project=. src/grafici.jl Brenta_caldonazzo_borgo
julia --project=. src/W_check.jl Brenta_caldonazzo_borgo
```

Tutti e quattro insieme:

```bash
for s in main_sezioni main_diversione grafici W_check; do julia --project=. src/$s.jl Brenta_caldonazzo_borgo || break; done
```

Gli ultimi tre hanno `Brenta_caldonazzo_borgo` come cartella predefinita e si possono lanciare
anche senza argomento; `main_sezioni.jl` no (il suo caso predefinito è `S_0050_6000_q2`). Il
calcolo delle tabelle richiede alcuni minuti, gli altri script pochi secondi più la compilazione
di Julia.

---

## Indice

1. [Calcolo delle tabelle (`main_sezioni.jl`)](#calcolo-delle-tabelle-main_sezionijl)
2. [Grandezze tabellate](#grandezze-tabellate)
3. [Input](#input)
4. [Output](#output)
5. [Come lanciarlo](#come-lanciarlo)
6. [Parametri di configurazione](#parametri-di-configurazione)
7. [Controlli e avvisi](#controlli-e-avvisi)
8. [Diversione Brenta nuovo / Brenta vecchio (`main_diversione.jl`)](#diversione-brenta-nuovo--brenta-vecchio-main_diversionejl)
9. [Grafici di controllo (`grafici.jl`)](#grafici-di-controllo-graficijl)
10. [Controllo della larghezza W (`W_check.jl`)](#controllo-della-larghezza-w-w_checkjl)
11. [Struttura della cartella](#struttura-della-cartella)
12. [Note e limitazioni](#note-e-limitazioni)

---

## Calcolo delle tabelle (`main_sezioni.jl`)

Il calcolo delle tabelle è in [src/main_sezioni.jl](src/main_sezioni.jl). In ordine:

1. **Controllo dei file di input.** Verifica che esistano `centerpoints.csv`, `bankpoints.csv`
   e un file per ogni sezione elencata in `centerpoints.csv` (più mesh e risultati 2D se β è
   calcolato dal 2D). Se manca qualcosa elenca tutti i file mancanti e si ferma. Fanno eccezione i
   file delle sezioni con `SaltaSezioniSenzaFile = true`: le sezioni senza file sono escluse con un
   avviso giallo che le elenca (in `Brenta_caldonazzo_borgo`, per esempio, B17700–B33800 e V00000,
   che la simulazione non ha salvato).
2. **Riconoscimento delle aste.** Le aste si ricavano dal nome delle sezioni, che è
   `<asta><numero>`: per esempio `N00050` → asta `N`, `B17650` → asta `B`, `S0400` → asta `S`.
   Ogni asta è elaborata come un **tratto indipendente**, con la sua ascissa `s_m` e le sue
   tabelle (vedi [Più aste](#più-aste)).
3. **Lettura delle sezioni** di ogni asta, ordinate da monte verso valle secondo l'ascissa
   curvilinea `s_m`. Tutte le sezioni devono avere lo stesso numero di righe (gli stessi scalini)
   e, se c'è la colonna del tempo, gli stessi istanti.
4. **Pulizia delle righe asciutte.** BASEMENT scrive 0 in tutte le colonne quando la sezione è
   asciutta. Quote di riferimento, fondo e WSE di quelle righe vengono riportate alla quota di
   riferimento della sezione.
5. **Pendenze**, per ogni sezione e per ogni scalino, dal profilo lungo l'asta allo stesso
   scalino (vedi [Pendenza dell'energia](#pendenza-dellenergia)):
   - `iWSE`, pendenza della superficie libera, solo informativa;
   - `J`, **pendenza dell'energia**, usata per Ks e jden;
   - entrambe alle differenze centrate se le due sezioni adiacenti sono bagnate, con la
     differenza unilaterale agli estremi o se una delle adiacenti è asciutta;
   - le J nulle o negative sono interpolate in WSE fra i livelli validi della stessa sezione, e
     ogni sostituzione viene segnalata.
6. **Differenza di portata ΔQ** con la sezione a monte, allo stesso scalino. Rappresenta lo
   scambio con la piana (esondazione/rientro) nel tratto fra le due sezioni.
7. **Coefficiente di quantità di moto β**: vale 1 oppure è calcolato dal campo di velocità 2D
   (vedi [β dal modello 2D](#β-dal-modello-2d)).
8. **Ricostruzione delle righe asciutte.** I primi scalini possono essere asciutti, perché il 2D
   non è stabile a portate troppo basse. La prima riga è ricostruita con tirante `h_dry`,
   portata `Q_dry` e area ricavata da Manning imponendo `Ks_dry`. Le eventuali righe successive
   sono interpolate linearmente fino al primo scalino bagnato. Il modello 1D è pensato per le
   piene, quindi a queste righe basta essere ragionevoli.
9. **Livelli monotoni.** Il modello 1D vuole WSE strettamente crescenti: le righe con livello non
   superiore al massimo dei precedenti vengono scartate e segnalate.
10. **Ultima sezione (condizione di uscita).** Per l'ultima sezione di ogni asta `J` è
   sostituita dalla pendenza media del fondo sull'intera asta. Così il moto uniforme che il modello 1D impone in
   uscita restituisce la portata degli scalini 2D. Vicino al bordo di valle del dominio 2D la
   pendenza del pelo libero è falsata dal richiamo dell'uscita.
11. **Calcolo della tabella** di ogni sezione (vedi sotto) e controlli di plausibilità.
12. **Scrittura** della tabella unica per il modello 1D, del file con la geometria longitudinale e
    di quello con la topologia (numero progressivo → nome della sezione). Con
    `CartelleSezioni = true` scrive anche una cartella per sezione con la tabella completa e i
    grafici. Con più aste e `AsteUnite = true` i tre file sono unici per tutte le aste (vedi
    [Più aste](#più-aste)).

---

## Grandezze tabellate

Ogni tabella ha **una riga per scalino di portata** (livello), con il livello crescente.
Le colonne con suffisso `_s` sono la **derivata rispetto al livello**, d(·)/dWSE, calcolata alle
differenze in avanti (l'ultima riga ripete la penultima).

| Colonna | Simbolo | Unità | Descrizione |
|---|---|---|---|
| `WSE` | η | m s.l.m. | Livello medio della superficie libera (*Mean wse*) |
| `Ω` | Ω | m² | Area bagnata (*Wetted area*) |
| `Ω_s` | dΩ/dη | m | Derivata dell'area rispetto al livello |
| `β` | β | – | Coefficiente di quantità di moto (Boussinesq) |
| `β_s` | dβ/dη | 1/m | Derivata di β |
| `W` | W | m | Larghezza bagnata (*Wetted geometric length*, contorno bagnato) |
| `W_s` | dW/dη | – | Derivata di W |
| `RH` | R_H = Ω/W | m | Raggio idraulico |
| `RH_s` | dR_H/dη | – | Derivata di R_H |
| `iWSE` | i | – | Pendenza della superficie libera (positiva se il livello scende verso valle); solo informativa |
| `J` | J | – | Pendenza dell'energia (positiva se l'energia scende verso valle), usata per Ks e jden |
| `Ks` | K_s | m^(1/3)/s | Coefficiente di Gauckler-Strickler equivalente dell'intera sezione |
| `jden` | – | m | Denominatore del termine d'attrito del modello 1D |
| `jden_s` | – | – | Derivata di `jden` |
| `ΔQ` | ΔQ | m³/s | Portata della sezione meno quella della sezione a monte (negativa se l'acqua esce dall'alveo) |
| `ΔQ_s` | – | m²/s | Derivata di ΔQ |

Formule usate:

- Ks equivalente, da Manning con la portata locale della sezione e la pendenza dell'energia:

  K_s = Q / ( Ω · R_H^(2/3) · √J )

- Denominatore dell'attrito, tale che il termine d'attrito del modello 1D sia j = u·|u| / (g·jden):

  jden = K_s² · R_H^(4/3) / g = u² / (g·J),  con u = Q/Ω

  Allo stato del 2D, quindi, il modello 1D ha come pendenza d'attrito proprio J.

- ΔQ vale 0 nella prima sezione e quando una delle due sezioni è asciutta.

### Pendenza dell'energia

Ks e jden si ricavano dalla pendenza dell'energia

J = −d/ds ( WSE + u²/2g )

e non da quella del pelo libero `iWSE`. Manning/Strickler lega l'attrito alla dissipazione di
energia, e le due pendenze coincidono solo se la velocità non cambia lungo l'asta. Dove la sezione
si allarga la corrente rallenta e il pelo libero può risalire anche se l'energia scende: lì
`iWSE` è nulla o negativa e non darebbe un Ks. Il termine cinetico usa α = 1 (velocità uniforme
sulla sezione).

J può comunque risultare nulla o negativa in pochi punti, per esempio dove la velocità è molto
non uniforme sulla sezione o dove l'acqua esce o rientra lateralmente. In quei punti Ks e jden
non sarebbero definiti, quindi J viene sostituita così:
1. interpolata linearmente in WSE fra il livello valido (bagnato, con J > 0) più vicino sotto e
   quello più vicino sopra, nella stessa sezione;
2. se c'è un solo livello valido, sotto o sopra, si prende quello;
3. se la sezione non ne ha nessuno, si usa la pendenza media del fondo dell'asta.

jden si calcola poi con la velocità del livello, quindi resta coerente con la sua portata. Ogni
sostituzione compare fra gli avvisi, con la sezione, la riga e il valore originale.

### β dal modello 2D

Se `β_dal_2D = true`, β è calcolato come nel metodo di Engelund

β = Ω / Q² · ∫ u² Y dy

ma usando le velocità del 2D al posto di quelle di moto uniforme. Le "strisce" sono i triangoli
della mesh appena **a monte** di ogni lato della nodestring. Per ciascuno si prende il tirante
Y = WSE − z_b, la portata specifica normale alla sezione q_n = **q**·**n** e la lunghezza del
lato ΔL, da cui:

Ω = Σ Y·ΔL,  Q = Σ q_n·ΔL,  β = Ω · Σ(q_n²/Y · ΔL) / Q²

Gli elementi con tirante non superiore a `MinWaterDepth` del modello 2D sono considerati
asciutti. Lo script verifica che l'area ricostruita dagli elementi coincida con quella del file
della sezione e segnala gli scalini in cui la portata degli elementi si discosta di più di
`Qrel_lim` da quella della sezione.

Con `β_dal_2D = false` si pone β = 1: mesh e risultati HDF5 non servono e non vengono letti.

---

## Input

Gli input di un caso stanno in una sottocartella di `input/`. La cartella da usare si passa
dalla riga di comando, vedi [Come lanciarlo](#come-lanciarlo).

```
input/<caso>/
├── centerpoints.csv                              # elenco delle sezioni da elaborare
├── bankpoints.csv                                # quote degli argini
├── <prefisso><nome>.csv                          # un file per sezione (output BASEMENT)
├── <qualsiasi nome>.2dm                          # mesh 2D      (solo se β_dal_2D)
└── <qualsiasi nome>results.h5                    # risultati 2D (solo se β_dal_2D)
```

I nomi dei file non vanno scritti nello script, perché vengono riconosciuti da soli:
- il **prefisso** dei file delle sezioni (es. `ns_50_scalini_`) si ricava dai file che finiscono
  con il nome di una sezione di `centerpoints.csv`, e deve essere uguale per tutte le sezioni;
- **mesh e risultati** sono l'unico file `*.2dm` e l'unico file `*results.h5` della cartella.
  Se ce ne sono zero o più di uno, lo script si ferma.

In tutti i file CSV le colonne sono cercate **per nome**, non per posizione. Il nome deve essere
esatto, maiuscole comprese. Le colonne in più, per esempio altri attributi esportati da QGIS,
vengono ignorate e l'ordine delle colonne non conta. Se manca una colonna richiesta, lo script si
ferma e stampa l'elenco delle colonne mancanti, di quelle trovate e un esempio di file valido.

### `centerpoints.csv`: elenco delle sezioni

Un punto per sezione sull'asse del canale (tipicamente esportato da QGIS). **Solo le sezioni
elencate qui vengono elaborate**, nell'ordine di `s_m` (l'ordine delle righe non conta).

| Colonna | Tipo | Descrizione |
|---|---|---|
| `fid` | intero | Identificativo GIS (non usato) |
| `nome` | testo | Nome della sezione, `<asta><numero>`, es. `S0400` o `N00050`. Deve coincidere con il nome della nodestring |
| `s_m` | reale | Ascissa curvilinea lungo l'asse della propria asta [m], crescente verso valle |
| `dist_next` | reale | Distanza dalla sezione successiva [m] (non usato) |

### `bankpoints.csv`: quote degli argini

Due punti per sezione, sulla sommità degli argini.

| Colonna | Tipo | Descrizione |
|---|---|---|
| `fid` | intero | Identificativo GIS (non usato) |
| `nome` | testo | `<nome sezione>_sx` o `<nome sezione>_dx` (sinistra/destra idrografica) |
| `s_m` | reale | Ascissa curvilinea [m] (non usata) |
| `z_dtm` | reale | Quota dell'argine dal DTM [m s.l.m.] |
| `dist_cp` | reale | Distanza dal centerpoint [m] (non usata) |

Ogni sezione di `centerpoints.csv` deve avere sia `_sx` che `_dx`.

### File delle sezioni: `<prefisso><nome>.csv`

È l'output per nodestring di BASEMENT 2D (*NodeStrg* → CSV), un file per sezione, ad esempio
`ns_50_scalini_S0400.csv`. **Ogni riga è un'uscita della simulazione**: la prima è la
condizione iniziale (tutta a zero), le successive sono gli scalini di portata a moto permanente.

Colonne usate (i nomi devono essere esattamente questi):

| Colonna | Uso |
|---|---|
| `Mean wse [m]` | livello WSE |
| `Discharge [m3/s]` | portata Q |
| `Wetted area [m2]` | area Ω; se è 0 la riga è considerata asciutta |
| `Reference elevation [m]` | quota di riferimento della sezione (talweg), letta dalla prima riga bagnata |
| `Mean bottom elevation [m]` | riparata nelle righe asciutte |
| `Wetted geometric length [m]` | larghezza/contorno bagnato W |
| `t [s]` | facoltativa; se c'è, tutte le sezioni devono avere gli stessi istanti |

Le altre colonne (volumi, trasporto solido, …) vengono ignorate.

Requisiti:
- tutte le sezioni hanno lo **stesso numero di righe** (gli stessi scalini);
- una volta bagnata, una sezione deve restare bagnata negli scalini successivi;
- una sezione deve essere bagnata almeno in una riga.

### Mesh `.2dm` e risultati `.h5` (solo con `β_dal_2D = true`)

- **Mesh** in formato SMS 2dm: nodi `ND`, triangoli `E3T` numerati in ordine e nodestring
  `NS n1 n2 … -nk nome`. Ogni sezione di `centerpoints.csv` deve essere una nodestring con lo
  stesso nome.
- **Risultati** HDF5 di BASEMENT 3.x, con i gruppi:
  - `CellsAll/Topology`, `CellsAll/BottomEl`: devono corrispondere alla mesh;
  - `Parameters/MinWaterDepth`;
  - `.config/model`: da qui si legge l'`upstream_direction` di ogni nodestring;
  - `NodeStrg/StrNames`, `RESULTS/NodeStrg/StateVar/*`: usati per verificare che i CSV delle
    sezioni corrispondano alle uscite del file;
  - `RESULTS/CellsAll/HydState/*`: livello, q_x, q_y per cella.

  Il numero di uscite deve essere uguale al numero di righe dei CSV delle sezioni.

### Casi presenti in `input/`

| Cartella | Sezioni | Scalini | Note |
|---|---|---|---|
| `Brenta_caldonazzo_borgo` | 738 in `centerpoints.csv` su due aste, B (B00000–B33800) e V (V00000–V03000), ogni 50 m; **414 con risultati** (B00000–B17650, V00050–V03000), prefisso `scalini_` | 21 + riga iniziale | caso attuale; contiene la diversione Brenta nuovo / Brenta vecchio; niente mesh/`.h5` |
| `S_0050_6000_q2` | 120, ogni 50 m (S0050–S6000), prefisso `ns_50_scalini_` | 52 + riga iniziale, scalini di ~2 m³/s | caso predefinito di `main_sezioni.jl`; c'è la mesh ma non il `.h5` (serve solo con β dal 2D) |
| `Brenta_minimum` | 735 su tre aste: N (321), V (60), B (354), ogni ~50 m (prefisso `minimum_scalini_`) | 21 + riga iniziale | niente mesh/`.h5` |

I file delle sezioni di `Brenta_caldonazzo_borgo` sono estratti da
`BM_results/inputs/scalini_results.h5` con lo script Python `BM_results/py_src/BMv41nsalldata.py`.

### Più aste

Se i nomi delle sezioni hanno prefissi diversi (es. `N00050`, `V00050`, `B00000`), ogni prefisso è
un'asta e viene elaborata **indipendentemente** dalle altre:
- ha la sua ascissa `s_m`, misurata lungo la propria asta, che può ripartire da 0;
- le pendenze, ΔQ e la condizione di uscita si calcolano solo con le sezioni della stessa asta;
- la prima sezione di ogni asta ha ΔQ = 0 e l'ultima è trattata come sezione di uscita.

Come vengono scritti i risultati dipende da `AsteUnite`:
- `true` (predefinito): **un solo** `tabella_sezioni.csv`, `localslope_xcoord.csv` e
  `topologia_sezioni.csv` per tutte le aste, con le aste una dopo l'altra in ordine alfabetico
  (B prima di V) e il numero progressivo che continua da un'asta all'altra. In
  `Brenta_caldonazzo_borgo`: B = sezioni 1–354, V = sezioni 355–414. `x_coord_[m]` riparte invece
  dall'inizio di ogni asta, e `topologia_sezioni.csv` dice a quale asta appartiene ogni numero;
- `false`: una sottocartella per asta, `output/<caso>/<asta>/`, con il numero progressivo che
  riparte da 1 in ogni asta.

Lo script non sa come le aste sono collegate fra loro (confluenze, diramazioni): i collegamenti
vanno gestiti nel modello 1D. Un'asta deve avere almeno due sezioni.

---

## Output

Gli output vanno in `output/<caso>/` (cartella `OutputFolder`; con β dal 2D il nome ha il
suffisso `_beta`). Contenuto completo, con i file di tutti gli script:

```
output/<caso>/
├── tabella_sezioni.csv        # main_sezioni.jl: tabella unica per il modello 1D
├── localslope_xcoord.csv      # main_sezioni.jl: geometria longitudinale
├── topologia_sezioni.csv      # main_sezioni.jl: numero progressivo → nome e asta della sezione
├── scala_diversione.csv       # main_diversione.jl
├── grafici/                   # grafici.jl
├── W_check/                   # W_check.jl
└── <nome sezione>/            # main_sezioni.jl, solo con CartelleSezioni = true
    ├── tabella.csv            # tabella completa della sezione
    ├── area_bagnata.svg
    ├── larghezza_bagnata.svg
    ├── raggio.idraulico.svg
    ├── pendenza_wse.svg
    ├── pendenza_energia.svg
    ├── ks.svg
    ├── jden.svg
    └── beta.svg               # solo con β_dal_2D = true
```

Con più aste e `AsteUnite = false` i file di `main_sezioni.jl` stanno invece in
`output/<caso>/<asta>/` (vedi [Più aste](#più-aste)). Gli altri script leggono solo la struttura
con le aste unite.

### `tabella_sezioni.csv`: input del modello 1D

Tutte le sezioni impilate, da monte verso valle, nel formato letto da `READ_VALUES!` del
modello 1D:

```
sezione,WSE,Ω,Ω_s,β,β_s,jden,jden_s,W,W_s,RH,RH_s,ΔQ,ΔQ_s
1,430.4703...,...
1,430.7386...,...
...
120,...
```

- `sezione` è il numero progressivo (1 = più a monte), non il nome;
- il modello legge **solo le prime 12 colonne**, in questo ordine fisso;
- `ΔQ` e `ΔQ_s` (colonne 13–14) per ora non sono usate dal modello;
- `iWSE`, `J` e `Ks` non compaiono qui: sono nelle tabelle delle singole sezioni.

### `localslope_xcoord.csv`: geometria longitudinale

Una riga per sezione:

| Colonna | Descrizione |
|---|---|
| `section` | numero progressivo, come in `tabella_sezioni.csv` |
| `local_slope` | pendenza locale del fondo, **negativa se il fondo scende verso valle** (differenze centrate sul talweg). Nell'ultima sezione è la pendenza media dell'intera asta |
| `x_coord_[m]` | ascissa curvilinea `s_m` |
| `talweg_[m]` | quota di riferimento della sezione (`Reference elevation`) |
| `argine_dx_[m]`, `argine_sx_[m]` | quote degli argini da `bankpoints.csv` |

### `topologia_sezioni.csv`: nomi delle sezioni

Il modello 1D numera le sezioni con un intero progressivo; questo file dice a quale sezione reale
corrisponde ogni numero. Una riga per sezione:

| Colonna | Descrizione |
|---|---|
| `section` | numero progressivo, come in `tabella_sezioni.csv` e `localslope_xcoord.csv` |
| `nome` | nome della sezione, es. `B13200` |
| `asta` | asta, ricavata dal prefisso del nome (`B`, `V`, …) |
| `x_coord_[m]` | ascissa curvilinea lungo la propria asta |

### `<nome sezione>/tabella.csv` (solo con `CartelleSezioni = true`)

Tutte le colonne della [tabella delle grandezze](#grandezze-tabellate), compresi `iWSE`, `J` e `Ks`.
Tabella e grafici di una sezione vengono **riscritti solo se la tabella è cambiata** rispetto a
quella già salvata, così i file non cambiano inutilmente.

### Output a terminale

Una riga di riepilogo per sezione:

```
Sezione     n WSE             J [-]    Ks [m^(1/3)/s]             β [-]
S0050          53  0.00xxx  0.00xxx    xx.xx    xx.xx    1.000    1.000
...
```

(numero di livelli; min e max di J, Ks e β), seguita dall'elenco degli
[avvisi](#controlli-e-avvisi).

### Output presenti in `output/`

| Cartella | Origine |
|---|---|
| `Brenta_caldonazzo_borgo` | caso `Brenta_caldonazzo_borgo`, β = 1, aste unite, senza cartelle per sezione; contiene anche gli output di `main_diversione.jl`, `grafici.jl` e `W_check.jl` |
| `S_50m_q2` | caso `S_0050_6000_q2`, β = 1, con una cartella per sezione |

---

## Come lanciarlo

Serve **Julia ≥ 1.12** (il `Manifest.toml` è generato con Julia 1.12.6). Le dipendenze sono in
[Project.toml](Project.toml): CSV, DataFrames, HDF5, Plots, Printf, Infiltrator.

La prima volta, dalla cartella del progetto:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Comando pronto da copiare per il caso attuale:

```bash
julia --project=. src/main_sezioni.jl Brenta_caldonazzo_borgo
```

In generale il caso si sceglie dalla riga di comando, senza modificare lo script:

```bash
julia --project=. src/main_sezioni.jl <cartella di input> [<cartella di output>]
```

| Argomento | Significato |
|---|---|
| `<cartella di input>` | sottocartella di `input/`. Va bene sia `Brenta_caldonazzo_borgo` sia `input/Brenta_caldonazzo_borgo/`, che si ottiene con il completamento automatico del terminale |
| `<cartella di output>` | facoltativo: sottocartella di `output/`. Se manca, ha lo stesso nome della cartella di input |

Esempi:

```bash
julia --project=. src/main_sezioni.jl Brenta_caldonazzo_borgo     # → output/Brenta_caldonazzo_borgo
julia --project=. src/main_sezioni.jl Brenta_minimum              # → output/Brenta_minimum
julia --project=. src/main_sezioni.jl S_0050_6000_q2 S_50m_q2     # → output/S_50m_q2
julia --project=. src/main_sezioni.jl                             # caso predefinito: S_0050_6000_q2 → output/S_50m_q2
julia --project=. src/main_sezioni.jl --help                      # uso e cartelle disponibili in input/
```

All'avvio lo script stampa la cartella di input, il prefisso dei file delle sezioni che ha
riconosciuto e la cartella di output. Se la cartella di input non esiste, si ferma ed elenca
quelle disponibili.

I percorsi sono risolti rispetto alla posizione dello script, quindi la cartella di lavoro non
conta. Anche `--project=.` è facoltativo: se manca, lo script attiva da solo l'ambiente Julia
del progetto. Le dipendenze vanno comunque installate una volta con `Pkg.instantiate()`.

In alternativa, dal REPL (utile per avere i grafici a schermo). Il caso si sceglie riempiendo
`ARGS` prima di `include`:

```julia
julia> ]activate .
julia> append!(empty!(ARGS), ["Brenta_caldonazzo_borgo"])   # facoltativo: senza, caso predefinito
julia> include("src/main_sezioni.jl")                      # esegue già main(; dosave=true)
julia> tabelle, tabella_unica, geometria = main(; doplot=true)   # rilancia mostrando i grafici
```

Per cambiare caso dal REPL bisogna riavviare Julia: la cartella viene fissata al primo `include`.

`main` accetta:

| Argomento | Default | Effetto |
|---|---|---|
| `dosave` | `false` | scrive tabelle, geometria e grafici su disco |
| `doplot` | `false` | mostra i grafici di ogni sezione |

e restituisce `(tabelle, tabella_unica, geometria)`: un `Dict` nome → `DataFrame` per sezione,
la tabella unica e la tabella della geometria. Con più aste restituisce invece un `Dict`
asta → `(tabelle, tabella_unica, geometria)`.

> L'ultima riga di `main_sezioni.jl` è `main(; dosave=true)`: con `include` o lanciandolo da terminale,
> l'elaborazione parte subito e **salva**.

Con `β_dal_2D = true` viene letto anche il file HDF5 (~130 MB per `q2`): il calcolo è più lento
e richiede più memoria.

---

## Parametri di configurazione

Sono le costanti in testa a [src/main_sezioni.jl](src/main_sezioni.jl). Il caso da elaborare **non** si cambia
qui ma dalla riga di comando (vedi [Come lanciarlo](#come-lanciarlo)).

| Costante | Significato |
|---|---|
| `InputDefault`, `OutputDefault` | caso elaborato quando si lancia lo script senza argomenti |
| `β_dal_2D` | `true`: β dal campo di velocità 2D; `false`: β = 1 |
| `CartelleSezioni` | `true`: scrive anche una cartella per sezione con tabella completa e grafici `.svg`; `false` (predefinito): solo i file riassuntivi, utile con centinaia di sezioni |
| `AsteUnite` | con più aste, `true` (predefinito): un unico file per tutte le aste; `false`: una sottocartella per asta (vedi [Più aste](#più-aste)) |
| `SaltaSezioniSenzaFile` | `true` (predefinito): le sezioni di `centerpoints.csv` senza file dei risultati sono escluse con un avviso; `false`: l'elaborazione si ferma. Le sezioni escluse devono stare in testa o in coda a un'asta, non in mezzo: lì le pendenze verrebbero calcolate su una distanza più lunga |
| `OutputFolder` | `output/<cartella di output>` con β = 1, `output/<cartella di output>_beta` con β dal 2D, così i due casi non si sovrascrivono |
| `Format` | estensione dei file tabellari (`.csv`) |
| `g` | accelerazione di gravità [m/s²] |
| `J_lim_rel` | intervallo di plausibilità di `J`, in multipli della pendenza media del fondo dell'asta (0,25–2,5), usato **solo** per gli avvisi |
| `Ks_lim`, `β_lim` | intervalli di plausibilità di Ks e β, usati **solo** per gli avvisi |
| `Qrel_lim` | scarto relativo massimo tollerato tra la portata degli elementi e quella della sezione (calcolo di β) |
| `h_dry`, `Q_dry`, `Ks_dry` | tirante, portata e Ks usati per ricostruire la prima riga asciutta |

Di norma si usa **β = 1** (`β_dal_2D = false`): con β dal 2D i risultati del modello 1D cambiano
poco. Il calcolo di β resta disponibile per confronto e scrive in una cartella separata.

---

## Controlli e avvisi

Gli **errori** fermano l'elaborazione con un messaggio che indica il file e cosa correggere:
- cartella o file di input mancanti, colonne mancanti nei CSV;
- in `centerpoints.csv`: nomi ripetuti, nomi che non finiscono con un numero, valori mancanti,
  ascisse ripetute nella stessa asta, aste con una sola sezione;
- punti `_sx`/`_dx` mancanti in `bankpoints.csv`;
- sezioni con numero di righe o istanti diversi;
- sezioni sempre asciutte, o asciutte dopo il primo scalino bagnato;
- sezioni senza nessuna pendenza dell'energia positiva in un'asta con fondo non discendente;
- (con β dal 2D) mesh e risultati che non corrispondono, nodestring mancanti, CSV delle sezioni
  diversi dalle uscite dell'HDF5, area degli elementi diversa da quella della sezione.

Gli **avvisi** non modificano i dati: sono stampati in fondo, con il loro numero
(`Controlli: N valori fuori dai limiti`):
- portata non crescente con gli scalini (normale dove c'è esondazione a monte);
- `J`, `Ks` o `β` fuori dagli intervalli di plausibilità;
- talweg che sale da monte a valle (probabile `s_m` orientata al contrario);
- pendenze dell'energia non positive sostituite (interpolate);
- righe scartate perché il livello non cresce;
- (con β dal 2D) portata degli elementi lontana da quella della sezione.

Molti avvisi sono attesi: per esempio i Ks bassi ai primi scalini o le pendenze elevate vicino
all'uscita del dominio 2D. Vanno letti, non per forza eliminati.

---

## Diversione Brenta nuovo / Brenta vecchio (`main_diversione.jl`)

Il Brenta è a canale singolo fino alla sezione `B13200`, dove si divide nel **Brenta nuovo**
(`B13250`, canale principale) e nel **Brenta vecchio** (`V00050`, prima sezione dell'asta V con
risultati: `V00000` non è stata salvata). Nel modello 1D le due aste sono risolte in modo
indipendente: il Brenta nuovo ha una portata uscente verso il Brenta vecchio, che va imposta come
funzione della portata in arrivo da monte. [src/main_diversione.jl](src/main_diversione.jl)
ricava questa funzione dalla simulazione a scalini.

Comando pronto da copiare:

```bash
julia --project=. src/main_diversione.jl Brenta_caldonazzo_borgo
```

Forma generale:

```bash
julia --project=. src/main_diversione.jl [<cartella di input>]   # predefinita: Brenta_caldonazzo_borgo
```

L'argomento è una sottocartella di `input/`; l'output va nella sottocartella di `output/` con lo
stesso nome. Come gli altri script, attiva da solo l'ambiente del progetto.

### Come funziona

1. Legge la colonna `Discharge [m3/s]` dei file delle tre sezioni (`*<nome>.csv` nella cartella di
   input; il prefisso non conta, ma il file deve essere uno solo). Ogni riga è un salvataggio della
   simulazione, quindi le tre sezioni devono avere lo stesso numero di righe.
2. **Bilancio**: per ogni salvataggio calcola `Q_monte − Q_nuovo − Q_vecchio` e avvisa se supera
   `Bilancio_lim` (2 %) di `Q_monte`. Uno scarto piccolo e positivo alle portate alte è atteso:
   parte dell'acqua passa fuori dalle nodestring (in golena) fra `B13200` e la biforcazione.
3. **Monotonia**: la scala deve avere `Q_monte` crescente, quindi i salvataggi in cui `Q_monte` non
   supera il massimo dei precedenti sono scartati con un avviso.
4. Calcola la frazione `Q_vecchio / Q_monte` (0 quando `Q_monte` = 0), stampa la tabella a terminale
   seguita dagli avvisi e scrive `output/<cartella>/scala_diversione.csv`.

### `scala_diversione.csv`

| Colonna | Significato |
|---|---|
| `salvataggio` | riga dei file delle sezioni (1 = istante iniziale, asciutto) |
| `Q_monte_[m3/s]` | portata in arrivo da monte, a `B13200` |
| `Q_brenta_vecchio_[m3/s]` | portata verso il Brenta vecchio, a `V00050` |
| `Q_brenta_nuovo_[m3/s]` | portata nel Brenta nuovo, a `B13250` (solo controllo) |
| `frazione_brenta_vecchio_[-]` | `Q_brenta_vecchio / Q_monte` |
| `residuo_bilancio_[m3/s]` | `Q_monte − Q_nuovo − Q_vecchio` |

La scala da usare nel modello 1D è la coppia `Q_monte_[m3/s]` → `Q_brenta_vecchio_[m3/s]`,
interpolata linearmente. Con `Brenta_caldonazzo_borgo` ha 22 punti, da 0 a 211 m³/s: la quota verso
il Brenta vecchio scende da circa il 19 % a 2 m³/s all'8 % a 211 m³/s (16,6 m³/s); il bilancio chiude
entro 1,5 m³/s (0,7 %). Il grafico è in `grafici/scala_diversione.png` (vedi
[Grafici di controllo](#grafici-di-controllo-graficijl)).

### Parametri

Costanti in testa allo script:

| Costante | Valore | Significato |
|---|---|---|
| `InputDefault` | `Brenta_caldonazzo_borgo` | cartella usata senza argomenti |
| `SezMonte` | `B13200` | ultima sezione a canale singolo, a monte della biforcazione |
| `SezNuovo` | `B13250` | prima sezione del Brenta nuovo, usata solo per il bilancio |
| `SezVecchio` | `V00050` | prima sezione del Brenta vecchio con risultati |
| `Bilancio_lim` | 0,02 | scarto massimo del bilancio, relativo a `Q_monte`, oltre il quale si avvisa |

---

## Grafici di controllo (`grafici.jl`)

[src/grafici.jl](src/grafici.jl) disegna la situazione del caso e aiuta a trovare le sezioni in cui
le tabelle hanno problemi. Legge **solo** i file già salvati in `output/<cartella>`:
`localslope_xcoord.csv` e `tabella_sezioni.csv` (obbligatori), `topologia_sezioni.csv` (se manca, le
sezioni sono indicate con il numero e trattate come un'asta sola) e `scala_diversione.csv` (se manca,
il grafico della diversione è saltato). Salva i PNG in `output/<cartella>/grafici/`, sovrascrivendo
quelli già presenti.

Comando pronto da copiare:

```bash
julia --project=. src/grafici.jl Brenta_caldonazzo_borgo
```

Forma generale:

```bash
julia --project=. src/grafici.jl [<cartella di output>]   # predefinita: Brenta_caldonazzo_borgo
```

L'argomento è una sottocartella di `output/`. I grafici sono solo salvati, non mostrati a schermo
(GR senza finestra), quindi lo script gira anche senza interfaccia grafica.

### Grafici prodotti

| File | Contenuto | A cosa serve |
|---|---|---|
| `profilo_<asta>.png` | talweg (con il fondo campito), argine dx e sx e livello massimo in tabella lungo l'asta; linea verticale alla diversione. In arancione le sezioni in cui il livello massimo supera l'argine più basso | vedere dove le tabelle arrivano oltre l'alveo inciso, cioè dove il 2D esonda |
| `pendenza_<asta>.png` | pendenza locale del fondo (`local_slope`, in ‰; negativa se il fondo scende verso valle) con le contropendenze (> 0) in arancione | trovare contropendenze e salti del talweg |
| `controlli_<asta>.png` | cinque pannelli lungo l'asta: tirante massimo in tabella e altezza dell'argine più basso sul talweg; W al livello massimo; jden massimo (scala log); ΔQ massimo e minimo; punteggio di anomalia | individuare a colpo d'occhio le sezioni fuori dall'andamento delle vicine |
| `mappa_jden.png` | log₁₀(jden) su una mappa sezione (ascissa) × riga della tabella (ordinata, livello crescente); linea bianca fra un'asta e l'altra | vedere dove jden salta, per sezione e per livello |
| `mappa_DQ.png` | ΔQ sulla stessa mappa, scala divergente limitata al 99° percentile di \|ΔQ\| (blu: la portata aumenta rispetto alla sezione a monte, rosso: diminuisce) | localizzare esondazioni e rientri |
| `scala_diversione.png` | Q Brenta vecchio e Q Brenta nuovo in funzione di Q monte; quota verso il Brenta vecchio in %; residuo del bilancio | controllare la scala di `main_diversione.jl` |
| `sezioni/<nome>.png` | per una sezione, sei pannelli in funzione del livello: Ω, W, RH, jden (log), ΔQ e Ω_s, quest'ultima confrontata con la differenza finita ΔΩ/ΔWSE; righe orizzontali a quota talweg e argini | esaminare le sezioni sospette |
| `sezioni_selezionate.csv` | sezioni con grafico dedicato, motivo della scelta e indicatori di anomalia | sapere perché una sezione è stata scelta |

Le sezioni con grafico dedicato sono la prima e l'ultima di ogni asta, quelle della diversione
(`SezioniDiversione`) e le `NSezioniAnomale` con il punteggio di anomalia più alto. La cartella
`sezioni/` non viene svuotata: se cambia la scelta, i grafici delle sezioni non più selezionate
restano.

### Punteggio di anomalia

Serve solo a scegliere quali sezioni guardare, non è una verifica. Per ogni sezione è la somma di:

| Indicatore (colonna di `sezioni_selezionate.csv`) | Peso nel punteggio |
|---|---|
| `geometria`: righe in cui Ω o W calano al crescere del livello, o Ω_s o W_s sono negative | 1 per riga |
| `salto_RH`: righe in cui RH cala di oltre il 20 % rispetto alla riga precedente | 0,5 per riga |
| `jden_rel`: jden massimo della sezione diviso per la mediana dei jden massimi dell'asta | 2 · log₁₀, solo se > 1 |
| `err_Ωs`: massimo scarto relativo fra Ω_s e la differenza finita ΔΩ/ΔWSE | fino a 5 |
| `ΔQ_max`: \|ΔQ\| massimo [m³/s] | 1 ogni 10 m³/s |
| `righe`: righe in meno rispetto alla sezione con più righe | 0,5 per riga |

I pesi sono nella funzione `anomalie()`. Il punteggio è stampato a terminale per le sezioni più
anomale.

### Parametri

| Costante | Valore | Significato |
|---|---|---|
| `OutputDefault` | `Brenta_caldonazzo_borgo` | cartella usata senza argomenti |
| `SezioniDiversione` | `["B13200", "B13250", "V00050"]` | sezioni della diversione, sempre disegnate e segnate sul profilo; vettore vuoto se non c'è diversione |
| `NSezioniAnomale` | 12 | quante sezioni più anomale disegnare |

Le funzioni di `grafici.jl` (lettura degli output, `plot_sezione`, …) sono riusate da `W_check.jl`,
che include il file: per questo `main()` parte solo quando `grafici.jl` è lanciato direttamente.

---

## Controllo della larghezza W (`W_check.jl`)

Nel modello 1D la larghezza W e la sua derivata W_s sono interpolate in funzione del livello. Se in
una sezione W sale e poi torna indietro in pochi centimetri di livello, W_s diventa enorme e cambia
segno, e il modello 1D ne risente. Per esempio in `B13100` W passa da 10,82 a 12,82 m e torna a
10,82 m in 5 cm di livello, con W_s = −368. [src/W_check.jl](src/W_check.jl) cerca questi casi in
tutte le sezioni.

Comando pronto da copiare:

```bash
julia --project=. src/W_check.jl Brenta_caldonazzo_borgo
```

Forma generale:

```bash
julia --project=. src/W_check.jl [<cartella di output>]   # predefinita: Brenta_caldonazzo_borgo
```

Come `grafici.jl`, legge solo gli output già salvati in `output/<cartella>` (lo include e ne usa le
funzioni) e non modifica le tabelle.

### Come funziona

Per ogni sezione scorre le righe della tabella per livello crescente e considera ogni **calo di W**
fra una riga e la successiva più grande della tolleranza max(`TolAbs`, `TolRel` · W). Ogni calo
diventa un evento, classificato così:

| Tipo | Quando | Righe anomale |
|---|---|---|
| **picco** | prima del calo W era salita da un valore non superiore (entro la tolleranza) a quello dopo il calo | le righe del picco |
| **buco** | non è un picco, e più in alto W risale almeno al valore prima del calo | le righe del buco |
| **calo** | non è né un picco né un buco: W scende e non risale più | dalla riga dopo il calo a fine tabella |
| **W_s<0** | righe con W_s negativa non già comprese in un evento (calo di W entro la tolleranza) | quelle righe |

Gli aumenti di W, anche bruschi, non sono segnalati (vedi i limiti sotto).

### Output in `output/<cartella>/W_check/`

- `report_W_check.md`: riepilogo (tolleranza, numero di sezioni per tipo di evento) e, per ogni
  sezione con anomalie, una riga di descrizione per evento e il grafico della sezione. Esempio:

  > **picco** (righe 5-5): W sale da 10.82 a 12.82 m tra WSE 404.679 e 404.724 e torna a 10.82 m a
  > WSE 404.730 (spessore 0.051 m); W_s minima -368.4

- `W_check.csv`: un evento per riga, con le colonne:

  | Colonna | Significato |
  |---|---|
  | `section`, `nome`, `asta` | sezione |
  | `tipo` | `picco`, `buco`, `calo` o `W_s<0` |
  | `righe` | righe anomale della tabella della sezione, `prima-ultima` (1 = livello più basso) |
  | `WSE_da`, `WSE_a` | livelli della prima e dell'ultima riga anomala [m] |
  | `W_prima`, `W_anomala`, `W_dopo` | W prima dell'evento, valore estremo (massimo del picco, minimo del buco) e W dopo [m] |
  | `ampiezza` | \|W_anomala − W_dopo\| [m] |
  | `W_s_min` | W_s minima nelle righe dell'evento e in quelle adiacenti |
  | `descrizione` | la stessa frase del report |

- `sezioni/<nome>.png`: per ogni sezione con anomalie, lo stesso grafico di `grafici/sezioni/`,
  con le righe anomale in arancione in tutti i pannelli. La cartella è svuotata a ogni lancio, così
  contiene solo le sezioni dell'ultimo controllo.

A terminale stampa il numero di sezioni con anomalie e la tabella degli eventi.

Con `Brenta_caldonazzo_borgo` trova 5 picchi, in `B05900`, `B06950`, `B13100`, `B13250` e
`B16850`; il peggiore è `B16850`, con W_s = −2528 perché il picco si chiude in 1 mm di livello.

### Parametri e limiti

| Costante | Valore | Significato |
|---|---|---|
| `TolAbs` | 0,10 m | calo minimo di W considerato |
| `TolRel` | 0,02 | calo minimo di W relativo a W |

La cartella predefinita è quella di `grafici.jl` (`OutputDefault`). Il controllo non segnala:
- i **salti verso l'alto** di W fra due livelli molto vicini (es. `B13500`, dove W passa da 16 a 19 m
  in pochi mm): W_s diventa molto grande ma resta positiva;
- i problemi delle altre grandezze (Ω, RH, jden), che si vedono nei grafici di `grafici.jl`.

---

## Struttura della cartella

```
calcolo_tabelle_2Dto1D/
├── README.md
├── Project.toml          # dipendenze Julia
├── Manifest.toml         # versioni esatte (non versionato, vedi .gitignore)
├── src/
│   ├── main_sezioni.jl   # tabelle delle sezioni per il modello 1D
│   ├── main_diversione.jl # scala di deflusso della diversione Brenta nuovo / Brenta vecchio
│   ├── grafici.jl        # grafici di controllo dagli output (incluso da W_check.jl)
│   └── W_check.jl        # controllo dei picchi della larghezza W
├── input/                # un caso per sottocartella (vedi "Input")
├── output/               # un caso per sottocartella (vedi "Output")
└── BM_results/           # estrazione dei risultati BASEMENT (.h5 → CSV delle sezioni), progetto a parte
```

I file dei risultati `.h5`, sopra i 50 MB, sono gestiti con **Git LFS** (vedi
`.gitattributes`). Dopo un clone serve `git lfs pull` per scaricarli.

---

## Note e limitazioni

- **Sezioni = canale principale.** Le nodestring tagliano solo l'alveo inciso, non la piana. Con
  l'esondazione la portata cambia lungo l'asta, e lo scambio con la piana è descritto da `ΔQ`.
- **Moto permanente a scalini.** Ogni riga corrisponde a uno scalino a regime. Si assume che a
  ogni livello di una sezione corrisponda una sola portata, cioè niente isteresi.
- **Ks equivalente.** È un Ks "di sezione" che riassume l'attrito 2D in un valore 1D. Valori
  plausibili (10–50 m^(1/3)/s) indicano che l'approccio è fisicamente coerente.
- **Righe asciutte ricostruite.** Sono artificiali: vanno bene solo perché il modello 1D lavora
  a portate ben superiori (> 1 m³/s).
- **Unicode nei nomi delle colonne.** Le tabelle usano `Ω`, `β`, `ΔQ`. Chi le legge deve
  usare UTF-8.
- **Infiltrator** è una dipendenza di debug: `@infiltrate false` in fondo a `main` non ha
  effetto finché non lo si cambia in `true`.
