# calcolo_tabelle_2Dto1D

Script in Julia che trasforma i risultati di una simulazione **BASEMENT 2D** in **tabelle
geometrico-idrauliche per sezione**, da usare come input di un modello 1D (moto vario, sezioni
tabellate).

L'idea è semplice. Nel modello 2D si fa crescere la portata a **scalini**, raggiungendo il
moto permanente a ogni scalino, e lungo l'asta si registrano le grandezze integrali su una
serie di **nodestring** (le sezioni). Ogni scalino dà, per ogni sezione, un punto
"livello → area, larghezza, …". Mettendo in fila gli scalini si ottiene, per ogni sezione, una
tabella in funzione del livello idrometrico (WSE) che il modello 1D può interpolare.

---

## Indice

1. [Cosa fa lo script](#cosa-fa-lo-script)
2. [Grandezze tabellate](#grandezze-tabellate)
3. [Input](#input)
4. [Output](#output)
5. [Come lanciarlo](#come-lanciarlo)
6. [Parametri di configurazione](#parametri-di-configurazione)
7. [Controlli e avvisi](#controlli-e-avvisi)
8. [Diversione Brenta nuovo / Brenta vecchio](#diversione-brenta-nuovo--brenta-vecchio)
9. [Struttura della cartella](#struttura-della-cartella)
10. [Note e limitazioni](#note-e-limitazioni)

---

## Cosa fa lo script

Il calcolo delle tabelle è in [src/main_sezioni.jl](src/main_sezioni.jl) (la scala di deflusso della diversione è a parte, vedi [Diversione](#diversione-brenta-nuovo--brenta-vecchio)). In ordine:

1. **Controllo dei file di input.** Verifica che esistano `centerpoints.csv`, `bankpoints.csv`
   e un file per ogni sezione elencata in `centerpoints.csv` (più mesh e risultati 2D se β è
   calcolato dal 2D). Se manca qualcosa elenca tutti i file mancanti e si ferma.
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
12. **Scrittura** delle tabelle per sezione, dei grafici, della tabella unica per il modello 1D e
    del file con la geometria longitudinale.

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
| `S_0050_6000_q2` | 120, ogni 50 m (S0050–S6000) | 52 + riga iniziale, scalini di ~2 m³/s | caso attuale; mesh e `.h5` presenti |
| `S_0050_6000_q5` | 120, ogni 50 m | 22 + riga iniziale, scalini di ~5 m³/s | mesh e `.h5` presenti |
| `S_0200_6000` | 30, ogni 200 m (prefisso `ns_channel_`) | 22 + riga iniziale | non c'è `centerpoints.csv` (c'è `sections_dx.csv`); niente mesh/`.h5` |
| `Brenta_minimum` | 735 su tre aste: N (321), V (60), B (354), ogni ~50 m (prefisso `minimum_scalini_`) | 21 + riga iniziale | niente mesh/`.h5` |

### Più aste

Se i nomi delle sezioni hanno prefissi diversi (es. `N00050`, `V00050`, `B00000`), ogni prefisso è
un'asta e viene elaborata **indipendentemente** dalle altre:
- ha la sua ascissa `s_m`, misurata lungo la propria asta, che può ripartire da 0;
- le pendenze, ΔQ e la condizione di uscita si calcolano solo con le sezioni della stessa asta;
- la prima sezione di ogni asta ha ΔQ = 0 e l'ultima è trattata come sezione di uscita;
- i risultati vanno in una sottocartella per asta (vedi [Output](#output)).

Lo script non sa come le aste sono collegate fra loro (confluenze, diramazioni): i collegamenti
vanno gestiti nel modello 1D. Un'asta deve avere almeno due sezioni.

---

## Output

Gli output vanno in `output/<caso>/` (cartella `OutputFolder`; con β dal 2D il nome ha il
suffisso `_beta`):

Con più aste, la stessa struttura si ripete in una sottocartella per asta:
`output/<caso>/<asta>/tabella_sezioni.csv`, `output/<caso>/<asta>/<nome sezione>/…`. Il numero
progressivo delle sezioni riparte da 1 in ogni asta.

```
output/<caso>/
├── tabella_sezioni.csv        # tabella unica per il modello 1D
├── localslope_xcoord.csv      # geometria longitudinale
└── <nome sezione>/            # una cartella per sezione
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

### `<nome sezione>/tabella.csv`

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
| `S_50m_q2` | caso `S_0050_6000_q2`, β = 1 (configurazione attuale dello script) |
| `S_50m` | sezioni ogni 50 m, β = 1, 22 scalini; versione precedente dello script (senza `ΔQ`) |
| `S_50m_beta` | come `S_50m`, con β dal 2D |
| `S_200m` | caso `S_0200_6000`; versione precedente dello script |

---

## Come lanciarlo

Serve **Julia ≥ 1.12** (il `Manifest.toml` è generato con Julia 1.12.6). Le dipendenze sono in
[Project.toml](Project.toml): CSV, DataFrames, HDF5, Plots, Printf, Infiltrator.

La prima volta, dalla cartella del progetto:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Poi si sceglie il caso dalla riga di comando, senza modificare lo script:

```bash
julia --project=. src/main_sezioni.jl <cartella di input> [<cartella di output>]
```

| Argomento | Significato |
|---|---|
| `<cartella di input>` | sottocartella di `input/`. Va bene sia `S_0050_6000_q5` sia `input/S_0050_6000_q5/`, che si ottiene con il completamento automatico del terminale |
| `<cartella di output>` | facoltativo: sottocartella di `output/`. Se manca, ha lo stesso nome della cartella di input |

Esempi:

```bash
julia --project=. src/main_sezioni.jl S_0050_6000_q5              # → output/S_0050_6000_q5
julia --project=. src/main_sezioni.jl S_0050_6000_q5 S_50m_q5     # → output/S_50m_q5
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
julia> append!(empty!(ARGS), ["S_0050_6000_q5"])   # facoltativo: senza, caso predefinito
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

## Diversione Brenta nuovo / Brenta vecchio

Il Brenta è a canale singolo fino alla sezione `B13200`, dove si divide nel **Brenta nuovo**
(`B13250`, canale principale) e nel **Brenta vecchio** (`V00050`, prima sezione dell'asta V con
risultati). Nel modello 1D le due aste sono risolte in modo indipendente: il Brenta nuovo ha una
portata uscente verso il Brenta vecchio, ricavata dalla simulazione a scalini con
[src/main_diversione.jl](src/main_diversione.jl):

```bash
julia --project=. src/main_diversione.jl [<cartella di input>]   # predefinita: Brenta_caldonazzo_borgo
```

Per ogni salvataggio legge la portata delle tre sezioni e scrive in
`output/<cartella>/scala_diversione.csv` la scala di deflusso portata-portata:

| Colonna | Significato |
|---|---|
| `salvataggio` | riga dei file delle sezioni (1 = istante iniziale, asciutto) |
| `Q_monte_[m3/s]` | portata in arrivo da monte, a `B13200` |
| `Q_brenta_vecchio_[m3/s]` | portata verso il Brenta vecchio, a `V00050` |
| `Q_brenta_nuovo_[m3/s]` | portata nel Brenta nuovo, a `B13250` (solo controllo) |
| `frazione_brenta_vecchio_[-]` | `Q_brenta_vecchio / Q_monte` |
| `residuo_bilancio_[m3/s]` | `Q_monte − Q_nuovo − Q_vecchio` |

I salvataggi con `Q_monte` non crescente sono scartati, e viene segnalato quando il bilancio non
chiude entro il 2% di `Q_monte`. Le sezioni sono costanti in testa allo script.

---

## Struttura della cartella

```
calcolo_tabelle_2Dto1D/
├── README.md
├── Project.toml          # dipendenze Julia
├── Manifest.toml         # versioni esatte (non versionato, vedi .gitignore)
├── src/
│   ├── main_sezioni.jl   # tabelle delle sezioni per il modello 1D
│   └── main_diversione.jl # scala di deflusso della diversione Brenta nuovo / Brenta vecchio
├── input/                # un caso per sottocartella (vedi "Input")
└── output/               # un caso per sottocartella (vedi "Output")
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
