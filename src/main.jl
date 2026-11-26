# Usa sempre l'ambiente Julia del progetto (Project.toml), anche se lo script è lanciato senza
# --project=. ; la prima volta servono le dipendenze: julia --project=. -e 'using Pkg; Pkg.instantiate()'
import Pkg
let progetto = normpath(joinpath(@__DIR__, ".."))
    normpath(dirname(Base.active_project())) == progetto || Pkg.activate(progetto; io=devnull)
end
using Plots, Plots.Measures, Printf, CSV, DataFrames, Infiltrator, HDF5

#= CASO DA ELABORARE
    Si sceglie dalla riga di comando, senza modificare questo file:
        julia --project=. src/main.jl <cartella di input> [<cartella di output>]
    - <cartella di input>: nome di una sottocartella di input/ (es. S_0050_6000_q5; va bene anche
      "input/S_0050_6000_q5/");
    - <cartella di output>: nome della sottocartella di output/; se manca è uguale a quella di input.
    Senza argomenti si elabora il caso predefinito (InputDefault → OutputDefault).
    I percorsi sono relativi alla cartella del progetto (calcolo_tabelle_2Dto1D), quindi lo script si
    può lanciare da qualsiasi cartella.
    La cartella di input contiene centerpoints.csv, bankpoints.csv e i file delle sezioni
    (<prefisso><nome>.csv: il prefisso è riconosciuto da solo, vedi find_file_prefix) e, solo se β
    viene dal 2D, la mesh (*.2dm) e i risultati (*results.h5) del modello 2D.
=#
const InputDefault = "S_0050_6000_q2"
const OutputDefault = "S_50m_q2"
# β calcolato dal modello 2D (true) oppure β = 1 (false); con β = 1 mesh e risultati 2D non servono.
# Normalmente si usa β = 1; le tabelle con β dal 2D vanno in una cartella separata (suffisso "_beta")
# per non sovrascrivere quelle con β = 1.
const β_dal_2D = false
const ProjectFolder = normpath(joinpath(@__DIR__, ".."))

if any(in(("-h", "--help")), ARGS) || length(ARGS) > 2
    println("""
        Uso: julia --project=. src/main.jl [<cartella di input> [<cartella di output>]]
          <cartella di input>   sottocartella di input/ (predefinita: $InputDefault)
          <cartella di output>  sottocartella di output/ (predefinita: uguale a quella di input;
                                senza argomenti: $OutputDefault)
        Cartelle disponibili in input/: $(join(filter(d -> isdir(joinpath(ProjectFolder, "input", d)), readdir(joinpath(ProjectFolder, "input"))), ", "))""")
    exit(length(ARGS) > 2 ? 1 : 0)
end
const InputName = isempty(ARGS) ? InputDefault : basename(rstrip(ARGS[1], '/'))
const OutputName = length(ARGS) >= 2 ? basename(rstrip(ARGS[2], '/')) : (isempty(ARGS) ? OutputDefault : InputName)
const InputFolder = joinpath(ProjectFolder, "input", InputName)
const OutputFolder = joinpath(ProjectFolder, "output", β_dal_2D ? OutputName*"_beta" : OutputName)
const Format = ".csv"
const g = 9.81      # [m/s^2]

# Limiti dei valori ritenuti ragionevoli, usati solo per i controlli
# [-], limiti della pendenza dell'energia J rispetto alla pendenza media del fondo dell'asta
# (con 0.002 dà 0.0005-0.005)
const J_lim_rel = (0.25, 2.5)
const Ks_lim = (10.0, 50.0)         # [m^(1/3)/s]
const β_lim = (1.0, 1.5)            # [-]
const Qrel_lim = 0.10               # scarto relativo massimo tra la portata degli elementi e quella della sezione

# Valori usati per ricostruire le righe asciutte
const h_dry = 0.01      # [m], tirante
const Q_dry = 0.1       # [m3/s], portata
const Ks_dry = 20.0     # [m^(1/3)/s], Ks imposto

#= COLONNE RICHIESTE NEI FILE DI INPUT
    Le colonne sono cercate per nome (esatto, maiuscole comprese), mai per posizione: i file possono
    avere altre colonne, in qualsiasi ordine, che vengono ignorate. Per ogni file: colonne richieste
    con la loro descrizione ed esempio di file valido, usati nel messaggio d'errore.
=#
const COLONNE_CENTERPOINTS = [
    "nome" => "nome della sezione (es. S0400), uguale al nome della nodestring e al suffisso del file della sezione",
    "s_m"  => "ascissa curvilinea lungo l'asse [m], crescente verso valle",
]
const ESEMPIO_CENTERPOINTS = """
    fid,nome,s_m,dist_next
    1,S0050,50.25,50
    2,S0100,100,50
    Il nome è <asta><numero>: con più aste (es. N00050, V00050, B00000) ognuna è elaborata a parte
    e s_m è l'ascissa lungo la propria asta."""

const COLONNE_BANKPOINTS = [
    "nome"  => "<nome sezione>_sx o <nome sezione>_dx (argine sinistro/destro), due righe per sezione",
    "z_dtm" => "quota della sommità dell'argine [m s.l.m.]",
]
const ESEMPIO_BANKPOINTS = """
    fid,nome,s_m,z_dtm,dist_cp
    1,S0050_sx,50.25,433.46,9.2
    2,S0050_dx,50.25,433.52,6.7"""

const COLONNE_SEZIONE = [
    "Mean wse [m]"                => "livello medio della superficie libera",
    "Discharge [m3/s]"            => "portata",
    "Wetted area [m2]"            => "area bagnata (0 se la sezione è asciutta)",
    "Mean bottom elevation [m]"   => "quota media del fondo",
    "Reference elevation [m]"     => "quota di riferimento (talweg)",
    "Wetted geometric length [m]" => "larghezza (contorno) bagnata",
]
const ESEMPIO_SEZIONE = """
    È l'output per nodestring di BASEMENT 2D, una riga per uscita (la prima è la condizione iniziale):
    Mean wse [m],Discharge [m3/s],Wetted area [m2],Mean bottom elevation [m],Reference elevation [m],Wetted geometric length [m],...
    0.0,0.0,0.0,0.0,0.0,0.0,...
    430.73,1.0,2.40,430.48,430.46,9.32,..."""

# Legge un file CSV e controlla che abbia le colonne richieste; se ne manca qualcuna si ferma
# spiegando come deve essere organizzato il file
function read_input_table(path, colonne, esempio)
    data = CSV.read(path, DataFrame)
    mancanti = [c for (c, _) in colonne if !(c in names(data))]
    isempty(mancanti) && return data
    larghezza = maximum(length(c) for (c, _) in colonne)
    msg = "ERRORE nel file $path\n" *
          "  Colonne mancanti: " * join(("\"$c\"" for c in mancanti), ", ") * "\n" *
          "  Colonne trovate:  " * join(("\"$c\"" for c in names(data)), ", ") * "\n\n" *
          "  Il file deve essere un CSV separato da virgole, con una riga di intestazione e almeno\n" *
          "  queste colonne (il nome deve essere esatto, maiuscole comprese):\n" *
          join(("    " * rpad("\"$c\"", larghezza + 2) * "  " * d for (c, d) in colonne), "\n") * "\n" *
          "  Le altre colonne e l'ordine delle colonne non contano.\n\n" *
          "  Esempio:\n" * join(("    " * r for r in split(esempio, "\n")), "\n")
    printstyled(msg, "\n"; color=:red, bold=true)
    error("Colonne mancanti in $(basename(path)): elaborazione interrotta")
end

function stop_with(msg, breve)
    printstyled(msg, "\n"; color=:red, bold=true)
    error(breve)
end

#= PREFISSO DEI FILE DELLE SEZIONI
    I file delle sezioni si chiamano <prefisso><nome>.csv (es. ns_50_scalini_S0400.csv). Il prefisso
    è ricavato dai file della cartella che finiscono con il nome di una sezione di centerpoints.csv:
    deve essere lo stesso per tutte le sezioni trovate.
=#
function find_file_prefix(nomi)
    files = readdir(InputFolder)
    prefissi = Set{String}()
    for nome in nomi, f in files
        endswith(f, nome*Format) && push!(prefissi, f[begin:end-length(nome*Format)])
    end
    length(prefissi) == 1 && return only(prefissi)
    isempty(prefissi) && stop_with("ERRORE: in $InputFolder non c'è nessun file delle sezioni elencate in " *
        "centerpoints.csv.\n  I file devono chiamarsi <prefisso><nome>$Format, es. ns_50_scalini_$(first(nomi))$Format.",
        "File delle sezioni non trovati: elaborazione interrotta")
    stop_with("ERRORE: in $InputFolder i file delle sezioni hanno prefissi diversi: " *
        join(("\"$p\"" for p in prefissi), ", ") * "\n  Tutti i file delle sezioni devono avere lo stesso prefisso.",
        "Prefissi dei file delle sezioni diversi: elaborazione interrotta")
end

# Restituisce l'unico file della cartella di input che finisce con `fine` (es. ".2dm")
function find_unique_file(fine, descrizione)
    trovati = filter(f -> endswith(f, fine), readdir(InputFolder))
    length(trovati) == 1 && return only(trovati)
    stop_with("ERRORE: in $InputFolder serve un solo file $descrizione (*$fine), trovati $(length(trovati))" *
        (isempty(trovati) ? "" : ": " * join(trovati, ", ")),
        "File $descrizione non trovato o non univoco: elaborazione interrotta")
end

#= CONTROLLO DEI FILE DI INPUT
    Prima di elaborare controlla che ci siano tutti i file: quelli fissi della cartella di input (mesh e
    risultati 2D solo se β viene dal 2D) e quelli delle sezioni elencate in centerpoints.csv. Se ne manca
    qualcuno li elenca tutti e si ferma. Restituisce il prefisso dei file delle sezioni e i nomi dei file
    di mesh e risultati 2D (nothing se β = 1).
=#
function check_input_files()
    if !isdir(InputFolder)
        disponibili = filter(d -> isdir(joinpath(ProjectFolder, "input", d)), readdir(joinpath(ProjectFolder, "input")))
        stop_with("ERRORE: la cartella di input $InputFolder non esiste.\n  Cartelle disponibili in input/: " *
            join(disponibili, ", "), "Cartella di input inesistente: elaborazione interrotta")
    end
    mancanti = String[]
    for f in ("centerpoints.csv", "bankpoints.csv")
        isfile(joinpath(InputFolder, f)) || push!(mancanti, joinpath(InputFolder, f))
    end
    sezioni = joinpath(InputFolder, "centerpoints.csv")
    FilePrefix = ""
    if isfile(sezioni)
        # prima controlla il contenuto di centerpoints.csv (nomi, ascisse, aste), poi i file delle sezioni
        nomi = reduce(vcat, [String.(sez[!, "nome"]) for (_, sez) in read_reaches()])
        FilePrefix = find_file_prefix(nomi)
        for nome in nomi
            f = joinpath(InputFolder, FilePrefix*nome*Format)
            isfile(f) || push!(mancanti, f)
        end
    end
    isempty(mancanti) || stop_with("AVVISO: mancano $(length(mancanti)) file di input:\n  " * join(mancanti, "\n  "),
                                   "File di input mancanti: elaborazione interrotta")
    MeshFile = β_dal_2D ? find_unique_file(".2dm", "di mesh") : nothing
    ResultsFile = β_dal_2D ? find_unique_file("results.h5", "dei risultati 2D") : nothing
    return FilePrefix, MeshFile, ResultsFile
end

function slope!(slope, y, x)
    @. slope[begin:end-1] = (y[begin+1:end] - y[begin:end-1]) / (x[begin+1:end] - x[begin:end-1])
    slope[end] = slope[end-1]
end

# Legge i dati di una sezione e "ripara" i valori artificialmente messi a 0 quando è asciutta
function load_section(InputPath)
    data = read_input_table(InputPath, COLONNE_SEZIONE, ESEMPIO_SEZIONE)
    wet =data[!, "Wetted area [m2]"] .> 0.0
    # quota di riferimento dalla prima riga bagnata: in quelle asciutte BASEMENT scrive 0
    k1 = findfirst(wet)
    k1 === nothing && stop_with("ERRORE nel file $InputPath\n  La sezione è asciutta (area bagnata nulla) in tutte " *
        "le righe: toglierla da centerpoints.csv.", "Sezione sempre asciutta: elaborazione interrotta")
    zref = data[k1, "Reference elevation [m]"]
    data[.!wet, "Reference elevation [m]"] .= zref
    data[.!wet, "Mean bottom elevation [m]"] .= zref
    data[.!wet, "Mean wse [m]"] .= zref
    return data, wet
end

#= RICOSTRUZIONE DELLE RIGHE ASCIUTTE
    La simulazione 2D non può usare portate e livelli troppo bassi per questioni di stabilità,
    quindi i primi scalini possono essere asciutti. Il modello 1D serve per le piene e non vede
    mai portate sotto 1 m3/s, quindi basta che questi valori siano ragionevoli.
    Prima riga: tirante h_dry, portata Q_dry, larghezza e pendenza dell'energia J del primo livello
    bagnato. L'area viene da Manning imponendo Ks_dry, con RH = Ω/W:
        Q = Ks * Ω * (Ω/W)^(2/3) * sqrt(J)  =>  Ω = (Q * W^(2/3) / (Ks * sqrt(J)))^(3/5)
    L'area non è W*h_dry, altrimenti RH sarebbe h_dry e il Ks sarebbe irragionevole.
    Eventuali altre righe asciutte sono interpolate fino al primo livello bagnato.
    Larghezza, pendenze (J e iWSE) e β sono quelli del primo livello bagnato.
=#
function rebuild_dry!(data, wet, J, iWSE, β, nome)
    k1 = findfirst(wet)
    kd = findfirst(!, wet[k1:end])
    kd === nothing || stop_with("ERRORE: la sezione $nome è asciutta alla riga $(k1+kd-1), dopo essere stata bagnata " *
        "alla riga $k1.\n  Una sezione deve restare bagnata negli scalini di portata successivi: toglierla da " *
        "centerpoints.csv o controllare la simulazione 2D.", "Sezione asciutta dopo il primo livello bagnato: elaborazione interrotta")
    WSE, Q, Ω, W = data[!, "Mean wse [m]"], data[!, "Discharge [m3/s]"], data[!, "Wetted area [m2]"], data[!, "Wetted geometric length [m]"]
    zref = data[begin, "Reference elevation [m]"]
    WSE[begin] = zref + h_dry
    Q[begin]   = Q_dry
    Ω[begin]   = (Q_dry * cbrt(W[k1])^2 / (Ks_dry * sqrt(J[k1])))^(3/5)
    for k in 2:k1-1
        θ = (k - 1) / (k1 - 1)
        WSE[k] = WSE[begin] + θ * (WSE[k1] - WSE[begin])
        Q[k]   = Q[begin] + θ * (Q[k1] - Q[begin])
        Ω[k]   = Ω[begin] + θ * (Ω[k1] - Ω[begin])
    end
    W[begin:k1-1] .= W[k1]
    J[begin:k1-1] .= J[k1]
    iWSE[begin:k1-1] .= iWSE[k1]
    β[begin:k1-1] .= β[k1]
    return data, J, iWSE, β
end

#= PENDENZA LUNGO L'ASTA (pelo libero o energia)
    WSE e wet hanno una riga per ogni scalino di portata e una colonna per ogni sezione,
    con le sezioni ordinate da monte verso valle secondo l'ascissa s. A parità di riga i
    livelli appartengono allo stesso istante, quindi formano il profilo di superficie libera.
    La stessa funzione, applicata al carico totale H = WSE + u²/2g, dà la pendenza dell'energia.
    La pendenza è positiva quando il livello scende verso valle ed è calcolata:
    - alle differenze centrate, se le sezioni adiacenti sono entrambe bagnate;
    - con la differenza unilaterale, alle estremità o se una delle adiacenti è asciutta;
    - NaN, se la sezione è asciutta o se lo sono entrambe le adiacenti.
=#
function free_surface_slope(WSE, wet, s)
    nlevels, nsections = size(WSE)
    iWSE = fill(NaN, nlevels, nsections)
    for j in 1:nsections, k in 1:nlevels
        wet[k, j] || continue
        up = (j > 1 && wet[k, j-1]) ? j - 1 : j
        dn = (j < nsections && wet[k, j+1]) ? j + 1 : j
        up == dn && continue
        iWSE[k, j] = (WSE[k, up] - WSE[k, dn]) / (s[dn] - s[up])
    end
    return iWSE
end

#= PENDENZA DELL'ENERGIA E PENDENZE NON POSITIVE
    Ks e jden si ricavano dalla pendenza dell'energia J = -d/ds (WSE + u²/2g), con u = Q/Ω, e non da
    quella del pelo libero: Manning/Strickler lega l'attrito alla dissipazione di energia, e le due
    pendenze coincidono solo se la velocità non cambia lungo l'asta. Dove la sezione si allarga la
    corrente rallenta e il pelo libero può risalire anche se l'energia scende.
    Il termine cinetico usa α = 1 (velocità uniforme sulla sezione), coerente con β = 1.
    Con la portata locale di ogni sezione jden = Ks² R_H^(4/3)/g = u²/(g J): a ogni livello della
    tabella il modello 1D ha come pendenza d'attrito proprio J.
=#
function energy_slope(WSE, Q, Ω, wet, s)
    H = @. ifelse(wet, WSE + (Q / Ω)^2 / (2g), WSE)
    return free_surface_slope(H, wet, s)
end

#= J NON POSITIVA
    Resta possibile in pochi punti (es. velocità molto non uniforme sulla sezione, scambi laterali con
    la piana): lì Ks e jden non sono definiti. J viene interpolata linearmente in WSE fra i livelli
    validi (bagnati, con J > 0) più vicini, sotto e sopra, della stessa sezione; se ce n'è uno solo si
    prende quello; se la sezione non ne ha nessuno si usa la pendenza media del fondo dell'asta.
    jden è poi calcolato con la velocità del livello, quindi resta coerente con la sua portata.
    Restituisce l'elenco delle sostituzioni, da segnalare tra i controlli.
=#
function fix_slope!(J, WSE, wet, nomi, i_fondo)
    avvisi = String[]
    nlevels, nsections = size(J)
    originale = copy(J)
    for j in 1:nsections, k in 1:nlevels
        (wet[k, j] && !(originale[k, j] > 0.0)) || continue
        ok(m) = wet[m, j] && originale[m, j] > 0.0
        b = findfirst(ok, k+1:nlevels)
        b = b === nothing ? nothing : k + b
        a = findlast(ok, 1:k-1)
        if a !== nothing && b !== nothing
            # interpolazione in WSE se il livello è compreso fra i due, altrimenti in base alla riga
            θ = WSE[a, j] < WSE[k, j] < WSE[b, j] ? (WSE[k, j] - WSE[a, j]) / (WSE[b, j] - WSE[a, j]) : (k - a) / (b - a)
            nuova, da = originale[a, j] + θ * (originale[b, j] - originale[a, j]), "interpolata fra le righe $a e $b"
        elseif a !== nothing || b !== nothing
            m = something(a, b)
            nuova, da = originale[m, j], "presa dalla riga $m"
        elseif i_fondo > 0.0
            nuova, da = i_fondo, "pendenza media del fondo"
        else
            stop_with("ERRORE: la sezione $(nomi[j]) non ha nessuna pendenza dell'energia positiva e anche la " *
                "pendenza media del fondo dell'asta non è positiva.\n  Controllare che s_m cresca verso valle.",
                "Pendenze non definite: elaborazione interrotta")
        end
        push!(avvisi, @sprintf("%-6s riga %d  WSE=%8.3f  J=%.6f non positiva, sostituita con %.6f (%s)",
            nomi[j], k, WSE[k, j], originale[k, j], nuova, da))
        J[k, j] = nuova
    end
    return avvisi
end

#= LETTURA DELLA MESH (.2dm)
    Restituisce le coordinate dei nodi, i triangoli (nodi in ordine, numerati come nel file) e le
    nodestring (elenco dei nodi da un estremo all'altro). Nel file una nodestring è una riga
    "NS n1 n2 ... -nk nome": l'ultimo nodo ha il segno meno ed è seguito dal nome.
=#
function read_mesh(path)
    nodes = Dict{Int,NTuple{2,Float64}}()
    tri = NTuple{3,Int}[]
    strings = Dict{String,Vector{Int}}()
    cur = Int[]
    for line in eachline(path)
        t = split(line)
        isempty(t) && continue
        if t[1] == "ND"
            nodes[parse(Int, t[2])] = (parse(Float64, t[3]), parse(Float64, t[4]))
        elseif t[1] == "E3T"
            parse(Int, t[2]) == length(tri) + 1 || error("$path: triangoli non numerati in ordine")
            push!(tri, (parse(Int, t[3]), parse(Int, t[4]), parse(Int, t[5])))
        elseif t[1] == "NS"
            for v in t[2:end]
                n = tryparse(Int, v)
                if n === nothing
                    strings[v] = cur
                    cur = Int[]
                else
                    push!(cur, abs(n))
                end
            end
        end
    end
    return nodes, tri, strings
end

#= ELEMENTI A MONTE DI UNA NODESTRING
    Ogni tratto della nodestring è un lato della mesh, condiviso da due triangoli: si prende quello
    dalla parte di monte. Nel modello 2D ogni nodestring ha "upstream_direction" (right o left), cioè
    da che parte sta monte percorrendo la nodestring dal primo all'ultimo nodo. Per ogni tratto
    restituisce il triangolo, la lunghezza del lato e la normale unitaria verso valle.
=#
function upstream_elements(ns, nodes, tri, edges, σ)
    elementi = Tuple{Int,Float64,NTuple{2,Float64}}[]
    for (p, q) in zip(ns[begin:end-1], ns[begin+1:end])
        (xp, yp), (xq, yq) = nodes[p], nodes[q]
        dx, dy = xq - xp, yq - yp
        L = hypot(dx, dy)
        lato = get(edges, minmax(p, q), Int[])
        # σ = +1 se monte è a destra: il baricentro del triangolo di monte ha prodotto vettoriale negativo
        monte = filter(lato) do e
            xc = sum(nodes[k][1] for k in tri[e]) / 3
            yc = sum(nodes[k][2] for k in tri[e]) / 3
            σ * (dx * (yc - yp) - dy * (xc - xp)) < 0
        end
        length(monte) == 1 || error("Nodestring: nessun triangolo (o più di uno) a monte del lato $p-$q")
        push!(elementi, (monte[1], L, (-σ * dy / L, σ * dx / L)))
    end
    return elementi
end

#= COEFFICIENTE DI QUANTITÀ DI MOTO β DAL MODELLO 2D
    Come nel metodo di Engelund (eq. 13 delle slide):
        β = Ω/Q² ∫ u² Y dy
    ma con le velocità del modello 2D al posto di quelle di moto uniforme striscia per striscia.
    Le strisce sono i triangoli appena a monte dei lati della nodestring: tirante Y = quota del pelo
    libero - quota del fondo dell'elemento, portata specifica normale alla sezione qn = q·n, larghezza
    ΔL = lunghezza del lato. Quindi u² Y = qn²/Y e
        Ω = Σ Y ΔL,   Q = Σ qn ΔL,   β = Ω Σ(qn²/Y ΔL) / Q²
    Ω coincide con l'area bagnata dei file delle sezioni: BASEMENT la calcola con gli stessi elementi.
    Q invece no, perché BASEMENT la calcola dai flussi sui lati e non dai valori nelle celle
    (in mediana gli elementi danno il 2% in più). Per β si usa la Q degli elementi, così con
    velocità uniforme vale esattamente 1 e non è mai minore di 1.
    Gli elementi con tirante non superiore al minimo del modello 2D sono asciutti (portata nulla).
    Restituisce β, una riga per scalino e una colonna per sezione (NaN se la sezione è asciutta),
    e gli avvisi sulla portata degli elementi.
=#
function beta_2D(nomi, dati, wet_all, MeshFile, ResultsFile)
    nodes, tri, strings = read_mesh(joinpath(InputFolder, MeshFile))
    edges = Dict{Tuple{Int,Int},Vector{Int}}()
    for (e, t) in enumerate(tri), (p, q) in ((t[1], t[2]), (t[2], t[3]), (t[3], t[1]))
        push!(get!(edges, minmax(p, q), Int[]), e)
    end
    avvisi = String[]
    h5open(joinpath(InputFolder, ResultsFile)) do f
        # i triangoli devono essere le celle del file dei risultati (che numera i nodi da 0)
        topo = read(f["CellsAll/Topology"])
        size(topo, 2) == length(tri) && all(topo[:, e] == collect(tri[e]) .- 1 for e in eachindex(tri)) ||
            error("La mesh $MeshFile non corrisponde a quella di $ResultsFile")
        zb = vec(read(f["CellsAll/BottomEl"]))
        hmin = only(read(f["Parameters/MinWaterDepth"]))
        model = read(f[".config/model"])
        ns_h5 = read(f["NodeStrg/StrNames"])
        uscite = sort(keys(f["RESULTS/CellsAll/HydState"]))
        nlevels = nrow(dati[first(nomi)][1])
        length(uscite) == nlevels || error("$ResultsFile ha $(length(uscite)) uscite, i file delle sezioni $nlevels righe")
        H  = [read(f["RESULTS/CellsAll/HydState/"*u]) for u in uscite]   # quota, qx, qy per cella
        SV = [read(f["RESULTS/NodeStrg/StateVar/"*u]) for u in uscite]   # come i file delle sezioni
        β = fill(NaN, nlevels, length(nomi))
        for (j, nome) in enumerate(nomi)
            haskey(strings, nome) || error("La nodestring $nome non è in $MeshFile")
            js = findfirst(==(nome), ns_h5)
            js === nothing && error("La nodestring $nome non è in $ResultsFile")
            m = match(Regex("\"name\":\"$nome\",\"upstream_direction\":\"(right|left)\""), model)
            m === nothing && error("Direzione di monte della nodestring $nome non trovata in $ResultsFile")
            elementi = upstream_elements(strings[nome], nodes, tri, edges, m[1] == "right" ? 1 : -1)
            data = dati[nome][1]
            for k in 1:nlevels
                # le righe del file della sezione devono essere le uscite del file dei risultati (portata,
                # area e, se bagnata, livello: in quelle asciutte load_section lo ha cambiato)
                stessa = isapprox(SV[k][2, js], data[k, "Discharge [m3/s]"]; rtol=1e-6, atol=1e-6) &&
                         isapprox(SV[k][3, js], data[k, "Wetted area [m2]"]; rtol=1e-6, atol=1e-6) &&
                         (!wet_all[k, j] || isapprox(SV[k][1, js], data[k, "Mean wse [m]"]; rtol=1e-9))
                stessa || error("$nome riga $k: il file della sezione non corrisponde all'uscita $(uscite[k]) di $ResultsFile")
                wet_all[k, j] || continue
                Ω = Qe = M = 0.0
                for (e, L, (nx, ny)) in elementi
                    Y = H[k][1, e] - zb[e]
                    Y > hmin || continue
                    qn = H[k][2, e] * nx + H[k][3, e] * ny
                    Ω += Y * L
                    Qe += qn * L
                    M += qn^2 / Y * L
                end
                # tolleranza per gli arrotondamenti di BASEMENT: a 50 m gli scarti arrivano a 1,4e-4, un
                # elemento sbagliato cambierebbe l'area di qualche percento
                isapprox(Ω, data[k, "Wetted area [m2]"]; rtol=1e-3, atol=1e-3) ||
                    error(@sprintf("%s riga %d: area degli elementi %.4f diversa da quella della sezione %.4f",
                                   nome, k, Ω, data[k, "Wetted area [m2]"]))
                Qe > 0.0 || error("$nome riga $k: portata degli elementi non positiva")
                Qs = data[k, "Discharge [m3/s]"]
                abs(Qe / Qs - 1) <= Qrel_lim || push!(avvisi, @sprintf("%-6s Q=%7.3f  portata degli elementi %7.3f (%+.1f%%)",
                                                                       nome, Qs, Qe, 100 * (Qe / Qs - 1)))
                β[k, j] = Ω * M / Qe^2
            end
        end
        return β, avvisi
    end
end

#= CALCOLO DI KS
    Non è banale. La pendenza dell'energia va calcolata per ogni sezione e per
    ogni valore di WSE. Per semplicità assumo che sia unica per tutte le portate, cioè
    che a ogni WSE corrisponda la sola portata del relativo scalino.
    Il Ks comunque deve cambiare a seconda delle diverse regioni bagnate.
    Valori normali di Ks indicano un approccio fisicamente consistente.
    Se la pendenza non è definita o non è positiva il Ks non è definito (NaN).
    La portata è quella locale della sezione, perché quella che resta nel canale
    principale dopo l'esondazione.
=#
#= DIFFERENZA DI PORTATA TRA SEZIONI ADIACENTI
    Q e wet hanno una riga per ogni scalino e una colonna per ogni sezione, da monte verso valle.
    ΔQ = Q della sezione − Q della sezione precedente, allo stesso scalino (moto permanente): è lo
    scambio con la piana nel tratto tra le due sezioni, negativo se l'acqua esce e positivo se rientra.
    Vale 0 nella prima sezione, che non ha una sezione a monte, e se una delle due sezioni è asciutta
    (le righe asciutte ricostruite hanno quindi ΔQ = 0).
=#
function discharge_difference(Q, wet)
    ΔQ = zeros(size(Q))
    for j in 2:size(Q, 2), k in axes(Q, 1)
        (wet[k, j] && wet[k, j-1]) && (ΔQ[k, j] = Q[k, j] - Q[k, j-1])
    end
    return ΔQ
end

function Ks_equivalent(Q, Ω, RH, J)
    (J > 0.0 && Ω > 0.0) || return NaN
    return sqrt( Q^2 / ((Ω * cbrt(RH)^2 )^2 * J) )
end

# Calcola le caratteristiche geometriche/idrauliche di una sezione, una riga per ogni WSE.
# J è la pendenza dell'energia (usata per Ks e jden), iWSE quella del pelo libero (solo informativa).
function section_table(data, J, iWSE, β, ΔQ)
    Q   = data[!, "Discharge [m3/s]"]
    WSE = data[!, "Mean wse [m]"]

    # compila l'area
    Ω   = data[!, "Wetted area [m2]"]
    Ω_s = similar(Ω)
    slope!(Ω_s, Ω, WSE)

    # compila il coefficiente di quantità di moto: dal modello 2D (beta_2D) o 1, secondo β_dal_2D
    β_s = similar(β)
    slope!(β_s, β, WSE)

    # compila la larghezza
    W   = data[!, "Wetted geometric length [m]"]
    W_s = similar(W)
    slope!(W_s, W, WSE)

    # compila il raggio idraulico
    RH   = Ω ./ W
    RH_s = similar(RH)
    slope!(RH_s, RH, WSE)

    # calcola il Ks equivalente dell'intera sezione
    Ks = Ks_equivalent.(Q, Ω, RH, J)

    # denominatore del termine di attrito del modello 1D, γ = |u|/jden
    jden   = @. Ks^2 * cbrt(RH)^4 / g
    jden_s = similar(jden)
    slope!(jden_s, jden, WSE)

    # differenza di portata rispetto alla sezione precedente (vedi discharge_difference)
    ΔQ_s = similar(ΔQ)
    slope!(ΔQ_s, ΔQ, WSE)

    return DataFrame(; WSE, Ω, Ω_s, β, β_s, W, W_s, RH, RH_s, iWSE, J, Ks, jden, jden_s, ΔQ, ΔQ_s)
end

#= Prepara la tabella di una sezione nel formato letto da READ_VALUES! nel modello 1D: le prime 12
    colonne sono fisse e il modello legge solo quelle. In seconda colonna c'è il livello, non il
    tirante. In fondo ci sono ΔQ e ΔQ_s (colonne 13 e 14), per ora non usate dal modello.
=#
function model_table(table, sezione)
    colonne = ["WSE", "Ω", "Ω_s", "β", "β_s", "jden", "jden_s", "W", "W_s", "RH", "RH_s", "ΔQ", "ΔQ_s"]
    return insertcols(table[!, colonne], 1, :sezione => sezione)
end

#= CONTROLLI SUI RISULTATI
    La simulazione prevede un'esondazione e le sezioni riguardano solo il canale principale,
    quindi la portata può variare lungo il tratto: se ne terrà conto nel modello 1D con un
    termine di portata laterale uscente. Per questo la portata di una sezione è confrontata
    solo con i propri scalini precedenti, e deve essere crescente.
    I valori fuori dai limiti sono solo segnalati, non modificati.
=#
function check_section(nome, table, Q, i_fondo)
    J_lim = J_lim_rel .* i_fondo
    avvisi = String[]
    WSE, J, Ks, β = table[!, "WSE"], table[!, "J"], table[!, "Ks"], table[!, "β"]
    for k in 2:nrow(table)
        Q[k] > Q[k-1] || push!(avvisi, @sprintf("%-6s WSE=%8.3f  Q non crescente: %.3f -> %.3f", nome, WSE[k], Q[k-1], Q[k]))
        # β vale almeno 1: sotto, di poco, solo per arrotondamento (un solo elemento bagnato)
        β_lim[1] - 1e-9 <= β[k] <= β_lim[2] || push!(avvisi, @sprintf("%-6s WSE=%8.3f  Q=%7.3f  β=%.3f", nome, WSE[k], Q[k], β[k]))
        J_lim[1] <= J[k] <= J_lim[2] || push!(avvisi, @sprintf("%-6s WSE=%8.3f  Q=%7.3f  J=%.5f", nome, WSE[k], Q[k], J[k]))
        Ks_lim[1] <= Ks[k] <= Ks_lim[2] || push!(avvisi, @sprintf("%-6s WSE=%8.3f  Q=%7.3f  Ks=%.2f", nome, WSE[k], Q[k], Ks[k]))
    end
    return avvisi
end

# Crea i grafici di una sezione, mostrandoli e/o salvandoli nella cartella della sezione
function plot_section(table, nome, Folder; dosave::Bool=false, doplot::Bool=false)
    grafici = [
        (:Ω,    "area_bagnata",        "Area bagnata",           "Area [m^2]"),
        (:W,    "larghezza_bagnata",   "Larghezza bagnata",      "Larghezza [m]"),
        (:RH,   "raggio.idraulico",    "Raggio idraulico",       "Raggio [m]"),
        (:iWSE, "pendenza_wse",        "Pendenza sup. libera",   "Pendenza [-]"),
        (:J,    "pendenza_energia",    "Pendenza dell'energia",  "Pendenza [-]"),
        (:Ks,   "ks",                  "Ks equivalente",         "Ks [m^(1/3)/s]"),
        (:β,    "beta",                "Coefficiente β",         "β [-]"),
        (:jden, "jden",                "jden",                   "jden [m]"),
    ]
    β_dal_2D || filter!(gr -> gr[1] != :β, grafici)    # con β = 1 il grafico di β non serve
    for (colonna, file, titolo, xlabel) in grafici
        plt = plot(table[!, colonna], table[!, "WSE"];
            title="$titolo - $nome",
            xlabel=xlabel,
            ylabel="Livello idrometrico [m]",
            legend=false,
            size=(450, 600)
        )
        doplot && display(plt)
        dosave && savefig(plt, joinpath(Folder, file*".svg"))
    end
    return nothing
end

#= SEZIONI E ASTE
    Legge centerpoints.csv e controlla nomi e ascisse. Le aste sono riconosciute dal nome delle sezioni,
    che deve essere <prefisso><numero>: il prefisso è il nome dell'asta (es. N00050 → asta N, B17650 →
    asta B, S0400 → asta S). Ogni asta ha la sua ascissa s_m ed è elaborata come un tratto indipendente,
    ognuna nella sua sottocartella di output se le aste sono più di una. Restituisce un vettore di
    (nome asta, sezioni ordinate da monte verso valle secondo s_m).
=#
reach_name(nome) = (m = match(r"^(.*?)\d+$", String(nome)); m === nothing ? nothing : String(m[1]))

function read_reaches()
    path = joinpath(InputFolder, "centerpoints.csv")
    cp = read_input_table(path, COLONNE_CENTERPOINTS, ESEMPIO_CENTERPOINTS)
    problemi = String[]
    any(ismissing, cp[!, "nome"]) && push!(problemi, "ci sono righe senza \"nome\"")
    any(ismissing, cp[!, "s_m"]) && push!(problemi, "ci sono righe senza \"s_m\"")
    eltype(skipmissing(cp[!, "s_m"])) <: Real || push!(problemi, "la colonna \"s_m\" contiene valori non numerici")
    if isempty(problemi)
        cp[!, "nome"] = String.(strip.(string.(cp[!, "nome"])))
        doppi = [n for (n, c) in countmap_(cp[!, "nome"]) if c > 1]
        isempty(doppi) || push!(problemi, "nomi ripetuti: " * join(sort(doppi), ", "))
        senza_numero = filter(n -> reach_name(n) === nothing, cp[!, "nome"])
        isempty(senza_numero) || push!(problemi, "nomi che non finiscono con un numero: " * join(senza_numero, ", ") *
            "\n  Il nome di una sezione deve essere <asta><numero>, es. N00050 (asta N) o S0400 (asta S).")
    end
    isempty(problemi) || stop_with("ERRORE nel file $path\n  " * join(problemi, "\n  "),
                                   "centerpoints.csv non valido: elaborazione interrotta")
    cp[!, "asta"] = reach_name.(cp[!, "nome"])
    reaches = Tuple{String,DataFrame}[]
    for asta in sort(unique(cp[!, "asta"]))
        sez = sort(cp[cp[!, "asta"] .== asta, :], "s_m")
        nrow(sez) >= 2 || stop_with("ERRORE nel file $path\n  L'asta \"$asta\" ha una sola sezione ($(only(sez.nome))): " *
            "ne servono almeno due per calcolare le pendenze.\n  Il nome dell'asta è la parte del nome della sezione " *
            "prima del numero finale.", "Asta con una sola sezione: elaborazione interrotta")
        s = sez[!, "s_m"]
        k = findfirst(i -> s[i+1] <= s[i], 1:length(s)-1)
        k === nothing || stop_with("ERRORE nel file $path\n  Asta \"$asta\": le sezioni $(sez.nome[k]) e $(sez.nome[k+1]) " *
            "hanno la stessa ascissa s_m = $(s[k]).", "Ascisse ripetute: elaborazione interrotta")
        push!(reaches, (asta, sez))
    end
    return reaches
end

countmap_(v) = (d = Dict{eltype(v),Int}(); foreach(x -> d[x] = get(d, x, 0) + 1, v); d)

#= ELABORAZIONE DI UN'ASTA
    Calcola tabelle e geometria delle sezioni di un'asta (ordinate da monte verso valle) e, se dosave,
    le scrive in Folder. Restituisce tabelle per sezione, tabella unica, geometria e avvisi.
=#
function process_reach(sections_all, Folder, FilePrefix, MeshFile, ResultsFile, z_argine; dosave, doplot)
    # leggi i dati delle sezioni, cioè quelle dell'asta elencate in centerpoints.csv
    nomi = String.(sections_all[!, "nome"])
    dati = Dict(nome => load_section(joinpath(InputFolder, FilePrefix*nome*Format)) for nome in nomi)
    # le righe di sezioni diverse devono essere gli stessi scalini: stesso numero di righe e, se c'è
    # la colonna del tempo, stessi istanti
    ref = dati[first(nomi)][1]
    for nome in nomi
        d = dati[nome][1]
        nrow(d) == nrow(ref) || stop_with("ERRORE: la sezione $nome ha $(nrow(d)) righe, la sezione $(first(nomi)) " *
            "ne ha $(nrow(ref)).\n  Tutti i file delle sezioni devono avere una riga per ogni uscita della simulazione 2D.",
            "Sezioni con numero di righe diverso: elaborazione interrotta")
        "t [s]" in names(d) && "t [s]" in names(ref) && (d[!, "t [s]"] == ref[!, "t [s]"] ||
            stop_with("ERRORE: la sezione $nome ha istanti diversi dalla sezione $(first(nomi)).",
                      "Sezioni con istanti diversi: elaborazione interrotta"))
    end

    # geometria delle sezioni: talweg e pendenza locale del fondo
    s = sections_all[!, "s_m"]
    talweg = [dati[nome][1][begin, "Reference elevation [m]"] for nome in nomi]
    avvisi = String[]
    talweg[end] < talweg[begin] || push!(avvisi, @sprintf("%-6s il talweg sale da monte (%.2f m) a valle (%.2f m): controllare che s_m cresca verso valle",
        nomi[end], talweg[begin], talweg[end]))
    # pendenza media del fondo sull'intera asta (positiva se il fondo scende verso valle)
    i_media = (talweg[begin] - talweg[end]) / (s[end] - s[begin])

    # pendenze per ogni sezione e per ogni WSE: del pelo libero (solo informativa) e dell'energia,
    # usata per Ks e jden; le J non positive vengono sostituite (vedi fix_slope!)
    WSE_all = reduce(hcat, [dati[nome][1][!, "Mean wse [m]"] for nome in nomi])
    Q_all   = reduce(hcat, [dati[nome][1][!, "Discharge [m3/s]"] for nome in nomi])
    Ω_all   = reduce(hcat, [dati[nome][1][!, "Wetted area [m2]"] for nome in nomi])
    wet_all = reduce(hcat, [dati[nome][2] for nome in nomi])
    iWSE_all = free_surface_slope(WSE_all, wet_all, s)
    J_all = energy_slope(WSE_all, Q_all, Ω_all, wet_all, s)
    append!(avvisi, fix_slope!(J_all, WSE_all, wet_all, nomi, abs(i_media)))

    # differenza di portata rispetto alla sezione precedente, dai dati 2D originali (prima di
    # ricostruire le righe asciutte)
    ΔQ_all = discharge_difference(Q_all, wet_all)

    # calcola β per ogni sezione e per ogni WSE dai risultati del modello 2D, oppure β = 1
    if β_dal_2D
        β_all, avvisi_β = beta_2D(nomi, dati, wet_all, MeshFile, ResultsFile)
        append!(avvisi, avvisi_β)
    else
        β_all = ones(size(WSE_all))
    end

    # pendenza locale del fondo: calcolata come quella della superficie libera, ma con il segno usato
    # dal modello 1D, negativo quando il fondo scende verso valle
    local_slope = -vec(free_surface_slope(permutedims(talweg), trues(1, length(nomi)), s))
    # l'ultima sezione può essere in contropendenza, ma il modello 1D ne usa la pendenza per la portata
    # di uscita: uso la pendenza media del fondo sull'intera asta
    local_slope[end] = -i_media

    # ciclo sulle sezioni, da monte verso valle
    tabelle = Dict{String,DataFrame}()
    @printf("%-10s %6s %17s %17s %17s\n", "Sezione", "n WSE", "J [-]", "Ks [m^(1/3)/s]", "β [-]")
    for (j, nome) in enumerate(nomi)
        data, wet = dati[nome]
        iWSE = iWSE_all[:, j]
        J = J_all[:, j]
        β = β_all[:, j]
        ΔQ = ΔQ_all[:, j]

        # sostituisci le righe asciutte con valori ricostruiti
        rebuild_dry!(data, wet, J, iWSE, β, nome)

        # il modello 1D vuole livelli crescenti: scarta le righe con livello non superiore al massimo
        # dei precedenti (es. S0050 ogni 50 m, dove negli ultimi scalini il livello scende)
        WSE = data[!, "Mean wse [m]"]
        tieni = [k == 1 || WSE[k] > maximum(WSE[1:k-1]) for k in eachindex(WSE)]
        for k in findall(.!tieni)
            push!(avvisi, @sprintf("%-6s WSE=%8.3f  Q=%7.3f  livello non crescente: riga scartata", nome, WSE[k], data[k, "Discharge [m3/s]"]))
        end
        data = data[tieni, :]
        iWSE = iWSE[tieni]
        J = J[tieni]
        β = β[tieni]
        ΔQ = ΔQ[tieni]

        #= ULTIMA SEZIONE: CONDIZIONE DI USCITA
            Il modello 1D calcola la portata in uscita come moto uniforme, Q = Ω·√(g·jden·S), con il jden
            della tabella dell'ultima sezione e la sua pendenza locale S (la media del fondo, vedi sopra).
            Qui quindi Ks e jden sono ricavati con quella pendenza invece che con la pendenza dell'energia: così, a ogni livello della tabella, il modello dà in uscita la portata degli scalini 2D.
            La pendenza del pelo libero nel 2D tra le ultime due sezioni è molto più ripida (fino a 0,0075
            contro 0,00213), probabilmente per il richiamo dell'uscita del dominio 2D: con quella e il Ks 33
            di settings.txt il tirante in uscita era 10-13 cm più alto del 2D.
        =#
        if j == length(nomi)
            J = fill(abs(local_slope[end]), length(J))
        end

        table = section_table(data, J, iWSE, β, ΔQ)
        tabelle[nome] = table

        i_ok = table[!, "J"]
        K_ok = table[!, "Ks"]
        @printf("%-10s %6d %8.5f %8.5f %8.2f %8.2f %8.3f %8.3f\n", nome, nrow(table),
            minimum(i_ok), maximum(i_ok), minimum(K_ok), maximum(K_ok), minimum(β), maximum(β))
        append!(avvisi, check_section(nome, table, data[!, "Discharge [m3/s]"], abs(i_media)))

        # salva tabella e grafici solo se la tabella è diversa da quella già salvata
        SecFolder = joinpath(Folder, nome)
        TablePath = joinpath(SecFolder, "tabella"*Format)
        salva = dosave && !(isfile(TablePath) && isequal(CSV.read(TablePath, DataFrame), table))
        if salva
            mkpath(SecFolder)
            CSV.write(TablePath, table)
        end
        (salva || doplot) && plot_section(table, nome, SecFolder; dosave=salva, doplot)
    end

    # tabella unica, con il numero progressivo della sezione da monte verso valle
    tabella_unica = reduce(vcat, [model_table(tabelle[nome], j) for (j, nome) in enumerate(nomi)])

    # quote degli argini destro e sinistro di ogni sezione, dai punti in bankpoints.csv
    argine_dx = [z_argine[nome*"_dx"] for nome in nomi]
    argine_sx = [z_argine[nome*"_sx"] for nome in nomi]
    geometria = DataFrame("section" => eachindex(nomi), "local_slope" => local_slope,
        "x_coord_[m]" => s, "talweg_[m]" => talweg,
        "argine_dx_[m]" => argine_dx, "argine_sx_[m]" => argine_sx)
    if dosave
        mkpath(Folder)
        CSV.write(joinpath(Folder, "tabella_sezioni"*Format), tabella_unica)
        CSV.write(joinpath(Folder, "localslope_xcoord"*Format), geometria)
        println("Tabelle e geometria scritte in $Folder")
    end
    return tabelle, tabella_unica, geometria, avvisi
end

# Quote degli argini da bankpoints.csv: controlla che ci siano entrambi i punti di ogni sezione
function read_banks(nomi)
    BankPath = joinpath(InputFolder, "bankpoints"*Format)
    bankpoints = read_input_table(BankPath, COLONNE_BANKPOINTS, ESEMPIO_BANKPOINTS)
    z_argine = Dict(strip(string(r["nome"])) => r["z_dtm"] for r in eachrow(bankpoints)
                    if !ismissing(r["nome"]) && !ismissing(r["z_dtm"]))
    senza_argine = [nome*lato for nome in nomi for lato in ("_sx", "_dx") if !haskey(z_argine, nome*lato)]
    isempty(senza_argine) || stop_with("ERRORE nel file $BankPath\n  Mancano (o non hanno z_dtm) i punti: " *
        join(senza_argine, ", ") * "\n  Ogni sezione di centerpoints.csv deve avere due righe, <nome>_sx e <nome>_dx.",
        "Punti mancanti in bankpoints.csv: elaborazione interrotta")
    return z_argine
end

#= Elabora tutte le aste. Con un'unica asta (tutte le sezioni con lo stesso prefisso) scrive in
    OutputFolder e restituisce (tabelle, tabella_unica, geometria); con più aste scrive ogni asta in
    OutputFolder/<asta> e restituisce un Dict asta => (tabelle, tabella_unica, geometria).
=#
function main(; dosave::Bool=false, doplot::Bool=false)
    FilePrefix, MeshFile, ResultsFile = check_input_files()
    reaches = read_reaches()
    multi = length(reaches) > 1
    println("Input:  $InputFolder (file delle sezioni: $(FilePrefix)<nome>$Format)")
    println("Output: $OutputFolder")
    multi && println("Aste riconosciute dai nomi delle sezioni: " *
                     join(("$a ($(nrow(sez)) sezioni)" for (a, sez) in reaches), ", "))
    z_argine = read_banks(reduce(vcat, [String.(sez[!, "nome"]) for (_, sez) in reaches]))

    risultati = Dict{String,Any}()
    avvisi = String[]
    for (asta, sez) in reaches
        multi && printstyled("\n== Asta $asta ==\n"; bold=true)
        Folder = multi ? joinpath(OutputFolder, asta) : OutputFolder
        tabelle, tabella_unica, geometria, avvisi_asta =
            process_reach(sez, Folder, FilePrefix, MeshFile, ResultsFile, z_argine; dosave, doplot)
        risultati[asta] = (tabelle, tabella_unica, geometria)
        append!(avvisi, multi ? ("[$asta] " .* avvisi_asta) : avvisi_asta)
    end

    println("\nControlli: $(length(avvisi)) valori fuori dai limiti")
    foreach(println, avvisi)

    @infiltrate false

    return multi ? risultati : only(values(risultati))
end

main(; dosave=true)
