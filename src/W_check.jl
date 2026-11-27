#= CONTROLLO DELLA LARGHEZZA W DELLE SEZIONI
    Cerca nelle tabelle già salvate (output/<cartella>/tabella_sezioni.csv) le righe in cui la larghezza
    bagnata W non cresce regolarmente con il livello WSE. Un picco di W (W sale e poi torna indietro in
    pochi cm di livello, come in B13100) dà derivate W_s enormi e di segno opposto, che creano problemi
    nel modello 1D.
    Per ogni sezione, scorrendo le righe per livello crescente, ogni calo di W oltre la tolleranza
    (max(TolAbs, TolRel·W)) è classificato come:
      - picco: W era salita da un valore pari a quello dopo il calo (righe anomale: quelle del picco)
      - buco:  W cala e poi risale almeno al valore precedente (righe anomale: quelle del buco)
      - calo:  W cala e non risale più (righe anomale: quelle dopo il calo, fino a fine tabella)
    Sono segnalate a parte anche le righe con derivata W_s < 0, anche senza calo di W oltre la tolleranza.

    Output in output/<cartella>/W_check/:
      W_check.csv         un evento per riga
      report_W_check.md   elenco delle sezioni con descrizione breve degli eventi
      sezioni/<nome>.png  grafici come in grafici/sezioni, con le righe anomale in arancione

    Uso:  julia --project=. src/W_check.jl [<cartella di output>]   (predefinita: Brenta_caldonazzo_borgo)
=#
include(joinpath(@__DIR__, "grafici.jl"))     # carica(), righe(), plot_sezione(), costanti e cartelle

const TolAbs = 0.10     # [m], calo minimo di W considerato
const TolRel = 0.02     # [-], calo minimo di W relativo a W
const CheckFolder = joinpath(Folder, "W_check")

tol(w) = max(TolAbs, TolRel * w)

# Eventi di una sezione: vettore di NamedTuple con tipo, righe anomale e valori di W
function eventi_W(t)
    W, η = t.W, t.WSE
    n = length(W)
    ev = []
    for j in 1:n-1
        W[j] - W[j+1] > tol(W[j]) || continue
        i = findlast(k -> W[k] <= W[j+1] + tol(W[j+1]), 1:j-1)
        m = findfirst(k -> W[k] >= W[j] - tol(W[j]), j+2:n)
        if i !== nothing
            tipo, r = "picco", i+1:j
            Wp, Wpk, Wd = W[i], maximum(W[r]), W[j+1]
            descr = @sprintf("W sale da %.2f a %.2f m tra WSE %.3f e %.3f e torna a %.2f m a WSE %.3f (spessore %.3f m)",
                Wp, Wpk, η[i], η[j], Wd, η[j+1], η[j+1] - η[i])
        elseif m !== nothing
            m = j + 1 + m
            tipo, r = "buco", j+1:m-1
            Wp, Wpk, Wd = W[j], minimum(W[r]), W[m]
            descr = @sprintf("W scende da %.2f a %.2f m tra WSE %.3f e %.3f e risale a %.2f m a WSE %.3f (spessore %.3f m)",
                Wp, Wpk, η[j], η[j+1], Wd, η[m], η[m] - η[j])
        else
            tipo, r = "calo", j+1:n
            Wp, Wpk, Wd = W[j], W[j+1], W[end]
            descr = @sprintf("W scende da %.2f a %.2f m tra WSE %.3f e %.3f e non risale (%.2f m a fine tabella)",
                Wp, Wpk, η[j], η[j+1], Wd)
        end
        push!(ev, (tipo=tipo, righe=r, WSE_da=η[first(r)], WSE_a=η[last(r)], W_prima=Wp, W_anomala=Wpk,
            W_dopo=Wd, ampiezza=abs(Wpk - Wd), W_s_min=minimum(t.W_s[max(first(r)-1, 1):min(last(r)+1, n)]),
            descr=descr))
    end
    # derivata negativa in righe non già coperte da un evento
    coperte = isempty(ev) ? Int[] : reduce(vcat, collect.(getfield.(ev, :righe)))
    neg = setdiff(findall(<(0), t.W_s), coperte)
    if !isempty(neg)
        push!(ev, (tipo="W_s<0", righe=neg, WSE_da=η[first(neg)], WSE_a=η[last(neg)], W_prima=NaN,
            W_anomala=NaN, W_dopo=NaN, ampiezza=0.0, W_s_min=minimum(t.W_s[neg]),
            descr=@sprintf("derivata W_s negativa in %d righe (minimo %.2f) senza calo di W oltre la tolleranza",
                length(neg), minimum(t.W_s[neg]))))
    end
    return ev
end

#= Controllo e scrittura di report, CSV e grafici in `cartella`.
    sezioni = nothing: tutte le sezioni, con report e grafico solo per quelle con anomalie.
    sezioni = vettore di nomi: solo quelle, sempre con report e grafico (anche senza anomalie) e con la
              tabella di W riga per riga (usato da W_check_sezioni.jl).
=#
function check(; sezioni=nothing, cartella=CheckFolder)
    geo, tab = carica()
    scelte = sezioni !== nothing
    if scelte
        ignote = setdiff(sezioni, geo.nome)
        isempty(ignote) || error("Sezioni non presenti in $(joinpath(Folder, "topologia_sezioni.csv")): " *
            join(ignote, ", ") * "\n  Il nome va scritto come in topologia_sezioni.csv, es. $(first(geo.nome))")
        geo = geo[in.(geo.nome, Ref(Set(sezioni))), :]
    end
    rm(joinpath(cartella, "sezioni"); force=true, recursive=true)   # niente grafici di sezioni vecchie
    mkpath(joinpath(cartella, "sezioni"))
    out = DataFrame(section=Int[], nome=String[], asta=String[], tipo=String[], righe=String[],
        WSE_da=Float64[], WSE_a=Float64[], W_prima=Float64[], W_anomala=Float64[], W_dopo=Float64[],
        ampiezza=Float64[], W_s_min=Float64[], descrizione=String[])
    md = ["# Controllo della larghezza W - $OutputName", "",
        @sprintf("Tolleranza sul calo di W: max(%.2f m, %.0f %% di W). Sezioni controllate: %d.", TolAbs, 100TolRel, nrow(geo)), ""]
    nsez = 0
    for gs in eachrow(sort(geo, :section))
        t = righe(tab, gs.section)
        ev = eventi_W(t)
        if isempty(ev)
            scelte || continue
            push!(md, "## $(gs.nome) (sezione $(gs.section), x = $(@sprintf("%.0f", gs["x_coord_[m]"])) m): nessuna anomalia", "")
            append!(md, tabella_W(t, Int[]))
            push!(md, "", "![$(gs.nome)](sezioni/$(gs.nome).png)", "")
            savefig(plot_sezione(t, gs, "controllo W: nessuna anomalia"), joinpath(cartella, "sezioni", gs.nome*".png"))
            continue
        end
        nsez += 1
        for e in ev
            push!(out, (gs.section, gs.nome, string(gs.asta), e.tipo, "$(first(e.righe))-$(last(e.righe))",
                e.WSE_da, e.WSE_a, e.W_prima, e.W_anomala, e.W_dopo, e.ampiezza, e.W_s_min, e.descr))
        end
        tipi = join(unique(getfield.(ev, :tipo)), ", ")
        push!(md, "## $(gs.nome) (sezione $(gs.section), x = $(@sprintf("%.0f", gs["x_coord_[m]"])) m): $tipi", "")
        append!(md, ["- **$(e.tipo)** (righe $(first(e.righe))-$(last(e.righe))): $(e.descr); W_s minima $(@sprintf("%.1f", e.W_s_min))" for e in ev])
        evid = sort(unique(reduce(vcat, collect.(getfield.(ev, :righe)))))
        scelte && (push!(md, ""); append!(md, tabella_W(t, evid)))
        push!(md, "", "![$(gs.nome)](sezioni/$(gs.nome).png)", "")
        savefig(plot_sezione(t, gs, "controllo W: $tipi"; evid=evid), joinpath(cartella, "sezioni", gs.nome*".png"))
    end
    insert!(md, 4, "Sezioni con anomalie: **$nsez** ($(count(==("picco"), out.tipo)) picchi, " *
        "$(count(==("buco"), out.tipo)) buchi, $(count(==("calo"), out.tipo)) cali, " *
        "$(count(==("W_s<0"), out.tipo)) con sola derivata negativa).\n")
    CSV.write(joinpath(cartella, "W_check.csv"), out)
    write(joinpath(cartella, "report_W_check.md"), join(md, "\n"))
    println("Controllo W: $nsez sezioni con anomalie su $(nrow(geo)), report in $cartella")
    show(select(out, :nome, :tipo, :righe, :WSE_da, :WSE_a, :W_prima, :W_anomala, :W_dopo, :W_s_min),
        allrows=true, summary=false, eltypes=false)
    println()
    return out
end

# Tabella markdown di W riga per riga, con le righe anomale in grassetto
function tabella_W(t, evid)
    md = ["| riga | WSE [m] | W [m] | W_s [-] |", "|---|---|---|---|"]
    for k in 1:nrow(t)
        c = @sprintf("%d | %.3f | %.3f | %.2f", k, t.WSE[k], t.W[k], t.W_s[k])
        push!(md, k in evid ? "| **" * replace(c, " | " => "** | **") * "** |" : "| $c |")
    end
    return md
end

# eseguito solo se lanciato direttamente: W_check_sezioni.jl include questo file
if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 1 || error("Uso: julia --project=. src/W_check.jl [<cartella di output>]")
    check()
end
