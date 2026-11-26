# Usa sempre l'ambiente Julia del progetto (Project.toml), anche se lo script è lanciato senza --project=.
import Pkg
let progetto = normpath(joinpath(@__DIR__, ".."))
    normpath(dirname(Base.active_project())) == progetto || Pkg.activate(progetto; io=devnull)
end
ENV["GKSwstype"] = "100"    # GR senza finestra: i grafici sono solo salvati
using Plots, Plots.Measures, Printf, CSV, DataFrames

#= GRAFICI DI CONTROLLO DELLE TABELLE
    Legge SOLO i file già salvati in output/<cartella> (da main_sezioni.jl e main_diversione.jl):
        localslope_xcoord.csv, topologia_sezioni.csv, tabella_sezioni.csv, scala_diversione.csv
    e salva i grafici in output/<cartella>/grafici:
        profilo_<asta>.png          talweg, argini e livello massimo in tabella lungo l'asta
        pendenza_<asta>.png         pendenza locale lungo l'asta
        controlli_<asta>.png        grandezze lungo l'asta utili a trovare sezioni anomale
        mappa_jden.png, mappa_DQ.png  mappe sezione × riga della tabella
        scala_diversione.png        portata verso il Brenta vecchio in funzione di quella da monte
        sezioni/<nome>.png          curve di una sezione in funzione del livello, per le sezioni scelte
        sezioni_selezionate.csv     sezioni scelte e motivo
    Le sezioni con grafico dedicato sono: prima e ultima di ogni asta, quelle della diversione e le
    NSezioniAnomale con il punteggio di anomalia più alto (vedi anomalie()).

    Uso:  julia --project=. src/grafici.jl [<cartella di output>]   (predefinita: Brenta_caldonazzo_borgo)
=#
const OutputDefault = "Brenta_caldonazzo_borgo"
const SezioniDiversione = ["B13200", "B13250", "V00050"]   # lasciare vuoto se non c'è diversione
const NSezioniAnomale = 12
# quantile per interpolazione lineare (evita la dipendenza da Statistics)
function quantile(v, p)
    x = sort(collect(v)); h = 1 + (length(x) - 1) * p; i = floor(Int, h)
    return i >= length(x) ? x[end] : x[i] + (h - i) * (x[i+1] - x[i])
end
const ProjectFolder = normpath(joinpath(@__DIR__, ".."))

length(ARGS) <= 1 || error("Uso: julia --project=. src/grafici.jl [<cartella di output>]")
const OutputName = isempty(ARGS) ? OutputDefault : basename(rstrip(ARGS[1], '/'))
const Folder = joinpath(ProjectFolder, "output", OutputName)
const PlotFolder = joinpath(Folder, "grafici")

default(fontfamily="sans-serif", dpi=150, framestyle=:box, grid=true, gridalpha=0.25,
    titlefontsize=10, guidefontsize=9, tickfontsize=8, legendfontsize=8)
const C_talweg, C_dx, C_sx, C_wse, C_evid = "#5b4636", "#1f77b4", "#d62728", "#17becf", "#ff7f0e"

leggi(file) = (p = joinpath(Folder, file); isfile(p) ? CSV.read(p, DataFrame) : nothing)
salva(plt, nome) = (savefig(plt, joinpath(PlotFolder, nome)); println("  ", nome))

# Unisce geometria e topologia; con un file di sola topologia (una asta) l'asta vale "".
function carica()
    geo = leggi("localslope_xcoord.csv")
    tab = leggi("tabella_sezioni.csv")
    (geo === nothing || tab === nothing) && error("In $Folder servono localslope_xcoord.csv e tabella_sezioni.csv")
    top = leggi("topologia_sezioni.csv")
    if top === nothing
        top = DataFrame(section=geo.section, nome=string.(geo.section), asta=fill("", nrow(geo)))
    end
    geo = innerjoin(geo, select(top, :section, :nome, :asta), on=:section)
    geo.x_km = geo[!, "x_coord_[m]"] ./ 1000
    return geo, tab
end

# Righe della tabella di una sezione (ordinate per livello)
righe(tab, s) = sort(tab[tab.sezione .== s, :], :WSE)

#= Punteggio di anomalia per sezione: somma di indicatori normalizzati
    - righe in cui Ω, W non crescono col livello o Ω_s, W_s < 0 (geometria incoerente)
    - righe in cui RH cala bruscamente (salto > 20 % tra righe successive)
    - jden massimo rispetto alla mediana dell'asta (attrito anomalo)
    - derivata Ω_s lontana dalla differenza finita dΩ/dWSE
    - |ΔQ| massimo rispetto alla portata (scambi laterali grandi)
    - righe mancanti rispetto al massimo
=#
function anomalie(geo, tab)
    nmax = maximum(nrow(righe(tab, s)) for s in geo.section)
    out = DataFrame(section=Int[], geometria=Int[], salto_RH=Int[], jden_max=Float64[],
        err_Ωs=Float64[], ΔQ_max=Float64[], righe=Int[])
    for s in geo.section
        t = righe(tab, s)
        dΩ, dW, dRH = diff(t.Ω), diff(t.W), diff(t.RH)
        geom = count(<(0), dΩ) + count(<(0), dW) + count(<(0), t.Ω_s) + count(<(0), t.W_s)
        salto = count(i -> t.RH[i+1] < 0.8 * t.RH[i], 1:nrow(t)-1)
        fd = dΩ ./ max.(diff(t.WSE), 1e-6)
        ref = (t.Ω_s[1:end-1] .+ t.Ω_s[2:end]) ./ 2
        errΩ = isempty(fd) ? 0.0 : maximum(abs.(fd .- ref) ./ max.(abs.(ref), 1.0))
        push!(out, (s, geom, salto, maximum(t.jden), errΩ, maximum(abs.(t.ΔQ)), nmax - nrow(t)))
    end
    out = innerjoin(out, select(geo, :section, :nome, :asta), on=:section)
    out.jden_rel = similar(out.jden_max)
    for asta in unique(out.asta)
        idx = out.asta .== asta
        out.jden_rel[idx] .= out.jden_max[idx] ./ quantile(out.jden_max[idx], 0.5)
    end
    out.punteggio = out.geometria .+ out.salto_RH ./ 2 .+ max.(log10.(out.jden_rel), 0) .* 2 .+
        min.(out.err_Ωs, 5) .+ out.ΔQ_max ./ 10 .+ out.righe ./ 2
    return sort(out, :punteggio, rev=true)
end

function plot_profilo(g, tab, asta)
    lmax = [maximum(righe(tab, s).WSE) for s in g.section]
    plt = plot(g.x_km, g[!, "talweg_[m]"], lw=1.5, c=C_talweg, label="talweg",
        fillrange=minimum(g[!, "talweg_[m]"]) - 2, fillalpha=0.15, fillcolor=C_talweg,
        title="Profilo longitudinale - asta $asta", xlabel="x [km]", ylabel="quota [m s.l.m.]",
        size=(1300, 550), left_margin=5mm, bottom_margin=5mm, legend=:topright)
    plot!(plt, g.x_km, g[!, "argine_dx_[m]"], lw=1, c=C_dx, label="argine dx")
    plot!(plt, g.x_km, g[!, "argine_sx_[m]"], lw=1, c=C_sx, label="argine sx")
    plot!(plt, g.x_km, lmax, lw=1, ls=:dash, c=C_wse, label="livello massimo in tabella")
    # sezioni in cui il livello massimo supera l'argine più basso: lì la tabella va oltre l'alveo
    sopra = lmax .> min.(g[!, "argine_dx_[m]"], g[!, "argine_sx_[m]"])
    any(sopra) && scatter!(plt, g.x_km[sopra], lmax[sopra], ms=1.8, msw=0, c=C_evid,
        label="livello > argine più basso ($(count(sopra)))")
    div = intersect(SezioniDiversione, g.nome)
    if !isempty(div)
        x = g.x_km[findfirst(==(div[1]), g.nome)]
        vline!(plt, [x], c=:gray, ls=:dot, label="")
        annotate!(plt, x, maximum(g[!, "argine_sx_[m]"]), text(" diversione ("*join(div, ", ")*")", 7, :left, :gray))
    end
    return plt
end

function plot_pendenza(g, asta)
    i = g.local_slope
    plt = plot(g.x_km, i .* 1000, lw=1, marker=:circle, ms=1.5, msw=0, c=C_talweg, label="pendenza locale",
        title="Pendenza locale - asta $asta", xlabel="x [km]", ylabel="pendenza [‰]  (< 0: fondo che scende verso valle)",
        size=(1300, 400), left_margin=5mm, bottom_margin=5mm)
    hline!(plt, [0], c=:black, lw=0.5, label="")
    # segno del modello 1D: negativa se il fondo scende verso valle, positiva in contropendenza
    contro = i .> 0
    any(contro) && scatter!(plt, g.x_km[contro], i[contro] .* 1000, c=C_evid, ms=3, msw=0,
        label="contropendenza, > 0 ($(count(contro)))")
    return plt
end

function plot_controlli(g, tab, an, asta)
    a = an[in.(an.section, Ref(Set(g.section))), :]
    a = a[sortperm(indexin(a.section, g.section)), :]
    ultima = [righe(tab, s)[end, :] for s in g.section]
    tw = g[!, "talweg_[m]"]
    p1 = plot(g.x_km, [r.WSE for r in ultima] .- tw, c=C_wse, label="tirante massimo in tabella",
        ylabel="[m]", title="asta $asta")
    plot!(p1, g.x_km, min.(g[!, "argine_dx_[m]"], g[!, "argine_sx_[m]"]) .- tw, c=C_sx, label="argine più basso - talweg")
    p2 = plot(g.x_km, [r.W for r in ultima], c=C_dx, label="W al livello massimo", ylabel="W [m]")
    p3 = plot(g.x_km, a.jden_max, yscale=:log10, c=C_talweg, label="jden massimo", ylabel="jden [m]")
    p4 = plot(g.x_km, [maximum(righe(tab, s).ΔQ) for s in g.section], c=C_dx, label="ΔQ massimo", ylabel="ΔQ [m³/s]")
    plot!(p4, g.x_km, [minimum(righe(tab, s).ΔQ) for s in g.section], c=C_sx, label="ΔQ minimo")
    p5 = bar(g.x_km, a.punteggio, c=C_evid, lw=0, label="punteggio di anomalia", ylabel="[-]", xlabel="x [km]")
    return plot(p1, p2, p3, p4, p5, layout=(5, 1), size=(1300, 1300), link=:x, left_margin=8mm, legend=:topleft)
end

# Mappa sezione × riga della tabella (righe dal livello più basso al più alto)
function plot_mappa(geo, tab, col, titolo; trasf=identity, clims=:auto, cmap=:viridis)
    nmax = maximum(nrow(righe(tab, s)) for s in geo.section)
    M = fill(NaN, nmax, nrow(geo))
    for (j, s) in enumerate(geo.section)
        v = righe(tab, s)[!, col]
        M[1:length(v), j] .= trasf.(v)
    end
    plt = heatmap(1:nrow(geo), 1:nmax, M, c=cmap, clims=clims, title=titolo,
        xlabel="sezione (numero progressivo)", ylabel="riga della tabella (livello crescente)",
        size=(1300, 450), left_margin=5mm, bottom_margin=5mm)
    # confine fra le aste
    for k in findall(diff(geo.asta .== geo.asta[1]) .!= 0)
        vline!(plt, [k + 0.5], c=:white, lw=1.5, label="")
    end
    return plt
end

function plot_sezione(t, gs, motivo)
    tw, adx, asx = gs["talweg_[m]"], gs["argine_dx_[m]"], gs["argine_sx_[m]"]
    function pannello(x, xl; kw...)
        p = plot(x, t.WSE; marker=:circle, ms=2, msw=0, c=C_dx, label="", xlabel=xl, ylabel="WSE [m]", kw...)
        hline!(p, [tw], c=C_talweg, ls=:dash, label="talweg")
        hline!(p, [adx], c=C_dx, ls=:dot, label="argine dx")
        hline!(p, [asx], c=C_sx, ls=:dot, label="argine sx")
        return p
    end
    fd = diff(t.Ω) ./ max.(diff(t.WSE), 1e-6)
    wm = (t.WSE[1:end-1] .+ t.WSE[2:end]) ./ 2
    p1 = pannello(t.Ω, "Ω [m²]"; legend=:topleft)
    p2 = pannello(t.W, "W [m]"; legend=false)
    p3 = pannello(t.RH, "RH [m]"; legend=false)
    p4 = pannello(t.jden, "jden [m]"; legend=false, xscale=:log10)
    p5 = pannello(t.ΔQ, "ΔQ [m³/s]"; legend=false)
    vline!(p5, [0], c=:black, lw=0.5, label="")
    p6 = pannello(t.Ω_s, "dΩ/dWSE [m]"; legend=false)
    plot!(p6, fd, wm, c=C_evid, ls=:dash, marker=:x, ms=3, label="differenza finita", legend=:bottomright)
    plot(p1, p2, p3, p4, p5, p6, layout=(2, 3), size=(1200, 800), left_margin=5mm, bottom_margin=5mm,
        plot_title="$(gs.nome) (sezione $(gs.section), x = $(@sprintf("%.0f", gs["x_coord_[m]"])) m) - $motivo",
        plot_titlefontsize=11)
end

function plot_diversione(sc)
    qm, qv, fr = sc[!, "Q_monte_[m3/s]"], sc[!, "Q_brenta_vecchio_[m3/s]"], sc[!, "frazione_brenta_vecchio_[-]"]
    p1 = plot(qm, qv, marker=:circle, ms=3, c=C_dx, label="Q Brenta vecchio",
        xlabel="Q monte [m³/s]", ylabel="Q [m³/s]", title="Scala di deflusso della diversione", legend=:topleft)
    plot!(p1, qm, sc[!, "Q_brenta_nuovo_[m3/s]"], marker=:circle, ms=3, c=C_sx, label="Q Brenta nuovo")
    p2 = plot(qm[qm .> 0], 100 .* fr[qm .> 0], marker=:circle, ms=3, c=C_dx, label="",
        xlabel="Q monte [m³/s]", ylabel="quota verso il Brenta vecchio [%]", title="Ripartizione")
    p3 = plot(qm, sc[!, "residuo_bilancio_[m3/s]"], marker=:circle, ms=3, c=C_talweg, label="",
        xlabel="Q monte [m³/s]", ylabel="Q monte - Q nuovo - Q vecchio [m³/s]", title="Bilancio")
    hline!(p3, [0], c=:black, lw=0.5, label="")
    plot(p1, p2, p3, layout=(1, 3), size=(1400, 450), left_margin=6mm, bottom_margin=6mm)
end

function main()
    geo, tab = carica()
    mkpath(joinpath(PlotFolder, "sezioni"))
    println("Grafici in $PlotFolder")
    an = anomalie(geo, tab)
    aste = unique(geo.asta)
    for asta in aste
        g = sort(geo[geo.asta .== asta, :], :section)
        suff = asta == "" ? "" : "_"*asta
        salva(plot_profilo(g, tab, asta), "profilo$suff.png")
        salva(plot_pendenza(g, asta), "pendenza$suff.png")
        salva(plot_controlli(g, tab, an, asta), "controlli$suff.png")
    end
    salva(plot_mappa(geo, tab, :jden, "log10(jden [m])"; trasf=log10), "mappa_jden.png")
    qlim = quantile(abs.(tab.ΔQ), 0.99)
    salva(plot_mappa(geo, tab, :ΔQ, "ΔQ [m³/s] (scala limitata al 99° percentile di |ΔQ|)";
        clims=(-qlim, qlim), cmap=:RdBu), "mappa_DQ.png")
    sc = leggi("scala_diversione.csv")
    sc === nothing || salva(plot_diversione(sc), "scala_diversione.png")

    # sezioni con grafico dedicato
    scelte = Dict{String,Vector{String}}()
    aggiungi(n, m) = push!(get!(scelte, n, String[]), m)
    for asta in aste
        g = sort(geo[geo.asta .== asta, :], :section)
        aggiungi(g.nome[1], "prima dell'asta"); aggiungi(g.nome[end], "ultima dell'asta")
    end
    foreach(n -> n in geo.nome && aggiungi(n, "diversione"), SezioniDiversione)
    for r in eachrow(first(an, NSezioniAnomale))
        aggiungi(r.nome, @sprintf("anomalia %.1f", r.punteggio))
    end
    elenco = DataFrame(section=Int[], nome=String[], motivo=String[])
    for (n, motivi) in scelte
        i = findfirst(==(n), geo.nome)
        push!(elenco, (geo.section[i], n, join(motivi, "; ")))
    end
    sort!(elenco, :section)
    for r in eachrow(elenco)
        salva(plot_sezione(righe(tab, r.section), geo[findfirst(==(r.section), geo.section), :], r.motivo),
            joinpath("sezioni", r.nome*".png"))
    end
    CSV.write(joinpath(PlotFolder, "sezioni_selezionate.csv"),
        innerjoin(elenco, select(an, Not([:nome, :asta])), on=:section))
    println("\nSezioni più anomale:")
    show(first(select(an, :nome, :punteggio, :geometria, :salto_RH, :jden_rel, :err_Ωs, :ΔQ_max, :righe), NSezioniAnomale),
        allrows=true, allcols=true, summary=false, eltypes=false)
    println()
end

main()
