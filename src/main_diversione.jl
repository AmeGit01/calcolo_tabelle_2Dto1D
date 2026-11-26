# Usa sempre l'ambiente Julia del progetto (Project.toml), anche se lo script è lanciato senza --project=.
import Pkg
let progetto = normpath(joinpath(@__DIR__, ".."))
    normpath(dirname(Base.active_project())) == progetto || Pkg.activate(progetto; io=devnull)
end
using Printf, CSV, DataFrames

#= SCALA DI DEFLUSSO DELLA DIVERSIONE BRENTA NUOVO / BRENTA VECCHIO
    Il Brenta è a canale singolo fino alla sezione B13200, dove si divide nel Brenta nuovo (B13250, il
    canale principale) e nel Brenta vecchio (V00050, prima sezione disponibile dell'asta V). Nel
    modello 1D le due aste sono risolte in modo indipendente: il Brenta nuovo ha una portata uscente
    verso il Brenta vecchio, che qui è ricavata dalla simulazione 2D a scalini.
    Per ogni salvataggio della simulazione (una riga dei file delle sezioni) legge la portata in arrivo
    da monte (B13200) e quella che va nel Brenta vecchio (V00050), e scrive la scala di deflusso
    portata-portata in output/<cartella>/scala_diversione.csv.

    Uso:  julia --project=. src/main_diversione.jl [<cartella di input>]
          (predefinita: Brenta_caldonazzo_borgo; l'output va nella sottocartella di output/ con lo
          stesso nome)
=#
const InputDefault = "Brenta_caldonazzo_borgo"
const SezMonte = "B13200"      # ultima sezione a canale singolo, a monte della divergenza
const SezNuovo = "B13250"      # prima sezione del Brenta nuovo, usata solo per il bilancio
const SezVecchio = "V00050"    # prima sezione del Brenta vecchio con risultati (V00000 non c'è)
# scarto massimo del bilancio Q_monte - Q_nuovo - Q_vecchio, relativo a Q_monte, oltre il quale si avvisa
const Bilancio_lim = 0.02
const Format = ".csv"
const ProjectFolder = normpath(joinpath(@__DIR__, ".."))

length(ARGS) <= 1 || error("Uso: julia --project=. src/main_diversione.jl [<cartella di input>]")
const InputName = isempty(ARGS) ? InputDefault : basename(rstrip(ARGS[1], '/'))
const InputFolder = joinpath(ProjectFolder, "input", InputName)
const OutputFolder = joinpath(ProjectFolder, "output", InputName)

# Portate di una sezione, una per salvataggio, dal file <prefisso><nome>.csv della cartella di input
function read_discharge(nome)
    trovati = filter(f -> endswith(f, nome*Format), readdir(InputFolder))
    length(trovati) == 1 || error("In $InputFolder serve un solo file *$nome$Format, trovati $(length(trovati))")
    data = CSV.read(joinpath(InputFolder, only(trovati)), DataFrame)
    "Discharge [m3/s]" in names(data) || error("Nel file $(only(trovati)) manca la colonna \"Discharge [m3/s]\"")
    return Float64.(data[!, "Discharge [m3/s]"])
end

function main(; dosave::Bool=true)
    Q_monte, Q_nuovo, Q_vecchio = read_discharge.((SezMonte, SezNuovo, SezVecchio))
    length(Q_monte) == length(Q_nuovo) == length(Q_vecchio) ||
        error("Le sezioni $SezMonte, $SezNuovo e $SezVecchio hanno un numero di salvataggi diverso")

    residuo = Q_monte .- Q_nuovo .- Q_vecchio
    avvisi = String[]
    for k in eachindex(Q_monte)
        Q_monte[k] > 0 && abs(residuo[k]) > Bilancio_lim * Q_monte[k] && push!(avvisi,
            @sprintf("salvataggio %3d  Q_monte=%8.3f  bilancio non chiuso: Q_monte - Q_nuovo - Q_vecchio = %.3f", k, Q_monte[k], residuo[k]))
    end

    # la scala di deflusso deve avere Q_monte crescente: scarta i salvataggi in cui non cresce
    tieni = [k == 1 || Q_monte[k] > maximum(Q_monte[1:k-1]) for k in eachindex(Q_monte)]
    for k in findall(.!tieni)
        push!(avvisi, @sprintf("salvataggio %3d  Q_monte=%8.3f  portata non crescente: riga scartata", k, Q_monte[k]))
    end

    frazione = [q > 0 ? v / q : 0.0 for (q, v) in zip(Q_monte, Q_vecchio)]
    scala = DataFrame("salvataggio" => findall(tieni),
        "Q_monte_[m3/s]" => Q_monte[tieni], "Q_brenta_vecchio_[m3/s]" => Q_vecchio[tieni],
        "Q_brenta_nuovo_[m3/s]" => Q_nuovo[tieni], "frazione_brenta_vecchio_[-]" => frazione[tieni],
        "residuo_bilancio_[m3/s]" => residuo[tieni])

    println("Diversione: monte $SezMonte, Brenta nuovo $SezNuovo, Brenta vecchio $SezVecchio")
    @printf("%5s %12s %12s %12s %10s %10s\n", "salv.", "Q_monte", "Q_vecchio", "Q_nuovo", "frazione", "residuo")
    for r in eachrow(scala)
        @printf("%5d %12.3f %12.3f %12.3f %10.4f %10.3f\n", r[1], r[2], r[3], r[4], r[5], r[6])
    end
    println("\nControlli: $(length(avvisi)) avvisi")
    foreach(println, avvisi)

    if dosave
        mkpath(OutputFolder)
        path = joinpath(OutputFolder, "scala_diversione"*Format)
        CSV.write(path, scala)
        println("Scala di deflusso scritta in $path")
    end
    return scala
end

main()
