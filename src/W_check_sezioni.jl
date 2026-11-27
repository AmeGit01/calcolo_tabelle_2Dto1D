#= CONTROLLO DELLA LARGHEZZA W DI SEZIONI SCELTE
    Come W_check.jl, ma controlla solo le sezioni passate come argomenti invece di cercare da solo
    quelle con anomalie. Utile per ricontrollare le sezioni corrette a mano.
    Per ogni sezione scrive sempre il grafico e la descrizione nel report, anche se non ha anomalie,
    insieme alla tabella di W riga per riga (righe anomale in grassetto).

    Output in output/<cartella>/W_check_sezioni/ (svuotata a ogni lancio, separata da quella di W_check.jl):
      W_check.csv         un evento per riga (solo sezioni con anomalie)
      report_W_check.md   tutte le sezioni richieste, con esito, tabella di W e grafico
      sezioni/<nome>.png  grafico di ogni sezione richiesta, righe anomale in arancione

    Uso:  julia --project=. src/W_check_sezioni.jl <cartella di output> <sezione> [<sezione> ...]
          es.  julia --project=. src/W_check_sezioni.jl Brenta_completo B13100 B32300
    Le sezioni si indicano con il nome (come in topologia_sezioni.csv).
=#
if length(ARGS) < 2 || any(in(("-h", "--help")), ARGS)
    println("Uso: julia --project=. src/W_check_sezioni.jl <cartella di output> <sezione> [<sezione> ...]\n" *
            "  es. julia --project=. src/W_check_sezioni.jl Brenta_completo B13100 B32300")
    exit(length(ARGS) < 2 ? 1 : 0)
end

include(joinpath(@__DIR__, "W_check.jl"))    # eventi_W(), check(), tolleranze; legge la cartella da ARGS[1]

check(; sezioni=unique(String.(ARGS[2:end])), cartella=joinpath(Folder, "W_check_sezioni"))
