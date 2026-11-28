# Independent phase-count and enthalpy check for issue #243.
include("generate_vlle.jl")

VERSION == v"1.11.7" || error("PH references pin Julia 1.11.7")
pkgversion(Clapeyron) == v"0.6.30" || error("PH references pin Clapeyron 0.6.30")

function check_ph_bracket()
    ids = ["nitrogen", "methanol", "water"]
    z = [0.628421052631579, 0.003157894736842105, 0.3684210526315789]
    p = 1e6
    for variant in ("PR", "PR-BM")
        model = matched_model(ids, variant, [zeros(3) for _ in ids])
        for (T, expected_count) in ((200.0, 3), (300.0, 2))
            beta, x, phases = if T == 200.0
                r = Clapeyron.tp_flash2(model, p, T, z, MultiPhaseTPFlash(full_tpd=true))
                pr_candidate(model, p, T, z, Float64.(r.fractions),
                             [Float64.(w) for w in r.compositions], Float64.(r.volumes))
            else
                # Independent two-phase candidate: most nitrogen in the vapor
                # and nearly all water in the liquid are only initial guesses.
                vapor_seed = z .* [0.99, 0.8, 0.005]
                roots = [:vapour, :liquid]
                amounts, compositions = pr_solve(model, p, T, z, roots,
                    permutedims(hcat(vapor_seed, z - vapor_seed)))
                amounts, compositions, roots
            end
            @assert length(beta) == expected_count (variant, T, beta)
            @assert all(>(0.0), beta) && isapprox(sum(beta), 1.0; atol=1e-12)
            ref = argmax(beta)
            d = pr_lnf(model, p, T, x[ref], phases[ref])
            @assert isempty(negative_minima(pr_grid(model, p, T, d); limit=1))
            for k in eachindex(beta)
                @assert isapprox(volume(model, p, T, x[k]), volume(model, p, T, x[k]; phase=phases[k]); rtol=1e-9)
            end
            h = sum(beta[k] * phase_properties(model, p, T, x[k], phases[k])["h"] for k in eachindex(beta))
            println(variant, " at ", T, " K: ", length(beta), " phases, h=", h, " J/mol, fractions ", beta)
        end
    end
end

check_ph_bracket()
