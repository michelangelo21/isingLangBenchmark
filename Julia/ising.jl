using Statistics

function ising(L::Int, T::Number, next_I, prev_I; MCS=230_000::Int)
    lattice = ones(Int8, L,L)
    
    boltzmann_factor = map(ΔE -> exp(-ΔE / T), -8:4:8)

    #initial arrays (returned)
    magnetizations = zeros(0)
    energies = zeros(0)

    for k=1 : MCS
        for i=1:L, j=1:L
        #metropolis algorithm
            @inbounds ΔE::Int = 2 * lattice[i,j] * (lattice[next_I[i],j] + lattice[prev_I[i],j] + lattice[i,next_I[j]] + lattice[i,prev_I[j]])
            if (ΔE ≤ 0) || (rand() ≤ boltzmann_factor[ΔE÷4 + 3])
                @inbounds lattice[i,j] *= -1
            end
        end

        if (30_000 < k && k%100 == 0)
            m_k::Float64 = abs(mean(lattice))
            push!(magnetizations, m_k)
            E::Float64 = 0.0
            for i=1:L, j=1:L
                # E += 0.5 * lattice[i,j] * (lattice[next_I[i],j] + lattice[prev_I[i],j] + lattice[i,next_I[j]] + lattice[i,prev_I[j]])
                @inbounds E += lattice[i,j] * (lattice[next_I[i],j] + lattice[i,next_I[j]])
            end
            push!(energies, E)
        end
    end
    return magnetizations, energies
end

function run_MC(Ls::Array{Int}, Ts)
    result = Dict{Int, Dict{Char, Dict{Float64,Float64}}}(l => Dict(c => Dict() for c in ['m','χ','C']) for l in Ls)
    for L ∈ Ls
        next_I = collect(2:L+1)
        prev_I = collect(0:L-1)
        next_I[end] = 1 
        prev_I[1] = L
        Threads.@threads for T ∈ Ts
            magnetizations, energies = ising(L,T,next_I, prev_I)

            result[L]['m'][T] = mean(magnetizations)                                                    # mean magnetization        <m>
            result[L]['χ'][T] = L^2 / T * var(magnetizations, mean=result[L]['m'][T], corrected=false)  # magnetic susceptibility    χ
            result[L]['C'][T] = var(energies, corrected=false)   / L^2 / T^2                            # heat capacity              C
        end
    end
    return result
end

function main()
    Ls = [10]
    # Ts = cat(0.5:0.05:1.75, 1.8:0.005:2.7, 2.75:0.05:3.5, dims=1)
    Ts = 0.5:0.1:3.5
    result_ret = @time run_MC(Ls,Ts)
end
main()