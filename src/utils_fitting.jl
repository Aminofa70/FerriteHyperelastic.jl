function calc_init_shear_mod(initial_mu, initial_alpha)
    shearMod = 0.5 * (initial_alpha .^ 2)' * (2 .* initial_mu)
    return shearMod
end

######################################################

function calc_mat_constants_hyper_biaxial(modelType::AbstractString,
                                          strainExp::AbstractVector,
                                          Sexp::AbstractVector)

    # detect zero stresses and filter
    detected_non_zero = Sexp .!= 0
    Sexp      = Sexp[detected_non_zero]
    strainExp = strainExp[detected_non_zero]

    # weights W = 1 ./ Sexp  (FLOATING division)
    W = 1.0 ./ Sexp

    λ = 1.0 .+ strainExp

    x = nothing
    if modelType == "neo-hookean"
        x = zeros(length(strainExp), 1)
        for i in eachindex(strainExp)
            xx = λ[i]
            x[i, 1] = 2 * (xx - 1 / xx^5)
        end

    elseif modelType == "mooney-rivlin"
        x = zeros(length(strainExp), 2)
        for i in eachindex(strainExp)
            xx = λ[i]
            x[i, 1] = 2 * (xx - 1 / xx^5)
            x[i, 2] = 2 * (xx^3 - 1 / xx^3)
        end

    elseif modelType == "yeoh"
        x = zeros(length(strainExp), 3)
        for i in eachindex(strainExp)
            xx = λ[i]
            Ibar1 = 2 * xx^2 + (1.0 / xx^2)^2
            x[i, 1] = 2 * (xx - 1 / xx^5)
            x[i, 2] = 4 * (xx - 1 / xx^5) * (Ibar1 - 3)
            x[i, 3] = 6 * (xx - 1 / xx^5) * (Ibar1 - 3)^2
        end

    else
        error("Unknown modelType: $modelType (expected 'neo-hookean', 'mooney-rivlin', or 'yeoh').")
    end

    xnew  = Diagonal(W) * x
    SeNew = Diagonal(W) * Sexp

    # Solve normal equations
    mat_cons = (xnew' * xnew) \ (xnew' * SeNew)

    return mat_cons
end

######################################################
function calc_mat_constants_hyper_ogden_linear(data_type::AbstractString,
    strainExp::AbstractVector,
    Sexp::AbstractVector)

    # pick c based on data_type
    c = data_type == "uniaxial" ? -0.5 :
        data_type == "biaxial" ? -2.0 :
        data_type == "shear" ? -1.0 :
        error("Unknown data_type: $data_type (expected 'uniaxial', 'biaxial', or 'shear').")

    initial_alpha = [2.0; 4.0; -2.0]
    maxIter = 100
    n = length(strainExp)
    λ = 1 .+ strainExp
    initial_mu = ones(3)
    W = 1.0 ./ Sexp          # weights

    x = zeros(n, 3)

    # inner function (captures c)
    function calc_stress(lambda::AbstractVector, mu::AbstractVector, alpha::AbstractVector)
        T = zeros(length(lambda))
        for p in eachindex(lambda)
            lam = lambda[p]
            s = 0.0
            @inbounds for o in 1:3
                mu_T = mu[o]
                alpha_T = alpha[o]
                s += (2 * mu_T / alpha_T) * (lam^(alpha_T - 1) - lam^(c * alpha_T - 1))
            end
            T[p] = s
        end
        return T
    end

    # Track best params seen (in case no improvement triggers early break)
    best_params = copy(initial_mu)
    best_err = Inf

    for iter in 1:maxIter
        T = calc_stress(λ, initial_mu, initial_alpha)
        residual = (1 .- T ./ Sexp) .^ 2
        current_error_norm = norm(residual)

        # build linear system matrix x
        for k in 1:n
            xx = λ[k]
            @inbounds for j in 1:3
                mu = initial_mu[j]
                alpha = initial_alpha[j]
                x[k, j] = (2 * mu / alpha) * (xx^(alpha - 1) - xx^(c * alpha - 1))
            end
        end

        # weighted least squares using normal equations (matches MATLAB)
        xnew = Diagonal(W) * x
        SeNew = Diagonal(W) * Sexp
        params = (xnew' * xnew) \ (xnew' * SeNew)

        shearMod = calc_init_shear_mod(params, initial_alpha)

        bad = any((params .== 0) .| isnan.(params)) ||
              any(params .<= 0) ||
              !(shearMod > 0)

        if bad
            # Change initial guess for the third alpha (random in [-6,-1])
            initial_alpha = [2.0; 4.0; float(rand(-6:-1))]
            continue
        end

        T2 = calc_stress(λ, params, initial_alpha)
        residual2 = (1 .- T2 ./ Sexp) .^ 2

        if norm(residual2) < current_error_norm
            # accept and return
            return params, initial_alpha
        else
            # keep searching; remember best so far
            if current_error_norm < best_err
                best_err = current_error_norm
                best_params = copy(params)
            end
            # "Changing the initial guess"
            initial_alpha = [2.0; 4.0; float(rand(-6:-1))]
        end
    end

    # If we reach here, no early break—return the best we saw
    return best_params, initial_alpha
end



function calc_residual_and_jacobian(data_type::AbstractString,
    strainExp::AbstractVector,
    Sexp::AbstractVector,
    params::AbstractVector,
    normalization_factor)

    # set c by data type
    c = data_type == "uniaxial" ? -0.5 :
        data_type == "biaxial" ? -2.0 :
        data_type == "shear" ? -1.0 :
        error("Unknown data_type: $data_type (expected 'uniaxial', 'biaxial', or 'shear').")

    n = length(strainExp)
    jacobian = zeros(eltype(float(one(eltype(strainExp)))), n, 6)
    residual = zeros(eltype(jacobian), n)

    λ = 1 .+ strainExp

    for k in 1:n
        lam = λ[k]
        L = log(lam)

        μ1, a1 = params[1], params[2]
        μ2, a2 = params[3], params[4]
        μ3, a3 = params[5], params[6]

        # derivatives wrt μ and α for each pair
        dSdmu1 = (2 / a1) * (lam^(a1 - 1) - lam^(c * a1 - 1))
        dSda1 = (2 * μ1 * (lam^(a1 - 1) * L - lam^(a1 * c - 1) * c * L)) / a1 +
                (2 * μ1 * (lam^(a1 * c - 1) - lam^(a1 - 1))) / a1^2

        dSdmu2 = (2 / a2) * (lam^(a2 - 1) - lam^(c * a2 - 1))
        dSda2 = (2 * μ2 * (lam^(a2 - 1) * L - lam^(a2 * c - 1) * c * L)) / a2 +
                (2 * μ2 * (lam^(a2 * c - 1) - lam^(a2 - 1))) / a2^2

        dSdmu3 = (2 / a3) * (lam^(a3 - 1) - lam^(c * a3 - 1))
        dSda3 = (2 * μ3 * (lam^(a3 - 1) * L - lam^(a3 * c - 1) * c * L)) / a3 +
                (2 * μ3 * (lam^(a3 * c - 1) - lam^(a3 - 1))) / a3^2

        jacobian[k, 1] = (-1 / Sexp[k]) * dSdmu1
        jacobian[k, 2] = (-1 / Sexp[k]) * dSda1
        jacobian[k, 3] = (-1 / Sexp[k]) * dSdmu2
        jacobian[k, 4] = (-1 / Sexp[k]) * dSda2
        jacobian[k, 5] = (-1 / Sexp[k]) * dSdmu3
        jacobian[k, 6] = (-1 / Sexp[k]) * dSda3
        # If you want the normalization-factor variant used in your comments:
        # jacobian[k,1] = -dSdmu1 / normalization_factor
        # jacobian[k,2] = -dSda1  / normalization_factor
        # jacobian[k,3] = -dSdmu2 / normalization_factor
        # jacobian[k,4] = -dSda2  / normalization_factor
        # jacobian[k,5] = -dSdmu3 / normalization_factor
        # jacobian[k,6] = -dSda3  / normalization_factor

        T1 = (2 * μ1 / a1) * (lam^(a1 - 1) - lam^(c * a1 - 1))
        T2 = (2 * μ2 / a2) * (lam^(a2 - 1) - lam^(c * a2 - 1))
        T3 = (2 * μ3 / a3) * (lam^(a3 - 1) - lam^(c * a3 - 1))
        Ttotal = T1 + T2 + T3

        residual[k] = (Sexp[k] - Ttotal) / Sexp[k]
        # or, using your commented line:
        # residual[k] = (Sexp[k] - Ttotal) / normalization_factor
    end

    return residual, jacobian
end

######################################################
function calc_mat_constants_hyper_ogden(dataType::AbstractString,
    strainExp::AbstractVector,
    Sexp::AbstractVector)

    # --- detect zero stress & filter ---
    keep = Sexp .!= 0
    Sexp = Sexp[keep]
    strainExp = strainExp[keep]

    # normalization factor (safe even if Sexp is empty or all ~0)
    normalization_factor =
        isempty(Sexp) ? 1.0 :
        max(maximum(abs.(Sexp)), 1e-12)  # avoid exactly 0; MATLAB used 1.0 as safety
    if normalization_factor == 0
        normalization_factor = 1.0
    end

    # --- initial guess from linear fit ---
    mu_params_linear, initial_alpha_linear =
        calc_mat_constants_hyper_ogden_linear(dataType, strainExp, Sexp)

    params = zeros(6)
    params[1:2:end] .= mu_params_linear
    params[2:2:end] .= initial_alpha_linear

    # --- restart settings (placeholders for future perturbations) ---
    params_base = copy(params)
    perturbation_factor = 5 / 100
    perturbation_rate = 3
    max_restarts = 3
    restart_count = 0

    # --- LM settings ---
    toler = 1e-3
    damping_increase = 10.0
    damping_decrease = 3.0
    maxIter = 1000

    fit_successful = false

    while !fit_successful && restart_count < max_restarts
        fit_stalled = false
        restart_count += 1

        if restart_count == 1
            params .= params_base
        else
            # future: perturbation = perturbation_factor .* params_base .* randn(length(params_base))
            perturbation = 0.0
            params .= params_base .+ perturbation
        end

        _, J0 = calc_residual_and_jacobian(dataType, strainExp, Sexp, params, normalization_factor)
        JtJ0 = J0' * J0
        τ = 1e-6 * maximum(diag(JtJ0))
        if !(τ > 0)   # covers τ == 0 and NaN
            τ = 1e-6
        end

        for _ in 1:maxIter
            r, J = calc_residual_and_jacobian(dataType, strainExp, Sexp, params, normalization_factor)
            curr_norm = norm(r)

            JtJ = J' * J
            Jtr = J' * r

            if curr_norm < toler
                # println("A good tolerance is reached")
                fit_successful = true
                break
            end

            # LM inner loop: adjust damping until improvement
            while true
                if τ > 1e20
                    fit_stalled = true
                    break
                end

                # (JᵀJ + τI) δ = Jᵀ r
                δ = (JtJ + τ * I) \ Jtr
                params_new = params .- δ

                r_new, _ = calc_residual_and_jacobian(dataType, strainExp, Sexp, params_new, normalization_factor)
                new_norm = norm(r_new)

                if new_norm < curr_norm
                    τ /= damping_decrease
                    params .= params_new
                    break
                else
                    τ *= damping_increase
                end
            end

            if fit_stalled
                perturbation_factor *= perturbation_rate
                break
            end
        end

        if fit_successful
            break
        end
    end

    return params
end

function calc_mat_constants_hyper_shear(modelType::AbstractString,
    strainExp::AbstractVector,
    Sexp::AbstractVector)

    # detect zero stress and filter
    keep = Sexp .!= 0
    Sexp = Sexp[keep]
    strainExp = strainExp[keep]

    W = 1.0 ./ Sexp
    λ = 1 .+ strainExp

    x = nothing
    if modelType == "neo-hookean"
        x = zeros(length(strainExp), 1)
        for i in eachindex(strainExp)
            xx = λ[i]
            x[i, 1] = 2 * (xx - 1 / xx^3)
        end

    elseif modelType == "mooney-rivlin"
        x = zeros(length(strainExp), 2)
        for i in eachindex(strainExp)
            xx = λ[i]
            x[i, 1] = 2 * (xx - 1 / xx^3)
            x[i, 2] = 2 * (xx - 1 / xx^3)   # matches your MATLAB code
        end

    elseif modelType == "yeoh"
        x = zeros(length(strainExp), 3)
        for i in eachindex(strainExp)
            xx = λ[i]
            Ibar1 = xx^2 + 1 + (1.0 / xx)^2
            base = 2 * (xx - 1 / xx^3)
            x[i, 1] = base
            x[i, 2] = 2 * base * (Ibar1 - 3)
            x[i, 3] = 3 * base * (Ibar1 - 3)^2
        end

    else
        error("Unknown modelType: $modelType (expected 'neo-hookean', 'mooney-rivlin', or 'yeoh').")
    end

    xnew = Diagonal(W) * x
    SeNew = Diagonal(W) * Sexp

    mat_cons = (xnew' * xnew) \ (xnew' * SeNew)
    return mat_cons
end


function calc_mat_constants_hyper_uniaxial(modelType::AbstractString,
    strainExp::AbstractVector,
    Sexp::AbstractVector)

    # detect zero stress and filter
    keep = Sexp .!= 0
    Sexp = Sexp[keep]
    strainExp = strainExp[keep]

    W = 1.0 ./ Sexp
    λ = 1 .+ strainExp

    x = nothing
    if modelType == "neo-hookean"
        x = zeros(length(strainExp), 1)
        for i in eachindex(strainExp)
            xx = λ[i]
            x[i, 1] = 2 * (xx - 1 / xx^2)
        end

    elseif modelType == "mooney-rivlin"
        x = zeros(length(strainExp), 2)
        for i in eachindex(strainExp)
            xx = λ[i]
            x[i, 1] = 2 * (xx - 1 / xx^2)
            x[i, 2] = 2 * (1 - 1 / xx^3)
        end

    elseif modelType == "yeoh"
        x = zeros(length(strainExp), 3)
        for i in eachindex(strainExp)
            xx = λ[i]
            Ibar1 = xx^2 + 2 * (1.0 / sqrt(xx))^2  # equals xx^2 + 2/xx
            base = 2 * (xx - 1 / xx^2)
            x[i, 1] = base
            x[i, 2] = 2 * base * (Ibar1 - 3)
            x[i, 3] = 3 * base * (Ibar1 - 3)^2
        end

    else
        error("Unknown modelType: $modelType (expected 'neo-hookean', 'mooney-rivlin', or 'yeoh').")
    end

    xnew = Diagonal(W) * x
    SeNew = Diagonal(W) * Sexp

    mat_cons = (xnew' * xnew) \ (xnew' * SeNew)
    return mat_cons
end

function drucker_stability_biaxial(modelType::AbstractString,
    mat_cons,
    λ::Real)

    # Helper to read coefficients whether mat_cons is a Vector or a Matrix
    getc(i) = ndims(mat_cons) == 1 ? mat_cons[i] : mat_cons[i, 1]

    λ1 = λ
    λ2 = λ1
    λ3 = λ1^(-2)

    if modelType == "neo-hookean"
        a = getc(1)
        D11 = 4 * (λ1^2 + λ3^2) * a
        D22 = 4 * (λ2^2 + λ3^2) * a
        D12 = 4 * λ3^2 * a
        term1 = D11 + D22
        term2 = D11 * D22 - D12^2
        dtdλ = 4 * a * (1 + 5 * λ1^(-6))
        return (term1 <= 0 || term2 <= 0 || dtdλ <= 0) ? 0 : 1

    elseif modelType == "mooney-rivlin"
        a = getc(1)
        b = getc(2)
        D11 = 4 * (λ1^2 + λ3^2) * (a + b * λ2^2)
        D22 = 4 * (λ2^2 + λ3^2) * (a + b * λ1^2)
        D12 = 4 * λ3^2 * a + 4 * b * λ3^(-2)
        term1 = D11 + D22
        term2 = D11 * D22 - D12^2
        return (term1 <= 0 || term2 <= 0) ? 0 : 1

    elseif modelType == "yeoh"
        a = getc(1)
        b = getc(2)
        c = getc(3)
        J = λ1 * λ2 * λ3
        I1 = J^(-2 / 3) * (λ1 + λ2 + λ3)  # mirrors your MATLAB (note: not squared terms)
        dWdI1 = a + 2 * b * (I1 - 3) + 3 * c * (I1 - 3)^2
        d2WdI1 = 2 * b + 6 * c * (I1 - 3)

        D11 = 4 * (λ1^2 + λ3^2) * dWdI1 + 4 * (λ1^2 - λ3^2)^2 * d2WdI1
        D22 = 4 * (λ2^2 + λ3^2) * dWdI1 + 4 * (λ2^2 - λ3^2)^2 * d2WdI1
        D12 = 4 * λ3^2 * dWdI1 + 4 * (λ1^2 - λ3^2) * (λ2^2 - λ3^2) * d2WdI1
        term1 = D11 + D22
        term2 = D11 * D22 - D12^2
        return (term1 <= 0 || term2 <= 0) ? 0 : 1

    else
        error("Unknown modelType: $modelType (expected 'neo-hookean', 'mooney-rivlin', or 'yeoh').")
    end
end

function drucker_stability_hyper(modelType::AbstractString, mat_cons,)
    # expects the following helpers to exist:
    #   drucker_stability_uniaxial(modelType, mat_cons, λ)::Int
    #   drucker_stability_biaxial(modelType, mat_cons, λ)::Int
    #   drucker_stability_shear(modelType,   mat_cons, λ)::Int

    flag_solver = 1

    lambdaArray_pos = 1.01:0.001:8.1
    lambdaArray_neg = 0.99:-0.01:0.5

    noNeedToCalc_uni = false
    noNeedToCalc_biax = false
    noNeedToCalc_shear = false

    # placeholders (to avoid undefined if never tripped)
    flag_uni = 1
    flag_biax = 1
    flag_sh = 1
    lambda_fail_uni = NaN
    lambda_fail_biax = NaN
    lambda_fail_sh = NaN

    # -------- positive lambda (tension) --------
    type = " tension"
    for λ in lambdaArray_pos
        if !noNeedToCalc_uni
            flag_uni = drucker_stability_uniaxial(modelType, mat_cons, λ)
            if flag_uni == 0
                noNeedToCalc_uni = true
                lambda_fail_uni = λ
            end
        end
        if !noNeedToCalc_biax
            flag_biax = drucker_stability_biaxial(modelType, mat_cons, λ)
            if flag_biax == 0
                noNeedToCalc_biax = true
                lambda_fail_biax = λ
            end
        end
        if !noNeedToCalc_shear
            flag_sh = drucker_stability_shear(modelType, mat_cons, λ)
            if flag_sh == 0
                noNeedToCalc_shear = true
                lambda_fail_sh = λ
            end
        end
    end

    if flag_uni == 0
        analysisType = "uniaxial"
        println("Stability margin in ", analysisType, type, " is ", -1 + lambda_fail_uni)
        flag_solver = 0
    elseif flag_uni == 1
        analysisType = "uniaxial"
        println("The material is stable for all strains in ", analysisType, type)
    end

    if flag_biax == 0
        analysisType = "biaxial"
        println("Stability margin in ", analysisType, type, " is ", -1 + lambda_fail_biax)
        flag_solver = 0
    elseif flag_biax == 1
        analysisType = "biaxial"
        println("The material is stable for all strains in ", analysisType, type)
    end

    if flag_sh == 0
        analysisType = "shear"
        println("Stability margin in ", analysisType, type, " is ", -1 + lambda_fail_sh)
        flag_solver = 0
    elseif flag_sh == 1
        analysisType = "shear"
        println("The material is stable for all strains in ", analysisType, type)
    end

    # -------- negative lambda (compression) --------
    noNeedToCalc_uni = false
    noNeedToCalc_biax = false
    noNeedToCalc_shear = false
    type = " compression"

    for λ in lambdaArray_neg
        if !noNeedToCalc_uni
            flag_uni = drucker_stability_uniaxial(modelType, mat_cons, λ)
            if flag_uni == 0
                noNeedToCalc_uni = true
                lambda_fail_uni = λ
            end
        end
        if !noNeedToCalc_biax
            flag_biax = drucker_stability_biaxial(modelType, mat_cons, λ)
            if flag_biax == 0
                noNeedToCalc_biax = true
                lambda_fail_biax = λ
            end
        end
        if !noNeedToCalc_shear
            flag_sh = drucker_stability_shear(modelType, mat_cons, λ)
            if flag_sh == 0
                noNeedToCalc_shear = true
                lambda_fail_sh = λ
            end
        end
    end

    if flag_uni == 0
        analysisType = "uniaxial"
        println("Stability margin in ", analysisType, type, " is ", -1 + lambda_fail_uni)
        flag_solver = 0
    elseif flag_uni == 1
        analysisType = "uniaxial"
        println("The material is stable for all strains in ", analysisType, type)
    end

    if flag_biax == 0
        analysisType = "biaxial"
        println("Stability margin in ", analysisType, type, " is ", -1 + lambda_fail_biax)
        flag_solver = 0
    elseif flag_biax == 1
        analysisType = "biaxial"
        println("The material is stable for all strains in ", analysisType, type)
    end

    if flag_sh == 0
        analysisType = "shear"
        println("Stability margin in ", analysisType, type, " is ", -1 + lambda_fail_sh)
        flag_solver = 0
    elseif flag_sh == 1
        analysisType = "shear"
        println("The material is stable for all strains in ", analysisType, type)
    end

    return flag_solver
end

function drucker_stability_ogden(data_type::AbstractString, mat_cons)
    # c by loading mode
    c = data_type == "uniaxial" ? -0.5 :
        data_type == "biaxial" ? -2.0 :
        data_type == "shear" ? -1.0 :
        error("Unknown data_type: $data_type (expected 'uniaxial', 'biaxial', or 'shear').")

    # helper to read params whether mat_cons is Vector or a column Matrix
    getc(i) = ndims(mat_cons) == 1 ? mat_cons[i] : mat_cons[i, 1]

    lambdaArray_pos = 1.01:0.001:8.1
    lambdaArray_neg = 0.99:-0.01:0.5

    # ---- positive lambda (tension) ----
    flag_pos = 1
    lam_pos = last(lambdaArray_pos)
    for λ in lambdaArray_pos
        λ1 = λ
        λ2 = λ1^c

        D11 = 0.0
        D22 = 0.0
        D12 = 0.0
        @inbounds for i in 1:3
            μ = getc(2i - 1)
            α = getc(2i)
            fac = 2 * μ * (λ1^(-α)) * (λ2^(-c * α))
            D11 += fac * (λ1^(2α) * λ2^(c * α) + 1)
            D22 += fac * (λ1^(α) * λ2^(2c * α) + 1)
            D12 += fac
        end

        term1 = D11 + D22
        term2 = D11 * D22 - D12^2
        if term1 <= 0 || term2 <= 0
            flag_pos = 0
            lam_pos = λ
            break
        end
    end

    # ---- negative lambda (compression) ----
    flag_neg = 1
    lam_neg = last(lambdaArray_neg)
    for λ in lambdaArray_neg
        λ1 = λ
        λ2 = λ1^c

        D11 = 0.0
        D22 = 0.0
        D12 = 0.0
        @inbounds for i in 1:3
            μ = getc(2i - 1)
            α = getc(2i)
            fac = 2 * μ * (λ1^(-α)) * (λ2^(-c * α))
            D11 += fac * (λ1^(2α) * λ2^(c * α) + 1)
            D22 += fac * (λ1^(α) * λ2^(2c * α) + 1)
            D12 += fac
        end

        term1 = D11 + D22
        term2 = D11 * D22 - D12^2
        if term1 <= 0 || term2 <= 0
            flag_neg = 0
            lam_neg = λ
            break
        end
    end

    return flag_pos, lam_pos, flag_neg, lam_neg
end

function drucker_stability_shear(modelType::AbstractString, mat_cons, λ::Real)
    # Read coefficients whether mat_cons is a Vector or 1-column Matrix
    getc(i) = ndims(mat_cons) == 1 ? mat_cons[i] : mat_cons[i, 1]

    λ1 = λ
    λ2 = 1.0
    λ3 = λ1^(-1)

    if modelType == "neo-hookean"
        a = getc(1)
        D11 = 4 * (λ1^2 + λ3^2) * a
        D22 = 4 * (λ2^2 + λ3^2) * a
        D12 = 4 * λ3^2 * a
        term1 = D11 + D22
        term2 = D11 * D22 - D12^2
        return (term1 <= 0 || term2 <= 0 || a <= 0) ? 0 : 1

    elseif modelType == "mooney-rivlin"
        a = getc(1)
        b = getc(2)
        D11 = 4 * (λ1^2 + λ3^2) * (a + b * λ2^2)
        D22 = 4 * (λ2^2 + λ3^2) * (a + b * λ1^2)
        D12 = 4 * λ3^2 * a + 4 * b * λ3^(-2)
        term1 = D11 + D22
        term2 = D11 * D22 - D12^2
        return (term1 <= 0 || term2 <= 0) ? 0 : 1

    elseif modelType == "yeoh"
        a = getc(1)
        b = getc(2)
        c = getc(3)
        J = λ1 * λ2 * λ3
        I1 = J^(-2 / 3) * (λ1 + λ2 + λ3)   # mirrors your MATLAB expression
        dWdI1 = a + 2 * b * (I1 - 3) + 3 * c * (I1 - 3)^2
        d2WdI1 = 2 * b + 6 * c * (I1 - 3)

        D11 = 4 * (λ1^2 + λ3^2) * dWdI1 + 4 * (λ1^2 - λ3^2)^2 * d2WdI1
        D22 = 4 * (λ2^2 + λ3^2) * dWdI1 + 4 * (λ2^2 - λ3^2)^2 * d2WdI1
        D12 = 4 * λ3^2 * dWdI1 + 4 * (λ1^2 - λ3^2) * (λ2^2 - λ3^2) * d2WdI1
        term1 = D11 + D22
        term2 = D11 * D22 - D12^2
        return (term1 <= 0 || term2 <= 0) ? 0 : 1

    else
        error("Unknown modelType: $modelType (expected 'neo-hookean', 'mooney-rivlin', or 'yeoh').")
    end
end

function drucker_stability_uniaxial(modelType::AbstractString, mat_cons, λ::Real)
    # read coefficients whether mat_cons is a Vector or a 1-column Matrix
    getc(i) = ndims(mat_cons) == 1 ? mat_cons[i] : mat_cons[i, 1]

    λ1 = λ
    λ2 = λ1^(-0.5)
    λ3 = λ1^(-0.5)

    if modelType == "neo-hookean"
        a = getc(1)
        D11 = 4 * (λ1^2 + λ3^2) * a
        D22 = 4 * (λ2^2 + λ3^2) * a
        D12 = 4 * λ3^2 * a
        term1 = D11 + D22
        term2 = D11 * D22 - D12^2
        dtdλ = 2 * a * (1 + 2 * λ1^(-3))
        return (term1 <= 0 || term2 <= 0 || dtdλ <= 0) ? 0 : 1

    elseif modelType == "mooney-rivlin"
        a = getc(1)
        b = getc(2)
        D11 = 4 * (λ1^2 + λ3^2) * (a + b * λ2^2)
        D22 = 4 * (λ2^2 + λ3^2) * (a + b * λ1^2)
        D12 = 4 * λ3^2 * a + 4 * b * λ3^(-2)
        term1 = D11 + D22
        term2 = D11 * D22 - D12^2
        return (term1 <= 0 || term2 <= 0) ? 0 : 1

    elseif modelType == "yeoh"
        a = getc(1)
        b = getc(2)
        c = getc(3)
        J = λ1 * λ2 * λ3
        I1 = J^(-2 / 3) * (λ1 + λ2 + λ3)  # mirrors your MATLAB expression
        dWdI1 = a + 2 * b * (I1 - 3) + 3 * c * (I1 - 3)^2
        d2WdI1 = 2 * b + 6 * c * (I1 - 3)

        D11 = 4 * (λ1^2 + λ3^2) * dWdI1 + 4 * (λ1^2 - λ3^2)^2 * d2WdI1
        D22 = 4 * (λ2^2 + λ3^2) * dWdI1 + 4 * (λ2^2 - λ3^2)^2 * d2WdI1
        D12 = 4 * λ3^2 * dWdI1 + 4 * (λ1^2 - λ3^2) * (λ2^2 - λ3^2) * d2WdI1
        term1 = D11 + D22
        term2 = D11 * D22 - D12^2
        return (term1 <= 0 || term2 <= 0) ? 0 : 1

    else
        error("Unknown modelType: $modelType (expected 'neo-hookean', 'mooney-rivlin', or 'yeoh').")
    end
end

function solver_constants_hyper(data_type::AbstractString,
    modelType::AbstractString,
    strainExp::AbstractVector,
    Sexp::AbstractVector)

    # Guard rails
    if data_type == "shear" && modelType == "mooney-rivlin"
        error("****Cannot evaluate Mooney-Rivlin model for pure shear. Use Yeoh model or etc.***")
    elseif data_type == "shear" && modelType == "ogden"
        error("****Cannot evaluate ogden model for pure shear. Use Yeoh model or etc.***")
    end

    # Fit material constants
    mat_cons =
        if modelType == "ogden"
            calc_mat_constants_hyper_ogden(data_type, strainExp, Sexp)
        elseif data_type == "uniaxial"
            calc_mat_constants_hyper_uniaxial(modelType, strainExp, Sexp)
        elseif data_type == "biaxial"
            calc_mat_constants_hyper_biaxial(modelType, strainExp, Sexp)
        elseif data_type == "shear"
            calc_mat_constants_hyper_shear(modelType, strainExp, Sexp)
        else
            error("Unknown data_type: $data_type")
        end

    # ---- ensure this is defined before any prints/returns ----
    mat_cons_solver = ""

    if modelType == "ogden"
        flag_pos, lam_pos, flag_neg, lam_neg = drucker_stability_ogden(data_type, mat_cons)
        if flag_pos == 0
            println("Instability margin for positive stretches is: ", lam_pos)
            mat_cons_solver = "Some instabilities exist"
        elseif flag_neg == 0
            println("Instability margin for negative stretches is: ", lam_neg)
            mat_cons_solver = "Some instabilities exist"
        else
            mat_cons_solver = "Successful: stable material constants"
        end
    else
        flag = drucker_stability_hyper(modelType, mat_cons)
        mat_cons_solver = (flag == 1) ? "Successful" : "Some instabilities exist"
    end

    # println("Material constants:")
    # println(permutedims(mat_cons))  # prints as a row, similar to MATLAB mat_cons'
    println(mat_cons_solver)

    return mat_cons
end
