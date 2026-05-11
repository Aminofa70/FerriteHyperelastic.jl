import{_ as s,o as n,c as e,aB as p}from"./chunks/framework.EAKmaEEb.js";const h=JSON.parse('{"title":"This tutorial is for curve fitting in hyperelastic models.","description":"","frontmatter":{},"headers":[],"relativePath":"tutorials/tutorial_curve_fitting.md","filePath":"tutorials/tutorial_curve_fitting.md","lastUpdated":null}'),i={name:"tutorials/tutorial_curve_fitting.md"};function t(l,a,o,c,r,u){return n(),e("div",null,[...a[0]||(a[0]=[p(`<h1 id="This-tutorial-is-for-curve-fitting-in-hyperelastic-models." tabindex="-1">This tutorial is for curve fitting in hyperelastic models. <a class="header-anchor" href="#This-tutorial-is-for-curve-fitting-in-hyperelastic-models." aria-label="Permalink to &quot;This tutorial is for curve fitting in hyperelastic models. {#This-tutorial-is-for-curve-fitting-in-hyperelastic-models.}&quot;">​</a></h1><div class="tip custom-block"><p class="custom-block-title">Hint</p><p>The package supports</p><ul><li><p>neo-Hookean model</p></li><li><p>Mooney-Rivlin model</p></li><li><p>Ogden model</p></li><li><p>Yeoh models model</p></li></ul></div><div class="tip custom-block"><p class="custom-block-title">Note</p><p>The package also supports <em>the Drucker stability criterion</em>.</p></div><p>In order to perform curve fitting, experimental data are required. Some experimental data are available in <a href="https://github.com/Aminofa70/HyperData.jl" target="_blank" rel="noreferrer">HyperData.jl</a> repository.</p><p>The experimental data can be <em>Strain</em> or <em>Stretch</em> for <code>x</code> axis and the nominal stress for <code>y</code>axis.</p><div class="tip custom-block"><p class="custom-block-title">Steps TODO Curve Fitting</p><ul><li><p>Load Package</p></li><li><p>Read Experimental Data</p></li><li><p>Curve Fitting Inputs</p></li><li><p>Call curve fitting function</p></li><li><p>Check the accuracy of the obtained material constant(s).</p></li></ul></div><h3 id="Load-Package" tabindex="-1">Load Package <a class="header-anchor" href="#Load-Package" aria-label="Permalink to &quot;Load Package {#Load-Package}&quot;">​</a></h3><p>Load the required packages for curve fitting as</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>using Revise          # to automatically reload changed code</span></span>
<span class="line"><span>using FerriteHyperelastic</span></span>
<span class="line"><span>using GLMakie</span></span>
<span class="line"><span>using HyperData</span></span></code></pre></div><h3 id="Read-Experimental-Data" tabindex="-1">Read Experimental Data <a class="header-anchor" href="#Read-Experimental-Data" aria-label="Permalink to &quot;Read Experimental Data {#Read-Experimental-Data}&quot;">​</a></h3><p>Read the experimental data as</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>λ_exp, P_exp = read_data!(Treloar_1944, uniaxial)</span></span></code></pre></div><h3 id="Curve-Fitting-Inputs" tabindex="-1">Curve Fitting Inputs <a class="header-anchor" href="#Curve-Fitting-Inputs" aria-label="Permalink to &quot;Curve Fitting Inputs {#Curve-Fitting-Inputs}&quot;">​</a></h3><p>Set inputs of curve fitting function</p><div class="warning custom-block"><p class="custom-block-title">Warning</p><p>The input MUST be Strain and the Nominal Stress, for this reasion if the data is different user must convert them to starin and nominal stress</p></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span># Input parameters</span></span>
<span class="line"><span>inputDataType = &quot;lambda&quot;</span></span>
<span class="line"><span>data_type = &quot;uniaxial&quot;</span></span>
<span class="line"><span>modelType = &quot;neo-hookean&quot;</span></span>
<span class="line"><span></span></span>
<span class="line"><span>Sexp = P_exp</span></span>
<span class="line"><span>strainExp = λ_exp .- 1</span></span></code></pre></div><div class="tip custom-block"><p class="custom-block-title">Note</p><p>Please use lowercase letter for all model types (this example is &quot;neo-hookean&quot;)</p></div><h3 id="Call-curve-fitting-function" tabindex="-1">Call curve fitting function <a class="header-anchor" href="#Call-curve-fitting-function" aria-label="Permalink to &quot;Call curve fitting function {#Call-curve-fitting-function}&quot;">​</a></h3><p>The curve-fitting function is now called. For neo-Hookean, the material constant is C1.</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span># Get material constants by fitting (assumed function)</span></span>
<span class="line"><span>mat_cons_solver = solver_constants_hyper(data_type, modelType, strainExp, Sexp)</span></span>
<span class="line"><span>C1 = mat_cons_solver[1]</span></span></code></pre></div><h3 id="Check-the-accuracy-of-the-obtained-material-constants" tabindex="-1">Check the accuracy of the obtained material constant(s) <a class="header-anchor" href="#Check-the-accuracy-of-the-obtained-material-constants" aria-label="Permalink to &quot;Check the accuracy of the obtained material constant(s) {#Check-the-accuracy-of-the-obtained-material-constants}&quot;">​</a></h3><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span># Generate strain values for smooth curve</span></span>
<span class="line"><span>ϵ = collect(LinRange(strainExp[1], strainExp[end], 100))</span></span>
<span class="line"><span></span></span>
<span class="line"><span># Calculate stretches for incompressible uniaxial tension</span></span>
<span class="line"><span>λ1 = @. ϵ + 1</span></span>
<span class="line"><span>λ = λ1</span></span>
<span class="line"><span>λ2 = @. 1 / sqrt(λ)</span></span>
<span class="line"><span>λ3 = λ2   </span></span>
<span class="line"><span></span></span>
<span class="line"><span># neo-Hookean stress (uniaxial nominal stress)</span></span>
<span class="line"><span>P_model = @. 2 * C1 * (λ - 1 / λ^2)</span></span>
<span class="line"><span></span></span>
<span class="line"><span># Plot results</span></span>
<span class="line"><span>GLMakie.closeall()</span></span>
<span class="line"><span></span></span>
<span class="line"><span>fig = Figure(size=(800, 600), fontsize=26)</span></span>
<span class="line"><span>ax = Axis(fig[1, 1], xlabel= L&quot;\\mathscr{ε}&quot;, ylabel=L&quot;P&quot;,  xgridvisible=false, ygridvisible=false)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>lines!(ax, ϵ, P_model, color = :black, label=&quot;Fit, neo-Hookean&quot;)</span></span>
<span class="line"><span>scatter!(ax, strainExp, P_exp, marker=:circle, color=:red, label=&quot;Uniaxial experiment&quot;)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>axislegend(ax, position=:lt, backgroundcolor=(:white, 0.7), framecolor=:gray)</span></span>
<span class="line"><span>display(fig)</span></span></code></pre></div><h2 id="Plain-program" tabindex="-1">Plain program <a class="header-anchor" href="#Plain-program" aria-label="Permalink to &quot;Plain program {#Plain-program}&quot;">​</a></h2><p>The version of the code without comments.</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>using Revise</span></span>
<span class="line"><span>using FerriteHyperelastic</span></span>
<span class="line"><span>using GLMakie</span></span>
<span class="line"><span>using HyperData</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>λ_exp, P_exp = read_data!(Treloar_1944, uniaxial)</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>inputDataType = &quot;lambda&quot;</span></span>
<span class="line"><span>data_type = &quot;uniaxial&quot;</span></span>
<span class="line"><span>modelType = &quot;neo-hookean&quot;</span></span>
<span class="line"><span></span></span>
<span class="line"><span>Sexp = P_exp</span></span>
<span class="line"><span>strainExp = λ_exp .- 1</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>mat_cons_solver = solver_constants_hyper(data_type, modelType, strainExp, Sexp)</span></span>
<span class="line"><span>C1 = mat_cons_solver[1]</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>ϵ = collect(LinRange(strainExp[1], strainExp[end], 100))</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>λ1 = @. ϵ + 1</span></span>
<span class="line"><span>λ = λ1</span></span>
<span class="line"><span>λ2 = @. 1 / sqrt(λ)</span></span>
<span class="line"><span>λ3 = λ2   </span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>P_model = @. 2 * C1 * (λ - 1 / λ^2)</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>GLMakie.closeall()</span></span>
<span class="line"><span></span></span>
<span class="line"><span>fig = Figure(size=(800, 600), fontsize=26)</span></span>
<span class="line"><span>ax = Axis(fig[1, 1], xlabel= L&quot;\\mathscr{ε}&quot;, ylabel=L&quot;P&quot;,  xgridvisible=false, ygridvisible=false)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>lines!(ax, ϵ, P_model, color = :black, label=&quot;Fit, neo-Hookean&quot;)</span></span>
<span class="line"><span>scatter!(ax, strainExp, P_exp, marker=:circle, color=:red, label=&quot;Uniaxial experiment&quot;)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>axislegend(ax, position=:lt, backgroundcolor=(:white, 0.7), framecolor=:gray)</span></span>
<span class="line"><span>display(fig)</span></span></code></pre></div>`,25)])])}const m=s(i,[["render",t]]);export{h as __pageData,m as default};
