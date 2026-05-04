import{_ as n,o as e,c as a,aB as t}from"./chunks/framework.CAlbP-O0.js";const h=JSON.parse('{"title":"This is the instruction to use the package","description":"","frontmatter":{},"headers":[],"relativePath":"instruction.md","filePath":"instruction.md","lastUpdated":null}'),p={name:"instruction.md"};function i(l,s,r,c,o,d){return e(),a("div",null,[...s[0]||(s[0]=[t(`<h1 id="This-is-the-instruction-to-use-the-package" tabindex="-1">This is the instruction to use the package <a class="header-anchor" href="#This-is-the-instruction-to-use-the-package" aria-label="Permalink to &quot;This is the instruction to use the package {#This-is-the-instruction-to-use-the-package}&quot;">​</a></h1><h2 id="Guid-to-use-Finite-Elementfem" tabindex="-1">Guid to use Finite Element(fem) <a class="header-anchor" href="#Guid-to-use-Finite-Elementfem" aria-label="Permalink to &quot;Guid to use Finite Element(fem) {#Guid-to-use-Finite-Elementfem}&quot;">​</a></h2><p>First we need to activate a dynamic structure to have all finite element parameters in it. We do it using</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>input = InputStruct()</span></span>
<span class="line"><span></span></span>
<span class="line"><span>\`\`\`\`</span></span>
<span class="line"><span>Here \`\`\`input\`\`\` is the a structure that will contains all fem parameters. </span></span>
<span class="line"><span></span></span>
<span class="line"><span>Now using [Ferrite.jl](https://ferrite-fem.github.io/Ferrite.jl/stable/) we define geometry, mesh, fem parameters for interpolation and numerical integration and define degree-of-freedom. More details for this parts can be found in [Ferrite.jl](https://ferrite-fem.github.io/Ferrite.jl/stable/). </span></span>
<span class="line"><span></span></span>
<span class="line"><span>As an example, we do it for a two dimensional problem.</span></span>
<span class="line"><span>Our element here is quad4 and then we define mesh and geometry as</span></span></code></pre></div><p>&quot;&quot;&quot; Lx : Lenght in x-direction Ly : Lenght in y-direction Ferrite.Quadrilateral: Define the type of elements. nx : Number of element in x-dir ny : Number of element in y-dir &quot;&quot;&quot; function create_grid(Lx, Ly, nx, ny) corners = [ Ferrite.Vec{2}((0.0, 0.0)), Ferrite.Vec{2}((Lx, 0.0)), Ferrite.Vec{2}((Lx, Ly)), Ferrite.Vec{2}((0.0, Ly)) ] grid = Ferrite.generate_grid(Ferrite.Quadrilateral, (nx, ny), corners) addnodeset!(grid, &quot;support_1&quot;, x -&gt; x[1] ≈ 0.0) addfacetset!(grid, &quot;pressure&quot;, x -&gt; x[1] ≈ Lx) return grid end</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>It is noted that</span></span></code></pre></div><p>addnodeset!(grid, &quot;support_1&quot;, x -&gt; x[1] ≈ 0.0) addfacetset!(grid, &quot;pressure&quot;, x -&gt; x[1] ≈ Lx)</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>are used to define our boundary conditions.</span></span>
<span class="line"><span></span></span>
<span class="line"><span>If you wan to chose other element types like tri3 and so on please refer to [generate_grid](https://github.com/Ferrite-FEM/Ferrite.jl/blob/f1d1d0deef7bdaf019bd63ce9e8d959b6ebc8c4d/src/Grid/grid_generators.jl#L1-L7) .</span></span>
<span class="line"><span></span></span>
<span class="line"><span>Next step is defining the fem values (interpolation and numerical integration)</span></span></code></pre></div><p>function create_values() dim, order = 2, 1 ip = Ferrite.Lagrange{Ferrite.RefQuadrilateral,order}()^dim qr = Ferrite.QuadratureRule{Ferrite.RefQuadrilateral}(2) qr_face = Ferrite.FacetQuadratureRule{Ferrite.RefQuadrilateral}(2) return Ferrite.CellValues(qr, ip), Ferrite.FacetValues(qr_face, ip) end</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Very important to know is that the type here should be the same as the element type. For example here we have</span></span></code></pre></div><p>Ferrite.Quadrilateral <code>and</code> Ferrite.RefQuadrilateral \`\`\` .</p><p>Then defining the dof (here is 2D, <code>ux</code>,<code>uy</code>)</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>function create_dofhandler(grid)</span></span>
<span class="line"><span>    dh = Ferrite.DofHandler(grid)</span></span>
<span class="line"><span>    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefQuadrilateral,1}()^2)</span></span>
<span class="line"><span>    Ferrite.close!(dh)</span></span>
<span class="line"><span>    return dh</span></span>
<span class="line"><span>end</span></span></code></pre></div><p>Now, we define the Dirichlet boundary condition</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>function create_bc(dh)</span></span>
<span class="line"><span>    ch = Ferrite.ConstraintHandler(dh)</span></span>
<span class="line"><span>    Ferrite.add!(ch, Ferrite.Dirichlet(:u, Ferrite.getnodeset(dh.grid, &quot;support_1&quot;), (x, t) -&gt; [0.0, 0.0], [1, 2]))</span></span>
<span class="line"><span>    Ferrite.close!(ch)</span></span>
<span class="line"><span>    return ch</span></span>
<span class="line"><span>end</span></span></code></pre></div><p>If we want to plot force-displacement curve, we need to define another boundary condition to get dof of the displacment for the plot. Displacement where the force is applied.</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>function create_bc_force(dh)</span></span>
<span class="line"><span>    dbc = Ferrite.ConstraintHandler(dh)</span></span>
<span class="line"><span>    Ferrite.add!(dbc, Ferrite.Dirichlet(:u, getfacetset(grid, &quot;pressure&quot;), (x, t) -&gt; 0*x))</span></span>
<span class="line"><span>    Ferrite.close!(dbc)</span></span>
<span class="line"><span>    return dbc</span></span>
<span class="line"><span>end</span></span></code></pre></div><p>To now we defined the required functions for fem. These parts are from the package <a href="https://ferrite-fem.github.io/Ferrite.jl/stable/" target="_blank" rel="noreferrer">Ferrite.jl</a>. More details can be seen in this package.</p><hr><p>The hyperelastic strain enenrgy function and stress are now defined. We use the neo-Hookean function, but other strain energy functions can also be employed</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>function Ψ(C, C10, D1)</span></span>
<span class="line"><span>    J = sqrt(det(C))</span></span>
<span class="line"><span>    I1 = tr(C)</span></span>
<span class="line"><span>    I1_bar = I1 * J^(-2 / 3)</span></span>
<span class="line"><span>    return C10 * (I1_bar - 3) + (1 / D1) * (J - 1)^2</span></span>
<span class="line"><span>end</span></span></code></pre></div><p>The second Piola kirchhoff stress and its differentiation are then define</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>function constitutive_driver(C, C10, D1)</span></span>
<span class="line"><span>    ∂²Ψ∂C², ∂Ψ∂C = Tensors.hessian(y -&gt; Ψ(y, C10, D1), C, :all)</span></span>
<span class="line"><span>    S = 2.0 * ∂Ψ∂C</span></span>
<span class="line"><span>    ∂S∂C = 2.0 * ∂²Ψ∂C²</span></span>
<span class="line"><span>    return S, ∂S∂C</span></span>
<span class="line"><span>end</span></span></code></pre></div><p>Now, we define a driver for the strain energy in order to inlcude it in the input structure</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>function make_constitutive_driver(C10, D1)</span></span>
<span class="line"><span>    return C -&gt; constitutive_driver(C, C10, D1)</span></span>
<span class="line"><span>end</span></span></code></pre></div><p>We should tell the solve our problem is 2D or 3D. Also we should define we have applied traction or displacement. Then we should script</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>input.model_type = :plane_strain   # or :plane_strain; :plane_stress; :threeD</span></span>
<span class="line"><span>input.load_type = :traction</span></span></code></pre></div><p>Now we should assign parameters for the materials ans also call the fem function to include them in the input structure.</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span></span></span>
<span class="line"><span>input.E , input.ν = 3.35, 0.45</span></span>
<span class="line"><span>E = input.E</span></span>
<span class="line"><span>ν = input.ν</span></span>
<span class="line"><span>C10 = E / (4 * (1 + ν))</span></span>
<span class="line"><span>D1 = 6.0 * (1.0 - 2.0 * ν) / E</span></span>
<span class="line"><span>input.material = make_constitutive_driver(C10, D1)</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>Lx, Ly = 3.17, 1.73  # Plate dimensions</span></span>
<span class="line"><span>nx, ny = 10, 10   # Number of elements along x and y</span></span>
<span class="line"><span>grid = create_grid(Lx, Ly, nx, ny)  # Generate the grid</span></span>
<span class="line"><span></span></span>
<span class="line"><span>input.grid = grid</span></span>
<span class="line"><span>input.dh = create_dofhandler(grid)</span></span>
<span class="line"><span>input.ch = create_bc(input.dh )</span></span>
<span class="line"><span># Create CellValues and FacetValues</span></span>
<span class="line"><span>input.cell_values, input.facet_values = create_values()</span></span></code></pre></div><p>Apply the traction (it can be multiple tractions)</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>input.ΓN = getfacetset(grid, &quot;pressure&quot;)</span></span>
<span class="line"><span>input.facetsets = [input.ΓN]</span></span>
<span class="line"><span>input.traction = [2.2, 0.0]</span></span>
<span class="line"><span>input.tractions = Dict(1 =&gt; input.traction)</span></span></code></pre></div><p>If the aim is also plot force-displacement, we shoul also find the dof of reaction force and displacement (if not leave them empty)</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span></span></span>
<span class="line"><span>dof_F_x = input.ch.prescribed_dofs[1:2:end]</span></span>
<span class="line"><span>input.dof_F = dof_F_x;</span></span>
<span class="line"><span></span></span>
<span class="line"><span>dbc= create_bc_force(input.dh)</span></span>
<span class="line"><span>dof_U_x = dbc.prescribed_dofs[1:2:end] </span></span>
<span class="line"><span>input.dof_U = dof_U_x</span></span>
<span class="line"><span># input.dof_F = []</span></span>
<span class="line"><span># input.dof_U = []</span></span></code></pre></div><p>Define the tolernce for solver and also parameters for time integration. Because the problem is nonlinear we need to solve it incrementally (time step here), then</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>input.tol = 1e-6</span></span>
<span class="line"><span></span></span>
<span class="line"><span>## default</span></span>
<span class="line"><span>#maxIterPerInc,totalTime,initInc,minInc,maxInc,totalInc = initialize_solver()</span></span>
<span class="line"><span></span></span>
<span class="line"><span># change like the following if you need</span></span>
<span class="line"><span>maxIterPerInc,totalTime,initInc,minInc,maxInc,totalInc = initialize_solver(500,1.0,1e-3,1e-15,0.8,1000)</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>input.maxIterPerInc = maxIterPerInc</span></span>
<span class="line"><span>input.totalTime = totalTime</span></span>
<span class="line"><span>input.initInc = initInc</span></span>
<span class="line"><span>input.minInc = minInc</span></span>
<span class="line"><span>input.maxInc = maxInc</span></span>
<span class="line"><span>input.totalInc = totalInc</span></span></code></pre></div><p>The solver saves the displacement of the last time step in vtu file for analysis of the results. So solver needs a name and dir to save it</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>input.filename = &quot;2D_Hyper&quot;</span></span>
<span class="line"><span>input.output_dir= &quot;/Users/aminalibakhshi/Desktop/vtu_geo/&quot;</span></span></code></pre></div><p>Howover, the solver saves the results for each time steps and returns all corresponding displacement vector and use can plot it using Plot.jl or GLMakie.jl</p><p>Now, calling the solve</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>sol  = run_fem(input)</span></span></code></pre></div><p>Solver is a structure that returns <code>U = sol.U_steps[end]</code> the displacement for each step and <code>sol.F_effect</code> the sum of reaction force and <code>sol.U_effect</code> average displacement in applied force (or displacement).</p><p>The force-displacement can be plotted using GLMakie.jl as</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>GLMakie.closeall()</span></span>
<span class="line"><span>fig = Figure(size=(800, 600))</span></span>
<span class="line"><span>ax = Axis(fig[1, 1], xlabel=&quot;Displacement&quot;, ylabel=&quot;Force&quot;, title=&quot;force-displacement&quot;, xgridvisible = false, ygridvisible = false)</span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span></span></span>
<span class="line"><span>lines!((sol.U_effect), abs.(sol.F_effect), color = :black)</span></span>
<span class="line"><span>scatter!((sol.U_effect), abs.(sol.F_effect), marker = :circle , color = :red)</span></span>
<span class="line"><span>display(fig)</span></span></code></pre></div>`,43)])])}const g=n(p,[["render",i]]);export{h as __pageData,g as default};
