### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ c7ab89d0-638e-4121-842e-b14a7d58a4b2
using PlutoUI

# ╔═╡ 720128f9-b311-4883-9444-f7d564008f66
using Statistics

# ╔═╡ 84fb4b5c-8fbd-11ed-165d-b758c56cbf89
using Plots 

# ╔═╡ 47c2a121-b626-4f8b-9722-cefa925ba9a5
md"""
# Interactive post-processing of capsule simulations data
"""

# ╔═╡ 2012189e-7f9c-4a59-85cc-c6dd4b04d129
md"""
## Some definitions
"""

# ╔═╡ 9d776581-c80d-4760-9ac2-e01294a60c5e
md"""
### Coordinate structure

We first create a strucutre representing the coordinates of a Lagrangian node. It has three attributes: the node's x, y and z coordinates.
"""

# ╔═╡ 6e5461ea-fae3-4c21-abf2-d2aae9dc64ba
begin
	mutable struct Coord
		x::Float64
		y::Float64
		z::Float64
	end
	Coord() = Coord(0, 0, 0)
end

# ╔═╡ 4cb4b887-a081-43e2-94eb-bdac56cc405a
md"""
The function below parses a Coord from a string assuming the coordinates are separated by spaces.
"""

# ╔═╡ 56710fed-1eff-4263-a48d-6dacc371bb86
function Base.parse(::Type{Coord}, str)
	c = split(str, ' ', limit=3)
	return Coord(parse(Float64, c[1]), parse(Float64, c[2]), parse(Float64, c[3]))
end

# ╔═╡ c74548f1-ac38-4b5f-a990-df5dbf01b7a8
md"""
It will be handy to be able to add Coord together, as well as multiply them by a scalar. As such, we overload below the +, -, *, / operators.
"""

# ╔═╡ e5e71a51-fe18-42fe-9acb-49f1147845f6
Base.:+(a::Coord, b::Coord) = Coord(a.x + b.x, a.y + b.y, a.z + b.z)

# ╔═╡ 111b0b5b-f40b-41f1-a8b8-bc8fcac65c4e
Base.:-(a::Coord, b::Coord) = Coord(a.x - b.x, a.y - b.y, a.z - b.z)

# ╔═╡ 0366ffe5-9f64-46f0-9aad-7d19c779a2e3
Base.:/(a::Coord, b::Number) = Coord(a.x/b, a.y/b, a.z/b)

# ╔═╡ 5695e259-55fc-40c7-af3d-04447486f71b
Base.:*(a::Coord, b::Number) = Coord(a.x*b, a.y*b, a.z*b)

# ╔═╡ 853ca926-d374-46a7-bdad-9141e506cd79
md"""
The method below returns the centroid of an array of Coord, by simply averaging them.
"""

# ╔═╡ febb23c9-4924-496b-a8e5-c881f26c7d8c
function compute_centroid(nodes::Vector{Coord})
	centroid = Coord(0, 0, 0)
	N = length(nodes)
	for i in 1:N
		centroid += nodes[i]
	end
	return centroid/N
end

# ╔═╡ cfaa17a3-5eaf-479e-bae3-c55608be8d20
md"""
### Triangle structure

The membrane of our capsule is made of triangles which vertices are the Lagrangian nodes we introduced above. We create the `Triangle` structure below, which attributes are the indices of its vertices, as well as the triangle area.
"""

# ╔═╡ 1cfe7a07-411b-4a05-a489-a1739d9a0702
begin
	mutable struct Triangle
		i::Integer
		j::Integer
		k::Integer
		A::Float64
	end
	Triangle() = Triangle(0, 0, 0, 0)
end

# ╔═╡ 7b25dd66-e499-4a66-b6bd-ac7b9663a603
md"""
The function below parses a string as a trianble, assuming its four attributes are separated by spaces.
"""

# ╔═╡ 78619fd5-866c-4153-a40e-0ec9f9ed8f4a
function Base.parse(::Type{Triangle}, str)
	c = split(str, ' ', limit=4)
	return Triangle(parse(Int64, c[1]), parse(Int64, c[2]), parse(Int64, c[3]),
		parse(Float64, c[4]))
end

# ╔═╡ 32a4b2f1-4299-4549-afa9-3be80ce19ffd
md"""
### Capsule position at a given instant

We define below the structure `CapsFrame` that stores relevant information about the capsule at a given time instant. Its attributes are:

* `time`, the physical time corresponding to this snapshot of the capsule
* `nb_nodes`, the number of Lagrangian nodes on the capsule
* `nodes`, an array of Coord which store information about the capsule's nodes
* `centroid`, the centroid of the capsule
* `nb_tri`, the number of triangles tiling the capsule's surface
* `triangles`, an array of Triangles representing the triangulation of the capsule
* `area`, the total discretized area of the capsule
"""

# ╔═╡ 81e03b17-9741-4c66-8047-cce319f3900e
begin
	mutable struct CapsFrame
		time::Float64
		nb_nodes::Integer
		nodes::Vector{Coord}
		centroid::Coord
		nb_tri::Integer
		triangles::Vector{Triangle}
		area::Float64
	end
	CapsFrame() = CapsFrame(0, 0, [Coord()], Coord(), 0, [Triangle()], 0)
end

# ╔═╡ 17e704a2-39ec-412a-a552-abcd04fc1a26
md"""
The function below parses a `CapsFrame` from a string containing all the nodes informations. It is assumed that this string corresponds to a line in the CSV data file storing the nodes information. Indeed, the data file containing the positions of the $n$ membrane nodes for the time interval $[t_0 : t_{end}]$ is strucured as follows:

```
t0,x1 y1 z1,x1 y1 z1,...,xn yn zn
t1,x1 y1 z1,x1 y1 z1,...,xn yn zn
.
.
.
tend,x1 y1 z1,x1 y1 z1,...,xn yn zn
```
and the function below parses one line of it to a `CapsFrame` structure. The information about the capsule's triangles will be added in a separate method, because the nodes and triangles informations are stored in two separate CSV files.

"""

# ╔═╡ 91daca3f-575b-4f71-83a7-8abde506e1fe
function Base.parse(::Type{CapsFrame}, nodeStr::AbstractString)
	e = split(nodeStr, ',')
	time = parse(Float64, e[1])
	nb_nodes = length(e) - 1
	nodes = Vector{Coord}(undef, nb_nodes)
	for i in 1:nb_nodes
		nodes[i] = parse(Coord, e[i+1])
	end
	centroid = compute_centroid(nodes)
	return CapsFrame(time, nb_nodes, nodes, centroid, 0, [Triangle()], 0)
end

# ╔═╡ f1f1f88f-4a8d-44d3-a516-39e6fa40e44f
md"""
### Plot information structure
When plotting various quantitative capsule data, it will prove useful to automatize the rescaling process of the graph, such that the plots for various conditions all properly superimpose and are ready to be analyzed. 
To this end, we create a structure `PlotInfo` below
"""

# ╔═╡ 2b3d71a8-ac2a-4b64-b547-661646343995
begin
	mutable struct PlotInfo
		liter::Integer
		veq::Float64
		vmin::Float64
		ivmin::Integer
		t0::Float64
		times::Vector{Float64}
	end
	PlotInfo() = PlotInfo(0, 0, 0, 0, 0, [0])
end

# ╔═╡ 91be20f2-79dd-4641-a4f5-04306dc341fd
md"""
### Capsule structure

Below we define the `Capsule` structure that aggregates the information of a capsule at several time instants. Its attributes are:

* `nb_times`, the number of snapshots of the capsule
* `cf`, an array of `CapsFrame` structures
* `cvel`, an array containing the centroid velocity as a Coord
* `vel`, an array containing the _norm_ of the centroid velocity as a Float
* `perm`, a permutation array that transforms the array of nodes in `cts` such that the first $n_{pp}$ nodes in the plane $z=0$ are at the beginning of the array and ordered counter-clockwise. It is useful for plotting the capsule outline in the plane $z=0$
* `npp`, the number of nodes in the plane $z=0$
* `area`, an array of areas of the capsule, useful to plot the area against time
"""

# ╔═╡ cbadf507-94a9-47bc-93b6-8512235ea715
begin
	mutable struct Capsule
		nb_times::Integer
		cf::Vector{CapsFrame}
		cvel::Vector{Coord}
		vel::Vector{Float64}
		perm::Vector{Int}
		npp::Integer
		area::Vector{Float64}
		pi::PlotInfo
	end
	Capsule() = Capsule(0, [CapsFrame()], [Coord()], [0], [0], 0, [0], PlotInfo())
end

# ╔═╡ 1c81d5f7-6506-40d8-a9c9-19b11cb901dc
md"""
The method below computes the velocity of the capsule centroid by differentiating its position. It fills the attributes `cvel` and `vel` of the `Capsule` structure.
"""

# ╔═╡ 5d7d1025-fcad-4cd9-9bb6-a2a3f033772c
function compute_cvel!(caps::Capsule)
	nbt = caps.nb_times
	caps.cvel = [Coord() for _ in 1:nbt-1]
	caps.vel = [0. for _ in 1:nbt-1]
	for i in 1:nbt-1
		dt = caps.cf[i+1].time - caps.cf[i].time
		caps.cvel[i] = (caps.cf[i+1].centroid - caps.cf[i].centroid)/dt
		caps.vel[i] = sqrt(caps.cvel[i].x^2 + caps.cvel[i].y^2 + caps.cvel[i].z^2)
	end
end

# ╔═╡ 3ee686a6-81d5-419a-9d0c-41719c4fe2bd
md"""
The function below computes the permutation array that transforms the array of nodes in `cts` such that the first $n_{pp}$ nodes in the plane $z=0$ are at the beginning of the array and ordered counter-clockwise. It is useful for plotting the capsule outline in the plane $z=0$.

This method takes a `Capsule` structure as an argument and modified its `perm` and `npp` attributes.

The first loop creates a permutation array aimed at gathering all the indices of the nodes which are in the $z=0$ plane at the beginning of the array.

The second and third loops computes the angle of the nodes contained in the $z=0$ plane, and attributes a large dummy angle to the rest of the nodes.

Sorting the array of angles computed previously and getting the corresponding indices permutation gives the desired permutation array.
"""

# ╔═╡ 1fcaa2aa-6876-44c0-8db2-35d00211bff4
function compute_permutations!(caps::Capsule)
	N = caps.cf[1].nb_nodes
	caps.perm = [0 for i in 1:N]
	theta = Array{Float64}(undef, N)
	perm1 = collect(1:N)
	perm2 = Array{Int}(undef, N)
	npp = 0 # number of points in the xy plane
	for i in 1:N
		if abs(caps.cf[1].nodes[i].z) < .05
			npp += 1
			perm1[npp] = i
			perm1[i] = npp
		end
	end
	for i in 1:npp
		theta[i] = atan(caps.cf[1].nodes[perm1[i]].y - caps.cf[1].centroid.y, 
			caps.cf[1].nodes[perm1[i]].x - caps.cf[1].centroid.x)
	end
	for i in npp+1:N
		theta[i] = 3*pi
	end
	perm2 = sortperm(theta)
	caps.perm .= perm1[perm2]
	caps.npp = npp
end

# ╔═╡ e118297a-fcb8-409c-bb8c-e1992833fe85
md"""
The function below takes a Capsule structure and an CSV file containing the triangle informations, and fills the `nb_tri`, `triangles` and `area` attribute of each CapsTimeStap frame, as well as the `area` array in the Capsule strucutre.

It is assumed that the CSV data file storing the triangles information has the following structure:
```
t0,i1 j1 k1 a1,i2 j2 k2 a2,...,in jn kn an
t1,i1 j1 k1 a1,i2 j2 k2 a2,...,in jn kn an
.
.
.
tend,i1 j1 k1 a1,i2 j2 k2 a2,...,in jn kn an
```
with `i`, `j`, `k` the indices of the vertices of the `n` triangles and `a` their areas.

In some cases, one more snapshot of the triangles information was saved at an additional time, earlier than what is recorded in the nodes datafile. If this occurs, we skip the first line and throw a warning.
"""

# ╔═╡ 0c86edd5-fb43-4949-9ace-06f923daec2f
function read_triangles!(c::Capsule, filename::AbstractString)
	f = open(filename, "r")
	nl = 1
	c.area = [0 for i in 1:c.nb_times]
	offset = 0
	for line in readlines(f)
		e = split(line, ',')
		t0 = parse(Float64, e[1])	
		initial_line = 1
		nb_tri = length(e) - 1
		c.cf[nl - offset].triangles = [Triangle() for i in 1:nb_tri]
		c.area[nl - offset] = 0
		for i in initial_line:nb_tri
			c.cf[nl - offset].triangles[i] = parse(Triangle, e[i+1])
			c.area[nl - offset] += c.cf[nl - offset].triangles[i].A
		end
		if (nl == 1 && abs(t0 - c.cf[1].time) > 0.01)
			println("--- Warning: initial time in "*filename*" not matching the initial time in the nodes data file. Skipping line ", offset + 1,".")
			offset += 1
			println("Info from nodes data file: t0 = ", c.cf[1].time)
			println("Info from triangles data file: t0 = ", t0)
		end
		nl += 1
	end
end

# ╔═╡ d582cf21-202e-4fb0-af4e-cf73f44a4fa0
function generate_plot_info!(c::Capsule)
	c.pi.liter = c.nb_times
	c.pi.veq = c.vel[c.pi.liter-1]
	c.pi.vmin, c.pi.ivmin = findmin(c.vel[2:c.pi.liter-1])
	c.pi.t0 = c.cf[2+c.pi.ivmin].time
	c.pi.times = [c.cf[i].time for i in 1:c.nb_times]
end

# ╔═╡ 74454cc5-8be1-4e6d-b2d7-e160b785e91b
md"""
The function below fills the a `Capsule` structure from a data file containing nodes information, and optionnally from a data file containing triangles information.
"""

# ╔═╡ f86282ba-646f-4d89-92c5-c6fbd4af4748
function read_capsule!(c::Capsule, filename1::AbstractString, 
	filename2::AbstractString="")
	c.nb_times = 0
	c.cf = Vector{CapsFrame}()
	f = open(filename1, "r")
	for line in readlines(f)
		c.nb_times += 1
		push!(c.cf, parse(CapsFrame, line))
	end
	compute_cvel!(c)
	compute_permutations!(c)
	if length(filename2) > 0
		read_triangles!(c, filename2)
	end
	generate_plot_info!(c);
end

# ╔═╡ bb2fecf9-2bb7-415f-ad41-2a510a3e6daf
md"""
The function below plots the outline of a capsule, passed as the first argument, at a given frame, the index of which is passed as the second argument.
The function returns a plot containing the capsule outline in the plane $z = 0$ as well as the position of its centroid.
"""

# ╔═╡ 5e0094a6-b590-4a96-b48a-10290a295982
function plot_caps_outline(caps::Capsule, iter::Int)
	N = caps.cf[iter].nb_nodes
	npp = caps.npp
	xplane = Array{Float64}(undef, npp + 1)
	yplane = Array{Float64}(undef, npp + 1)
	xplane = [[caps.cf[iter].nodes[caps.perm[i]].x for i in 1:npp]; 
		caps.cf[iter].nodes[caps.perm[1]].x]
	yplane = [[caps.cf[iter].nodes[caps.perm[i]].y for i in 1:npp];
		caps.cf[iter].nodes[caps.perm[1]].y]

	p = plot(xplane, yplane, aspect_ratio=1.0, label="")
	p = plot!(scatter!([caps.cf[iter].centroid.x], [caps.cf[iter].centroid.y], 
		color="red", markersize=5.0, label="Centroid"))
	return p
end

# ╔═╡ ee492782-8c43-4d37-a626-5ba91ff2285a
md"""
## Single capsule analysis

### Reading the data of the capsules
"""

# ╔═╡ db900dca-28db-4676-ae2b-a3dc0ba7787b
md"""
First, we specify the simulations we want to visualize: which Capillary numbers are we interested in? Which Reynolds numbers? Alter the two cells below accordingly.
"""

# ╔═╡ 9ec3fff3-b7d3-4290-b68e-49418f420b34
Ca = [0.025, 0.05, 0.083, 0.12]

# ╔═╡ 257ae4c8-bf3e-4401-a708-5d2a92bafe64
Re = [0.01, 1.0, 25, 50]

# ╔═╡ 3f3705b2-d93a-4b4a-bd61-030b134335f7
md"""
We now read the position of all the nodes of the capsules we are interested in post-processing.
"""

# ╔═╡ b28bfb6b-1082-4345-9f44-8cb450c376b0
nca = length(Ca) 

# ╔═╡ ff00d3f0-15e3-4423-aa97-034cdd07dd2f
nre = length(Re)

# ╔═╡ f03efe73-0791-4813-95c8-e7835a384bcc
md"""
The function below ensures that parsing our `Float` Capillary and Reynolds number as Strings will match the names of the corresponding data files, e.g. that there is no trailing '.0'.
"""

# ╔═╡ 4047d010-3fb1-4a6c-964c-558efe816051
function format_name(s::AbstractString)
	if length(s) > 1
		if s[1] == '0'
			s_result = s[2:end]
		elseif s[end-1:end] == ".0"
			s_result = s[1:end-2]
		end
	end
	return s_result
end

# ╔═╡ 3e88e22e-7c57-46e0-8847-9375dd726d93
md"""
The data files are read for each combinations of Capillary and Reynolds numbers specified in the array `Ca` and `Re`. The resulting `Capsule` structures are stored in an array of arrays of `Capsules` (for some unknown reason I wasn't able to create a two-dimensional arrays of capsules: I don't know how to declare it without initializing it. This is left for future works :)).
"""

# ╔═╡ ad01eab3-0951-4975-9427-ba21da747aa9
mbs = [[Capsule() for i in 1:nre] for j in 1:nca];

# ╔═╡ ff974dab-dd03-4732-aa0c-df9fbeb6fe68
custom_filepath = @bind filepath TextField(default="data/");

# ╔═╡ 996e25ff-bd66-4776-b9c9-c4664176e49e
md"""
Default path of the folder in which the _folders_ containing the data files of the simulations are contained: $custom_filepath
"""

# ╔═╡ ad630b5b-c5df-4977-8568-7f106ac4c9e0
md"""
From the various capsule data that was just processed, we find the largest time duration of a simulation and we create a time array out of it. It will be our x-axis for most of the plots below
"""

# ╔═╡ 3310df40-a1ad-4065-a9c6-b0d99720d16d
ntimes = Matrix{Integer}(undef, nca, nre);

# ╔═╡ 710d4e0e-e7d7-4f12-964b-207d04f3e1e2
for i in 1:nca
	for j in 1:nre
		sre = string(Re[j])
		sre = format_name(sre)
		basename = filepath*"Cca"*string(Ca[i])[2:end]*"R"*sre*"/"
		read_capsule!(mbs[i][j],basename*"/mb_pos.csv")
		ntimes[i,j] = mbs[i][j].nb_times
	end
end

# ╔═╡ 4fc07215-aae7-47ca-9388-227d5ee748dc
maxiter = maximum(ntimes);

# ╔═╡ 9747ff93-ae58-4720-8afe-7260500d2e63
imaxiter, jmaxiter = convert(Tuple, argmax(ntimes));

# ╔═╡ eb136683-99af-4fc5-8485-0cac92d74c5b
times = [mbs[imaxiter][jmaxiter].cf[i].time for i in 1:maxiter];

# ╔═╡ 18ea582a-565c-4c90-aae0-afe13d74ce54
md"""
## Generating the outline of the capsule
"""

# ╔═╡ 01dacf9e-5ff9-444e-ae0a-4fc6ef00c9cb
md"""
The outline of the capsule in the plane $z = 0$ is shown only for one simulation, corresponding to the Capillary and Reynolds numbers located at indices ```ica``` and ```ire``` in their respective lists.
"""

# ╔═╡ c2ed42aa-0420-4e3f-be43-8d221d803dbc
custom_ica = @bind cica TextField(default="4");

# ╔═╡ f0595f36-4f5d-4620-8e66-272ff3741492
md"""
Index of the Capillary number of the capsule to display the outline of: $custom_ica
"""

# ╔═╡ 2dfd67ad-6dc6-4117-bc9f-9cfe5acfdfd7
ica = parse(Int64, cica)

# ╔═╡ 5a36ae91-a91f-46c3-9beb-00f08c6a7c3e
custom_ire = @bind cire TextField(default="1");

# ╔═╡ 87bb7986-f96e-4b9c-a043-4848d7b30cc8
md"""
Index of the Reynolds number of the capsule to display the outline of: $custom_ire
"""

# ╔═╡ 3e59eac6-bb05-4860-a520-90fb99c3e18b
ire = parse(Int64, cire)

# ╔═╡ 12858436-cf62-43ec-ad9c-767123022962
md"""
The last thing we need to do before plotting the capsule outline is to decide which time frame we want to visualize it at. We leave it as a(n interactive) choice for the user of this notebook, by defining the slider below.
"""

# ╔═╡ 1db76df8-817b-4aed-9286-8da5d3502b1e
iteration_slider = @bind iter Slider(1:mbs[ica][ire].nb_times, default=200);

# ╔═╡ 999d36c2-2b85-4ace-bb13-04498e1fb07b
Markdown.parse("""
The capsule outline for Ca=`$(Ca[ica])` and Re=`$(Re[ire])` is plotted below
""")

# ╔═╡ 9dd12fcf-834e-4792-a17e-f920662e4835
md"""
Control display frame: $iteration_slider
"""

# ╔═╡ 4ee1047c-ff72-4076-8103-f4c9cd0317db
p1 = plot_caps_outline(mbs[ica][ire], iter)

# ╔═╡ 7abca575-87a5-48fa-9a43-b7ac9e576a4d
md"""
## Plotting the normalized capsule velocity
"""

# ╔═╡ 0659e7d5-85d3-477d-aa98-73240f589b52
md"""
We plot the deviation of the normalized capsule velocity for the previously selected Capillary and Reynolds numbers.
"""

# ╔═╡ b6c49954-247f-4172-8d21-96b1746eb176
md"""
### Graphs for fixed Reynolds
"""

# ╔═╡ 22778a45-f773-44a3-bc01-15278d52e90f
md"""
The function below plots the normalized velocity curves of the centroid of a single capsule for all Capilarry numbers, and for a fixed Reynolds number.
"""

# ╔═╡ 9fcaa230-00ad-4ea9-aacd-09ae0568c8ef
function plot_vel_fixed_re(m, re, labels=["" for i in 1:4], title="", xlims=(-3, 9), ylims=(.85, 1.15))
	p = plot(bg=:white)
	n = size(m, 1) # number of curves
	for i in 1:n
		plot!(p, m[i][re].pi.times[2:end].-m[i][re].pi.t0,
			m[i][re].vel./m[i][re].pi.veq,
			label="Ca="*string(labels[i]),
			title=title)
	end
	xlims!(p, xlims)
	ylims!(p, ylims)
	return p
end

# ╔═╡ 6e3889e6-9aaf-4db6-a41b-6625f509c641
plot_vel_fixed_re(mbs, 1, Ca, "Re="*string(Re[1]))

# ╔═╡ 0ebdc692-3ee8-41f3-86ba-4eef96e25f8e
if nre > 1
	plot_vel_fixed_re(mbs, 2, Ca, "Re="*string(Re[2]))
end

# ╔═╡ 12dc155c-94d9-4f2b-91c4-10192bfde8d5
if nre > 2
	plot_vel_fixed_re(mbs, 3, Ca, "Re="*string(Re[3]))
end

# ╔═╡ 3f518649-1692-4f45-9df8-e4be84cdd116
if nre > 3
	plot_vel_fixed_re(mbs, 4, Ca, "Re="*string(Re[4]))
end

# ╔═╡ 9c8cccaa-8883-4dc2-9992-d4b592070c9e
md"""
### Graphs for fixed Capillary
"""

# ╔═╡ 9a1d0587-a2ff-475e-a33c-4d22f5083318
md"""
The function below plots the normalized velocity curves of the centroid of a single capsule for all Reynolds numbers, and for a fixed Capillary number.
"""

# ╔═╡ dac6695e-0196-4054-9629-a2445093063e
function plot_vel_fixed_ca(m, ca, labels=["" for i in 1:4], title="", xlims=(-3, 9), ylims=(.85, 1.15))
	p = plot(bg=:white)
	n = size(m[1], 1) # number of curves
	for i in 1:n
		plot!(p, m[ca][i].pi.times[2:end].-m[ca][i].pi.t0,
			m[ca][i].vel./m[ca][i].pi.veq,
			label="Re="*string(labels[i]),
			title=title)
	end
	xlims!(p, xlims)
	ylims!(p, ylims)
	return p
end

# ╔═╡ aff87583-e545-4293-a4a3-745b4c5dbe2d
plot_vel_fixed_ca(mbs, 1, Re, "Ca="*string(Ca[1]))

# ╔═╡ f67eefc4-b2a6-43fd-9580-4400944493c8
if nca > 1
	plot_vel_fixed_ca(mbs, 2, Re, "Ca="*string(Ca[2]))
end

# ╔═╡ c92b8978-be23-47bf-b345-0cc611eb2104
if nca > 2
	plot_vel_fixed_ca(mbs, 3, Re, "Ca="*string(Ca[3]))
end

# ╔═╡ b46984df-1eaa-413a-a3d8-ca74b961becb
if nca > 3
	plot_vel_fixed_ca(mbs, 4, Re, "Ca="*string(Ca[3]))
end

# ╔═╡ e47b2a95-75c7-4f94-ba23-035c68dc79ad
md"""
## Space convergence study
"""

# ╔═╡ 0ef1eabb-2fb3-4da8-86d1-7430d052ee24
l11 = Capsule();

# ╔═╡ 0adee9d3-73df-48e7-bad1-bacaeeaefbae
read_capsule!(l11, "data/Cca.12R50l11/mb_pos.csv");

# ╔═╡ c0c1c111-e7bc-48b6-adfb-04e3b4965663
l9 = Capsule();

# ╔═╡ 72d0aac0-0609-4810-864f-1a742d86e41e
read_capsule!(l9, "data/Cca.12R50l9lag3/mb_pos.csv");

# ╔═╡ d9ad599f-f0d5-41f9-ad87-c7b47d45b23f
begin
	pconv = plot(bg=:white)
	plot!(pconv, mbs[4][4].pi.times[2:end].-mbs[4][4].pi.t0,
			mbs[4][4].vel./mbs[4][4].pi.veq, label="level 10")
	plot!(pconv, l11.pi.times[2:end].-l11.pi.t0,
			l11.vel./l11.pi.veq, label="level 11")
	xlims!(pconv, (-3, 5.4))
	ylims!(pconv, (.6, 1.4))
	xlabel!(pconv, "time")
	ylabel!(pconv, "Centroid velocity deviation")
	title!(pconv, "Convergence study, Ca=0.12, Re=50")
end

# ╔═╡ 388d3a92-2064-460e-b4e7-cfa805b89b3f
md"""
## Surface area dilatation
"""

# ╔═╡ 620b0843-6389-4e32-ab6c-0a8ca1c498e0
begin
	a = Capsule()
	read_capsule!(a, "data/Cca.12R.01/mb_pos.csv", "data/Cca.12R.01/mb_tri.csv")
	
	b = Capsule()
	read_capsule!(b, "data/Cca.12R1/mb_pos.csv", "data/Cca.12R1/mb_tri.csv")
end

# ╔═╡ 383f5087-f619-46c2-bf1a-26c02f7a49ad
begin	
	parea = plot(bg=:white)
	plot!(parea, a.pi.times[1:end-1].-a.pi.t0, a.area[1:end-1]./(4*pi),
		label = "Re = 0.01")
	plot!(parea, b.pi.times[1:end-1].-b.pi.t0, b.area[1:end-1]./(4*pi),
		label = "Re = 1")
	parea
end

# ╔═╡ 7146f184-13aa-43cb-967b-d99d53e51c04
md"""
## Binary capsule analysis
"""

# ╔═╡ e7ca3688-2586-4d38-b881-1ed83db1b03d
md"""
Coming soon...
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
Plots = "~1.38.4"
PlutoUI = "~0.7.49"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "f433cc62030fa516dae9d8996e9de991da8944d1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "844b061c104c408b24537482469400af6075aae4"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.5"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random", "SnoopPrecompile"]
git-tree-sha1 = "aa3edc8f8dea6cbfa176ee12f7c2fc82f0608ed3"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.20.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "61fdd77467a5c3ad071ef8277ac6bd6af7dd4c04"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "9e23bd6bb3eb4300cb567bdf63e2c14e5d2ffdbc"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.71.5"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "aa23c9f9b7c0ba6baeabe966ea1c7d2c7487ef90"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.71.5+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "37e4657cd56b11abe3d10cd4a1ec5fbdb4180263"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.7.4"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "2422f47b34d4b127720a18f86fa7b1aa2e141f29"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.18"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "45b288af6956e67e621c5cbb2d75a261ab58300b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.20"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "6503b77492fd7fcb9379bf73cd31035670e3c509"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.3.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6e9dba33f9f2c44e08a020b0caf6903be540004"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.19+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "8175fc2b118a3755113c8e68084dc1a9e63c61ee"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.3"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "c95373e73290cf50a8a22c3375e4625ded5c5280"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.4"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SnoopPrecompile", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "87036ff7d1277aa624ce4d211ddd8720116f80bf"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.38.4"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eadad7b14cf046de6eb41f13c9275e5aa2711ab6"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.49"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "261dddd3b862bd2c940cf6ca4d1c8fe593e457c8"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.3"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase", "SnoopPrecompile"]
git-tree-sha1 = "e974477be88cb5e3040009f3767611bc6357846f"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.11"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "94f38103c984f89cf77c402f2a68dbd870f8165f"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.11"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "ac00576f90d8a259f2c9d823e91d1de3fd44d348"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "ed8d92d9774b077c53e1da50fd81a36af3744c1c"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╟─47c2a121-b626-4f8b-9722-cefa925ba9a5
# ╠═c7ab89d0-638e-4121-842e-b14a7d58a4b2
# ╠═720128f9-b311-4883-9444-f7d564008f66
# ╠═84fb4b5c-8fbd-11ed-165d-b758c56cbf89
# ╟─2012189e-7f9c-4a59-85cc-c6dd4b04d129
# ╟─9d776581-c80d-4760-9ac2-e01294a60c5e
# ╠═6e5461ea-fae3-4c21-abf2-d2aae9dc64ba
# ╟─4cb4b887-a081-43e2-94eb-bdac56cc405a
# ╠═56710fed-1eff-4263-a48d-6dacc371bb86
# ╟─c74548f1-ac38-4b5f-a990-df5dbf01b7a8
# ╠═e5e71a51-fe18-42fe-9acb-49f1147845f6
# ╠═111b0b5b-f40b-41f1-a8b8-bc8fcac65c4e
# ╠═0366ffe5-9f64-46f0-9aad-7d19c779a2e3
# ╠═5695e259-55fc-40c7-af3d-04447486f71b
# ╟─853ca926-d374-46a7-bdad-9141e506cd79
# ╠═febb23c9-4924-496b-a8e5-c881f26c7d8c
# ╟─cfaa17a3-5eaf-479e-bae3-c55608be8d20
# ╠═1cfe7a07-411b-4a05-a489-a1739d9a0702
# ╟─7b25dd66-e499-4a66-b6bd-ac7b9663a603
# ╠═78619fd5-866c-4153-a40e-0ec9f9ed8f4a
# ╟─32a4b2f1-4299-4549-afa9-3be80ce19ffd
# ╠═81e03b17-9741-4c66-8047-cce319f3900e
# ╟─17e704a2-39ec-412a-a552-abcd04fc1a26
# ╠═91daca3f-575b-4f71-83a7-8abde506e1fe
# ╟─f1f1f88f-4a8d-44d3-a516-39e6fa40e44f
# ╠═2b3d71a8-ac2a-4b64-b547-661646343995
# ╟─91be20f2-79dd-4641-a4f5-04306dc341fd
# ╠═cbadf507-94a9-47bc-93b6-8512235ea715
# ╟─1c81d5f7-6506-40d8-a9c9-19b11cb901dc
# ╠═5d7d1025-fcad-4cd9-9bb6-a2a3f033772c
# ╟─3ee686a6-81d5-419a-9d0c-41719c4fe2bd
# ╠═1fcaa2aa-6876-44c0-8db2-35d00211bff4
# ╟─e118297a-fcb8-409c-bb8c-e1992833fe85
# ╠═0c86edd5-fb43-4949-9ace-06f923daec2f
# ╠═d582cf21-202e-4fb0-af4e-cf73f44a4fa0
# ╟─74454cc5-8be1-4e6d-b2d7-e160b785e91b
# ╠═f86282ba-646f-4d89-92c5-c6fbd4af4748
# ╟─bb2fecf9-2bb7-415f-ad41-2a510a3e6daf
# ╠═5e0094a6-b590-4a96-b48a-10290a295982
# ╟─ee492782-8c43-4d37-a626-5ba91ff2285a
# ╟─db900dca-28db-4676-ae2b-a3dc0ba7787b
# ╠═9ec3fff3-b7d3-4290-b68e-49418f420b34
# ╠═257ae4c8-bf3e-4401-a708-5d2a92bafe64
# ╟─3f3705b2-d93a-4b4a-bd61-030b134335f7
# ╠═b28bfb6b-1082-4345-9f44-8cb450c376b0
# ╠═ff00d3f0-15e3-4423-aa97-034cdd07dd2f
# ╟─f03efe73-0791-4813-95c8-e7835a384bcc
# ╠═4047d010-3fb1-4a6c-964c-558efe816051
# ╟─3e88e22e-7c57-46e0-8847-9375dd726d93
# ╠═ad01eab3-0951-4975-9427-ba21da747aa9
# ╠═ff974dab-dd03-4732-aa0c-df9fbeb6fe68
# ╟─996e25ff-bd66-4776-b9c9-c4664176e49e
# ╠═710d4e0e-e7d7-4f12-964b-207d04f3e1e2
# ╟─ad630b5b-c5df-4977-8568-7f106ac4c9e0
# ╠═3310df40-a1ad-4065-a9c6-b0d99720d16d
# ╠═4fc07215-aae7-47ca-9388-227d5ee748dc
# ╠═9747ff93-ae58-4720-8afe-7260500d2e63
# ╠═eb136683-99af-4fc5-8485-0cac92d74c5b
# ╟─18ea582a-565c-4c90-aae0-afe13d74ce54
# ╟─01dacf9e-5ff9-444e-ae0a-4fc6ef00c9cb
# ╠═c2ed42aa-0420-4e3f-be43-8d221d803dbc
# ╟─f0595f36-4f5d-4620-8e66-272ff3741492
# ╠═2dfd67ad-6dc6-4117-bc9f-9cfe5acfdfd7
# ╠═5a36ae91-a91f-46c3-9beb-00f08c6a7c3e
# ╟─87bb7986-f96e-4b9c-a043-4848d7b30cc8
# ╠═3e59eac6-bb05-4860-a520-90fb99c3e18b
# ╟─12858436-cf62-43ec-ad9c-767123022962
# ╠═1db76df8-817b-4aed-9286-8da5d3502b1e
# ╟─999d36c2-2b85-4ace-bb13-04498e1fb07b
# ╟─9dd12fcf-834e-4792-a17e-f920662e4835
# ╠═4ee1047c-ff72-4076-8103-f4c9cd0317db
# ╟─7abca575-87a5-48fa-9a43-b7ac9e576a4d
# ╟─0659e7d5-85d3-477d-aa98-73240f589b52
# ╟─b6c49954-247f-4172-8d21-96b1746eb176
# ╟─22778a45-f773-44a3-bc01-15278d52e90f
# ╠═9fcaa230-00ad-4ea9-aacd-09ae0568c8ef
# ╠═6e3889e6-9aaf-4db6-a41b-6625f509c641
# ╠═0ebdc692-3ee8-41f3-86ba-4eef96e25f8e
# ╠═12dc155c-94d9-4f2b-91c4-10192bfde8d5
# ╠═3f518649-1692-4f45-9df8-e4be84cdd116
# ╟─9c8cccaa-8883-4dc2-9992-d4b592070c9e
# ╟─9a1d0587-a2ff-475e-a33c-4d22f5083318
# ╠═dac6695e-0196-4054-9629-a2445093063e
# ╠═aff87583-e545-4293-a4a3-745b4c5dbe2d
# ╠═f67eefc4-b2a6-43fd-9580-4400944493c8
# ╠═c92b8978-be23-47bf-b345-0cc611eb2104
# ╠═b46984df-1eaa-413a-a3d8-ca74b961becb
# ╟─e47b2a95-75c7-4f94-ba23-035c68dc79ad
# ╠═0ef1eabb-2fb3-4da8-86d1-7430d052ee24
# ╠═0adee9d3-73df-48e7-bad1-bacaeeaefbae
# ╠═c0c1c111-e7bc-48b6-adfb-04e3b4965663
# ╠═72d0aac0-0609-4810-864f-1a742d86e41e
# ╠═d9ad599f-f0d5-41f9-ad87-c7b47d45b23f
# ╟─388d3a92-2064-460e-b4e7-cfa805b89b3f
# ╠═620b0843-6389-4e32-ab6c-0a8ca1c498e0
# ╠═383f5087-f619-46c2-bf1a-26c02f7a49ad
# ╟─7146f184-13aa-43cb-967b-d99d53e51c04
# ╟─e7ca3688-2586-4d38-b881-1ed83db1b03d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
