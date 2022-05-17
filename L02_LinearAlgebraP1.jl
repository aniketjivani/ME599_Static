### A Pluto.jl notebook ###
# v0.17.7

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

# ╔═╡ 51a640f1-7007-44ac-9751-655083073fab
begin
	# We activate our project environment, which disables the notebook's built in package manager.
	using Pkg
	Pkg.activate("../Project.toml")
	using LinearAlgebra
	using PlutoUI
	using Colors
	using Plots
	gr()
	nbsp = html"&nbsp"
	lbreak = html"<br>"

end

# ╔═╡ 8a09a85c-cfcd-4478-ae21-7852227b3f61
begin
	using Images
end

# ╔═╡ e4277174-73eb-4925-a104-21fa937ae3b6
using Random

# ╔═╡ a876703f-6e85-488b-adb2-c566073f304d
html"""
	<article class="title">
	<p>
	MECHENG 599-001 Computational & Data-Driven Methods in Engineering
	</p>
	</article>
	
	<article class="instructor">
	<p>
	Instructor: Xun Huan
	</p>
	</article>
	

	<article class="lecture">
	<p>
	Lecture 02: Linear Algebra Review 1
	</p>
	</article>


	<style>
	
	article.title p{
		
		color:cyan;
		font-size: 30px;
		font-weight: bold;
		text-align: center;
	}

	article.instructor p{
		color:grey;
		font-size: 18px;
		font-weight: bold;
		text-align: center;
	}


	article.lecture p{
		color:grey;
		font-size: 25px;
		font-weight: bold;
		text-align: center;
	}

	</style>
"""

# ╔═╡ ae321ee5-c21f-4d5a-b1b7-4b2f0e4d8cd9
begin
	md"""
	**Delete this cell afterwards**: there are blank spaces created throughout this notebook with the following code:
	
	```julia
	md"
	$lbreak
	$lbreak
	$lbreak
	$lbreak
	$lbreak
	$lbreak
	$lbreak
	"
	```
	where `lbreak = html"<br>"` - so that if this notebook is exported as a pdf, that space can be used for drawing and explaining. If more blank spaces are needed, the above code can be copy pasted in a new cell, with cell visibility set to hidden. (click the eye that appears in the top left when hovering over a cell)
	
	"""
end

# ╔═╡ 6f622d7e-0684-11ec-0e37-dbc7eccb7a73
md"""
# Matrices and vectors

Key points are summarized below:

"""

# ╔═╡ 80fd4f57-dfd8-4e1d-a5c8-67bf6cbfb281
begin
	md"""
	
	```julia
	julia> y = 10
	10
	
	julia> typeof(y)
	Int64
	```
	"""
end

# ╔═╡ afadb596-8200-47c1-933d-b4d60dced054
md""" Let us define a variable as follows - """

# ╔═╡ 5fa5bbfb-27ee-4621-b66b-7f05ab37fb2d
y = 10

# ╔═╡ d5b726b0-7aef-4c35-8fcf-55be19ec8ffa
md"""
Hitting `Ctrl + Enter` (or `Command + Enter` on a Mac) executes the cell you are currently in and adds a new cell below. To see the full list of shortcuts in Pluto, hit `F1` on your keyboard. 
"""

# ╔═╡ d347889f-15bf-4459-8411-cfc00058afc3
typeof(y)

# ╔═╡ 1f40c451-c984-46da-9a35-337b97e254ce
im # `im` represents the imaginary unit i for complex numbers

# ╔═╡ bcf0eee2-0a56-43c5-836e-80a29a552e01
z = 3 + 4im

# ╔═╡ 40e0d8b2-6f34-4dcd-bf3c-6e2ca30422f5
typeof(z)

# ╔═╡ b7a25f65-92cc-446f-8f76-5d1d1774a83d
md"""
* The complex conjugate of a scalar $z$ is written as $\bar{z}$ or $z^{*}$, obtained by negating its imaginary part. If $z$ is a real number, then $\bar{z} = z$.
"""

# ╔═╡ e67dffdf-cc8f-4ee3-93cc-a863a61bc3f3
z'

# ╔═╡ 74cad985-f418-41b5-8f78-84e7f5475838
md"""
Note that using `transpose(z)` will not return the complex conjugate - it is safer to use `'` everywhere.
"""

# ╔═╡ 2026b715-c4c9-4b91-82da-96eed87f3d39
transpose(z)

# ╔═╡ 58ff8911-289a-4d31-ad87-c6e94525dbb0
y' == y

# ╔═╡ d9d3e037-dc13-45e5-a97f-408122d512a0
md"""
## Notation

* We will use upper case to denote a matrix, and lower case to denote a (column) vector, e.g., $A \in \mathbb{C}^{m\times n}$ and $x \in\mathbb{C}^{n}$. Each column of $A$ is denoted as $a_j, j=1,\ldots,n$; each row of $A$ is $a^{\ast}_i,i=1,\ldots,m$; and each element of $A$ is denoted as $a_{ij},i=1,\ldots,m$ and $j=1,\ldots,n$. Each element of the vector $x$ is denoted by $x_j,j=1,\ldots,n$. Context will often help. Using this notation, we can write

$$\begin{align}
 A = \left[\begin{array}{c|c|c|c} &&& \\ &&& \\
a_1 & a_2 & \cdots & a_n \\ &&& \\ &&& \end{array}\right]
= \left[\begin{array}{ccccc} & & a_1^{\ast} & & \\\hline && a_2^{\ast} && \\\hline && \vdots && \\\hline && a_m^{\ast} && \end{array}\right], 
\qquad x = \begin{bmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{bmatrix}.
\end{align}$$


"""

# ╔═╡ e2e120fe-778e-44e1-8337-37c9e57d1951
md"""
* The matrix-vector product $b=Ax$ is then written as:

$$\begin{align}
    b_i = \sum_{j=1}^{n} a_{ij} x_j, \quad i=1,\ldots,m.
\end{align}$$

"""

# ╔═╡ 8ef28a47-417f-45a8-acc7-76547a63da9a


# ╔═╡ b975fb26-094c-4587-9b03-1eef425b4cff
md"""

Equivalently, it can also be expressed as

$$\begin{align}
    b = \sum_{j=1}^{n} a_{j} x_j
    = x_1 \begin{bmatrix} \\ \\ a_1 \\ \\ \\ \end{bmatrix} 
    + x_2 \begin{bmatrix} \\ \\ a_2 \\ \\ \\ \end{bmatrix}
    + \ldots
    + x_n \begin{bmatrix} \\ \\ a_n \\ \\ \\ \end{bmatrix}.
\end{align}$$


This interpretation is very useful:
$b$ is being expressed as a linear combination of the columns $a_j$ (analogy: $a_j$ are the ingredients, and $x_j$ is the recipe telling us what proportions to combine these ingredients in order to produce the final dish $b$). Instead of viewing $Ax=b$ as $A$ operating on $x$ to produce $b$, can think as $x$ acting on $A$ to produce $b$.

"""

# ╔═╡ 894ed0d4-941e-46a2-8b9a-357bbd30fde9
md"""
We can create our own dish below by changing the ingredients in $x$ 😉
"""

# ╔═╡ f98ba8f4-2ce6-47d8-9231-d55b178385a7
begin
	a₁ = [3;
		  -1]

	a₂ = [-4;
		  3]

	a₃ = [0;
		  2]

	a₄ = [1;
		  1]

end

# ╔═╡ 9e15631a-a0bb-4e1a-b622-1f4d2d8b789a
md"""

x₁ = $(@bind x₁ Slider(1:15, default=5, show_value=true)) 
$nbsp 
$nbsp 
x₂ = $(@bind x₂ Slider(3:20, default=10, show_value=true))

$lbreak
$lbreak

x₃ = $(@bind x₃ Slider(-5:10, default=-3, show_value=true)) 
$nbsp 
$nbsp 
x₄ = $(@bind x₄ Slider(2:8, default=6, show_value=true))


"""

# ╔═╡ 12e6e9e0-629d-4fef-abc5-1ce917f5cf6d
b = x₁ * a₁ + x₂ * a₂ + x₃ * a₃ + x₄ * a₄

# ╔═╡ a36c5035-194c-4c98-b120-5481ac10fc25
md"""
* Interpretation of matrix-matrix multiplication $B=AC$  (here $A$ is $\ell\times m$, $C$ is $m\times n$, and $B$ is $\ell \times n$):

$$\begin{align}
    b_{ij} = \sum_{k=1}^{m} a_{ik} c_{kj}.
\end{align}$$

Break this down into $n$ matrix-vector multiplications $b_j=Ac_j$ (analogy: each column of $B$ is a dish formed by following a recipe $c_j$ that tells us what proportions to combine together the ingredients which are the columns of $A$).

"""

# ╔═╡ ce3c3a23-7edd-4697-918d-222a78291538
md"
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
"

# ╔═╡ d69540d9-2647-48f3-9b80-3a6f72f19035
md"""
# The fundamental subspaces

* The *range* of $A$, written $\textrm{range}(A)$, is the set of vectors that can be expressed as $Ax$ for some $x$ (analogy: $\textrm{range}(A)$ is the set of all possible dishes that can be produced by different combinations of the ingredients $a_j$). $\textrm{range}(A)$ is the space spanned by the columns of $A$, thus also called the *column space* of $A$.

"""

# ╔═╡ 8c6db3b1-a11a-42ff-809d-21e2dff6b1ba
md"""
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
"""

# ╔═╡ e29cd04a-c020-4c7a-8fda-39a0126fb83f
md"""
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
"""

# ╔═╡ 802b8714-3f01-4d80-a6ed-b28d432839c4
# Simple exercise - look at the columns of A and check if a particular b can belong to the column space of A. If not, what would be the column space?

# ╔═╡ 1903e8b6-f8d3-4ed5-93ba-4aaca44608c9
A = [3 0;
	 0 -2] 

# ╔═╡ 4e3801f4-cb8f-4d45-8e8c-d2556bdbcbc8
bb = [5
	 3]

# ╔═╡ a48392f9-2eb2-46f2-af6a-b032816390c3
md"""
Does $bb \in \textrm{range}(A)$? ⟹ Is there any possible linear combination of $A$ that can produce $bb$? 
"""

# ╔═╡ 1b439146-1fe0-4756-ae3d-d23982aa3e01
# Play around with values of x1 and x2 till btil coincides with b

# ╔═╡ a89bc125-a20c-40d6-87c3-4b9ba6754665
# PSA: Slightly buggy, type negative sign by typing number and moving cursor to left
# or avoid entering value altogether by using up and down arrows in box
# change values of x1 and x2 till the box below turns green.

# ╔═╡ 77c7cc1e-cd25-4dcd-86e5-8c615bde9017
md"""
x1 = $(@bind x1 NumberField(-100:0.01:100, default=1))
$lbreak
$lbreak
x2 = $(@bind x2 NumberField(-100:0.01:100, default=2))
"""

# ╔═╡ 16f0cd08-6cb9-446f-876d-977430fd6705
x = [x1;
	 x2]

# ╔═╡ 514b975e-0fb3-4943-8596-37dc99be3c10
btil = A*x

# ╔═╡ 25ae9da1-c85e-490b-a0da-693cf9cd5499
md"""
!!! warning "Example"
    What if we now change $A(2, 2)$ to $0$? Would $bb$ still $\in \textrm{range}(A)$? What would be the new range?
"""

# ╔═╡ 0f370f2c-0371-40d0-a4d7-6979c81dfbed
# Correct solution for A = [3 0; 0 -2]
xx = A\bb

# ╔═╡ 920ed5b7-cda9-49fc-bac6-91a2cf36f79d
begin
	if isapprox(x, xx, atol=0.01)
		md"""
!!! correct "You got this!"
    Since we are able to find an $x$ that satisfies $Ax = bb$ we can say that $bb \in \textrm{range}(A)$
		"""
	else
		md"""
!!! danger "Close enough..."
    Need to try again :)
		"""
	end
	
end

# ╔═╡ 3e7052ea-743b-4c02-b266-b90b402d87c2
md"""
* The _nullspace_ of $A$, written $\textrm{null}(A)$, is the set of vectors satisfying $Ax=0$.
"""

# ╔═╡ f5ec28f0-b9ac-478a-adda-bf641f868efa
md"""
* A _nonsingular_ or _invertible_ matrix is a square matrix with full rank. For an invertible $m\times m$ matrix $A$, its $m$ columns span the entire $\mathbb{C}^m$. Therefore, we can expand out the canonical unit vector 
$$e_j=Az_j = \sum_{i=1}^{m} z_{ij} a_j$$

through some "recipe" vector $z_j$. Doing this for all unit vectors $j=1,\ldots,m$ and collecting the $z_j$ vectors, then $Z=[z_1,\ldots,z_m]$ is the inverse matrix of $A$.
"""

# ╔═╡ b2466928-15a6-4626-9acb-58f1b3ff7e74
md"""
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
"""

# ╔═╡ 0220b07d-abc1-4c38-acd7-eaaa3148bd5d
md"""
* Consider $x=A^{-1}b$. Earlier we already established $A^{-1}b$ (which equals $x$ here) is the vector of coefficients of the expansion of $b$ (the dish) in the basis of columns of $A$ (the ingredients). Hence, multiplication by $A^{-1}$ is a _change of basis_ operation representing a vector from the "standard" coordinates to the "$A$-columns" coordinates! Please study and fully understand the diagram on page 9 of Trefethen and Bau (1997).
"""

# ╔═╡ 8c43733c-2786-4d31-94f4-0d190fcd3dca
PlutoUI.Resource("https://i.imgur.com/HTO18vV.jpg")

# ╔═╡ 36e25a97-78b1-437d-ab98-295166daba5b
md"""
# Rank

* The _column rank_ of $A$ is the dimension of its columns space; the _row rank_ is the dimension of the space spanned by its rows. Column rank and row rank are always equal, often simply called the _rank_ of the matrix. A matrix of full rank means it has the highest possible rank (i.e. the lesser of $m$ and $n$). If $m\geq n$, then that matrix maps no two distinct vectors to the same vector (one-to-one).

"""

# ╔═╡ 02b49d9c-098b-45ec-8837-347ec3e2e9f0
md"""
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
"""

# ╔═╡ 0eafbb01-c078-4ea0-914f-0412ab6c890a
md"""
As a matter of fact, for simple, structured arrays, we can even identify the rank by means of visual inspection. Gilbert Strang's "Introduction to Linear Algebra" features a pretty cool section devoted to finding the ranks of flags. Let's look at a few!!🇩🇰🇧🇯
"""

# ╔═╡ fc2770f0-425b-4b07-8315-070f43ea30cf
PlutoUI.Resource("https://i.imgur.com/E1D0q7S.jpg")

# ╔═╡ 9cd171cd-82b8-4208-b6fa-94d298a6e173


# ╔═╡ e304923f-9ab2-4a43-abd5-95ebb84d023c
md"""
!!! correct "Images"
    A crude algorithm for image compression relies on us being able to use Singular Value Decomposition (SVD) to approximate our original image by a low-rank version. While it is not the most effective algorithm in practice for the task. However, it goes to show how central the idea of rank is to the methods we are interested in. 

In fact, when operating on multi-dimensional arrays ("tensors"), a major challenge is determining the rank of the array uniquely, which limits our ability to generalize algorithms that are so well established for matrices. Read more about tensor decompositions here!! [Tensor Decompositions and Applications (Kolda)](https://www.kolda.net/publication/TensorReview.pdf)

""" 

# ╔═╡ b0832388-5c01-40c0-ad79-f368b646abc9
md"""
# Orthogonality

Please read Lecture 2 of Trefethen and Bau (1997). Key points are summarized below:
"""

# ╔═╡ 31ab7c4f-335d-46a4-9688-cbda0feec3ec
md"""
* The _hermitian conjugate_ or _adjoint_ of an $m\times n $ matrix $A$, written as $A^{*}$, has $i,j$ entry equal to the complex conjugate of the $j,i$ entry of $A$. For example
$$A = \begin{bmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \\ a_{31} & a_{32} \end{bmatrix},
    \qquad
    A^{*} = \begin{bmatrix} \bar{a}_{11} & \bar{a}_{21} & \bar{a}_{31} \\ \bar{a}_{12} & \bar{a}_{22} & \bar{a}_{32} \end{bmatrix}$$
If the matrix is real, then the hermitian is simply the matrix transpose, written $A^{T}$.
"""

# ╔═╡ 98b1b9dc-0f37-4f3e-91bf-30be8e957f34
[2 + 3im 4 - im;
 1 + 5im -6]'

# ╔═╡ 2a9d5631-8176-44f0-beb1-3fda3058d56f
md"""
* The _inner product_ of two vectors $x,y \in \mathcal{C}^{m}$ is
$$x^{* }y = \sum_{i=1}^{m} \bar{x}_i y_i.$$

The Euclidean length of $x$, written as $\|x\|$, is

$$\|x\| = \sqrt{x^{\ast}x} = \left(\sum_{i=1}^{m} |x_i|^2\right)^{\frac{1}{2}}.$$
The inner product is bilinear (linear in each vector separately).
"""

# ╔═╡ b6c98f01-51b0-47e3-889d-d2b0eb7e767f
a = [2, 3, 5]

# ╔═╡ ddee31eb-56a4-4240-9c90-a484d2dea11f
typeof(a)

# ╔═╡ 810e9ddc-41a0-4c90-80bf-84e978ba1bb9
Array{Int64, 1} == Vector{Int64} # Vector is an alias for array of dim 1

# ╔═╡ 3a77c857-8a57-49ab-bcf5-79c1dece416e
c = [3, -1, 0]

# ╔═╡ d811dd3d-c419-491e-97f4-35f7159f6a6a
a'c # or a' * c

# ╔═╡ 0c1814e3-4b62-4434-8868-598eb7127836
a'a

# ╔═╡ e6588752-18c3-4503-a427-f0b4b8246c0d
norm(a) == sqrt(a'a)

# ╔═╡ c41d4e8e-58b8-4eff-95ce-88aa50951059
md"""
We will learn that the Euclidean length is a kind of norm in the next section. The function `norm` accepts two arguments. Looking it up in Live Docs shows us the following: 
```julia
norm(A, p::Real=2)
```
where `p` takes the value 2 if the second argument is not supplied. 
"""

# ╔═╡ e9ab39f4-21aa-4929-8608-3a1c577136e9
md"""
* Some useful properties:
$$\begin{align} 
(AB)^{\ast} &= B^{\ast} A^{\ast}\\ 
(AB)^{-1} &= B^{-1} A^{-1} \\ 
A^{-\ast} &:= (A^{\ast})^{-1} = (A^{-1})^{\ast}
\end{align}$$
"""

# ╔═╡ 0b20a6b4-7a38-4f3c-9c4a-e2268e417898
md"""
* Two vectors $x$ and $y$ are _orthogonal_ if $x^{\ast}y=0$. Two sets of vectors $X$ and $Y$ are orthogonal if every $x\in X$ is orthogonal to every $y\in Y$. A set of nonzero vectors $S$ is said to be orthogonal if for $x,y\in S,x\neq y,\Rightarrow x^{\ast}y=0$. The set is _orthonormal_ if every $x\in S$ has $\|x\|=1$.
"""

# ╔═╡ 23cf0c46-1b15-4659-a6fb-17c036bb5e5c
aa = [-1, 1, 0]

# ╔═╡ c47ea5f8-2282-4fa6-b387-c4116cff0b85
cc = [1, 1, 0]

# ╔═╡ 8d3635ea-f227-41e1-b524-9a80213f2fdf
norm(aa)

# ╔═╡ 441b7b2f-0b4a-409f-9672-b5bf103eef44
aa'cc

# ╔═╡ 2b37dac7-7dbf-40c8-826d-cde1800ac896
md"""
* Vectors in an orthogonal set $S$ are linearly independent.
"""

# ╔═╡ 890bfe38-bf37-45f9-a45f-3a665152747a
md"""
* Inner products can be used to decompose arbitrary vectors into orthogonal components. Let $\{q_1,q_2,\ldots,q_m\}$ be an orthonormal set (i.e. $q_i^{\ast} q_i=1$ and $q_i^\ast q_j=0$ for $i \neq j$) that spans $\mathbb{C}^{m}$, then an arbitrary vector $v$ can be decomposed into $m$ orthogonal components 
$$\begin{align}
    v = \sum_{i=1}^{m} (q_i^{\ast} v) q_i = \sum_{i=1}^{m} (q_i q_i^{\ast}) v.
\end{align}$$
"""

# ╔═╡ 3913c8f5-e3a4-404d-989a-024291f3912f
md"""
In the first form, since $q_i^{\ast}v$ is a scalar, we can interpret $v$ being expanded as a scalar combination of vectors $q_i$. In the second form, we can interpret $q_i q_i^{\ast}$ (a rank-one matrix) as a projection operation of $v$ onto the direction $q_i$.
"""

# ╔═╡ e54eb17d-d4d3-4d35-bebf-b0a488b9b56a
md"""
* A square matrix $Q\in \mathbb{C}^{m\times m}$ is _unitary_ if $Q^{\ast}=Q^{-1}$, or equivalently $Q^{\ast}Q=I$. If $Q$ is a real matrix, then it is known as an orthogonal matrix.
"""

# ╔═╡ fbb084c7-6daf-46f8-851b-cc6bd7dc236f
md"""
  $Qx$ & $Q^{-1}b$ have the same "change-of-basis" interpretation as the previous section. Furthermore, multiplication by unitary matrix (or its adjoint) preserves (Euclidean) geometric structure because the inner product is unchanged under such a change of basis: $(Qx)^{\ast} (Qy) = x^{\ast}y$. Consequently, the norm is also preserved: $\|Qx\| = \|x\|$. If $Q$ is real, it corresponds to a rigid rotation (rotation matrix) or reflection of the vector space.
"""

# ╔═╡ d6d3ce82-4211-483d-8353-0fd3a7756321
md"""

# Concept of distances and norms

A (vector) _norm_ is a function $\|\cdot\|: \mathbb{C}^{m} \mapsto \mathbb{R}$ that assigns a real-valued length to each vector. A norm must satisfy the following three properties:

-  $\|x\|\geq 0$, and $\|x\|=0$ only if $x=0$
-  $\|x + y\| \leq \|x\| + \|y\|$
-  $\|\alpha x\| = |\alpha| \|x\|$.


Note that in a normed metric space $(X, \|\cdot\|)$, $d : X \times X \mapsto \mathbb{R}$ defined as $d(x, y) = \|x-y\|$ is the _distance_ metric on X. 

There are a variety of norms and distances defined in the literature, which prove to be useful in different settings. `Distances.jl`, a Julia package, supports the following!
"""

# ╔═╡ fa8fed50-c898-44f6-aa38-fb211e89e62a
PlutoUI.Resource("https://i.imgur.com/RC0Tn6K.jpg")

# ╔═╡ f66c61bf-eca1-47d1-b8f3-928459daab08
md"""
!!! warning "Note"
    Not all the metrics above fall into the definition of "distance" and "norm" that we have explained above. Divergence in particular, while not technically a distance, is still very useful, as we will get to know when we dive into obtaining probability distributions in future problems!
"""

# ╔═╡ a0f85f74-a6d8-4efd-b295-3437be3dbf78
md"""
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
"""

# ╔═╡ 9ddeb1b1-7436-4d02-b772-9ed1ee5c513e
md"""
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
"""

# ╔═╡ 97c03b72-f8a5-40e6-841a-d0ebd4a40c30
md"""
_Example:_ Vector $p$-norm

$$\begin{align}
\|x\|_1 &= \sum_{i=1}^{m} |x_i| \\
\|x\|_2 &= \(\sum_{i=1}^{m} |x_i|^2 \)^{\frac{1}{2}} \\
\|x\|_{\infty} &= \max_{1\leq i \leq m} |x_i| \\
\|x\|_p &= \(\sum_{i=1}^{m} |x_i|^p \)^{\frac{1}{p}}, \quad 1\leq p < \infty.
\end{align}$$
"""

# ╔═╡ 52700936-6738-4731-b643-335f73fa2111
md"""
_Example:_ Weighted-norm (of a 2-norm / Euclidean norm) weighed by a (here assuming diagonal) matrix $W$:

$$\begin{align}
\left\|x\right\|_{W} = \left\|Wx\right\| = \left(\sum_{i=1}^{m} |w_i x_i|^2 \right)^{\frac{1}{2}}.
\end{align}$$

"""

# ╔═╡ eeedaee2-7f68-4431-aaab-f3f0962455f9
md"""
**Review**:
Weighted norms do have their uses, especially when the magnitudes of individual entries of $x$ fluctuate wildly, and a particular entry can disproportionately affect the value of the norm. In that case, a diagonal matrix like the one described above can help "normalize" the entry. Here's a helpful post that elaborates on this! [Weighted Inner Product](https://math.stackexchange.com/a/396845)

Using a weighted norm will turn out to be useful in a clustering application too, down the line!

"""

# ╔═╡ dced6ea7-55e8-43f6-8652-69bdabecfa6b
md"""
## Matrix norms

A _matrix norm_ for a $m\times n$ matrix satisfies the same three conditions as the vector norm but in the $mn$-dimensional vector space. 

_Example:_ Frobenius norm is essentially the 2-norm but viewing the matrix as a $mn$-dimensional vector:
$$\begin{align}
    \norm{A}{F} = \left( \sum_{i=1}^{m} \sum_{j=1}^{n} |a_{ij}|^2 \right)^{\frac{1}{2}} = \sqrt{\textrm{tr}(A^{\ast}A)} = \sqrt{\textrm{tr}(AA^{\ast})}.
\end{align}$$

"""

# ╔═╡ 11135275-4514-4bf4-9b46-f6c36a5b98f0
md"""
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
"""

# ╔═╡ ab8fcb4c-d2c5-4334-838c-d65a60fdf7ec
md"""
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
"""

# ╔═╡ 6d8c90ca-14b2-4c55-a3e9-d621d420e816
md"""
_Example:_ Induced matrix norms describe the behavior of a matrix as an operator between its domain and range spaces. Let $\left\|\cdot\right\|_{(n)}$
and $\left\|\cdot \right\|_{(m)}$ be vector norms on the domain and range of $A \in \mathbb{C}^{m\times n}$, respectively. 

Then the norm $\left\|A\right\|_{(m, n)}$ induced by the vector norms $\left\|\cdot\right\|_{(n)}$
and $\left\|\cdot\right\|_{(m)}$, is the smallest number $C$ satisfying
$$\begin{align}
\left\|Ax\right\|_{(m)} \leq C \left\|x\right\|_{(n)}\end{align}$$

for all vector $x\in \mathbb{C}^{n}$. In other words,

$$\begin{align}
\left\|A\right\|_{(m, n)} = \sup_{x\in \mathbb{C}^{n}, x\neq 0} \frac{\left\|Ax\right\|_{(m)}}{\left\|x\right\|_{(n)}} = \sup_{x\in \mathbb{C}^{n}, \left\|x\right\|_{(n)}=1} \left\|Ax\right\|_{(m)}.
\end{align}$$

Supremum means "least upper bound". Under this definition, we can show the 1-norm and $\infty$-norm of a matrix are respectively the "maximum column sum" and the "maximum row sum":

$$\begin{align}
\left\|A\right\|_{1} = \max_{1\leq j \leq n} \left\|a_j\right\|_1, \\
\left\|A\right\|_{\infty} = \max_{1\leq i \leq m} \left\|a_i^{\ast}\right\|_{1}.
\end{align}$$
"""

# ╔═╡ d26a3567-015a-4793-8fb2-fc251f3a1d3f
md"""
* For any $A\in\mathbb{C}^{m\times n}$ and unitary $Q\in \mathbb{C}^{m\times m}$, the 2-norm and Frobenius norms are invariant under unitary multiplication: $\left\|QA\right\|_{2} = \left\|A\right\|_{2}$ and $\left\|QA\right\|_{F} = \left\|A\right\|_{2}$.
"""

# ╔═╡ 7f0c9a43-cb5c-4844-935b-526fa4f58342
md"""
**Note**: For more about unitarily invariant norms and links to the singular value decomposition, this post by Nick Higham is a good starting point! [Nick Higham's blog](https://nhigham.com/2021/02/02/what-is-a-unitarily-invariant-norm/)
"""

# ╔═╡ 7847aa02-3f6b-4bd1-bb96-93e67eac889d
md"""
## Visualizing the norms

One way to compare different norms is to examine the locus of points with equal distances for the different norms! We will do so below for 1 norm and 2 norm, and we can change the number of samples to make the shape more apparent:

"""

# ╔═╡ e5f59ca2-cd27-4bc7-b0a5-cf073445a658
@bind nSamples Slider(6:500, default=200, show_value=true)

# ╔═╡ 51241efc-e41e-4958-bfc3-dc1558e45dd3
samples = randn(nSamples, 2)

# ╔═╡ bb304dad-6d6f-4628-a21a-ffdabeab5909
normalized_n2_samples = samples ./ [norm(samples[i, :]) for i in 1:nSamples] # geneerate points at same l2 distance from the origin

# ╔═╡ ea1d9a5f-e314-4ea4-9373-6e833305af4c
normalized_n1_samples = samples ./ [norm(samples[i, :], 1) for i in 1:nSamples]

# ╔═╡ 79739c96-3de6-433d-b593-74e74471c6f6
md"""
It is also fun to look at other norms this way. Let us do so for `p = 1.5` and `p = Inf`:
"""

# ╔═╡ 6f8fd801-ccb8-4fb2-ba15-3faaa64221d5
PlutoUI.Resource("https://i.imgur.com/dqB1yRs.jpg")

# ╔═╡ 90e70e76-daaa-4525-98ac-00118e965ea5
md"""
The idea of using different norms, as well as the significance of these shapes will become apparent when we dive into linear regression problems in subsequent chapters.😃◈◉
"""

# ╔═╡ ed417cf7-db23-4c9e-b542-fbaae9cf8194
md"""
# Singular Value Decomposition (SVD)

"""

# ╔═╡ 3be12500-601d-4426-a4ec-8ae7e845988e
md"""
!!! quote "Quote by Lloyd Trefethen"
    *Many problems of linear algebra can be better understood if we first ask the question: what if we take the SVD?*
"""

# ╔═╡ 32dda940-dea1-4c78-8c95-c1507bcac801
md"""
The SVD used to be regarded as an obscure theoretical concept in linear algebra right until the 1960s, when Gene Golub came along and published the first practical algorithm to compute it, and in fact a modified version of that is still in use today in several programming environments ! Today, it is a workhorse of so many different methods, and any discussion of matrix factorization would be incomplete without diving into the SVD, and its closely related cousin, the PCA. 

**Fun Fact**: Here's the actual license plate of Gene Golub, bearing his most famous contribution. 
"""

# ╔═╡ 7158f8c8-d46c-413b-981e-6a3139d2e9c2
PlutoUI.Resource("https://i.imgur.com/aU7CIoe.jpeg")

# ╔═╡ 06aadb33-d629-411c-93e9-cfb5b51040f1
md"""
* Given $A\in \mathbb{C}^{m\times n}$, a singular value decomposition (SVD) of $A$ is a factorization
$$\begin{align}
    A = U\Sigma V^{\ast},
\end{align}$$

where $U\in \mathbb{C}^{m\times m}$ is unitary, $V\in \mathbb{C}^{n\times n}$ is unitary, and $\Sigma \in \mathbb{R}^{m\times n}$ is diagonal. The diagonal entries $\sigma_j$ of $\Sigma$ are nonnegative and in nonincreasing order $\sigma_1 \geq \sigma_2 \geq \ldots \geq \sigma_p \,\,\geq\,\, 0$, where $p=\min(m,n)$. 
"""

# ╔═╡ 5db08a2b-557c-4c00-a3b6-c27d282e9117
md"""
* Every matrix $A\in \mathbb{C}^{m\times n}$ has a SVD. The singular values $\{\sigma_j\}$ are uniquely determined.
"""

# ╔═╡ 4bde7f51-8d8f-4dfc-8735-5838449aac36
md"""
* Geometry interpretation: the image of a unit sphere under any $m\times n$ matrix is a hyperellipse. That is, let $S$ be a unit sphere in $\mathbb{R}^n$, and given any $A\in \mathbb{R}^{m \times n}$, then $AS$ is a hyperellipse in $\mathbb{R}^m$.  

$\{\sigma_i u_i\}_{i=1}^{m}$ denote the _principal semiaxes_ of the hyperellipse, with magnitude $\sigma_i$ along orthogonal (unit) directions $u_i$. 
If $m\geq n$ for example, then $n$ of the $\sigma_i$ will be non-zero. 
"""

# ╔═╡ f3d40a46-ff76-465b-8516-41dfd52e1525
md"""
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
"""

# ╔═╡ 3b15e60d-a86e-481e-a311-190837820d6d
md"""
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
"""

# ╔═╡ 25905881-dac2-44fb-a923-fcf38f5d32fe
md"""
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
"""

# ╔═╡ b25f85c0-1a52-46a1-a1c4-341bef17166f
md"""
A _reduced SVD_ for $m\times n$ matrix $A$ when $m\geq n$ is a factorization
$$\begin{align}
A = \hat{U}\hat{\Sigma}V^{\ast},
\end{align}$$

where $\hat{\Sigma}$ is an $n\times n$ diagonal matrix with positive real entries (if $A$ full rank), $\hat{U}$ is $m\times n$ with orthonormal columns, and $V$ is $n\times n$ with orthonormal columns. $V$ is therefore unitary; $\hat{U}$ is not unitary (it is not square).
"""

# ╔═╡ 1e8fde01-11a4-4a93-81aa-19b6104dcfc0
md"""
* The rank of $A$ is $r$, the number of nonzero singular values. 
"""

# ╔═╡ 89e4f2b1-978e-4483-ab66-513cda8d398c
md"""
*  $\left\|A\right\|_{2}=\sigma_1$, and $\left\|A\right\|_{F}=\sqrt{\sigma_1^2+\sigma_2^2+\ldots+\sigma_r^2}$.
"""

# ╔═╡ b5782f0d-653e-40c0-a810-353707dc4063
A

# ╔═╡ e8565946-2bba-4242-b991-21edb8686a54
opnorm(A)

# ╔═╡ 66f74b82-c99e-437e-b333-26660a51b67e
U, s, V = svd(A)

# ╔═╡ 80e205fd-5e75-435e-8811-8fbbdb6af24d
s[1]

# ╔═╡ c85dafce-2923-4998-8501-9a96687761de
isapprox(opnorm(A), s[1])

# ╔═╡ 693772f0-a030-40df-9707-139a06ab7cc8
# Exercise: Try writing out the factorization for B = [0 -2
# 													   3 0]

# ╔═╡ 0896d8b7-ec44-4326-aed4-e63d2a2dee48
# Try: Write out the SVD for B. Would the singular values change from A? What about the singular vectors?
B = [0 -2   
	 3 0] 

# ╔═╡ f8f86fc8-0933-4ea8-941d-e80a5621c275
md"""
* A square matrix $A\in \mathbb{C}^{m\times m}$ has the _eigenvalue decomposition_
$$\begin{align}
    A = X\Lambda X^{-1},
\end{align}$$

where the columns of $X\in \mathbb{C}^{m\times m}$ contain eigenvectors of $A$,
and $\Lambda$ is a diagonal matrix whose entries $\lambda_j$ are the eigenvalues of $A$. 

Individually, we can write the more familiar form of $Ax_j=\lambda_j x_j$. Lecture 24 of Trefethen & Bau (1997) has details on eigenvalue problems.
"""

# ╔═╡ b8a46009-b47a-403d-a263-ca510515b4be
md"""
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
"""

# ╔═╡ b7c937c2-e872-44fb-8eb0-97b259b587e3

md"""
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
"""

# ╔═╡ f34fb827-c834-4397-88f7-a8fe11f7b7e6
md"""
* Nonzero singular values of $A$ are the square roots of nonzero eigenvalues of $A^{\ast}A$ or $AA^{\ast}$. 
"""

# ╔═╡ 9b6b6f7b-166c-426c-b564-42e5ab359638
md"""
* **very cool property!** If $A=A^{\ast}$, then singular values of $A$ are the absolute values of the eigenvalues of $A$.
"""

# ╔═╡ 7bd5aea7-793f-4bde-92e7-5e55cba8f4b1
md"""
* For $A\in\mathbb{C}^{m\times m}$, $|\det(A)|=\prod_{i=1}^{m} \sigma_i$.
"""

# ╔═╡ fb3be364-a167-44da-bf3f-30b35564c50f


# ╔═╡ 364fa99f-7e4a-4f0b-97c3-1537b7006088
md"""
## Low rank approximation 

* It is easy to see from the SVD that $A$ is the sum of $r$ rank-one matrices:
$$\begin{align}
    A = \sum_{j=1}^{r} \sigma_j u_j v_j^{\ast}.
\end{align}$$
"""

# ╔═╡ fa6baeb0-2614-41ae-a59c-02ad86f862f0
md"""
We can thus truncate this expansion to arrive at a $\nu$-term low-rank approximation ($0\leq \nu \leq r$) 
(also known as a truncated SVD):

$$\begin{align}
    A_{\nu} = \sum_{j=1}^{\nu} \sigma_j u_j v_j^{\ast}
\end{align}$$

then we can show $\left\|A-A_{\nu}\right\|_2 = \sigma_{\nu+1}$, and also $\left\|A-A_{\nu}\right\|_{F} = \sqrt{\sigma_{\nu+1}^2+\ldots+\sigma_{r}^2}$ (if $\nu=\min(m,n)$ already, then define $\sigma_{\nu+1}=0$). In other words, the $\nu$-term partial sum captures as much of the energy of $A$ as possible, since it retains the $\nu$ \emph{largest} singular values. In this sense, this truncated approximation is optimal.
"""

# ╔═╡ 2a5254b2-4613-4cbd-bde8-d966eab0684c
md"""
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
"""

# ╔═╡ a4fdc854-28a2-44e7-b170-586159e062c2
md"""
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
$lbreak
"""

# ╔═╡ d473ee59-9526-43ad-b083-42204e9e41b3
#= md"""
## Food for thought

A small part of our notebooks where we can reflect on what we learnt, and see what other interesting questions we can ask about it!

""" =#

# ╔═╡ ea85bc4e-72c4-46b4-8852-552734c4dc8e
md"""

**Delete Cell Afterwards**: All of this stuff for next notebook and / or include interesting problems based on this stuff as assignments???

#### Going beyond the classics: Example of Randomized Methods - using sketching for solving a large linear systems


#### Condition numbers for linear system and for eigenvalue problem


#### Higher Order SVD - generalization of matrix SVD


#### Other transforms formulated as matrix problems - for example, using Fourier, wavelet basis for image compression - superior methods over SVD!

"""

# ╔═╡ c1143820-ca59-4cc2-bef6-52198a669962


# ╔═╡ 487a7420-c4dc-403b-a89a-f264335aa67c


# ╔═╡ 3d424b6b-5c82-444b-93a0-e0db512c0488
md"""
# Appendix
Miscellaneous code required to keep the notebook running and polish it up!
"""

# ╔═╡ 9072fcb2-2fcf-4dc1-b5b7-96dc35a56f7d
# add table of contents
	PlutoUI.TableOfContents(title="Lecture 02: Linear Algebra Review 1", indent=true, depth=3)
	

# ╔═╡ ab28b878-5083-4ffc-935c-e00e23a9e272
logocolors = Colors.JULIA_LOGO_COLORS

# ╔═╡ c3c9f4b1-6b88-409f-9208-703e56520fe3
begin
		# Trying the polygon method for adding vectors - albeit super crude!
		
		quiver([0], [0],
			   quiver=([x₁ * a₁[1]], [x₁ * a₁[2]]),
			   line=(logocolors.red, 2.5, :dash),
			   label="x₁a₁")
	
		quiver!([x₁ * a₁[1]], [x₁ * a₁[2]],
				quiver=([x₂ * a₂[1]], [x₂ * a₂[2]]),
				line=(logocolors.green, 2.5, :dash),
				label="x₂a₂")		

		quiver!([x₂ * a₂[1] + x₁ * a₁[1]], [x₂ * a₂[2] + x₁ * a₁[2]],
		   		quiver=([x₃ * a₃[1]], [x₃ * a₃[2]]),
		   		line=(logocolors.purple, 2.5, :dash),
				label="x₃a₃")
		
		quiver!([x₃ * a₃[1] + x₂ * a₂[1] + x₁ * a₁[1]], [x₃ * a₃[2] + x₂ * a₂[2] + x₁ * a₁[2]],
		   		quiver=([x₄ * a₄[1]], [x₄ * a₄[2]]),
		   		line=(logocolors.blue, 2.5, :dash),
				label="x₄a₄")
	
		quiver!([0], [0],
		   		quiver=([b[1]], [b[2]]),
		   		line=(:black, 3),
				label="b (Resultant)")
		
		annotate!(x₁ * a₁[1], b[2], text("x₁a₁", logocolors.red, 14, :center, :bold))
		annotate!(x₁ * a₁[1], b[2] - 2, text("x₂a₂", logocolors.green, 14, :center, :bold))
		annotate!(x₁ * a₁[1], b[2] - 4, text("x₃a₃", logocolors.purple, 14, :center, :bold))
		annotate!(x₁ * a₁[1], b[2] - 6, text("x₄a₄", logocolors.blue, 14, :center, :bold))
		annotate!(b[1], b[2] + 1, text("b", :black, :center, 14, :bold))

		plot!(title = "x1 = $(x₁), x2 = $(x₂), x3 = $(x₃), x4 = $(x₄)")
		plot!(aspect_ratio=:equal,
			  size = (600, 600))
		plot!(legend=:topright) # For some reason, this won't work for quiver with gr()!! Using text annotations instead
		plot!(framestyle=:origin)
end

# ╔═╡ 973f5717-d2b3-46f2-b7b3-01f07311efae
begin	
		quiver([0], [0],
			   quiver=([btil[1]], [btil[2]]),
			   line=(logocolors.red, 2.5, :dash),
			   label="Ax")
		quiver!([0], [0],
			   quiver=([bb[1]], [bb[2]]),
			   line=(:black, 3),
			   label="b")
		plot!(framestyle=:origin)
end

# ╔═╡ 5c9195c1-02f1-4226-924a-dfb55aa02fe1
begin
	p = plot()
	for i in 1:nSamples
		scatter!(p, [normalized_n2_samples[i, 1]], [normalized_n2_samples[i, 2]], 
				 marker=(:circle, logocolors.red), label="")
	end
	
	
	for i in 1:nSamples
		scatter!(p, [normalized_n1_samples[i, 1]], 
				 [normalized_n1_samples[i, 2]],
				 marker=(:circle, logocolors.purple), label="")
	end
	
	plot!(framestyle=:origin, aspect_ratio=:equal)
end

# ╔═╡ 1e91a35f-f6c8-4e9e-a065-2e568bcf73ac
# Download flags
begin
	flag_denmark = mktemp() do fn,f
		download("https://cdn.britannica.com/07/8007-004-8CF0B1A9/Flag-Denmark.jpg", fn)
		load(fn)
	end

	flag_colombia = mktemp() do fn,f
		download("https://cdn.britannica.com/68/7668-050-9304EBB7/Flag-Colombia.jpg", fn)
		load(fn)
	end
	
	flag_benin = mktemp() do fn,f
		download("https://cdn.britannica.com/37/5037-050-D5FCC732/Flag-Benin.jpg", fn)
		load(fn)
	end
	
	md""" We are done downloading the flags!"""
	
end

# ╔═╡ 3615c495-69a7-4bff-9da9-22bfb6a98f64
colombia = Gray.(flag_colombia)

# ╔═╡ 9b8a847f-6ebd-4dbc-8aa5-f7df36355e87
colombia2 = convert(Array{Float64}, colombia)

# ╔═╡ 49271f41-0633-4faf-89f1-b494e66fbcfc
size(colombia)

# ╔═╡ c308d131-e5dd-4993-988d-aa7248cc59d2
rank(colombia)

# ╔═╡ d7cfd955-49cd-406b-9527-5b8525f6e02e
benin = Gray.(flag_benin)

# ╔═╡ 706efb94-e834-4a3c-bb79-b4dd85931853
size(benin)

# ╔═╡ f8c6ebd3-528e-484f-951a-e448399a0869
rank(benin)

# ╔═╡ 45781360-ae14-4a1a-8ee9-ccdeb2e2347f
# begin
# 		# Vectors in their original positions! We shifted them to show the polygon method of adding vectors.
		
# 		quiver([0], [0],
# 			   quiver=([x₁ * a₁[1]], [x₁ * a₁[2]]),
# 			   line=(logocolors.red, 2))
		
		
# 		quiver!([0], [0],
# 			   quiver=([x₂ * a₂[1]], [x₂ * a₂[2]]),
# 			   line=(logocolors.green, 2))
			
# 		quiver!([0], [0],
# 		   		quiver=([x₃ * a₃[1]], [x₃ * a₃[2]]),
# 		   		line=(logocolors.purple, 2))
		
		
# 		quiver!([0], [0],
# 		   		quiver=([x₄ * a₄[1]], [x₄ * a₄[2]]),
# 		   		line=(logocolors.blue, 2))
	
	
# 		quiver!([0], [0],
# 		   		quiver=([b[1]], [b[2]]),
# 		   		line=(:black, 3))
		
# 		plot!(title = "x1 = $(x₁), x2 = $(x₂), x3 = $(x₃), x4 = $(x₄)")
# 		plot!(aspect_ratio=:equal,
# 			  size = (500, 500))
# end

# ╔═╡ Cell order:
# ╟─a876703f-6e85-488b-adb2-c566073f304d
# ╠═ae321ee5-c21f-4d5a-b1b7-4b2f0e4d8cd9
# ╟─6f622d7e-0684-11ec-0e37-dbc7eccb7a73
# ╟─80fd4f57-dfd8-4e1d-a5c8-67bf6cbfb281
# ╟─afadb596-8200-47c1-933d-b4d60dced054
# ╠═5fa5bbfb-27ee-4621-b66b-7f05ab37fb2d
# ╟─d5b726b0-7aef-4c35-8fcf-55be19ec8ffa
# ╠═d347889f-15bf-4459-8411-cfc00058afc3
# ╠═1f40c451-c984-46da-9a35-337b97e254ce
# ╠═bcf0eee2-0a56-43c5-836e-80a29a552e01
# ╠═40e0d8b2-6f34-4dcd-bf3c-6e2ca30422f5
# ╟─b7a25f65-92cc-446f-8f76-5d1d1774a83d
# ╠═e67dffdf-cc8f-4ee3-93cc-a863a61bc3f3
# ╟─74cad985-f418-41b5-8f78-84e7f5475838
# ╠═2026b715-c4c9-4b91-82da-96eed87f3d39
# ╠═58ff8911-289a-4d31-ad87-c6e94525dbb0
# ╟─d9d3e037-dc13-45e5-a97f-408122d512a0
# ╟─e2e120fe-778e-44e1-8337-37c9e57d1951
# ╠═8ef28a47-417f-45a8-acc7-76547a63da9a
# ╟─b975fb26-094c-4587-9b03-1eef425b4cff
# ╟─894ed0d4-941e-46a2-8b9a-357bbd30fde9
# ╠═f98ba8f4-2ce6-47d8-9231-d55b178385a7
# ╠═12e6e9e0-629d-4fef-abc5-1ce917f5cf6d
# ╟─9e15631a-a0bb-4e1a-b622-1f4d2d8b789a
# ╟─c3c9f4b1-6b88-409f-9208-703e56520fe3
# ╟─a36c5035-194c-4c98-b120-5481ac10fc25
# ╟─ce3c3a23-7edd-4697-918d-222a78291538
# ╟─d69540d9-2647-48f3-9b80-3a6f72f19035
# ╟─8c6db3b1-a11a-42ff-809d-21e2dff6b1ba
# ╟─e29cd04a-c020-4c7a-8fda-39a0126fb83f
# ╠═802b8714-3f01-4d80-a6ed-b28d432839c4
# ╠═1903e8b6-f8d3-4ed5-93ba-4aaca44608c9
# ╟─4e3801f4-cb8f-4d45-8e8c-d2556bdbcbc8
# ╟─514b975e-0fb3-4943-8596-37dc99be3c10
# ╟─a48392f9-2eb2-46f2-af6a-b032816390c3
# ╠═16f0cd08-6cb9-446f-876d-977430fd6705
# ╠═1b439146-1fe0-4756-ae3d-d23982aa3e01
# ╠═a89bc125-a20c-40d6-87c3-4b9ba6754665
# ╟─77c7cc1e-cd25-4dcd-86e5-8c615bde9017
# ╟─973f5717-d2b3-46f2-b7b3-01f07311efae
# ╟─920ed5b7-cda9-49fc-bac6-91a2cf36f79d
# ╟─25ae9da1-c85e-490b-a0da-693cf9cd5499
# ╠═0f370f2c-0371-40d0-a4d7-6979c81dfbed
# ╟─3e7052ea-743b-4c02-b266-b90b402d87c2
# ╟─f5ec28f0-b9ac-478a-adda-bf641f868efa
# ╟─b2466928-15a6-4626-9acb-58f1b3ff7e74
# ╟─0220b07d-abc1-4c38-acd7-eaaa3148bd5d
# ╟─8c43733c-2786-4d31-94f4-0d190fcd3dca
# ╟─36e25a97-78b1-437d-ab98-295166daba5b
# ╟─02b49d9c-098b-45ec-8837-347ec3e2e9f0
# ╟─0eafbb01-c078-4ea0-914f-0412ab6c890a
# ╟─fc2770f0-425b-4b07-8315-070f43ea30cf
# ╠═3615c495-69a7-4bff-9da9-22bfb6a98f64
# ╠═9b8a847f-6ebd-4dbc-8aa5-f7df36355e87
# ╠═49271f41-0633-4faf-89f1-b494e66fbcfc
# ╠═c308d131-e5dd-4993-988d-aa7248cc59d2
# ╠═d7cfd955-49cd-406b-9527-5b8525f6e02e
# ╠═706efb94-e834-4a3c-bb79-b4dd85931853
# ╠═f8c6ebd3-528e-484f-951a-e448399a0869
# ╠═9cd171cd-82b8-4208-b6fa-94d298a6e173
# ╟─e304923f-9ab2-4a43-abd5-95ebb84d023c
# ╟─b0832388-5c01-40c0-ad79-f368b646abc9
# ╟─31ab7c4f-335d-46a4-9688-cbda0feec3ec
# ╠═98b1b9dc-0f37-4f3e-91bf-30be8e957f34
# ╟─2a9d5631-8176-44f0-beb1-3fda3058d56f
# ╠═b6c98f01-51b0-47e3-889d-d2b0eb7e767f
# ╠═ddee31eb-56a4-4240-9c90-a484d2dea11f
# ╠═810e9ddc-41a0-4c90-80bf-84e978ba1bb9
# ╠═3a77c857-8a57-49ab-bcf5-79c1dece416e
# ╠═d811dd3d-c419-491e-97f4-35f7159f6a6a
# ╠═0c1814e3-4b62-4434-8868-598eb7127836
# ╠═e6588752-18c3-4503-a427-f0b4b8246c0d
# ╟─c41d4e8e-58b8-4eff-95ce-88aa50951059
# ╟─e9ab39f4-21aa-4929-8608-3a1c577136e9
# ╟─0b20a6b4-7a38-4f3c-9c4a-e2268e417898
# ╠═23cf0c46-1b15-4659-a6fb-17c036bb5e5c
# ╠═c47ea5f8-2282-4fa6-b387-c4116cff0b85
# ╠═8d3635ea-f227-41e1-b524-9a80213f2fdf
# ╠═441b7b2f-0b4a-409f-9672-b5bf103eef44
# ╟─2b37dac7-7dbf-40c8-826d-cde1800ac896
# ╟─890bfe38-bf37-45f9-a45f-3a665152747a
# ╟─3913c8f5-e3a4-404d-989a-024291f3912f
# ╟─e54eb17d-d4d3-4d35-bebf-b0a488b9b56a
# ╟─fbb084c7-6daf-46f8-851b-cc6bd7dc236f
# ╟─d6d3ce82-4211-483d-8353-0fd3a7756321
# ╟─fa8fed50-c898-44f6-aa38-fb211e89e62a
# ╟─f66c61bf-eca1-47d1-b8f3-928459daab08
# ╟─a0f85f74-a6d8-4efd-b295-3437be3dbf78
# ╟─9ddeb1b1-7436-4d02-b772-9ed1ee5c513e
# ╟─97c03b72-f8a5-40e6-841a-d0ebd4a40c30
# ╟─52700936-6738-4731-b643-335f73fa2111
# ╟─eeedaee2-7f68-4431-aaab-f3f0962455f9
# ╟─dced6ea7-55e8-43f6-8652-69bdabecfa6b
# ╟─11135275-4514-4bf4-9b46-f6c36a5b98f0
# ╟─ab8fcb4c-d2c5-4334-838c-d65a60fdf7ec
# ╟─6d8c90ca-14b2-4c55-a3e9-d621d420e816
# ╟─d26a3567-015a-4793-8fb2-fc251f3a1d3f
# ╟─7f0c9a43-cb5c-4844-935b-526fa4f58342
# ╟─7847aa02-3f6b-4bd1-bb96-93e67eac889d
# ╟─e5f59ca2-cd27-4bc7-b0a5-cf073445a658
# ╠═51241efc-e41e-4958-bfc3-dc1558e45dd3
# ╠═bb304dad-6d6f-4628-a21a-ffdabeab5909
# ╟─5c9195c1-02f1-4226-924a-dfb55aa02fe1
# ╠═ea1d9a5f-e314-4ea4-9373-6e833305af4c
# ╟─79739c96-3de6-433d-b593-74e74471c6f6
# ╟─6f8fd801-ccb8-4fb2-ba15-3faaa64221d5
# ╟─90e70e76-daaa-4525-98ac-00118e965ea5
# ╟─ed417cf7-db23-4c9e-b542-fbaae9cf8194
# ╟─3be12500-601d-4426-a4ec-8ae7e845988e
# ╟─32dda940-dea1-4c78-8c95-c1507bcac801
# ╟─7158f8c8-d46c-413b-981e-6a3139d2e9c2
# ╟─06aadb33-d629-411c-93e9-cfb5b51040f1
# ╟─5db08a2b-557c-4c00-a3b6-c27d282e9117
# ╟─4bde7f51-8d8f-4dfc-8735-5838449aac36
# ╟─f3d40a46-ff76-465b-8516-41dfd52e1525
# ╟─3b15e60d-a86e-481e-a311-190837820d6d
# ╟─25905881-dac2-44fb-a923-fcf38f5d32fe
# ╟─b25f85c0-1a52-46a1-a1c4-341bef17166f
# ╟─1e8fde01-11a4-4a93-81aa-19b6104dcfc0
# ╟─89e4f2b1-978e-4483-ab66-513cda8d398c
# ╠═b5782f0d-653e-40c0-a810-353707dc4063
# ╠═e8565946-2bba-4242-b991-21edb8686a54
# ╠═66f74b82-c99e-437e-b333-26660a51b67e
# ╠═80e205fd-5e75-435e-8811-8fbbdb6af24d
# ╠═c85dafce-2923-4998-8501-9a96687761de
# ╠═693772f0-a030-40df-9707-139a06ab7cc8
# ╠═0896d8b7-ec44-4326-aed4-e63d2a2dee48
# ╟─f8f86fc8-0933-4ea8-941d-e80a5621c275
# ╟─b8a46009-b47a-403d-a263-ca510515b4be
# ╟─b7c937c2-e872-44fb-8eb0-97b259b587e3
# ╟─f34fb827-c834-4397-88f7-a8fe11f7b7e6
# ╟─9b6b6f7b-166c-426c-b564-42e5ab359638
# ╟─7bd5aea7-793f-4bde-92e7-5e55cba8f4b1
# ╠═fb3be364-a167-44da-bf3f-30b35564c50f
# ╟─364fa99f-7e4a-4f0b-97c3-1537b7006088
# ╟─fa6baeb0-2614-41ae-a59c-02ad86f862f0
# ╟─2a5254b2-4613-4cbd-bde8-d966eab0684c
# ╟─a4fdc854-28a2-44e7-b170-586159e062c2
# ╠═d473ee59-9526-43ad-b083-42204e9e41b3
# ╟─ea85bc4e-72c4-46b4-8852-552734c4dc8e
# ╠═c1143820-ca59-4cc2-bef6-52198a669962
# ╠═487a7420-c4dc-403b-a89a-f264335aa67c
# ╟─3d424b6b-5c82-444b-93a0-e0db512c0488
# ╠═51a640f1-7007-44ac-9751-655083073fab
# ╠═8a09a85c-cfcd-4478-ae21-7852227b3f61
# ╠═e4277174-73eb-4925-a104-21fa937ae3b6
# ╠═9072fcb2-2fcf-4dc1-b5b7-96dc35a56f7d
# ╠═ab28b878-5083-4ffc-935c-e00e23a9e272
# ╠═1e91a35f-f6c8-4e9e-a065-2e568bcf73ac
# ╠═45781360-ae14-4a1a-8ee9-ccdeb2e2347f
