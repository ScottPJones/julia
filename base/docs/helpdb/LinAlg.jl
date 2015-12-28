# This file is a part of Julia. License is MIT: http://julialang.org/license

"""
    svdvals(A)

Returns the singular values of `A`.
"""
svdvals(A)

"""
    svdvals(A, B)

Return only the singular values from the generalized singular value decomposition of `A` and `B`.
"""
svdvals(A, B)

"""
    svd(A, [thin=true]) -> U, S, V

Wrapper around `svdfact` extracting all parts the factorization to a tuple. Direct use of
`svdfact` is therefore generally more efficient. Computes the SVD of `A`, returning `U`,
vector `S`, and `V` such that `A == U*diagm(S)*V'`. If `thin` is `true`, an economy mode
decomposition is returned. The default is to produce a thin decomposition.
"""
svd

"""
    svd(A, B) -> U, V, Q, D1, D2, R0

Wrapper around `svdfact` extracting all parts the factorization to a tuple. Direct use of
`svdfact` is therefore generally more efficient. The function returns the generalized SVD of
`A` and `B`, returning `U`, `V`, `Q`, `D1`, `D2`, and `R0` such that `A = U*D1*R0*Q'` and `B =
V*D2*R0*Q'`.
    """
svd(A::AbstractMatrix, B::AbstractMatrix)

"""
    svdvals!(A)

Returns the singular values of `A`, while saving space by overwriting the input.
"""
svdvals!

"""
    svdfact(A, [thin=true]) -> SVD

Compute the Singular Value Decomposition (SVD) of `A` and return an `SVD` object. `U`, `S`,
`V` and `Vt` can be obtained from the factorization `F` with `F[:U]`, `F[:S]`, `F[:V]` and
`F[:Vt]`, such that `A = U*diagm(S)*Vt`. If `thin` is `true`, an economy mode decomposition
is returned. The algorithm produces `Vt` and hence `Vt` is more efficient to extract than
`V`. The default is to produce a thin decomposition.
"""
svdfact(A)

"""
    svdfact(A, B) -> GeneralizedSVD

Compute the generalized SVD of `A` and `B`, returning a `GeneralizedSVD` Factorization
object `F`, such that `A = F[:U]*F[:D1]*F[:R0]*F[:Q]'` and `B = F[:V]*F[:D2]*F[:R0]*F[:Q]'`.
"""
svdfact(A, B)

"""
    svdfact!(A, [thin=true]) -> SVD

`svdfact!` is the same as [`svdfact`](:func:`svdfact`), but saves space by overwriting the
input `A`, instead of creating a copy. If `thin` is `true`, an economy mode decomposition is
returned. The default is to produce a thin decomposition.
"""
svdfact!

"""
    eigfact!(A, [B])

Same as [`eigfact`](:func:`eigfact`), but saves space by overwriting the input `A` (and
`B`), instead of creating a copy.
"""
eigfact!

"""
    eigmax(A)

Returns the largest eigenvalue of `A`.
"""
eigmax

"""
    triu!(M)

Upper triangle of a matrix, overwriting `M` in the process.
"""
triu!(M)

"""
    triu!(M, k)

Returns the upper triangle of `M` starting from the `k`th superdiagonal, overwriting `M` in
the process.
"""
triu!(M, k)

"""
    diag(M[, k])

The `k`th diagonal of a matrix, as a vector. Use `diagm` to construct a diagonal matrix.
"""
diag

"""
    istril(A) -> Bool

Test whether a matrix is lower triangular.
"""
istril

"""
    nullspace(M)

Basis for nullspace of `M`.
"""
nullspace

"""
    peakflops(n; parallel=false)

`peakflops` computes the peak flop rate of the computer by using double precision
[`gemm!`](:func:`Base.LinAlg.BLAS.gemm!`). By default, if no arguments are specified, it
multiplies a matrix of size `n x n`, where `n = 2000`. If the underlying BLAS is using
multiple threads, higher flop rates are realized. The number of BLAS threads can be set with
`blas_set_num_threads(n)`.

If the keyword argument `parallel` is set to `true`, `peakflops` is run in parallel on all
the worker processors. The flop rate of the entire parallel computer is returned. When
running in parallel, only 1 BLAS thread is used. The argument `n` still refers to the size
of the problem that is solved on each processor.
"""
peakflops

"""
    isposdef!(A) -> Bool

Test whether a matrix is positive definite, overwriting `A` in the processes.
"""
isposdef!

"""
    tril!(M)

Lower triangle of a matrix, overwriting `M` in the process.
"""
tril!(M)

"""
    tril!(M, k)

Returns the lower triangle of `M` starting from the `k`th superdiagonal, overwriting `M` in
the process.
"""
tril!(M, k)

"""
    lu(A) -> L, U, p

Compute the LU factorization of `A`, such that `A[p,:] = L*U`.
"""
lu

"""
    tril(M)

Lower triangle of a matrix.
"""
tril(M)

"""
    tril(M, k)

Returns the lower triangle of `M` starting from the `k`th superdiagonal.
"""
tril(M,k)

"""
    triu(M)

Upper triangle of a matrix.
"""
triu(M)

"""
    triu(M, k)

Returns the upper triangle of `M` starting from the `k`th superdiagonal.
"""
triu(M, k)

"""
    SymTridiagonal(d, du)

Construct a real symmetric tridiagonal matrix from the diagonal and upper diagonal,
respectively. The result is of type `SymTridiagonal` and provides efficient specialized
eigensolvers, but may be converted into a regular matrix with [`full`](:func:`full`).
"""
SymTridiagonal

"""
    rank(M)

Compute the rank of a matrix.
"""
rank

"""
    sylvester(A, B, C)

Computes the solution `X` to the Sylvester equation `AX + XB + C = 0`, where `A`, `B` and
`C` have compatible dimensions and `A` and `-B` have no eigenvalues with equal real part.
"""
sylvester

"""
    cross(x, y)
    ×(x,y)

Compute the cross product of two 3-vectors.
"""
cross

"""
    lufact(A [,pivot=Val{true}]) -> F

Compute the LU factorization of `A`. The return type of `F` depends on the type of `A`. In
most cases, if `A` is a subtype `S` of AbstractMatrix with an element type `T` supporting `+`, `-`, `*`
and `/` the return type is `LU{T,S{T}}`. If pivoting is chosen (default) the element type
should also support `abs` and `<`. When `A` is sparse and have element of type `Float32`,
`Float64`, `Complex{Float32}`, or `Complex{Float64}` the return type is `UmfpackLU`. Some
examples are shown in the table below.

| Type of input `A`                              | Type of output `F`     | Relationship between `F` and `A`             |
|:-----------------------------------------------|:-----------------------|:---------------------------------------------|
| [`Matrix`](:func:`Matrix`)                     | `LU`                   | `F[:L]*F[:U] == A[F[:p], :]`                 |
| [`Tridiagonal`](:func:`Tridiagonal`)           | `LU{T,Tridiagonal{T}}` | `F[:L]*F[:U] == A[F[:p], :]`                 |
| [`SparseMatrixCSC`](:func:`SparseMatrixCSC`)   | `UmfpackLU`            | `F[:L]*F[:U] == (F[:Rs] .* A)[F[:p], F[:q]]` |

The individual components of the factorization `F` can be accessed by indexing:

| Component | Description                         | `LU` | `LU{T,Tridiagonal{T}}` | `UmfpackLU` |
|:----------|:------------------------------------|:-----|:-----------------------|:------------|
| `F[:L]`   | `L` (lower triangular) part of `LU` | ✓    | ✓                      | ✓           |
| `F[:U]`   | `U` (upper triangular) part of `LU` | ✓    | ✓                      | ✓           |
| `F[:p]`   | (right) permutation `Vector`        | ✓    | ✓                      | ✓           |
| `F[:P]`   | (right) permutation `Matrix`        | ✓    | ✓                      |             |
| `F[:q]`   | left permutation `Vector`           |      |                        | ✓           |
| `F[:Rs]`  | `Vector` of scaling factors         |      |                        | ✓           |
| `F[:(:)]` | `(L,U,p,q,Rs)` components           |      |                        | ✓           |

| Supported function | `LU` | `LU{T,Tridiagonal{T}}` | `UmfpackLU` |
|:-------------------|:-----|:-----------------------|:------------|
| `/`                | ✓    |                        |             |
| `\\`               | ✓    | ✓                      | ✓           |
| `cond`             | ✓    |                        | ✓           |
| `det`              | ✓    | ✓                      | ✓           |
| `logdet`           | ✓    | ✓                      |             |
| `logabsdet`        | ✓    | ✓                      |             |
| `size`             | ✓    | ✓                      |             |
"""
lufact

"""
    bkfact!(A) -> BunchKaufman

`bkfact!` is the same as [`bkfact`](:func:`bkfact`), but saves space by overwriting the
input `A`, instead of creating a copy.
"""
bkfact!

"""
    hessfact(A)

Compute the Hessenberg decomposition of `A` and return a `Hessenberg` object. If `F` is the
factorization object, the unitary matrix can be accessed with `F[:Q]` and the Hessenberg
matrix with `F[:H]`. When `Q` is extracted, the resulting type is the `HessenbergQ` object,
and may be converted to a regular matrix with [`full`](:func:`full`).
"""
hessfact

"""
    eigmin(A)

Returns the smallest eigenvalue of `A`.
"""
eigmin

"""
    diagind(M[, k])

A `Range` giving the indices of the `k`th diagonal of the matrix `M`.
"""
diagind

"""
    A_mul_B!(Y, A, B) -> Y

Calculates the matrix-matrix or matrix-vector product ``A⋅B`` and stores the result in `Y`,
overwriting the existing value of `Y`. Note that `Y` must not be aliased with either `A` or
`B`.

```jldoctest
julia> A=[1.0 2.0; 3.0 4.0]; B=[1.0 1.0; 1.0 1.0]; Y = similar(B); A_mul_B!(Y, A, B);

julia> Y
2x2 Array{Float64,2}:
 3.0  3.0
 7.0  7.0
```
"""
A_mul_B!

"""
    Bidiagonal(dv, ev, isupper)

Constructs an upper (`isupper=true`) or lower (`isupper=false`) bidiagonal matrix using the
given diagonal (`dv`) and off-diagonal (`ev`) vectors.  The result is of type `Bidiagonal`
and provides efficient specialized linear solvers, but may be converted into a regular
matrix with [`full`](:func:`full`).
"""
Bidiagonal

"""
    lufact!(A) -> LU

`lufact!` is the same as [`lufact`](:func:`lufact`), but saves space by overwriting the
input `A`, instead of creating a copy.  For sparse `A` the `nzval` field is not overwritten
but the index fields, `colptr` and `rowval` are decremented in place, converting from
1-based indices to 0-based indices.
"""
lufact!

"""
    isposdef(A) -> Bool

Test whether a matrix is positive definite.
"""
isposdef

"""
    pinv(M[, tol])

Computes the Moore-Penrose pseudoinverse.

For matrices `M` with floating point elements, it is convenient to compute
the pseudoinverse by inverting only singular values above a given threshold,
`tol`.

The optimal choice of `tol` varies both with the value of `M` and the intended application
of the pseudoinverse. The default value of `tol` is
`eps(real(float(one(eltype(M)))))*maximum(size(A))`, which is essentially machine epsilon
for the real part of a matrix element multiplied by the larger matrix dimension. For
inverting dense ill-conditioned matrices in a least-squares sense,
`tol = sqrt(eps(real(float(one(eltype(M))))))` is recommended.

For more information, see [^issue8859], [^B96], [^S84], [^KY88].

[^issue8859]: Issue 8859, "Fix least squares", https://github.com/JuliaLang/julia/pull/8859

[^B96]: Åke Björck, "Numerical Methods for Least Squares Problems",  SIAM Press, Philadelphia, 1996, "Other Titles in Applied Mathematics", Vol. 51. [doi:10.1137/1.9781611971484](http://epubs.siam.org/doi/book/10.1137/1.9781611971484)

[^S84]: G. W. Stewart, "Rank Degeneracy", SIAM Journal on Scientific and Statistical Computing, 5(2), 1984, 403-413. [doi:10.1137/0905030](http://epubs.siam.org/doi/abs/10.1137/0905030)

[^KY88]: Konstantinos Konstantinides and Kung Yao, "Statistical analysis of effective singular values in matrix rank determination", IEEE Transactions on Acoustics, Speech and Signal Processing, 36(5), 1988, 757-763. [doi:10.1109/29.1585](http://dx.doi.org/10.1109/29.1585)
"""
pinv

"""
    isdiag(A) -> Bool

Test whether a matrix is diagonal.
"""
isdiag

"""
    ishermitian(A) -> Bool

Test whether a matrix is Hermitian.
"""
ishermitian

"""
    issymmetric(A) -> Bool

Test whether a matrix is symmetric.
"""
issymmetric

"""
    svds(A; nsv=6, ritzvec=true, tol=0.0, maxiter=1000) -> (left_sv, s, right_sv, nconv, niter, nmult, resid)

`svds` computes largest singular values `s` of `A` using Lanczos or Arnoldi iterations. Uses
[`eigs`](:func:`eigs`) underneath.

Inputs are:

* `A`: Linear operator. It can either subtype of `AbstractArray` (e.g., sparse matrix) or
  duck typed. For duck typing `A` has to support `size(A)`, `eltype(A)`, `A * vector` and
  `A' * vector`.
* `nsv`: Number of singular values.
* `ritzvec`: Whether to return the left and right singular vectors `left_sv` and `right_sv`,
  default is `true`. If `false` the singular vectors are omitted from the output.
* `tol`: tolerance, see [`eigs`](:func:`eigs`).
* `maxiter`: Maximum number of iterations, see [`eigs`](:func:`eigs`).

**Example**

```julia
X = sprand(10, 5, 0.2)
svds(X, nsv = 2)
```
"""
svds

"""
    linreg(x, y) -> a, b

Perform linear regression. Returns `a` and `b` such that `a + b*x` is the closest straight
line to the given points `(x, y)`, i.e., such that the squared error between `y` and `a +
b*x` is minimized.

**Example**:

    using PyPlot
    x = [1.0:12.0;]
    y = [5.5, 6.3, 7.6, 8.8, 10.9, 11.79, 13.48, 15.02, 17.77, 20.81, 22.0, 22.99]
    a, b = linreg(x, y)          # Linear regression
    plot(x, y, "o")              # Plot (x, y) points
    plot(x, [a+b*i for i in x])  # Plot line determined by linear regression
"""
linreg(x,y)

"""
    linreg(x, y, w)

Weighted least-squares linear regression.
"""
linreg(x,y,w)

"""
    sqrtm(A)

If `A` has no negative real eigenvalues, compute the principal matrix square root of `A`,
that is the unique matrix ``X`` with eigenvalues having positive real part such that
``X^2 = A``. Otherwise, a nonprincipal square root is returned.

If `A` is symmetric or Hermitian, its eigendecomposition ([`eigfact`](:func:`eigfact`)) is
used to compute the square root. Otherwise, the square root is determined by means of the
Björck-Hammarling method, which computes the complex Schur form ([`schur`](:func:`schur`))
and then the complex square root of the triangular factor.

[^BH83]: Åke Björck and Sven Hammarling, "A Schur method for the square root of a matrix", Linear Algebra and its Applications, 52-53, 1983, 127-140. [doi:10.1016/0024-3795(83)80010-X](http://dx.doi.org/10.1016/0024-3795(83)80010-X)

"""
sqrtm

"""
    expm(A)

Compute the matrix exponential of `A`, defined by

```math
e^A = \\sum_{n=0}^{\\infty} \\frac{A^n}{n!}.
```

For symmetric or Hermitian `A`, an eigendecomposition ([`eigfact`](:func:`eigfact`)) is
used, otherwise the scaling and squaring algorithm (see [^H05]) is chosen.

[^H05]: Nicholas J. Higham, "The squaring and scaling method for the matrix exponential revisited", SIAM Journal on Matrix Analysis and Applications, 26(4), 2005, 1179-1193. [doi:10.1137/090768539](http://dx.doi.org/10.1137/090768539)

"""
expm

"""
    hessfact!(A)

`hessfact!` is the same as [`hessfact`](:func:`hessfact`), but saves space by overwriting
the input `A`, instead of creating a copy.
"""
hessfact!

"""
    diagm(v[, k])

Construct a diagonal matrix and place `v` on the `k`th diagonal.
"""
diagm

"""
    logabsdet(M)

Log of absolute value of determinant of real matrix. Equivalent to `(log(abs(det(M))), sign(det(M)))`,
but may provide increased accuracy and/or speed.
"""
logabsdet

"""
    norm(A, [p])

Compute the `p`-norm of a vector or the operator norm of a matrix `A`, defaulting to the `p=2`-norm.

For vectors, `p` can assume any numeric value (even though not all values produce a
mathematically valid vector norm). In particular, `norm(A, Inf)` returns the largest value
in `abs(A)`, whereas `norm(A, -Inf)` returns the smallest.

For matrices, the matrix norm induced by the vector `p`-norm is used, where valid values of
`p` are `1`, `2`, or `Inf`. (Note that for sparse matrices, `p=2` is currently not
implemented.) Use [`vecnorm`](:func:`vecnorm`) to compute the Frobenius norm.
"""
norm

"""
    lyap(A, C)

Computes the solution `X` to the continuous Lyapunov equation `AX + XA' + C = 0`, where no
eigenvalue of `A` has a zero real part and no two eigenvalues are negative complex
conjugates of each other.
"""
lyap

"""
    condskeel(M, [x, p])

```math
\\kappa_S(M, p) & = \\left\\Vert \\left\\vert M \\right\\vert \\left\\vert M^{-1} \\right\\vert  \\right\\Vert_p \\\\
\\kappa_S(M, x, p) & = \\left\\Vert \\left\\vert M \\right\\vert \\left\\vert M^{-1} \\right\\vert \\left\\vert x \\right\\vert \\right\\Vert_p
```

Skeel condition number ``\\kappa_S`` of the matrix `M`, optionally with respect to the
vector `x`, as computed using the operator `p`-norm. `p` is `Inf` by default, if not
provided. Valid values for `p` are `1`, `2`, or `Inf`.

This quantity is also known in the literature as the Bauer condition number, relative
condition number, or componentwise relative condition number.
"""
condskeel

"""
    det(M)

Matrix determinant.
"""
det

"""
    A_rdiv_Bt(A, B)

For matrices or vectors ``A`` and ``B``, calculates ``A / Bᵀ``.
"""
A_rdiv_Bt

"""
    qr(A [,pivot=Val{false}][;thin=true]) -> Q, R, [p]

Compute the (pivoted) QR factorization of `A` such that either `A = Q*R` or `A[:,p] = Q*R`.
Also see `qrfact`. The default is to compute a thin factorization. Note that `R` is not
extended with zeros when the full `Q` is requested.
"""
qr

"""
    Tridiagonal(dl, d, du)

Construct a tridiagonal matrix from the lower diagonal, diagonal, and upper diagonal,
respectively.  The result is of type `Tridiagonal` and provides efficient specialized linear
solvers, but may be converted into a regular matrix with [`full`](:func:`full`).
"""
Tridiagonal

"""
    bkfact(A) -> BunchKaufman

Compute the Bunch-Kaufman [^Bunch1977] factorization of a real symmetric or complex Hermitian
matrix `A` and return a `BunchKaufman` object. The following functions are available for
`BunchKaufman` objects: `size`, `\\`, `inv`, `issymmetric`, `ishermitian`.

[^Bunch1977]: J R Bunch and L Kaufman, Some stable methods for calculating inertia and solving symmetric linear systems, Mathematics of Computation 31:137 (1977), 163-179. [url](http://www.ams.org/journals/mcom/1977-31-137/S0025-5718-1977-0428694-0).

"""
bkfact

"""
    diff(A, [dim])

Finite difference operator of matrix or vector.
"""
diff

"""
    gradient(F, [h])

Compute differences along vector `F`, using `h` as the spacing between points. The default
spacing is one.
"""
gradient

"""
    logdet(M)

Log of matrix determinant. Equivalent to `log(det(M))`, but may provide increased accuracy
and/or speed.
"""
logdet

"""
    iseltype(A,T)

Tests whether `A` or its elements are of type `T`.
"""
iseltype

"""
    eigvecs(A, [eigvals,][permute=true,][scale=true]) -> Matrix

Returns a matrix `M` whose columns are the eigenvectors of `A`. (The `k`th eigenvector can
be obtained from the slice `M[:, k]`.) The `permute` and `scale` keywords are the same as
for [`eigfact`](:func:`eigfact`).

For [`SymTridiagonal`](:class:`SymTridiagonal`) matrices, if the optional vector of
eigenvalues `eigvals` is specified, returns the specific corresponding eigenvectors.
"""
eigvecs

"""
    scale(A, b)
    scale(b, A)

Scale an array `A` by a scalar `b`, returning a new array.

If `A` is a matrix and `b` is a vector, then `scale(A,b)` scales each column `i` of `A` by
`b[i]` (similar to `A*diagm(b)`), while `scale(b,A)` scales each row `i` of `A` by `b[i]`
(similar to `diagm(b)*A`), returning a new array.

Note: for large `A`, `scale` can be much faster than `A .* b` or `b .* A`, due to the use of BLAS.
"""
scale

"""
    scale!(A, b)
    scale!(b, A)

Scale an array `A` by a scalar `b`, similar to [`scale`](:func:`scale`) but overwriting `A`
in-place.

If `A` is a matrix and `b` is a vector, then `scale!(A,b)` scales each column `i` of `A` by
`b[i]` (similar to `A*diagm(b)`), while `scale!(b,A)` scales each row `i` of `A` by `b[i]`
(similar to `diagm(b)*A`), again operating in-place on `A`.

"""
scale!

"""
    factorize(A)

Compute a convenient factorization (including LU, Cholesky, Bunch-Kaufman, LowerTriangular,
UpperTriangular) of `A`, based upon the type of the input matrix. The return value can then
be reused for efficient solving of multiple systems. For example: `A=factorize(A); x=A\\b; y=A\\C`.
"""
factorize

"""
    vecdot(x, y)

For any iterable containers `x` and `y` (including arrays of any dimension) of numbers (or
any element type for which `dot` is defined), compute the Euclidean dot product (the sum of
`dot(x[i],y[i])`) as if they were vectors.
"""
vecdot

"""
    Ac_mul_B(A, B)

For matrices or vectors ``A`` and ``B``, calculates ``Aᴴ⋅B``.
"""
Ac_mul_B

"""
    qrfact!(A [,pivot=Val{false}])

`qrfact!` is the same as [`qrfact`](:func:`qrfact`) when `A` is a subtype of
`StridedMatrix`, but saves space by overwriting the input `A`, instead of creating a copy.
"""
qrfact!

"""
    At_rdiv_B(A, B)

For matrices or vectors ``A`` and ``B``, calculates ``Aᵀ / B``.
"""
At_rdiv_B

"""
    A_ldiv_Bt(A, B)

For matrices or vectors ``A`` and ``B``, calculates ``A`` \\ ``Bᵀ``.
"""
A_ldiv_Bt

"""

    eigvals(A,[irange,][vl,][vu]) -> values

Returns the eigenvalues of `A`. If `A` is `Symmetric`, `Hermitian` or `SymTridiagonal`,
it is possible to calculate only a subset of the eigenvalues by specifying either a
`UnitRange` `irange` covering indices of the sorted eigenvalues, or a pair `vl` and `vu`
for the lower and upper boundaries of the eigenvalues.

For general non-symmetric matrices it is possible to specify how the matrix is balanced
before the eigenvector calculation. The option `permute=true` permutes the matrix to
become closer to upper triangular, and `scale=true` scales the matrix by its diagonal
elements to make rows and columns moreequal in norm. The default is `true` for both
options.
"""
    eigvals

"""
    A_ldiv_Bc(A, B)

For matrices or vectors ``A`` and ``B``, calculates ``A`` \\ ``Bᴴ``.
"""
A_ldiv_Bc

"""
    qrfact(A [,pivot=Val{false}]) -> F

Computes the QR factorization of `A`. The return type of `F` depends on the element type of
`A` and whether pivoting is specified (with `pivot==Val{true}`).

| Return type   | `eltype(A)`     | `pivot`      | Relationship between `F` and `A` |
|:--------------|:----------------|:-------------|:---------------------------------|
| `QR`          | not `BlasFloat` | either       | `A==F[:Q]*F[:R]`                 |
| `QRCompactWY` | `BlasFloat`     | `Val{false}` | `A==F[:Q]*F[:R]`                 |
| `QRPivoted`   | `BlasFloat`     | `Val{true}`  | `A[:,F[:p]]==F[:Q]*F[:R]`        |

`BlasFloat` refers to any of: `Float32`, `Float64`, `Complex64` or `Complex128`.

The individual components of the factorization `F` can be accessed by indexing:

| Component | Description                               | `QR`            | `QRCompactWY`      | `QRPivoted`     |
|:----------|:------------------------------------------|:----------------|:-------------------|:----------------|
| `F[:Q]`   | `Q` (orthogonal/unitary) part of `QR`     | ✓ (`QRPackedQ`) | ✓ (`QRCompactWYQ`) | ✓ (`QRPackedQ`) |
| `F[:R]`   | `R` (upper right triangular) part of `QR` | ✓               | ✓                  | ✓               |
| `F[:p]`   | pivot `Vector`                            |                 |                    | ✓               |
| `F[:P]`   | (pivot) permutation `Matrix`              |                 |                    | ✓               |

The following functions are available for the `QR` objects: `size`, `\\`. When `A` is
rectangular, `\\` will return a least squares solution and if the solution is not unique,
the one with smallest norm is returned.

Multiplication with respect to either thin or full `Q` is allowed, i.e. both `F[:Q]*F[:R]`
and `F[:Q]*A` are supported. A `Q` matrix can be converted into a regular matrix with
[`full`](:func:`full`) which has a named argument `thin`.

**note**

`qrfact` returns multiple types because LAPACK uses several representations that minimize
the memory storage requirements of products of Householder elementary reflectors, so that
the `Q` and `R` matrices can be stored compactly rather as two separate dense matrices.

The data contained in `QR` or `QRPivoted` can be used to construct the `QRPackedQ` type,
which is a compact representation of the rotation matrix:

```math
Q = \\prod_{i=1}^{\\min(m,n)} (I - \\tau_i v_i v_i^T)
```

where ``\\tau_i`` is the scale factor and ``v_i`` is the projection vector associated with
the ``i^{th}`` Householder elementary reflector.

The data contained in `QRCompactWY` can be used to construct the `QRCompactWYQ` type,
which is a compact representation of the rotation matrix

```math
Q = I + Y T Y^T
```

where `Y` is ``m \\times r`` lower trapezoidal and `T` is ``r \\times r`` upper
triangular. The *compact WY* representation [^Schreiber1989] is not to be confused with the
older, *WY* representation [^Bischof1987]. (The LAPACK documentation uses `V` in lieu of `Y`.)

[^Bischof1987]: C Bischof and C Van Loan, "The WY representation for products of Householder matrices", SIAM J Sci Stat Comput 8 (1987), s2-s13. [doi:10.1137/0908009](http://dx.doi.org/10.1137/0908009)

[^Schreiber1989]: R Schreiber and C Van Loan, "A storage-efficient WY representation for products of Householder transformations", SIAM J Sci Stat Comput 10 (1989), 53-57. [doi:10.1137/0910005](http://dx.doi.org/10.1137/0910005)

"""
qrfact(A,?)


"""
    qrfact(A) -> SPQR.Factorization

Compute the QR factorization of a sparse matrix `A`. A fill-reducing permutation is used.
The main application of this type is to solve least squares problems with `\\`. The function
calls the C library SPQR and a few additional functions from the library are wrapped but not
exported.
"""
qrfact(A)

"""
    istriu(A) -> Bool

Test whether a matrix is upper triangular.
"""
istriu

"""
    A_mul_Bt(A, B)

For matrices or vectors ``A`` and ``B``, calculates ``A⋅Bᵀ``.
"""
A_mul_Bt

"""
    vecnorm(A, [p])

For any iterable container `A` (including arrays of any dimension) of numbers (or any
element type for which `norm` is defined), compute the `p`-norm (defaulting to `p=2`) as if
`A` were a vector of the corresponding length.

For example, if `A` is a matrix and `p=2`, then this is equivalent to the Frobenius norm.
"""
vecnorm

"""
    trace(M)

Matrix trace.
"""
trace

"""
    eigfact(A,[irange,][vl,][vu,][permute=true,][scale=true]) -> Eigen

Computes the eigenvalue decomposition of `A`, returning an `Eigen` factorization object `F`
which contains the eigenvalues in `F[:values]` and the eigenvectors in the columns of the
matrix `F[:vectors]`. (The `k`th eigenvector can be obtained from the slice `F[:vectors][:, k]`.)

The following functions are available for `Eigen` objects: `inv`, `det`.

If `A` is [`Symmetric`](:class:`Symmetric`), [`Hermitian`](:class:`Hermitian`) or
[`SymTridiagonal`](:class:`SymTridiagonal`), it is possible to calculate only a subset of
the eigenvalues by specifying either a [`UnitRange`](:class:`UnitRange`) `irange` covering
indices of the sorted eigenvalues or a pair `vl` and `vu` for the lower and upper boundaries
of the eigenvalues.

For general nonsymmetric matrices it is possible to specify how the matrix is balanced
before the eigenvector calculation. The option `permute=true` permutes the matrix to become
closer to upper triangular, and `scale=true` scales the matrix by its diagonal elements to
make rows and columns more equal in norm. The default is `true` for both options.
"""
eigfact(A,?,?,?,?)


"""
    eigfact(A, B) -> GeneralizedEigen

Computes the generalized eigenvalue decomposition of `A` and `B`, returning a
`GeneralizedEigen` factorization object `F` which contains the generalized eigenvalues in
`F[:values]` and the generalized eigenvectors in the columns of the matrix `F[:vectors]`.
(The `k`th generalized eigenvector can be obtained from the slice `F[:vectors][:, k]`.)
"""
eigfact(A,B)

"""
    eig(A,[irange,][vl,][vu,][permute=true,][scale=true]) -> D, V

Computes eigenvalues and eigenvectors of `A`. See [`eigfact`](:func:`eigfact`) for details
on the `permute` and `scale` keyword arguments. The eigenvectors are returned columnwise.

```jldoctest
julia> eig([1.0 0.0 0.0; 0.0 3.0 0.0; 0.0 0.0 18.0])
([1.0,3.0,18.0],
3x3 Array{Float64,2}:
 1.0  0.0  0.0
 0.0  1.0  0.0
 0.0  0.0  1.0)
```

`eig` is a wrapper around [`eigfact`](:func:`eigfact`), extracting all parts of the
factorization to a tuple; where possible, using [`eigfact`](:func:`eigfact`) is recommended.
"""
eig(A,?,?,?)

"""
    eig(A, B) -> D, V

Computes generalized eigenvalues and vectors of `A` with respect to `B`.

`eig` is a wrapper around [`eigfact`](:func:`eigfact`), extracting all parts of the
factorization to a tuple; where possible, using [`eigfact`](:func:`eigfact`) is recommended.
"""
eig(A,B)

"""
    At_ldiv_B(A, B)

For matrices or vectors ``A`` and ``B``, calculates ``Aᵀ`` \\ ``B``.
"""
At_ldiv_B

"""
    dot(x, y)
    ⋅(x,y)

Compute the dot product. For complex vectors, the first vector is conjugated.
"""
dot

"""
    cond(M, [p])

Condition number of the matrix `M`, computed using the operator `p`-norm. Valid values for
`p` are `1`, `2` (default), or `Inf`.
"""
cond


