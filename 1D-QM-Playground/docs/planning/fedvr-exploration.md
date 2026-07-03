# FEDVR: A Future Alternative to B-Splines
*Design exploration — 2026-07-03*
*Status: deferred. To be revisited once the B-spline code is working.*

---

## 1. What Is FEDVR?

**Finite-Element Discrete Variable Representation (FEDVR)** is a basis for 1D quantum mechanics that can be thought of as a special limiting case of the B-spline basis. It inherits B-splines' compact support and piecewise-polynomial structure, but is constructed so that two particularly useful properties emerge:

- The **potential energy matrix** is diagonal (not merely banded).
- The **overlap matrix is the identity** — the basis is orthonormal.

These properties open up significantly more efficient algorithms, particularly for time-dependent problems.

---

## 2. Relationship to B-Splines

Recall that a B-spline basis of order $k$ on an interval is built from a knot vector — a sequence of points on the $x$-axis with associated multiplicities. The multiplicity of an interior knot controls how smooth the B-spline functions are across that point: multiplicity 1 gives $C^{k-2}$ continuity (maximum smoothness), while multiplicity $k-1$ gives $C^0$ continuity (only continuous, not differentiable).

**FEDVR is the B-spline basis in the limit where every interior knot has multiplicity $k-1$.**

This maximum degeneracy of interior knots partitions the interval into independent "finite elements." Within each element $[x_i, x_{i+1}]$, there are exactly $k$ B-spline basis functions, and they are decoupled from those of neighboring elements except at the shared boundary knot (which connects adjacent elements through a single bridging function).

The result is a basis that looks like many small, locally-polynomial patches stitched together at their endpoints — a finite element method in the traditional sense, but with a very particular choice of local basis within each element.

---

## 3. The DVR Transformation: Diagonalizing Position

The "DVR" part of FEDVR comes from **Discrete Variable Representation**, a technique from molecular physics. The idea is to transform the local B-spline basis within each finite element into a new basis that diagonalizes the position operator $\hat{x}$.

Within a single element, the matrix representative of $\hat{x}$ is:

$$X_{ij} = \langle B_i \mid x \mid B_j \rangle = \int_{x_{\text{elem}}} B_i(x)\, x\, B_j(x)\, dx$$

Diagonalizing this matrix yields eigenvalues $\{x_\alpha\}$ and eigenvectors $\{|\delta_\alpha\rangle\}$. The eigenvalues are special quadrature points — specifically, the **Gauss-Lobatto points** of the element, which include both endpoints and a set of interior points that are optimal for numerical integration. The eigenfunctions $\delta_\alpha(x)$ are the FEDVR basis functions.

Because $|\delta_\alpha\rangle$ is an eigenstate of $\hat{x}$ with eigenvalue $x_\alpha$, each FEDVR function is, in a precise sense, "localized at" the point $x_\alpha$. This is the DVR property: the basis function $\delta_\alpha(x)$ is maximally concentrated around the single representative point $x_\alpha$ within the element.

---

## 4. Key Properties

### 4.1 Orthonormality

FEDVR basis functions satisfy:

$$\langle \delta_\alpha \mid \delta_\beta \rangle = \delta_{\alpha\beta}$$

This is in contrast to B-splines, where $\langle B_i \mid B_j \rangle = S_{ij}$ is a banded but non-trivial overlap matrix. With FEDVR, the generalized eigenvalue problem

$$\mathbf{H}\mathbf{c} = \mathbf{S}\mathbf{c}E$$

reduces to a standard eigenvalue problem

$$\mathbf{H}\mathbf{c} = \mathbf{c}E$$

because $\mathbf{S} = \mathbf{I}$.

### 4.2 Diagonal Potential Matrix

For a local (multiplicative) potential $V(x)$, the matrix element between two DVR functions is:

$$\langle \delta_\alpha \mid V(x) \mid \delta_\beta \rangle \approx V(x_\alpha)\, \delta_{\alpha\beta}$$

The approximation symbol $\approx$ appears because this equality holds exactly only when $V(x)$ is a polynomial of degree low enough to be integrated exactly by the Gauss-Lobatto quadrature rule. For smooth potentials, the error is exponentially small in the number of quadrature points per element and is negligible in practice.

The practical consequence is profound: **the potential energy matrix is diagonal.** Computing $V|\psi\rangle$ requires only $N$ multiplications (one per basis function), not a matrix-vector product.

### 4.3 Banded Kinetic Energy Matrix

The kinetic energy matrix $T_{\alpha\beta} = \langle \delta_\alpha \mid -\frac{1}{2}\frac{d^2}{dx^2} \mid \delta_\beta \rangle$ is not diagonal — the kinetic energy operator involves derivatives and couples basis functions. However, it remains **banded**: functions within the same finite element are coupled to each other, and functions in adjacent elements are coupled through the shared boundary DVR point. Functions in non-adjacent elements are exactly zero.

The Hamiltonian matrix in FEDVR is therefore:

$$H = T + V$$

where $T$ is banded and $V$ is diagonal — so $H$ itself is banded and can be stored and applied efficiently.

### 4.4 Summary of Matrix Structure

| Operator | B-spline matrix | FEDVR matrix |
|---|---|---|
| Overlap $S$ | Banded (bandwidth $k$) | Identity $\mathbf{I}$ |
| Potential $V(x)$ | Banded (bandwidth $k$) | **Diagonal** |
| Kinetic energy $T$ | Banded (bandwidth $k$) | Banded (block-tridiagonal) |
| Hamiltonian $H$ | Banded (bandwidth $k$) | Banded (block-tridiagonal) |

---

## 5. Algorithmic Advantages

The diagonal potential and orthonormal overlap unlock two classes of efficiency gains.

**For the TISE:** The standard eigenvalue problem (no $\mathbf{S}^{-1}$ factor) is solved by more efficient routines. The diagonal potential means the full Hamiltonian matrix is easier to assemble and apply.

**For the TDSE:** Time propagation via eigenstate expansion works identically to the B-spline case (diagonalize $H$, multiply each eigencoefficient by a phase). But FEDVR also opens the door to explicit propagation methods — algorithms that apply $H$ to $|\psi\rangle$ repeatedly without ever diagonalizing. Because $V|\psi\rangle$ is a diagonal scaling and $T|\psi\rangle$ is a banded matrix-vector product, each application of $H$ is $\mathcal{O}(N)$. Chebyshev or Krylov-based propagators (such as the Arnoldi/Lanczos method) become efficient for large systems where full diagonalization is prohibitive.

---

## 6. FEDVR as a Side Product of the B-Spline Code

The stakeholder noted that **FEDVR could emerge as a side product** once the B-spline solver is in place. This is because:

- FEDVR is constructed entirely within B-spline space — it is a linear transformation of the B-spline basis within each element.
- The knot vector and element structure are the same; only the local basis within each element changes (from the canonical B-spline functions to the position eigenstates).
- The assembly routines (quadrature, local matrix construction) are shared. The DVR transformation is a post-processing step applied to the local matrices.

In practice, the migration path from B-splines to FEDVR requires:
1. Choosing interior knot multiplicities of $k-1$ (instead of 1) in the knot vector.
2. Computing the position matrix $X_{ij}$ within each element and diagonalizing it.
3. Transforming all local matrices (kinetic energy, potential) from the B-spline basis to the DVR basis using the eigenvector matrix.
4. Assembling the global system in the DVR basis, respecting the shared boundary DVR points that connect adjacent elements.

The global solver (LAPACK eigenvalue routines, output routines, controller, analysis) is otherwise unchanged.

---

## 7. Status and Next Steps

This exploration is **deferred**. The current plan is:

1. Complete and validate the B-spline TISE and TDSE solvers.
2. Revisit FEDVR as an extension once the B-spline code is confirmed working.
3. Implement FEDVR as a configurable basis option alongside B-splines (e.g., a `basis.type: fedvr` flag in `config.yaml`).

The mathematical connection between the two bases means the effort to add FEDVR support should be incremental once the B-spline infrastructure is solid.
