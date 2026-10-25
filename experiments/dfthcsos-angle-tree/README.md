# DFTHC-SOS angle-tree screening

This is a classical screening experiment for the Givens-tree representation. It does not run the quantum verifier and does not claim that any candidate is a valid submission.

## Mathematical proxy

For one target unit vector (u), an ordered rooted Givens tree distributes the norm of a parent across its child subtrees. If the children of node (p) are processed in order (c_1,\ldots,c_k), define

\[
m_v=\sum_{j\in\operatorname{subtree}(v)} a_j,
\qquad
a_j=\sum_i \omega_i u_{ij}^2.
\]

Immediately before the rotation on edge \((p,c_\ell)\), the squared amplitude still held by the parent is

\[
s_{p,c_\ell}^2=u_p^2+\sum_{j=\ell}^k
\|u_{\operatorname{subtree}(c_j)}\|^2.
\]

The derivative of the final vector with respect to that angle has norm \(s_{p,c_\ell}\), because all later Givens rotations are orthogonal. With independent, zero-mean rounding errors of variance \(\Delta^2/12\), the first-order expected projector error is proportional to

\[
J(T)=\sum_p\left(k_p a_p+
\sum_{\ell=1}^{k_p}\ell m_{c_\ell}\right).
\]

For a fixed topology, the rearrangement inequality proves that this proxy is minimized by processing children in **descending subtree mass** \(m_c\). The script uses this ordering to produce a valid parent-before-child edge list.

For the two-body part, it uses the local sensitivity matrix

\[
X_{rc}=\sum_b w_{rbc}u_{rb}u_{rb}^{\mathsf T},
\qquad
Y_{rb}=\sum_c w_{rbc}X_{rc}.
\]

The perturbation of the Gram tensor caused by a projector perturbation \(D\) in row \((r,b)\) is \(Y_{rb}\otimes D+D\otimes Y_{rb}\). Its squared norm is at most \(8\|Y_{rb}\|_F^2 s^2\), motivating the proxy weight \(4\|Y_{rb}\|_F^2\). The 76 one-body rows are balanced to the same total proxy weight at this screening stage; their exact eigenvalue-weighted error is checked later by the official verifier.

## Search procedure

1. Decode the current angle table to obtain target vectors.
2. Score every single-edge tree swap using \(J(T)\); no FCIDUMP parsing or quantum simulation is used in this pass.
3. For the best proxy candidates, re-encode all target vectors on each candidate tree at the requested angle precision and compute the deterministic weighted projector distortion.
4. Save the candidate factorization files for subsequent physical scoring and circuit compilation.

The proxy is an expected-error / sensitivity bound, not a prediction of the official score. Deterministic rounding, interactions among rows, the one-body eigenvector assignment, the SOS gap, normalization, and actual ANGLE-circuit cost can change the ranking. Only candidates that pass the unchanged official verifier count as valid.

The current prototype re-encodes shortlisted trees in the signed angle gauge, so at brot=8 the top two word bits are equal. It still does not compile ANGLE circuits or optimize the global topological ordering of edges; both affect lookup cost and belong to the next stage.

## Run

From the repository root:

```bash
python experiments/dfthcsos-angle-tree/tree_proxy.py \
  --source-brot 12 --target-brot 8 --top-k 12 \
  --out /tmp/dfthcsos-angle-tree
```

The script uses NumPy only. It writes a JSON report and one NPZ per shortlisted candidate. The NPZs are **screening artifacts**, not submissions until matching ANGLE circuits are built and the full verifier passes.

## Cheap exact contraction for the two-body sensitivity

The sensitivity norms are computed without materializing the large X matrices or any N^4 tensor. Define K_r = W_r W_r^T, G_r = U_r U_r^T, and H_r = G_r elementwise-multiplied by itself. Then

\[
\|Y_{rb}\|_F^2 = (K_r H_r K_r)_{bb}.
\]

This follows from \(\langle u_bu_b^T,u_{b'}u_{b'}^T\rangle=(u_b^T u_{b'})^2\). The contraction costs O(R B^3 + R B^2 N), so the screening stage avoids the much larger O(R B C N^2) matrix construction.
