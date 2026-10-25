#!/usr/bin/env python3
"""Screen single-edge Givens-tree swaps with a quantization-sensitivity proxy.

This script is deliberately classical and does not invoke the official verifier.
Candidate NPZs are research artifacts until their matching ANGLE circuits pass it.
"""
from __future__ import annotations

import argparse
import heapq
import json
from pathlib import Path
from typing import Iterable

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = ROOT / "google-verifiable-quantum-circuit-compilation/sample_submission/dfthcsos/submission.npz"


def canonical_tree(edges: np.ndarray | Iterable[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(tuple(sorted((int(a), int(b)))) for a, b in np.asarray(list(edges), dtype=np.int64)))


def decode_words(words: np.ndarray, edges: np.ndarray, brot: int) -> np.ndarray:
    """Decode unsigned angle words; negative signed words are represented modulo 2**brot."""
    words = np.asarray(words, dtype=np.int64)
    edges = np.asarray(edges, dtype=np.int64)
    vectors = np.zeros((words.shape[0], words.shape[1] + 1), dtype=np.float64)
    vectors[:, 0] = 1.0
    for j, (parent, child) in enumerate(edges):
        vp = vectors[:, parent].copy()
        vc = vectors[:, child].copy()
        theta = np.pi * words[:, j] / (1 << (brot - 1))
        co, si = np.cos(theta), np.sin(theta)
        vectors[:, parent] = co * vp - si * vc
        vectors[:, child] = si * vp + co * vc
    return vectors


def encode_vectors(vectors: np.ndarray, ordered_edges: np.ndarray, brot: int) -> np.ndarray:
    """Inverse Givens elimination with a signed gauge and modular integer words."""
    work = np.asarray(vectors, dtype=np.float64).copy()
    edges = np.asarray(ordered_edges, dtype=np.int64)
    num_angles = work.shape[1] - 1
    if edges.shape != (num_angles, 2):
        raise ValueError(f"expected {(num_angles, 2)} ordered edges, got {edges.shape}")
    step = np.pi / (1 << (brot - 1))
    half = 1 << (brot - 2)
    signed_words = np.zeros((work.shape[0], num_angles), dtype=np.int64)
    for j in reversed(range(num_angles)):
        parent, child = edges[j]
        vp = work[:, parent].copy()
        vc = work[:, child].copy()
        theta = np.arctan2(vc, vp)
        flip = np.ones(work.shape[0], dtype=np.float64)

        # Shift angles outside [-pi/2, pi/2) by pi; the retained parent
        # amplitude changes sign, which is a gauge transformation of the
        # same final projector.
        high = theta >= np.pi / 2
        low = theta < -np.pi / 2
        theta[high] -= np.pi
        theta[low] += np.pi
        flip[high | low] = -1.0

        q = np.rint(theta / step).astype(np.int64)
        # +half has no signed (brot-1)-bit representation. Use -half and
        # flip the retained amplitude, the exact boundary gauge correction.
        boundary = q >= half
        q[boundary] = -half
        flip[boundary] *= -1.0
        q = np.clip(q, -half, half - 1)
        signed_words[:, j] = q
        work[:, parent] = flip * np.hypot(vp, vc)
        work[:, child] = 0.0

    words = signed_words % (1 << brot)
    if not np.all(((words >> (brot - 1)) & 1) == ((words >> (brot - 2)) & 1)):
        raise AssertionError("signed-gauge encoding must make the top two word bits equal")
    return words


def tree_order_and_proxy(
    tree: tuple[tuple[int, int], ...], orbital_weights: np.ndarray
) -> tuple[float, np.ndarray, np.ndarray]:
    """Root a tree at 0, sort siblings by descending subtree mass, and compute J(T)."""
    d = np.asarray(orbital_weights, dtype=np.float64)
    n = d.size
    if len(tree) != n - 1:
        raise ValueError(f"a tree on {n} vertices must have {n - 1} edges")
    adjacency: list[list[int]] = [[] for _ in range(n)]
    for a, b in tree:
        if a == b or not (0 <= a < n and 0 <= b < n):
            raise ValueError("invalid tree edge")
        adjacency[a].append(b)
        adjacency[b].append(a)

    parent = np.full(n, -2, dtype=np.int64)
    parent[0] = -1
    traversal = [0]
    for p in traversal:
        for c in adjacency[p]:
            if c == parent[p]:
                continue
            if parent[c] != -2:
                raise ValueError("edges contain a cycle")
            parent[c] = p
            traversal.append(c)
    if len(traversal) != n:
        raise ValueError("edges are disconnected")

    children: list[list[int]] = [[] for _ in range(n)]
    for v in range(1, n):
        children[int(parent[v])].append(v)

    mass = d.copy()
    for v in reversed(traversal[1:]):
        mass[int(parent[v])] += mass[v]
    for p in range(n):
        children[p].sort(key=lambda c: (-mass[c], c))

    proxy = 0.0
    for p in range(n):
        proxy += len(children[p]) * d[p]
        proxy += sum((rank + 1) * mass[c] for rank, c in enumerate(children[p]))

    ordered: list[tuple[int, int]] = []

    def visit(p: int) -> None:
        for c in children[p]:
            ordered.append((p, c))
            visit(c)

    visit(0)
    return float(proxy), np.asarray(ordered, dtype=np.int64), mass


def proxy_row_weights(vectors: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Build exact two-body first-order sensitivity weights without forming X or N**4 tensors.

    Let K = W W^T and H[b,b'] = (u_b dot u_b')**2. Then
    ||Y_b||_F**2 = (K H K)[b,b], where Y_b = sum_c w_bc X_c.
    One-body rows are balanced to the same total proxy weight; exact eigenvalue
    weights require the later physical evaluation.
    """
    n = vectors.shape[1]
    rank_r, bases_b, copies_c = w.shape
    if vectors.shape[0] != n + rank_r * bases_b:
        raise ValueError("angle table row count is inconsistent with W shape")
    u = vectors[n:].reshape(rank_r, bases_b, n)
    k = np.einsum("rbc,rdc->rbd", w, w, optimize=True)
    gram = np.einsum("rbp,rdp->rbd", u, u, optimize=True)
    projector_gram = np.square(gram)
    y_norm2 = np.einsum("rbi,rij,rbj->rb", k, projector_gram, k, optimize=True)
    two_body = 4.0 * np.maximum(y_norm2, 0.0).reshape(-1)
    total = float(two_body.sum())
    if not np.isfinite(total) or total <= 1e-30:
        two_body = np.ones(rank_r * bases_b, dtype=np.float64)
        total = float(two_body.sum())
    one_body = np.full(n, total / n, dtype=np.float64)
    return np.concatenate((one_body, two_body))


def one_edge_swaps(tree: tuple[tuple[int, int], ...], n: int):
    """Yield each tree obtainable by deleting one edge and adding one cross-cut edge."""
    base = set(tree)
    all_vertices = set(range(n))
    for removed in sorted(base):
        adjacency: list[list[int]] = [[] for _ in range(n)]
        for a, b in base:
            if (a, b) == removed:
                continue
            adjacency[a].append(b)
            adjacency[b].append(a)
        component = {removed[0]}
        stack = [removed[0]]
        while stack:
            p = stack.pop()
            for c in adjacency[p]:
                if c not in component:
                    component.add(c)
                    stack.append(c)
        other = all_vertices - component
        reduced = base - {removed}
        for a in sorted(component):
            for b in sorted(other):
                added = tuple(sorted((a, b)))
                if added in base:
                    continue
                candidate = tuple(sorted(reduced | {added}))
                yield candidate, removed, added


def weighted_projector_distortion(
    target: np.ndarray, decoded: np.ndarray, row_weights: np.ndarray
) -> float:
    """Sum_i weight_i ||u_i u_i^T - v_i v_i^T||_F^2, without forming projectors."""
    overlap = np.sum(target * decoded, axis=1)
    return float(np.sum(row_weights * 2.0 * np.maximum(0.0, 1.0 - np.square(overlap))))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--source-brot", type=int, default=12)
    parser.add_argument("--target-brot", type=int, default=8)
    parser.add_argument("--top-k", type=int, default=12)
    parser.add_argument("--max-swaps", type=int, default=0, help="0 means examine all single-edge swaps")
    parser.add_argument("--out", type=Path, default=Path("/tmp/dfthcsos-angle-tree"))
    args = parser.parse_args()
    if not 2 <= args.source_brot <= 16 or not 2 <= args.target_brot <= 16:
        raise ValueError("brot must be between 2 and 16")
    if args.top_k < 1:
        raise ValueError("top-k must be positive")

    with np.load(args.source, allow_pickle=False) as loaded:
        data = {key: np.array(loaded[key]) for key in loaded.files}
    if "angle_words" not in data or "tree_edges" not in data or "W" not in data:
        raise ValueError("source NPZ must contain angle_words, tree_edges, and W")

    source_words = data["angle_words"].astype(np.int64)
    if np.any(source_words < 0) or np.any(source_words >= (1 << args.source_brot)):
        raise ValueError("source angle words do not fit source-brot")
    source_tree = canonical_tree(data["tree_edges"])
    if "angle_tree_target_vectors" in data:
        target_vectors = np.asarray(data["angle_tree_target_vectors"], dtype=np.float64)
        expected_shape = (source_words.shape[0], source_words.shape[1] + 1)
        if target_vectors.shape != expected_shape:
            raise ValueError(f"saved target vectors must have shape {expected_shape}")
    else:
        target_vectors = decode_words(source_words, data["tree_edges"], args.source_brot)
    row_weights = proxy_row_weights(target_vectors, np.asarray(data["W"], dtype=np.float64))
    orbital_weights = np.sum(row_weights[:, None] * np.square(target_vectors), axis=0)

    base_proxy, base_order, _ = tree_order_and_proxy(source_tree, orbital_weights)
    base_words = encode_vectors(target_vectors, base_order, args.target_brot)
    base_decoded = decode_words(base_words, base_order, args.target_brot)
    base_distortion = weighted_projector_distortion(target_vectors, base_decoded, row_weights)

    heap: list[tuple[float, tuple[tuple[int, int], ...], tuple[int, int], tuple[int, int]]] = []
    examined = 0
    for candidate, removed, added in one_edge_swaps(source_tree, target_vectors.shape[1]):
        score, _, _ = tree_order_and_proxy(candidate, orbital_weights)
        examined += 1
        item = (-score, candidate, removed, added)
        if len(heap) < args.top_k:
            heapq.heappush(heap, item)
        elif score < -heap[0][0]:
            heapq.heapreplace(heap, item)
        if args.max_swaps and examined >= args.max_swaps:
            break

    shortlist = sorted(
        [(-neg_score, candidate, removed, added) for neg_score, candidate, removed, added in heap],
        key=lambda row: row[0],
    )
    args.out.mkdir(parents=True, exist_ok=True)
    screened = []
    for proxy_rank, (proxy, candidate, removed, added) in enumerate(shortlist, start=1):
        _, ordered_edges, _ = tree_order_and_proxy(candidate, orbital_weights)
        words = encode_vectors(target_vectors, ordered_edges, args.target_brot)
        decoded = decode_words(words, ordered_edges, args.target_brot)
        distortion = weighted_projector_distortion(target_vectors, decoded, row_weights)
        screened.append(
            {
                "proxy_rank": proxy_rank,
                "proxy_J": proxy,
                "proxy_ratio_vs_reordered_baseline": proxy / base_proxy if base_proxy else float("nan"),
                "weighted_projector_distortion": distortion,
                "removed_edge": list(removed),
                "added_edge": list(added),
                "tree_edges": ordered_edges,
                "angle_words": words,
                "decoded": decoded,
            }
        )

    screened.sort(key=lambda row: row["weighted_projector_distortion"])
    results = []
    for rank, row in enumerate(screened, start=1):
        candidate_data = dict(data)
        candidate_data["tree_edges"] = row["tree_edges"]
        candidate_data["angle_words"] = row["angle_words"]
        candidate_data["angle_tree_target_vectors"] = target_vectors
        rank_r, bases_b, _ = np.asarray(data["W"]).shape
        num_orbitals = target_vectors.shape[1]
        candidate_data["U"] = row["decoded"][num_orbitals:].reshape(rank_r, bases_b, num_orbitals)
        artifact = f"candidate_{rank:03d}.npz"
        np.savez_compressed(args.out / artifact, **candidate_data)
        results.append(
            {
                "rank_by_deterministic_distortion": rank,
                "rank_by_proxy": row["proxy_rank"],
                "proxy_J": row["proxy_J"],
                "proxy_ratio_vs_reordered_baseline": row["proxy_ratio_vs_reordered_baseline"],
                "weighted_projector_distortion": row["weighted_projector_distortion"],
                "removed_edge": row["removed_edge"],
                "added_edge": row["added_edge"],
                "tree_edges": row["tree_edges"].tolist(),
                "artifact": artifact,
            }
        )

    report = {
        "status": "SCREEN_ONLY",
        "source": str(args.source),
        "source_brot": args.source_brot,
        "target_brot": args.target_brot,
        "vectors": int(target_vectors.shape[0]),
        "orbitals": int(target_vectors.shape[1]),
        "single_edge_swaps_examined": examined,
        "baseline_proxy_J": base_proxy,
        "baseline_reordered_projector_distortion": base_distortion,
        "top_candidates": results,
        "warning": "No physical metrics, ANGLE circuit costs, or official verification were run.",
    }
    (args.out / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
