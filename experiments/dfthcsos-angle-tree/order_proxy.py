#!/usr/bin/env python3
"""Reorder commuting Givens edges to reduce the signed-QROM active-bit proxy."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def active_bit_masks(words: np.ndarray, brot: int, batch_size: int):
    """Precompute the exact active difference-bit mask for every edge-column pair."""
    payload_bits = brot - 1
    payload_mask = (1 << payload_bits) - 1
    payload = np.asarray(words, dtype=np.int64) & payload_mask
    if payload.shape[0] % 2:
        payload = np.vstack((payload, payload[-1:]))

    pair_diff = payload[0::2] ^ payload[1::2]
    num_angles = payload.shape[1]
    column_mask = np.bitwise_or.reduce(pair_diff, axis=0)
    pair_mask = np.zeros((num_angles, num_angles), dtype=np.int64)
    for i in range(num_angles):
        pair_mask[i] = np.bitwise_or.reduce(pair_diff[:, i, None] ^ pair_diff, axis=0)
    return column_mask, pair_mask


def active_counts(
    permutation: np.ndarray, num_angles: int, batch_size: int,
    column_mask: np.ndarray, pair_mask: np.ndarray,
) -> tuple[int, ...]:
    """Count active packed payload bits in each XOR-delta batch.

    For batch k>0, the row-pair difference at a slot is the XOR of the
    pair-difference masks for the current and previous-batch columns.
    """
    num_batches = (num_angles + batch_size - 1) // batch_size
    padded = np.full(num_batches * batch_size, -1, dtype=np.int64)
    padded[:num_angles] = permutation
    counts = []
    for k in range(num_batches):
        count = 0
        for slot in range(batch_size):
            pos = k * batch_size + slot
            current = int(padded[pos])
            if k == 0:
                mask = 0 if current < 0 else int(column_mask[current])
            else:
                previous = int(padded[pos - batch_size])
                if current < 0 and previous < 0:
                    mask = 0
                elif current < 0:
                    mask = int(column_mask[previous])
                elif previous < 0:
                    mask = int(column_mask[current])
                else:
                    mask = int(pair_mask[current, previous])
            count += mask.bit_count()
        counts.append(count)
    return tuple(counts)


def score_counts(counts: tuple[int, ...]) -> tuple[int, int, int]:
    """Lexicographic objective: minimize peak workspace, then total active bits."""
    return max(counts, default=0), sum(counts), sum(x * x for x in counts)


def rooted_parent(edges: np.ndarray, num_orbitals: int) -> np.ndarray:
    parent = np.full(num_orbitals, -1, dtype=np.int64)
    seen = {0}
    for p, c in np.asarray(edges, dtype=np.int64):
        p, c = int(p), int(c)
        if p not in seen or c in seen:
            raise ValueError("tree_edges must be in parent-before-child order")
        parent[c] = p
        seen.add(c)
    if len(seen) != num_orbitals:
        raise ValueError("tree_edges do not span all orbitals")
    return parent


def is_ancestor(ancestor: int, node: int, parent: np.ndarray) -> bool:
    current = node
    while current >= 0:
        if current == ancestor:
            return True
        current = int(parent[current])
    return False


def can_swap_adjacent(
    edge_a: np.ndarray, edge_b: np.ndarray, parent: np.ndarray
) -> bool:
    p, c = map(int, edge_a)
    q, d = map(int, edge_b)
    if len({p, c, q, d}) != 4:
        return False
    # Disjoint rotations commute, but an ancestor edge must remain before
    # every descendant edge to satisfy the verifier's tree-order contract.
    if is_ancestor(c, q, parent) or is_ancestor(d, p, parent):
        return False
    return True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--brot", type=int, default=8)
    parser.add_argument("--batch-size", type=int, default=7)
    parser.add_argument("--max-steps", type=int, default=100)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.brot < 2 or args.batch_size < 1 or args.max_steps < 0:
        raise ValueError("invalid brot, batch-size, or max-steps")

    with np.load(args.source, allow_pickle=False) as loaded:
        data = {key: np.array(loaded[key]) for key in loaded.files}
    words = np.asarray(data["angle_words"], dtype=np.int64)
    edges = np.asarray(data["tree_edges"], dtype=np.int64)
    if words.ndim != 2 or edges.shape != (words.shape[1], 2):
        raise ValueError("angle_words and tree_edges shapes do not match")
    if np.any(words < 0) or np.any(words >= (1 << args.brot)):
        raise ValueError("angle words do not fit brot")
    if not np.all(((words >> (args.brot - 1)) & 1) == ((words >> (args.brot - 2)) & 1)):
        raise ValueError("expected signed-gauge words with equal top two bits")

    num_angles = words.shape[1]
    parent = rooted_parent(edges, words.shape[1] + 1)
    column_mask, pair_mask = active_bit_masks(words, args.brot, args.batch_size)
    permutation = np.arange(num_angles, dtype=np.int64)
    counts = active_counts(permutation, num_angles, args.batch_size, column_mask, pair_mask)
    initial_score = score_counts(counts)
    accepted = []

    for _ in range(args.max_steps):
        best_score = score_counts(counts)
        best_position = None
        best_counts = counts
        for pos in range(num_angles - 1):
            if not can_swap_adjacent(edges[permutation[pos]], edges[permutation[pos + 1]], parent):
                continue
            trial = permutation.copy()
            trial[pos], trial[pos + 1] = trial[pos + 1], trial[pos]
            trial_counts = active_counts(trial, num_angles, args.batch_size, column_mask, pair_mask)
            trial_score = score_counts(trial_counts)
            if trial_score < best_score:
                best_score = trial_score
                best_position = pos
                best_counts = trial_counts
        if best_position is None:
            break
        pos = best_position
        a, b = int(permutation[pos]), int(permutation[pos + 1])
        permutation[pos], permutation[pos + 1] = permutation[pos + 1], permutation[pos]
        counts = best_counts
        accepted.append({"position": pos, "swapped_columns": [a, b], "score": list(best_score)})

    reordered = dict(data)
    reordered["tree_edges"] = edges[permutation]
    reordered["angle_words"] = words[:, permutation]
    args.out.mkdir(parents=True, exist_ok=True)
    output_npz = args.out / "ordered_candidate.npz"
    np.savez_compressed(output_npz, **reordered)
    report = {
        "status": "ORDER_PROXY_SCREEN_ONLY",
        "source": str(args.source),
        "brot": args.brot,
        "batch_size": args.batch_size,
        "active_bits_before": list(initial_score),
        "per_batch_active_bits_before": list(active_counts(np.arange(num_angles), num_angles, args.batch_size, column_mask, pair_mask)),
        "active_bits_after": list(score_counts(counts)),
        "per_batch_active_bits_after": list(counts),
        "accepted_adjacent_swaps": len(accepted),
        "swaps": accepted,
        "artifact": output_npz.name,
        "warning": "This minimizes a select-swap QROM active-bit proxy, not exact Toffoli cost. Physical vectors should be unchanged because only disjoint commuting edges are swapped. Run the unchanged official verifier after building matching ANGLE circuits.",
    }
    (args.out / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
