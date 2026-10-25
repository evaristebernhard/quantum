#!/usr/bin/env python3
"""Search valid Givens-edge orderings with a signed-QROM cost proxy."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def active_bit_masks(words: np.ndarray, brot: int):
    """Precompute row-pair bit masks for every angle column and column pair."""
    payload_mask = (1 << (brot - 1)) - 1
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


def active_counts(permutation, num_angles, batch_size, column_mask, pair_mask):
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


def qrom_node_cost(words, permutation, brot, batch_size, lowbits):
    """Count recursive QROM split/leaf Toffoli proxies, excluding fixed controls.

    Equality of compressed table rows is exactly equality of the packed first
    rows and packed XOR differences. Prefix sums of adjacent signature changes
    make each subtree-constant query O(1).
    """
    payload_bits = brot - 1
    if batch_size * payload_bits > 63:
        raise ValueError("this proxy currently requires batch_size * (brot - 1) <= 63")
    num_rows, num_angles = words.shape
    num_batches = (num_angles + batch_size - 1) // batch_size
    padded = np.zeros((num_rows, num_batches * batch_size), dtype=np.int64)
    padded[:, :num_angles] = words[:, permutation]
    batches = padded.reshape(num_rows, num_batches, batch_size).transpose(1, 0, 2)
    deltas = batches.copy()
    deltas[1:] ^= batches[:-1]

    total_cost = 0
    per_batch = []
    payload_mask = (1 << payload_bits) - 1
    for k in range(num_batches):
        delta = deltas[k]
        if delta.shape[0] % 2:
            delta = np.vstack((delta, delta[-1:]))
        first = delta[0::2]
        diff = first ^ delta[1::2]
        pair_count = first.shape[0]
        base_pack = np.zeros(pair_count, dtype=np.uint64)
        diff_pack = np.zeros(pair_count, dtype=np.uint64)
        for slot in range(batch_size):
            shift = slot * payload_bits
            base_pack |= (first[:, slot] & payload_mask).astype(np.uint64) << shift
            diff_pack |= (diff[:, slot] & payload_mask).astype(np.uint64) << shift

        change = (base_pack[:-1] != base_pack[1:]) | (diff_pack[:-1] != diff_pack[1:])
        prefix = np.concatenate(([0], np.cumsum(change, dtype=np.int64)))
        table_rows = pair_count

        def constant_interval(lo: int, end: int) -> bool:
            return end - lo <= 1 or prefix[end - 1] - prefix[lo] == 0

        n = (table_rows - 1).bit_length()
        low = min(lowbits, n)
        high_bits = n - low
        cost = 0

        def node(lo: int, hi: int, offset: int) -> None:
            nonlocal cost
            if lo + offset >= table_rows:
                return
            end = min(hi + offset, table_rows)
            if hi - lo == 1 or constant_interval(lo + offset, end):
                return
            mid = (lo + hi) // 2
            if mid + offset >= table_rows:
                node(lo, mid, offset)
                return
            if hi - lo == 2:
                cost += 1
                return
            cost += 2
            node(lo, mid, offset)
            node(mid, hi, offset)

        for high in range(1 << high_bits):
            node(0, 1 << low, high << low)
        per_batch.append(cost)
        total_cost += cost
    return int(total_cost), tuple(int(x) for x in per_batch)


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


def build_precedence(edges: np.ndarray, parent: np.ndarray):
    """Parent edges precede child edges; sibling edges keep their existing order."""
    edge_count = len(edges)
    edge_index_by_child = {int(child): i for i, (_, child) in enumerate(edges)}
    successors = [set() for _ in range(edge_count)]
    indegree = np.zeros(edge_count, dtype=np.int64)

    def add(a: int, b: int) -> None:
        if b not in successors[a]:
            successors[a].add(b)
            indegree[b] += 1

    for i, (p, c) in enumerate(edges):
        p = int(p)
        if p != 0:
            add(edge_index_by_child[p], i)

    by_parent: dict[int, list[int]] = {}
    for i, (p, _) in enumerate(edges):
        by_parent.setdefault(int(p), []).append(i)
    for indices in by_parent.values():
        for a, b in zip(indices, indices[1:]):
            add(a, b)
    return successors, indegree


def random_topological_order(successors, base_indegree, rng):
    indegree = base_indegree.copy()
    ready = [i for i, degree in enumerate(indegree) if degree == 0]
    order = []
    while ready:
        pick = int(rng.integers(len(ready)))
        node = ready.pop(pick)
        order.append(node)
        for child in successors[node]:
            indegree[child] -= 1
            if indegree[child] == 0:
                ready.append(child)
    if len(order) != len(indegree):
        raise ValueError("edge precedence constraints contain a cycle")
    return np.asarray(order, dtype=np.int64)


def score_order(permutation, words, brot, batch_size, lowbits, column_mask, pair_mask):
    active = active_counts(permutation, words.shape[1], batch_size, column_mask, pair_mask)
    qrom_cost, per_batch_qrom = qrom_node_cost(words, permutation, brot, batch_size, lowbits)
    # Peak width first, then active select-swap bits, then recursive QROM Toffolis.
    score = (max(active, default=0), sum(active), qrom_cost, sum(x * x for x in active))
    return score, active, per_batch_qrom


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--brot", type=int, default=8)
    parser.add_argument("--batch-size", type=int, default=7)
    parser.add_argument("--lowbits", type=int, default=8)
    parser.add_argument("--samples", type=int, default=128)
    parser.add_argument("--seed", type=int, default=20261025)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.brot < 2 or args.batch_size < 1 or args.samples < 1:
        raise ValueError("invalid brot, batch-size, or samples")
    if args.batch_size * (args.brot - 1) > 63:
        raise ValueError("batch-size * (brot - 1) must be <= 63")

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

    parent = rooted_parent(edges, words.shape[1] + 1)
    successors, indegree = build_precedence(edges, parent)
    column_mask, pair_mask = active_bit_masks(words, args.brot)
    identity = np.arange(words.shape[1], dtype=np.int64)
    initial_score, initial_active, initial_qrom = score_order(
        identity, words, args.brot, args.batch_size, args.lowbits, column_mask, pair_mask
    )

    rng = np.random.default_rng(args.seed)
    seen = {tuple(identity.tolist())}
    best_perm = identity.copy()
    best_score, best_active, best_qrom = initial_score, initial_active, initial_qrom
    sampled = 0
    for _ in range(args.samples):
        permutation = random_topological_order(successors, indegree, rng)
        key = tuple(permutation.tolist())
        if key in seen:
            continue
        seen.add(key)
        sampled += 1
        score, active, qrom = score_order(
            permutation, words, args.brot, args.batch_size, args.lowbits, column_mask, pair_mask
        )
        if score < best_score:
            best_perm, best_score, best_active, best_qrom = permutation, score, active, qrom

    reordered = dict(data)
    reordered["tree_edges"] = edges[best_perm]
    reordered["angle_words"] = words[:, best_perm]
    args.out.mkdir(parents=True, exist_ok=True)
    output_npz = args.out / "ordered_candidate.npz"
    np.savez_compressed(output_npz, **reordered)
    report = {
        "status": "ORDER_PROXY_SCREEN_ONLY",
        "source": str(args.source),
        "brot": args.brot,
        "batch_size": args.batch_size,
        "lowbits": args.lowbits,
        "sampled_topological_orders": sampled,
        "initial_score": list(initial_score),
        "best_score": list(best_score),
        "active_bits_before": list(initial_active),
        "active_bits_after": list(best_active),
        "qrom_node_toffoli_proxy_before": int(sum(initial_qrom)),
        "qrom_node_toffoli_proxy_after": int(sum(best_qrom)),
        "per_batch_qrom_proxy_before": list(initial_qrom),
        "per_batch_qrom_proxy_after": list(best_qrom),
        "permutation": best_perm.tolist(),
        "artifact": output_npz.name,
        "warning": "This is a structural QROM proxy, not exact circuit cost or official score. Only disjoint commuting edges are reordered; parent-before-child and sibling order are preserved.",
    }
    (args.out / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
