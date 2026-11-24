#!/usr/bin/env python3
"""
Standalone drop-in replacement for ROSE_bamToGFF.py
Powered by pysam (htslib) + multiprocessing + atomic file writes.
Zero dependency on ROSE_utils or rose2 package.
"""

import argparse
import logging
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pysam

logger = logging.getLogger("fast_ROSE_bamToGFF")


def _parse_gff(gff_path: str) -> List[List[str]]:
    rows: List[List[str]] = []
    with open(gff_path, "r") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) >= 7:
                rows.append(cols[:9])
    return rows


def _resolve_bam_chrom(orig_chrom: str, ref_lengths: Dict[str, int]) -> Optional[str]:
    """Match GFF chromosome name against BAM header safely (handles chr / non-chr / custom contigs)."""
    if orig_chrom in ref_lengths:
        return orig_chrom
    if orig_chrom.startswith("chr"):
        alt = orig_chrom[3:]
        if alt in ref_lengths:
            return alt
    else:
        alt = "chr" + orig_chrom
        if alt in ref_lengths:
            return alt
    return None


def _compute_bins_floor0(
    bam: pysam.AlignmentFile,
    chrom: str,
    chrom_len: int,
    gff_start: int,
    gff_end: int,
    gff_strand: str,
    sense: str,
    extension: int,
    n_bins: int,
    bin_size: float,
    mmr: float,
) -> List[float]:
    if sense in ("both", "."):
        want_reverse: Optional[bool] = None
    else:
        gff_is_fwd = gff_strand in ("+", ".")
        want_reverse = (not gff_is_fwd) if sense == "+" else gff_is_fwd

    search_start = max(0, gff_start - extension - 1)
    search_end = min(chrom_len, gff_end + extension)
    if search_start >= search_end:
        return [0.0] * n_bins

    if n_bins == 1:
        total_cov = 0
        for read in bam.fetch(chrom, search_start, search_end):
            if read.is_unmapped:
                continue
            is_rev = read.is_reverse
            if want_reverse is not None and is_rev != want_reverse:
                continue

            r_start = read.reference_start + 1
            r_end = read.reference_end
            if r_end is None:
                continue

            if is_rev:
                r_start -= extension
            else:
                r_end += extension

            ov_s = r_start if r_start > gff_start else gff_start
            ov_e = r_end if r_end < gff_end else gff_end
            if ov_s <= ov_e:
                total_cov += ov_e - ov_s + 1

        return [round((total_cov / bin_size) / mmr, 4)]

    bin_bounds: List[Tuple[int, int]] = []
    if gff_strand in ("+", ".", "both"):
        cur = float(gff_start)
        for _ in range(n_bins):
            nxt = cur + bin_size
            bin_bounds.append((int(cur), int(nxt)))
            cur = nxt
        reverse_bins = False
    else:
        cur = float(gff_end)
        for _ in range(n_bins):
            prv = cur - bin_size
            bin_bounds.append((int(prv), int(cur)))
            cur = prv
        bin_bounds.reverse()
        reverse_bins = True

    bin_covs = [0] * n_bins
    inv_bin_size = 1.0 / bin_size

    for read in bam.fetch(chrom, search_start, search_end):
        if read.is_unmapped:
            continue
        is_rev = read.is_reverse
        if want_reverse is not None and is_rev != want_reverse:
            continue

        r_start = read.reference_start + 1
        r_end = read.reference_end
        if r_end is None:
            continue

        if is_rev:
            r_start -= extension
        else:
            r_end += extension

        ov_s = r_start if r_start > gff_start else gff_start
        ov_e = r_end if r_end < gff_end else gff_end
        if ov_s > ov_e:
            continue

        b0 = max(0, int((ov_s - gff_start) * inv_bin_size) - 1)
        b1 = min(n_bins - 1, int((ov_e - gff_start) * inv_bin_size) + 1)
        for b in range(b0, b1 + 1):
            bs, be = bin_bounds[b]
            s = ov_s if ov_s > bs else bs
            e = ov_e if ov_e < be else be
            if s <= e:
                bin_covs[b] += e - s + 1

    if reverse_bins:
        bin_covs.reverse()

    denom = bin_size * mmr
    return [round(cov / denom, 4) for cov in bin_covs]


def _compute_bins_with_floor(
    bam: pysam.AlignmentFile,
    chrom: str,
    chrom_len: int,
    gff_start: int,
    gff_end: int,
    gff_strand: str,
    sense: str,
    extension: int,
    floor: int,
    n_bins: int,
    bin_size: float,
    mmr: float,
) -> List[float]:
    if sense in ("both", "."):
        want_reverse: Optional[bool] = None
    else:
        gff_is_fwd = gff_strand in ("+", ".")
        want_reverse = (not gff_is_fwd) if sense == "+" else gff_is_fwd

    search_start = max(0, gff_start - extension - 1)
    search_end = min(chrom_len, gff_end + extension)
    if search_start >= search_end:
        return [0.0] * n_bins

    events: List[Tuple[int, int]] = []
    for read in bam.fetch(chrom, search_start, search_end):
        if read.is_unmapped:
            continue
        is_rev = read.is_reverse
        if want_reverse is not None and is_rev != want_reverse:
            continue

        r_start = read.reference_start + 1
        r_end = read.reference_end
        if r_end is None:
            continue

        if is_rev:
            r_start -= extension
        else:
            r_end += extension

        ov_s = r_start if r_start > gff_start else gff_start
        ov_e = r_end if r_end < gff_end else gff_end
        if ov_s <= ov_e:
            events.append((ov_s, 1))
            events.append((ov_e + 1, -1))

    if not events:
        return [0.0] * n_bins

    events.sort()

    intervals: List[Tuple[int, int, int]] = []
    depth = 0
    prev_pos: Optional[int] = None
    idx = 0
    n_ev = len(events)

    while idx < n_ev:
        pos = events[idx][0]
        if prev_pos is not None and depth > floor and pos > prev_pos:
            intervals.append((prev_pos, pos - 1, depth))
        while idx < n_ev and events[idx][0] == pos:
            depth += events[idx][1]
            idx += 1
        prev_pos = pos

    if n_bins == 1:
        total_cov = sum((iv_e - iv_s + 1) * d for iv_s, iv_e, d in intervals)
        return [round((total_cov / bin_size) / mmr, 4)]

    bin_bounds: List[Tuple[int, int]] = []
    if gff_strand in ("+", ".", "both"):
        cur = float(gff_start)
        for _ in range(n_bins):
            nxt = cur + bin_size
            bin_bounds.append((int(cur), int(nxt)))
            cur = nxt
        reverse_bins = False
    else:
        cur = float(gff_end)
        for _ in range(n_bins):
            prv = cur - bin_size
            bin_bounds.append((int(prv), int(cur)))
            cur = prv
        bin_bounds.reverse()
        reverse_bins = True

    bin_covs = [0] * n_bins
    inv_bin_size = 1.0 / bin_size

    for iv_s, iv_e, d in intervals:
        b0 = max(0, int((iv_s - gff_start) * inv_bin_size) - 1)
        b1 = min(n_bins - 1, int((iv_e - gff_start) * inv_bin_size) + 1)
        for b in range(b0, b1 + 1):
            bs, be = bin_bounds[b]
            s = iv_s if iv_s > bs else bs
            e = iv_e if iv_e < be else be
            if s <= e:
                bin_covs[b] += (e - s + 1) * d

    if reverse_bins:
        bin_covs.reverse()

    denom = bin_size * mmr
    return [round(cov / denom, 4) for cov in bin_covs]


def _process_chunk(
    args: Tuple[str, List[List[str]], str, int, int, int, float]
) -> List[List[Any]]:
    bam_file, gff_chunk, sense, extension, floor, n_bins, mmr = args
    results: List[List[Any]] = []

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        ref_lengths = dict(zip(bam.references, bam.lengths))

        for line in gff_chunk:
            orig_chrom = line[0]
            gff_id = line[1]
            gff_start = int(line[3])
            gff_end = int(line[4])
            gff_strand = line[6]

            # Preserve original GFF chromosome name so ROSE_main.py mapCollection matches 100%
            locus_str = f"{orig_chrom}({gff_strand}):{gff_start}-{gff_end}"
            cluster_line: List[Any] = [gff_id, locus_str]

            if n_bins > 0:
                bin_size = (gff_end - gff_start) / n_bins
                bam_chrom = _resolve_bam_chrom(orig_chrom, ref_lengths)

                if bin_size <= 0:
                    cluster_line.extend(["NA"] * n_bins)
                elif bam_chrom is None:
                    cluster_line.extend([0.0] * n_bins)
                elif floor > 0:
                    cluster_line.extend(
                        _compute_bins_with_floor(
                            bam, bam_chrom, ref_lengths[bam_chrom],
                            gff_start, gff_end, gff_strand,
                            sense, extension, floor, n_bins, bin_size, mmr,
                        )
                    )
                else:
                    cluster_line.extend(
                        _compute_bins_floor0(
                            bam, bam_chrom, ref_lengths[bam_chrom],
                            gff_start, gff_end, gff_strand,
                            sense, extension, n_bins, bin_size, mmr,
                        )
                    )

            results.append(cluster_line)

    return results


def main(argv: Optional[List[str]] = None) -> int:
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    parser = argparse.ArgumentParser(description="Fast pysam drop-in for ROSE_bamToGFF.py")
    parser.add_argument("-b", "--bam", required=True, help="Sorted and indexed BAM file")
    parser.add_argument("-i", "--input", required=True, help="Input .gff file")
    parser.add_argument("-o", "--output", default=None, help="Output filename")
    parser.add_argument("-s", "--sense", choices=["+", "-", ".", "both"], default="both")
    parser.add_argument("-f", "--floor", type=int, default=0)
    parser.add_argument("-e", "--extension", type=int, default=200)
    parser.add_argument("-r", "--rpm", action="store_true", default=False)
    parser.add_argument("-m", "--matrix", type=int, default=None)
    parser.add_argument("-p", "--threads", type=int, default=None)

    args = parser.parse_args(argv)

    bam_path = Path(args.bam).resolve()
    if not bam_path.exists():
        logger.error(f"BAM file not found: {args.bam}")
        return 1

    # Support both .bai and .csi indices
    has_idx = any(
        Path(cand).exists()
        for cand in (
            f"{bam_path}.bai",
            f"{bam_path}.csi",
            str(bam_path.with_suffix(".bai")),
            str(bam_path.with_suffix(".csi")),
        )
    )
    if not has_idx:
        logger.error("No .bai or .csi index found for BAM file.")
        return 1

    if not Path(args.input).exists():
        logger.error(f"Input GFF not found: {args.input}")
        return 1

    output = args.output or os.path.join(os.getcwd(), Path(args.input).name + ".mapped")
    n_bins = int(args.matrix) if args.matrix else 0

    # Optional escape hatch: set ROSE_IGNORE_FLOOR=1 to match pyrose2 v2.0 floor=0 behavior
    floor = 0 if os.environ.get("ROSE_IGNORE_FLOOR") == "1" else int(args.floor)

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        mmr = round(float(bam.mapped) / 1_000_000, 4) if args.rpm else 1.0

    gff_rows = _parse_gff(args.input)

    workers = args.threads if args.threads is not None else int(os.environ.get("ROSE_WORKERS", "1"))
    workers = max(1, min(workers, max(1, len(gff_rows))))

    logger.info(
        f"[fast-pysam] BAM={bam_path.name} MMR={mmr} regions={len(gff_rows)} "
        f"floor={floor} ext={args.extension} bins={n_bins} workers={workers}"
    )

    out_rows: List[List[Any]] = []
    if n_bins > 0:
        out_rows.append(
            ["GENE_ID", "locusLine"] + [f"bin_{n}_{bam_path.name}" for n in range(1, n_bins + 1)]
        )

    if gff_rows:
        if workers == 1:
            out_rows.extend(
                _process_chunk((str(bam_path), gff_rows, args.sense, args.extension, floor, n_bins, mmr))
            )
        else:
            n_chunks = min(len(gff_rows), workers * 4)
            chunk_size = (len(gff_rows) + n_chunks - 1) // n_chunks
            tasks = [
                (str(bam_path), gff_rows[i : i + chunk_size], args.sense, args.extension, floor, n_bins, mmr)
                for i in range(0, len(gff_rows), chunk_size)
            ]
            with ProcessPoolExecutor(max_workers=workers) as ex:
                for chunk_res in ex.map(_process_chunk, tasks):
                    out_rows.extend(chunk_res)

    # Atomic write: prevent ROSE_main.py polling loop from reading a half-written file
    out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_name(f".{out_path.name}.tmp.{os.getpid()}")

    with open(tmp_path, "w") as out_fh:
        for row in out_rows:
            out_fh.write("\t".join(map(str, row)) + "\n")

    os.replace(tmp_path, out_path)
    logger.info(f"[fast-pysam] Done -> {output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
