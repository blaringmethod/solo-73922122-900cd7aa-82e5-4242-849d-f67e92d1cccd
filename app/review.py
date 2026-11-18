"""SNV 复核引擎：按比对位置从双端 FASTQ 结果提取位点证据。

设计约束（与需求逐条对应）：

* 只处理至多 ``REVIEW_MAX_PAIRS`` 个互异读号、每端不超过
  ``REVIEW_MAX_READ_LENGTH`` 碱基的结果；
* 完整流式扫描两份结果 FASTQ，后段损坏/缺件/重复或未知读号/非法比对
  一律整次报错（``ReviewError``），绝不返回部分统计；
* 结果文件以只读方式打开，任何路径都不写入；
* 内存中只保留被引用读号的记录与位点级计数，符合 256 MiB 进程上限。

证据规则：

* 反向链（reverse）先对序列反向互补、对质量串整体反转，再按 CIGAR 取碱基；
* CIGAR 仅允许 M/I/D/S，S 只能出现在首尾；M 推进读与参考、I 仅推进读、
  D 仅推进参考、S 仅推进读且不覆盖任何参考位点；
* 只有覆盖目标位点且 Phred+33 质量值达到阈值的碱基才算有效证据；
* 同一对读段的两端都有效时：碱基一致只计一次，碱基相异计一次冲突，
  不选一端、不重复计数；一端有效时按该端碱基分类，两端均无效则计入
  “无有效证据”。
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from . import config, db

# 单条 CIGAR：操作码仅 M/I/D/S，长度为正整数（首位不能为 0）。
_CIGAR_TOKEN_RE = re.compile(rb"([1-9]\d*)([MIDS])")
_COMPLEMENT = bytes.maketrans(b"ACGTacgtNn", b"TGCAtgcaNn")


class ReviewError(ValueError):
    """复核请求无法完成；整次失败，不产生任何部分统计。"""

    def __init__(self, message: str, *, status_code: int = 400) -> None:
        super().__init__(message)
        self.status_code = status_code


@dataclass(frozen=True)
class SiteSpec:
    position: int           # 从 1 开始的参考坐标
    ref: str                # 参考碱基 ACGT
    alt: str                # 候选 ALT 碱基 ACGT，且与 REF 不同


@dataclass(frozen=True)
class AlignmentSpec:
    read_id: str
    end: int                # 1 或 2
    ref_start: int          # CIGAR 第一个参考消耗操作覆盖的参考坐标（>=1）
    strand: str             # "forward" / "reverse"
    cigar: str


@dataclass(frozen=True)
class ReviewSpec:
    task_id: str
    min_quality: int
    sites: tuple[SiteSpec, ...]
    alignments: tuple[AlignmentSpec, ...]


@dataclass
class _ParsedAlignment:
    spec: AlignmentSpec
    # (长度, 操作码字节) 元组，操作码为 b"M"/b"I"/b"D"/b"S"
    ops: tuple[tuple[int, bytes], ...]
    query_consumed: int
    ref_span: int


@dataclass
class _SiteCounter:
    position: int
    ref: str
    alt: str
    ref_count: int = 0
    alt_count: int = 0
    other_count: int = 0
    conflict_count: int = 0
    no_evidence_count: int = 0


# --------------------------------------------------------------------------- #
# CIGAR 解析与校验
# --------------------------------------------------------------------------- #

def parse_cigar(cigar: str) -> tuple[tuple[tuple[int, bytes], ...], int, int]:
    """解析 CIGAR 字符串。

    返回 (操作序列, 读消耗总长, 参考消耗总长)。非法 CIGAR 抛 ReviewError。
    仅允许 M/I/D/S，S 只能出现在首尾，且必须至少含一个 M。
    """
    if not cigar:
        raise ReviewError("empty CIGAR")
    raw = cigar.encode("ascii", errors="strict") if cigar.isascii() else b""
    pos = 0
    ops: list[tuple[int, bytes]] = []
    query_consumed = 0
    ref_span = 0
    while pos < len(raw):
        match = _CIGAR_TOKEN_RE.match(raw, pos)
        if match is None or match.start() != pos:
            raise ReviewError(f"illegal CIGAR {cigar!r}: only M/I/D/S with "
                              "positive integer lengths are allowed")
        length = int(match.group(1))
        op = match.group(2)
        ops.append((length, op))
        if op in (b"M", b"I", b"S"):
            query_consumed += length
        if op in (b"M", b"D"):
            ref_span += length
        pos = match.end()
    if not ops:
        raise ReviewError(f"illegal CIGAR {cigar!r}: no operations")
    if all(op != b"M" for _, op in ops):
        raise ReviewError(f"illegal CIGAR {cigar!r}: at least one M operation is required")
    # S 只允许出现在首尾。
    interior = ops[1:-1] if len(ops) > 2 else []
    if any(op == b"S" for _, op in interior):
        raise ReviewError(f"illegal CIGAR {cigar!r}: soft clipping (S) is only "
                          "allowed at the leading or trailing edge")
    if query_consumed <= 0:
        raise ReviewError(f"illegal CIGAR {cigar!r}: consumes no query bases")
    return tuple(ops), query_consumed, ref_span


# --------------------------------------------------------------------------- #
# FASTQ 只读流式扫描
# --------------------------------------------------------------------------- #

def _iter_fastq_records(path: Path, label: str):
    """逐条产出 (读ID, 序列 bytes, 质量 bytes)，完整校验四行式结构。

    质量字符必须落在 Phred+33 可表示范围（ASCII 33..126）。任何结构性
    损坏（含文件后段截断）都抛 ReviewError。文件始终以只读方式打开。
    """
    with path.open("rb") as fh:
        record_no = 0
        while True:
            lines = [fh.readline(), fh.readline(), fh.readline(), fh.readline()]
            if not lines[0]:
                if any(lines[1:]):
                    raise ReviewError(
                        f"{label}: truncated final record (expected four lines)"
                    )
                if record_no == 0:
                    raise ReviewError(f"{label}: result FASTQ contains no records")
                return
            record_no += 1
            if any(not line for line in lines):
                raise ReviewError(
                    f"{label}: record {record_no} is missing lines "
                    "(FASTQ records must span exactly four lines)"
                )
            header, sequence, plus, quality = lines
            if not header.startswith(b"@"):
                raise ReviewError(
                    f"{label}: record {record_no} header does not start with '@'"
                )
            if not plus.startswith(b"+"):
                raise ReviewError(
                    f"{label}: record {record_no} third line does not start with '+'"
                )
            seq = sequence.rstrip(b"\r\n")
            qual = quality.rstrip(b"\r\n")
            if not seq or not qual:
                raise ReviewError(
                    f"{label}: record {record_no} has an empty sequence or quality string"
                )
            if len(seq) != len(qual):
                raise ReviewError(
                    f"{label}: record {record_no} sequence length {len(seq)} "
                    f"!= quality length {len(qual)}"
                )
            # Phred+33：合法质量字符为 ASCII 33..126；min/max 由 C 实现，扫描廉价。
            if qual and (qual[0] < 33 or max(qual) > 126):
                raise ReviewError(
                    f"{label}: record {record_no} contains a quality byte outside "
                    "the Phred+33 range (ASCII 33..126)"
                )
            first_field = header[1:].split(None, 1)[0] if len(header) > 1 else b""
            if not first_field:
                raise ReviewError(f"{label}: record {record_no} has an empty read id")
            if len(first_field) >= 2 and first_field[-2:-1] == b"/" and first_field[-1:] in (
                b"1", b"2",
            ):
                first_field = first_field[:-2]
            if not first_field:
                raise ReviewError(f"{label}: record {record_no} has an empty read id")
            yield first_field.decode("latin-1"), seq, qual


def _scan_side(
    path: Path,
    label: str,
    wanted: dict[str, int],
) -> dict[str, tuple[bytes, bytes]]:
    """扫描一侧结果 FASTQ，只保留 ``wanted`` 中的读号记录。

    ``wanted`` 为 {读ID: 该端允许的最大读长}；同时检查该文件内部读号重复。
    """
    found: dict[str, tuple[bytes, bytes]] = {}
    seen: set[str] = set()
    for read_id, seq, qual in _iter_fastq_records(path, label):
        if read_id in seen:
            raise ReviewError(f"{label}: duplicate read id {read_id!r} in result FASTQ")
        seen.add(read_id)
        if read_id in wanted:
            if len(seq) > wanted[read_id]:
                raise ReviewError(
                    f"{label}: read {read_id!r} is {len(seq)} bases long; the review "
                    f"limit is {config.REVIEW_MAX_READ_LENGTH} bases per end"
                )
            found[read_id] = (seq, qual)
    return found


# --------------------------------------------------------------------------- #
# 单条比对 -> 位点证据
# --------------------------------------------------------------------------- #

def _reverse_complement(seq: bytes) -> bytes:
    return seq.translate(_COMPLEMENT)[::-1]


def extract_evidence(
    alignment: _ParsedAlignment,
    seq: bytes,
    qual: bytes,
    site_positions: set[int],
    min_quality: int,
) -> dict[int, str]:
    """按比对返回 {参考坐标: 达标碱基}，仅含覆盖且质量达标的位点。"""
    spec = alignment.spec
    if spec.strand == "reverse":
        # 反向链：反向互补序列，并整体反转质量串后再按 CIGAR 同步。
        seq = _reverse_complement(seq)
        qual = qual[::-1]

    evidence: dict[int, str] = {}
    query_pos = 0
    ref_pos = spec.ref_start
    for length, op in alignment.ops:
        if op == b"M":
            for _ in range(length):
                if ref_pos in site_positions and qual[query_pos] - 33 >= min_quality:
                    evidence[ref_pos] = chr(seq[query_pos]).upper()
                query_pos += 1
                ref_pos += 1
        elif op == b"I":
            query_pos += length
        elif op == b"D":
            ref_pos += length
        else:  # b"S"：软裁剪，不覆盖任何参考位点
            query_pos += length
    return evidence


# --------------------------------------------------------------------------- #
# 复核入口
# --------------------------------------------------------------------------- #

def run_review(spec: ReviewSpec) -> dict:
    """执行整次复核；任何输入/数据问题抛 ReviewError，成功才返回统计。"""
    # 1) 任务必须存在且成功，结果文件齐备且位于本任务目录内（只读）。
    conn = db.connect()
    try:
        row = db.get_task(conn, spec.task_id)
    finally:
        conn.close()
    if row is None:
        raise ReviewError("task not found", status_code=404)
    if row["status"] != "succeeded":
        raise ReviewError(
            f"task is {row['status']}; reviews are only available for succeeded tasks",
            status_code=409,
        )
    base = config.TASKS_DIR / spec.task_id
    base_resolved = base.resolve()
    paths: dict[int, Path] = {}
    for end, rel in ((1, row["result_r1"]), (2, row["result_r2"])):
        if not rel:
            raise ReviewError(
                f"result R{end} file is missing for task {spec.task_id!r}",
                status_code=404,
            )
        path = (base / rel).resolve()
        if base_resolved not in path.parents or not path.is_file():
            raise ReviewError(
                f"result R{end} file is missing for task {spec.task_id!r}",
                status_code=404,
            )
        paths[end] = path

    # 2) 比对表整体校验：读号数上限、端别、链向、起点与 CIGAR 合法性。
    parsed: list[_ParsedAlignment] = []
    pair_keys: set[str] = set()
    seen_pairs: set[tuple[str, int]] = set()
    per_end_wanted: dict[int, dict[str, int]] = {1: {}, 2: {}}
    for item in spec.alignments:
        key = (item.read_id, item.end)
        if key in seen_pairs:
            raise ReviewError(
                f"duplicate alignment for read {item.read_id!r} end R{item.end}"
            )
        seen_pairs.add(key)
        pair_keys.add(item.read_id)
        if item.end not in (1, 2):
            raise ReviewError(f"alignment end must be 1 or 2, got {item.end}")
        if item.strand not in ("forward", "reverse"):
            raise ReviewError(
                f"strand must be 'forward' or 'reverse', got {item.strand!r}"
            )
        if item.ref_start < 1:
            raise ReviewError(
                f"ref_start must be a 1-based positive coordinate, got {item.ref_start}"
            )
        ops, query_consumed, _ref_span = parse_cigar(item.cigar)
        parsed.append(_ParsedAlignment(item, ops, query_consumed, _ref_span))
        # 同一读号同一端只可能出现一次（上面已查重），直接登记长度上限。
        per_end_wanted[item.end][item.read_id] = config.REVIEW_MAX_READ_LENGTH

    if len(pair_keys) > config.REVIEW_MAX_PAIRS:
        raise ReviewError(
            f"the alignment table references {len(pair_keys)} distinct read ids; "
            f"at most {config.REVIEW_MAX_PAIRS} read pairs can be reviewed"
        )

    # 3) 完整流式扫描两份结果（含结构与 Phred+33 范围校验），仅留存被引用记录。
    records: dict[int, dict[str, tuple[bytes, bytes]]] = {}
    for end in (1, 2):
        records[end] = _scan_side(paths[end], f"R{end}", per_end_wanted[end])

    # 4) 未知读号与 CIGAR/读长不一致均属非法比对，整次报错。
    for alignment in parsed:
        item = alignment.spec
        try:
            seq, qual = records[item.end][item.read_id]
        except KeyError:
            raise ReviewError(
                f"unknown read id {item.read_id!r} (not present in result R{item.end})"
            ) from None
        if alignment.query_consumed != len(seq):
            raise ReviewError(
                f"CIGAR {item.cigar!r} of read {item.read_id!r} R{item.end} consumes "
                f"{alignment.query_consumed} query bases but the read is {len(seq)} "
                "bases long"
            )

    # 5) 逐位点汇总。每个读号（一对读段）在每个位点至多贡献一次分类。
    site_positions = {site.position for site in spec.sites}
    counters = {
        site.position: _SiteCounter(site.position, site.ref, site.alt)
        for site in spec.sites
    }
    ref_by_pos = {site.position: site.ref for site in spec.sites}
    alt_by_pos = {site.position: site.alt for site in spec.sites}

    for read_id in sorted(pair_keys):
        # 一条比对在一个位点最多产生一个碱基；同一端不可能重复（已查重）。
        per_site: dict[int, list[str]] = {}
        for alignment in parsed:
            if alignment.spec.read_id != read_id:
                continue
            seq, qual = records[alignment.spec.end][read_id]
            evidence = extract_evidence(
                alignment, seq, qual, site_positions, spec.min_quality
            )
            for position, base in evidence.items():
                per_site.setdefault(position, []).append(base)

        for site in spec.sites:
            bases = per_site.get(site.position)
            counter = counters[site.position]
            if not bases:
                counter.no_evidence_count += 1
                continue
            if len(bases) == 2 and bases[0] != bases[1]:
                # 两端都有效却相异：冲突，不选一端、不重复计数。
                counter.conflict_count += 1
                continue
            base = bases[0]
            if base == ref_by_pos[site.position]:
                counter.ref_count += 1
            elif base == alt_by_pos[site.position]:
                counter.alt_count += 1
            else:
                counter.other_count += 1

    return {
        "task_id": spec.task_id,
        "pairs_total": len(pair_keys),
        "min_quality": spec.min_quality,
        "sites": [
            {
                "position": counter.position,
                "ref": counter.ref,
                "alt": counter.alt,
                "ref_count": counter.ref_count,
                "alt_count": counter.alt_count,
                "other_count": counter.other_count,
                "conflict_count": counter.conflict_count,
                "no_evidence_count": counter.no_evidence_count,
            }
            for counter in (counters[site.position] for site in spec.sites)
        ],
    }
