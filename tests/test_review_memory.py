"""复核路径内存压力测试：打满 2000 对 × 每端 500 碱基的上限输入。

在子进程中对一份 succeeded 任务的结果执行 run_review（引用全部 2000 个
读号、每端 500 碱基，外加 8 个位点），从 /proc/<pid>/status 读取 VmHWM，
要求峰值 RSS 低于 230 MiB，并核对位点统计的读对总数。

用法：.venv/bin/python -I tests/test_review_memory.py
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAIRS = 2000
READ_LEN = 500
RSS_LIMIT_KIB = 230 * 1024

CHILD = textwrap.dedent("""
    import sys
    sys.path.insert(0, %r)
    from pathlib import Path
    from app import config, db
    from app.review import (
        AlignmentSpec, ReviewSpec, SiteSpec, run_review,
    )

    data_dir = Path(sys.argv[1])
    config.DATA_DIR = data_dir
    config.DB_PATH = data_dir / "tasks.db"
    config.TASKS_DIR = data_dir / "tasks"
    db._initialized = False
    db.init_db()

    alignments = []
    for i in range(%d):
        rid = "read-%%06d" %% i
        alignments.append(AlignmentSpec(
            read_id=rid, end=1, ref_start=1,
            strand="forward" if i %% 2 == 0 else "reverse",
            cigar="%dM"))
        alignments.append(AlignmentSpec(
            read_id=rid, end=2, ref_start=1,
            strand="reverse" if i %% 2 == 0 else "forward",
            cigar="%dM"))
    sites = [SiteSpec(position=p, ref="A", alt="C") for p in (1, 50, 100,
            200, 250, 300, 400, 500)]
    spec = ReviewSpec(
        task_id="reviewstress", min_quality=20,
        sites=tuple(sites), alignments=tuple(alignments))
    result = run_review(spec)
    print("PAIRS=" + str(result["pairs_total"]))
    for site in result["sites"]:
        total = sum(site[k] for k in (
            "ref_count", "alt_count", "other_count",
            "conflict_count", "no_evidence_count"))
        print("SITE %%d TOTAL %%d" %% (site["position"], total))
""") % (str(ROOT), PAIRS, READ_LEN, READ_LEN)


def read_vmhwm(pid: int) -> int:
    for line in Path(f"/proc/{pid}/status").read_text().splitlines():
        if line.startswith("VmHWM:"):
            return int(line.split()[1])
    raise RuntimeError("VmHWM not available")


def write_results(base: Path) -> None:
    """直接写出一份结构合法的 succeeded 任务结果（绕过上传/配对流程）。"""
    out = base / "tasks" / "reviewstress" / "output"
    out.mkdir(parents=True)
    # 交替碱基保证反向互补后仍为 ACGT；质量全部 Q30('?')。
    seq1 = (b"ACGT" * ((READ_LEN // 4) + 1))[:READ_LEN]
    seq2 = (b"TGCA" * ((READ_LEN // 4) + 1))[:READ_LEN]
    qual = b"?" * READ_LEN
    with (out / "R1.fastq").open("wb") as f1, (out / "R2.fastq").open("wb") as f2:
        buf1 = bytearray()
        buf2 = bytearray()
        for i in range(PAIRS):
            rid = f"read-{i:06d}".encode()
            buf1 += b"@" + rid + b"/1\n" + seq1 + b"\n+\n" + qual + b"\n"
            buf2 += b"@" + rid + b"/2\n" + seq2 + b"\n+\n" + qual + b"\n"
            if len(buf1) >= 1024 * 1024:
                f1.write(buf1)
                f2.write(buf2)
                buf1.clear()
                buf2.clear()
        f1.write(buf1)
        f2.write(buf2)


def main() -> int:
    scratch = Path(tempfile.mkdtemp(prefix="fastq-review-mem-",
                                    dir=str(ROOT / "tests" / "output")))
    data_dir = scratch / "data"
    data_dir.mkdir(parents=True)
    write_results(data_dir)

    # 登记 succeeded 任务元数据，result_r1/r2 指向输出文件。
    sys.path.insert(0, str(ROOT))
    from app import config, db  # noqa: E402
    config.DATA_DIR = data_dir
    config.DB_PATH = data_dir / "tasks.db"
    config.TASKS_DIR = data_dir / "tasks"
    db._initialized = False
    db.init_db()
    conn = db.connect()
    db.create_task(conn, "reviewstress", "R1.fastq", "R2.fastq", 0, 0)
    db.update_task(
        conn, "reviewstress", status="succeeded", stage=config.STAGE_DONE,
        progress=1.0, message="synthetic",
        result_r1="output/R1.fastq", result_r2="output/R2.fastq",
    )
    conn.close()

    env = dict(os.environ)
    env["PYTHONIOENCODING"] = "utf-8"
    proc = subprocess.Popen(
        [sys.executable, "-I", "-c", CHILD, str(data_dir)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        env=env, text=True,
    )
    peak = 0
    while proc.poll() is None:
        try:
            peak = max(peak, read_vmhwm(proc.pid))
        except (FileNotFoundError, ProcessLookupError):
            break
        time.sleep(0.05)
    out, err = proc.communicate()
    print(out.strip())
    if proc.returncode != 0:
        print(err[-4000:])
        print("RESULT: FAIL (child exited %d)" % proc.returncode)
        return 1

    pairs_line = next((l for l in out.splitlines() if l.startswith("PAIRS=")), "")
    ok = pairs_line == f"PAIRS={PAIRS}"
    site_lines = [l for l in out.splitlines() if l.startswith("SITE ")]
    ok = ok and len(site_lines) == 8 and all(
        l.endswith(f"TOTAL {PAIRS}") for l in site_lines
    )
    print(f"peak RSS: {peak / 1024:.1f} MiB "
          f"(limit {RSS_LIMIT_KIB / 1024:.0f} MiB)")
    if peak > RSS_LIMIT_KIB:
        print("RESULT: FAIL (memory cap exceeded)")
        ok = False
    print("RESULT: PASS" if ok else "RESULT: FAIL (count verification)")
    shutil.rmtree(scratch, ignore_errors=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
