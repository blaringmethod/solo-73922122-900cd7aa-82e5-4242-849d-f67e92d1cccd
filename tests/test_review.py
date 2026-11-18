"""SNV 复核 API 的功能自测。

通过真实的 /tasks 上传流程产出 succeeded 任务，再调用
POST /tasks/{id}/reviews 验证证据提取、整次失败语义与只读保证。
运行：.venv/bin/python -I tests/test_review.py
"""

from __future__ import annotations

import asyncio
import io
import os
import shutil
import sys
import tempfile
import time
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import httpx  # noqa: E402

from app import config, db, worker  # noqa: E402
from app.main import app  # noqa: E402

RESULTS: list[tuple[str, bool, str]] = []


def record(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" - {detail}" if detail else ""))


def make_record(rid: str, seq: bytes, quality: bytes,
                suffix: str = "") -> bytes:
    header = b"@" + rid.encode() + suffix.encode() + b"\n"
    return header + seq + b"\n+\n" + quality + b"\n"


def fastq(records: list[bytes]) -> io.BytesIO:
    return io.BytesIO(b"".join(records))


def q(chars: str) -> bytes:
    """质量串：'!' = Q0，'5' = Q20，'?' = Q30，'I' = Q40。"""
    return chars.encode("ascii")


async def wait_terminal(client: httpx.AsyncClient, task_id: str,
                        timeout: float = 30.0) -> dict:
    deadline = time.time() + timeout
    while time.time() < deadline:
        body = (await client.get(f"/tasks/{task_id}")).json()
        if body["status"] in ("succeeded", "failed", "cancelled"):
            return body
        await asyncio.sleep(0.05)
    raise AssertionError(f"task {task_id} did not finish in {timeout}s")


async def upload_pair(client: httpx.AsyncClient,
                      r1: list[bytes], r2: list[bytes]) -> str:
    resp = await client.post(
        "/tasks",
        files={"r1": ("r1.fastq", fastq(r1), "application/octet-stream"),
               "r2": ("r2.fastq", fastq(r2), "application/octet-stream")},
    )
    assert resp.status_code == 202, resp.text
    body = await wait_terminal(client, resp.json()["id"])
    assert body["status"] == "succeeded", body
    return body["id"]


def site_counter(payload: dict, position: int) -> dict:
    for site in payload["sites"]:
        if site["position"] == position:
            return site
    raise KeyError(position)


async def post_review(client: httpx.AsyncClient, task_id: str,
                      body: dict) -> httpx.Response:
    return await client.post(f"/tasks/{task_id}/reviews", json=body)


async def main_async() -> int:
    tmp_root = Path(tempfile.mkdtemp(prefix="fastq-review-",
                                     dir=str(ROOT / "tests" / "output")))
    data_dir = tmp_root / "data"
    os.environ["FASTQ_DATA_DIR"] = str(data_dir)
    config.DATA_DIR = data_dir
    config.DB_PATH = data_dir / "tasks.db"
    config.TASKS_DIR = data_dir / "tasks"

    db._initialized = False
    db.init_db()
    worker.recover_on_startup()

    transport = httpx.ASGITransport(app=app)
    async with httpx.AsyncClient(transport=transport, base_url="http://test",
                                 timeout=60) as client:
        # ------------------------------------------------------------------ #
        # 数据集：
        #  rA R1 ACGTACGT (全 Q30)，R2 TGCATGCA 全 Q30，反向链
        #  rB R1 acgtacgt 小写、质量 5,5,5,5,!,!,!,!（后四个 Q0）
        #     R2 GGGG.... 设计为与 R1 在位点 5 冲突
        #  rC R1 含插入/缺失/软裁剪的复杂 CIGAR
        #  rD R1 含 N 碱基；R2 缺比对（只给 R1 一条）
        # ------------------------------------------------------------------ #
        r1 = [
            make_record("rA", b"ACGTACGT", q("????????")),
            make_record("rB", b"acgtacgt", q("5555!!!!")),
            make_record("rC", b"TTACGTAGT", q("?????????")),
            make_record("rD", b"ACGTNCGT", q("????????")),
        ]
        r2 = [
            # rA 反向链：取反向互补后为 TGCATGCA（与存储序列一致设计）
            make_record("rA", b"TGCATGCA", q("????????")),
            # rB: 位点 5 的碱基为 G（R1 在 5 是 a=A），制造冲突
            make_record("rB", b"TTTTGTTT", q("????????")),
            make_record("rC", b"ACGTACGT", q("????????")),
            make_record("rD", b"ACGTACGT", q("????????")),
        ]
        tid = await upload_pair(client, r1, r2)
        result_dir = config.TASKS_DIR / tid / "output"
        r1_path = result_dir / "R1.fastq"
        r2_path = result_dir / "R2.fastq"
        r1_before, r2_before = r1_path.read_bytes(), r2_path.read_bytes()

        # ---------- 1. 基础正向链 + 质量阈值 + REF/ALT/其他/无证据 ---------- #
        body = {
            "min_quality": 20,
            "sites": [
                {"position": 5, "ref": "A", "alt": "T"},
                {"position": 7, "ref": "G", "alt": "A"},
            ],
            "alignments": [
                {"read_id": "rA", "end": 1, "ref_start": 1,
                 "strand": "forward", "cigar": "8M"},
                {"read_id": "rB", "end": 1, "ref_start": 1,
                 "strand": "forward", "cigar": "8M"},
                {"read_id": "rC", "end": 1, "ref_start": 1,
                 "strand": "forward", "cigar": "2S2M2D2M2I1M"},
                {"read_id": "rD", "end": 1, "ref_start": 1,
                 "strand": "forward", "cigar": "8M"},
            ],
        }
        resp = await post_review(client, tid, body)
        record("basic review returns 200", resp.status_code == 200, resp.text)
        data = resp.json()
        record("pairs_total counts distinct read ids",
               data["pairs_total"] == 4, str(data["pairs_total"]))

        # 位点 5（ref A / alt T，阈值 Q20）：
        #  rA: R1 位5=A 达标 -> REF
        #  rB: R1 位5 质量 Q0 -> 无有效证据
        #  rC: CIGAR 2S 2M(ref1-2) 2D(ref3-4) 2M(ref5-6) 2I 1M(ref7)，
        #      位5 取 query index 4 = C，既非 REF(A) 也非 ALT(T) -> other
        #  rD: 位5=N 达标 -> other
        c5 = site_counter(data, 5)
        record("site5 REF count", c5["ref_count"] == 1, str(c5))
        record("site5 ALT count zero", c5["alt_count"] == 0, str(c5))
        record("site5 other counts C and N", c5["other_count"] == 2, str(c5))
        record("site5 no-evidence counts low quality",
               c5["no_evidence_count"] == 1, str(c5))
        record("site5 conflict zero", c5["conflict_count"] == 0, str(c5))
        totals5 = sum(c5[k] for k in (
            "ref_count", "alt_count", "other_count",
            "conflict_count", "no_evidence_count"))
        record("site5 buckets sum to pairs_total", totals5 == 4, str(totals5))

        # 位点 7：rA 位7=G=REF；rB 位7 质量 Q0 无证据；
        # rC 在 2I 后由最后的 1M 覆盖 ref7，取 query index 8=T，ALT=A -> other；
        # rD 位7=G=REF。
        c7 = site_counter(data, 7)
        record("site7 insertion shifts query only",
               c7["other_count"] == 1 and c7["no_evidence_count"] == 1
               and c7["ref_count"] == 2 and c7["alt_count"] == 0, str(c7))

        # ---------- 2. 反向链：反向互补 + 质量反转，双端一致只计一次 ---------- #
        body = {
            "min_quality": 20,
            "sites": [{"position": 1, "ref": "T", "alt": "C"}],
            "alignments": [
                {"read_id": "rA", "end": 1, "ref_start": 1,
                 "strand": "forward", "cigar": "8M"},
                {"read_id": "rA", "end": 2, "ref_start": 1,
                 "strand": "reverse", "cigar": "8M"},
            ],
        }
        resp = await post_review(client, tid, body)
        record("reverse review 200", resp.status_code == 200, resp.text)
        c1 = site_counter(resp.json(), 1)
        # R1 位1=A（other）；R2 反向互补 TGCATGCA -> TGCATGCA 位1=T，与 R1 相异
        # => 冲突（两端都有效）
        record("reverse complement applied",
               c1["conflict_count"] == 1 and c1["other_count"] == 0, str(c1))

        # 把 R1 改成与反向链结果一致：rA R1 换成与 R2 rc 后同碱基的设计
        # 直接验证“双端一致只计一次”：两端都报 T(ref)
        body = {
            "min_quality": 20,
            "sites": [{"position": 1, "ref": "T", "alt": "C"}],
            "alignments": [
                {"read_id": "rA", "end": 2, "ref_start": 1,
                 "strand": "reverse", "cigar": "8M"},
            ],
        }
        resp = await post_review(client, tid, body)
        c1 = site_counter(resp.json(), 1)
        record("single-end reverse strand evidence classified as REF",
               c1["ref_count"] == 1, str(c1))

        # ---------- 3. 双端冲突：rB 位5 R1=A、R2=G，都达标即冲突 ---------- #
        body = {
            "min_quality": 0,
            "sites": [{"position": 5, "ref": "A", "alt": "G"}],
            "alignments": [
                {"read_id": "rB", "end": 1, "ref_start": 1,
                 "strand": "forward", "cigar": "8M"},
                {"read_id": "rB", "end": 2, "ref_start": 1,
                 "strand": "forward", "cigar": "8M"},
            ],
        }
        resp = await post_review(client, tid, body)
        c5 = site_counter(resp.json(), 5)
        record("discordant mates counted once as conflict",
               c5["conflict_count"] == 1 and c5["ref_count"] == 0
               and c5["alt_count"] == 0, str(c5))

        # 同阈值下只给 R1（A=REF），不得重复计数
        body["alignments"].pop()
        resp = await post_review(client, tid, body)
        c5 = site_counter(resp.json(), 5)
        record("one-end evidence counted once",
               c5["ref_count"] == 1 and c5["conflict_count"] == 0, str(c5))

        # ---------- 4. 质量阈值边界：Q20 字符 '5' 恰好达标 ---------- #
        body = {
            "min_quality": 20,
            "sites": [{"position": 1, "ref": "A", "alt": "C"}],
            "alignments": [
                {"read_id": "rB", "end": 1, "ref_start": 1,
                 "strand": "forward", "cigar": "8M"},
            ],
        }
        resp = await post_review(client, tid, body)
        c1 = site_counter(resp.json(), 1)
        record("Q20 boundary accepted", c1["ref_count"] == 1, str(c1))
        body["min_quality"] = 21
        resp = await post_review(client, tid, body)
        c1 = site_counter(resp.json(), 1)
        record("Q20 fails threshold 21 -> no evidence",
               c1["no_evidence_count"] == 1, str(c1))

        # ---------- 5. 未知读号 / 重复 (read_id,end) ---------- #
        body = {
            "min_quality": 10,
            "sites": [{"position": 1, "ref": "A", "alt": "C"}],
            "alignments": [
                {"read_id": "ghost", "end": 1, "ref_start": 1,
                 "strand": "forward", "cigar": "8M"},
            ],
        }
        resp = await post_review(client, tid, body)
        record("unknown read id -> 400 whole-request error",
               resp.status_code == 400 and "ghost" in resp.text, resp.text)
        record("no partial stats on error",
               "sites" not in resp.json(), resp.text)

        body = {
            "min_quality": 10,
            "sites": [{"position": 1, "ref": "A", "alt": "C"}],
            "alignments": [
                {"read_id": "rA", "end": 1, "ref_start": 1,
                 "strand": "forward", "cigar": "8M"},
                {"read_id": "rA", "end": 1, "ref_start": 2,
                 "strand": "forward", "cigar": "8M"},
            ],
        }
        resp = await post_review(client, tid, body)
        record("duplicate (read_id,end) -> 422",
               resp.status_code == 422, resp.text)

        # ---------- 6. 非法 CIGAR ---------- #
        for bad_cigar in ("8X", "0M", "M8", "2S2M2S3M", "8I", "8S", "3M2S2M"):
            body = {
                "min_quality": 10,
                "sites": [{"position": 1, "ref": "A", "alt": "C"}],
                "alignments": [
                    {"read_id": "rA", "end": 1, "ref_start": 1,
                     "strand": "forward", "cigar": bad_cigar},
                ],
            }
            resp = await post_review(client, tid, body)
            record(f"illegal CIGAR {bad_cigar!r} -> 422",
                   resp.status_code == 422, resp.text)

        # 词法合法但与读长不符 -> 引擎层 400
        body = {
            "min_quality": 10,
            "sites": [{"position": 1, "ref": "A", "alt": "C"}],
            "alignments": [
                {"read_id": "rA", "end": 1, "ref_start": 1,
                 "strand": "forward", "cigar": "7M"},
            ],
        }
        resp = await post_review(client, tid, body)
        record("CIGAR/read length mismatch -> 400",
               resp.status_code == 400, resp.text)

        # ---------- 7. schema 边界：位点数、质量值、allele、坐标重复 ---------- #
        bad_sites = {
            "min_quality": 10,
            "sites": [{"position": 1, "ref": "A", "alt": "A"}],
            "alignments": body["alignments"],
        }
        resp = await post_review(client, tid, bad_sites)
        record("ref==alt -> 422", resp.status_code == 422, resp.text)

        bad_sites = {
            "min_quality": 41,
            "sites": [{"position": 1, "ref": "A", "alt": "C"}],
            "alignments": body["alignments"],
        }
        resp = await post_review(client, tid, bad_sites)
        record("min_quality 41 -> 422", resp.status_code == 422, resp.text)

        bad_sites = {
            "min_quality": 10,
            "sites": [
                {"position": 1, "ref": "A", "alt": "C"},
                {"position": 1, "ref": "C", "alt": "A"},
            ],
            "alignments": body["alignments"],
        }
        resp = await post_review(client, tid, bad_sites)
        record("duplicate site position -> 422",
               resp.status_code == 422, resp.text)

        bad_sites = {
            "min_quality": 10,
            "sites": [{"position": 0, "ref": "A", "alt": "C"}],
            "alignments": body["alignments"],
        }
        resp = await post_review(client, tid, bad_sites)
        record("position 0 -> 422", resp.status_code == 422, resp.text)

        too_many = {
            "min_quality": 10,
            "sites": [
                {"position": i + 1, "ref": "A", "alt": "C"}
                for i in range(config.REVIEW_MAX_SITES + 1)
            ],
            "alignments": body["alignments"],
        }
        resp = await post_review(client, tid, too_many)
        record("more than 8 sites -> 422", resp.status_code == 422, resp.text)

        # ---------- 8. 2000 对上限 ---------- #
        over = config.REVIEW_MAX_PAIRS + 1
        body = {
            "min_quality": 10,
            "sites": [{"position": 1, "ref": "A", "alt": "C"}],
            "alignments": [
                {"read_id": f"n{i}", "end": 1, "ref_start": 1,
                 "strand": "forward", "cigar": "8M"}
                for i in range(over)
            ],
        }
        resp = await post_review(client, tid, body)
        record(">2000 distinct read ids -> 422",
               resp.status_code == 422, resp.text)

        # ---------- 9. 任务状态与缺件 ---------- #
        valid_body = {
            "min_quality": 10,
            "sites": [{"position": 1, "ref": "A", "alt": "C"}],
            "alignments": [
                {"read_id": "rA", "end": 1, "ref_start": 1,
                 "strand": "forward", "cigar": "8M"},
            ],
        }
        unknown = await post_review(
            client, uuid.uuid4().hex, valid_body)
        record("review unknown task -> 404",
               unknown.status_code == 404, unknown.text)

        # 未成功任务：直接构造一条 processing 任务
        stuck = uuid.uuid4().hex
        conn = db.connect()
        db.create_task(conn, stuck, "r1.fastq", "r2.fastq", 1, 1)
        db.update_task(conn, stuck, status="processing",
                       stage=config.STAGE_PARSE)
        conn.close()
        resp = await post_review(client, stuck, valid_body)
        record("review non-succeeded task -> 409",
               resp.status_code == 409, resp.text)

        r2_path.unlink()
        resp = await post_review(client, tid, valid_body)
        record("missing result file -> 404 whole-request error",
               resp.status_code == 404, resp.text)
        # 恢复文件以验证只读语义（由后续成功请求侧面证明）
        r2_path.write_bytes(r2_before)

        # ---------- 10. FASTQ 后段损坏：整次报错 ---------- #
        corrupted = r1_before + b"@tail\nACGT\n+\n"  # 截断的末条记录
        r1_path.write_bytes(corrupted)
        resp = await post_review(client, tid, valid_body)
        record("corrupted FASTQ tail -> 400, no partial stats",
               resp.status_code == 400 and "sites" not in resp.json(),
               resp.text)
        r1_path.write_bytes(r1_before)  # 还原，确保后续检查基于原始文件

        # ---------- 11. 只读保证：成功复核不改动结果文件 ---------- #
        good_body = {
            "min_quality": 20,
            "sites": [{"position": 5, "ref": "A", "alt": "T"}],
            "alignments": [
                {"read_id": "rA", "end": 1, "ref_start": 1,
                 "strand": "forward", "cigar": "8M"},
                {"read_id": "rA", "end": 2, "ref_start": 1,
                 "strand": "reverse", "cigar": "8M"},
            ],
        }
        resp = await post_review(client, tid, good_body)
        record("review after restore works",
               resp.status_code == 200, resp.text)
        record("R1.fastq byte-identical after reviews",
               r1_path.read_bytes() == r1_before)
        record("R2.fastq byte-identical after reviews",
               r2_path.read_bytes() == r2_before)

        # ---------- 12. 超 500 碱基的被引用读段 ---------- #
        long_seq = b"A" * (config.REVIEW_MAX_READ_LENGTH + 1)
        long_qual = b"?" * len(long_seq)
        lr1 = [make_record("long", long_seq, long_qual)]
        lr2 = [make_record("long", b"ACGT", q("????"))]
        ltid = await upload_pair(client, lr1, lr2)
        resp = await post_review(
            client, ltid,
            {"min_quality": 10,
             "sites": [{"position": 1, "ref": "A", "alt": "C"}],
             "alignments": [
                 {"read_id": "long", "end": 1, "ref_start": 1,
                  "strand": "forward",
                  "cigar": f"{config.REVIEW_MAX_READ_LENGTH + 1}M"}]})
        record("referenced read over 500 bases -> 400",
               resp.status_code == 400, resp.text)

        # 未被引用的长读段不影响只引用 R2 的复核
        resp = await post_review(
            client, ltid,
            {"min_quality": 10,
             "sites": [{"position": 1, "ref": "A", "alt": "C"}],
             "alignments": [
                 {"read_id": "long", "end": 2, "ref_start": 1,
                  "strand": "forward", "cigar": "4M"}]})
        record("unreferenced long end does not block review",
               resp.status_code == 200
               and site_counter(resp.json(), 1)["ref_count"] == 1, resp.text)

    shutil.rmtree(tmp_root, ignore_errors=True)
    failed = [name for name, ok, _ in RESULTS if not ok]
    print("\n%d/%d checks passed" % (len(RESULTS) - len(failed), len(RESULTS)))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(asyncio.run(main_async()))
