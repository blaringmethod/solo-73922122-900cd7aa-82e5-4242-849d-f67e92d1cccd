# Paired FASTQ Re-pairing Service

测序湖入湖工具：接收顺序被打乱的双端未压缩 FASTQ（R1/R2），按读段号重新
配对，按 R1 顺序交付两份文件，原文逐字节保留。FastAPI + SQLite 实现，单文件
上限 256 MiB，进程内存上限 256 MiB，可处理总量超过内存的输入。

## 配对规则

- 读段号取头行首个空白前的字段，去掉开头 `@`，再去掉末尾的 `/1` 或 `/2`
  （例如 `@read-7/1 sample=x` → `read-7`）。
- 输出按 **R1 的原始顺序**排列；R2 按相同读段号对齐到同一位置。
- 四行记录原文（含头行描述、换行风格）原样写入结果，不做任何改写。

## 失败即整次失败，不开放部分结果

任一情况发生，任务整体标记为 `failed`/`cancelled`，已写入的暂存全部删除，
下载接口返回 409：

- 同一文件内读段号重复（含 `/1`、`/2` 归一化后重复）
- 孤端（某读段号只出现在一侧）
- 格式错误：非四行式、头行不以 `@` 开头、第三行不以 `+` 开头、空文件
- 序列与质量字符串长度不等

## 内存与大输入

- 上传经 `SpooledTemporaryFile` 落盘（超过 1 MiB 即写临时文件），分块复制，
  单文件流式校验 256 MiB 上限，超限返回 413。
- 解析按 512 条/批写入 SQLite（WAL、小 page cache、临时表走文件
  `temp_store=FILE`、关闭 mmap），任意时刻内存里只有一个批次。
- 配对发布使用 `JOIN ... ORDER BY a.ord` 流式游标，分批 `fetchmany` 写出。
- 因此输入总量可以远超 256 MiB；自测以约 393 MiB 总量（两份各 196.5 MiB、
  各 200 万条）验证，峰值 RSS 约 19 MiB（见下文自测结果）。

## 任务与取消语义

- `POST /tasks` 上传后立即返回 `202`，后台守护线程处理。
- `GET /tasks/{id}` 查询状态、阶段（upload/parse/validate/publish/done）与
  0–1 进度。
- `POST /tasks/{id}/cancel`：置取消标记后**立即返回**，不等待配对结束；
  工作线程在下一个批次边界中止。取消与完成之间用条件 UPDATE 仲裁：
  - 取消先到 → 任务最终 `cancelled`，结果永不发布，暂存目录删除；
  - 完成先到 → 取消请求得到 `409`，明确告知任务已完成、不可取消，结果保留。
- 取消后可立刻创建新任务；任务目录以 UUID 隔离，清理只删除本任务目录，
  不会误删新任务文件。
- 服务重启（启动钩子）：所有非终态任务标记 `failed`（原因
  `interrupted by service restart`）并清理暂存与记录；`succeeded` 任务的
  结果文件保留，仍可下载。

## API

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| POST | `/tasks` | multipart 字段 `r1`、`r2`，返回任务对象（202） |
| GET  | `/tasks/{id}` | 状态/进度/消息/下载链接 |
| POST | `/tasks/{id}/cancel` | 请求取消（立即返回） |
| GET  | `/tasks/{id}/download/r1` | 下载配对后的 R1 |
| GET  | `/tasks/{id}/download/r2` | 下载配对后的 R2 |
| POST | `/tasks/{id}/reviews` | 候选 SNV 双端读段证据复核（见下节） |
| GET  | `/health` | 健康检查 |

启动后打开交互式文档：<http://localhost:8000/docs>，可直接在页面上演示
上传、取消、下载全流程；复核接口的请求体 schema 与示例同样可在 `/docs`
中直接调用。

## SNV 复核接口 `POST /tasks/{id}/reviews`

变异复核员据此判断候选 SNV 是否有可信的双端读段支持。仅接受
**已成功（succeeded）** 的任务号；请求体包含：

- `min_quality`：Phred+33 最低质量值，整数 0～40；
- `sites`：1～8 个位点，每个含从 1 开始的参考坐标 `position`、互异的
  ACGT `ref` 与 `alt`，位点坐标必须互异；
- `alignments`：比对表，每条含 `read_id`（读号）、`end`（端别 1/2）、
  `ref_start`（1 开始的参考起点）、`strand`（`forward`/`reverse`）、
  `cigar`（仅 M/I/D/S，S 只允许在首尾，至少含一个 M）。

处理上限（超出整次 422，不处理）：

- 比对表至多引用 **2000 个互异读号**（即 2000 对读段），同一
  `(read_id, end)` 不得重复；
- 被引用读段**每端不超过 500 碱基**（CIGAR 读消耗必须与读长一致）。

证据提取与计数规则：

- 按比对位置在两端结果 FASTQ 上取证据；`reverse` 链先对序列**反向互补**、
  对质量串**整体反转**，再按 CIGAR 同步取碱基；
- CIGAR 语义：M 同时推进读与参考、I 仅推进读、D 仅推进参考、S 仅推进读
  且不覆盖任何参考位点；
- 只有**覆盖目标位点且质量达标**的碱基才计入；质量不达标或未覆盖视为
  该端无证据；
- 同一对读段在每个位点至多贡献一次：两端都有效且碱基一致只计一次，
  两端都有效却相异计一次**冲突**（不选一端、不重复计数），仅一端有效按
  该端碱基分类，两端均无有效证据计入 `no_evidence_count`；
- 碱基分类：等于 REF 计 `ref_count`，等于 ALT 计 `alt_count`，其余有效
  碱基（含 N、小写归一后仍非 REF/ALT 的碱基）计 `other_count`。

响应逐位点返回五个桶的读对数，五者之和恒等于 `pairs_total`：

```json
{
  "task_id": "...",
  "pairs_total": 2,
  "min_quality": 20,
  "sites": [{
    "position": 1, "ref": "A", "alt": "T",
    "ref_count": 0, "alt_count": 1, "other_count": 0,
    "conflict_count": 1, "no_evidence_count": 0
  }]
}
```

### 整次失败，不给部分统计、不改原文件

下列情况整次请求报错，响应中不含任何位点统计；结果 FASTQ 始终以只读
方式打开，复核前后逐字节不变：

- 任务不存在 → 404；任务未成功 → 409；结果文件缺件 → 404；
- 未知读号、重复 `(read_id, end)`、非法 CIGAR、CIGAR 与读长不符、
  引用读段超过 500 碱基、比对表超过 2000 个互异读号 → 400/422；
- 结果 FASTQ 结构损坏（含**后段截断**）、质量字节超出 Phred+33
  （ASCII 33～126）范围 → 400。引擎会完整流式扫描两份结果文件后才出
  统计，因此后段损坏同样能被捕获。

### 复核的内存边界

- 扫描逐记录流式进行，内存里只保留被引用读号（≤2000）的序列/质量与
  位点级计数（≤8 位点），不加载整份 FASTQ；
- 自测以打满上限的 2000 对 × 每端 500 碱基（双端各 4000 条比对、8 位点）
  实测子进程峰值 RSS 约 21 MiB，远低于 256 MiB 上限。

## 运行（Docker Compose，内存上限 256 MiB）

```bash
docker compose up --build
# 打开 http://localhost:8000/docs
```

`docker-compose.yml` 中 `mem_limit: 256m` 即交付要求的内存上限；数据保存在
命名卷 `fastq-data`（容器内 `/data`，含 SQLite 与各任务目录），重启后成功
结果仍可下载。

## 本地直接运行（开发/自测）

```bash
python3 -m venv --without-pip .venv
.venv/bin/python tools/get-pip.py
.venv/bin/pip install -r requirements.txt
FASTQ_DATA_DIR=./data .venv/bin/uvicorn app.main:app --port 8000
```

## 演示脚本

```bash
./scripts/demo.sh                 # 需要一个运行在 localhost:8000 的服务
BASE_URL=http://localhost:8000 ./scripts/demo.sh
```

脚本会生成一对乱序 FASTQ，演示上传 → 轮询进度 → 下载校验，并额外演示一次
取消与一次错误输入。

## 自测

```bash
# 短测试（功能/取消/重启语义 + SNV 复核 + 取消与完成的仲裁）
bash scripts/run_tests.sh
# 或分别运行
.venv/bin/python -I tests/test_api.py          # 30 项：配对、整次失败、取消、重启恢复
.venv/bin/python -I tests/test_review.py       # 42 项：SNV 复核证据规则与整次失败
.venv/bin/python -I tests/test_cancel_race.py  # 9 项：取消/完成条件 UPDATE 仲裁

# 内存压力测试（一）：配对路径，两份各 196.5 MiB（合计约 393 MiB，超过
# 256 MiB 内存上限），子进程实测峰值 RSS 约 19 MiB，并逐条校验配对顺序
.venv/bin/python -I tests/test_memory.py 2000000

# 内存压力测试（二）：复核路径，2000 对 × 每端 500 碱基的上限输入，
# 4000 条比对 + 8 位点，子进程实测峰值 RSS 约 21 MiB，逐位点核对计数总和
.venv/bin/python -I tests/test_review_memory.py
```

## 目录结构

```
app/
  config.py    # 上限、路径、批大小等配置（环境变量可覆盖，含复核上限）
  db.py        # SQLite schema 与任务/记录访问
  fastq.py     # 流式解析、校验、JOIN 配对发布
  worker.py    # 后台线程、取消仲裁、启动恢复、任务目录清理
  schemas.py   # 复核接口的 Pydantic 请求/响应模型与结构性校验
  review.py    # 复核引擎：CIGAR 取证据、反向互补、逐位点读对计数
  main.py      # FastAPI 路由（上传/状态/取消/下载/复核）与有界流式上传
tests/         # 功能自测与内存压力测试（含复核专项）
scripts/demo.sh
Dockerfile, docker-compose.yml, requirements.txt
```

## 遗留项 / 边界说明

- 仅接受未压缩 FASTQ；`.gz` 等压缩格式不在范围内。
- 服务为单实例 + 单 SQLite 文件设计；多实例横向扩展需要换用共享存储/队列。
- 没有内置任务 TTL 清理，成功结果会一直保留在数据卷中，需外部运维定期清理。
- 开发沙箱内未安装 Docker，`docker compose up` 未实跑；compose/Dockerfile 仅做
  了 YAML 与结构校验，应用本身已用等价的本地 uvicorn 命令完整验证（含
  `/docs`、上传、取消、下载、重启恢复）。在装有 Docker 的环境执行
  `docker compose up --build` 即可，`mem_limit: 256m` 即内存上限。
