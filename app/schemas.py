"""Pydantic request/response models for the SNV review API.

结构性约束（数量、取值域、CIGAR 词法、互异性）在模型层校验，FastAPI 会
统一返回 422；需要读取结果数据才能判定的语义错误（未知读号、结果缺件、
FASTQ 后段损坏等）由引擎抛 ReviewError，路由层转为 400/404/409。
"""

from __future__ import annotations

import re
from typing import Literal

from pydantic import BaseModel, Field, field_validator, model_validator

from . import config

_CIGAR_RE = re.compile(r"^[1-9]\d*[MIDS]([1-9]\d*[MIDS])*$")
_ALLELES = frozenset("ACGT")


class SiteIn(BaseModel):
    position: int = Field(..., ge=1, description="从 1 开始的参考坐标")
    ref: str = Field(..., min_length=1, max_length=1, description="参考碱基（ACGT）")
    alt: str = Field(..., min_length=1, max_length=1, description="候选 ALT 碱基（ACGT）")

    @field_validator("ref", "alt")
    @classmethod
    def _acgt(cls, value: str) -> str:
        value = value.upper()
        if value not in _ALLELES:
            raise ValueError("only ACGT alleles are allowed")
        return value

    @model_validator(mode="after")
    def _distinct(self) -> "SiteIn":
        if self.ref == self.alt:
            raise ValueError("ref and alt must be different alleles")
        return self


class AlignmentIn(BaseModel):
    read_id: str = Field(..., min_length=1, max_length=1024, description="读段号")
    end: Literal[1, 2] = Field(..., description="端别：1=R1，2=R2")
    ref_start: int = Field(..., ge=1, description="比对在参考上的起点（1 开始）")
    strand: Literal["forward", "reverse"] = Field(
        ..., description="链向：forward 正向，reverse 反向"
    )
    cigar: str = Field(
        ..., min_length=1, max_length=4096,
        description="CIGAR 比对表，仅允许 M/I/D/S 操作且 S 只能在首尾",
    )

    @field_validator("cigar")
    @classmethod
    def _cigar_shape(cls, value: str) -> str:
        value = value.strip()
        if not _CIGAR_RE.fullmatch(value):
            raise ValueError("CIGAR must contain only positive-length M/I/D/S operations")
        tokens = re.findall(r"\d+([MIDS])", value)
        if "M" not in tokens:
            raise ValueError("CIGAR must contain at least one M operation")
        if "S" in tokens[1:-1]:
            raise ValueError("soft clipping (S) is only allowed at the edges")
        return value


class ReviewRequest(BaseModel):
    min_quality: int = Field(
        ..., ge=0, le=40, description="Phred+33 最低质量值（0～40）"
    )
    sites: list[SiteIn] = Field(
        ..., min_length=1, max_length=config.REVIEW_MAX_SITES,
        description=f"候选位点，至多 {config.REVIEW_MAX_SITES} 个，坐标互异",
    )
    alignments: list[AlignmentIn] = Field(
        ..., min_length=1, max_length=config.REVIEW_MAX_ALIGNMENTS,
        description=(
            "比对表；至多引用 "
            f"{config.REVIEW_MAX_PAIRS} 个互异读号，每端读长不超过 "
            f"{config.REVIEW_MAX_READ_LENGTH} 碱基"
        ),
    )

    @model_validator(mode="after")
    def _cross_checks(self) -> "ReviewRequest":
        positions = [site.position for site in self.sites]
        if len(set(positions)) != len(positions):
            raise ValueError("site positions must be unique")
        keys = [(item.read_id, item.end) for item in self.alignments]
        if len(set(keys)) != len(keys):
            raise ValueError("duplicate (read_id, end) entries in the alignment table")
        distinct_reads = {item.read_id for item in self.alignments}
        if len(distinct_reads) > config.REVIEW_MAX_PAIRS:
            raise ValueError(
                f"at most {config.REVIEW_MAX_PAIRS} distinct read ids are allowed"
            )
        return self


class SiteOut(BaseModel):
    position: int
    ref: str
    alt: str
    ref_count: int = Field(..., description="支持 REF 的读对数")
    alt_count: int = Field(..., description="支持 ALT 的读对数")
    other_count: int = Field(..., description="有效但既非 REF 也非 ALT 的读对数")
    conflict_count: int = Field(..., description="两端有效碱基相异的读对数")
    no_evidence_count: int = Field(..., description="无有效证据的读对数")


class ReviewResponse(BaseModel):
    task_id: str
    pairs_total: int = Field(..., description="本次复核覆盖的读对（互异读号）数")
    min_quality: int
    sites: list[SiteOut]
