"""Model-agnostic source provenance records."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any


@dataclass(frozen=True)
class SourceRef:
    """Reference connecting a model component to authoritative source code."""

    path: str
    symbol: str
    role: str
    repository: str = ""
    commit: str | None = None
    function: str | None = None
    units: str | None = None
    assumptions: tuple[str, ...] = ()
    approximation_level: str = "native-source"

    def __post_init__(self) -> None:
        for field_name in ("path", "symbol", "role"):
            if not getattr(self, field_name):
                raise ValueError(f"{field_name} must not be empty")
        if not self.approximation_level:
            raise ValueError("approximation_level must not be empty")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible representation."""
        result = asdict(self)
        result["assumptions"] = list(self.assumptions)
        return result


def provenance_manifest(
    *, model: str, approximation_level: str, sources: tuple[SourceRef, ...]
) -> dict[str, Any]:
    """Build a JSON-compatible provenance manifest for a model component."""
    if not model:
        raise ValueError("model must not be empty")
    if not approximation_level:
        raise ValueError("approximation_level must not be empty")
    return {
        "model": model,
        "approximation_level": approximation_level,
        "sources": [source.to_dict() for source in sources],
    }