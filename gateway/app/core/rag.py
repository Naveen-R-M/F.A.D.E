import re, difflib, textwrap
from pathlib import Path
from typing import List, Tuple

Passage = Tuple[str, str]  # (pid, passage)

class RAG:
    """Lightweight stdlib retrieval over project_context + sample_queries."""
    def __init__(self, base_dir: Path):
        self.base = base_dir
        self.corpus = self._load_corpus()

    def _read(self, p: Path) -> str:
        try:
            return p.read_text(encoding="utf-8", errors="ignore") if p.exists() else ""
        except Exception:
            return ""

    def _tokenize(self, text: str):
        return [t for t in re.findall(r"[A-Za-z0-9]+", text.lower()) if len(t) > 2]

    def _chunk(self, text: str, size: int = 700, overlap: int = 120):
        if not text: return []
        words = text.split(); out, i = [], 0
        while i < len(words):
            out.append(" ".join(words[i:i+size]))
            i += max(1, size - overlap)
        return out

    def _score(self, query: str, passage: str) -> float:
        qtok, ptok = set(self._tokenize(query)), set(self._tokenize(passage))
        jacc = len(qtok & ptok) / (len(qtok | ptok) + 1e-9)
        fuzz = difflib.SequenceMatcher(None, query.lower(), passage.lower()).ratio()
        return 0.6 * jacc + 0.4 * fuzz

    def _parse_samples(self, md: str):
        labeled, unlabeled, current = [], [], None
        for raw in md.splitlines():
            line = raw.strip()
            if not line: continue
            if re.match(r"^##\s*CHAT QUERIES", line, re.I): current = "CHAT"; continue
            if re.match(r"^##\s*JOB QUERIES", line, re.I): current = "JOB"; continue
            if re.match(r"^[-*]\s+", line) or re.match(r"^\d+\.\s+", line):
                text = re.sub(r"^[-*]\s+|\d+\.\s+", "", line).strip()
                if text:
                    (labeled if current in {"CHAT","JOB"} else unlabeled).append((current, text) if current else text)
        return labeled, unlabeled

    def _load_corpus(self):
        project_text = self._read(self.base / "project_context.txt")
        samples_md = self._read(self.base / "sample_queries.txt")

        labeled, unlabeled = self._parse_samples(samples_md)
        project_chunks: List[Passage] = [(f"project_context.txt:{i}", c) for i, c in enumerate(self._chunk(project_text))]
        sample_chunks: List[Passage]  = [(f"sample_queries.txt:{i}", c) for i, c in enumerate([t for _, t in labeled] + unlabeled)]

        return {
            "project_text": project_text,
            "project_chunks": project_chunks,
            "sample_labeled": labeled,
            "sample_unlabeled": unlabeled,
            "sample_chunks": sample_chunks
        }

    def retrieve(self, query: str, k: int = 6):
        passages: List[Passage] = []
        passages += self.corpus.get("project_chunks", [])
        passages += self.corpus.get("sample_chunks", [])
        if not passages: return []
        scored = [(self._score(query, p), pid, p) for pid, p in passages]
        scored.sort(key=lambda x: x[0], reverse=True)
        hits = scored[:k]
        # Compactify for tool result
        return [
            {"pid": pid, "score": float(score),
             "snippet": textwrap.shorten(p.replace("\n"," "), 400)}
            for (score, pid, p) in hits
        ]
