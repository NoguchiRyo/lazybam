from typing import Dict, Any, List, Tuple
from .error import DuplicateProgramIDError
import warnings


class BamHeader:
    """
    Manage a BAM header: parse from SAM-style header text bytes,
    edit program/reference entries, and serialize back to header text bytes.
    """

    def __init__(
        self, header: Dict[str, Any], refs: List[Tuple[str, int]] | None = None
    ):
        # header: SAM-style parsed dict, refs: list of (name, length)
        self.header = header
        self.refs = refs if refs is not None else []

    def get_ref_name(self, rid: int) -> str:
        """Get reference name by index."""
        if rid < 0 or rid >= len(self.refs):
            raise IndexError(f"Reference ID {rid} out of range")
        return self.refs[rid][0]

    @classmethod
    def parse_header_text(
        cls, text: str
    ) -> tuple[Dict[str, Any], list[tuple[str, int]]]:
        """Parse SAM-style header text into a dict structure."""
        h: Dict[str, Any] = {}
        refs: List[Tuple[str, int]] = []
        # Parse header lines
        for line in text.strip().split("\n"):
            if not line or not line.startswith("@"):
                continue
            tag, *fields = line.split("\t")
            if tag == "@HD":
                h[tag] = {}
                for f in fields:
                    key, value = f.split(":", 1)
                    h[tag][key] = value
            elif tag in ("@SQ", "@RG", "@PG"):
                rec: Dict[str, Any] = {}
                for f in fields:
                    key, value = f.split(":", 1)
                    rec[key] = value
                h.setdefault(tag, []).append(rec)
                if tag == "@SQ":
                    # Add reference name and length to refs
                    if "SN" in rec and "LN" in rec:
                        refs.append((rec["SN"], int(rec["LN"])))
            elif tag == "@CO":
                comment = fields[0] if fields else ""
                h.setdefault(tag, []).append(comment)
        return h, refs

    @classmethod
    def from_bytes(cls, data: bytes) -> "BamHeader":
        """
        Construct BamHeader by decoding header bytes (SAM-style) and parsing text.
        """
        text = data.decode("ascii")
        header_dict, refs = cls.parse_header_text(text)
        return cls(header_dict, refs)

    def get_header_text(self) -> str:
        """Reconstruct SAM-style header text from the internal dict."""
        lines: List[str] = []
        # HD
        hd = self.header.get("@HD", {})
        fields = [f"{k}:{v}" for k, v in hd.items()]
        lines.append("@HD" + ("\t" + "\t".join(fields) if fields else ""))
        # SQ, RG, PG
        for tag in ("@SQ", "@RG", "@PG"):
            for rec in self.header.get(tag, []):
                fields = [f"{k}:{v}" for k, v in rec.items()]
                lines.append(tag + ("\t" + "\t".join(fields) if fields else ""))
        # CO
        for comment in self.header.get("@CO", []):
            lines.append("@CO" + ("\t" + comment if comment else ""))
        return "\n".join(lines) + "\n"

    def to_bytes(self) -> bytes:
        """
        Serialize header to SAM-style text bytes.
        """
        return self.get_header_text().encode("ascii")

    def to_string(self) -> str:
        """
        Serialize header to SAM-style text string.
        """
        return self.get_header_text()

    def add_program(self, **fields: str) -> None:
        """Add a @PG program record. Fields must include at least ID."""
        if "ID" not in fields:
            raise ValueError("Program record must include 'ID' field")

        # check if ID already exists
        try:
            for rec in self.header.get("@PG", []):
                if rec.get("ID") == fields["ID"]:
                    raise DuplicateProgramIDError()
            self.header.setdefault("@PG", []).append(fields)
        except DuplicateProgramIDError:
            warnings.warn(f"Program ID '{fields['ID']}' already exists")
            warnings.warn("Trying again with a different ID.")
            # Generate a new ID
            if str.isdigit(fields["ID"][-1]):
                new_id = fields["ID"][:-1] + str(int(fields["ID"][-1]) + 1)
            else:
                new_id = fields["ID"] + "2"
            fields["ID"] = new_id
            self.add_program(**fields)

    def change_SO_tag(self, value: str) -> None:
        """Change the SO tag in the header."""
        if "@HD" not in self.header:
            self.header["@HD"] = {
                "VN": "1.6",
                "SO": value,
            }
        else:
            self.header["@HD"]["SO"] = value

    def remove_program(self, program_id: str) -> None:
        """Remove a @PG record by its ID."""
        self.header["@PG"] = [
            r for r in self.header.get("@PG", []) if r.get("ID") != program_id
        ]

    def replace_references(self, refs: List[Tuple[str, int]]) -> None:
        """
        Replace all references with new list. Each tuple is (name, length).
        This updates both self.refs and '@SQ' entries in header.
        """
        self.refs = refs
        self.header["@SQ"] = [{"SN": name, "LN": str(length)} for name, length in refs]

    def __repr__(self) -> str:
        """Return a string representation of the header."""
        return self.to_string()

    def __str__(self) -> str:
        """Return a string representation of the header."""
        return self.to_string()
