"""
molecule_plotter.py
"""

import io
import re
from typing import Optional

from rdkit.Chem import AllChem as rdkit_all_chem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageDraw, ImageFont


class MoleculePlotter:
    """
    Handles all 2D parameter visualizations for a parameterized molecule.

    Requires:
        mol    - RDKit Mol object (from GenMolecule.mol)
        struct - ParmEd Structure object (from GenMolecule.struct)
        name   - molecule name string

    Usage:
        plotter = MoleculePlotter(mol=gen.mol, struct=gen.struct, name=gen.name)
        png = plotter.plot_2d_parameters(kind="bond", out_png="bonds.png")
    """

    FONT_PATH = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"

    KIND_COLOR = {
        "bond": (180, 0, 0, 255),  # red
        "angle": (0, 0, 180, 255),  # blue
        "dihedral": (0, 140, 0, 255),  # green
        "charge": None,  # dynamic: red=+, blue=-
    }

    def __init__(self, mol, struct, name: str):
        self.mol = mol
        self.struct = struct
        self.name = name

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def plot_2d_parameters(
        self,
        kind: str = "bond",
        width: int = 900,
        height: int = 700,
        out_png: Optional[str] = None,
    ) -> bytes:
        """
        Render a 2D structure diagram annotated with force-field parameters.

        Parameters
        ----------
        kind : str
            One of 'bond', 'angle', 'dihedral', 'charge'.
        width, height : int
            Output image dimensions in pixels.
        out_png : str, optional
            If given, write PNG bytes to this path.

        Returns
        -------
        bytes
            Raw PNG image bytes.
        """
        self._validate(kind)
        drawer = self._make_drawer(width, height, add_atom_indices=(kind != "charge"))
        img, draw = self._png_to_pil(drawer)
        font, font_small = self._load_fonts()

        def label(txt, cx, cy, color, small=False):
            self._draw_label(
                draw, txt, cx, cy, color, font=font_small if small else font
            )

        def px(idx):
            return self._get_px(drawer, idx)

        if kind == "bond":
            for b in self.struct.bonds:
                if not b.type:
                    continue
                i, j = b.atom1.idx, b.atom2.idx
                x1, y1 = px(i)
                x2, y2 = px(j)
                txt = f"k={float(b.type.k):.2f}\nr0={float(b.type.req):.3f}"
                label(txt, 0.5 * (x1 + x2), 0.5 * (y1 + y2), self.KIND_COLOR["bond"])

        elif kind == "angle":
            for a in self.struct.angles:
                if not a.type:
                    continue
                i, j, k = a.atom1.idx, a.atom2.idx, a.atom3.idx
                x1, y1 = px(i)
                x2, y2 = px(j)
                x3, y3 = px(k)
                txt = f"k={float(a.type.k):.2f}\nθ0={float(a.type.theteq):.2f}"
                label(
                    txt,
                    (x1 + x2 + x3) / 3,
                    (y1 + y2 + y3) / 3,
                    self.KIND_COLOR["angle"],
                )

        elif kind == "dihedral":
            self._draw_dihedrals(draw, drawer, font_small)

        elif kind == "charge":
            self._draw_atom_indices(draw, drawer, font_small)
            self._draw_charges(draw, drawer, font_small)

        png = self._pil_to_png(img)
        if out_png:
            with open(out_png, "wb") as f:
                f.write(png)
        return png

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _validate(self, kind: str):
        if self.mol is None:
            raise RuntimeError("mol is None — call generate_3d_structure first.")
        if self.struct is None:
            raise RuntimeError("struct is None — call create_pmd first.")
        valid = set(self.KIND_COLOR)
        if kind not in valid:
            raise ValueError(f"kind must be one of {sorted(valid)}, got {kind!r}")

    def _make_drawer(
        self, width: int, height: int, add_atom_indices: bool = True
    ) -> rdMolDraw2D.MolDraw2DCairo:
        """Compute 2D coords, draw molecule, return finished drawer."""
        rdkit_all_chem.Compute2DCoords(self.mol)
        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
        drawer.drawOptions().addAtomIndices = add_atom_indices
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, self.mol)
        drawer.FinishDrawing()
        return drawer

    def _png_to_pil(self, drawer) -> tuple[Image.Image, ImageDraw.ImageDraw]:
        """Convert drawer PNG bytes to a Pillow Image + Draw pair."""
        img = Image.open(io.BytesIO(drawer.GetDrawingText())).convert("RGBA")
        return img, ImageDraw.Draw(img)

    def _pil_to_png(self, img: Image.Image) -> bytes:
        buf = io.BytesIO()
        img.convert("RGB").save(buf, format="PNG")
        return buf.getvalue()

    def _load_fonts(self) -> tuple:
        try:
            return (
                ImageFont.truetype(self.FONT_PATH, 13),
                ImageFont.truetype(self.FONT_PATH, 11),
            )
        except Exception:
            f = ImageFont.load_default()
            return f, f

    @staticmethod
    def _get_px(drawer, idx: int) -> tuple[float, float]:
        p = drawer.GetDrawCoords(idx)
        return p.x, p.y

    @staticmethod
    def _draw_label(
        draw: ImageDraw.ImageDraw,
        txt: str,
        cx: float,
        cy: float,
        color: tuple,
        font,
    ):
        """Draw a text label centered at (cx, cy), no background."""
        lines = txt.split("\n")
        line_h = 15
        total_h = line_h * len(lines)
        for n, line in enumerate(lines):
            lw = draw.textlength(line, font=font)
            draw.text(
                (cx - lw / 2, cy - total_h / 2 + n * line_h),
                line,
                font=font,
                fill=color,
            )

    def _draw_dihedrals(self, draw, drawer, font):
        """
        One label per unique atom quartet, placed perpendicular to the
        central j-k bond. All Fourier terms for that quartet are stacked.
        """
        seen: set = set()
        all_dihedrals = list(self.struct.dihedrals)

        for d in all_dihedrals:
            if not d.type:
                continue
            i, j, k, l = d.atom1.idx, d.atom2.idx, d.atom3.idx, d.atom4.idx
            key = (i, j, k, l)
            if key in seen:
                continue
            seen.add(key)

            x2, y2 = self._get_px(drawer, j)
            x3, y3 = self._get_px(drawer, k)
            cx, cy = 0.5 * (x2 + x3), 0.5 * (y2 + y3)

            # Perpendicular offset so label doesn't sit on the bond
            dx, dy = x3 - x2, y3 - y2
            length = max((dx**2 + dy**2) ** 0.5, 1e-6)
            cx += (-dy / length) * 28
            cy += (dx / length) * 28

            # Collect all Fourier terms sharing this quartet
            terms = [
                dd
                for dd in all_dihedrals
                if dd.type
                and (dd.atom1.idx, dd.atom2.idx, dd.atom3.idx, dd.atom4.idx) == key
            ]
            txt = "\n".join(
                f"k={float(t.type.phi_k):.2f} n={int(t.type.per)} φ={float(t.type.phase):.0f}°"
                for t in terms
            )
            self._draw_label(draw, txt, cx, cy, self.KIND_COLOR["dihedral"], font)

    def _draw_atom_indices(self, draw, drawer, font):
        """Draw atom indices offset from atom centers to avoid label overlap."""
        for atom in self.mol.GetAtoms():
            idx = atom.GetIdx()
            x, y = self._get_px(drawer, idx)
            draw.text(
                (x - 18, y - 18),
                str(idx),
                font=font,
                fill=(90, 90, 90, 255),
            )

    def _draw_charges(self, draw, drawer, font):
        """
        Partial charge label above each atom, color-coded red (+) / blue (-).
        """
        pmd_atoms = list(self.struct.atoms)  # Ensure we have the correct atom structure
        for atom in self.mol.GetAtoms():
            idx = atom.GetIdx()
            charge = float(pmd_atoms[idx].charge)
            x, y = self._get_px(drawer, idx)
            color = (180, 0, 0, 255) if charge >= 0 else (0, 0, 180, 255)
            self._draw_label(
                draw,
                f"q={charge:+.3f}",
                x,
                y + 16,
                color,
                font,
            )

    @staticmethod
    def plot_mol_file(
        path: str,
        type_labels: Optional[dict] = None,
        width: int = 900,
        height: int = 700,
        out_png: Optional[str] = None,
    ) -> bytes:
        """
        Parse a LAMMPS molecule template file and render a 2D diagram
        annotated with atom type labels and partial charges.

        Parameters
        ----------
        path : str
            Path to the LAMMPS .mol file.
        type_labels : dict, optional
            Mapping of integer type ID -> label string, e.g. {1: 'ca', 2: 'nv'}.
            If provided, labels are shown instead of integer IDs.
            Can also be a string->string mapping (for files that already use labels).
        width, height : int
            Output image size in pixels.
        out_png : str, optional
            If given, write the PNG to this path.

        Returns
        -------
        bytes
            Raw PNG image bytes.
        """
        atoms, bonds, charges, types = MoleculePlotter._parse_mol_file(path)

        # Build an RDKit RWMol from connectivity only (no element info needed)
        rwmol = rdkit_all_chem.RWMol()
        for _ in atoms:
            rwmol.AddAtom(rdkit_all_chem.Atom(6))  # placeholder carbon — layout only
        for a1, a2 in bonds:
            rwmol.AddBond(a1 - 1, a2 - 1, rdkit_all_chem.rdchem.BondType.SINGLE)
        mol = rwmol.GetMol()
        rdkit_all_chem.Compute2DCoords(mol)

        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
        drawer.drawOptions().addAtomIndices = False
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
        drawer.FinishDrawing()

        img = Image.open(io.BytesIO(drawer.GetDrawingText())).convert("RGBA")
        draw = ImageDraw.Draw(img)

        try:
            font_type = ImageFont.truetype(
                "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 15
            )
            font_charge = ImageFont.truetype(
                "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 13
            )
            font_idx = ImageFont.truetype(
                "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 12
            )
        except Exception:
            fallback = ImageFont.load_default()
            font_type = fallback
            font_charge = fallback
            font_idx = fallback

        tl = type_labels or {}
        occupied_boxes: list[tuple[float, float, float, float]] = []

        for atom_id in atoms:
            idx = atom_id - 1
            p = drawer.GetDrawCoords(idx)
            x, y = p.x, p.y

            raw_type = types.get(atom_id, "?")
            label_str = str(tl.get(raw_type, tl.get(str(raw_type), raw_type)))
            charge = charges.get(atom_id, 0.0)

            # Type label near atom center; shifted if collision occurs.
            MoleculePlotter._place_label_static(
                draw,
                label_str,
                x,
                y,
                (60, 60, 60, 255),
                font_type,
                occupied_boxes,
                candidates=[
                    (0, -16),
                    (18, -10),
                    (-18, -10),
                    (22, 8),
                    (-22, 8),
                    (0, 18),
                    (26, -24),
                    (-26, -24),
                ],
            )

            # Charge label, prefer below, then fan out if crowded.
            c_color = (180, 0, 0, 255) if charge >= 0 else (0, 0, 180, 255)
            MoleculePlotter._place_label_static(
                draw,
                f"q={charge:+.3f}",
                x,
                y,
                c_color,
                font_charge,
                occupied_boxes,
                candidates=[
                    (0, 18),
                    (22, 18),
                    (-22, 18),
                    (0, 30),
                    (24, 30),
                    (-24, 30),
                    (30, 4),
                    (-30, 4),
                ],
            )

        # Atom index labels, prefer top-left but avoid existing labels.
        for atom_id in atoms:
            idx = atom_id - 1
            p = drawer.GetDrawCoords(idx)
            MoleculePlotter._place_label_static(
                draw,
                str(atom_id),
                p.x,
                p.y,
                (100, 100, 100, 220),
                font_idx,
                occupied_boxes,
                candidates=[
                    (-18, -20),
                    (16, -20),
                    (-18, 20),
                    (16, 20),
                    (0, -28),
                    (0, 28),
                ],
            )

        png_buf = io.BytesIO()
        img.convert("RGB").save(png_buf, format="PNG")
        png = png_buf.getvalue()

        if out_png:
            with open(out_png, "wb") as fh:
                fh.write(png)
        return png

    @staticmethod
    def _parse_mol_file(path: str) -> tuple[list, list, dict, dict]:
        """
        Parse a LAMMPS molecule template file.

        Returns
        -------
        atoms   : list[int]           — 1-based atom IDs in order
        bonds   : list[tuple[int,int]]— (atom1, atom2) pairs
        charges : dict[int, float]    — atom_id -> charge
        types   : dict[int, str/int]  — atom_id -> type (int or str)
        """
        section = None
        atoms: list[int] = []
        bonds: list[tuple] = []
        charges: dict[int, float] = {}
        types: dict = {}

        with open(path, encoding="utf-8") as fh:
            for raw in fh:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue

                # Section headers (case-insensitive)
                low = line.lower()
                if low == "types":
                    section = "types"
                    continue
                if low == "charges":
                    section = "charges"
                    continue
                if low == "coords":
                    section = "coords"
                    continue
                if low == "bonds":
                    section = "bonds"
                    continue
                if re.match(r"^(angles|dihedrals|impropers)", low):
                    section = None
                    continue
                # header counts line e.g. "20 atoms"
                if re.search(r"\batoms\b", low) and re.match(r"^\d+", line):
                    section = None
                    continue
                if re.search(
                    r"\b(bonds|angles|dihedrals|impropers|types)\b", low
                ) and re.match(r"^\d+", line):
                    section = None
                    continue

                parts = line.split()
                if not parts:
                    continue

                try:
                    atom_id = int(parts[0])
                except ValueError:
                    continue

                if section == "types" and len(parts) >= 2:
                    # type may be int or str label
                    try:
                        types[atom_id] = int(parts[1])
                    except ValueError:
                        types[atom_id] = parts[1]
                    if atom_id not in [a for a in atoms]:
                        atoms.append(atom_id)

                elif section == "charges" and len(parts) >= 2:
                    charges[atom_id] = float(parts[1])

                elif section == "coords" and len(parts) >= 4:
                    if atom_id not in atoms:
                        atoms.append(atom_id)

                elif section == "bonds" and len(parts) >= 3:
                    # format: bond_id  type  atom1  atom2
                    # OR:     atom1  atom2  (if no type column)
                    if len(parts) >= 4:
                        bonds.append((int(parts[2]), int(parts[3])))
                    else:
                        bonds.append((int(parts[1]), int(parts[2])))

        if not atoms:
            atoms = sorted(types.keys())

        return atoms, bonds, charges, types

    @staticmethod
    def _draw_label_static(draw, txt, cx, cy, color, font):
        """Draw a centered text label, no background."""
        lines = txt.split("\n")
        line_h = 14
        total_h = line_h * len(lines)
        for n, line in enumerate(lines):
            lw = draw.textlength(line, font=font)
            draw.text(
                (cx - lw / 2, cy - total_h / 2 + n * line_h),
                line,
                font=font,
                fill=color,
            )

    @staticmethod
    def _place_label_static(
        draw,
        txt,
        anchor_x,
        anchor_y,
        color,
        font,
        occupied_boxes,
        candidates,
    ):
        """Place a centered text label at the first non-overlapping candidate offset."""
        for dx, dy in candidates:
            cx = anchor_x + dx
            cy = anchor_y + dy
            left, top, right, bottom = draw.textbbox((0, 0), txt, font=font)
            w = right - left
            h = bottom - top
            box = (cx - w / 2 - 1, cy - h / 2 - 1, cx + w / 2 + 1, cy + h / 2 + 1)
            if not any(MoleculePlotter._boxes_overlap(box, b) for b in occupied_boxes):
                MoleculePlotter._draw_label_static(draw, txt, cx, cy, color, font)
                occupied_boxes.append(box)
                return

        dx, dy = candidates[0]
        cx = anchor_x + dx
        cy = anchor_y + dy
        MoleculePlotter._draw_label_static(draw, txt, cx, cy, color, font)
        left, top, right, bottom = draw.textbbox((0, 0), txt, font=font)
        w = right - left
        h = bottom - top
        occupied_boxes.append(
            (cx - w / 2 - 1, cy - h / 2 - 1, cx + w / 2 + 1, cy + h / 2 + 1)
        )

    @staticmethod
    def _boxes_overlap(a, b) -> bool:
        return not (a[2] <= b[0] or a[0] >= b[2] or a[3] <= b[1] or a[1] >= b[3])
