"""
molecule_plotter.py
"""

import io
from typing import Optional

from rdkit.Chem import AllChem as rdkit_all_chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point2D
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
        drawer = self._make_drawer(width, height)
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

    def _make_drawer(self, width: int, height: int) -> rdMolDraw2D.MolDraw2DCairo:
        """Compute 2D coords, draw molecule, return finished drawer."""
        rdkit_all_chem.Compute2DCoords(self.mol)
        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
        drawer.drawOptions().addAtomIndices = True
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
        """Draw a text label with a white background box, centered at (cx, cy)."""
        lines = txt.split("\n")
        line_h = 15
        max_w = max(draw.textlength(l, font=font) for l in lines)
        total_h = line_h * len(lines)
        pad = 3
        draw.rectangle(
            [
                cx - max_w / 2 - pad,
                cy - total_h / 2 - pad,
                cx + max_w / 2 + pad,
                cy + total_h / 2 + pad,
            ],
            fill=(255, 255, 255, 210),
        )
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

    def _draw_charges(self, draw, drawer, font):
        """
        Partial charge label above each atom, color-coded red (+) / blue (-).
        """
        pmd_atoms = list(self.struct.atoms)
        for atom in self.mol.GetAtoms():
            idx = atom.GetIdx()
            charge = float(pmd_atoms[idx].charge)
            x, y = self._get_px(drawer, idx)
            color = (180, 0, 0, 255) if charge >= 0 else (0, 0, 180, 255)
            self._draw_label(draw, f"q={charge:+.3f}", x, y - 22, color, font)
