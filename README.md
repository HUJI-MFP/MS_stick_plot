# MS analysis code

Desktop Python tools for mass-spectrometry results of depsipeptide / polyester oligomer libraries (Moran's lab). The tools read MS peak lists from Excel, find which oligomers of the chosen monomers explain each peak, and either draw an editable stick plot or quantify the amino-acid content.

| File | What it does |
| --- | --- |
| `stick plot MS.py` | Interactive MS stick plot with automatic oligomer labels (main tool) |
| `amino acid quantification MS.py` | Amino-acid percentage of the identified oligomers, written to Excel |
| `plot actions.txt` | Keyboard shortcuts of the stick plot (same as the table below) |
| `MS_Stick_Plot.spec` | PyInstaller recipe for building `MS_Stick_Plot.exe` |
| `old codes/stick plot MS V1.py` | First version of the stick plot, kept for reference |

## Input

An Excel workbook with one sheet per sample, each holding a mass (*m/z*) column and an intensity column (% intensity preferred). The sample name is read from the sheet when it is there, otherwise the sheet name is used.

The ionisation mode is taken from the file name (`positive` / `negative`) and can be overridden when asked. Masses are matched as [M+H]⁺ in positive mode and [M−H]⁻ in negative mode (±1.00784 Da), within 0.5 Da.

Built-in monomers: the 20 standard amino acids plus ornithine, Dab, Dap, β-aminobutyric acid, β-homo-arginine and β-homo-lysine, and the hydroxy acids glycolic, lactic, phenyllactic, hydroxypropionic and hydroxybutyric acid. Labels use 3-letter codes, which you can rename in the plot (key **S**).

## Stick plot (`stick plot MS.py`)

1. Pick the Excel file.
2. Choose the mode:
   - **binary**: one amino acid + one hydroxy acid per sample;
   - **multi-monomer**: any set of monomers, with an optional custom adduct (e.g. +Na, −CO₂).
3. Analyse a single sheet or all sheets (batch).
4. Set the minimum intensity threshold (default 5).

Every peak above the threshold is drawn:
- **red**: amino-acid / depsipeptide oligomers;
- **green**: hydroxy-acid / polyester oligomers;
- **co-polymer** colour: mixed oligomers;
- **blue**: unidentified peaks.

Noncovalent dimers (2M+H / 2M−H) are found automatically. Use the matplotlib toolbar to zoom, pan and save (PNG, SVG, PDF, …).

### Plot actions (keyboard)

Click the plot window first, then press:

| Key | Action | Notes |
| --- | --- | --- |
| E | Edit the main plot title | Updates immediately |
| X | Edit the X-axis label | Updates immediately |
| Y | Edit the Y-axis label | Updates immediately |
| T | Edit axis/tick text style | X/Y title font size + bold, X/Y tick font size + bold |
| Z | Edit X/Y axis limits | Set min/max for both axes |
| B | Toggle top/right border (spines) | Hide/show the top and right outline |
| L | Legend editor | Rename/delete/add items, show/hide, move, title, font size |
| G | Grid menu | Show / remove / toggle grid |
| N | Remove grid (clean mode) | No major/minor grid, white background |
| F | Edit figure size | Width/height in inches |
| P | Peak label font options | Increase/decrease or set the size of all peak labels |
| S | Synonym editor | Change monomer abbreviations; all labels update live |
| C | Category colour editor | Amino/hydroxy/co-polymer colours (labels, sticks, connectors, legend) |
| A | Add a custom label | Press A, then click near a peak |
| R | Remove a label | Press R, then click near a label |
| Mouse drag | Move a label | Left-click and drag the label text |

## Amino-acid quantification (`amino acid quantification MS.py`)

For each sample it identifies the oligomers (amino-only, hydroxy-only and co-polymers) and weights each one's amino-acid count by its abundance:

contribution = n(amino acids) × abundance %

It writes one Excel file per sample with two sheets:
- **Oligomers**: formula, matched *m/z*, abundance %, amino-acid count and normalised contribution;
- **Summary**: number of oligomers identified and the total normalised amino-acid contribution.

## Running

Python 3.10+ with:

```
pip install pandas matplotlib openpyxl xlrd
```

`tkinter` is used for the file picker and dialogs (bundled with Python on Windows and macOS). Run it with:

```
python "stick plot MS.py"
python "amino acid quantification MS.py"
```

### Windows .exe

```
pip install pyinstaller
pyinstaller MS_Stick_Plot.spec
```

The executable is written to `dist/MS_Stick_Plot.exe`. Build output (`build/`, `dist/`) is not kept in the repository.
