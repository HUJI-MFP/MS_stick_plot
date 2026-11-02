"""
MS Stick Plot Analysis Tool

Features:
- Plots ALL mass peaks above the specified intensity threshold
- Color-codes identified oligomers/monomers (red for amino/depsipeptide, green for hydroxy/polyester)
- Automatically identifies non-covalent dimers (2M+H or 2M-H)
- Shows unidentified peaks in blue
- Fully interactive plot editing:
  * Pan/Zoom using matplotlib toolbar
  * Drag labels to reposition them
  * Press 'e' to edit plot title
  * Press 'x' to edit X-axis label
  * Press 'y' to edit Y-axis label
  * Press 'a' then click near a peak to add a custom label
  * Press 'r' then click near a label to remove it
  * Save plots in various formats (use toolbar save button)
"""
import os
# Ensure interactive matplotlib backend under EXE
os.environ.setdefault('MPLBACKEND', 'TkAgg')
try:
    import matplotlib
    matplotlib.use('TkAgg')
except Exception:
    pass
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from itertools import product
import shutil
import tempfile
from typing import Optional

# --- GUI prompt helpers (Tk first, console fallback) ---
def ask_yes_no_gui(title: str, question: str) -> Optional[bool]:
    try:
        import tkinter as tk
        from tkinter import messagebox
        root = tk.Tk()
        root.withdraw()
        root.attributes('-topmost', True)
        ans = messagebox.askyesno(title, question)
        root.destroy()
        return bool(ans)
    except Exception:
        return None

def prompt_yes_no(question: str, default: bool = False) -> bool:
    ans = ask_yes_no_gui("Question", question)
    if ans is not None:
        return ans
    try:
        reply = input(f"{question} [y/N]: ").strip().lower()
        return reply in ("y", "yes")
    except Exception:
        return default

def prompt_text(question: str, default: Optional[str] = None, title: str = "Input") -> Optional[str]:
    try:
        import tkinter as tk
        from tkinter import simpledialog
        root = tk.Tk()
        root.withdraw()
        root.attributes('-topmost', True)
        ans = simpledialog.askstring(title, question, initialvalue=default, parent=root)
        root.destroy()
        if ans is None:
            return default
        return ans
    except Exception:
        pass
    try:
        prompt = f"{question}"
        if default is not None:
            prompt += f" (default: {default})"
        prompt += ": "
        ans = input(prompt).strip()
        return ans if ans else default
    except Exception:
        return default

def prompt_float(question: str, default: Optional[float] = None, title: str = "Input") -> Optional[float]:
    s_default = None if default is None else str(default)
    try:
        import tkinter as tk
        from tkinter import simpledialog
        root = tk.Tk()
        root.withdraw()
        root.attributes('-topmost', True)
        ans = simpledialog.askstring(title, question, initialvalue=s_default, parent=root)
        root.destroy()
        if ans is None or ans.strip() == "":
            return default
        try:
            return float(ans)
        except ValueError:
            return default
    except Exception:
        pass
    try:
        ans = input(f"{question} (default: {s_default}): ").strip()
        return float(ans) if ans else default
    except Exception:
        return default

def choose_from_list(options, prompt):
    # Prefer a GUI list selection dialog
    sel = gui_select_from_list(options, title="Select", message=prompt)
    if sel is not None:
        return sel
    # Console fallback (rare in EXE)
    try:
        for i, item in enumerate(options, 1):
            print(f"{i}. {item}")
        idx = int(input(prompt)) - 1
        if 0 <= idx < len(options):
            return options[idx]
    except Exception:
        pass
    return options[0]

def gui_select_from_list(options, title: str, message: str) -> Optional[str]:
    try:
        import tkinter as tk
        from tkinter import ttk
        sel = {"value": None}
        root = tk.Tk()
        root.title(title)
        root.attributes('-topmost', True)
        root.resizable(False, False)
        frame = ttk.Frame(root, padding=10)
        frame.pack(fill='both', expand=True)
        ttk.Label(frame, text=message, wraplength=380).pack(anchor='w', pady=(0,8))
        lb = tk.Listbox(frame, height=min(12, max(4, len(options))), exportselection=False)
        for item in options:
            lb.insert('end', str(item))
        lb.selection_set(0)
        lb.pack(fill='both', expand=True)
        btns = ttk.Frame(frame)
        btns.pack(fill='x', pady=(8,0))
        def on_ok():
            try:
                idx = lb.curselection()
                if idx:
                    sel["value"] = options[idx[0]]
            except Exception:
                sel["value"] = None
            root.destroy()
        def on_cancel():
            sel["value"] = None
            root.destroy()
        ttk.Button(btns, text="OK", command=on_ok).pack(side='right')
        ttk.Button(btns, text="Cancel", command=on_cancel).pack(side='right', padx=(0,8))
        root.mainloop()
        return sel["value"]
    except Exception:
        return None

def gui_select_amino_hydroxy(amino_options, hydroxy_options, sample_label: str, default_amino: Optional[str] = None, default_hydroxy: Optional[str] = None):
    """Return (amino, hydroxy) or None if user clicks Skip/Cancel.
    Always forces an explicit selection to avoid carryover between sheets.
    """
    try:
        import tkinter as tk
        from tkinter import ttk, messagebox
        result = {"amino": None, "hydroxy": None, "skip": False}
        root = tk.Tk()
        root.title(f"Select for {sample_label}")
        root.attributes('-topmost', True)
        root.resizable(False, False)
        frm = ttk.Frame(root, padding=10)
        frm.pack(fill='both', expand=True)
        ttk.Label(frm, text=f"Sample: {sample_label}").grid(row=0, column=0, columnspan=2, sticky='w', pady=(0,8))
        ttk.Label(frm, text="Amino acid:").grid(row=1, column=0, sticky='w')
        ttk.Label(frm, text="Hydroxy acid:").grid(row=2, column=0, sticky='w', pady=(6,0))

        placeholder = "— Select —"
        amino_values = sorted(list(amino_options))
        hydroxy_values = sorted(list(hydroxy_options))
        amino_var = tk.StringVar()
        hydroxy_var = tk.StringVar()
        amino_cb = ttk.Combobox(frm, textvariable=amino_var, state='readonly', width=32, values=[placeholder] + amino_values)
        hydroxy_cb = ttk.Combobox(frm, textvariable=hydroxy_var, state='readonly', width=32, values=[placeholder] + hydroxy_values)
        amino_var.set(placeholder)
        hydroxy_var.set(placeholder)
        amino_cb.grid(row=1, column=1, sticky='w')
        hydroxy_cb.grid(row=2, column=1, sticky='w', pady=(6,0))

        btns = ttk.Frame(frm)
        btns.grid(row=3, column=0, columnspan=2, sticky='e', pady=(10,0))

        def ok():
            a = amino_var.get()
            h = hydroxy_var.get()
            if (a == placeholder) or (h == placeholder) or (a not in amino_values) or (h not in hydroxy_values):
                try:
                    messagebox.showwarning("Selection required", "Please select both amino and hydroxy acids.")
                except Exception:
                    pass
                return
            result["amino"], result["hydroxy"] = a, h
            root.destroy()

        def skip():
            result["skip"] = True
            root.destroy()

        ttk.Button(btns, text="Skip", command=skip).pack(side='right', padx=(0,8))
        ttk.Button(btns, text="OK", command=ok).pack(side='right')
        root.mainloop()
        if result["skip"]:
            return None
        return result["amino"], result["hydroxy"]
    except Exception:
        return None

# Minimal diagnostics to help when running as an EXE
try:
    _backend = getattr(__import__('matplotlib'), 'get_backend', lambda: 'unknown')()
    # Touch optional engines so PyInstaller includes them
    try:
        import tkinter  # noqa: F401
    except Exception:
        pass
    try:
        import openpyxl  # noqa: F401
    except Exception:
        pass
    try:
        import xlrd  # noqa: F401
    except Exception:
        pass
    print(f"[MS diag] backend={_backend}")
except Exception:
    pass

# --- Draggable Text with Connector Line ---
class DraggableTextWithLine:
    def __init__(self, text, line):
        self.text = text
        self.line = line
        self.press = None
        self.cid_press = text.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cid_release = text.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cid_motion = text.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_press(self, event):
        if event.inaxes != self.text.axes: return
        contains, _ = self.text.contains(event)
        if not contains: return
        x0, y0 = self.text.get_position()
        self.press = (x0, y0, event.xdata, event.ydata)

    def on_motion(self, event):
        if self.press is None or event.inaxes != self.text.axes: return
        x0, y0, xpress, ypress = self.press
        dx = event.xdata - xpress
        dy = event.ydata - ypress
        new_x, new_y = x0 + dx, y0 + dy
        self.text.set_position((new_x, new_y))
        self.line.set_data([self.line.get_xdata()[0], new_x], [self.line.get_ydata()[0], new_y])
        self.text.figure.canvas.draw()

    def on_release(self, event):
        self.press = None

# --- Interactive Plot Editor ---
class InteractivePlotEditor:
    """Allows interactive editing of plot elements with right-click context menu"""
    
    def __init__(self, fig, ax, draggables, filtered_data, labels_by_idx, ionization_mode, mass_adjustment, tolerance, stick_lines):
        self.fig = fig
        self.ax = ax
        self.draggables = draggables
        self.filtered_idx, self.filtered_mass, self.filtered_intensity = filtered_data
        self.labels_by_idx = labels_by_idx
        self.ionization_mode = ionization_mode
        self.mass_adjustment = mass_adjustment
        self.tolerance = tolerance
        self.stick_lines = stick_lines  # Dictionary mapping idx to LineCollection or vlines object
        
        # Connect keyboard events for editing
        self.cid_key = fig.canvas.mpl_connect('key_press_event', self.on_key_press)
        self.cid_click = fig.canvas.mpl_connect('button_press_event', self.on_click)
        
        # Initialize mode flags
        self.add_mode = False
        self.delete_mode = False
        
        print("\n" + "="*60)
        print("INTERACTIVE EDITING CONTROLS:")
        print("="*60)
        print("  'e'  - Edit plot title")
        print("  'x'  - Edit X-axis label")
        print("  'y'  - Edit Y-axis label")
        print("  'a'  - Add label to a peak (click near a peak)")
        print("       Stick color will change to match label color")
        print("  'r'  - Remove label (click near a label)")
        print("       Stick color will change back to blue (unidentified)")
        print("  Left-click and drag - Move existing labels")
        print("="*60)
        print("TIP: Make sure the plot window is focused (clicked) before pressing keys!")
        print("="*60 + "\n")
    
    def on_key_press(self, event):
        """Handle keyboard shortcuts for editing"""
        if event.key in ('e', 'E'):
            self.edit_title()
        elif event.key in ('x', 'X'):
            self.edit_xlabel()
        elif event.key in ('y', 'Y'):
            self.edit_ylabel()
        elif event.key in ('a', 'A'):
            print("ADD MODE: Click near a peak to add a label...")
            self.add_mode = True
        elif event.key in ('r', 'R'):
            print("REMOVE MODE: Click near a label to remove it...")
            self.delete_mode = True
    
    def on_click(self, event):
        """Handle mouse clicks for add/delete operations"""
        if event.inaxes != self.ax or event.button != 1:  # Left click only
            return
        
        if self.add_mode:
            self.add_label_at_position(event.xdata, event.ydata)
            self.add_mode = False
        elif self.delete_mode:
            self.delete_label_at_position(event.xdata, event.ydata)
            self.delete_mode = False
    
    def edit_title(self):
        """Edit the plot title"""
        try:
            import tkinter as tk
            from tkinter import simpledialog
            root = tk.Tk()
            root.withdraw()
            root.attributes('-topmost', True)
            current_title = self.ax.get_title()
            new_title = simpledialog.askstring("Edit Title", "Enter new plot title:", 
                                               initialvalue=current_title, parent=root)
            root.destroy()
            if new_title is not None:
                self.ax.set_title(new_title, fontsize=16, fontweight='bold')
                self.fig.canvas.draw()
        except Exception as e:
            print(f"Error editing title: {e}")
    
    def edit_xlabel(self):
        """Edit the X-axis label"""
        try:
            import tkinter as tk
            from tkinter import simpledialog
            root = tk.Tk()
            root.withdraw()
            root.attributes('-topmost', True)
            current_label = self.ax.get_xlabel()
            new_label = simpledialog.askstring("Edit X-axis Label", "Enter new X-axis label:", 
                                               initialvalue=current_label, parent=root)
            root.destroy()
            if new_label is not None:
                self.ax.set_xlabel(new_label, fontsize=14, fontweight='bold')
                self.fig.canvas.draw()
        except Exception as e:
            print(f"Error editing X-axis label: {e}")
    
    def edit_ylabel(self):
        """Edit the Y-axis label"""
        try:
            import tkinter as tk
            from tkinter import simpledialog
            root = tk.Tk()
            root.withdraw()
            root.attributes('-topmost', True)
            current_label = self.ax.get_ylabel()
            new_label = simpledialog.askstring("Edit Y-axis Label", "Enter new Y-axis label:", 
                                               initialvalue=current_label, parent=root)
            root.destroy()
            if new_label is not None:
                self.ax.set_ylabel(new_label, fontsize=14, fontweight='bold')
                self.fig.canvas.draw()
        except Exception as e:
            print(f"Error editing Y-axis label: {e}")
    
    def update_stick_color(self, idx, new_color):
        """Update the color of a stick line based on its label status"""
        if idx in self.stick_lines:
            stick = self.stick_lines[idx]
            try:
                stick.set_color(new_color)
            except Exception as e:
                print(f"Error updating stick color: {e}")
    
    def add_label_at_position(self, x, y):
        """Add a custom label at the clicked position"""
        # Find the nearest peak
        if not self.filtered_mass:
            print("No peaks to label!")
            return
        
        distances = [(abs(m - x), idx, m, i) for idx, m, i in zip(self.filtered_idx, self.filtered_mass, self.filtered_intensity)]
        distances.sort()
        _, nearest_idx, nearest_m, nearest_i = distances[0]
        
        try:
            import tkinter as tk
            from tkinter import simpledialog, colorchooser
            root = tk.Tk()
            root.withdraw()
            root.attributes('-topmost', True)
            
            # Ask for label text
            label_text = simpledialog.askstring("Add Label", 
                                                f"Enter label for peak at m/z={nearest_m:.2f}:", 
                                                parent=root)
            if not label_text:
                root.destroy()
                return
            
            # Ask for color
            color_choice = simpledialog.askstring("Label Color", 
                                                  "Enter color (red/green/blue/black) or leave blank for red:", 
                                                  initialvalue="red", parent=root)
            root.destroy()
            
            color = color_choice if color_choice else 'red'
            
            # Create the label
            txt = self.ax.text(nearest_m, nearest_i + 1, label_text, rotation=90, 
                             color=color, fontsize=9, fontweight='bold', 
                             ha='center', va='bottom')
            line = Line2D([nearest_m, nearest_m], [nearest_i, nearest_i + 1], 
                         color=color, linestyle='--', linewidth=0.8)
            self.ax.add_line(line)
            
            # Make it draggable
            draggable = DraggableTextWithLine(txt, line)
            self.draggables.append(draggable)
            
            # Update the labels dictionary
            self.labels_by_idx[nearest_idx] = (label_text, 'custom')
            
            # Update the stick color to match the label color
            self.update_stick_color(nearest_idx, color)
            
            self.fig.canvas.draw()
            
        except Exception as e:
            print(f"Error adding label: {e}")
    
    def delete_label_at_position(self, x, y):
        """Delete a label near the clicked position"""
        # Find labels near the click
        min_dist = float('inf')
        target_draggable = None
        target_idx = None
        
        for i, draggable in enumerate(self.draggables):
            label_x, label_y = draggable.text.get_position()
            dist = ((label_x - x)**2 + (label_y - y)**2)**0.5
            if dist < min_dist:
                min_dist = dist
                target_draggable = draggable
                target_idx = i
        
        # If close enough, delete it
        if target_draggable and min_dist < 50:  # threshold in data coordinates
            try:
                # Remove the text and line from the plot
                target_draggable.text.remove()
                target_draggable.line.remove()
                
                # Remove from draggables list
                self.draggables.pop(target_idx)
                
                # Find and remove from labels_by_idx, and update stick color
                label_pos_x = target_draggable.line.get_xdata()[0]
                removed_idx = None
                for idx, m in zip(self.filtered_idx, self.filtered_mass):
                    if abs(m - label_pos_x) < 0.01:
                        if idx in self.labels_by_idx:
                            del self.labels_by_idx[idx]
                            removed_idx = idx
                        break
                
                # Change stick color back to blue (unidentified)
                if removed_idx is not None:
                    self.update_stick_color(removed_idx, 'blue')
                
                self.fig.canvas.draw()
            except Exception as e:
                print(f"Error deleting label: {e}")

# --- Molecular Weights ---
molecular_weights = {
    "Glycine": 75.07, "Alanine": 89.09, "Valine": 117.15, "Leucine": 131.17,
    "Isoleucine": 131.17, "Proline": 115.13, "Phenylalanine": 165.19,
    "Tyrosine": 181.19, "Tryptophan": 204.23, "Serine": 105.09,
    "Threonine": 119.12, "Cysteine": 121.16, "Methionine": 149.21,
    "Asparagine": 132.12, "Glutamine": 146.15, "Aspartic Acid": 133.10,
    "Glutamic Acid": 147.13, "Lysine": 146.19, "Arginine": 174.20,
    "Histidine": 155.16, "Ornithine": 132.16, "Diaminobutyric acid": 118.13,
    "beta amino butyric acid": 103.12, "diamino propionic acid": 104.11,
    "Glycolic Acid": 76.05, "Lactic Acid": 90.08, "Phenyllactic Acid": 166.17,
    "beta homo-arginine": 188.13, "beta homo-lysine": 160.21, 
    "hydroxy propionic acid": 90.08, "hydroxybutyric acid": 104.10
}

# 3-letter codes for amino acids
amino_three_letter = {
    "Glycine": "Gly", "Alanine": "Ala", "Valine": "Val", "Leucine": "Leu",
    "Isoleucine": "Ile", "Proline": "Pro", "Phenylalanine": "Phe",
    "Tyrosine": "Tyr", "Tryptophan": "Trp", "Serine": "Ser",
    "Threonine": "Thr", "Cysteine": "Cys", "Methionine": "Met",
    "Asparagine": "Asn", "Glutamine": "Gln", "Aspartic Acid": "Asp",
    "Glutamic Acid": "Glu", "Lysine": "Lys", "Arginine": "Arg",
    "Histidine": "His", "Ornithine": "Orn", "Diaminobutyric acid": "Dab",
    "beta amino butyric acid": "bAba", "diamino propionic acid": "Dpr",
    "beta homo-arginine": "bhArg", "beta homo-lysine": "bhLys",
}

def choose_from_list(options, prompt):
    for i, item in enumerate(options, 1):
        print(f"{i}. {item}")
    idx = None
    while idx is None:
        try:
            idx = int(input(prompt)) - 1
            if not (0 <= idx < len(options)):
                print("Please enter a valid number from the list.")
                idx = None
        except ValueError:
            print("Please enter a number.")
    return options[idx]


def infer_ion_mode_from_filename(path: str) -> str:
    name = os.path.basename(path).lower()
    if "negative" in name:
        return "negative"
    if "positive" in name:
        return "positive"
    ans = prompt_text("Ionization mode not found in filename. Enter ionization mode (positive/negative)", default="positive", title="Ionization mode")
    ans = (ans or "positive").strip().lower()
    return ans if ans in ("positive", "negative") else "positive"


def get_mass_adjustment(mode: str) -> float:
    if mode == "positive":
        return 1.00784
    if mode == "negative":
        return -1.00784
    raise ValueError("Invalid ionization mode. Please choose 'positive' or 'negative'.")


def parse_ms_dataframe(df: pd.DataFrame):
    # Try common column names first (case-insensitive); prefer % intensity
    cols = list(df.columns)
    low_cols = [str(c).lower().strip() for c in cols]

    def find_contains(candidates):
        for pat in candidates:
            for orig, low in zip(cols, low_cols):
                if pat in low:
                    return orig
        return None

    # Mass column detection
    mass_col = find_contains(["m/z", "mz", "mass", "m over z", "m\u2215z", "m/z (da)"])

    # Intensity: prefer percentage-like columns
    percent_priority = ["%", "percent", "relative intensity (%)", "intensity (%)", "%intensity", "intensity %"]
    inten_col = find_contains(percent_priority)
    if inten_col is None:
        inten_col = find_contains(["relative intensity", "rel. intensity", "abundance", "intensity", "height"])

    # Identify numeric columns and compute heuristics
    numeric_cols = []
    numeric_series = {}
    for c in cols:
        s = pd.to_numeric(df[c], errors='coerce')
        if s.notna().sum() >= max(3, int(0.5 * len(s))):  # mostly numeric
            numeric_cols.append(c)
            numeric_series[c] = s

    # Prefer percent-like numeric column for intensity (0..100/110 range)
    percent_like = []
    for c in numeric_cols:
        s = numeric_series[c]
        mx = s.max(skipna=True)
        mn = s.min(skipna=True)
        if pd.notna(mx) and pd.notna(mn) and 0 <= mn <= mx <= 110:
            percent_like.append(c)

    # Heuristic for mass if still unknown
    if mass_col is None and numeric_cols:
        def mass_score(c):
            rng = numeric_series[c].max(skipna=True) - numeric_series[c].min(skipna=True)
            mx = numeric_series[c].max(skipna=True)
            # prefer columns with max > 150 (not percent-like) and large range
            return (mx if pd.notna(mx) else -1) + (rng if pd.notna(rng) else 0)
        mass_col = max(numeric_cols, key=mass_score)

    # Heuristic for intensity if still unknown
    if inten_col is None:
        if percent_like:
            # Choose the percent-like column with highest mean
            inten_col = max(percent_like, key=lambda c: numeric_series[c].mean(skipna=True))
        elif numeric_cols:
            # Fallback to a different numeric column than mass
            remain = [c for c in numeric_cols if c != mass_col]
            inten_col = remain[0] if remain else numeric_cols[0]

    # Last resort fallback
    if mass_col is None or inten_col is None:
        if len(cols) >= 2:
            mass_col = mass_col or cols[0]
            inten_col = inten_col or cols[1]
        else:
            raise ValueError("Could not identify mass and intensity columns.")

    # Convert to numeric
    mass = pd.to_numeric(df[mass_col], errors='coerce')
    intensity = pd.to_numeric(df[inten_col], errors='coerce')

    # Basic plausibility checks and auto-correction
    mass_max = mass.max(skipna=True)
    inten_max = intensity.max(skipna=True)

    # If mass looks percent-like and intensity looks mass-like, swap
    if pd.notna(mass_max) and pd.notna(inten_max):
        if (mass_max <= 150 and inten_max > 150) or (mass_max <= 110 and inten_max > 110):
            mass, intensity = intensity, mass
            mass_col, inten_col = inten_col, mass_col

    # Drop NaNs and sort by mass ascending to ensure proper x-axis
    valid = mass.notna() & intensity.notna()
    mass = mass[valid]
    intensity = intensity[valid]
    order = mass.sort_values().index
    mass = mass.loc[order].reset_index(drop=True)
    intensity = intensity.loc[order].reset_index(drop=True)

    # Optional: print chosen columns for transparency
    print(f"Using columns -> mass: '{mass_col}', intensity: '{inten_col}'")

    return mass, intensity


def load_ms_from_excel(xlsx_path: str, sheet_name: str):
    # First try label-block parsing (headers embedded in arbitrary cells)
    try:
        raw = pd.read_excel(xlsx_path, sheet_name=sheet_name, header=None)
    except PermissionError:
        # Copy to a temp file if locked and retry
        tmpdir = tempfile.mkdtemp(prefix="ms_tmp_")
        tmp_path = os.path.join(tmpdir, os.path.basename(xlsx_path))
        shutil.copy2(xlsx_path, tmp_path)
        try:
            raw = pd.read_excel(tmp_path, sheet_name=sheet_name, header=None)
        finally:
            # cleanup temp copy
            shutil.rmtree(tmpdir, ignore_errors=True)

    def try_parse_label_blocks(raw_df: pd.DataFrame):
        max_rows = min(500, len(raw_df))
        max_cols = min(200, len(raw_df.columns))

        def clean(v):
            if pd.isna(v):
                return None
            s = str(v).strip()
            return s if s and s.lower() != 'nan' else None

        def is_intensity_header(text: str) -> bool:
            t = text.lower()
            return ('abundance' in t and '%' in text) or ('%' in text) or ('intensity' in t)

        # find m/z cell
        r_mz = c_mz = None
        for r in range(max_rows):
            for c in range(max_cols):
                cell = clean(raw_df.iloc[r, c])
                if cell and ('m/z' in cell.lower() or cell.lower() == 'mz' or 'm\u2215z' in cell.lower()):
                    r_mz, c_mz = r, c
                    break
            if r_mz is not None:
                break

        if r_mz is None:
            return None

        # find intensity header: "abundance %" preferably two cells to the right
        r_int = r_mz
        c_int = None
        candidate = clean(raw_df.iloc[r_mz, c_mz + 2]) if c_mz + 2 < max_cols else None
        if candidate and is_intensity_header(candidate):
            c_int = c_mz + 2
        else:
            # search the same row for a header containing 'abundance' and/or '%'
            for c in range(c_mz + 1, max_cols):
                cell = clean(raw_df.iloc[r_mz, c])
                if cell and is_intensity_header(cell):
                    c_int = c
                    break
        if c_int is None:
            # broader search anywhere
            found = False
            for rr in range(max_rows):
                for cc in range(max_cols):
                    cell = clean(raw_df.iloc[rr, cc])
                    if cell and is_intensity_header(cell):
                        r_int, c_int = rr, cc
                        found = True
                        break
                if found:
                    break

        if c_int is None:
            return None

        # Extract columns below headers
        mass_series = pd.to_numeric(raw_df.iloc[r_mz + 1:, c_mz], errors='coerce').reset_index(drop=True)
        inten_series = pd.to_numeric(raw_df.iloc[r_int + 1:, c_int], errors='coerce').reset_index(drop=True)
        tmp = pd.DataFrame({'mass': mass_series, 'intensity': inten_series}).dropna()
        # filter out obviously wrong values
        tmp = tmp[(tmp['mass'] > 0)]
        # sort by mass
        tmp = tmp.sort_values('mass').reset_index(drop=True)
        if len(tmp) == 0:
            return None
        print(f"Using label-block columns -> mass at col {c_mz}, intensity at col {c_int}")
        return tmp['mass'], tmp['intensity']

    res = try_parse_label_blocks(raw)
    if res is not None:
        return res

    # Fallback to header-based parsing
    try:
        df = pd.read_excel(xlsx_path, sheet_name=sheet_name)
    except PermissionError:
        tmpdir = tempfile.mkdtemp(prefix="ms_tmp_")
        tmp_path = os.path.join(tmpdir, os.path.basename(xlsx_path))
        shutil.copy2(xlsx_path, tmp_path)
        try:
            df = pd.read_excel(tmp_path, sheet_name=sheet_name)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)
    return parse_ms_dataframe(df)


def extract_sample_name_from_df(df: pd.DataFrame) -> str | None:
    # df is expected to be read with header=None (raw cells)
    max_rows = min(50, len(df))
    max_cols = min(30, len(df.columns))

    def clean_cell(v):
        if pd.isna(v):
            return None
        s = str(v).strip()
        return s if s and s.lower() != 'nan' else None

    def first_non_empty_right(r, c, limit=15):
        for cc in range(c + 1, min(c + 1 + limit, max_cols)):
            val = clean_cell(df.iloc[r, cc])
            if val:
                return val
        return None

    def first_non_empty_down(r, c, limit=10):
        for rr in range(r + 1, min(r + 1 + limit, max_rows)):
            val = clean_cell(df.iloc[rr, c])
            if val:
                return val
        return None

    # Priority 1: look for a cell containing 'sample name'
    for r in range(max_rows):
        for c in range(max_cols):
            cell = clean_cell(df.iloc[r, c])
            if not cell:
                continue
            low = cell.lower()
            if 'sample name' in low or low == 'sample' or low.startswith('sample:'):
                # Prefer right neighbor(s), then below
                val = first_non_empty_right(r, c) or first_non_empty_down(r, c)
                if val:
                    return val

    # Priority 2: more generic 'sample' key-value pairs
    for r in range(max_rows):
        for c in range(max_cols):
            cell = clean_cell(df.iloc[r, c])
            if not cell:
                continue
            low = cell.lower()
            if 'sample' in low:
                val = first_non_empty_right(r, c) or first_non_empty_down(r, c)
                if val:
                    return val

    return None


def get_sheet_display_names(xlsx_path: str):
    sheet_names = None
    tmpdir_main = None
    try:
        xls = pd.ExcelFile(xlsx_path)
    except PermissionError:
        tmpdir_main = tempfile.mkdtemp(prefix="ms_tmp_")
        tmp_path = os.path.join(tmpdir_main, os.path.basename(xlsx_path))
        shutil.copy2(xlsx_path, tmp_path)
        xls = pd.ExcelFile(tmp_path)
    try:
        sheet_names = list(xls.sheet_names)
    finally:
        try:
            xls.close()
        except Exception:
            pass
        if tmpdir_main:
            shutil.rmtree(tmpdir_main, ignore_errors=True)
    items = []
    sample_names = []
    for s in sheet_names:
        try:
            # Read raw cells to scan labels
            try:
                df = pd.read_excel(xlsx_path, sheet_name=s, header=None, nrows=50)
            except PermissionError:
                tmpdir2 = tempfile.mkdtemp(prefix="ms_tmp_")
                tmp_path2 = os.path.join(tmpdir2, os.path.basename(xlsx_path))
                shutil.copy2(xlsx_path, tmp_path2)
                try:
                    df = pd.read_excel(tmp_path2, sheet_name=s, header=None, nrows=50)
                finally:
                    shutil.rmtree(tmpdir2, ignore_errors=True)
            sample_name = extract_sample_name_from_df(df)
        except Exception:
            sample_name = None
        sample_names.append(sample_name if sample_name else s)
        items.append((sample_name, s))

    # Disambiguate duplicates by appending sheet name
    display_items = []
    counts = {}
    for name, s in items:
        key = name if name else s
        counts[key] = counts.get(key, 0) + 1
    for name, s in items:
        base = name if name else s
        if counts[base] > 1:
            display = f"{base} (sheet: {s})"
        else:
            display = base
        display_items.append((display, s))
    return display_items

def pick_excel_file() -> str:
    # Try a GUI file dialog; fallback to manual input
    path: Optional[str] = None
    try:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        path = filedialog.askopenfilename(
            title="Select Excel file",
            filetypes=[("Excel files", "*.xlsx;*.xls")],
        )
        try:
            root.update()
        except Exception:
            pass
        try:
            root.destroy()
        except Exception:
            pass
    except Exception:
        path = None

    if not path:
        # GUI failed; last-resort console prompt
        try:
            path = input("Enter path to Excel file (.xlsx): ").strip().strip('"')
        except Exception:
            path = ""
    return path

def analyze_and_plot(
    excel_path: str,
    sheet_name: str,
    chosen_display: str,
    selected_amino: str,
    selected_hydroxy: str,
    ionization_mode: str,
    min_intensity: float,
    tolerance: float = 0.5,
    max_mass: float = 1500,
    block_show: bool = True,
    custom_title: Optional[str] = None,
):
    mass_adjustment = get_mass_adjustment(ionization_mode)
    # --- Load Data ---
    mass, intensity = load_ms_from_excel(excel_path, sheet_name)

    # --- Mass Calculations (amino/hydroxy homopolymers and co-polymers) ---
    try:
        amino_weight = molecular_weights[selected_amino]
    except KeyError:
        raise KeyError(f"Unknown amino acid selected: {selected_amino}. Available: {list(molecular_weights.keys())}")
    try:
        hydroxy_weight = molecular_weights[selected_hydroxy]
    except KeyError:
        raise KeyError(f"Unknown hydroxy acid selected: {selected_hydroxy}. Available: {list(molecular_weights.keys())}")
    condensation_products = []  # tuples: (mass, label, category)
    amino_code = amino_three_letter.get(selected_amino, selected_amino)
    print(f"[MS diag] using amino={selected_amino} (code={amino_code}, mw={amino_weight}) hydroxy={selected_hydroxy} (mw={hydroxy_weight}) mode={ionization_mode}")

    # co-polymers of amino + hydroxy
    for n in range(1, 21):
        for m in range(1, 21):
            mass_val = n * amino_weight + m * hydroxy_weight - (n + m - 1) * 18.015
            if mass_val <= max_mass:
                condensation_products.append((mass_val + mass_adjustment, f"{n}x{amino_code} + {m}x{selected_hydroxy}", "copoly"))

    # hydroxy-only oligomers
    for m in range(1, 21):
        mass_val = m * hydroxy_weight - (m - 1) * 18.015
        if mass_val <= max_mass:
            condensation_products.append((mass_val + mass_adjustment, f"{m}x{selected_hydroxy}", "hydroxy"))

    # amino-only oligomers
    for n in range(1, 21):
        mass_val = n * amino_weight - (n - 1) * 18.015
        if mass_val <= max_mass:
            condensation_products.append((mass_val + mass_adjustment, f"{n}x{amino_code}", "amino"))

    # --- Label Matching with de-duplication ---
    matches = []  # each: (idx, m, i, label, delta, category)
    mass_list = list(mass)
    inten_list = list(intensity)
    for idx, (m, i) in enumerate(zip(mass_list, inten_list)):
        is_isotope = False
        for other_mass, other_intensity in zip(mass_list, inten_list):
            if other_mass > m and abs((other_mass - m) - 1.003) <= 0.2 and other_intensity > i:
                is_isotope = True
                break
        if is_isotope:
            continue
        best = None  # (label, delta, category)
        for product_mass, label, cat in condensation_products:
            delta = abs(m - product_mass)
            if delta <= tolerance:
                if best is None or delta < best[1]:
                    best = (label, delta, cat)
        if best is not None:
            matches.append((idx, m, i, best[0], best[1], best[2]))

    # Group by mass bin to avoid duplicate labels at the same mass location
    best_by_bin = {}
    def mass_bin(x):
        return round(float(x), 3)  # 0.001 Da binning
    for idx, m, i, label, delta, cat in matches:
        key = mass_bin(m)
        cur = best_by_bin.get(key)
        if cur is None:
            best_by_bin[key] = (idx, m, i, label, delta, cat)
        else:
            _, _, i0, _, d0, _ = cur
            if (i > i0) or (i == i0 and delta < d0):
                best_by_bin[key] = (idx, m, i, label, delta, cat)

    labels_by_idx = {idx: (label, cat) for idx, m, i, label, delta, cat in best_by_bin.values()}

    # --- Identify Non-Covalent Dimers (2M+H or 2M-H) ---
    # For each identified monomer, look for masses that match 2M±H
    dimer_matches = []
    
    for idx, (label, cat) in labels_by_idx.items():
        # Get the monomer mass (already includes mass adjustment)
        monomer_mass = mass[idx]
        
        # Calculate expected non-covalent dimer mass
        # For non-covalent dimers: 2M + mass_adjustment (since we already have M + adjustment)
        # Essentially: (M + adjustment) * 2 - adjustment = 2M + adjustment
        expected_dimer = 2 * monomer_mass - mass_adjustment
        
        # Search for peaks matching the dimer mass
        for dimer_idx, (m, i) in enumerate(zip(mass, intensity)):
            if dimer_idx not in labels_by_idx:  # Only label unidentified peaks
                delta = abs(m - expected_dimer)
                if delta <= tolerance:
                    dimer_label = f"{label} (2M{'−H' if ionization_mode == 'negative' else '+H'})"
                    dimer_matches.append((dimer_idx, dimer_label, cat, delta))
    
    # Add best dimer matches to labels (avoid duplicates)
    for dimer_idx, dimer_label, cat, delta in dimer_matches:
        if dimer_idx not in labels_by_idx:
            labels_by_idx[dimer_idx] = (dimer_label, cat)

    # --- Filter Data (preserve indices) ---
    # Show ALL masses above the threshold, not just labeled ones
    filtered = [(idx, m, i) for idx, (m, i) in enumerate(zip(mass, intensity)) if i >= min_intensity]
    if filtered:
        filtered_idx, filtered_mass, filtered_intensity = zip(*filtered)
    else:
        filtered_idx, filtered_mass, filtered_intensity = [], [], []

    # --- Plotting ---
    fig, ax = plt.subplots(figsize=(11, 6))
    
    # Y limit - use ALL filtered intensities, not just labeled ones
    if filtered_intensity:
        y_max = max(150, max(filtered_intensity) + 10)
    else:
        y_max = 150
    ax.set_ylim(0, y_max)

    # X limit - use ALL filtered masses, not just labeled ones
    if filtered_mass:
        x_min = max(0, min(filtered_mass) - 10)
        x_max = max(filtered_mass) + 10
        ax.set_xlim(x_min, x_max)
    else:
        ax.set_xlim(0, 1200)

    draggables = []
    stick_lines = {}  # Store references to stick lines for color updates
    
    # Plot ALL filtered masses (above threshold), not just labeled ones
    for idx, m, i in zip(filtered_idx, filtered_mass, filtered_intensity):
        # Plot all stick lines - labeled ones will be colored according to their category
        if idx in labels_by_idx:
            label, cat = labels_by_idx[idx]
            color = 'red' if cat in ('amino', 'copoly') else 'green'
        else:
            # Unlabeled peaks shown in blue
            color = 'blue'
        
        # Create and store the stick line
        stick = ax.vlines(x=m, ymin=0, ymax=i, color=color, linewidth=1)
        stick_lines[idx] = stick
        
        # Add labels only for identified peaks
        if idx in labels_by_idx:
            label, cat = labels_by_idx[idx]
            color = 'red' if cat in ('amino', 'copoly') else 'green'
            txt = ax.text(m, i + 1, label, rotation=90, color=color, fontsize=9, fontweight='bold', ha='center', va='bottom')
            line = Line2D([m, m], [i, i + 1], color=color, linestyle='--', linewidth=0.8)
            ax.add_line(line)
            draggables.append(DraggableTextWithLine(txt, line))

    ax.set_xlabel("m/z (Mass to Charge Ratio)", fontsize=14, fontweight='bold')
    ax.set_ylabel("Intensity (%)", fontsize=14, fontweight='bold')
    title_text = (custom_title or "").strip() or f"MS spectrum of {amino_code} & {selected_hydroxy} ({ionization_mode.capitalize()})"
    ax.set_title(title_text, fontsize=16, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.4)

    # Legend for label colors
    legend_handles = [
        Line2D([0], [0], color='red', lw=2, label='Amino acid / depsipeptide'),
        Line2D([0], [0], color='green', lw=2, label='Hydroxy acid / polyester'),
        Line2D([0], [0], color='blue', lw=2, label='Unidentified peaks'),
    ]
    ax.legend(handles=legend_handles, loc='upper right', frameon=True, framealpha=0.85, fontsize=9)

    plt.tight_layout()
    
    # Create interactive editor for titles, labels, and annotations
    editor = InteractivePlotEditor(
        fig, ax, draggables,
        (filtered_idx, filtered_mass, filtered_intensity),
        labels_by_idx,
        ionization_mode,
        mass_adjustment,
        tolerance,
        stick_lines
    )
    
    # keep references so event callbacks remain alive
    try:
        fig._draggables = draggables  # type: ignore[attr-defined]
        fig._editor = editor  # type: ignore[attr-defined]
    except Exception:
        pass
    
    # Ensure the figure manager is created and window is interactive
    manager = plt.get_current_fig_manager()
    try:
        manager.window.state('normal')  # bring to front
        manager.window.attributes('-topmost', True)  # Make window top-most temporarily
        manager.window.attributes('-topmost', False)  # Then allow normal z-order
    except Exception as e:
        print(f"[Warning] Could not set window state: {e}")
    
    if block_show:
        plt.show(block=True)
    else:
        plt.show(block=False)
        try:
            plt.pause(0.1)
        except Exception:
            pass
    return fig

def main():
    # --- Choose Excel file (browse) ---
    excel_path = None
    while not excel_path or not os.path.exists(excel_path):
        excel_path = pick_excel_file()
        if not os.path.exists(excel_path):
            print(f"File not found: {excel_path}")

    # --- Choose mode: batch or single ---
    sheet_items = get_sheet_display_names(excel_path)
    batch = prompt_yes_no("Batch analyze all sheets?", default=False)

    # Ionization & threshold (global)
    ionization_mode = infer_ion_mode_from_filename(excel_path)
    override = prompt_text(f"Detected ionization mode: '{ionization_mode}'. Press Enter to accept, or type 'positive'/'negative' to override:", default="", title="Ionization mode")
    override = (override or "").strip().lower()
    if override in ("positive", "negative"):
        ionization_mode = override
    min_intensity = prompt_float("Enter the minimum intensity threshold (e.g., 5.0)", default=5.0, title="Intensity threshold") or 5.0

    amino_acids = {k: v for k, v in molecular_weights.items() if k not in ["Glycolic Acid", "Lactic Acid", "Phenyllactic Acid", "hydroxy propionic acid", "hydroxybutyric acid"]}
    hydroxy_acids = {k: v for k, v in molecular_weights.items() if k in ["Glycolic Acid", "Lactic Acid", "Phenyllactic Acid", "hydroxy propionic acid", "hydroxybutyric acid"]}

    if batch:
        # GUI-only per-sheet selections (with Skip). Sequential: wait for user to close each plot before next selection.
        amino_list = list(amino_acids.keys())
        hydroxy_list = list(hydroxy_acids.keys())
        for disp, sname in sheet_items:
            while True:
                sel = gui_select_amino_hydroxy(amino_list, hydroxy_list, disp)
                if not sel:
                    # Skip this sheet
                    break
                amino_choice, hydroxy_choice = sel
                if (amino_choice not in amino_acids) or (hydroxy_choice not in hydroxy_acids):
                    # Safety: invalid selection; re-prompt
                    continue
                print(f"[MS diag] sheet='{sname}' display='{disp}' chosen amino='{amino_choice}' hydroxy='{hydroxy_choice}'")
                # Optional custom title per sample
                custom_title = prompt_text("Optional: enter a custom plot title (leave blank for default)", default=disp, title="Plot title")
                analyze_and_plot(
                    excel_path,
                    sname,
                    disp,
                    amino_choice,
                    hydroxy_choice,
                    ionization_mode,
                    min_intensity,
                    block_show=True,
                    custom_title=custom_title,
                )
                break
    else:
        # Single-sample interactive via GUI
        display_list = [d for d, _ in sheet_items]
        chosen_display = gui_select_from_list(display_list, title="Select sample", message="Select a sample") or display_list[0]
        sheet_idx = display_list.index(chosen_display)
        sheet_name = sheet_items[sheet_idx][1]

        selected_amino = gui_select_from_list(list(amino_acids.keys()), title="Select amino acid", message="Select an amino acid") or list(amino_acids.keys())[0]
        selected_hydroxy = gui_select_from_list(list(hydroxy_acids.keys()), title="Select hydroxy acid", message="Select a hydroxy acid") or list(hydroxy_acids.keys())[0]

        custom_title = prompt_text("Optional: enter a custom plot title (leave blank for default)", default=chosen_display, title="Plot title")
        analyze_and_plot(
            excel_path,
            sheet_name,
            chosen_display,
            selected_amino,
            selected_hydroxy,
            ionization_mode,
            min_intensity,
            block_show=True,
            custom_title=custom_title,
        )


if __name__ == "__main__":
    main()
