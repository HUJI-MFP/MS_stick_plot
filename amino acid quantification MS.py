"""
MS Amino Acid Quantification Tool

Analyzes MS files to calculate the percentage of amino acids in identified oligomers,
normalized by oligomer abundance.

For each sample:
- Identifies oligomers (amino-only, hydroxy-only, and co-polymers)
- Calculates amino acid contribution: (n_amino × abundance%) for each oligomer
- Outputs an Excel file with:
  * Oligomer formula
  * m/z matched
  * Abundance (%)
  * Amino acid count
  * Normalized amino acid contribution
  * Total amino acid percentage across all oligomers
"""
import os
import pandas as pd
from itertools import product
import shutil
import tempfile
from typing import Optional, Tuple, List, Dict
from datetime import datetime

# --- GUI prompt helpers (Tk first, console fallback) ---
def prompt_yes_no(question: str, default: bool = False) -> bool:
    try:
        import tkinter as tk
        from tkinter import messagebox
        root = tk.Tk()
        root.withdraw()
        root.attributes('-topmost', True)
        ans = messagebox.askyesno("Question", question)
        root.destroy()
        return bool(ans)
    except Exception:
        pass
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
    sel = gui_select_from_list(options, title="Select", message=prompt)
    if sel is not None:
        return sel
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

def gui_select_amino_hydroxy(amino_options, hydroxy_options, sample_label: str):
    """Return (amino, hydroxy) or None if user clicks Skip/Cancel."""
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

# --- Molecular weights dictionary ---
molecular_weights = {
    "Glycine": 75.07, "Alanine": 89.09, "Valine": 117.15, "Leucine": 131.17,
    "Isoleucine": 131.17, "Proline": 115.13, "Phenylalanine": 165.19,
    "Tyrosine": 181.19, "Tryptophan": 204.23, "Serine": 105.09,
    "Threonine": 119.12, "Cysteine": 121.16, "Methionine": 149.21,
    "Asparagine": 132.12, "Glutamine": 146.15, "Aspartic Acid": 133.10,
    "Glutamic Acid": 147.13, "Lysine": 146.19, "Arginine": 174.20,
    "Histidine": 155.16, "Glycolic Acid": 76.05, "Lactic Acid": 90.08, "Phenyllactic Acid": 166.17,
    "beta homo-arginine": 188.13, "beta homo-lysine": 160.21, "hydroxy propionic acid": 90.08
}

# 3-letter codes for amino acids
amino_three_letter = {
    "Glycine": "Gly", "Alanine": "Ala", "Valine": "Val", "Leucine": "Leu",
    "Isoleucine": "Ile", "Proline": "Pro", "Phenylalanine": "Phe",
    "Tyrosine": "Tyr", "Tryptophan": "Trp", "Serine": "Ser",
    "Threonine": "Thr", "Cysteine": "Cys", "Methionine": "Met",
    "Asparagine": "Asn", "Glutamine": "Gln", "Aspartic Acid": "Asp",
    "Glutamic Acid": "Glu", "Lysine": "Lys", "Arginine": "Arg",
    "Histidine": "His", "beta homo-arginine": "bhArg", "beta homo-lysine": "bhLys",
}

def infer_ion_mode_from_filename(path: str) -> str:
    fname = os.path.basename(path).lower()
    if 'negative' in fname:
        return 'negative'
    elif 'positive' in fname:
        return 'positive'
    else:
        reply = input("Ionization mode (positive/negative): ").strip().lower()
        return reply if reply in ('positive', 'negative') else 'negative'

def get_mass_adjustment(mode: str) -> float:
    H = 1.00784
    return H if mode == 'positive' else -H

def parse_ms_dataframe(df: pd.DataFrame):
    """Extract mass and intensity from DataFrame with intelligent column detection."""
    # Try header-based detection with strong preference for % intensity
    mass_col, intensity_col = None, None
    abundance_col = None
    percent_col = None
    
    # First pass: find columns
    for col in df.columns:
        col_lower = str(col).lower()
        if 'm/z' in col_lower or 'mass' in col_lower:
            if mass_col is None:
                mass_col = col
        # Look for abundance column (without %)
        if 'abundance' in col_lower and '%' not in col_lower:
            abundance_col = col
        # Look for % abundance column
        if '%' in col_lower and 'abundance' in col_lower:
            percent_col = col
    
    # Strategy: if we found an abundance column, check one column to its right for % abundance
    if abundance_col is not None and percent_col is None:
        col_list = list(df.columns)
        try:
            abundance_idx = col_list.index(abundance_col)
            if abundance_idx + 1 < len(col_list):
                next_col = col_list[abundance_idx + 1]
                next_col_lower = str(next_col).lower()
                # Check if next column has % or looks like percentage data
                if '%' in next_col_lower or 'percent' in next_col_lower:
                    percent_col = next_col
                else:
                    # Check if values are in 0-100 range
                    sample = df[next_col].dropna()
                    if len(sample) > 0 and pd.api.types.is_numeric_dtype(df[next_col]):
                        if ((sample >= 0) & (sample <= 100)).mean() > 0.7:
                            percent_col = next_col
        except (ValueError, KeyError):
            pass
    
    # Use percent column if found, otherwise fall back
    if percent_col is not None:
        intensity_col = percent_col
    elif abundance_col is not None:
        intensity_col = abundance_col

    # Numeric heuristics as last resort
    numeric_cols = df.select_dtypes(include=['number']).columns.tolist()
    if len(numeric_cols) >= 2:
        if mass_col is None:
            mass_col = numeric_cols[0]
        if intensity_col is None:
            # Look for column with values mostly in 0-100 range
            for c in numeric_cols[1:]:
                if c != mass_col:
                    sample = df[c].dropna()
                    if len(sample) > 0 and ((sample >= 0) & (sample <= 100)).mean() > 0.5:
                        intensity_col = c
                        break
            if intensity_col is None:
                intensity_col = numeric_cols[1] if numeric_cols[1] != mass_col else numeric_cols[0]

    if mass_col is None or intensity_col is None:
        raise ValueError("Could not identify mass and intensity columns")

    # Auto-swap if reversed
    mass_sample = df[mass_col].dropna().iloc[:10] if len(df[mass_col].dropna()) >= 10 else df[mass_col].dropna()
    inten_sample = df[intensity_col].dropna().iloc[:10] if len(df[intensity_col].dropna()) >= 10 else df[intensity_col].dropna()
    if len(mass_sample) > 0 and len(inten_sample) > 0:
        if mass_sample.max() <= 100 and inten_sample.max() > 100:
            mass_col, intensity_col = intensity_col, mass_col

    mass = df[mass_col].dropna()
    intensity = df[intensity_col].dropna()
    common_idx = mass.index.intersection(intensity.index)
    mass = mass.loc[common_idx]
    intensity = intensity.loc[common_idx]
    
    # Sort by mass
    sorted_idx = mass.sort_values().index
    mass = mass.loc[sorted_idx]
    intensity = intensity.loc[sorted_idx]
    
    print(f"  Using columns: mass='{mass_col}', intensity='{intensity_col}'")
    return mass.values, intensity.values

def load_ms_from_excel(xlsx_path: str, sheet_name: str):
    """Load mass and intensity from an Excel sheet with robust parsing."""
    temp_dir = None
    work_path = xlsx_path
    try:
        xls = pd.ExcelFile(xlsx_path)
    except PermissionError:
        temp_dir = tempfile.mkdtemp()
        work_path = os.path.join(temp_dir, os.path.basename(xlsx_path))
        shutil.copy2(xlsx_path, work_path)
        xls = pd.ExcelFile(work_path)
    
    try:
        # Label-block parser
        df_raw = pd.read_excel(xls, sheet_name=sheet_name, header=None)
        mz_cell = None
        for r_idx in range(min(50, len(df_raw))):
            for c_idx in range(min(10, len(df_raw.columns))):
                val = df_raw.iloc[r_idx, c_idx]
                if isinstance(val, str) and val.strip().lower() == "m/z":
                    mz_cell = (r_idx, c_idx)
                    break
            if mz_cell:
                break
        
        if mz_cell:
            r, c = mz_cell
            abundance_col = None
            percent_abundance_col = None
            
            # Look for abundance columns - prioritize % abundance
            for offset in range(1, min(5, len(df_raw.columns) - c)):
                hdr = df_raw.iloc[r, c + offset]
                if isinstance(hdr, str):
                    hdr_lower = hdr.lower()
                    # Prioritize columns with % AND abundance
                    if '%' in hdr_lower and 'abundance' in hdr_lower:
                        percent_abundance_col = c + offset
                        break
                    elif 'abundance' in hdr_lower and abundance_col is None:
                        abundance_col = c + offset
            
            # If we found a regular abundance column, check the next column for % abundance
            if abundance_col is not None and percent_abundance_col is None:
                if abundance_col + 1 < len(df_raw.columns):
                    next_hdr = df_raw.iloc[r, abundance_col + 1]
                    if isinstance(next_hdr, str) and ('%' in next_hdr.lower() or 'percent' in next_hdr.lower()):
                        percent_abundance_col = abundance_col + 1
            
            # Use % abundance if found, otherwise fall back to abundance
            final_col = percent_abundance_col if percent_abundance_col is not None else abundance_col
            
            if final_col is not None:
                mass_data = df_raw.iloc[r+1:, c].dropna()
                intensity_data = df_raw.iloc[r+1:, final_col].dropna()
                mass_data = pd.to_numeric(mass_data, errors='coerce').dropna()
                intensity_data = pd.to_numeric(intensity_data, errors='coerce').dropna()
                common_idx = mass_data.index.intersection(intensity_data.index)
                if len(common_idx) > 10:
                    col_name = df_raw.iloc[r, final_col] if isinstance(df_raw.iloc[r, final_col], str) else f"column_{final_col}"
                    print(f"  Parsed via label-block method (using '{col_name}')")
                    mass = mass_data.loc[common_idx].sort_values()
                    intensity = intensity_data.loc[common_idx].loc[mass.index]
                    return mass.values, intensity.values
        
        # Fallback: header-based
        df = pd.read_excel(xls, sheet_name=sheet_name)
        return parse_ms_dataframe(df)
    finally:
        xls.close()
        if temp_dir and os.path.exists(temp_dir):
            try:
                shutil.rmtree(temp_dir)
            except Exception:
                pass

def get_sheet_display_names(xlsx_path: str):
    """Return list of (display_name, sheet_name) tuples. Optimized to load workbook once."""
    temp_dir = None
    work_path = xlsx_path
    
    # Get sheet names from pandas (fast)
    try:
        xls = pd.ExcelFile(xlsx_path)
        sheet_names = xls.sheet_names
        xls.close()
    except PermissionError:
        temp_dir = tempfile.mkdtemp()
        work_path = os.path.join(temp_dir, os.path.basename(xlsx_path))
        shutil.copy2(xlsx_path, work_path)
        xls = pd.ExcelFile(work_path)
        sheet_names = xls.sheet_names
        xls.close()
    
    # Extract all sample names in one pass with openpyxl (optimized)
    try:
        from openpyxl import load_workbook
        wb = load_workbook(work_path, data_only=True, read_only=True)  # read_only mode is faster
        
        sample_names = []
        for sname in sheet_names:
            try:
                ws = wb[sname]
                extracted = None
                
                # OPTIMIZED: Check column C (3), row 1 first (most common location)
                try:
                    val = ws.cell(1, 3).value  # Column C, Row 1
                    if val and str(val).strip() and str(val).strip().lower() not in ['sample name', 'signal name', '']:
                        extracted = str(val).strip()
                except Exception:
                    pass
                
                # If not found at C1, fallback to searching
                if not extracted:
                    for row_idx in range(1, min(21, ws.max_row + 1)):
                        for col_idx in range(1, min(8, ws.max_column + 1)):
                            cell = ws.cell(row_idx, col_idx)
                            if cell.value and isinstance(cell.value, str) and cell.value.lower().strip() == 'sample name':
                                # Found "Sample Name" - check +2 columns first (optimized)
                                merge_end_col = col_idx
                                for merged_range in ws.merged_cells.ranges:
                                    if (merged_range.min_row <= row_idx <= merged_range.max_row and
                                        merged_range.min_col <= col_idx <= merged_range.max_col):
                                        merge_end_col = merged_range.max_col
                                        break
                                
                                # Check +2 first (most common case)
                                check_col = merge_end_col + 2
                                if check_col <= ws.max_column:
                                    val = ws.cell(row_idx, check_col).value
                                    if val and str(val).strip() and str(val).strip().lower() not in ['sample name', 'signal name']:
                                        extracted = str(val).strip()
                                        break
                                
                                # Fallback to other offsets
                                if not extracted:
                                    for offset in [1, 3, 4, 5]:
                                        check_col = merge_end_col + offset
                                        if check_col <= ws.max_column:
                                            val = ws.cell(row_idx, check_col).value
                                            if val and str(val).strip() and str(val).strip().lower() not in ['sample name', 'signal name']:
                                                extracted = str(val).strip()
                                                break
                                break
                        if extracted:
                            break
                
                sample_names.append(extracted if extracted else sname)
            except Exception:
                sample_names.append(sname)
        
        wb.close()
        
    except Exception:
        # Fallback: use sheet names
        sample_names = sheet_names.copy()
    
    finally:
        if temp_dir and os.path.exists(temp_dir):
            try:
                shutil.rmtree(temp_dir)
            except Exception:
                pass
    
    # Disambiguate duplicate names
    name_counts = {}
    for name in sample_names:
        name_counts[name] = name_counts.get(name, 0) + 1
    
    seen = {}
    display_items = []
    for sname, extracted in zip(sheet_names, sample_names):
        if name_counts[extracted] > 1:
            display = f"{extracted} ({sname})"
        else:
            display = extracted
        if display in seen:
            display = f"{display} [{seen[display]}]"
            seen[display] = seen.get(display, 0) + 1
        else:
            seen[display] = 1
        display_items.append((display, sname))
    
    return display_items

def pick_excel_file() -> str:
    """GUI file picker with fallback to manual input."""
    try:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        root.attributes('-topmost', True)
        path = filedialog.askopenfilename(
            title="Select Excel file",
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")]
        )
        root.destroy()
        if path:
            return path
    except Exception:
        pass
    return input("Enter path to Excel file: ").strip().strip('"')

def analyze_sample(
    excel_path: str,
    sheet_name: str,
    sample_display: str,
    selected_amino: str,
    selected_hydroxy: str,
    ionization_mode: str,
    min_intensity: float,
    tolerance: float = 0.5,
    max_mass: float = 1500,
) -> pd.DataFrame:
    """
    Analyze one sample and return a DataFrame with oligomer details and amino acid contributions.
    
    Returns DataFrame with columns:
    - Oligomer: formula string
    - Type: 'Amino-only', 'Hydroxy-only', or 'Co-polymer'
    - m/z: matched peak mass
    - Abundance (%): intensity
    - Amino Count: number of amino acid units in oligomer
    - Hydroxy Count: number of hydroxy acid units in oligomer
    - Normalized Amino Contribution: (Amino Count × Abundance) / 100
    """
    mass_adjustment = get_mass_adjustment(ionization_mode)
    mass, intensity = load_ms_from_excel(excel_path, sheet_name)
    
    amino_weight = molecular_weights[selected_amino]
    hydroxy_weight = molecular_weights[selected_hydroxy]
    amino_code = amino_three_letter.get(selected_amino, selected_amino)
    
    # Generate all oligomers
    condensation_products = []  # (mass, label, type, n_amino, n_hydroxy)
    
    # Co-polymers
    for n in range(1, 21):
        for m in range(1, 21):
            mass_val = n * amino_weight + m * hydroxy_weight - (n + m - 1) * 18.015
            if mass_val <= max_mass:
                label = f"{n}x{amino_code} + {m}x{selected_hydroxy}"
                condensation_products.append((mass_val + mass_adjustment, label, "Co-polymer", n, m))
    
    # Hydroxy-only
    for m in range(1, 21):
        mass_val = m * hydroxy_weight - (m - 1) * 18.015
        if mass_val <= max_mass:
            label = f"{m}x{selected_hydroxy}"
            condensation_products.append((mass_val + mass_adjustment, label, "Hydroxy-only", 0, m))
    
    # Amino-only
    for n in range(1, 21):
        mass_val = n * amino_weight - (n - 1) * 18.015
        if mass_val <= max_mass:
            label = f"{n}x{amino_code}"
            condensation_products.append((mass_val + mass_adjustment, label, "Amino-only", n, 0))
    
    # Match peaks to oligomers
    matches = []  # (idx, m, i, label, delta, type, n_amino, n_hydroxy)
    mass_list = list(mass)
    inten_list = list(intensity)
    
    for idx, (m, i) in enumerate(zip(mass_list, inten_list)):
        if i < min_intensity:
            continue
        
        # Isotope filtering
        is_isotope = False
        for other_mass, other_intensity in zip(mass_list, inten_list):
            if other_mass > m and abs((other_mass - m) - 1.003) <= 0.2 and other_intensity > i:
                is_isotope = True
                break
        if is_isotope:
            continue
        
        best = None  # (label, delta, type, n_amino, n_hydroxy)
        for product_mass, label, oligo_type, n_amino, n_hydroxy in condensation_products:
            delta = abs(m - product_mass)
            if delta <= tolerance:
                if best is None or delta < best[1]:
                    best = (label, delta, oligo_type, n_amino, n_hydroxy)
        
        if best is not None:
            matches.append((idx, m, i, best[0], best[1], best[2], best[3], best[4]))
    
    # Deduplicate by mass bin
    def mass_bin(x):
        return round(float(x), 3)
    
    best_by_bin = {}
    for idx, m, i, label, delta, oligo_type, n_amino, n_hydroxy in matches:
        key = mass_bin(m)
        cur = best_by_bin.get(key)
        if cur is None:
            best_by_bin[key] = (idx, m, i, label, delta, oligo_type, n_amino, n_hydroxy)
        else:
            _, _, i0, _, d0, _, _, _ = cur
            if (i > i0) or (i == i0 and delta < d0):
                best_by_bin[key] = (idx, m, i, label, delta, oligo_type, n_amino, n_hydroxy)
    
    # Build results dataframe
    results = []
    for idx, m, i, label, delta, oligo_type, n_amino, n_hydroxy in best_by_bin.values():
        normalized_contrib = (n_amino * i) / 100.0
        results.append({
            'Oligomer': label,
            'Type': oligo_type,
            'm/z': round(m, 4),
            'Abundance (%)': round(i, 2),
            'Amino Count': n_amino,
            'Hydroxy Count': n_hydroxy,
            'Normalized Amino Contribution': round(normalized_contrib, 4)
        })
    
    df = pd.DataFrame(results)
    if len(df) > 0:
        df = df.sort_values('Abundance (%)', ascending=False)
    
    return df

def main():
    print("\n" + "="*60)
    print("MS AMINO ACID QUANTIFICATION TOOL")
    print("="*60)
    
    # Choose Excel file
    excel_path = None
    while not excel_path or not os.path.exists(excel_path):
        excel_path = pick_excel_file()
        if not os.path.exists(excel_path):
            print(f"File not found: {excel_path}")
    
    print(f"\nSelected file: {excel_path}")
    
    # Get sheets
    sheet_items = get_sheet_display_names(excel_path)
    
    # Batch or single?
    batch = prompt_yes_no("\nBatch analyze all sheets?", default=False)
    
    # Global parameters - ionization mode (GUI prompt)
    ionization_mode = infer_ion_mode_from_filename(excel_path)
    try:
        import tkinter as tk
        from tkinter import messagebox
        root = tk.Tk()
        root.withdraw()
        root.attributes('-topmost', True)
        override = messagebox.askyesnocancel(
            "Ionization Mode",
            f"Detected ionization mode: '{ionization_mode}'\n\nClick 'Yes' for Positive, 'No' for Negative, or 'Cancel' to use detected mode."
        )
        root.destroy()
        if override is True:
            ionization_mode = 'positive'
        elif override is False:
            ionization_mode = 'negative'
        # If None (Cancel), keep detected mode
    except Exception:
        # Fallback to console if GUI fails
        print(f"\nDetected ionization mode: '{ionization_mode}'")
        override = input("Press Enter to accept, or type 'positive'/'negative' to override: ").strip().lower()
        if override in ("positive", "negative"):
            ionization_mode = override
    
    print(f"Using ionization mode: {ionization_mode}")
    
    min_intensity = prompt_float("\nEnter the minimum intensity threshold (e.g., 5.0)", default=5.0)
    if min_intensity is None:
        min_intensity = 5.0
    
    tolerance = prompt_float("\nEnter mass tolerance in Da (default: 0.5)", default=0.5)
    if tolerance is None:
        tolerance = 0.5
    
    amino_acids = {k: v for k, v in molecular_weights.items() if k not in ["Glycolic Acid", "Lactic Acid", "Phenyllactic Acid", "hydroxy propionic acid"]}
    hydroxy_acids = {k: v for k, v in molecular_weights.items() if k in ["Glycolic Acid", "Lactic Acid", "Phenyllactic Acid", "hydroxy propionic acid"]}
    
    # Output directory
    output_dir = os.path.join(os.path.dirname(excel_path), f"MS_Quantification_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
    os.makedirs(output_dir, exist_ok=True)
    print(f"\nOutput directory: {output_dir}")
    
    if batch:
        print("\nBatch mode: analyzing all sheets...")
        for disp, sname in sheet_items:
            print(f"\n{'='*60}")
            print(f"Sample: {disp}")
            print('='*60)
            
            # GUI selection for amino/hydroxy
            selection = gui_select_amino_hydroxy(list(amino_acids.keys()), list(hydroxy_acids.keys()), disp)
            if selection is None:
                print("  Skipped.")
                continue
            
            amino_choice, hydroxy_choice = selection
            print(f"  Amino: {amino_choice}, Hydroxy: {hydroxy_choice}")
            
            # Analyze
            try:
                df_results = analyze_sample(
                    excel_path,
                    sname,
                    disp,
                    amino_choice,
                    hydroxy_choice,
                    ionization_mode,
                    min_intensity,
                    tolerance=tolerance,
                )
                
                if len(df_results) == 0:
                    print("  No oligomers identified.")
                    continue
                
                # Add summary row
                total_amino_contrib = df_results['Normalized Amino Contribution'].sum()
                print(f"  Total normalized amino acid contribution: {total_amino_contrib:.2f}")
                
                # Save to Excel
                safe_name = "".join(c if c.isalnum() or c in (' ', '_', '-') else '_' for c in disp)
                output_path = os.path.join(output_dir, f"{safe_name}.xlsx")
                
                with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
                    # Write results with total row
                    df_results.to_excel(writer, sheet_name='Oligomers', index=False, startrow=0)
                    
                    # Add total row
                    ws = writer.sheets['Oligomers']
                    total_row = len(df_results) + 2  # +1 for header, +1 for next row
                    ws.cell(total_row, 1, "TOTAL")
                    ws.cell(total_row, 7, round(total_amino_contrib, 4))  # Column G (Normalized Amino Contribution)
                    
                    # Bold the total row
                    from openpyxl.styles import Font
                    for col in range(1, 8):
                        ws.cell(total_row, col).font = Font(bold=True)
                    
                    # Write summary
                    summary_df = pd.DataFrame({
                        'Parameter': ['Sample', 'Amino Acid', 'Hydroxy Acid', 'Ionization Mode', 
                                     'Min Intensity (%)', 'Mass Tolerance (Da)', 
                                     'Total Oligomers Identified', 'Total Normalized Amino Contribution'],
                        'Value': [disp, amino_choice, hydroxy_choice, ionization_mode.capitalize(),
                                 min_intensity, tolerance, len(df_results), round(total_amino_contrib, 4)]
                    })
                    summary_df.to_excel(writer, sheet_name='Summary', index=False)
                
                print(f"  Saved: {output_path}")
                
            except Exception as e:
                print(f"  Error analyzing sample: {e}")
                continue
    
    else:
        # Single sample mode
        print("\nAvailable samples:")
        display_list = [d for d, _ in sheet_items]
        chosen_display = choose_from_list(display_list, "Select a sample by number: ")
        sheet_idx = display_list.index(chosen_display)
        sheet_name = sheet_items[sheet_idx][1]
        
        print("\nAvailable amino acids:")
        selected_amino = choose_from_list(list(amino_acids.keys()), "Select an amino acid by number: ")
        print("\nAvailable hydroxy acids:")
        selected_hydroxy = choose_from_list(list(hydroxy_acids.keys()), "Select a hydroxy acid by number: ")
        
        print(f"\nAnalyzing: {chosen_display}")
        print(f"  Amino: {selected_amino}, Hydroxy: {selected_hydroxy}")
        
        try:
            df_results = analyze_sample(
                excel_path,
                sheet_name,
                chosen_display,
                selected_amino,
                selected_hydroxy,
                ionization_mode,
                min_intensity,
                tolerance=tolerance,
            )
            
            if len(df_results) == 0:
                print("  No oligomers identified.")
                return
            
            # Display results
            print("\n" + "="*60)
            print("RESULTS:")
            print("="*60)
            print(df_results.to_string(index=False))
            
            total_amino_contrib = df_results['Normalized Amino Contribution'].sum()
            print("\n" + "="*60)
            print(f"Total normalized amino acid contribution: {total_amino_contrib:.2f}")
            print("="*60)
            
            # Save to Excel
            safe_name = "".join(c if c.isalnum() or c in (' ', '_', '-') else '_' for c in chosen_display)
            output_path = os.path.join(output_dir, f"{safe_name}.xlsx")
            
            with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
                # Write results with total row
                df_results.to_excel(writer, sheet_name='Oligomers', index=False, startrow=0)
                
                # Add total row
                ws = writer.sheets['Oligomers']
                total_row = len(df_results) + 2  # +1 for header, +1 for next row
                ws.cell(total_row, 1, "TOTAL")
                ws.cell(total_row, 7, round(total_amino_contrib, 4))  # Column G (Normalized Amino Contribution)
                
                # Bold the total row
                from openpyxl.styles import Font
                for col in range(1, 8):
                    ws.cell(total_row, col).font = Font(bold=True)
                
                summary_df = pd.DataFrame({
                    'Parameter': ['Sample', 'Amino Acid', 'Hydroxy Acid', 'Ionization Mode', 
                                 'Min Intensity (%)', 'Mass Tolerance (Da)', 
                                 'Total Oligomers Identified', 'Total Normalized Amino Contribution'],
                    'Value': [chosen_display, selected_amino, selected_hydroxy, ionization_mode.capitalize(),
                             min_intensity, tolerance, len(df_results), round(total_amino_contrib, 4)]
                })
                summary_df.to_excel(writer, sheet_name='Summary', index=False)
            
            print(f"\nSaved: {output_path}")
            
        except Exception as e:
            print(f"  Error analyzing sample: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "="*60)
    print("Analysis complete!")
    print("="*60)

if __name__ == "__main__":
    main()
