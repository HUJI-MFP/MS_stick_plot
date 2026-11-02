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
    try:
        # Prevent default matplotlib shortcuts from conflicting with custom editor keys
        matplotlib.rcParams['keymap.yscale'] = []
        matplotlib.rcParams['keymap.xscale'] = []
        matplotlib.rcParams['keymap.grid'] = []
        matplotlib.rcParams['keymap.grid_minor'] = []
        matplotlib.rcParams['keymap.fullscreen'] = []
        matplotlib.rcParams['keymap.all_axes'] = []
    except Exception:
        pass
except Exception:
    pass
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import colors as mcolors
import shutil
import tempfile
from typing import Optional, List

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

def gui_select_from_list(options, title: str, message: str, parent=None) -> Optional[str]:
    try:
        import tkinter as tk
        from tkinter import ttk
        sel = {"value": None}
        owns_root = parent is None
        root = tk.Tk() if owns_root else tk.Toplevel(parent)
        root.title(title)
        try:
            root.attributes('-topmost', True)
        except Exception:
            pass
        root.resizable(False, False)
        if parent is not None:
            try:
                root.transient(parent)
                root.grab_set()
            except Exception:
                pass

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

        if owns_root:
            root.mainloop()
        else:
            root.wait_window()

        if owns_root:
            try:
                root.destroy()
            except Exception:
                pass
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

def gui_select_multiple_monomers(all_monomers, sample_label: str):
    """
    Allow user to select multiple monomers for analysis.
    Returns: list of selected monomer names, or None if Skip/Cancel
    """
    try:
        import tkinter as tk
        from tkinter import ttk, messagebox
        
        result = {"monomers": [], "skip": False}
        root = tk.Tk()
        root.title(f"Select Monomers for {sample_label}")
        root.attributes('-topmost', True)
        root.resizable(False, False)
        
        frm = ttk.Frame(root, padding=10)
        frm.pack(fill='both', expand=True)
        
        ttk.Label(frm, text=f"Sample: {sample_label}", font=('Arial', 10, 'bold')).pack(anchor='w', pady=(0,8))
        ttk.Label(frm, text="Select monomers in this mixture:", font=('Arial', 9)).pack(anchor='w', pady=(0,4))
        
        # Create scrollable listbox with checkboxes
        list_frame = ttk.Frame(frm)
        list_frame.pack(fill='both', expand=True, pady=(0,10))
        
        scrollbar = ttk.Scrollbar(list_frame)
        scrollbar.pack(side='right', fill='y')
        
        listbox = tk.Listbox(list_frame, selectmode='multiple', height=15, width=40, 
                            yscrollcommand=scrollbar.set, exportselection=False)
        listbox.pack(side='left', fill='both', expand=True)
        scrollbar.config(command=listbox.yview)
        
        sorted_monomers = sorted(all_monomers)
        for monomer in sorted_monomers:
            listbox.insert('end', monomer)
        
        def ok():
            selected_indices = listbox.curselection()
            if len(selected_indices) == 0:
                try:
                    messagebox.showwarning("Selection required", "Please select at least one monomer.")
                except Exception:
                    pass
                return
            result["monomers"] = [sorted_monomers[i] for i in selected_indices]
            root.destroy()
        
        def skip():
            result["skip"] = True
            root.destroy()
        
        btns = ttk.Frame(frm)
        btns.pack(fill='x')
        ttk.Button(btns, text="Skip", command=skip).pack(side='right', padx=(8,0))
        ttk.Button(btns, text="OK", command=ok).pack(side='right')
        
        root.mainloop()
        
        if result["skip"]:
            return None
        return result["monomers"]
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
        if event.inaxes != self.text.axes:
            return
        contains, _ = self.text.contains(event)
        if not contains:
            return
        x0, y0 = self.text.get_position()
        self.press = (x0, y0, event.xdata, event.ydata)

    def on_motion(self, event):
        if self.press is None or event.inaxes != self.text.axes:
            return
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
    
    def __init__(self, fig, ax, draggables, filtered_data, labels_by_idx, ionization_mode, mass_adjustment, tolerance, stick_lines, legend_items=None, legend_config=None, synonym_map=None, category_colors=None):
        self.fig = fig
        self.ax = ax
        self.draggables = draggables
        self.filtered_idx, self.filtered_mass, self.filtered_intensity = filtered_data
        self.labels_by_idx = labels_by_idx
        self.ionization_mode = ionization_mode
        self.mass_adjustment = mass_adjustment
        self.tolerance = tolerance
        self.stick_lines = stick_lines  # Dictionary mapping idx to LineCollection or vlines object
        self.legend_items = [dict(item) for item in (legend_items or [])]
        self.default_legend_items = [dict(item) for item in (legend_items or [])]
        self.legend_config = dict(legend_config or {})
        self.legend_visible = True
        self.synonym_map = dict(synonym_map or {})
        self.default_synonym_map = dict(synonym_map or {})
        self.category_colors = dict(category_colors or {})
        
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
        print("  't'  - Edit axis/tick text style (size/bold)")
        print("  'z'  - Edit X/Y axis limits")
        print("  'b'  - Toggle top/right border (spines)")
        print("  'l'  - Edit legend (rename/delete/add/hide)")
        print("  'g'  - Grid options (show/hide/toggle)")
        print("  'n'  - Remove grid (clean look)")
        print("  'f'  - Edit figure size")
        print("  'p'  - Peak label font size options")
        print("  's'  - Edit monomer synonyms in labels")
        print("  'c'  - Edit amino/hydroxy label colors")
        print("  'a'  - Add label to a peak (click near a peak)")
        print("       Stick color will change to match label color")
        print("  'r'  - Remove label (click near a label)")
        print("       Stick color will change back to light gray (unidentified)")
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
        elif event.key in ('t', 'T'):
            self.edit_axis_text_style()
        elif event.key in ('z', 'Z'):
            self.edit_axis_limits()
        elif event.key in ('b', 'B'):
            self.toggle_top_right_spines()
        elif event.key in ('l', 'L'):
            self.edit_legend()
        elif event.key in ('g', 'G'):
            self.edit_grid()
        elif event.key in ('n', 'N'):
            self.set_grid_visible(False)
        elif event.key in ('f', 'F'):
            self.edit_figure_size()
        elif event.key in ('p', 'P'):
            self.edit_peak_label_font()
        elif event.key in ('s', 'S'):
            self.edit_synonyms()
        elif event.key in ('c', 'C'):
            self.edit_category_colors()
        elif event.key in ('a', 'A'):
            print("ADD MODE: Click near a peak to add a label...")
            self.add_mode = True
        elif event.key in ('r', 'R'):
            print("REMOVE MODE: Click near a label to remove it...")
            self.delete_mode = True

    def _apply_synonyms_to_text(self, text: str) -> str:
        if not text:
            return text
        updated = str(text)
        for full_name in sorted(self.synonym_map.keys(), key=len, reverse=True):
            code = self.synonym_map.get(full_name, full_name)
            updated = updated.replace(full_name, code)
        return updated

    def _refresh_all_label_synonyms(self):
        # Update canonical labels storage
        for idx, (label, cat) in list(self.labels_by_idx.items()):
            self.labels_by_idx[idx] = (self._apply_synonyms_to_text(label), cat)

        # Update visible text objects
        for draggable in self.draggables:
            try:
                x_pos = draggable.line.get_xdata()[0]
                matched_idx = None
                for idx, m in zip(self.filtered_idx, self.filtered_mass):
                    if abs(m - x_pos) < 0.01:
                        matched_idx = idx
                        break
                if matched_idx is not None and matched_idx in self.labels_by_idx:
                    draggable.text.set_text(self.labels_by_idx[matched_idx][0])
                else:
                    draggable.text.set_text(self._apply_synonyms_to_text(draggable.text.get_text()))
            except Exception:
                pass
        self._redraw()

    def _replace_synonym_code_in_labels(self, old_code: str, new_code: str):
        old_code = str(old_code or "").strip()
        new_code = str(new_code or "").strip()
        if not old_code or not new_code or old_code == new_code:
            return

        # Update canonical labels storage
        for idx, (label, cat) in list(self.labels_by_idx.items()):
            self.labels_by_idx[idx] = (str(label).replace(old_code, new_code), cat)

        # Update visible text objects
        for draggable in self.draggables:
            try:
                current_text = draggable.text.get_text()
                draggable.text.set_text(str(current_text).replace(old_code, new_code))
            except Exception:
                pass

        self._redraw()

    def _category_group_for_label(self, cat: str) -> str:
        c = str(cat or '').strip().lower()
        if c in ('amino', 'depsipeptide'):
            return 'amino'
        if c == 'hydroxy':
            return 'hydroxy'
        if c == 'copoly':
            return 'copoly'
        if is_hydroxy_monomer(c):
            return 'hydroxy'
        return 'other'

    def _find_idx_for_x_position(self, x_pos: float) -> Optional[int]:
        for idx, m in zip(self.filtered_idx, self.filtered_mass):
            try:
                if abs(float(m) - float(x_pos)) < 0.01:
                    return idx
            except Exception:
                pass
        return None

    def _sync_legend_category_colors(self):
        changed = False
        for item in self.legend_items:
            group = item.get('color_group')
            if group and group in self.category_colors:
                item['color'] = self.category_colors[group]
                changed = True
        if changed:
            self.apply_legend()
        else:
            self._redraw()

    def _apply_category_colors(self):
        for idx, (_, cat) in self.labels_by_idx.items():
            group = self._category_group_for_label(cat)
            color = self.category_colors.get(group)
            if color:
                self.update_stick_color(idx, color)

        for draggable in self.draggables:
            try:
                x_pos = draggable.line.get_xdata()[0]
                idx = self._find_idx_for_x_position(x_pos)
                if idx is None or idx not in self.labels_by_idx:
                    continue
                _, cat = self.labels_by_idx[idx]
                group = self._category_group_for_label(cat)
                color = self.category_colors.get(group)
                if color:
                    draggable.text.set_color(color)
                    draggable.line.set_color(color)
            except Exception:
                pass

        self._sync_legend_category_colors()

    def edit_category_colors(self):
        try:
            while True:
                action = self._select_from_list(
                    [
                        "Set amino color",
                        "Set hydroxy color",
                        "Set copoly color",
                        "Reset default colors",
                        "Done",
                    ],
                    title="Label Colors",
                    message="Change category colors in real time"
                )
                if not action or action == "Done":
                    break

                if action == "Reset default colors":
                    self.category_colors['amino'] = 'red'
                    self.category_colors['hydroxy'] = 'blue'
                    if 'copoly' not in self.category_colors:
                        self.category_colors['copoly'] = 'black'
                    self._apply_category_colors()
                    continue

                if action == "Set amino color":
                    group = 'amino'
                    default = self.category_colors.get(group, 'red')
                elif action == "Set hydroxy color":
                    group = 'hydroxy'
                    default = self.category_colors.get(group, 'blue')
                else:
                    group = 'copoly'
                    default = self.category_colors.get(group, 'black')

                new_color = self._ask_text("Category Color", f"Enter color for {group} labels (name or #hex)", default=default)
                if not new_color:
                    continue
                new_color = new_color.strip()
                if not mcolors.is_color_like(new_color):
                    print(f"Invalid color '{new_color}'.")
                    continue
                self.category_colors[group] = new_color
                self._apply_category_colors()
        except Exception as e:
            print(f"Error editing category colors: {e}")

    def _redraw(self):
        try:
            self.fig.canvas.draw_idle()
        except Exception:
            self.fig.canvas.draw()

    def _get_tk_parent(self):
        try:
            manager = getattr(self.fig.canvas, 'manager', None)
            return getattr(manager, 'window', None)
        except Exception:
            return None

    def _ask_text(self, title: str, question: str, default: Optional[str] = None) -> Optional[str]:
        try:
            from tkinter import simpledialog
            parent = self._get_tk_parent()
            if parent is not None:
                ans = simpledialog.askstring(title, question, initialvalue=default, parent=parent)
                if ans is None:
                    return default
                return ans
        except Exception:
            pass
        return prompt_text(question, default=default, title=title)

    def _ask_float(self, title: str, question: str, default: Optional[float] = None) -> Optional[float]:
        s_default = None if default is None else str(default)
        try:
            from tkinter import simpledialog
            parent = self._get_tk_parent()
            if parent is not None:
                ans = simpledialog.askstring(title, question, initialvalue=s_default, parent=parent)
                if ans is None or ans.strip() == "":
                    return default
                try:
                    return float(ans)
                except ValueError:
                    return default
        except Exception:
            pass
        return prompt_float(question, default=default, title=title)

    def _select_from_list(self, options, title: str, message: str) -> Optional[str]:
        return gui_select_from_list(options, title=title, message=message, parent=self._get_tk_parent())

    def _apply_tick_style(self, axis: str, size: float, is_bold: bool):
        weight = 'bold' if is_bold else 'normal'
        if axis in ('x', 'both'):
            for lbl in self.ax.get_xticklabels():
                lbl.set_fontsize(size)
                lbl.set_fontweight(weight)
        if axis in ('y', 'both'):
            for lbl in self.ax.get_yticklabels():
                lbl.set_fontsize(size)
                lbl.set_fontweight(weight)

    def _current_tick_size(self, axis: str) -> float:
        labels = self.ax.get_xticklabels() if axis == 'x' else self.ax.get_yticklabels()
        for lbl in labels:
            txt = lbl.get_text()
            if txt is not None and str(txt).strip() != '':
                return float(lbl.get_fontsize())
        return 12.0

    def _current_tick_is_bold(self, axis: str) -> bool:
        labels = self.ax.get_xticklabels() if axis == 'x' else self.ax.get_yticklabels()
        for lbl in labels:
            txt = lbl.get_text()
            if txt is not None and str(txt).strip() != '':
                fw = str(lbl.get_fontweight()).lower()
                return fw in ('bold', 'heavy', 'semibold', 'demibold', 'black', '700', '800', '900')
        return False

    def _build_legend_handles(self):
        handles = []
        for item in self.legend_items:
            if not item.get('enabled', True):
                continue
            color = item.get('color', 'black')
            label = item.get('label', 'Item')
            handles.append(Line2D([0], [0], color=color, lw=2, label=label))
        return handles

    def apply_legend(self):
        if self.ax.get_legend() is not None:
            self.ax.get_legend().remove()
        if not self.legend_visible:
            self._redraw()
            return
        handles = self._build_legend_handles()
        if handles:
            self.ax.legend(
                handles=handles,
                loc=self.legend_config.get('loc', 'upper right'),
                frameon=self.legend_config.get('frameon', True),
                framealpha=self.legend_config.get('framealpha', 0.85),
                fontsize=self.legend_config.get('fontsize', 9),
                title=self.legend_config.get('title', None),
            )
        self._redraw()

    def _get_peak_labels(self):
        labels = []
        for draggable in self.draggables:
            try:
                if draggable.text is not None:
                    labels.append(draggable.text)
            except Exception:
                pass
        return labels

    def _apply_peak_label_fontsize(self, size: float):
        labels = self._get_peak_labels()
        for txt in labels:
            txt.set_fontsize(size)
        self._redraw()

    def edit_peak_label_font(self):
        try:
            labels = self._get_peak_labels()
            if not labels:
                print("No peak labels available to resize.")
                return
            current_size = float(labels[0].get_fontsize())
            action = self._select_from_list(
                [
                    "Increase label font (+1)",
                    "Decrease label font (-1)",
                    "Set exact label font size",
                ],
                title="Peak Label Font",
                message=f"Current label font size: {current_size:.1f}"
            )
            if not action:
                return
            if action == "Increase label font (+1)":
                self._apply_peak_label_fontsize(current_size + 1.0)
            elif action == "Decrease label font (-1)":
                new_size = max(1.0, current_size - 1.0)
                self._apply_peak_label_fontsize(new_size)
            else:
                new_size = self._ask_float("Peak Label Font", "Enter label font size", default=current_size)
                if new_size is None or float(new_size) <= 0:
                    return
                self._apply_peak_label_fontsize(float(new_size))
        except Exception as e:
            print(f"Error editing peak label font: {e}")

    def edit_synonyms(self):
        try:
            if not self.synonym_map:
                print("No synonyms available to edit.")
                return

            while True:
                action = self._select_from_list(
                    [
                        "Change one synonym",
                        "Reset to default synonyms",
                        "Done",
                    ],
                    title="Synonym Editor",
                    message="Choose synonym action"
                )
                if not action or action == "Done":
                    break

                if action == "Reset to default synonyms":
                    # Replace current codes in-place so labels update immediately
                    for monomer, current_code in list(self.synonym_map.items()):
                        default_code = self.default_synonym_map.get(monomer, current_code)
                        self._replace_synonym_code_in_labels(current_code, default_code)
                    self.synonym_map = dict(self.default_synonym_map)
                    continue

                monomer = self._select_from_list(
                    sorted(self.synonym_map.keys()),
                    title="Select Monomer",
                    message="Pick a monomer name to change its synonym"
                )
                if not monomer:
                    continue
                current = self.synonym_map.get(monomer, monomer)
                new_code = self._ask_text(
                    "Edit Synonym",
                    f"Enter new synonym for '{monomer}'",
                    default=current
                )
                if new_code is None:
                    continue
                new_code = new_code.strip()
                if not new_code:
                    print("Synonym cannot be empty.")
                    continue
                old_code = self.synonym_map.get(monomer, monomer)
                self.synonym_map[monomer] = new_code
                self._replace_synonym_code_in_labels(old_code, new_code)
        except Exception as e:
            print(f"Error editing synonyms: {e}")
    
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
            current_title = self.ax.get_title()
            new_title = self._ask_text("Edit Title", "Enter new plot title:", default=current_title)
            if new_title is not None:
                self.ax.set_title(new_title, fontsize=16, fontweight='bold')
                self._redraw()
        except Exception as e:
            print(f"Error editing title: {e}")
    
    def edit_xlabel(self):
        """Edit the X-axis label"""
        try:
            current_label = self.ax.get_xlabel()
            new_label = self._ask_text("Edit X-axis Label", "Enter new X-axis label:", default=current_label)
            if new_label is not None:
                self.ax.set_xlabel(new_label, fontsize=14, fontweight='bold')
                self._redraw()
        except Exception as e:
            print(f"Error editing X-axis label: {e}")
    
    def edit_ylabel(self):
        """Edit the Y-axis label"""
        try:
            current_label = self.ax.get_ylabel()
            new_label = self._ask_text("Edit Y-axis Label", "Enter new Y-axis label:", default=current_label)
            if new_label is not None:
                self.ax.set_ylabel(new_label, fontsize=14, fontweight='bold')
                self._redraw()
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
            from tkinter import simpledialog
            temp_root = None
            parent = self._get_tk_parent()
            if parent is None:
                temp_root = tk.Tk()
                temp_root.withdraw()
                try:
                    temp_root.attributes('-topmost', True)
                except Exception:
                    pass
                parent = temp_root
            
            # Ask for label text
            label_text = simpledialog.askstring(
                "Add Label",
                f"Enter label for peak at m/z={nearest_m:.2f}:",
                parent=parent,
            )
            if not label_text:
                if temp_root is not None:
                    temp_root.destroy()
                return
            
            # Ask for color
            color_choice = simpledialog.askstring(
                "Label Color",
                "Enter color (red/green/blue/black) or leave blank for red:",
                initialvalue="red",
                parent=parent,
            )
            if temp_root is not None:
                temp_root.destroy()
            
            color = color_choice if color_choice else 'red'
            if not mcolors.is_color_like(color):
                print(f"Invalid color '{color}'. Using red.")
                color = 'red'
            
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
            
            self._redraw()
            
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
        if target_draggable:
            try:
                xlim = self.ax.get_xlim()
                ylim = self.ax.get_ylim()
                x_span = max(1e-9, abs(xlim[1] - xlim[0]))
                y_span = max(1e-9, abs(ylim[1] - ylim[0]))
                norm_dist = ((target_draggable.text.get_position()[0] - x) / x_span) ** 2 + ((target_draggable.text.get_position()[1] - y) / y_span) ** 2
                if norm_dist > 0.01:
                    print("No nearby label found to remove.")
                    return

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
                
                # Change stick color back to lightgray (unidentified)
                if removed_idx is not None:
                    self.update_stick_color(removed_idx, 'lightgray')
                
                self._redraw()
            except Exception as e:
                print(f"Error deleting label: {e}")

    def edit_axis_limits(self):
        try:
            x0, x1 = self.ax.get_xlim()
            y0, y1 = self.ax.get_ylim()
            new_x0 = self._ask_float("Axis Limits", "X-axis minimum", default=float(x0))
            new_x1 = self._ask_float("Axis Limits", "X-axis maximum", default=float(x1))
            new_y0 = self._ask_float("Axis Limits", "Y-axis minimum", default=float(y0))
            new_y1 = self._ask_float("Axis Limits", "Y-axis maximum", default=float(y1))
            if new_x0 is None or new_x1 is None or new_y0 is None or new_y1 is None:
                return
            x_min = float(new_x0)
            x_max = float(new_x1)
            y_min = float(new_y0)
            y_max = float(new_y1)
            if x_min >= x_max or y_min >= y_max:
                print("Invalid axis limits: min must be smaller than max.")
                return
            self.ax.set_xlim(x_min, x_max)
            self.ax.set_ylim(y_min, y_max)
            self._redraw()
        except Exception as e:
            print(f"Error editing axis limits: {e}")

    def edit_axis_text_style(self):
        try:
            while True:
                action = self._select_from_list(
                    [
                        "X-axis title style",
                        "Y-axis title style",
                        "X ticks style",
                        "Y ticks style",
                        "Both ticks style",
                        "Done",
                    ],
                    title="Axis and Tick Text Style",
                    message="Choose what to edit"
                )

                if not action or action == "Done":
                    break

                if action == "X-axis title style":
                    current_size = float(self.ax.xaxis.label.get_fontsize())
                    current_weight = str(self.ax.xaxis.label.get_fontweight()).lower()
                    current_is_bold = current_weight in ('bold', 'heavy', 'semibold', 'demibold', 'black', '700', '800', '900')
                    size = self._ask_float("X-axis Title Style", "Font size", default=current_size)
                    if size is None or float(size) <= 0:
                        continue
                    bold_choice = self._select_from_list(
                        ['Normal', 'Bold'],
                        title="X-axis Title Style",
                        message="Font weight"
                    )
                    if not bold_choice:
                        continue
                    self.ax.xaxis.label.set_fontsize(float(size))
                    self.ax.xaxis.label.set_fontweight('bold' if bold_choice == 'Bold' else 'normal')
                    self._redraw()

                elif action == "Y-axis title style":
                    current_size = float(self.ax.yaxis.label.get_fontsize())
                    current_weight = str(self.ax.yaxis.label.get_fontweight()).lower()
                    current_is_bold = current_weight in ('bold', 'heavy', 'semibold', 'demibold', 'black', '700', '800', '900')
                    size = self._ask_float("Y-axis Title Style", "Font size", default=current_size)
                    if size is None or float(size) <= 0:
                        continue
                    bold_choice = self._select_from_list(
                        ['Normal', 'Bold'],
                        title="Y-axis Title Style",
                        message="Font weight"
                    )
                    if not bold_choice:
                        continue
                    self.ax.yaxis.label.set_fontsize(float(size))
                    self.ax.yaxis.label.set_fontweight('bold' if bold_choice == 'Bold' else 'normal')
                    self._redraw()

                elif action == "X ticks style":
                    current_size = self._current_tick_size('x')
                    current_is_bold = self._current_tick_is_bold('x')
                    size = self._ask_float("X Ticks Style", "Tick font size", default=current_size)
                    if size is None or float(size) <= 0:
                        continue
                    bold_default = 'Bold' if current_is_bold else 'Normal'
                    bold_choice = self._select_from_list(
                        ['Normal', 'Bold'],
                        title="X Ticks Style",
                        message=f"Tick font weight (current: {bold_default})"
                    )
                    if not bold_choice:
                        continue
                    self._apply_tick_style('x', float(size), bold_choice == 'Bold')
                    self._redraw()

                elif action == "Y ticks style":
                    current_size = self._current_tick_size('y')
                    current_is_bold = self._current_tick_is_bold('y')
                    size = self._ask_float("Y Ticks Style", "Tick font size", default=current_size)
                    if size is None or float(size) <= 0:
                        continue
                    bold_default = 'Bold' if current_is_bold else 'Normal'
                    bold_choice = self._select_from_list(
                        ['Normal', 'Bold'],
                        title="Y Ticks Style",
                        message=f"Tick font weight (current: {bold_default})"
                    )
                    if not bold_choice:
                        continue
                    self._apply_tick_style('y', float(size), bold_choice == 'Bold')
                    self._redraw()

                elif action == "Both ticks style":
                    current_size = (self._current_tick_size('x') + self._current_tick_size('y')) / 2.0
                    current_is_bold = self._current_tick_is_bold('x') and self._current_tick_is_bold('y')
                    size = self._ask_float("Ticks Style", "Tick font size", default=current_size)
                    if size is None or float(size) <= 0:
                        continue
                    bold_default = 'Bold' if current_is_bold else 'Normal'
                    bold_choice = self._select_from_list(
                        ['Normal', 'Bold'],
                        title="Ticks Style",
                        message=f"Tick font weight (current: {bold_default})"
                    )
                    if not bold_choice:
                        continue
                    self._apply_tick_style('both', float(size), bold_choice == 'Bold')
                    self._redraw()
        except Exception as e:
            print(f"Error editing axis/tick text style: {e}")

    def toggle_top_right_spines(self):
        try:
            top_visible = self.ax.spines['top'].get_visible()
            right_visible = self.ax.spines['right'].get_visible()
            make_visible = not (top_visible and right_visible)
            self.ax.spines['top'].set_visible(make_visible)
            self.ax.spines['right'].set_visible(make_visible)
            self._redraw()
            if make_visible:
                print("Top/right border shown.")
            else:
                print("Top/right border hidden.")
        except Exception as e:
            print(f"Error toggling top/right border: {e}")

    def set_grid_visible(self, visible: bool):
        try:
            if visible:
                self.ax.grid(True, which='major', axis='both', linestyle='--', alpha=0.4)
            else:
                # Fully remove grid lines (major + minor) and keep a clean white background
                self.ax.grid(False, which='both', axis='both')
                try:
                    self.ax.minorticks_off()
                except Exception:
                    pass
                for line in self.ax.get_xgridlines() + self.ax.get_ygridlines():
                    try:
                        line.set_visible(False)
                    except Exception:
                        pass
                self.ax.set_facecolor('white')
                self.fig.patch.set_facecolor('white')
            self._redraw()
            if visible:
                print("Grid shown.")
            else:
                print("Grid removed.")
        except Exception as e:
            print(f"Error setting grid visibility: {e}")

    def toggle_grid(self):
        try:
            grid_visible = any(gl.get_visible() for gl in self.ax.get_xgridlines() + self.ax.get_ygridlines())
            self.set_grid_visible(not grid_visible)
        except Exception as e:
            print(f"Error toggling grid: {e}")

    def edit_grid(self):
        try:
            action = self._select_from_list(
                ["Show grid", "Remove grid", "Toggle grid"],
                title="Grid",
                message="Choose grid action"
            )
            if not action:
                return
            if action == "Show grid":
                self.set_grid_visible(True)
            elif action == "Remove grid":
                self.set_grid_visible(False)
            else:
                self.toggle_grid()
        except Exception as e:
            print(f"Error editing grid: {e}")

    def edit_figure_size(self):
        try:
            w, h = self.fig.get_size_inches()
            new_w = self._ask_float("Figure Size", "Figure width (inches)", default=float(w))
            new_h = self._ask_float("Figure Size", "Figure height (inches)", default=float(h))
            if new_w is None or new_h is None:
                return
            width = float(new_w)
            height = float(new_h)
            if width <= 0 or height <= 0:
                print("Figure size must be positive.")
                return
            self.fig.set_size_inches(width, height, forward=True)
            self.fig.tight_layout()
            self._redraw()
        except Exception as e:
            print(f"Error editing figure size: {e}")

    def edit_legend(self):
        try:
            while True:
                action = gui_select_from_list(
                    [
                        "Rename legend item",
                        "Delete legend item",
                        "Add legend item",
                        "Increase legend font (+1)",
                        "Decrease legend font (-1)",
                        "Set legend font size",
                        "Toggle legend show/hide",
                        "Move legend position",
                        "Set legend title",
                        "Restore default legend",
                        "Done",
                    ],
                    title="Legend Editor",
                    message="Choose legend edit action",
                    parent=self._get_tk_parent()
                )

                if not action or action == "Done":
                    break

                if action == "Rename legend item":
                    visible_items = [it for it in self.legend_items if it.get('enabled', True)]
                    if not visible_items:
                        print("No legend items to rename.")
                        continue
                    labels = [it.get('label', 'Item') for it in visible_items]
                    current = self._select_from_list(labels, title="Rename", message="Select item to rename")
                    if not current:
                        continue
                    target = next((it for it in visible_items if it.get('label', 'Item') == current), None)
                    if target is None:
                        continue
                    new_name = self._ask_text("Rename Legend Item", "New legend text", default=current)
                    if new_name is not None and new_name.strip() != "":
                        target['label'] = new_name.strip()
                        self.apply_legend()

                elif action == "Delete legend item":
                    visible_items = [it for it in self.legend_items if it.get('enabled', True)]
                    if not visible_items:
                        print("No legend items to delete.")
                        continue
                    labels = [it.get('label', 'Item') for it in visible_items]
                    current = self._select_from_list(labels, title="Delete", message="Select item to delete")
                    if not current:
                        continue
                    for item in self.legend_items:
                        if item.get('label', 'Item') == current and item.get('enabled', True):
                            item['enabled'] = False
                            break
                    self.apply_legend()

                elif action == "Add legend item":
                    new_label = self._ask_text("Add Legend Item", "Legend text", default="Custom item")
                    if not new_label or not new_label.strip():
                        continue
                    new_color = self._ask_text("Add Legend Item", "Color (e.g., red, #1f77b4)", default="black")
                    new_color = (new_color or "black").strip()
                    if not mcolors.is_color_like(new_color):
                        print(f"Invalid color '{new_color}'. Using black.")
                        new_color = 'black'
                    self.legend_items.append({'label': new_label.strip(), 'color': new_color, 'enabled': True})
                    self.apply_legend()

                elif action == "Toggle legend show/hide":
                    self.legend_visible = not self.legend_visible
                    self.apply_legend()

                elif action == "Increase legend font (+1)":
                    current = float(self.legend_config.get('fontsize', 9))
                    self.legend_config['fontsize'] = current + 1.0
                    self.apply_legend()

                elif action == "Decrease legend font (-1)":
                    current = float(self.legend_config.get('fontsize', 9))
                    self.legend_config['fontsize'] = max(1.0, current - 1.0)
                    self.apply_legend()

                elif action == "Set legend font size":
                    current = float(self.legend_config.get('fontsize', 9))
                    new_size = self._ask_float("Legend Font", "Legend font size", default=current)
                    if new_size is None or float(new_size) <= 0:
                        continue
                    self.legend_config['fontsize'] = float(new_size)
                    self.apply_legend()

                elif action == "Move legend position":
                    positions = ['best', 'upper right', 'upper left', 'lower right', 'lower left', 'center right', 'center left', 'upper center', 'lower center', 'center']
                    pos = self._select_from_list(positions, title="Legend Position", message="Choose legend position")
                    if pos:
                        self.legend_config['loc'] = pos
                        self.apply_legend()

                elif action == "Set legend title":
                    current_title = self.legend_config.get('title', '') or ''
                    new_title = self._ask_text("Legend Title", "Legend title (leave blank to clear)", default=current_title)
                    self.legend_config['title'] = (new_title or '').strip() or None
                    self.apply_legend()

                elif action == "Restore default legend":
                    self.legend_items = [dict(item) for item in self.default_legend_items]
                    self.legend_visible = True
                    self.apply_legend()
        except Exception as e:
            print(f"Error editing legend: {e}")

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

HYDROXY_ACID_NAMES = {
    "glycolic acid",
    "lactic acid",
    "phenyllactic acid",
    "hydroxy propionic acid",
    "hydroxybutyric acid",
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
    "Glycolic Acid": "Glc", "Lactic Acid": "Lac", "Phenyllactic Acid": "PLA",
    "hydroxy propionic acid": "HPA", "Hydroxy propionic acid": "HPA",
    "hydroxybutyric acid": "HBA", "Hydroxybutyric acid": "HBA",
}


def get_monomer_code(name: str) -> str:
    code = amino_three_letter.get(name)
    if code:
        return code
    # fallback: compact 3-letter uppercase code from words
    tokens = [t for t in str(name).replace('-', ' ').split() if t]
    if len(tokens) >= 2:
        return "".join(t[0].upper() for t in tokens[:3])
    return str(name)[:3].upper()


def build_default_synonym_map() -> dict:
    return {name: get_monomer_code(name) for name in molecular_weights.keys()}


def is_hydroxy_monomer(name: str) -> bool:
    return str(name).strip().lower() in HYDROXY_ACID_NAMES

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
        if c_mz is None:
            return None
        c_mz_idx = int(c_mz)

        # find intensity header: "abundance %" preferably two cells to the right
        r_int = r_mz
        c_int = None
        candidate = clean(raw_df.iloc[r_mz, c_mz_idx + 2]) if c_mz_idx + 2 < max_cols else None
        if candidate and is_intensity_header(candidate):
            c_int = c_mz_idx + 2
        else:
            # search the same row for a header containing 'abundance' and/or '%'
            for c in range(c_mz_idx + 1, max_cols):
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
        c_int_idx = int(c_int)

        # Extract columns below headers
        mass_series = pd.to_numeric(raw_df.iloc[r_mz + 1:, c_mz_idx], errors='coerce').reset_index(drop=True)
        inten_series = pd.to_numeric(raw_df.iloc[r_int + 1:, c_int_idx], errors='coerce').reset_index(drop=True)
        tmp = pd.DataFrame({'mass': mass_series, 'intensity': inten_series}).dropna()
        # filter out obviously wrong values
        tmp = tmp[(tmp['mass'] > 0)]
        # sort by mass
        tmp = tmp.sort_values('mass').reset_index(drop=True)
        if len(tmp) == 0:
            return None
        print(f"Using label-block columns -> mass at col {c_mz_idx}, intensity at col {c_int_idx}")
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
    amino_code = get_monomer_code(selected_amino)
    hydroxy_code = get_monomer_code(selected_hydroxy)
    print(f"[MS diag] using amino={selected_amino} (code={amino_code}, mw={amino_weight}) hydroxy={selected_hydroxy} (code={hydroxy_code}, mw={hydroxy_weight}) mode={ionization_mode}")

    # co-polymers of amino + hydroxy
    for n in range(1, 21):
        for m in range(1, 21):
            mass_val = n * amino_weight + m * hydroxy_weight - (n + m - 1) * 18.015
            if mass_val <= max_mass:
                condensation_products.append((mass_val + mass_adjustment, f"{n}x{amino_code} + {m}x{hydroxy_code}", "copoly"))

    # hydroxy-only oligomers
    for m in range(1, 21):
        mass_val = m * hydroxy_weight - (m - 1) * 18.015
        if mass_val <= max_mass:
            condensation_products.append((mass_val + mass_adjustment, f"{m}x{hydroxy_code}", "hydroxy"))

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
            color = 'red' if cat in ('amino', 'copoly') else 'blue'
        else:
            # Unlabeled peaks shown in light gray
            color = 'lightgray'
        
        # Create and store the stick line
        stick = ax.vlines(x=m, ymin=0, ymax=i, color=color, linewidth=1)
        stick_lines[idx] = stick
        
        # Add labels only for identified peaks
        if idx in labels_by_idx:
            label, cat = labels_by_idx[idx]
            color = 'red' if cat in ('amino', 'copoly') else 'blue'
            txt = ax.text(m, i + 1, label, rotation=90, color=color, fontsize=11, fontweight='bold', ha='center', va='bottom')
            line = Line2D([m, m], [i, i + 1], color=color, linestyle='--', linewidth=0.8)
            ax.add_line(line)
            draggables.append(DraggableTextWithLine(txt, line))

    ax.set_xlabel("m/z (Mass to Charge Ratio)", fontsize=13, fontweight='bold')
    ax.set_ylabel("Intensity (%)", fontsize=13, fontweight='bold')
    title_text = (custom_title or "").strip() or f"MS spectrum of {amino_code} & {hydroxy_code} ({ionization_mode.capitalize()})"
    ax.set_title(title_text, fontsize=14, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.4)

    # Legend for label colors
    legend_items = [
        {'label': 'Amino acid / depsipeptide', 'color': 'red', 'enabled': True, 'color_group': 'amino'},
        {'label': 'Hydroxy acid / polyester', 'color': 'blue', 'enabled': True, 'color_group': 'hydroxy'},
        {'label': 'Unidentified peaks', 'color': 'lightgray', 'enabled': True},
    ]
    legend_handles = [Line2D([0], [0], color=item['color'], lw=2, label=item['label']) for item in legend_items if item.get('enabled', True)]
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
        stick_lines,
        legend_items=legend_items,
        legend_config={'loc': 'upper right', 'frameon': True, 'framealpha': 0.85, 'fontsize': 9},
        synonym_map=build_default_synonym_map(),
        category_colors={'amino': 'red', 'hydroxy': 'blue', 'copoly': 'red'}
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
        win = getattr(manager, 'window', None)
        if win is not None:
            win.state('normal')  # bring to front
            win.attributes('-topmost', True)  # Make window top-most temporarily
            win.attributes('-topmost', False)  # Then allow normal z-order
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

def analyze_and_plot_multi_monomer(
    excel_path: str,
    sheet_name: str,
    chosen_display: str,
    selected_monomers: List[str],
    ionization_mode: str,
    min_intensity: float,
    tolerance: float = 0.5,
    max_mass: float = 1500,
    block_show: bool = True,
    custom_title: Optional[str] = None,
    custom_adduct_mass: Optional[float] = None,
    custom_adduct_name: Optional[str] = None,
):
    """
    Analyze and plot MS data with multiple monomers (not just binary amino+hydroxy).
    Generates all possible oligomer combinations from the selected monomers.
    """
    mass_adjustment = get_mass_adjustment(ionization_mode)
    
    # Load Data
    mass, intensity = load_ms_from_excel(excel_path, sheet_name)
    
    # Get molecular weights and codes for selected monomers
    monomer_weights = {}
    monomer_codes = {}
    for mon in selected_monomers:
        if mon not in molecular_weights:
            raise KeyError(f"Unknown monomer: {mon}. Available: {list(molecular_weights.keys())}")
        monomer_weights[mon] = molecular_weights[mon]
        monomer_codes[mon] = get_monomer_code(mon)
    
    print(f"[MS diag multi] sample='{chosen_display}' monomers={selected_monomers} mode={ionization_mode}")
    
    # Generate all possible oligomers
    condensation_products = []  # tuples: (mass, label, category)
    
    # Single monomers (homopolymers up to 20-mers)
    for monomer in selected_monomers:
        weight = monomer_weights[monomer]
        code = monomer_codes[monomer]
        for n in range(1, 21):
            mass_val = n * weight - (n - 1) * 18.015
            if mass_val <= max_mass:
                label = f"{n}x{code}" if n > 1 else code
                condensation_products.append((mass_val + mass_adjustment, label, monomer))
    
    # Co-polymers (all combinations of 2+ different monomers)
    if len(selected_monomers) >= 2:
        # Generate all combinations of oligomers with multiple monomer types
        for total_units in range(2, 21):  # Total oligomer size
            # For each way to distribute units among monomers
            from itertools import product
            for composition in product(range(total_units + 1), repeat=len(selected_monomers)):
                if sum(composition) == total_units and sum(1 for c in composition if c > 0) >= 2:
                    # Calculate mass and create label
                    mass_val = -18.015 * (total_units - 1)  # Condensation water loss
                    label_parts = []
                    category_parts = []
                    
                    for count, monomer in zip(composition, selected_monomers):
                        if count > 0:
                            mass_val += count * monomer_weights[monomer]
                            code = monomer_codes[monomer]
                            label_parts.append(f"{count}x{code}" if count > 1 else code)
                            category_parts.append(monomer)
                    
                    if mass_val <= max_mass:
                        label = " + ".join(label_parts)
                        category = "copoly"
                        condensation_products.append((mass_val + mass_adjustment, label, category))
    
    # Label Matching with de-duplication
    matches = []
    mass_list = list(mass)
    inten_list = list(intensity)
    
    for idx, (m, i) in enumerate(zip(mass_list, inten_list)):
        # Check if this is an isotope peak
        is_isotope = False
        for other_mass, other_intensity in zip(mass_list, inten_list):
            if other_mass > m and abs((other_mass - m) - 1.003) <= 0.2 and other_intensity > i:
                is_isotope = True
                break
        
        if is_isotope:
            continue
        
        # Find matching oligomers
        for calc_mass, label, cat in condensation_products:
            delta = abs(m - calc_mass)
            if delta <= tolerance:
                matches.append((idx, m, i, label, delta, cat))
    
    # De-duplicate: keep best match per peak
    matches.sort(key=lambda x: x[4])  # Sort by delta
    labels_by_idx = {}
    for idx, m, i, label, delta, cat in matches:
        if idx not in labels_by_idx:
            labels_by_idx[idx] = (label, cat)
    
    # Track which indices already have adducts to prevent stacking
    adduct_indices = set()
    
    # Check for non-covalent dimers
    dimer_matches = []
    for idx, (label, cat) in labels_by_idx.items():
        monomer_mass = mass[idx]
        expected_dimer = 2 * monomer_mass - mass_adjustment
        
        for dimer_idx, (m, i) in enumerate(zip(mass, intensity)):
            if dimer_idx not in labels_by_idx:
                delta = abs(m - expected_dimer)
                if delta <= tolerance:
                    dimer_label = f"{label} (2M{'−H' if ionization_mode == 'negative' else '+H'})"
                    dimer_matches.append((dimer_idx, dimer_label, cat, delta))
    
    for dimer_idx, dimer_label, cat, delta in dimer_matches:
        if dimer_idx not in labels_by_idx:
            labels_by_idx[dimer_idx] = (dimer_label, cat)
            adduct_indices.add(dimer_idx)  # Mark as having an adduct
    
    # Check for water adducts (polymer mass + 18) - only for oligomers without other adducts
    water_adduct_matches = []
    for idx, (label, cat) in labels_by_idx.items():
        # Skip if this peak already has an adduct (like 2M-H)
        if idx in adduct_indices:
            continue
        
        polymer_mass = mass[idx]
        # Expected water adduct: polymer mass + 18.015 (water)
        expected_water_adduct = polymer_mass + 18.015
        
        for adduct_idx, (m, i) in enumerate(zip(mass, intensity)):
            if adduct_idx not in labels_by_idx and adduct_idx not in adduct_indices:
                delta = abs(m - expected_water_adduct)
                if delta <= tolerance:
                    adduct_label = f"{label} + H₂O"
                    water_adduct_matches.append((adduct_idx, adduct_label, cat, delta))
    
    for adduct_idx, adduct_label, cat, delta in water_adduct_matches:
        if adduct_idx not in labels_by_idx:
            labels_by_idx[adduct_idx] = (adduct_label, cat)
            adduct_indices.add(adduct_idx)  # Mark as having an adduct
    
    # Check for custom adducts if provided
    if custom_adduct_mass is not None and custom_adduct_name:
        custom_adduct_matches = []
        for idx, (label, cat) in labels_by_idx.items():
            # Skip if this peak already has an adduct
            if idx in adduct_indices:
                continue
            
            polymer_mass = mass[idx]
            expected_custom_adduct = polymer_mass + custom_adduct_mass
            
            for adduct_idx, (m, i) in enumerate(zip(mass, intensity)):
                if adduct_idx not in labels_by_idx and adduct_idx not in adduct_indices:
                    delta = abs(m - expected_custom_adduct)
                    if delta <= tolerance:
                        adduct_label = f"{label} + {custom_adduct_name}"
                        custom_adduct_matches.append((adduct_idx, adduct_label, cat, delta))
        
        for adduct_idx, adduct_label, cat, delta in custom_adduct_matches:
            if adduct_idx not in labels_by_idx:
                labels_by_idx[adduct_idx] = (adduct_label, cat)
                adduct_indices.add(adduct_idx)  # Mark as having an adduct
    
    # Filter Data
    filtered = [(idx, m, i) for idx, (m, i) in enumerate(zip(mass, intensity)) if i >= min_intensity]
    if filtered:
        filtered_idx, filtered_mass, filtered_intensity = zip(*filtered)
    else:
        filtered_idx, filtered_mass, filtered_intensity = [], [], []
    
    # Plotting
    fig, ax = plt.subplots(figsize=(11, 6))
    
    if filtered_intensity:
        y_max = max(150, max(filtered_intensity) + 10)
    else:
        y_max = 150
    ax.set_ylim(0, y_max)
    
    if filtered_mass:
        x_min = max(0, min(filtered_mass) - 10)
        x_max = max(filtered_mass) + 10
        ax.set_xlim(x_min, x_max)
    else:
        ax.set_xlim(0, 1200)
    
    draggables = []
    stick_lines = {}
    
    # Color mapping by monomer identity (fixed): amino=red, hydroxy=blue
    color_map = {}
    for monomer in selected_monomers:
        color_map[monomer] = 'blue' if is_hydroxy_monomer(monomer) else 'red'
    color_map['copoly'] = 'black'
    
    # Plot all filtered masses
    for idx, m, i in zip(filtered_idx, filtered_mass, filtered_intensity):
        if idx in labels_by_idx:
            label, cat = labels_by_idx[idx]
            color = color_map.get(cat, 'black')
        else:
            label = None
            color = 'lightgray'
        
        line, = ax.plot([m, m], [0, i], color=color, linewidth=1.5, alpha=0.8)
        stick_lines[idx] = line
        
        if label:
            txt = ax.text(m, i + 2, label, ha='center', va='bottom', fontsize=10,
                         rotation=90, color=color, picker=True, weight='bold')
            connector = ax.plot([m, m], [i, i + 2], color=color, linewidth=0.5, alpha=0.5)[0]
            draggables.append(DraggableTextWithLine(txt, connector))
    
    # Create legend
    legend_elements = []
    legend_items = []
    has_amino = any(not is_hydroxy_monomer(m) for m in selected_monomers)
    has_hydroxy = any(is_hydroxy_monomer(m) for m in selected_monomers)

    if has_amino:
        legend_elements.append(Line2D([0], [0], color='red', linewidth=2, label='Amino acid oligomers'))
        legend_items.append({'label': 'Amino acid oligomers', 'color': 'red', 'enabled': True, 'color_group': 'amino'})

    if has_hydroxy:
        legend_elements.append(Line2D([0], [0], color='blue', linewidth=2, label='Hydroxy acid oligomers'))
        legend_items.append({'label': 'Hydroxy acid oligomers', 'color': 'blue', 'enabled': True, 'color_group': 'hydroxy'})
    
    if len(selected_monomers) >= 2:
        legend_elements.append(Line2D([0], [0], color='black', linewidth=2, label='Co-polymers'))
        legend_items.append({'label': 'Co-polymers', 'color': 'black', 'enabled': True, 'color_group': 'copoly'})
    
    legend_elements.append(Line2D([0], [0], color='lightgray', linewidth=2, label='Unidentified'))
    legend_items.append({'label': 'Unidentified', 'color': 'lightgray', 'enabled': True})
    
    ax.legend(handles=legend_elements, loc='upper right', fontsize=9, framealpha=0.9)
    
    # Set title and labels with bigger, bold font
    if custom_title:
        title = custom_title
    else:
        monomer_str = ", ".join([monomer_codes[m] for m in selected_monomers])
        title = f"{chosen_display} - Monomers: {monomer_str}"
    ax.set_title(title, fontsize=14, weight='bold')
    ax.set_xlabel("m/z", fontsize=13, weight='bold')
    ax.set_ylabel("Intensity (%)", fontsize=13, weight='bold')
    ax.grid(True, linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    
    # Create interactive editor for titles, labels, and annotations
    editor = InteractivePlotEditor(
        fig, ax, draggables,
        (filtered_idx, filtered_mass, filtered_intensity),
        labels_by_idx,
        ionization_mode,
        mass_adjustment,
        tolerance,
        stick_lines,
        legend_items=legend_items,
        legend_config={'loc': 'upper right', 'frameon': True, 'framealpha': 0.9, 'fontsize': 9},
        synonym_map=build_default_synonym_map(),
        category_colors={'amino': 'red', 'hydroxy': 'blue', 'copoly': 'black'}
    )
    
    # Keep references so event callbacks remain alive
    try:
        fig._draggables = draggables  # type: ignore[attr-defined]
        fig._editor = editor  # type: ignore[attr-defined]
    except Exception:
        pass
    
    # Ensure the figure manager is created and window is interactive
    manager = plt.get_current_fig_manager()
    try:
        win = getattr(manager, 'window', None)
        if win is not None:
            win.state('normal')  # bring to front
            win.attributes('-topmost', True)  # Make window top-most temporarily
            win.attributes('-topmost', False)  # Then allow normal z-order
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

    # --- Choose analysis mode: binary (amino+hydroxy) or multi-monomer ---
    use_multi_monomer = prompt_yes_no("Use multi-monomer mode? (Yes = select multiple monomers, No = traditional amino+hydroxy binary mode)", default=False)
    
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

    # Ask about custom adduct (only for multi-monomer mode)
    custom_adduct_mass = None
    custom_adduct_name = None
    if use_multi_monomer:
        add_custom_adduct = prompt_yes_no("Add a custom adduct to search for?", default=False)
        if add_custom_adduct:
            custom_adduct_mass = prompt_float("Enter the adduct mass (Da, always positive number):", default=None, title="Custom Adduct Mass")
            if custom_adduct_mass is not None:
                # Ask if adding or subtracting
                add_or_subtract = prompt_yes_no("Add this mass? (Yes = add mass like +Na, No = subtract mass like -CO₂)", default=True)
                if not add_or_subtract:
                    custom_adduct_mass = -custom_adduct_mass  # Make it negative for subtraction
                
                custom_adduct_name = prompt_text("Enter a short name for this adduct (e.g., 'Na', 'K', '-CO₂'):", default="Adduct", title="Adduct Name")
                if not custom_adduct_name:
                    custom_adduct_name = "Adduct"

    amino_acids = {k: v for k, v in molecular_weights.items() if k not in ["Glycolic Acid", "Lactic Acid", "Phenyllactic Acid", "hydroxy propionic acid", "hydroxybutyric acid"]}
    hydroxy_acids = {k: v for k, v in molecular_weights.items() if k in ["Glycolic Acid", "Lactic Acid", "Phenyllactic Acid", "hydroxy propionic acid", "hydroxybutyric acid"]}
    all_monomers = list(molecular_weights.keys())

    if batch:
        # Batch mode: analyze all sheets
        if use_multi_monomer:
            # Multi-monomer mode
            for disp, sname in sheet_items:
                while True:
                    selected_monomers = gui_select_multiple_monomers(all_monomers, disp)
                    if not selected_monomers:
                        # Skip this sheet
                        break
                    print(f"[MS diag multi] sheet='{sname}' display='{disp}' monomers={selected_monomers}")
                    custom_title = prompt_text("Optional: enter a custom plot title (leave blank for default)", default=disp, title="Plot title")
                    analyze_and_plot_multi_monomer(
                        excel_path,
                        sname,
                        disp,
                        selected_monomers,
                        ionization_mode,
                        min_intensity,
                        block_show=True,
                        custom_title=custom_title,
                        custom_adduct_mass=custom_adduct_mass,
                        custom_adduct_name=custom_adduct_name,
                    )
                    break
        else:
            # Binary mode (traditional amino+hydroxy)
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
        # Single-sample mode
        display_list = [d for d, _ in sheet_items]
        chosen_display = gui_select_from_list(display_list, title="Select sample", message="Select a sample") or display_list[0]
        sheet_idx = display_list.index(chosen_display)
        sheet_name = sheet_items[sheet_idx][1]

        if use_multi_monomer:
            # Multi-monomer mode
            selected_monomers = gui_select_multiple_monomers(all_monomers, chosen_display)
            if selected_monomers:
                custom_title = prompt_text("Optional: enter a custom plot title (leave blank for default)", default=chosen_display, title="Plot title")
                analyze_and_plot_multi_monomer(
                    excel_path,
                    sheet_name,
                    chosen_display,
                    selected_monomers,
                    ionization_mode,
                    min_intensity,
                    block_show=True,
                    custom_title=custom_title,
                    custom_adduct_mass=custom_adduct_mass,
                    custom_adduct_name=custom_adduct_name,
                )
        else:
            # Binary mode (traditional amino+hydroxy)
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
