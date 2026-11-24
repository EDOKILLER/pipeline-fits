#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PIPELINE FITS - dall'immagine grezza all'analisi scientifica
============================================================
Programma con finestra per Ubuntu/Linux e Windows 11, per i FITS di
telescopi spaziali e a terra (UV, visibile, infrarosso, raggi X, ...).

Sezioni (collassabili), una per ogni passo:
  1. Cartella e file      cartella madre, rinomina dei FITS (f090.fits ...), memoria
  2. Riduzione rumore     anteprime affiancate, scelta con un clic, per ogni file
  3. Sottrazione fondo    idem, con la riduzione del rumore scelta
  4. Stretching           idem; ricostruzione zone rovinate o scarto del file
  5. Stack bianco e nero  riproiezione WCS e combinazione dei filtri
  6. Stack RGB            pesi R/G/B per ogni filtro, anteprima in tempo reale
  7. Analisi scientifiche griglia WCS, SNR, gradiente, strutture, isofote,
                          sorgenti, istogramma, profilo, oggetti con nome (SIMBAD)

Tutto viene salvato nella cartella madre (scelte: pipeline_stato.json,
log: pipeline_log.txt, risultati: stack_*, analisi_*).

Opzioni: --autotest (prove complete su dati sintetici) e --diagnostica (SIMBAD).
"""

import json
import math
import os
import queue
import sys
import threading
import time
import traceback
from datetime import datetime
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.patches import Rectangle

import tkinter as tk
import tkinter.font as tkfont
from tkinter import ttk, filedialog, messagebox

import motore_pipeline as M

NOME = "Pipeline FITS"
TITOLO = f"{NOME} {M.VERSIONE}  —  {M.COPYRIGHT}"
FILE_CONFIG = Path.home() / ".pipeline_fits.json"
ZOOM_PREDEFINITO = 1.1 if sys.platform == "win32" else 1.5
ZOOM_MIN, ZOOM_MAX, ZOOM_PASSO = 0.8, 2.8, 0.1
MAX_LATO_PANNELLO = 1400
MODI_MEMORIA = {"Ibrida (consigliata)": "ibrido", "Solo RAM (più veloce)": "ram",
                "Solo disco (poca RAM)": "disco"}
METODI_STACK = {"Media pesata": "media", "Mediana": "mediana", "Massimo": "massimo"}
GRIGLIE = {"Griglia del file di riferimento": "riferimento",
           "Mosaico di tutti i campi, orientato come le immagini": "mosaico"}


def carica_config():
    try:
        c = json.loads(FILE_CONFIG.read_text(encoding="utf-8"))
        return c if isinstance(c, dict) else {}
    except (OSError, ValueError):
        return {}


def salva_config(c):
    try:
        FILE_CONFIG.write_text(json.dumps(c, indent=1, ensure_ascii=False), encoding="utf-8")
    except OSError:
        pass


def fmt_valore(v):
    return "—" if v is None else f"{float(v):g}"


# =============================================================================
# WIDGET GRANDI E CHIARI
# =============================================================================
class CasellaGrande(ttk.Button):
    """Casella sì/no grande: ☑ / ☐ seguiti dal testo, si clicca su tutta la riga."""

    def __init__(self, master, testo, valore=False, comando=None, **kw):
        self.var = tk.BooleanVar(value=valore)
        self.testo, self.comando = testo, comando
        super().__init__(master, style="Check.TButton", command=self._cambia, **kw)
        self._aggiorna()

    def _cambia(self):
        self.var.set(not self.var.get())
        self._aggiorna()
        if self.comando:
            self.comando(self.var.get())

    def _aggiorna(self):
        self.configure(text=("☑  " if self.var.get() else "☐  ") + self.testo)

    def get(self):
        return self.var.get()

    def set(self, v):
        self.var.set(bool(v))
        self._aggiorna()


class NumeroGrande(ttk.Frame):
    """Campo numerico con pulsanti − e + grandi al posto delle freccette."""

    def __init__(self, master, valore, passo=1, minimo=None, massimo=None, intero=False,
                 larghezza=7, comando=None):
        super().__init__(master)
        self.passo, self.minimo, self.massimo, self.intero = passo, minimo, massimo, intero
        self.comando = comando
        self.var = tk.StringVar(value=self._fmt(valore))
        ttk.Button(self, text="−", width=3, style="Piccolo.TButton", command=lambda: self._incr(-1)).pack(side=tk.LEFT)
        e = ttk.Entry(self, textvariable=self.var, width=larghezza, justify="center")
        e.pack(side=tk.LEFT, padx=3, fill=tk.Y)
        e.bind("<Return>", lambda ev: self._notifica())
        e.bind("<FocusOut>", lambda ev: self._notifica())
        ttk.Button(self, text="+", width=3, style="Piccolo.TButton", command=lambda: self._incr(+1)).pack(side=tk.LEFT)

    def _fmt(self, v):
        return str(int(round(v))) if self.intero else f"{v:g}"

    def _limita(self, v):
        if self.minimo is not None:
            v = max(self.minimo, v)
        if self.massimo is not None:
            v = min(self.massimo, v)
        return int(round(v)) if self.intero else v

    def get(self):
        try:
            return self._limita(float(self.var.get().replace(",", ".")))
        except ValueError:
            return self._limita(self.minimo or 0)

    def set(self, v):
        self.var.set(self._fmt(self._limita(v)))

    def _incr(self, s):
        self.set(round(self.get() + s * self.passo, 6))
        self._notifica()

    def _notifica(self):
        self.set(self.get())
        if self.comando:
            self.comando(self.get())


class EditorValori(ttk.Frame):
    """Lista di valori modificabile: aggiungi, modifica, rimuovi, predefiniti."""

    def __init__(self, master, titolo, comando=None):
        super().__init__(master)
        self.comando = comando
        self.intero, self.minimo, self.predefiniti = False, None, []
        self.lbl = ttk.Label(self, text=titolo, style="Campo.TLabel")
        self.lbl.grid(row=0, column=0, columnspan=3, sticky="w")
        self.lista = tk.Listbox(self, height=5, exportselection=False, font="TkDefaultFont",
                                activestyle="none", selectbackground="#3d7bd9")
        self.lista.grid(row=1, column=0, rowspan=4, sticky="nsew")
        self.lista.bind("<<ListboxSelect>>", self._seleziona)
        self.var = tk.StringVar()
        self.entry = ttk.Entry(self, textvariable=self.var, width=9, justify="center")
        self.entry.grid(row=1, column=1, columnspan=2, sticky="ew", padx=(6, 0))
        self.entry.bind("<Return>", lambda e: self._aggiungi())
        self.btn = []
        for i, (t, f) in enumerate((("➕ Aggiungi", self._aggiungi), ("✏ Modifica", self._modifica),
                                    ("🗑 Rimuovi", self._rimuovi), ("↺ Predefiniti", self._predef))):
            b = ttk.Button(self, text=t, command=f, style="Piccolo.TButton")
            b.grid(row=2 + i // 2, column=1 + i % 2, sticky="ew", padx=(6 if i % 2 == 0 else 3, 0), pady=2)
            self.btn.append(b)
        self.columnconfigure(0, weight=1)
        self.nota = ttk.Label(self, text="", style="Nota.TLabel", wraplength=300)
        self.nota.grid(row=5, column=0, columnspan=3, sticky="w")

    def configura(self, valori, predefiniti, intero=False, minimo=None, abilitato=True, nota=""):
        self.intero, self.minimo, self.predefiniti = intero, minimo, list(predefiniti)
        self.lista.delete(0, tk.END)
        for v in valori:
            self.lista.insert(tk.END, fmt_valore(v))
        stato = tk.NORMAL if abilitato else tk.DISABLED
        for w in self.btn + [self.entry]:
            w.configure(state=stato)
        self.lista.configure(state=stato)
        self.nota.configure(text=nota)

    def valori(self):
        out = []
        for t in self.lista.get(0, tk.END):
            out.append(None if t == "—" else (int(float(t)) if self.intero else float(t)))
        return out

    def _leggi(self):
        try:
            v = float(self.var.get().replace(",", "."))
        except ValueError:
            messagebox.showerror("Valore non valido", "Scrivi un numero (decimali con il punto o la virgola).")
            return None
        if self.minimo is not None and v < self.minimo:
            messagebox.showerror("Valore non valido", f"Il valore minimo è {self.minimo}.")
            return None
        return int(round(v)) if self.intero else v

    def _imposta(self, valori):
        valori = sorted(set(valori), key=lambda x: (x is None, x or 0))
        self.lista.delete(0, tk.END)
        for v in valori:
            self.lista.insert(tk.END, fmt_valore(v))
        if self.comando:
            self.comando(self.valori())

    def _seleziona(self, _e):
        s = self.lista.curselection()
        if s:
            self.var.set(self.lista.get(s[0]))

    def _aggiungi(self):
        v = self._leggi()
        if v is not None:
            self._imposta(self.valori() + [v])

    def _modifica(self):
        s = self.lista.curselection()
        v = self._leggi()
        if not s:
            messagebox.showinfo("Modifica", "Seleziona prima il valore da modificare nella lista.")
            return
        if v is not None:
            vals = self.valori()
            vals[s[0]] = v
            self._imposta(vals)

    def _rimuovi(self):
        s = self.lista.curselection()
        if s:
            vals = self.valori()
            vals.pop(s[0])
            if not vals:
                messagebox.showinfo("Rimuovi", "Deve restare almeno un valore.")
                return
            self._imposta(vals)

    def _predef(self):
        self._imposta(self.predefiniti)


class Sezione(ttk.Frame):
    """Sezione collassabile con titolo grande e stato a destra."""

    def __init__(self, master, titolo, al_apri):
        super().__init__(master)
        self.titolo, self.al_apri = titolo, al_apri
        self.aperta = False
        self.testa = ttk.Button(self, style="Sezione.TButton", command=self.commuta)
        self.testa.pack(fill=tk.X)
        self.corpo = ttk.Frame(self, padding=(12, 6, 6, 10))
        self.stato = ""
        self._aggiorna()

    def _aggiorna(self):
        freccia = "▼" if self.aperta else "▶"
        self.testa.configure(text=f"{freccia}  {self.titolo}" + (f"     [{self.stato}]" if self.stato else ""))

    def commuta(self):
        self.apri() if not self.aperta else self.chiudi()

    def apri(self):
        if not self.aperta:
            self.aperta = True
            self.corpo.pack(fill=tk.X)
            self._aggiorna()
            self.al_apri(self)

    def chiudi(self):
        if self.aperta:
            self.aperta = False
            self.corpo.pack_forget()
            self._aggiorna()

    def imposta_stato(self, testo):
        self.stato = testo
        self._aggiorna()


class PannelloScorrevole(ttk.Frame):
    """Contenitore verticale con barra di scorrimento grande e rotellina del mouse."""

    def __init__(self, master):
        super().__init__(master)
        self.canvas = tk.Canvas(self, highlightthickness=0, borderwidth=0)
        self.sb = ttk.Scrollbar(self, orient=tk.VERTICAL, command=self.canvas.yview)
        self.interno = ttk.Frame(self.canvas)
        self.finestra = self.canvas.create_window((0, 0), window=self.interno, anchor="nw")
        self.canvas.configure(yscrollcommand=self.sb.set)
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.interno.bind("<Configure>", lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        self.canvas.bind("<Configure>", lambda e: self.canvas.itemconfigure(self.finestra, width=e.width))
        self.bind("<Enter>", lambda e: self._ruota(True))
        self.bind("<Leave>", lambda e: self._ruota(False))

    def _ruota(self, attiva):
        if attiva:
            self.bind_all("<MouseWheel>", lambda e: self.canvas.yview_scroll(int(-e.delta / 120) or (-1 if e.delta > 0 else 1), "units"))
            self.bind_all("<Button-4>", lambda e: self.canvas.yview_scroll(-3, "units"))
            self.bind_all("<Button-5>", lambda e: self.canvas.yview_scroll(3, "units"))
        else:
            for ev in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
                self.unbind_all(ev)


# =============================================================================
# VISORE DELLE IMMAGINI
# =============================================================================
class Visore(ttk.Frame):
    def __init__(self, master, app):
        super().__init__(master, padding=6)
        self.app = app
        self.lbl_titolo = ttk.Label(self, text="Scegli la cartella madre nella sezione 1.", style="Titolo.TLabel")
        self.lbl_titolo.pack(fill=tk.X)
        self.bind("<Configure>", lambda e: self.lbl_titolo.configure(wraplength=max(300, e.width - 30)))
        self.fig = Figure(layout="constrained")
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        barra = ttk.Frame(self)
        barra.pack(fill=tk.X)
        self.toolbar = NavigationToolbar2Tk(self.canvas, barra, pack_toolbar=False)
        self.toolbar.update()
        self.toolbar.pack(side=tk.LEFT)
        ttk.Label(barra, text="  (lente = zoom, mano = sposta, casetta = vista intera; "
                              "disattivale per cliccare sulle immagini)", style="Nota.TLabel").pack(side=tk.LEFT)
        basso = ttk.Frame(self)
        basso.pack(side=tk.BOTTOM, fill=tk.X, pady=(6, 0))
        self.btn_prec = ttk.Button(basso, text="◀  Precedente", command=lambda: self._nav(-1), style="Grande.TButton")
        self.btn_prec.pack(side=tk.LEFT)
        self.btn_succ = ttk.Button(basso, text="Successivo  ▶", command=lambda: self._nav(+1), style="Grande.TButton")
        self.btn_succ.pack(side=tk.LEFT, padx=6)
        self.lbl_info = ttk.Label(basso, text="", style="Scelta.TLabel")
        self.lbl_info.pack(side=tk.LEFT, padx=10)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.canvas.mpl_connect("button_press_event", self._clic)
        self.nav = None
        self.clic = None
        self.assi = {}
        self.img_rgb = None
        self.selettore = None
        self.imposta_navigazione(None)

    def imposta_navigazione(self, funzione, etichette=("◀  Precedente", "Successivo  ▶")):
        self.nav = funzione
        stato = tk.NORMAL if funzione else tk.DISABLED
        self.btn_prec.configure(state=stato, text=etichette[0])
        self.btn_succ.configure(state=stato, text=etichette[1])

    def _nav(self, d):
        if self.nav:
            self.nav(d)

    def info(self, testo, colore="black"):
        self.lbl_info.configure(text=testo, foreground=colore)

    def _pulisci(self, titolo):
        if self.selettore is not None:
            try:
                self.selettore.set_active(False)
            except Exception:
                pass
            self.selettore = None
        self.fig.clear()
        self.assi, self.clic, self.img_rgb = {}, None, None
        self.lbl_titolo.configure(text=titolo)

    def messaggio(self, titolo, testo=""):
        self._pulisci(titolo)
        ax = self.fig.add_subplot(1, 1, 1)
        ax.axis("off")
        ax.text(0.5, 0.5, testo, ha="center", va="center", fontsize=12 * self.app.zoom, wrap=True,
                transform=ax.transAxes, color="#555")
        self.canvas.draw_idle()

    def griglia(self, pannelli, titolo, clic=None, scelto=None):
        """pannelli: lista di dict {img, titolo, dati}; clic(dati) quando si clicca."""
        self._pulisci(titolo)
        n = len(pannelli)
        ncol = 1 if n == 1 else 2 if n <= 4 else 3 if n <= 9 else 4 if n <= 16 else 5
        nrow = math.ceil(n / ncol)
        assi = self.fig.subplots(nrow, ncol, squeeze=False, sharex=True, sharey=True)
        fs = (10 if n <= 6 else 8.5) * self.app.zoom * 0.8
        for ax, p in zip(assi.flat, pannelli):
            img = p["img"]
            ax.imshow(np.nan_to_num(img, nan=0), origin="lower", cmap="gray", vmin=0, vmax=1,
                      interpolation="nearest")
            ax.set_title(p["titolo"], fontsize=fs)
            ax.set_xticks([]); ax.set_yticks([])
            self.assi[ax] = p
            if scelto is not None and p.get("dati") == scelto:
                self._bordo(ax, True)
        for ax in assi.flat[n:]:
            ax.axis("off")
        self.clic = clic
        self.canvas.draw_idle()

    def _bordo(self, ax_scelto, solo=False):
        if not solo:
            for ax in self.assi:
                for sp in ax.spines.values():
                    sp.set_edgecolor("black"); sp.set_linewidth(0.8)
        for sp in ax_scelto.spines.values():
            sp.set_edgecolor("red"); sp.set_linewidth(5)

    def evidenzia(self, dati):
        for ax, p in self.assi.items():
            if p.get("dati") == dati:
                self._bordo(ax)
        self.canvas.draw_idle()

    def immagine(self, img, titolo, clic=None, rgb=False, cmap="gray", estensione=None, zone=()):
        self._pulisci(titolo)
        ax = self.fig.add_subplot(1, 1, 1)
        kw = dict(origin="lower", interpolation="nearest")
        if estensione is not None:
            kw["extent"] = estensione
        if rgb:
            self.img_rgb = ax.imshow(np.clip(np.nan_to_num(img), 0, 1), **kw)
        else:
            ax.imshow(np.nan_to_num(img, nan=0), cmap=cmap, vmin=0, vmax=1, **kw)
        for (y0, y1, x0, x1, col) in zone:
            ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, edgecolor=col, linewidth=1.5))
        ax.set_xticks([]); ax.set_yticks([])
        self.assi = {ax: {"dati": None}}
        self.clic = clic
        self.canvas.draw_idle()
        return ax

    def aggiorna_rgb(self, img):
        if self.img_rgb is None:
            return False
        self.img_rgb.set_data(np.clip(np.nan_to_num(img), 0, 1))
        self.canvas.draw_idle()
        return True

    def figura_rgba(self, rgba, titolo):
        self._pulisci(titolo)
        ax = self.fig.add_subplot(1, 1, 1)
        ax.imshow(rgba, interpolation="antialiased")
        ax.axis("off")
        self.canvas.draw_idle()

    def rettangolo(self, img, titolo, estensione, zone, al_disegno):
        from matplotlib.widgets import RectangleSelector
        ax = self.immagine(img, titolo, estensione=estensione, zone=zone)
        self.clic = None

        def fatto(e1, e2):
            x0, x1 = sorted([e1.xdata, e2.xdata])
            y0, y1 = sorted([e1.ydata, e2.ydata])
            if x1 - x0 > 2 and y1 - y0 > 2:
                al_disegno([int(y0), int(math.ceil(y1)), int(x0), int(math.ceil(x1))])

        self.selettore = RectangleSelector(ax, fatto, useblit=False, button=[1], interactive=False,
                                           props=dict(edgecolor="yellow", fill=False, linewidth=2))

    def _clic(self, ev):
        if ev.inaxes is None or self.clic is None or ev.button != 1:
            return
        modo = getattr(self.toolbar.mode, "value", self.toolbar.mode)
        if modo:
            return
        p = self.assi.get(ev.inaxes)
        if p is None:
            return
        if p.get("dati") is not None or len(self.assi) > 1:
            self._bordo(ev.inaxes)
            self.canvas.draw_idle()
        self.clic(p.get("dati"), ev)


# =============================================================================
# APPLICAZIONE
# =============================================================================
class App(tk.Tk):
    FONT_SISTEMA = ("TkDefaultFont", "TkTextFont", "TkFixedFont", "TkMenuFont", "TkHeadingFont",
                    "TkCaptionFont", "TkSmallCaptionFont", "TkIconFont", "TkTooltipFont")

    def __init__(self):
        super().__init__()
        self.title(TITOLO)
        self.geometry("1440x880")
        self.minsize(1000, 640)
        try:
            self.state("zoomed")
        except tk.TclError:
            try:
                self.attributes("-zoomed", True)
            except tk.TclError:
                pass
        self.cfg = carica_config()
        self.zoom = min(max(float(self.cfg.get("zoom", ZOOM_PREDEFINITO)), ZOOM_MIN), ZOOM_MAX)
        self.coda = queue.Queue()
        self.occupato = False
        self.madre = None
        self.stato = None
        self.pipe = None
        self.nomi = []
        self._immagini = {}
        self.riproiettate = {}         # nome -> array sulla griglia dello stack
        self.griglia_stack = None      # (chiave, wcs, shape)
        self.stack = None              # (array, wcs)
        self.rgb_anteprime = {}
        self.pesi_rgb = {}
        self.risultati_analisi = []
        self.indice_analisi = 0
        self.sezioni = []
        self.sezione_attiva = None
        self.valori_nr = self.cfg.get("valori_nr", {k: v["valori"] for k, v in M.ALGORITMI_NR.items()})
        self.valori_fondo = self.cfg.get("valori_fondo", list(M.VALORI_FONDO))
        self.valori_st = self.cfg.get("valori_st", {k: v["valori"] for k, v in M.ALGORITMI_STRETCH.items()})

        self._stile()
        self._costruisci()
        self._applica_zoom()
        for t in ("<Control-plus>", "<Control-equal>", "<Control-KP_Add>"):
            self.bind_all(t, lambda e: self._cambia_zoom(+ZOOM_PASSO))
        for t in ("<Control-minus>", "<Control-KP_Subtract>"):
            self.bind_all(t, lambda e: self._cambia_zoom(-ZOOM_PASSO))
        self.bind_all("<Control-0>", lambda e: self._cambia_zoom(None))
        self.protocol("WM_DELETE_WINDOW", self._chiudi)
        self.after(100, self._controlla_coda)
        if self.cfg.get("madre") and Path(self.cfg["madre"]).is_dir():
            self.after(300, lambda: self._apri_madre(Path(self.cfg["madre"])))
        self.sezioni[0].apri()

    # ================================================================ stile
    def _stile(self):
        self._font_base = {}
        for n in self.FONT_SISTEMA:
            try:
                self._font_base[n] = tkfont.nametofont(n).cget("size") or 10
            except tk.TclError:
                pass
        fam = tkfont.nametofont("TkDefaultFont").actual()["family"]
        self.f_titolo = tkfont.Font(family=fam, weight="bold")
        self.f_sezione = tkfont.Font(family=fam, weight="bold")
        self.f_grassetto = tkfont.Font(family=fam, weight="bold")
        self.f_nota = tkfont.Font(family=fam, slant="italic")
        self.f_log = tkfont.Font(family=tkfont.nametofont("TkFixedFont").actual()["family"])
        self.st = ttk.Style(self)
        try:
            self.st.theme_use("clam")
        except tk.TclError:
            pass
        self.option_add("*TCombobox*Listbox.font", "TkDefaultFont")

    @staticmethod
    def _scala(d, f):
        n = max(1, round(abs(d) * f))
        return -n if d < 0 else n

    def _applica_zoom(self):
        z = self.zoom
        for n, b in self._font_base.items():
            tkfont.nametofont(n).configure(size=self._scala(b, z))
        b = self._font_base.get("TkDefaultFont", 10)
        self.f_titolo.configure(size=self._scala(b, 1.25 * z))
        self.f_sezione.configure(size=self._scala(b, 1.2 * z))
        self.f_grassetto.configure(size=self._scala(b, 1.05 * z))
        self.f_nota.configure(size=self._scala(b, 0.92 * z))
        self.f_log.configure(size=self._scala(self._font_base.get("TkFixedFont", b), 0.9 * z))
        p = lambda v: max(1, round(v * z))
        s = self.st
        s.configure("TButton", padding=(p(8), p(5)))
        s.configure("Grande.TButton", font=self.f_grassetto, padding=(p(14), p(8)))
        s.configure("Piccolo.TButton", padding=(p(6), p(3)))
        s.configure("Sezione.TButton", font=self.f_sezione, anchor="w", padding=(p(10), p(9)),
                    background="#dfe6ef")
        s.map("Sezione.TButton", background=[("active", "#c9d6e8")])
        s.configure("Check.TButton", anchor="w", padding=(p(6), p(4)), relief="flat")
        s.configure("Titolo.TLabel", font=self.f_titolo)
        s.configure("Campo.TLabel", font=self.f_grassetto)
        s.configure("Nota.TLabel", font=self.f_nota, foreground="#555555")
        s.configure("Scelta.TLabel", font=self.f_grassetto)
        s.configure("TEntry", padding=p(4))
        s.configure("TCombobox", padding=p(4), arrowsize=p(20))
        s.configure("Treeview", rowheight=tkfont.nametofont("TkDefaultFont").metrics("linespace") + p(8))
        s.configure("Treeview.Heading", font=self.f_grassetto, padding=p(4))
        for o in ("Vertical", "Horizontal"):
            s.configure(f"{o}.TScrollbar", arrowsize=p(22), width=p(22))
        s.configure("TProgressbar", thickness=p(18))
        self.lbl_zoom.configure(text=f"Testo {round(z * 100)}%")
        larg = tkfont.nametofont("TkDefaultFont").measure("0")
        for k, l in zip(("inc", "file", "filtro", "lambda", "nr", "bg", "st"), (3, 11, 8, 8, 7, 6, 10)):
            self.tab_file.column(k, width=l * larg + p(12))
        for sc in getattr(self, "_scale_rgb", []):
            sc.configure(width=p(18), sliderlength=p(28))
        for w in getattr(self, "_note_wrap", []):
            w.configure(wraplength=p(380))

    def _cambia_zoom(self, d):
        self.zoom = ZOOM_PREDEFINITO if d is None else min(max(round(self.zoom + d, 2), ZOOM_MIN), ZOOM_MAX)
        self._applica_zoom()
        self.cfg["zoom"] = self.zoom
        salva_config(self.cfg)
        return "break"

    def _nota(self, master, testo):
        l = ttk.Label(master, text=testo, style="Nota.TLabel", wraplength=380, justify="left")
        self._note_wrap = getattr(self, "_note_wrap", []) + [l]
        return l

    # ================================================================ layout
    def _costruisci(self):
        alto = ttk.Frame(self, padding=(8, 6))
        alto.pack(fill=tk.X)
        ttk.Button(alto, text="A−", width=4, command=lambda: self._cambia_zoom(-ZOOM_PASSO)).pack(side=tk.LEFT)
        ttk.Button(alto, text="A+", width=4, command=lambda: self._cambia_zoom(+ZOOM_PASSO)).pack(side=tk.LEFT, padx=4)
        self.lbl_zoom = ttk.Label(alto, text="")
        self.lbl_zoom.pack(side=tk.LEFT, padx=(4, 16))
        self.lbl_madre = ttk.Label(alto, text="Nessuna cartella madre", style="Campo.TLabel")
        self.lbl_madre.pack(side=tk.LEFT)
        ttk.Button(alto, text="ℹ Informazioni", command=self._informazioni).pack(side=tk.RIGHT, padx=(8, 0))
        self.chk_firma = CasellaGrande(alto, "Firma © sulle immagini salvate", self.cfg.get("firma", True),
                                       lambda v: self._salva_cfg("firma", v))
        self.chk_firma.pack(side=tk.RIGHT, padx=(8, 0))
        self.lbl_ram = ttk.Label(alto, text="", style="Nota.TLabel")
        self.lbl_ram.pack(side=tk.RIGHT)

        pw = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        pw.pack(fill=tk.BOTH, expand=True)
        sx = ttk.Frame(pw)
        pw.add(sx, weight=0)
        self.scorri = PannelloScorrevole(sx)
        self.scorri.pack(fill=tk.BOTH, expand=True)
        dx = ttk.Frame(pw)
        pw.add(dx, weight=1)

        # log e progresso sotto il visore
        bas = ttk.Frame(dx, padding=(6, 0, 6, 6))
        bas.pack(side=tk.BOTTOM, fill=tk.X)
        # filigrana: nome dell'autore, tenue, sotto il registro
        ttk.Label(bas, text=f"{NOME} {M.VERSIONE}  ·  {M.COPYRIGHT}", style="Nota.TLabel",
                  anchor="e").pack(side=tk.BOTTOM, fill=tk.X)
        self.lbl_prog = ttk.Label(bas, text="")
        self.lbl_prog.pack(fill=tk.X)
        self.prog = ttk.Progressbar(bas, mode="determinate", maximum=1000)
        self.prog.pack(fill=tk.X)
        fl = ttk.Frame(bas)
        fl.pack(fill=tk.X, pady=(4, 0))
        self.log = tk.Text(fl, height=4, wrap="word", font=self.f_log)
        sbl = ttk.Scrollbar(fl, orient=tk.VERTICAL, command=self.log.yview)
        self.log.configure(yscrollcommand=sbl.set)
        sbl.pack(side=tk.RIGHT, fill=tk.Y)
        self.log.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.visore = Visore(dx, self)
        self.visore.pack(fill=tk.BOTH, expand=True)

        c = self.scorri.interno
        self._sez_file(c)
        self._sez_rumore(c)
        self._sez_fondo(c)
        self._sez_stretch(c)
        self._sez_stack(c)
        self._sez_rgb(c)
        self._sez_analisi(c)
        self.after(2000, self._aggiorna_ram)

    def _nuova_sezione(self, master, titolo, al_apri):
        s = Sezione(master, titolo, lambda sez: self._apertura(sez, al_apri))
        s.pack(fill=tk.X, pady=(0, 4))
        self.sezioni.append(s)
        return s

    def _apertura(self, sez, al_apri):
        for s in self.sezioni:
            if s is not sez:
                s.chiudi()
        self.sezione_attiva = sez
        self.visore.imposta_navigazione(None)
        self.visore.info("")
        al_apri()
        self.after(50, lambda: self.scorri.canvas.yview_moveto(
            max(0, sez.winfo_y() / max(1, self.scorri.interno.winfo_height()))))

    def _navigatore(self, master, comando):
        f = ttk.Frame(master)
        ttk.Label(f, text="File:", style="Campo.TLabel").pack(side=tk.LEFT)
        ttk.Button(f, text="◀", width=3, command=lambda: comando(-1)).pack(side=tk.LEFT, padx=(6, 2))
        var = tk.StringVar()
        cb = ttk.Combobox(f, textvariable=var, state="readonly", width=14)
        cb.pack(side=tk.LEFT)
        cb.bind("<<ComboboxSelected>>", lambda e: comando(0))
        ttk.Button(f, text="▶", width=3, command=lambda: comando(+1)).pack(side=tk.LEFT, padx=2)
        return f, var, cb

    # ---------------------------------------------------------------- sezione 1
    def _sez_file(self, c):
        s = self._nuova_sezione(c, "1. Cartella madre e file FITS", self._apri_sez_file)
        b = s.corpo
        ttk.Button(b, text="📁  Scegli la cartella madre", style="Grande.TButton",
                   command=self._scegli_madre).pack(fill=tk.X)
        self._nota(b, "La cartella che contiene i FITS (anche dentro sottocartelle). "
                      "Tutti i risultati vengono salvati qui.").pack(fill=tk.X, pady=(2, 8))
        ttk.Button(b, text="🔤  Rinomina e copia i FITS delle sottocartelle",
                   command=self._rinomina).pack(fill=tk.X)
        self._nota(b, "Ogni file prende il nome del suo filtro (f090.fits, f1130.fits; doppioni: f090a, f090b) "
                      "e viene copiato nella cartella madre. Gli originali si cancellano solo dopo conferma."
                   ).pack(fill=tk.X, pady=(2, 8))
        ttk.Label(b, text="File della pipeline (clic = includi / escludi):", style="Campo.TLabel").pack(anchor="w")
        ft = ttk.Frame(b)
        ft.pack(fill=tk.X)
        col = ("inc", "file", "filtro", "lambda", "nr", "bg", "st")
        tit = ("✔", "File", "Filtro", "λ", "Rumore", "Fondo", "Stretch")
        lar = (3, 9, 8, 8, 7, 7, 9)
        self.tab_file = ttk.Treeview(ft, columns=col, show="headings", height=7, selectmode="browse")
        for k, t, l in zip(col, tit, lar):
            self.tab_file.heading(k, text=t)
            self.tab_file.column(k, width=l * 11, anchor="center" if k != "file" else "w", stretch=k == "file")
        sb = ttk.Scrollbar(ft, orient=tk.VERTICAL, command=self.tab_file.yview)
        self.tab_file.configure(yscrollcommand=sb.set)
        self.tab_file.pack(side=tk.LEFT, fill=tk.X, expand=True)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.tab_file.bind("<Button-1>", self._clic_tab_file)

        g = ttk.LabelFrame(b, text="Impostazioni", padding=8)
        g.pack(fill=tk.X, pady=8)
        ttk.Label(g, text="Uso della memoria:", style="Campo.TLabel").grid(row=0, column=0, sticky="w")
        self.var_mem = tk.StringVar(value=self.cfg.get("memoria", "Ibrida (consigliata)"))
        cb = ttk.Combobox(g, textvariable=self.var_mem, values=list(MODI_MEMORIA), state="readonly", width=22)
        cb.grid(row=0, column=1, sticky="w", padx=6)
        cb.bind("<<ComboboxSelected>>", lambda e: self._cambia_memoria())
        ttk.Label(g, text="Limite RAM (GB):", style="Campo.TLabel").grid(row=1, column=0, sticky="w", pady=4)
        self.num_limite = NumeroGrande(g, self.cfg.get("limite_gb", 4.0), passo=0.5, minimo=0.5, massimo=256,
                                       comando=lambda v: self._cambia_memoria())
        self.num_limite.grid(row=1, column=1, sticky="w", padx=6)
        self._nota(g, "Ibrida: usa la RAM fino al limite, poi il disco. Solo disco: per PC con poca RAM "
                      "(più lento). I file temporanei stanno in _pipeline_cache.").grid(row=2, column=0, columnspan=2, sticky="w")
        ttk.Label(g, text="Lato ritaglio anteprime (px):", style="Campo.TLabel").grid(row=3, column=0, sticky="w", pady=(8, 4))
        self.num_lato = NumeroGrande(g, self.cfg.get("lato", 1500), passo=250, minimo=250, massimo=8000, intero=True,
                                     comando=lambda v: self._salva_cfg("lato", v))
        self.num_lato.grid(row=3, column=1, sticky="w", padx=6, pady=(8, 4))
        ttk.Button(g, text="📍  Scegli la zona delle anteprime", command=self._scegli_zona).grid(
            row=4, column=0, columnspan=2, sticky="ew")
        self._nota(g, "Le anteprime usano un riquadro a piena risoluzione (veloce anche con file da 1 GB). "
                      "I risultati finali usano sempre l'immagine intera.").grid(row=5, column=0, columnspan=2, sticky="w")
        ttk.Label(g, text="Calcolo dell'SNR medio:", style="Campo.TLabel").grid(row=6, column=0, sticky="w", pady=(8, 0))
        self.var_snr = tk.StringVar(value=self.cfg.get("snr", M.METODI_SNR["segnale"]))
        ttk.Combobox(g, textvariable=self.var_snr, values=list(M.METODI_SNR.values()), state="readonly",
                     width=34).grid(row=7, column=0, columnspan=2, sticky="ew")
        self.chk_auto = CasellaGrande(g, "Dopo ogni scelta passa da solo al file successivo",
                                      self.cfg.get("auto", True), lambda v: self._salva_cfg("auto", v))
        self.chk_auto.grid(row=8, column=0, columnspan=2, sticky="ew", pady=(8, 0))
        self.sez_file = s

    # ---------------------------------------------------------------- sezione 2
    def _sez_rumore(self, c):
        s = self._nuova_sezione(c, "2. Riduzione del rumore", self._apri_sez_rumore)
        b = s.corpo
        ttk.Label(b, text="Algoritmo:", style="Campo.TLabel").pack(anchor="w")
        self.var_nr_algo = tk.StringVar(value=M.ALGORITMI_NR[self.cfg.get("nr_algo", "gaussiano")]["etichetta"])
        cb = ttk.Combobox(b, textvariable=self.var_nr_algo, state="readonly",
                          values=[v["etichetta"] for v in M.ALGORITMI_NR.values()])
        cb.pack(fill=tk.X)
        cb.bind("<<ComboboxSelected>>", lambda e: self._cambia_algo_nr())
        self.lbl_nr_descr = self._nota(b, "")
        self.lbl_nr_descr.pack(fill=tk.X, pady=(2, 6))
        self.ed_nr = EditorValori(b, "Valori da provare:", comando=self._salva_valori_nr)
        self.ed_nr.pack(fill=tk.X)
        ttk.Label(b, text="Stretching per valutare il risultato:", style="Campo.TLabel").pack(anchor="w", pady=(8, 0))
        f = ttk.Frame(b)
        f.pack(fill=tk.X)
        self.var_ev_algo = tk.StringVar(value=M.ALGORITMI_STRETCH[self.cfg.get("ev_algo", "asinh")]["etichetta"])
        ttk.Combobox(f, textvariable=self.var_ev_algo, state="readonly", width=16,
                     values=[v["etichetta"] for v in M.ALGORITMI_STRETCH.values()]).pack(side=tk.LEFT)
        ttk.Label(f, text="  valore:").pack(side=tk.LEFT)
        self.var_ev_val = tk.StringVar(value=str(self.cfg.get("ev_val", 0.1)))
        ttk.Entry(f, textvariable=self.var_ev_val, width=8, justify="center").pack(side=tk.LEFT)
        self._nota(b, "Oggetti diversi rendono meglio con stretching diversi (asinh per nebulose, "
                      "log per galassie...). Il valore è ignorato se lo stretching non ha parametri.").pack(fill=tk.X)
        nav, self.var_nr_file, self.cb_nr_file = self._navigatore(b, lambda d: self._vai_file("nr", d))
        nav.pack(fill=tk.X, pady=(10, 4))
        ttk.Button(b, text="▶  Calcola le anteprime", style="Grande.TButton",
                   command=lambda: self._anteprime("nr")).pack(fill=tk.X)
        self._nota(b, "La prima immagine ha solo lo stretching (riferimento). Clicca su quella migliore: "
                      "la scelta viene salvata.").pack(fill=tk.X, pady=(2, 0))
        self.sez_nr = s

    # ---------------------------------------------------------------- sezione 3
    def _sez_fondo(self, c):
        s = self._nuova_sezione(c, "3. Sottrazione del fondo", self._apri_sez_fondo)
        b = s.corpo
        ttk.Label(b, text="Stimatore del fondo (photutils Background2D):", style="Campo.TLabel").pack(anchor="w")
        self.var_bg_algo = tk.StringVar(value=M.ALGORITMI_FONDO[self.cfg.get("bg_algo", "mediana")]["etichetta"])
        cb = ttk.Combobox(b, textvariable=self.var_bg_algo, state="readonly",
                          values=[v["etichetta"] for v in M.ALGORITMI_FONDO.values()])
        cb.pack(fill=tk.X)
        cb.bind("<<ComboboxSelected>>", lambda e: self._descr_bg())
        self.lbl_bg_descr = self._nota(b, "")
        self.lbl_bg_descr.pack(fill=tk.X, pady=(2, 6))
        self.ed_bg = EditorValori(b, "Lato dei riquadri (box) in pixel:", comando=self._salva_valori_bg)
        self.ed_bg.pack(fill=tk.X)
        self._nota(b, "Box piccoli seguono meglio il fondo ma possono 'mangiare' la nebulosa; box grandi "
                      "danno un fondo più liscio. Stelle e nebulosa vengono escluse dalla stima.").pack(fill=tk.X)
        nav, self.var_bg_file, self.cb_bg_file = self._navigatore(b, lambda d: self._vai_file("bg", d))
        nav.pack(fill=tk.X, pady=(10, 4))
        ttk.Button(b, text="▶  Calcola le anteprime", style="Grande.TButton",
                   command=lambda: self._anteprime("bg")).pack(fill=tk.X)
        self._nota(b, "La prima immagine ha solo la riduzione del rumore: cliccala se NON vuoi sottrarre il fondo.").pack(fill=tk.X)
        self.sez_bg = s

    # ---------------------------------------------------------------- sezione 4
    def _sez_stretch(self, c):
        s = self._nuova_sezione(c, "4. Stretching e pixel rovinati", self._apri_sez_stretch)
        b = s.corpo
        ttk.Label(b, text="Tipo di stretching:", style="Campo.TLabel").pack(anchor="w")
        self.var_st_algo = tk.StringVar(value=self.cfg.get("st_algo_etichetta", "Tutti i tipi"))
        cb = ttk.Combobox(b, textvariable=self.var_st_algo, state="readonly",
                          values=["Tutti i tipi"] + [v["etichetta"] for v in M.ALGORITMI_STRETCH.values()])
        cb.pack(fill=tk.X)
        cb.bind("<<ComboboxSelected>>", lambda e: self._cambia_algo_st())
        self.ed_st = EditorValori(b, "Valori da provare:", comando=self._salva_valori_st)
        self.ed_st.pack(fill=tk.X, pady=(6, 0))
        self.chk_usa_bg = CasellaGrande(b, "Applica la sottrazione del fondo prima dello stretching", True,
                                        self._cambia_usa_bg)
        self.chk_usa_bg.pack(fill=tk.X, pady=(8, 0))
        self._nota(b, "Togli la spunta se la sottrazione del fondo rovina l'immagine (vale per il file corrente).").pack(fill=tk.X)
        rp = ttk.LabelFrame(b, text="Pixel rovinati", padding=6)
        rp.pack(fill=tk.X, pady=8)
        self.chk_ripara = CasellaGrande(rp, "🔧 Ricostruisci automaticamente le zone nere interne", False,
                                        self._cambia_ripara)
        self.chk_ripara.pack(fill=tk.X)
        ttk.Button(rp, text="▭  Disegna sull'immagine una zona da ricostruire", command=self._disegna_zona).pack(fill=tk.X, pady=(4, 0))
        ttk.Button(rp, text="Cancella le zone disegnate", command=self._cancella_zone).pack(fill=tk.X, pady=(4, 0))
        self.lbl_zone = self._nota(rp, "")
        self.lbl_zone.pack(fill=tk.X)
        self.chk_scarta = CasellaGrande(rp, "🚫 Scarta questo file dalle fasi successive", False, self._cambia_scarta)
        self.chk_scarta.pack(fill=tk.X, pady=(6, 0))
        nav, self.var_st_file, self.cb_st_file = self._navigatore(b, lambda d: self._vai_file("st", d))
        nav.pack(fill=tk.X, pady=(6, 4))
        ttk.Button(b, text="▶  Calcola le anteprime", style="Grande.TButton",
                   command=lambda: self._anteprime("st")).pack(fill=tk.X)
        ttk.Button(b, text="↩  Rivedi il fondo di questo file", command=self._torna_fondo).pack(fill=tk.X, pady=(4, 0))
        self.sez_st = s

    # ---------------------------------------------------------------- sezione 5
    def _sez_stack(self, c):
        s = self._nuova_sezione(c, "5. Stack in bianco e nero", self._apri_sez_stack)
        b = s.corpo
        ttk.Label(b, text="Griglia dell'immagine finale:", style="Campo.TLabel").pack(anchor="w")
        g = self.cfg.get("griglia")
        self.var_griglia = tk.StringVar(value=g if g in GRIGLIE else list(GRIGLIE)[0])
        cb = ttk.Combobox(b, textvariable=self.var_griglia, values=list(GRIGLIE), state="readonly")
        cb.pack(fill=tk.X)
        cb.bind("<<ComboboxSelected>>", lambda e: self._stima_griglia())
        ttk.Label(b, text="File di riferimento:", style="Campo.TLabel").pack(anchor="w", pady=(6, 0))
        self.var_rif = tk.StringVar()
        self.cb_rif = ttk.Combobox(b, textvariable=self.var_rif, state="readonly")
        self.cb_rif.pack(fill=tk.X)
        self.cb_rif.bind("<<ComboboxSelected>>", lambda e: (setattr(self, "_rif_utente", True), self._stima_griglia()))
        self._nota(b, "Proposto in automatico il file con i pixel più piccoli: è quello che dà la "
                      "risoluzione più alta.").pack(fill=tk.X)
        f = ttk.Frame(b)
        f.pack(fill=tk.X, pady=(6, 0))
        ttk.Label(f, text="Riduzione risoluzione:", style="Campo.TLabel").pack(side=tk.LEFT)
        self.num_rid = NumeroGrande(f, self.cfg.get("riduzione", 1), passo=1, minimo=1, massimo=32, intero=True,
                                    comando=lambda v: self._stima_griglia())
        self.num_rid.pack(side=tk.LEFT, padx=6)
        self.lbl_griglia = self._nota(b, "")
        self.lbl_griglia.pack(fill=tk.X)
        self._nota(b, "1 = piena risoluzione (massimo dettaglio). 2 = metà per lato (4 volte meno memoria "
                      "e tempo): usalo solo se manca memoria.").pack(fill=tk.X)
        ttk.Label(b, text="Come combinare i filtri:", style="Campo.TLabel").pack(anchor="w", pady=(6, 0))
        self.var_metodo = tk.StringVar(value=self.cfg.get("metodo", "Media pesata"))
        ttk.Combobox(b, textvariable=self.var_metodo, values=list(METODI_STACK), state="readonly").pack(fill=tk.X)
        ttk.Label(b, text="Peso di ogni filtro (media pesata):", style="Campo.TLabel").pack(anchor="w", pady=(6, 0))
        self.fr_pesi = ttk.Frame(b)
        self.fr_pesi.pack(fill=tk.X)
        self.num_pesi = {}
        ttk.Button(b, text="▶  Crea lo stack", style="Grande.TButton", command=self._crea_stack).pack(fill=tk.X, pady=(10, 4))
        rt = ttk.LabelFrame(b, text="Ritaglio dell'immagine finale", padding=6)
        rt.pack(fill=tk.X, pady=(6, 0))
        ttk.Button(rt, text="✂  Disegna il ritaglio sull'immagine", command=self._disegna_ritaglio).pack(fill=tk.X)
        self.chk_quadrato = CasellaGrande(rt, "Ritaglio quadrato", self.cfg.get("quadrato", False),
                                          lambda v: self._salva_cfg("quadrato", v))
        self.chk_quadrato.pack(fill=tk.X, pady=(4, 0))
        ttk.Button(rt, text="🪄  Ritaglio automatico: togli i bordi senza dati",
                   command=self._ritaglio_auto).pack(fill=tk.X, pady=(4, 0))
        ttk.Button(rt, text="↺  Annulla il ritaglio (immagine intera)", command=self._annulla_ritaglio).pack(fill=tk.X, pady=(4, 0))
        self.lbl_ritaglio = self._nota(rt, "Nessun ritaglio: si salva l'immagine intera.")
        self.lbl_ritaglio.pack(fill=tk.X)
        self._nota(rt, "Il ritaglio automatico tiene il rettangolo più grande in cui TUTTI i filtri hanno dati. "
                       "Il ritaglio vale anche per l'RGB e per le analisi; le coordinate WCS dei FITS "
                       "salvati vengono aggiornate automaticamente.").pack(fill=tk.X)
        ttk.Label(b, text="Formato dell'immagine:", style="Campo.TLabel").pack(anchor="w", pady=(6, 0))
        self.var_fmt = tk.StringVar(value=self.cfg.get("formato", "TIFF 16 bit"))
        ttk.Combobox(b, textvariable=self.var_fmt, values=list(M.FORMATI_IMMAGINE), state="readonly").pack(fill=tk.X)
        ttk.Button(b, text="💾  Salva stack (FITS + immagine)", style="Grande.TButton",
                   command=self._salva_stack).pack(fill=tk.X, pady=(10, 0))
        self._nota(b, "Salva stack_bn.fits (con coordinate WCS) e stack_bn nel formato scelto, nella cartella madre.").pack(fill=tk.X)
        self.sez_stack = s

    # ---------------------------------------------------------------- sezione 6
    def _sez_rgb(self, c):
        s = self._nuova_sezione(c, "6. Stack a colori RGB", self._apri_sez_rgb)
        b = s.corpo
        ttk.Button(b, text="▶  Prepara le immagini", command=self._prepara_rgb).pack(fill=tk.X)
        self._nota(b, "Usa la stessa griglia e le stesse scelte dello stack in bianco e nero.").pack(fill=tk.X)
        ttk.Button(b, text="🎨  Colorazione automatica (ordine cromatico)", style="Grande.TButton",
                   command=self._rgb_auto).pack(fill=tk.X, pady=(6, 0))
        self._nota(b, "Filtri ordinati per lunghezza d'onda: i più corti verso il blu, i più lunghi verso il rosso "
                      "(il metodo delle immagini JWST e Hubble).").pack(fill=tk.X)
        self.chk_bil = CasellaGrande(b, "Bilancia automaticamente i tre canali", True, lambda v: self._aggiorna_rgb())
        self.chk_bil.pack(fill=tk.X, pady=(4, 0))
        ttk.Label(b, text="Pesi rosso / verde / blu per ogni filtro:", style="Campo.TLabel").pack(anchor="w", pady=(8, 0))
        self.fr_rgb = ttk.Frame(b)
        self.fr_rgb.pack(fill=tk.X)
        self._scale_rgb = []
        ttk.Label(b, text="Formato dell'immagine:", style="Campo.TLabel").pack(anchor="w", pady=(6, 0))
        self.var_fmt_rgb = tk.StringVar(value=self.cfg.get("formato_rgb", "TIFF 16 bit"))
        ttk.Combobox(b, textvariable=self.var_fmt_rgb, values=list(M.FORMATI_IMMAGINE), state="readonly").pack(fill=tk.X)
        ttk.Button(b, text="💾  Salva RGB (FITS + immagine)", style="Grande.TButton",
                   command=self._salva_rgb).pack(fill=tk.X, pady=(10, 0))
        self.sez_rgb = s

    # ---------------------------------------------------------------- sezione 7
    def _sez_analisi(self, c):
        s = self._nuova_sezione(c, "7. Analisi scientifiche", self._apri_sez_analisi)
        b = s.corpo
        self._nota(b, "Usano stack_bn.fits della cartella madre (con le coordinate WCS).").pack(fill=tk.X)
        self.chk_analisi = {}
        scelte = self.cfg.get("analisi", {})
        for k, t in M.ANALISI.items():
            ch = CasellaGrande(b, t, scelte.get(k, True))
            ch.pack(fill=tk.X)
            self.chk_analisi[k] = ch
        f = ttk.Frame(b)
        f.pack(fill=tk.X, pady=(6, 0))
        ttk.Label(f, text="Nomi massimi sull'immagine:", style="Campo.TLabel").pack(side=tk.LEFT)
        self.num_nomi = NumeroGrande(f, self.cfg.get("max_nomi", 12), passo=1, minimo=1, massimo=60, intero=True)
        self.num_nomi.pack(side=tk.LEFT, padx=6)
        self._nota(b, "Vengono scelti i più famosi (più citati nella letteratura scientifica); "
                      "il CSV contiene tutti gli oggetti trovati. Serve Internet (SIMBAD).").pack(fill=tk.X)
        ttk.Label(b, text="Formato delle immagini:", style="Campo.TLabel").pack(anchor="w", pady=(6, 0))
        self.var_fmt_an = tk.StringVar(value=self.cfg.get("formato_an", "PNG 8 bit"))
        ttk.Combobox(b, textvariable=self.var_fmt_an, values=list(M.FORMATI_IMMAGINE), state="readonly").pack(fill=tk.X)
        ttk.Button(b, text="▶  Esegui le analisi", style="Grande.TButton", command=self._analisi).pack(fill=tk.X, pady=(10, 0))
        self.sez_an = s

    # ================================================================ utilità
    def _informazioni(self):
        """Finestra con nome del programma, versione e autore."""
        d = tk.Toplevel(self)
        d.title(f"Informazioni su {NOME}")
        d.transient(self)
        d.resizable(False, False)
        f = ttk.Frame(d, padding=24)
        f.pack(fill=tk.BOTH, expand=True)
        ttk.Label(f, text=f"{NOME}", style="Titolo.TLabel").pack(anchor="w")
        ttk.Label(f, text=f"Versione {M.VERSIONE}").pack(anchor="w", pady=(4, 0))
        ttk.Label(f, text=M.COPYRIGHT, style="Campo.TLabel").pack(anchor="w", pady=(12, 0))
        ttk.Label(f, text="Dall'immagine FITS grezza all'analisi scientifica: rumore, fondo, stretching,\n"
                          "stack in bianco e nero e RGB, analisi e oggetti con nome.", justify="left"
                  ).pack(anchor="w", pady=(12, 0))
        ttk.Label(f, text="Catalogo OpenNGC di Mattia Verga (licenza CC BY-SA 4.0).\n"
                          "Usa astropy, photutils, reproject, scikit-image, matplotlib.",
                  style="Nota.TLabel", justify="left").pack(anchor="w", pady=(12, 0))
        ttk.Button(f, text="Chiudi", command=d.destroy).pack(anchor="e", pady=(16, 0))
        d.grab_set()

    def _salva_cfg(self, k, v):
        self.cfg[k] = v
        salva_config(self.cfg)

    def scrivi(self, testo):
        self.log.insert(tk.END, f"{datetime.now():%H:%M:%S}  {testo}\n")
        self.log.see(tk.END)

    def _aggiorna_ram(self):
        g = M.ram_processo_gb()
        libera = M.ram_libera()
        if g == g:
            self.lbl_ram.configure(text=f"RAM programma {g:.1f} GB" +
                                   (f"  |  libera {libera / 1024**3:.1f} GB" if libera else ""))
        self.after(3000, self._aggiorna_ram)

    def esegui(self, descrizione, funzione, al_termine):
        """Calcolo in un thread separato: la finestra resta reattiva."""
        if self.occupato:
            messagebox.showinfo("Attendi", "C'è già un calcolo in corso.")
            return
        self.occupato = True
        self.lbl_prog.configure(text=descrizione + "...")
        self.prog.configure(mode="indeterminate")
        self.prog.start(15)
        log = lambda m: self.coda.put(("log", m))
        prog = lambda testo, k, t: self.coda.put(("prog", (testo, k, t)))

        def corpo():
            try:
                r = funzione(log, prog)
                self.coda.put(("fatto", (al_termine, r)))
            except MemoryError:
                self.coda.put(("errore", "Memoria esaurita: scegli 'Solo disco' o un limite più basso "
                                         "nella sezione 1, oppure aumenta la riduzione dello stack."))
            except Exception:
                self.coda.put(("errore", traceback.format_exc()))

        threading.Thread(target=corpo, daemon=True).start()

    def _controlla_coda(self):
        try:
            while True:
                tipo, d = self.coda.get_nowait()
                if tipo == "log":
                    self.scrivi(d)
                elif tipo == "prog":
                    testo, k, t = d
                    if str(self.prog.cget("mode")) != "determinate":
                        self.prog.stop()
                        self.prog.configure(mode="determinate")
                    self.prog.configure(value=1000 * k / max(t, 1))
                    self.lbl_prog.configure(text=f"{testo}: {100 * k / max(t, 1):.0f}%")
                elif tipo == "fatto":
                    self._fine_lavoro()
                    cb, r = d
                    try:
                        cb(r)
                    except Exception:
                        self.scrivi("ERRORE:\n" + traceback.format_exc())
                elif tipo == "errore":
                    self._fine_lavoro()
                    self.scrivi("ERRORE:\n" + d)
                    messagebox.showerror("Errore", d.strip().splitlines()[-1][:500])
        except queue.Empty:
            pass
        self.after(100, self._controlla_coda)

    def _fine_lavoro(self):
        self.occupato = False
        self.prog.stop()
        self.prog.configure(mode="determinate", value=0)
        self.lbl_prog.configure(text="")

    def immagine(self, nome):
        if nome not in self._immagini:
            self._immagini[nome] = M.ImmagineFITS(self.madre / nome)
        return self._immagini[nome]

    def scelte(self, nome):
        return self.stato.file(nome)

    def rumore(self, nome):
        s = self.scelte(nome)
        if not s.get("rumore"):
            s["rumore"] = list(M.stima_rumore(self.immagine(nome).dati))
            self.stato.salva()
        return tuple(s["rumore"])

    def centro(self):
        c = self.stato.globale("centro_anteprime") if self.stato else None
        return (c[0], c[1]) if c else (None, None)

    def metodo_snr(self):
        for k, v in M.METODI_SNR.items():
            if v == self.var_snr.get():
                return k
        return "segnale"

    def _chiave(self, diz, etichetta):
        for k, v in diz.items():
            if v["etichetta"] == etichetta:
                return k
        return next(iter(diz))

    def _richiedi_madre(self):
        if not self.madre or not self.nomi:
            messagebox.showinfo("Cartella madre", "Scegli prima la cartella madre con i FITS (sezione 1).")
            return False
        return True

    # ================================================================ sezione 1: azioni
    def _apri_sez_file(self):
        self.visore.messaggio("Cartella madre e file", "Scegli la cartella madre, poi rinomina i FITS se sono in "
                              "sottocartelle.\nClic su una riga della tabella = includi / escludi il file.")

    def _scegli_madre(self):
        if self.occupato:
            return
        d = filedialog.askdirectory(title="Cartella madre con i FITS",
                                    initialdir=self.cfg.get("madre") or str(Path.home()))
        if d:
            self._apri_madre(Path(d))

    def _apri_madre(self, p):
        for img in self._immagini.values():
            img.chiudi()
        self._immagini.clear()
        self.riproiettate.clear()
        self.griglia_stack = self.stack = None
        self.rgb_anteprime.clear()
        self.madre = p
        self.stato = M.Stato(p)
        self.pipe = M.Pipeline(p, M.Memoria(MODI_MEMORIA.get(self.var_mem.get(), "ibrido"), self.num_limite.get()))
        self._salva_cfg("madre", str(p))
        self.lbl_madre.configure(text=f"Cartella madre: {p}")
        self.scrivi(f"Cartella madre: {p}")
        self._aggiorna_file()
        if not self.nomi:
            sotto = [x for x in p.rglob("*") if x.is_file() and M.e_fits(x) and x.parent != p]
            if sotto:
                if messagebox.askyesno("FITS nelle sottocartelle",
                                       f"Nella cartella madre non ci sono FITS, ma ce ne sono {len(sotto)} nelle "
                                       f"sottocartelle.\nVuoi rinominarli e copiarli ora nella cartella madre?"):
                    self._rinomina()

    def _aggiorna_file(self):
        self.nomi = [x.name for x in M.elenca_file_pipeline(self.madre)] if self.madre else []
        self.tab_file.delete(*self.tab_file.get_children())
        for n in self.nomi:
            s = self.scelte(n)
            try:
                img = self.immagine(n)
                filtro, lam = img.filtro, M.descrivi_lunghezza(img.lambda_um)
                s.setdefault("filtro", filtro)
                s.setdefault("lambda_um", img.lambda_um)
            except Exception as e:
                filtro, lam = f"errore: {e}"[:20], ""
                s["incluso"] = False
            nr, bg, st = s.get("nr"), s.get("bg"), s.get("stretch")
            t_nr = f"{nr['valore']:g}" if nr and nr.get("valore") else ("no" if nr else "")
            t_bg = (f"{bg['box']}" if bg and bg.get("box") else ("no" if bg else ""))
            if not s.get("usa_bg", True):
                t_bg = "off"
            t_st = st["algo"] + (f" {st['valore']:g}" if st.get("valore") is not None else "") if st else ""
            self.tab_file.insert("", tk.END, iid=n, values=("☑" if s.get("incluso", True) else "☐", n, filtro,
                                                             lam, t_nr, t_bg, t_st))
        if self.stato:
            self.stato.salva()
        for cb in (self.cb_nr_file, self.cb_bg_file, self.cb_st_file):
            cb.configure(values=self.nomi)
        for var in (self.var_nr_file, self.var_bg_file, self.var_st_file):
            if var.get() not in self.nomi and self.nomi:
                var.set(self.nomi[0])
        inc = self.inclusi()
        self.cb_rif.configure(values=inc)
        if self.var_rif.get() not in inc and inc:
            self.var_rif.set(inc[0])
        n = len(self.nomi)
        self.sez_file.imposta_stato(f"{len(inc)}/{n} file inclusi" if n else "")
        for sez, k in ((self.sez_nr, "nr"), (self.sez_bg, "bg"), (self.sez_st, "stretch")):
            fatti = sum(1 for x in inc if self.scelte(x).get(k))
            sez.imposta_stato(f"{fatti}/{len(inc)} file scelti" if inc else "")

    def inclusi(self):
        return [n for n in self.nomi if self.scelte(n).get("incluso", True)]

    def _clic_tab_file(self, ev):
        if self.tab_file.identify_region(ev.x, ev.y) != "cell":
            return
        r = self.tab_file.identify_row(ev.y)
        if r:
            s = self.scelte(r)
            s["incluso"] = not s.get("incluso", True)
            self.stato.evento(f"{r}: {'incluso' if s['incluso'] else 'escluso'}")
            self._aggiorna_file()

    def _cambia_memoria(self):
        self._salva_cfg("memoria", self.var_mem.get())
        self._salva_cfg("limite_gb", self.num_limite.get())
        if self.pipe:
            self.pipe.mem.modo = MODI_MEMORIA.get(self.var_mem.get(), "ibrido")
            self.pipe.mem.limite = self.num_limite.get() * 1024**3
        self.scrivi(f"Memoria: {self.var_mem.get()}, limite {self.num_limite.get()} GB")

    def _rinomina(self):
        if not self.madre:
            messagebox.showinfo("Rinomina", "Scegli prima la cartella madre.")
            return
        self.esegui("Ricerca dei FITS nelle sottocartelle",
                    lambda log, prog: M.trova_fits_da_rinominare(self.madre, log), self._mostra_rinomina)

    def _mostra_rinomina(self, trovati):
        if not trovati:
            messagebox.showinfo("Rinomina", "Nelle sottocartelle non ci sono FITS da rinominare.")
            return
        d = tk.Toplevel(self)
        d.title("Rinomina e copia dei FITS")
        d.geometry("1100x600")
        d.transient(self)
        ttk.Label(d, text=f"{len(trovati)} file trovati. Controlla i nuovi nomi e premi 'Copia e rinomina'.",
                  style="Titolo.TLabel", padding=8).pack(fill=tk.X)
        f = ttk.Frame(d, padding=8)
        f.pack(fill=tk.BOTH, expand=True)
        col = ("orig", "filtro", "lam", "tel", "nuovo", "dim")
        tv = ttk.Treeview(f, columns=col, show="headings")
        for k, t, w in zip(col, ("File originale", "Filtro", "λ", "Telescopio / strumento", "Nuovo nome", "Dimensione"),
                           (420, 90, 90, 200, 120, 100)):
            tv.heading(k, text=t)
            tv.column(k, width=w, anchor="w")
        for t in trovati:
            tv.insert("", tk.END, values=(t["relativo"], t["filtro"], M.descrivi_lunghezza(t["lambda_um"]),
                                          f"{t['telescopio']} {t['strumento']}", t["nuovo_nome"],
                                          f"{t['dimensione'] / 1024**2:.0f} MB"))
        sb = ttk.Scrollbar(f, orient=tk.VERTICAL, command=tv.yview)
        tv.configure(yscrollcommand=sb.set)
        tv.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        totale = sum(t["dimensione"] for t in trovati)
        libero = M.shutil.disk_usage(self.madre).free
        bb = ttk.Frame(d, padding=8)
        bb.pack(fill=tk.X)
        ttk.Label(bb, text=f"Spazio necessario {totale / 1024**3:.1f} GB, libero {libero / 1024**3:.1f} GB").pack(side=tk.LEFT)

        def vai():
            d.destroy()
            if totale > libero * 0.98:
                messagebox.showerror("Spazio insufficiente", "Non c'è abbastanza spazio su disco per le copie.")
                return
            self.esegui("Copia e verifica dei file",
                        lambda log, prog: M.copia_e_verifica(trovati, self.madre,
                                                             lambda k, t, n: prog("Copia", k, t), log),
                        self._dopo_copia)

        ttk.Button(bb, text="Annulla", command=d.destroy).pack(side=tk.RIGHT)
        ttk.Button(bb, text="✔  Copia e rinomina", style="Grande.TButton", command=vai).pack(side=tk.RIGHT, padx=8)

    def _dopo_copia(self, copiati):
        self._aggiorna_file()
        if not copiati:
            return
        self.stato.evento(f"Rinominati e copiati {len(copiati)} file nella cartella madre")
        elenco = "\n".join(f"  {t['relativo']}" for t in copiati[:15])
        if len(copiati) > 15:
            elenco += f"\n  ... e altri {len(copiati) - 15}"
        if messagebox.askyesno("Cancellare gli originali?",
                               f"{len(copiati)} file copiati e verificati (copie identiche).\n\n"
                               f"Cancello gli originali nelle sottocartelle?\n\n{elenco}\n\n"
                               f"Le sottocartelle rimaste vuote verranno eliminate.", icon="warning"):
            M.cancella_originali(copiati, self.madre, self.scrivi)
            self.stato.evento(f"Cancellati {len(copiati)} file originali")
            self.scrivi(f"Cancellati {len(copiati)} file originali.")
        self.stato.salva()

    def _scegli_zona(self):
        if not self._richiedi_madre():
            return
        nome = self.var_nr_file.get() or self.nomi[0]

        def lavoro(log, prog):
            img = self.immagine(nome)
            rid, f = self.pipe.ridotta(img)
            return nome, M.applica_stretch(rid, "asinh", 0.1), img.shape

        def mostra(r):
            nome, ant, shape = r
            ny, nx = shape
            lato = self.num_lato.get()
            ra, dec = self.centro()
            img = self.immagine(nome)
            xc, yc = img.centro_pixel(ra, dec)
            zone = [(yc - lato / 2, yc + lato / 2, xc - lato / 2, xc + lato / 2, "yellow")]
            self.visore.immagine(ant, f"{nome}: clicca sul centro della zona da usare per le anteprime "
                                      f"(riquadro giallo = zona attuale)", clic=self._clic_zona,
                                 estensione=(0, nx, 0, ny), zone=zone)
            self._zona_nome = nome

        self.esegui("Preparazione della panoramica", lavoro, mostra)

    def _clic_zona(self, _dati, ev):
        img = self.immagine(self._zona_nome)
        x, y = ev.xdata, ev.ydata
        if img.wcs is not None:
            ra, dec = img.wcs.all_pix2world([[x, y]], 0)[0]
            self.stato.imposta_globale("centro_anteprime", [float(ra), float(dec)],
                                       f"Zona delle anteprime centrata su RA {ra:.5f}, Dec {dec:.5f}")
        self.scrivi(f"Zona delle anteprime: pixel ({x:.0f}, {y:.0f}) di {self._zona_nome}")
        self._scegli_zona()

    # ================================================================ sezioni 2-4: anteprime
    def _apri_sez_rumore(self):
        self._cambia_algo_nr(salva=False)
        self.visore.messaggio("Riduzione del rumore", "Scegli algoritmo e valori, poi premi 'Calcola le anteprime'.")
        self.visore.imposta_navigazione(lambda d: self._vai_file("nr", d), ("◀  File precedente", "File successivo  ▶"))

    def _cambia_algo_nr(self, salva=True):
        k = self._chiave(M.ALGORITMI_NR, self.var_nr_algo.get())
        a = M.ALGORITMI_NR[k]
        self.lbl_nr_descr.configure(text=a["descrizione"])
        self.ed_nr.lbl.configure(text=f"Valori da provare ({a['parametro']}):")
        self.ed_nr.configura(self.valori_nr.get(k, a["valori"]), a["valori"], a["intero"], a["minimo"])
        if salva:
            self._salva_cfg("nr_algo", k)

    def _salva_valori_nr(self, valori):
        k = self._chiave(M.ALGORITMI_NR, self.var_nr_algo.get())
        self.valori_nr[k] = valori
        self._salva_cfg("valori_nr", self.valori_nr)

    def _apri_sez_fondo(self):
        self._descr_bg()
        self.ed_bg.configura(self.valori_fondo, M.VALORI_FONDO, intero=True, minimo=8)
        self.visore.messaggio("Sottrazione del fondo", "Premi 'Calcola le anteprime' per il file selezionato.")
        self.visore.imposta_navigazione(lambda d: self._vai_file("bg", d), ("◀  File precedente", "File successivo  ▶"))

    def _descr_bg(self):
        k = self._chiave(M.ALGORITMI_FONDO, self.var_bg_algo.get())
        self.lbl_bg_descr.configure(text=M.ALGORITMI_FONDO[k]["descrizione"])
        self._salva_cfg("bg_algo", k)

    def _salva_valori_bg(self, valori):
        self.valori_fondo = valori
        self._salva_cfg("valori_fondo", valori)

    def _apri_sez_stretch(self):
        self._cambia_algo_st()
        self._mostra_stato_file_st()
        self.visore.messaggio("Stretching", "Premi 'Calcola le anteprime' per il file selezionato.")
        self.visore.imposta_navigazione(lambda d: self._vai_file("st", d), ("◀  File precedente", "File successivo  ▶"))

    def _algo_st(self):
        t = self.var_st_algo.get()
        return None if t == "Tutti i tipi" else self._chiave(M.ALGORITMI_STRETCH, t)

    def _cambia_algo_st(self):
        self._salva_cfg("st_algo_etichetta", self.var_st_algo.get())
        k = self._algo_st()
        if k is None:
            self.ed_st.configura([], [], abilitato=False,
                                 nota="Con 'Tutti i tipi' si usano i valori salvati per ciascuno stretching "
                                      "(scegli un tipo per modificarli).")
            return
        a = M.ALGORITMI_STRETCH[k]
        if a["parametro"] is None:
            self.ed_st.configura([None], [None], abilitato=False, nota="Questo stretching non ha parametri.")
        else:
            self.ed_st.lbl.configure(text=f"Valori da provare ({a['parametro']}):")
            self.ed_st.configura(self.valori_st.get(k, a["valori"]), a["valori"], minimo=1e-6)

    def _salva_valori_st(self, valori):
        k = self._algo_st()
        if k:
            self.valori_st[k] = valori
            self._salva_cfg("valori_st", self.valori_st)

    def _var_file(self, fase):
        return {"nr": self.var_nr_file, "bg": self.var_bg_file, "st": self.var_st_file}[fase]

    def _vai_file(self, fase, d):
        if not self.nomi:
            return
        var = self._var_file(fase)
        i = self.nomi.index(var.get()) if var.get() in self.nomi else 0
        i = (i + d) % len(self.nomi)
        var.set(self.nomi[i])
        if fase == "st":
            self._mostra_stato_file_st()
        self._anteprime(fase)

    def _mostra_stato_file_st(self):
        n = self.var_st_file.get()
        if not n or not self.stato:
            return
        s = self.scelte(n)
        bg = s.get("bg") or {}
        self.chk_usa_bg.set(s.get("usa_bg", bool(bg.get("box"))))
        rp = s.get("riparazione") or {}
        self.chk_ripara.set(rp.get("auto", False))
        self.chk_scarta.set(not s.get("incluso", True))
        t = []
        if rp.get("auto"):
            t.append(f"{len(rp.get('zone_auto', []))} zone trovate automaticamente")
        if rp.get("rettangoli"):
            t.append(f"{len(rp['rettangoli'])} zone disegnate")
        self.lbl_zone.configure(text="; ".join(t) or "Nessuna zona da ricostruire.")

    def _parametri_valutazione(self):
        k = self._chiave(M.ALGORITMI_STRETCH, self.var_ev_algo.get())
        try:
            v = float(self.var_ev_val.get().replace(",", "."))
        except ValueError:
            v = M.ALGORITMI_STRETCH[k]["valori"][0]
        if M.ALGORITMI_STRETCH[k]["parametro"] is None:
            v = None
        self._salva_cfg("ev_algo", k)
        if v is not None:
            self._salva_cfg("ev_val", v)
        return k, v

    @staticmethod
    def _decima(a):
        s = max(1, math.ceil(max(a.shape) / MAX_LATO_PANNELLO))
        return np.ascontiguousarray(a[::s, ::s])

    def _anteprime(self, fase):
        if not self._richiedi_madre():
            return
        nome = self._var_file(fase).get()
        if not nome:
            return
        s = self.scelte(nome)
        lato = self.num_lato.get()
        ra, dec = self.centro()
        metodo = self.metodo_snr()
        ev_algo, ev_val = self._parametri_valutazione()
        if fase in ("bg", "st") and not s.get("nr"):
            messagebox.showinfo("Manca la riduzione del rumore",
                                f"Per {nome} scegli prima la riduzione del rumore (sezione 2).")
            return
        if fase == "nr":
            algo = self._chiave(M.ALGORITMI_NR, self.var_nr_algo.get())
            valori = self.ed_nr.valori()
        elif fase == "bg":
            algo = self._chiave(M.ALGORITMI_FONDO, self.var_bg_algo.get())
            valori = self.ed_bg.valori()
        else:
            k = self._algo_st()
            combinazioni = []
            for a in ([k] if k else list(M.ALGORITMI_STRETCH)):
                vals = self.valori_st.get(a, M.ALGORITMI_STRETCH[a]["valori"]) if M.ALGORITMI_STRETCH[a]["parametro"] else [None]
                combinazioni += [(a, v) for v in vals]
        scelte = json.loads(json.dumps(s))

        def lavoro(log, prog):
            img = self.immagine(nome)
            med, sig = self.rumore(nome)
            crop, (y0, x0) = self.pipe.ritaglio_preparato(img, ra, dec, lato, scelte)
            pannelli = []
            st_ev = scelte.get("stretch_valutazione") or {"algo": ev_algo, "valore": ev_val}
            if fase == "nr":
                st_ev = {"algo": ev_algo, "valore": ev_val}
                lista = [None] + list(valori)
                masc = M.maschera_segnale(crop)
                for i, v in enumerate(lista):
                    prog("Anteprime", i, len(lista))
                    r = crop if v is None else M.riduci_rumore(crop, algo, v, med, sig)
                    snr = M.snr_medio(r, metodo, masc)
                    if v is None:
                        t = f"Solo stretching\nSNR medio = {snr:.2f}"
                        dati = {"algo": algo, "valore": None}
                    else:
                        t = f"{M.ALGORITMI_NR[algo]['etichetta']} {v:g}\nSNR medio = {snr:.2f}"
                        dati = {"algo": algo, "valore": v}
                    pannelli.append({"img": self._decima(M.applica_stretch(r, st_ev["algo"], st_ev["valore"])),
                                     "titolo": t, "dati": dati})
                return nome, pannelli, st_ev
            nr = scelte["nr"]
            base = M.riduci_rumore(crop, nr["algo"], nr["valore"], med, sig) if nr.get("valore") else crop
            if fase == "bg":
                lista = [None] + list(valori)
                masc = M.maschera_segnale(base)
                for i, box in enumerate(lista):
                    prog("Anteprime", i, len(lista))
                    if box is None:
                        r = base
                        t = "Fondo NON sottratto"
                    else:
                        log(f"{nome}: fondo con box {box}")
                        r = base - self.pipe.fondo_ritaglio(img, algo, box, y0, x0, base.shape, scelte)
                        t = f"Fondo {M.ALGORITMI_FONDO[algo]['etichetta']} box {box}"
                    snr = M.snr_medio(r, metodo, masc)
                    pannelli.append({"img": self._decima(M.applica_stretch(r, st_ev["algo"], st_ev["valore"])),
                                     "titolo": f"{t}\nSNR medio = {snr:.2f}", "dati": {"algo": algo, "box": box}})
                return nome, pannelli, None
            bg = scelte.get("bg") or {}
            usa = scelte.get("usa_bg", bool(bg.get("box")))
            if usa and bg.get("box"):
                base = base - self.pipe.fondo_ritaglio(img, bg["algo"], bg["box"], y0, x0, base.shape, scelte)
            limiti = M.limiti_zscale(base)
            for i, (a, v) in enumerate(combinazioni):
                prog("Anteprime", i, len(combinazioni))
                pannelli.append({"img": self._decima(M.applica_stretch(base, a, v, limiti)),
                                 "titolo": M.descrivi_stretch(a, v), "dati": {"algo": a, "valore": v}})
            return nome, pannelli, None

        self.esegui(f"Anteprime di {nome}", lavoro, lambda r: self._mostra_anteprime(fase, *r))

    def _mostra_anteprime(self, fase, nome, pannelli, st_ev):
        s = self.scelte(nome)
        chiave = {"nr": "nr", "bg": "bg", "st": "stretch"}[fase]
        scelto = s.get(chiave)
        if fase == "nr":
            self._st_ev_corrente = st_ev
        titoli = {"nr": "Riduzione del rumore", "bg": "Sottrazione del fondo", "st": "Stretching"}
        inc = "" if s.get("incluso", True) else "   (FILE SCARTATO)"
        self.visore.griglia(pannelli, f"{titoli[fase]} — {nome}{inc}: clicca sull'immagine migliore",
                            clic=lambda dati, ev: self._scelta(fase, nome, dati), scelto=scelto)
        self.visore.imposta_navigazione(lambda d: self._vai_file(fase, d), ("◀  File precedente", "File successivo  ▶"))
        self.visore.info(f"Scelta attuale: {self._testo_scelta(fase, scelto)}" if scelto else "Nessuna scelta per questo file")

    def _testo_scelta(self, fase, d):
        if not d:
            return "-"
        if fase == "nr":
            return "nessuna riduzione" if d.get("valore") is None else f"{M.ALGORITMI_NR[d['algo']]['etichetta']} {d['valore']:g}"
        if fase == "bg":
            return "nessuna sottrazione" if d.get("box") is None else f"{M.ALGORITMI_FONDO[d['algo']]['etichetta']} box {d['box']}"
        return M.descrivi_stretch(d["algo"], d.get("valore"))

    def _scelta(self, fase, nome, dati):
        chiave = {"nr": "nr", "bg": "bg", "st": "stretch"}[fase]
        self.stato.imposta(nome, chiave, dati, f"{ {'nr': 'riduzione rumore', 'bg': 'fondo', 'st': 'stretching'}[fase]} "
                                               f"= {self._testo_scelta(fase, dati)}")
        if fase == "nr":
            self.stato.imposta(nome, "stretch_valutazione", getattr(self, "_st_ev_corrente", None))
        if fase == "bg":
            self.stato.imposta(nome, "usa_bg", dati.get("box") is not None)
        if fase in ("nr", "bg"):
            self.scelte(nome).pop("cache_stretch", None)
        self.visore.info(f"✔ Salvato per {nome}: {self._testo_scelta(fase, dati)}", "#0a7a0a")
        self.scrivi(f"{nome}: scelta {self._testo_scelta(fase, dati)}")
        self._aggiorna_file()
        if self.chk_auto.get():
            inc = self.inclusi()
            chiave_fatto = chiave
            restanti = [n for n in inc if not self.scelte(n).get(chiave_fatto)]
            if restanti:
                self._var_file(fase).set(restanti[0])
                if fase == "st":
                    self._mostra_stato_file_st()
                self.after(400, lambda: self._anteprime(fase))
            else:
                successiva = {"nr": self.sez_bg, "bg": self.sez_st, "st": self.sez_stack}[fase]
                messagebox.showinfo("Fatto", "Tutti i file inclusi hanno una scelta. "
                                             "Passo alla sezione successiva.")
                successiva.apri()

    # ---------------------------------------------------------------- sezione 4: extra
    def _cambia_usa_bg(self, v):
        n = self.var_st_file.get()
        if n and self.stato:
            self.stato.imposta(n, "usa_bg", bool(v), f"fondo {'applicato' if v else 'NON applicato'} prima dello stretching")
            self._aggiorna_file()

    def _cambia_scarta(self, v):
        n = self.var_st_file.get()
        if n and self.stato:
            self.stato.imposta(n, "incluso", not v, "SCARTATO dalle fasi successive" if v else "di nuovo incluso")
            self._aggiorna_file()

    def _cambia_ripara(self, v):
        n = self.var_st_file.get()
        if not n or not self.stato:
            return
        rp = dict(self.scelte(n).get("riparazione") or {})
        rp["auto"] = bool(v)
        if not v:
            self.stato.imposta(n, "riparazione", rp, "ricostruzione automatica disattivata")
            self.pipe.dimentica_file(n)
            self._mostra_stato_file_st()
            return

        def lavoro(log, prog):
            return M.trova_zone_rovinate(self.immagine(n))

        def fatto(zone):
            rp["zone_auto"] = zone
            self.stato.imposta(n, "riparazione", rp, f"ricostruzione automatica: {len(zone)} zone")
            self.pipe.dimentica_file(n)
            self._mostra_stato_file_st()
            self.scrivi(f"{n}: trovate {len(zone)} zone rovinate interne")
            if zone:
                self._anteprime("st")
            else:
                messagebox.showinfo("Pixel rovinati", "Non ho trovato zone nere all'interno dell'immagine. "
                                                      "Se ne vedi qualcuna, disegnala con il rettangolo.")

        self.esegui(f"Ricerca delle zone rovinate in {n}", lavoro, fatto)

    def _disegna_zona(self):
        n = self.var_st_file.get()
        if not n or not self._richiedi_madre():
            return

        def lavoro(log, prog):
            img = self.immagine(n)
            rid, f = self.pipe.ridotta(img)
            return M.applica_stretch(rid, "asinh", 0.1), img.shape

        def mostra(r):
            ant, (ny, nx) = r
            rp = self.scelte(n).get("riparazione") or {}
            zone = [(z[0], z[1], z[2], z[3], "red") for z in rp.get("rettangoli", [])]
            zone += [(z[0], z[1], z[2], z[3], "orange") for z in (rp.get("zone_auto", []) if rp.get("auto") else [])]
            self.visore.rettangolo(ant, f"{n}: trascina con il mouse un rettangolo sulla zona rovinata "
                                        f"(rosso = zone disegnate, arancione = trovate in automatico)",
                                   (0, nx, 0, ny), zone, lambda rett: self._aggiungi_zona(n, rett))
            self.visore.info("Puoi disegnare più zone; poi premi 'Calcola le anteprime'.")

        self.esegui("Preparazione della panoramica", lavoro, mostra)

    def _aggiungi_zona(self, n, rett):
        rp = dict(self.scelte(n).get("riparazione") or {})
        rp.setdefault("rettangoli", []).append(rett)
        self.stato.imposta(n, "riparazione", rp, f"zona da ricostruire disegnata {rett}")
        self.pipe.dimentica_file(n)
        self._mostra_stato_file_st()
        self.scrivi(f"{n}: aggiunta zona y {rett[0]}-{rett[1]}, x {rett[2]}-{rett[3]}")
        self._disegna_zona()

    def _cancella_zone(self):
        n = self.var_st_file.get()
        if n and self.stato:
            rp = dict(self.scelte(n).get("riparazione") or {})
            rp["rettangoli"] = []
            self.stato.imposta(n, "riparazione", rp, "zone disegnate cancellate")
            self.pipe.dimentica_file(n)
            self._mostra_stato_file_st()

    def _torna_fondo(self):
        n = self.var_st_file.get()
        self.sez_bg.apri()
        if n:
            self.var_bg_file.set(n)
            self._anteprime("bg")

    # ================================================================ sezione 5: stack
    def _apri_sez_stack(self):
        for w in self.fr_pesi.winfo_children():
            w.destroy()
        self.num_pesi = {}
        pesi = (self.stato.globale("pesi_stack", {}) if self.stato else {})
        for i, n in enumerate(self.inclusi()):
            ttk.Label(self.fr_pesi, text=n).grid(row=i, column=0, sticky="w")
            num = NumeroGrande(self.fr_pesi, pesi.get(n, 1.0), passo=0.1, minimo=0, massimo=100)
            num.grid(row=i, column=1, sticky="w", padx=6, pady=1)
            self.num_pesi[n] = num
        inc = [n for n in self.inclusi() if self.immagine(n).wcs is not None]
        if inc and not getattr(self, "_rif_utente", False):
            self.var_rif.set(min(inc, key=lambda n: M.scala_pixel_deg(self.immagine(n).wcs)))
        self._stima_griglia()
        self._aggiorna_lbl_ritaglio()
        if self.stack is not None:
            self._mostra_stack()
        else:
            self.visore.messaggio("Stack in bianco e nero", "Premi 'Crea lo stack'.")

    def _parametri_griglia(self):
        return GRIGLIE.get(self.var_griglia.get(), "riferimento"), self.var_rif.get(), self.num_rid.get()

    def _calcola_griglia(self, nomi):
        tipo, rif, k = self._parametri_griglia()
        imgs = [self.immagine(n) for n in nomi]
        if tipo == "mosaico":
            return M.griglia_mosaico(imgs, k)
        r = self.immagine(rif if rif in nomi else nomi[0])
        return M.griglia_riferimento(r, k)

    def _stima_griglia(self):
        self._salva_cfg("griglia", self.var_griglia.get())
        self._salva_cfg("riduzione", self.num_rid.get())
        if not self.stato:
            return
        nomi = [n for n in self.inclusi() if self.immagine(n).wcs is not None]
        if not nomi:
            self.lbl_griglia.configure(text="Nessun file incluso con coordinate WCS.")
            return
        try:
            w, (ny, nx) = self._calcola_griglia(nomi)
            mb = ny * nx * 4 / 1024**2
            self.lbl_griglia.configure(text=f"Immagine finale: {nx} × {ny} pixel, {M.scala_pixel_deg(w) * 3600:.3f}″/pixel "
                                            f"(~{mb:.0f} MB per filtro, {len(nomi)} filtri)")
        except Exception as e:
            self.lbl_griglia.configure(text=f"Griglia non calcolabile: {e}")

    def _file_pronti(self):
        nomi = self.inclusi()
        mancano = [n for n in nomi if not (self.scelte(n).get("nr") and self.scelte(n).get("stretch"))]
        senza_wcs = [n for n in nomi if self.immagine(n).wcs is None]
        if mancano:
            messagebox.showinfo("Scelte mancanti", "Questi file non hanno ancora riduzione del rumore e stretching:\n"
                                + "\n".join(mancano) + "\n\nCompletali o scartali (sezioni 2-4).")
            return None
        if senza_wcs:
            messagebox.showinfo("Coordinate mancanti", "Questi file non hanno coordinate WCS e non possono essere "
                                                       "sovrapposti:\n" + "\n".join(senza_wcs) + "\n\nScartali nella sezione 4.")
            return None
        if not nomi:
            messagebox.showinfo("Nessun file", "Non ci sono file inclusi.")
            return None
        return nomi

    def _prepara_riproiezioni(self, nomi, log, prog):
        """Elabora ogni file a piena risoluzione e lo riproietta sulla griglia comune."""
        w, shape = self._calcola_griglia(nomi)
        chiave_g = (self._parametri_griglia(), tuple(nomi))
        if self.griglia_stack is None or self.griglia_stack[0] != chiave_g:
            for _ch, a in self.riproiettate.values():
                self.pipe.mem.libera(a)
            self.riproiettate = {}
        self.griglia_stack = (chiave_g, w, shape)
        scala = M.scala_pixel_deg(w)
        for i, n in enumerate(nomi):
            s = self.scelte(n)
            ch = M.chiave_scelte(self.immagine(n), s)
            if n in self.riproiettate and self.riproiettate[n][0] == ch:
                continue
            log(f"[{i + 1}/{len(nomi)}] {n}: elaborazione a piena risoluzione")
            img = self.immagine(n)
            arr, limiti = self.pipe.elabora(img, s, self.rumore(n),
                                            lambda fase, k, t: prog(f"{n} {fase}", k, t), log)
            arr_r, w_r = M.riduci_per_griglia(arr, img.wcs, scala, self.pipe.mem)
            log(f"{n}: riproiezione sulla griglia comune {shape[1]}×{shape[0]}")
            dst = self.pipe.mem.alloca(f"rip_{Path(n).stem}", shape)
            M.riproietta(arr_r, w_r, w, shape, dst, progresso=lambda k, t: prog(f"{n} riproiezione", k, t))
            self.riproiettate[n] = (ch, dst)
            del arr, arr_r
            M.restituisci_ram()
        return w, shape

    def _crea_stack(self):
        if not self._richiedi_madre():
            return
        nomi = self._file_pronti()
        if not nomi:
            return
        pesi = [self.num_pesi[n].get() if n in self.num_pesi else 1.0 for n in nomi]
        metodo = METODI_STACK.get(self.var_metodo.get(), "media")
        self._salva_cfg("metodo", self.var_metodo.get())
        self.stato.imposta_globale("pesi_stack", dict(zip(nomi, pesi)))

        def lavoro(log, prog):
            w, shape = self._prepara_riproiezioni(nomi, log, prog)
            log("Combinazione dei filtri")
            if self.stack is not None:
                self.pipe.mem.libera(self.stack[0])
            dst = self.pipe.mem.alloca("stack_bn", shape)
            M.combina_bn([self.riproiettate[n][1] for n in nomi], pesi, metodo, dst,
                         lambda k, t: prog("Combinazione", k, t))
            return dst, w, M.anteprima(dst)

        def fatto(r):
            dst, w, ant = r
            self.stack = (dst, w)
            self._stack_ant = ant
            self.stato.imposta_globale("stack", {"file": nomi, "pesi": pesi, "metodo": metodo,
                                                 "griglia": self.var_griglia.get(), "riferimento": self.var_rif.get(),
                                                 "riduzione": self.num_rid.get(),
                                                 "dimensione": f"{dst.shape[1]}x{dst.shape[0]}"},
                                       f"Stack B/N creato con {len(nomi)} filtri ({metodo})")
            self.sez_stack.imposta_stato("creato (da salvare)")
            self._mostra_stack()

        self.esegui("Creazione dello stack in bianco e nero", lavoro, fatto)

    # ---------------------------------------------------------------- ritaglio
    def ritaglio(self):
        """[y0, y1, x0, x1] del ritaglio valido per la griglia attuale, oppure None."""
        if not self.stato or self.griglia_stack is None:
            return None
        r = self.stato.globale("ritaglio")
        if r and list(r.get("forma", [])) == list(self.griglia_stack[2]):
            return r["box"]
        return None

    def _imposta_ritaglio(self, box, testo):
        ny, nx = self.griglia_stack[2]
        y0, y1, x0, x1 = box
        y0, x0 = max(0, int(y0)), max(0, int(x0))
        y1, x1 = min(ny, int(y1)), min(nx, int(x1))
        if y1 - y0 < 16 or x1 - x0 < 16:
            messagebox.showinfo("Ritaglio", "Il rettangolo è troppo piccolo.")
            return
        self.stato.imposta_globale("ritaglio", {"box": [y0, y1, x0, x1], "forma": [ny, nx]},
                                   f"Ritaglio {testo}: x {x0}-{x1}, y {y0}-{y1} ({x1 - x0}×{y1 - y0} pixel)")
        self._aggiorna_lbl_ritaglio()
        self.sez_stack.imposta_stato("ritagliato (da salvare)")
        self._mostra_stack()
        if self.rgb_anteprime:
            self._aggiorna_rgb(ricrea=True) if self.sezione_attiva is self.sez_rgb else None

    def _aggiorna_lbl_ritaglio(self):
        r = self.ritaglio()
        if r:
            y0, y1, x0, x1 = r
            self.lbl_ritaglio.configure(text=f"Ritaglio attivo: {x1 - x0} × {y1 - y0} pixel "
                                             f"(x {x0}–{x1}, y {y0}–{y1}).")
        else:
            self.lbl_ritaglio.configure(text="Nessun ritaglio: si salva l'immagine intera.")

    def _disegna_ritaglio(self):
        if self.stack is None:
            messagebox.showinfo("Ritaglio", "Crea prima lo stack.")
            return
        ny, nx = self.griglia_stack[2]
        zone = []
        r = self.ritaglio()
        if r:
            zone.append((r[0], r[1], r[2], r[3], "yellow"))
        quadrato = " (verrà reso quadrato)" if self.chk_quadrato.get() else ""
        self.visore.rettangolo(self._stack_ant, f"Trascina con il mouse il rettangolo da tenere{quadrato}. "
                                                f"Giallo = ritaglio attuale.", (0, nx, 0, ny), zone,
                               self._fine_disegno_ritaglio)
        self.visore.info("Dopo aver disegnato, l'anteprima mostra il risultato. Premi 'Salva' per salvarlo.")

    def _fine_disegno_ritaglio(self, rett):
        y0, y1, x0, x1 = rett
        if self.chk_quadrato.get():
            lato = min(y1 - y0, x1 - x0)
            cy, cx = (y0 + y1) // 2, (x0 + x1) // 2
            y0, y1, x0, x1 = cy - lato // 2, cy - lato // 2 + lato, cx - lato // 2, cx - lato // 2 + lato
        self._imposta_ritaglio([y0, y1, x0, x1], "disegnato" + (" (quadrato)" if self.chk_quadrato.get() else ""))

    def _ritaglio_auto(self):
        if self.stack is None:
            messagebox.showinfo("Ritaglio", "Crea prima lo stack.")
            return
        nomi = self.stato.globale("stack", {}).get("file", [])
        arrays = [self.riproiettate[n][1] for n in nomi if n in self.riproiettate]
        if not arrays:
            arrays = [self.stack[0]]
        quadrato = self.chk_quadrato.get()

        def lavoro(log, prog):
            box = M.ritaglio_automatico(arrays)
            if box and quadrato:
                y0, y1, x0, x1 = box
                lato = min(y1 - y0, x1 - x0)
                cy, cx = (y0 + y1) // 2, (x0 + x1) // 2
                box = [cy - lato // 2, cy - lato // 2 + lato, cx - lato // 2, cx - lato // 2 + lato]
            return box

        def fatto(box):
            if not box:
                messagebox.showinfo("Ritaglio", "Non ho trovato una zona in cui tutti i filtri hanno dati.")
                return
            self._imposta_ritaglio(box, "automatico" + (" (quadrato)" if quadrato else ""))

        self.esegui("Ricerca della zona con i dati di tutti i filtri", lavoro, fatto)

    def _annulla_ritaglio(self):
        if self.stato:
            self.stato.imposta_globale("ritaglio", None, "Ritaglio annullato")
            self._aggiorna_lbl_ritaglio()
            if self.stack is not None:
                self._mostra_stack()

    def _mostra_stack(self):
        dst, w = self.stack
        ant = self._stack_ant
        ny, nx = dst.shape
        r = self.ritaglio()
        if r:
            y0, y1, x0, x1 = r
            f = max(1, math.ceil(max(ny, nx) / max(ant.shape)))
            self.visore.immagine(ant[y0 // f:math.ceil(y1 / f), x0 // f:math.ceil(x1 / f)],
                                 f"Stack in bianco e nero, RITAGLIATO: {x1 - x0} × {y1 - y0} pixel "
                                 f"(dell'intero {nx} × {ny}). Premi 'Salva' per scriverlo su disco.")
        else:
            self.visore.immagine(ant, f"Stack in bianco e nero — {nx} × {ny} pixel "
                                      f"(anteprima ridotta). Premi 'Salva' per scriverlo su disco.")

    def _salva_stack(self):
        if self.stack is None:
            messagebox.showinfo("Stack", "Crea prima lo stack.")
            return
        dst, w = self.stack
        fmt = self.var_fmt.get()
        self._salva_cfg("formato", fmt)
        p_fits = self.madre / "stack_bn.fits"
        p_img = self.madre / f"stack_bn{M.FORMATI_IMMAGINE[fmt]}"
        if (p_fits.exists() or p_img.exists()) and not messagebox.askyesno("Sovrascrivere?",
                                                                            f"{p_fits.name} esiste già. Sovrascriverlo?"):
            return
        st = self.stato.globale("stack", {})
        storia = [f"Stack B/N di {', '.join(st.get('file', []))}", f"metodo {st.get('metodo')}"]
        for n in st.get("file", []):
            s = self.scelte(n)
            storia.append(f"{n}: NR {self._testo_scelta('nr', s.get('nr'))}; fondo "
                          f"{self._testo_scelta('bg', s.get('bg')) if s.get('usa_bg', True) else 'no'}; "
                          f"{self._testo_scelta('st', s.get('stretch'))}")

        r = self.ritaglio()
        dati, w_out = dst, w
        if r:
            y0, y1, x0, x1 = r
            dati, w_out = dst[y0:y1, x0:x1], M.ritaglia_wcs(w, y0, x0)
            storia.append(f"Ritaglio: x {x0}-{x1}, y {y0}-{y1} dello stack {dst.shape[1]}x{dst.shape[0]}")

        def lavoro(log, prog):
            M.salva_fits(dati, w_out, p_fits, storia)
            M.salva_immagine(dati, p_img, fmt, firma=self.chk_firma.get())
            return p_fits, p_img

        def fatto(r):
            self.stato.evento(f"Stack B/N salvato: {r[0].name}, {r[1].name}")
            self.stato.salva()
            self.sez_stack.imposta_stato("salvato")
            self.scrivi(f"Salvati {r[0]} e {r[1].name}")
            messagebox.showinfo("Salvato", f"Salvati nella cartella madre:\n{r[0].name}\n{r[1].name}")

        self.esegui("Salvataggio dello stack", lavoro, fatto)

    # ================================================================ sezione 6: RGB
    def _apri_sez_rgb(self):
        if self.rgb_anteprime:
            if self.pesi_rgb:                       # conserva i pesi già regolati
                self.pesi_iniziali = self._pesi_correnti()
            self._costruisci_barre_rgb()
            self._aggiorna_rgb(ricrea=True)
        else:
            self.visore.messaggio("Stack RGB", "Premi 'Prepara le immagini' (usa la griglia e le scelte dello stack B/N).")

    def _prepara_rgb(self):
        if not self._richiedi_madre():
            return
        nomi = self._file_pronti()
        if not nomi:
            return

        def lavoro(log, prog):
            self._prepara_riproiezioni(nomi, log, prog)
            return {n: M.anteprima(self.riproiettate[n][1], 1400) for n in nomi}

        def fatto(ant):
            self.rgb_anteprime = ant
            salvati = self.stato.globale("rgb", {}).get("pesi", {})
            if all(n in salvati for n in ant):
                self.pesi_iniziali = {n: salvati[n] for n in ant}
            else:
                self.pesi_iniziali = self._pesi_auto()
            self._costruisci_barre_rgb()
            self._aggiorna_rgb(ricrea=True)

        self.esegui("Preparazione delle immagini per l'RGB", lavoro, fatto)

    def _pesi_auto(self):
        nomi = list(self.rgb_anteprime)
        lam = [self.immagine(n).lambda_um for n in nomi]
        return dict(zip(nomi, M.pesi_cromatici(lam)))

    def _costruisci_barre_rgb(self):
        for w in self.fr_rgb.winfo_children():
            w.destroy()
        self._scale_rgb = []
        self.pesi_rgb = {}
        p = lambda v: max(1, round(v * self.zoom))
        for i, n in enumerate(self.rgb_anteprime):
            img = self.immagine(n)
            f = ttk.LabelFrame(self.fr_rgb, text=f"{n}  ·  {img.filtro}  ·  {M.descrivi_lunghezza(img.lambda_um)}", padding=4)
            f.pack(fill=tk.X, pady=3)
            vari = []
            for c, (nome_c, colore) in enumerate((("R", "#e04040"), ("G", "#30a030"), ("B", "#3060e0"))):
                v = tk.DoubleVar(value=float(self.pesi_iniziali.get(n, (1, 1, 1))[c]))
                riga = ttk.Frame(f)
                riga.pack(fill=tk.X)
                tk.Label(riga, text=nome_c, fg="white", bg=colore, width=3, font=self.f_grassetto).pack(side=tk.LEFT)
                sc = tk.Scale(riga, from_=0, to=1, resolution=0.01, orient=tk.HORIZONTAL, variable=v,
                              troughcolor=colore, width=p(18), sliderlength=p(28), showvalue=True,
                              font="TkDefaultFont", highlightthickness=0, command=lambda _v: self._aggiorna_rgb())
                sc.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=4)
                self._scale_rgb.append(sc)
                vari.append(v)
            self.pesi_rgb[n] = vari

    def _rgb_auto(self):
        if not self.rgb_anteprime:
            self._prepara_rgb()
            return
        auto = self._pesi_auto()
        for n, vari in self.pesi_rgb.items():
            for c in range(3):
                vari[c].set(auto[n][c])
        self._aggiorna_rgb()

    def _pesi_correnti(self):
        return {n: tuple(round(v.get(), 3) for v in vari) for n, vari in self.pesi_rgb.items()}

    def _aggiorna_rgb(self, ricrea=False):
        if not self.rgb_anteprime:
            return
        if getattr(self, "_rgb_timer", None):
            self.after_cancel(self._rgb_timer)

        def calcola():
            self._rgb_timer = None
            nomi = list(self.rgb_anteprime)
            pesi = self._pesi_correnti()
            arr = [self.rgb_anteprime[n] for n in nomi]
            r = self.ritaglio()
            if r and self.griglia_stack is not None:
                ny, nx = self.griglia_stack[2]
                f = max(1, math.ceil(max(ny, nx) / max(arr[0].shape)))
                y0, y1, x0, x1 = r
                arr = [a[y0 // f:math.ceil(y1 / f), x0 // f:math.ceil(x1 / f)] for a in arr]
            lista = [pesi[n] for n in nomi]
            img = M.combina_rgb(arr, lista)
            self._scala_rgb = M.bilanciamento(img) if self.chk_bil.get() else (1.0, 1.0, 1.0)
            img = np.clip(img * np.array(self._scala_rgb, dtype=np.float32), 0, 1)
            if ricrea or not self.visore.aggiorna_rgb(img):
                rit = " (con il ritaglio dello stack B/N)" if self.ritaglio() else ""
                self.visore.immagine(img, f"Stack RGB{rit} — muovi le barre a sinistra: l'anteprima si aggiorna subito",
                                     rgb=True)

        self._rgb_timer = self.after(60, calcola)

    def _salva_rgb(self):
        if not self.rgb_anteprime or self.griglia_stack is None:
            messagebox.showinfo("RGB", "Premi prima 'Prepara le immagini'.")
            return
        fmt = self.var_fmt_rgb.get()
        self._salva_cfg("formato_rgb", fmt)
        nomi = list(self.rgb_anteprime)
        pesi = self._pesi_correnti()
        scala = getattr(self, "_scala_rgb", (1.0, 1.0, 1.0))
        _, w, shape = self.griglia_stack
        r = self.ritaglio()
        ry0, ry1, rx0, rx1 = r if r else (0, shape[0], 0, shape[1])
        if r:
            w = M.ritaglia_wcs(w, ry0, rx0)
            shape = (ry1 - ry0, rx1 - rx0)
        p_fits = self.madre / "stack_rgb.fits"
        p_img = self.madre / f"stack_rgb{M.FORMATI_IMMAGINE[fmt]}"
        if p_fits.exists() and not messagebox.askyesno("Sovrascrivere?", f"{p_fits.name} esiste già. Sovrascriverlo?"):
            return
        self.stato.imposta_globale("rgb", {"pesi": pesi, "scala_canali": [round(x, 4) for x in scala]},
                                   "Pesi RGB: " + "; ".join(f"{n} R{p[0]} G{p[1]} B{p[2]}" for n, p in pesi.items()))

        def lavoro(log, prog):
            ny, nx = shape
            rgb = self.pipe.mem.alloca("rgb", (ny, nx, 3))
            arr = [self.riproiettate[n][1] for n in nomi]
            lista = [pesi[n] for n in nomi]
            for y0 in range(0, ny, 512):
                y1 = min(y0 + 512, ny)
                rgb[y0:y1] = M.combina_rgb([a[ry0 + y0:ry0 + y1, rx0:rx1] for a in arr], lista, scala)
                prog("Composizione RGB", y1, ny)
            storia = ["Stack RGB"] + [f"{n}: R {p[0]} G {p[1]} B {p[2]}" for n, p in pesi.items()]
            if r:
                storia.append(f"Ritaglio: x {rx0}-{rx1}, y {ry0}-{ry1}")
            log("Scrittura del FITS RGB (cubo a 3 piani)")
            M.salva_fits([rgb[..., 0], rgb[..., 1], rgb[..., 2]], w, p_fits, storia)
            log("Scrittura dell'immagine")
            M.salva_immagine(rgb, p_img, fmt, firma=self.chk_firma.get())
            self.pipe.mem.libera(rgb)
            del rgb
            M.restituisci_ram()
            return p_fits, p_img

        def fatto(r):
            self.stato.evento(f"RGB salvato: {r[0].name}, {r[1].name}")
            self.stato.salva()
            self.sez_rgb.imposta_stato("salvato")
            messagebox.showinfo("Salvato", f"Salvati nella cartella madre:\n{r[0].name}\n{r[1].name}")

        self.esegui("Salvataggio dell'immagine RGB", lavoro, fatto)

    # ================================================================ sezione 7: analisi
    def _apri_sez_analisi(self):
        if self.risultati_analisi:
            self._mostra_analisi(self.indice_analisi)
        else:
            self.visore.messaggio("Analisi scientifiche", "Scegli le analisi e premi 'Esegui le analisi'.")

    def _analisi(self):
        if not self.madre:
            messagebox.showinfo("Analisi", "Scegli prima la cartella madre.")
            return
        p = self.madre / "stack_bn.fits"
        if not p.exists():
            messagebox.showinfo("Analisi", "Manca stack_bn.fits nella cartella madre: crea e salva lo stack (sezione 5).")
            return
        scelte = [k for k, ch in self.chk_analisi.items() if ch.get()]
        self._salva_cfg("analisi", {k: ch.get() for k, ch in self.chk_analisi.items()})
        fmt = self.var_fmt_an.get()
        self._salva_cfg("formato_an", fmt)
        max_nomi = self.num_nomi.get()
        self._salva_cfg("max_nomi", max_nomi)

        def lavoro(log, prog):
            d, w, f = M.carica_stack(p)
            if f > 1:
                log(f"Stack ridotto di {f} volte per le analisi ({d.shape[1]}×{d.shape[0]})")
            risultati = []
            for i, k in enumerate(scelte):
                prog("Analisi", i, len(scelte))
                log(f"Analisi: {M.ANALISI[k]}")
                try:
                    if k == "oggetti":
                        fig, righe = M.analisi_oggetti(d, w, max_nomi=max_nomi, log=log)
                    else:
                        fig, righe = M.FUNZIONI_ANALISI[k](d, w)
                    pf = M.salva_figura(fig, self.madre / f"analisi_{k}", fmt, firma=self.chk_firma.get())
                    pc = M.scrivi_csv(self.madre / M.CSV_ANALISI[k], righe) if k in M.CSV_ANALISI and righe else None
                    risultati.append({"chiave": k, "titolo": M.ANALISI[k], "rgba": M.figura_in_rgba(fig),
                                      "file": pf, "csv": pc, "righe": len(righe) if righe else 0})
                    log(f"  salvato {pf.name}" + (f" e {pc.name} ({len(righe)} righe)" if pc else ""))
                except Exception as e:
                    log(f"  ERRORE in '{M.ANALISI[k]}': {type(e).__name__}: {e}")
                    risultati.append({"chiave": k, "titolo": M.ANALISI[k] + " — NON riuscita", "rgba": None,
                                      "errore": f"{type(e).__name__}: {e}"})
            return risultati

        def fatto(ris):
            self.risultati_analisi = ris
            self.stato.evento("Analisi eseguite: " + ", ".join(r["chiave"] for r in ris))
            self.stato.salva()
            ok = sum(1 for r in ris if r.get("rgba") is not None)
            self.sez_an.imposta_stato(f"{ok}/{len(ris)} completate")
            self._mostra_analisi(0)

        self.esegui("Analisi scientifiche", lavoro, fatto)

    def _mostra_analisi(self, i):
        if not self.risultati_analisi:
            return
        self.indice_analisi = i % len(self.risultati_analisi)
        r = self.risultati_analisi[self.indice_analisi]
        titolo = f"Analisi {self.indice_analisi + 1}/{len(self.risultati_analisi)}: {r['titolo']}"
        if r.get("rgba") is None:
            self.visore.messaggio(titolo, r.get("errore", ""))
        else:
            self.visore.figura_rgba(r["rgba"], titolo)
        info = f"File: {r['file'].name}" if r.get("file") else ""
        if r.get("csv"):
            info += f"   |   {r['csv'].name} ({r['righe']} righe)"
        self.visore.info(info)
        self.visore.imposta_navigazione(lambda d: self._mostra_analisi(self.indice_analisi + d),
                                        ("◀  Analisi precedente", "Analisi successiva  ▶"))

    # ================================================================ chiusura
    def _chiudi(self):
        if self.occupato and not messagebox.askyesno("Uscire?", "C'è un calcolo in corso. Uscire comunque?"):
            return
        for img in self._immagini.values():
            img.chiudi()
        if self.pipe:
            self.pipe.mem.pulisci_temporanei()
        self.destroy()


# =============================================================================
# AUTOTEST (dati sintetici, senza rete) E DIAGNOSTICA
# =============================================================================
def _fits_sintetico(percorso, filtro, rot_deg=0.0, scala_as=0.06, n=700, buco=True, seme=0):
    from astropy.io import fits
    from astropy.wcs import WCS
    rng = np.random.default_rng(seme)
    yy, xx = np.mgrid[:n, :n]
    d = (rng.normal(0, 1, (n, n)) + 10 + 0.003 * xx
         + 60 * np.exp(-((yy - n / 2) ** 2 + (xx - n / 2) ** 2) / (2 * (n / 8) ** 2)))
    for k in range(40):
        y, x = rng.integers(20, n - 20, 2)
        d += 200 * np.exp(-((yy - y) ** 2 + (xx - x) ** 2) / (2 * 1.5 ** 2))
    d[:15] = 0
    d[:, :10] = 0
    if buco:
        d[300:330, 400:440] = 0                 # zona nera interna
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [14.77, -72.18]
    w.wcs.crpix = [n / 2, n / 2]
    s = scala_as / 3600
    c, si = math.cos(math.radians(rot_deg)), math.sin(math.radians(rot_deg))
    w.wcs.cd = [[-s * c, s * si], [s * si, s * c]]
    hp = fits.Header()
    hp["TELESCOP"], hp["INSTRUME"], hp["FILTER"] = "JWST", "NIRCAM", filtro
    sci = fits.ImageHDU(d.astype(np.float32), header=w.to_header(), name="SCI")
    Path(percorso).parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList([fits.PrimaryHDU(header=hp), sci]).writeto(percorso, overwrite=True)


def autotest():
    import tempfile
    import shutil
    print(f"AUTOTEST {NOME} {M.VERSIONE} - {M.COPYRIGHT}")
    ok = True

    def v(cond, testo):
        nonlocal ok
        print(("[ok] " if cond else "[ERRORE] ") + testo, flush=True)
        ok &= bool(cond)

    madre = Path(tempfile.mkdtemp(prefix="autotest_pipeline_"))
    try:
        _fits_sintetico(madre / "mastDownload" / "jw01" / "jw01_nircam_f090w_i2d.fits", "F090W", 0, 0.03, seme=1)
        _fits_sintetico(madre / "mastDownload" / "jw02" / "jw02_nircam_f200w_i2d.fits", "F200W", 20, 0.06, seme=2)
        _fits_sintetico(madre / "altro" / "copia_f200w.fits", "F200W", 35, 0.06, buco=False, seme=3)
        trovati = M.trova_fits_da_rinominare(madre)
        nomi = sorted(t["nuovo_nome"] for t in trovati)
        v(nomi == ["f090.fits", "f200a.fits", "f200b.fits"], f"nomi dai filtri: {nomi}")
        copiati = M.copia_e_verifica(trovati, madre, log=lambda m: None)
        M.cancella_originali(copiati, madre, log=lambda m: None)
        v(not (madre / "mastDownload").exists() and len(M.elenca_file_pipeline(madre)) == 3,
          "copia verificata, originali e cartelle vuote cancellati")
        img = M.ImmagineFITS(madre / "f090.fits")
        v(abs(img.lambda_um - 0.90) < 1e-6 and img.wcs is not None, f"filtro {img.filtro}, λ {img.lambda_um} µm, WCS")
        v(M.lunghezza_onda_um("F656N", {"TELESCOP": "HST"}) == 0.656 and
          M.lunghezza_onda_um("F1130W", {"TELESCOP": "JWST"}) == 11.3, "lunghezze d'onda HST/JWST")
        med, sig = M.stima_rumore(img.dati)
        v(0.5 < sig < 2, f"stima del rumore σ = {sig:.2f}")
        crop, (y0, x0) = img.ritaglio(None, None, 400)
        for a in M.ALGORITMI_NR:
            t0 = time.time()
            r = M.riduci_rumore(crop, a, M.ALGORITMI_NR[a]["valori"][1], med, sig)
            masc = M.maschera_segnale(crop)
            s0, s1 = M.snr_medio(crop, "segnale", masc), M.snr_medio(r, "segnale", masc)
            v(r.shape == crop.shape and np.isfinite(r).mean() > 0.8 and s1 > s0,
              f"riduzione rumore '{a}' ({time.time() - t0:.1f} s), SNR {s0:.1f} -> {s1:.1f}")
        mem = M.Memoria("disco", 4)
        pipe = M.Pipeline(madre, mem)
        for a in M.ALGORITMI_FONDO:
            f = pipe.fondo_ritaglio(img, a, 128, y0, x0, crop.shape)
            v(f.shape == crop.shape and np.isfinite(f).all() and abs(np.median(f) - 10.6) < 3,
              f"fondo '{a}' mediana {np.median(f):.2f}")
        zone = M.trova_zone_rovinate(img)
        v(len(zone) >= 1, f"zone rovinate trovate: {len(zone)}")
        c2, (a0, b0) = img.ritaglio(None, None, 600)
        prima = np.isnan(c2[300 - a0:330 - a0, 400 - b0:440 - b0]).mean()
        M.ripara(c2, a0, b0, zone, [[100, 120, 100, 130]])
        v(prima > 0.9 and np.isfinite(c2[300 - a0:330 - a0, 400 - b0:440 - b0]).all(), "ricostruzione automatica + manuale")
        for a, info in M.ALGORITMI_STRETCH.items():
            r = M.applica_stretch(crop, a, info["valori"][0])
            v(np.nanmin(r) >= 0 and np.nanmax(r) <= 1, f"stretching '{a}'")
        scelte = {"nr": {"algo": "gaussiano", "valore": 1.0}, "bg": {"algo": "mediana", "box": 128},
                  "usa_bg": True, "stretch": {"algo": "asinh", "valore": 0.1},
                  "riparazione": {"auto": True, "zone_auto": zone, "rettangoli": []}}
        risultati = {}
        for modo in ("ram", "disco"):
            pipe.mem.modo = modo
            for p in pipe.cache.glob("elab_*"):
                p.unlink()
            arr, lim = pipe.elabora(img, scelte, (med, sig), log=lambda m: None)
            risultati[modo] = np.array(arr)
        v(np.allclose(risultati["ram"], risultati["disco"], equal_nan=True), "elaborazione completa RAM = disco")
        v(np.isfinite(risultati["ram"][310, 420]), "zona rovinata ricostruita nell'elaborazione completa")
        imgs = [M.ImmagineFITS(madre / n) for n in ("f090.fits", "f200a.fits", "f200b.fits")]
        st = []
        for im in imgs:
            a, _ = pipe.elabora(im, scelte, M.stima_rumore(im.dati), log=lambda m: None)
            st.append(a)
        # tre file con la stessa rotazione: il mosaico orientato deve essere quasi tutto pieno
        ruotati = []
        for k_, f_ in enumerate(("r1", "r2")):
            _fits_sintetico(madre / "ruotati" / f"{f_}.fits", "F200W", 30, 0.06, buco=False, seme=10 + k_)
            ruotati.append(M.ImmagineFITS(madre / "ruotati" / f"{f_}.fits"))
        wo, so = M.griglia_mosaico(ruotati, 1)
        wn, sn_ = M.griglia_mosaico(ruotati, 1, nord_in_alto=True)
        v(so[0] * so[1] < 0.8 * sn_[0] * sn_[1] and abs(so[0] - 700) < 20 and abs(so[1] - 700) < 20,
          f"mosaico orientato {so[1]}×{so[0]} (Nord in alto {sn_[1]}×{sn_[0]})")
        d_ = mem.alloca("ro", so)
        M.riproietta(np.array(ruotati[0].blocco(0, 700, 0, 700)), ruotati[0].wcs, wo, so, d_)
        v(np.isfinite(d_).mean() > 0.9, f"mosaico orientato senza angoli vuoti (pieno al {100 * np.isfinite(d_).mean():.0f}%)")
        for im in ruotati:
            im.chiudi()
        for tipo in ("riferimento", "mosaico"):
            w, shape = M.griglia_riferimento(imgs[0], 2) if tipo == "riferimento" else M.griglia_mosaico(imgs, 2)
            rip = []
            for im, a in zip(imgs, st):
                ar, wr = M.riduci_per_griglia(a, im.wcs, M.scala_pixel_deg(w), mem)
                d = mem.alloca("r", shape)
                M.riproietta(ar, wr, w, shape, d, righe=97)
                rip.append(d)
            stack = mem.alloca("s", shape)
            M.combina_bn(rip, [1, 1, 1], "media", stack)
            cen = w.wcs_world2pix([[14.77, -72.18]], 0)[0]
            yc, xc = int(cen[1]), int(cen[0])
            v(np.isfinite(stack[yc, xc]) and stack[yc, xc] > 0.5,
              f"stack su griglia '{tipo}' {shape[1]}×{shape[0]}, centro allineato ({stack[yc, xc]:.2f})")
        box = M.ritaglio_automatico(rip)
        y0r, y1r, x0r, x1r = box
        v(all(np.isfinite(np.asarray(a)[y0r:y1r, x0r:x1r]).all() for a in rip),
          f"ritaglio automatico {x1r - x0r}×{y1r - y0r} su {shape[1]}×{shape[0]}: dati di tutti i filtri ovunque")
        wr = M.ritaglia_wcs(w, y0r, x0r)
        a_ = w.wcs_pix2world([[x0r + 10, y0r + 20]], 0)[0]
        b_ = wr.wcs_pix2world([[10, 20]], 0)[0]
        v(np.allclose(a_, b_, atol=1e-9), "coordinate WCS del ritaglio corrette")
        for metodo in ("mediana", "massimo"):
            M.combina_bn(rip, [1, 1, 1], metodo, stack)
            v(np.nanmax(stack) <= 1, f"combinazione '{metodo}'")
        p = M.salva_fits(stack, w, madre / "stack_bn.fits", ["autotest: λ 0,9 µm, 350×350, è già → ok"])
        d, w2, f = M.carica_stack(p)
        v(d.shape == stack.shape and w2.has_celestial, "FITS dello stack con WCS")
        hp = M.fits.getheader(p)
        v(hp.get("AUTHOR") == M.AUTORE and "(C) 2026" in str(hp.get("COMMENT")),
          f"autore e copyright nell'intestazione del FITS: {M.AUTORE}")
        for fmt in M.FORMATI_IMMAGINE:
            q = M.salva_immagine(stack, madre / f"stack_bn{M.FORMATI_IMMAGINE[fmt]}", fmt, firma=True)
            v(q.stat().st_size > 1000, f"immagine {fmt} con firma")
        prova = np.zeros((400, 600), np.float32)
        q = M.salva_immagine(prova, madre / "firma.png", "PNG 8 bit", firma=True)
        from PIL import Image
        fa = np.asarray(Image.open(q))
        v(fa[-60:, -300:].max() > 100 and fa[:200, :300].max() == 0, "firma visibile solo nell'angolo in basso a destra")
        pesi = M.pesi_cromatici([im.lambda_um for im in imgs])
        v(pesi[0][2] == 1.0 and max(pesi[1][0], pesi[2][0]) > 0.5, f"pesi cromatici {pesi}")
        rgb = M.combina_rgb(rip, pesi)
        M.salva_fits([rgb[..., 0], rgb[..., 1], rgb[..., 2]], w, madre / "stack_rgb.fits", ["rgb"])
        M.salva_immagine(rgb, madre / "stack_rgb.tif", "TIFF 16 bit")
        v((madre / "stack_rgb.tif").stat().st_size > 1000, "RGB salvato (FITS a 3 piani + TIFF 16 bit)")
        for k in M.ANALISI:
            if k == "oggetti":
                fig, righe = M.analisi_oggetti(d, w2, oggetti=[dict(nome="NGC 346", nome_simbad="NGC 346",
                    altri_nomi="", ra=14.77, dec=-72.18, x=d.shape[1] / 2, y=d.shape[0] / 2, tipo="Cl*",
                    riferimenti=900, catalogo=True)])
            else:
                fig, righe = M.FUNZIONI_ANALISI[k](d, w2)
            q = M.salva_figura(fig, madre / f"analisi_{k}", "PNG 8 bit")
            if righe:
                M.scrivi_csv(madre / M.CSV_ANALISI.get(k, f"analisi_{k}.csv"), righe)
            v(q.stat().st_size > 5000, f"analisi '{k}'" + (f" ({len(righe)} righe CSV)" if righe else ""))
        ogg = M.oggetti_openngc(w2, d.shape)
        v(any(o["nome"] == "NGC 346" for o in ogg), f"catalogo OpenNGC incluso: {len(ogg)} oggetti nel campo, "
                                                    f"tra cui {', '.join(o['nome'] for o in ogg[:4])}")
        fig, righe = M.analisi_oggetti(d, w2, oggetti=ogg)
        v(righe and righe[0]["fonte"] == "OpenNGC", "analisi oggetti con OpenNGC (senza Internet)")
        stato = M.Stato(madre)
        stato.imposta("f090.fits", "nr", scelte["nr"], "prova")
        v((madre / "pipeline_log.txt").read_text(encoding="utf-8").count("f090.fits") >= 1, "stato e log")
        mem.pulisci_temporanei()
        for im in imgs + [img]:
            im.chiudi()
        r = tk.Tk()
        r.withdraw()
        ttk.Style(r).theme_use("clam")
        r.destroy()
        v(True, "finestra Tk")
    except Exception:
        traceback.print_exc()
        ok = False
    finally:
        import gc
        gc.collect()
        shutil.rmtree(madre, ignore_errors=True)
    print("AUTOTEST SUPERATO" if ok else "AUTOTEST FALLITO")
    return 0 if ok else 1


def diagnostica():
    from astropy.wcs import WCS
    print(f"DIAGNOSTICA {NOME} - {datetime.now():%Y-%m-%d %H:%M}")
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [14.77083, -72.17583]
    w.wcs.crpix = [500, 500]
    w.wcs.cdelt = [-0.0002, 0.0002]
    try:
        t0 = time.time()
        ogg = M.oggetti_simbad(w, (1000, 1000))
        print(f"[ok] SIMBAD: {len(ogg)} oggetti nel campo di NGC 346 ({time.time() - t0:.0f} s)")
        for o in ogg[:8]:
            print(f"       {o['nome']:<25} {o['tipo']:<8} rif={o['riferimenti']}")
    except Exception as e:
        print(f"[ERRORE] SIMBAD: {type(e).__name__}: {e}")
    try:
        ogg = M.oggetti_openngc(w, (1000, 1000))
        print(f"[ok] catalogo OpenNGC incluso: {len(ogg)} oggetti nel campo di NGC 346: "
              + ", ".join(o["nome"] for o in ogg[:6]))
    except Exception as e:
        print(f"[ERRORE] catalogo OpenNGC: {type(e).__name__}: {e}")
    return 0


def _uscita_testo(nome):
    if sys.stdout is None or not hasattr(sys.stdout, "write"):
        sys.stdout = sys.stderr = open(Path.cwd() / f"{nome}_risultato.txt", "w", encoding="utf-8", buffering=1)


def main():
    if "--autotest" in sys.argv:
        _uscita_testo("autotest")
        sys.exit(autotest())
    if "--diagnostica" in sys.argv:
        _uscita_testo("diagnostica")
        sys.exit(diagnostica())
    if getattr(sys, "frozen", False) or sys.stderr is None:
        try:
            log = open(Path.home() / ".pipeline_fits.log", "a", encoding="utf-8", buffering=1)
            log.write(f"\n--- avvio {datetime.now():%Y-%m-%d %H:%M:%S} ---\n")
            sys.stdout = sys.stderr = log
        except OSError:
            sys.stdout = sys.stderr = open(os.devnull, "w")
    if sys.platform == "win32":
        try:
            import ctypes
            ctypes.windll.shcore.SetProcessDpiAwareness(1)
        except Exception:
            pass
    App().mainloop()


if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()
    main()
