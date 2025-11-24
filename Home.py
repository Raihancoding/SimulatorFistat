import streamlit as st
from math import factorial, log
from collections import Counter
import itertools
import math
import os
import matplotlib.pyplot as plt
import pandas as pd

# Diasumsikan Konstanta Boltzmann k_B = 1 untuk perbandingan nilai
KB = 1.0 

# ----------------------------------------
# FUNGSI-FUNGSI UTILITY & PERHITUNGAN (TIDAK BERUBAH)
# ----------------------------------------
def nCr(n, r):
    """Menghitung nCr (n choose r) secara aman."""
    if r < 0 or r > n: return 0
    if hasattr(math, 'comb'): return math.comb(n, r)
    return factorial(n) // (factorial(r) * factorial(n - r))

# Microstate (W) dan Macrostate (M) untuk 1 level g (Digunakan di Mode Umum)
def microstates_MB(N, g): return g**N
def microstates_BE(N, g): return nCr(N + g - 1, N)
def microstates_FD(N, g): return 0 if N > g else nCr(g, N)
def macrostates_MB(N, g): return nCr(N + g - 1, N) 
def macrostates_BE(N, g): return microstates_BE(N, g)
def macrostates_FD(N, g): return microstates_FD(N, g)

MAX_SHOW = 5000

# Fungsi untuk membuat daftar Microstate (digunakan di Mode Umum)
def list_MB(N, g):
    total = g**N
    if total > MAX_SHOW: return None, total
    return list(itertools.product(range(0, g), repeat=N)), total

def list_BE(N, g):
    total = nCr(N + g - 1, N)
    if total > MAX_SHOW: return None, total
    return list(itertools.combinations_with_replacement(range(0, g), N)), total

def list_FD(N, g):
    total = 0 if N > g else nCr(g, N)
    if total > MAX_SHOW: return None, total
    return list(itertools.combinations(range(0, g), N)), total

def state_distribution(E_list):
    c = Counter(E_list)
    levels = sorted(c.keys())
    distrib = [c[l] for l in levels]
    return levels, distrib

def find_macrostates(N, E_levels, E_total):
    """Mencari semua Macrostate (n1, n2, ...) yang memenuhi N dan E_total (Menggunakan Backtracking)."""
    results = []
    L = len(E_levels)
    def backtrack(i, remN, remE, current):
        if i == L:
            if remN == 0 and remE == 0:
                results.append(current.copy())
            return
        
        max_n = remN
        
        for n in range(max_n + 1):
            neededE = n * E_levels[i]
            if neededE > remE:
                break
            
            current.append(n)
            backtrack(i + 1, remN - n, remE - neededE, current)
            current.pop()
            
    backtrack(0, N, E_total, [])
    return results

# NEW FUNCTION: Menganalisis semua energi yang mungkin
def analyze_all_energy_states(N, E_levels_int, g_levels_search, stat_mode):
    if not E_levels_int or N <= 0:
        return None, 0, 0
    
    E_min = N * min(E_levels_int)
    E_max = N * max(E_levels_int)
    
    energy_density = {} 
    
    for E_total in range(E_min, E_max + 1):
        macrostates = find_macrostates(N, E_levels_int, E_total)
        
        total_W = 0
        for macro in macrostates:
            g_used = g_levels_search
            if len(g_used) < len(macro): g_used = (g_levels_search + [1] * (len(macro) - len(g_levels_search)))
            elif len(g_used) > len(macro): g_used = g_used[:len(macro)]
            
            W_macro = 0
            
            # Hitung W Macrostate sesuai statistik yang dipilih
            if stat_mode == "Maxwell–Boltzmann":
                numerator = factorial(N)
                denom = 1
                g_factor = 1
                for ni, gi in zip(macro, g_used):
                    denom *= factorial(ni)
                    g_factor *= (gi ** ni)
                W_macro = (numerator // denom) * g_factor
            elif stat_mode == "Bose–Einstein":
                W_macro = 1
                for ni, gi in zip(macro, g_used):
                    W_macro *= nCr(ni + gi - 1, ni)
            else: # Fermi-Dirac
                valid_FD = all(ni <= gi for ni, gi in zip(macro, g_used))
                if valid_FD:
                    W_macro = 1
                    for ni, gi in zip(macro, g_used):
                        W_macro *= nCr(gi, ni)
                else:
                    W_macro = 0
            
            total_W += W_macro
        
        if total_W > 0:
            energy_density[E_total] = total_W

    return energy_density, E_min, E_max

# ATOM_DATA 
ATOM_DATA = {
    "H": [(0, 2)], "He": [(0, 2), (2, 4)], 
    "Uranium": [(0, 2), (1, 2), (2, 4), (3, 4), (4, 6), (5, 6), (6, 8), (7, 8)]
}

# Helper parse functions
def parse_energy_text(txt):
    try:
        return [float(x.strip()) for x in txt.split(",") if x.strip() != ""]
    except: return None
def parse_int_list(txt):
    try:
        return [int(x.strip()) for x in txt.split(",") if x.strip() != ""]
    except: return None

# ----------------------------------------
# Streamlit UI
# ----------------------------------------
st.set_page_config(layout="wide")
st.title(" Simulator Statistik Terpadu")

# ====== BLOK PANDUAN SELAMAT DATANG BARU ======

# Memisahkan title dari konten
if not st.session_state.get('app_started'):
    st.markdown("""
        Selamat datang di **Simulator Statistik Terpadu**!
        
        Aplikasi ini dirancang untuk memvisualisasikan Macrostate, Microstate, Entropi ($S$),
        dan Probabilitas ($P$) untuk sistem kuantum dasar (MB, BE, FD).
    """)
    
    st.subheader("Petunjuk Penggunaan Cepat:")
    
    st.markdown("""
    1.  **Konfigurasi Sistem (Sidebar Kiri):** Pilih Statistik, masukkan **N** (jumlah partikel), dan tentukan **Level Energi Unik** serta **Degenerasi ($g_i$)**.
    2.  **Pilih Mode Analisis:**
        * **Mode Umum:** Biarkan semua kotak centang di sidebar nonaktif untuk mendapatkan $\Omega$ dan $M$ total berdasarkan $N$ dan $g$.
        * **Makrostate Otomatis:** Aktifkan kotak centang ini dan masukkan $E_{\\text{total}}$ untuk mencari semua solusi $n_i$ yang mungkin.
        * **Plot Kerapatan Keadaan:** Aktifkan kotak centang ini untuk memvisualisasikan Microstate total ($\Omega$) untuk **semua** kemungkinan energi.
    3.  **Mulai:** Klik tombol **"Hitung 🚀"** di bagian bawah sidebar untuk melihat hasil di area ini.
    """)
    
    st.warning("Pastikan jumlah level energi dan degenerasi ($g_i$) yang dimasukkan konsisten.")

# Mengatur state untuk menghilangkan panduan setelah tombol ditekan
if 'app_started' not in st.session_state:
    st.session_state['app_started'] = False

# ====== END OF BLOK PANDUAN SELAMAT DATANG BARU ======


# --- SIDEBAR (KONFIGURASI) ---
st.sidebar.header("1. Konfigurasi Sistem ⚙️")

stat = st.sidebar.selectbox("Mode Statistik", ["Maxwell–Boltzmann", "Bose–Einstein", "Fermi–Dirac"])

# --- PENAMBAHAN RUMUS STATISTIK DI SIDEBAR ---
if stat == "Maxwell–Boltzmann":
    st.sidebar.markdown("**MB: Partikel Dapat Dibedakan**")
    st.sidebar.latex(r"\Omega_{\text{Multi-level}} = N! \prod_i \frac{g_i^{n_i}}{n_i!}")
    st.sidebar.latex(r"M = \binom{N+L-1}{N} \quad \text{(L = Jumlah Level)}")
elif stat == "Bose–Einstein":
    st.sidebar.markdown("**BE: Partikel Identik (Boson)**")
    st.sidebar.latex(r"\Omega_{\text{Multi-level}} = \prod_i \binom{n_i + g_i - 1}{n_i}")
    st.sidebar.latex(r"M = \binom{N+L-1}{N} \quad \text{(L = Jumlah Level)}")
else: # Fermi-Dirac
    st.sidebar.markdown("**FD: Partikel Identik (Fermion)**")
    st.sidebar.latex(r"\Omega_{\text{Multi-level}} = \prod_i \binom{g_i}{n_i} \quad \text{(Jika } n_i \leq g_i \text{)}")
    st.sidebar.latex(r"M = \binom{N+L-1}{N} \quad \text{(L = Jumlah Level)}")

st.sidebar.subheader("Pilih Atom / Manual")
atom_choice = st.sidebar.selectbox("Atom", ["Manual"] + list(ATOM_DATA.keys()))

deg_option = st.sidebar.radio("Tipe degenerasi ($g$):", ["Sama", "Per-level"])

# Logic Atom/Manual
if atom_choice != "Manual":
    level_data = ATOM_DATA[atom_choice]
    st.sidebar.markdown("**Level Energi (E, g)**")
    text_block = "".join(f"E={E}, g={gg}\n" for E, gg in level_data)
    st.sidebar.text(text_block)
    E_list_from_atom = [float(E) for E, gg in level_data for _ in range(gg)]
    use_atom = st.sidebar.checkbox("Gunakan data atom ini", value=True)
else:
    E_list_from_atom = []
    use_atom = False

st.sidebar.markdown("**Input Manual**")
manual_energy_text = st.sidebar.text_area("Level Energi Unik (pisah koma)", "0,1,2")
manual_N = st.sidebar.number_input("Jumlah partikel N", min_value=1, value=3)

if deg_option == "Per-level":
    deg_text = st.sidebar.text_input("Degenerasi per-level (pisah koma)", "1,1,1")
else:
    deg_text = None
    g_manual = st.sidebar.number_input("g (total keadaan energi)", min_value=1, value=5)

# Mode Khusus Nx
st.sidebar.markdown("---")
st.sidebar.subheader("Mode Khusus: Distribusi $N_i$")
mode_nx = st.sidebar.checkbox("Aktifkan mode $N_i$ (Hitung $\Omega$ dari Macrostate tunggal)", value=False)
if mode_nx:
    energy_levels_text = st.sidebar.text_area("Energi Level ($E_i$) [Unik]", "1,2,3,4")
    degeneracy_levels_text = st.sidebar.text_input("Degenerasi Level ($g_i$)", "1,3,3,5")
    nx_text = st.sidebar.text_input("Populasi per Level ($N_i$)", "1,1,2,1")

# Makrostate Otomatis
st.sidebar.markdown("---")
st.sidebar.subheader("Pencarian Makrostate Otomatis")
use_auto_macro = st.sidebar.checkbox("Aktifkan Makrostate Otomatis (Cari $n_i$ berdasarkan $E_{total}$)", value=False)
E_total_input = st.sidebar.number_input("Total Energi ($E_{total}$)", min_value=0, value=12) if use_auto_macro else None

# NEW FEATURE INPUT
st.sidebar.markdown("---")
st.sidebar.subheader("Analisis Kerapatan Keadaan (Density of States)")
use_all_energy_plot = st.sidebar.checkbox("Plot $\Omega$ vs Semua $E_{\text{total}}$ yang mungkin", value=False)


do = st.sidebar.button("Hitung 🚀")

# ----------------------------------------
# MAIN COMPUTATION & OUTPUT
# ----------------------------------------
if do:
    # Set state agar panduan selamat datang hilang
    st.session_state['app_started'] = True 
    
    # --- Mode Nx ---
    if mode_nx:
        E_levels = parse_energy_text(energy_levels_text)
        g_levels = parse_int_list(degeneracy_levels_text)
        Nx_list = parse_int_list(nx_text)

        if E_levels is None or g_levels is None or Nx_list is None or not (len(E_levels) == len(g_levels) == len(Nx_list)):
            st.error("Input mode Nx salah atau panjangnya tidak sama.")
            st.stop()
            
        N_total = sum(Nx_list)
        
        # Hitung W_nx_mb, W_nx_be, W_nx_fd (Rumus Macrostate Multi-Level)
        numerator = factorial(N_total)
        denominator = 1
        g_term_mb = 1
        W_be_term = 1
        W_fd_term = 1
        is_fd_valid = True
        
        for Nx, gval in zip(Nx_list, g_levels):
            denominator *= factorial(Nx)
            g_term_mb *= (gval ** Nx)
            W_be_term *= nCr(Nx + gval - 1, Nx)
            if Nx > gval: is_fd_valid = False
            W_fd_term *= nCr(gval, Nx)

        W_nx_mb = (numerator // denominator) * g_term_mb
        W_nx_be = W_be_term
        W_nx_fd = W_fd_term if is_fd_valid else 0
        
        W_display = W_nx_mb if stat == "Maxwell–Boltzmann" else W_nx_be if stat == "Bose–Einstein" else W_nx_fd
        
        # Hitung Entropi
        S = KB * log(W_display) if W_display > 0 else 0

        st.header("✅ Hasil Mode Khusus ($N_i$)")
        # Menggunakan sintaks f-string yang benar
        st.info(f"Dihitung untuk Macrostate $\\mathbf{{({', '.join(map(str, Nx_list))})}}$ dengan $N={N_total}$")
        
        col1, col2, col3 = st.columns(3)
        col1.metric("Microstate (MB)", f"{W_nx_mb:,}")
        col2.metric("Microstate (BE)", f"{W_nx_be:,}")
        col3.metric("Microstate (FD)", f"{W_nx_fd:,}")
        
        st.metric("Entropi ($S/k_B$)", f"{S/KB:.4f}")

        st.success(f"Microstate $\Omega$ Sesuai Statistik ({stat}): **{W_display:,}**")

        st.subheader("Grafik Populasi")
        fig_nx, ax_nx = plt.subplots(figsize=(8, 4))
        ax_nx.bar([f"E={e}" for e in E_levels], Nx_list, color='darkorange')
        ax_nx.set_xlabel("Level Energi")
        ax_nx.set_ylabel("Jumlah Partikel $N_i$")
        ax_nx.set_title("Distribusi Partikel per Level ($N_i$)")
        st.pyplot(fig_nx)

        st.stop()
    
    # --- PERSIAPAN UNTUK PERHITUNGAN UMUM ATAU MAKROSTATE OTOMATIS/PLOT ENERGI ---
    
    # Logic penentuan N, g, E_list, E_levels_search, g_levels_search
    if use_atom and atom_choice != "Manual":
        E_list = E_list_from_atom.copy()
        N = len(E_list)
        atom_level_pairs = ATOM_DATA[atom_choice]
        E_levels_search = [int(pair[0]) for pair in atom_level_pairs]
        g_levels_search = [int(pair[1]) for pair in atom_level_pairs]
        g = sum(g_levels_search)
    else:
        parsed = parse_energy_text(manual_energy_text)
        if parsed is None: st.error("Gagal mem-parse energi manual."); st.stop()
        E_list = parsed
        N = int(manual_N)
        E_levels_search = sorted(list(set(E_list)))
        if deg_option == "Per-level":
            parsed_deg = parse_int_list(deg_text)
            if parsed_deg is None: st.error("Gagal mem-parse degenerasi per-level."); st.stop()
            g = sum(parsed_deg)
            if len(parsed_deg) == len(E_levels_search): g_levels_search = parsed_deg
            else: g_levels_search = [1] * len(E_levels_search)
        else:
            g = int(g_manual)
            g_levels_search = [g] * len(E_levels_search)
            
    if N <= 0 or g <= 0: st.error("N dan g harus positif."); st.stop()
    if not use_auto_macro and not use_all_energy_plot and len(E_list) != N: st.error(f"Jumlah energi ({len(E_list)}) harus sama dengan N ({N}) untuk mode umum."); st.stop()
    
    E_levels_int = [int(x) for x in E_levels_search]

    # --- TATA LETAK KOLOM HASIL UTAMA ---
    col_visual, col_summary = st.columns([3, 1])

    # ----------------------------------------
    # KOLOM KANAN: RINGKASAN & RUMUS
    # ----------------------------------------
    with col_summary:
        st.subheader("📌 Ringkasan Hasil")

        # ----------------------------------------
        # 4. PLOT SEMUA ENERGI (OVERRIDE DISPLAY UTAMA)
        # ----------------------------------------
        if use_all_energy_plot:
            energy_density, E_min, E_max = analyze_all_energy_states(N, E_levels_int, g_levels_search, stat)
            
            if not energy_density:
                st.error("Tidak ada Microstate yang ditemukan dalam rentang energi ini.")
                st.stop()
            
            E_values = list(energy_density.keys())
            Omega_values = list(energy_density.values())
            
            W = sum(Omega_values)
            S = KB * log(W) if W > 0 else 0
            
            st.metric("Microstate Total ($\Omega$)", f"{W:,}")
            st.metric("Entropi ($S/k_B$)", f"{S/KB:.4f}")
            st.metric("Rentang Energi", f"{E_min} - {E_max}")
            st.metric("Total Partikel ($N$)", N)
            
            formula = r"\Omega_{\text{total}} = \sum_{E} \Omega(E)"
            
            st.markdown("---")
            with st.expander("Lihat Rumus Perhitungan Detail"):
                st.latex(f"Statistik: {stat}")
                st.latex(r"S = k_B \ln \Omega \quad \text{(Entropi Boltzmann)}")
                st.latex(formula)
                st.write(f"Level Energi Unik: {E_levels_search}")
                st.warning("Perhitungan umum dan makrostate otomatis dilewati.")
            
            # Pindah ke Kolom Visual
            with col_visual:
                st.header("📈 Kerapatan Keadaan Diskrit ($\Omega$ vs $E_{\text{total}}$)")
                # Baris ini diperbaiki untuk sintaks yang bersih:
                st.info(f"Menghitung semua Macrostate untuk $E_{{\\text{{total}}}}$ dari {E_min} hingga {E_max}.") 
                
                fig_dos, ax_dos = plt.subplots(figsize=(9, 5))
                ax_dos.bar([str(e) for e in E_values], Omega_values, color='purple', alpha=0.7)
                ax_dos.set_xlabel("Energi Total ($E_{\\text{total}}$)")
                ax_dos.set_ylabel(f"Total Microstate $\Omega(E)$ ({stat})")
                ax_dos.set_title("Density of States (Kerapatan Keadaan)")
                st.pyplot(fig_dos)

            st.stop()


        # ----------------------------------------
        # 2. RUN MAKROSTATE OTOMATIS (Jika diaktifkan)
        # ----------------------------------------
        if use_auto_macro:
            macrostates = find_macrostates(N, E_levels_int, int(E_total_input))
            
            if len(macrostates) == 0: 
                st.error("Tidak ada makrostate yang memenuhi E_total."); st.stop()

            macro_results = []
            total_W_MB, total_W_BE, total_W_FD = 0, 0, 0
            
            for macro in macrostates:
                g_used = g_levels_search
                if len(g_used) < len(macro): g_used = (g_levels_search + [1] * (len(macro) - len(g_levels_search)))
                elif len(g_used) > len(macro): g_used = g_used[:len(macro)]

                # Hitung W untuk ketiga statistik
                numerator = factorial(N)
                denom = 1
                g_factor = 1
                W_BE_macro = 1
                W_FD_macro = 1
                valid_FD = True
                
                for ni, gi in zip(macro, g_used):
                    denom *= factorial(ni)
                    g_factor *= (gi ** ni)
                    W_BE_macro *= nCr(ni + gi - 1, ni)
                    if ni > gi: valid_FD = False
                    W_FD_macro *= nCr(gi, ni)

                W_MB_macro = (numerator // denom) * g_factor
                W_FD_macro = W_FD_macro if valid_FD else 0
                
                total_W_MB += W_MB_macro
                total_W_BE += W_BE_macro
                total_W_FD += W_FD_macro

                macro_results.append({
                    "Konfigurasi ($n_1,n_2,...$)": str(tuple(macro)),
                    "W_MB": W_MB_macro, "W_BE": W_BE_macro, "W_FD": W_FD_macro,
                    "P_MB": 0.0, "P_BE": 0.0, "P_FD": 0.0,
                })
            
            W = total_W_MB if stat == "Maxwell–Boltzmann" else total_W_BE if stat == "Bose–Einstein" else total_W_FD
            M = len(macrostates)
            
            # Perhitungan Probabilitas Final
            for res in macro_results:
                if total_W_MB > 0: res['P_MB'] = res['W_MB'] / total_W_MB
                if total_W_BE > 0: res['P_BE'] = res['W_BE'] / total_W_BE
                if total_W_FD > 0: res['P_FD'] = res['W_FD'] / total_W_FD
            
            S = KB * log(W) if W > 0 else 0
            
            st.metric("Microstate Total ($\Omega$)", f"{W:,}")
            st.metric("Macrostate Total ($M$)", f"{M:,}")
            st.metric("Entropi ($S/k_B$)", f"{S/KB:.4f}")
            st.metric("Total Partikel ($N$)", N)
            
            formula = r"\Omega_{\text{total}} = \sum_{n_i} W(n_i)"
            
        # ----------------------------------------
        # 3. PERHITUNGAN UMUM (NON-MAKROSTATE OTOMATIS & NON-PLOT)
        # ----------------------------------------
        else:
            if stat == "Maxwell–Boltzmann":
                W = microstates_MB(N, g); M = macrostates_MB(N, g)
                micro_list, total = list_MB(N, g)
                formula = r"W = g^N,\quad M = \binom{N+g-1}{N}"
            elif stat == "Bose–Einstein":
                W = microstates_BE(N, g); M = macrostates_BE(N, g)
                micro_list, total = list_BE(N, g)
                formula = r"W = \binom{N+g-1}{N},\quad M=W"
            else:
                W = microstates_FD(N, g); M = macrostates_FD(N, g)
                micro_list, total = list_FD(N, g)
                formula = r"W = \binom{g}{N},\quad M=W"
            
            S = KB * log(W) if W > 0 else 0

            st.metric("Microstate ($W$)", f"{W:,}")
            st.metric("Macrostate ($M$)", f"{M:,}")
            st.metric("Entropi ($S/k_B$)", f"{S/KB:.4f}")
            st.metric("Total Partikel ($N$)", N)
            
        # Rumus dan Ringkasan Teknis
        st.markdown("---")
        with st.expander("Lihat Rumus Perhitungan Detail"):
            st.latex(f"Statistik: {stat}")
            st.latex(r"S = k_B \ln \Omega \quad \text{(Entropi Boltzmann)}")
            st.latex(formula)
            st.write(f"Degenerasi Total $g$: {g}")


    # ----------------------------------------
    # KOLOM KIRI: VISUALISASI & TABEL DATA
    # ----------------------------------------
    with col_visual:
        st.header("📊 Analisis Data")
        
        # Penentuan Distribusi untuk Grafik
        levels, distrib = state_distribution(E_list)
        if use_auto_macro:
             st.info(f"Pencarian Selesai: Ditemukan {M} makrostate untuk $E_{{total}}={E_total_input}$.")
        
        # --- TAB GROUP: VISUALISASI & DETAIL ---
        tab_viz, tab_macro, tab_micro = st.tabs(["Visualisasi & Distribusi", "Tabel Makrostate ($N_i$)", "Daftar Mikrostate"])

        # TAB 1: VISUALISASI DAN DISTRIBUSI
        with tab_viz:
            st.subheader("Grafik Distribusi Energi")
            st.info(f"Level Energi Unik ($E_i$): {levels} | Distribusi partikel per level ($n_i$): $\\mathbf{{{tuple(distrib)}}}$")

            fig, ax = plt.subplots(figsize=(8, 4))
            ax.bar([str(lv) for lv in levels], distrib, color='skyblue')
            ax.set_xlabel("Level Energi $E_i$")
            ax.set_ylabel("Jumlah Partikel $n_i$")
            ax.set_title("Populasi Partikel per Level Energi")
            st.pyplot(fig)
            
        # TAB 2: TABEL MAKROSTATE OTOMATIS
        with tab_macro:
            if use_auto_macro:
                st.subheader(f"Detail Solusi Makrostate ($N_i$)")
                df_macro = pd.DataFrame(macro_results)
                
                # Memilih kolom W dan P yang relevan
                W_col = 'W_MB' if stat == "Maxwell–Boltzmann" else 'W_BE' if stat == "Bose–Einstein" else 'W_FD'
                P_col = 'P_MB' if stat == "Maxwell–Boltzmann" else 'P_BE' if stat == "Bose–Einstein" else 'P_FD'
                
                df_display = df_macro[['Konfigurasi ($n_1,n_2,...$)', W_col, P_col]].copy()
                df_display = df_display.rename(columns={W_col: 'Microstate (W)', P_col: 'Probabilitas (P)'})

                st.dataframe(df_display, height=400, use_container_width=True, 
                             column_config={"Probabilitas (P)": st.column_config.ProgressColumn("Probabilitas (P)", format="%.4f")})
            else:
                st.warning("Mode ini menampilkan daftar semua solusi $N_i$. Aktifkan **Makrostate Otomatis** di sidebar.")
            
        # TAB 3: DAFTAR MICROSTATE
        with tab_micro:
            st.subheader("Daftar Konfigurasi Mikro")
            # Pastikan micro_list ada dan mode Makrostate Otomatis TIDAK aktif
            if 'micro_list' in locals() and not use_auto_macro and not use_all_energy_plot and micro_list is not None:
                st.info(f"Menampilkan {len(micro_list):,} mikrostate (total = {total:,})")
                df_micro = pd.DataFrame({
                    "No.": range(1, len(micro_list) + 1),
                    "Konfigurasi": [str(m) for m in micro_list]
                })
                st.data_editor(df_micro, use_container_width=True, hide_index=True)
            else:
                 st.warning("Daftar Mikrostate terlalu besar (batas 5000), mode Otomatis/Plot Energi sedang aktif.")