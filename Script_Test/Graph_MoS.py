import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import spm1d
from ezc3d import c3d

# === Chemin de ton fichier events_CTL_27.csv ===
file_path = "/Users/bastin/events_CTL_27.csv"  # adapte si besoin

# === Chargement des données ===
df = pd.read_csv(file_path)
df.columns = df.columns.str.strip()
df = df[df['Sujet'] == 'CTL_27']

# Séparation par surface
df_plat = df[df["Surface"] == "Plat"]
df_mixte = df[df["Surface"] == "Mixte"]

# Chemins des fichiers C3D selon ton organisation
files_plat = {i: f"/Users/bastin/Desktop/CTL_27/CTL_27_Plat_lunette_0{i}.c3d" for i in range(1, 5)}
files_mixte = {i: f"/Users/bastin/Desktop/CTL_27/CTL_27_Mixte_0{i}.c3d" for i in range(1, 5)}

# Fonction utilitaire
def extract_marker_positions(file_path):
    c3d_data = c3d(file_path)
    labels = c3d_data["parameters"]["POINT"]["LABELS"]["value"]
    data = c3d_data["data"]["points"]
    return {label: np.vstack((data[0, i], data[1, i], data[2, i])).T for i, label in enumerate(labels)}

# Initialisation
mos_ap_plat, mos_ml_plat = [], []
mos_ap_mixte, mos_ml_mixte = [], []
x_norm = np.linspace(0, 1, 100)

# === Extraction des cycles et interpolation ===
for df_cond, files, mos_ap_list, mos_ml_list in [
    (df_plat, files_plat, mos_ap_plat, mos_ml_plat),
    (df_mixte, files_mixte, mos_ap_mixte, mos_ml_mixte)]:

    for _, row in df_cond.iterrows():
        try:
            essai = int(row["Essai"])
            hs = int(row["Pose Talon"])
            ho = int(row["Décollement Talon"])
            to = int(row["Décollement Pied"])
            lat = row["Latéralité"]

            file = files[essai]
            pos = extract_marker_positions(file)

            side = "L" if lat == "Gauche" else "R"
            oppside = "R" if side == "L" else "L"

            COM_x = np.mean([pos[m][:, 0] for m in ['LPSI', 'RPSI', 'LASI', 'RASI']], axis=0)
            COM_y = np.mean([pos[m][:, 1] for m in ['LPSI', 'RPSI', 'LASI', 'RASI']], axis=0)
            COM_z = np.mean([pos[m][:, 2] for m in ['LPSI', 'RPSI', 'LASI', 'RASI']], axis=0)

            vx = np.diff(COM_x) * 100
            vy = np.diff(COM_y) * 100
            lz = np.abs(pos[f"{side}ANK"][:-1, 2] - COM_z[:-1])
            k = np.sqrt(9810 / (lz + 1e-6))
            xCOM_x = COM_x[:-1] + vx / k
            xCOM_y = COM_y[:-1] + vy / k

            dirML = -1 if pos[f"{side}M5"][0, 0] - pos[f"{oppside}M5"][0, 0] < 0 else 1

            mos_ap = np.concatenate([
                (pos[f"{side}HEE"][hs:ho, 1] - xCOM_y[hs:ho]),
                (pos[f"{side}TOE"][ho:to, 1] - xCOM_y[ho:to])
            ])
            mos_ml = np.concatenate([
                (pos[f"{side}ANK"][hs:ho, 0] - xCOM_x[hs:ho]) * dirML,
                (pos[f"{side}M5"][ho:to, 0] - xCOM_x[ho:to]) * dirML
            ])

            mos_ap_list.append(np.interp(x_norm, np.linspace(0, 1, len(mos_ap)), mos_ap))
            mos_ml_list.append(np.interp(x_norm, np.linspace(0, 1, len(mos_ml)), mos_ml))

        except Exception as e:
            print(f"Erreur essai {essai} {lat}: {e}")

# Conversion en arrays
ctl27_mos_ap_plat = np.array(mos_ap_plat)
ctl27_mos_ap_mixte = np.array(mos_ap_mixte)
ctl27_mos_ml_plat = np.array(mos_ml_plat)
ctl27_mos_ml_mixte = np.array(mos_ml_mixte)

# === SPM et tracés MoS AP ===
ti_ap = spm1d.stats.ttest2(ctl27_mos_ap_plat, ctl27_mos_ap_mixte)
spmi_ap = ti_ap.inference(alpha=0.05, two_tailed=True, interp=True)

fig, axs = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
x = x_norm * 100

axs[0].plot(x, ctl27_mos_ap_plat.mean(axis=0), color='blue', label='Plat')
axs[0].fill_between(x,
                    ctl27_mos_ap_plat.mean(axis=0) - ctl27_mos_ap_plat.std(axis=0),
                    ctl27_mos_ap_plat.mean(axis=0) + ctl27_mos_ap_plat.std(axis=0),
                    color='blue', alpha=0.3)
axs[0].plot(x, ctl27_mos_ap_mixte.mean(axis=0), color='red', label='Mixte')
axs[0].fill_between(x,
                    ctl27_mos_ap_mixte.mean(axis=0) - ctl27_mos_ap_mixte.std(axis=0),
                    ctl27_mos_ap_mixte.mean(axis=0) + ctl27_mos_ap_mixte.std(axis=0),
                    color='red', alpha=0.3)
axs[0].set_ylabel("MoS AP (mm)")
axs[0].set_title("MoS Antéro-Postérieur CTL_27 (Moyenne ± ET)")
axs[0].legend()
axs[0].grid(True)

spmi_ap.plot(ax=axs[1])
spmi_ap.plot_threshold_label(ax=axs[1])
axs[1].set_xlabel("% du cycle d'appui")
axs[1].set_ylabel("t-value")
axs[1].set_title("SPM t-indépendant : MoS AP Plat vs Mixte")
axs[1].grid(True)
plt.tight_layout()
plt.show()

# === SPM et tracés MoS ML ===
ti_ml = spm1d.stats.ttest2(ctl27_mos_ml_plat, ctl27_mos_ml_mixte)
spmi_ml = ti_ml.inference(alpha=0.05, two_tailed=True, interp=True)

fig, axs = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

axs[0].plot(x, ctl27_mos_ml_plat.mean(axis=0), color='blue', label='Plat')
axs[0].fill_between(x,
                    ctl27_mos_ml_plat.mean(axis=0) - ctl27_mos_ml_plat.std(axis=0),
                    ctl27_mos_ml_plat.mean(axis=0) + ctl27_mos_ml_plat.std(axis=0),
                    color='blue', alpha=0.3)
axs[0].plot(x, ctl27_mos_ml_mixte.mean(axis=0), color='red', label='Mixte')
axs[0].fill_between(x,
                    ctl27_mos_ml_mixte.mean(axis=0) - ctl27_mos_ml_mixte.std(axis=0),
                    ctl27_mos_ml_mixte.mean(axis=0) + ctl27_mos_ml_mixte.std(axis=0),
                    color='red', alpha=0.3)
axs[0].set_ylabel("MoS ML (mm)")
axs[0].set_title("MoS Medio-Latéral CTL_27 (Moyenne ± ET)")
axs[0].legend()
axs[0].grid(True)

spmi_ml.plot(ax=axs[1])
spmi_ml.plot_threshold_label(ax=axs[1])
axs[1].set_xlabel("% du cycle d'appui")
axs[1].set_ylabel("t-value")
axs[1].set_title("SPM t-indépendant : MoS ML Plat vs Mixte")
axs[1].grid(True)
plt.tight_layout()
plt.show()

