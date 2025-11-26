import pandas as pd
import numpy as np
from ezc3d import c3d
import os

# === Fonctions utilitaires ===
def extract_marker_positions(file_path):
    c3d_data = c3d(file_path)
    marker_names = c3d_data["parameters"]["POINT"]["LABELS"]["value"]
    marker_data = np.array(c3d_data["data"]["points"])
    return {
        marker_names[i]: np.vstack((marker_data[0, i, :], marker_data[1, i, :], marker_data[2, i, :])).T
        for i in range(len(marker_names))
    }

def get_marker_names(file_path):
    return c3d(file_path)["parameters"]["POINT"]["LABELS"]["value"]

def calculate_MoS(file_path, heel_strike_frame, heel_off_frame, toe_off_frame):
    positions = extract_marker_positions(file_path)
    marker_names = get_marker_names(file_path)

    if "RM5" not in marker_names and "RM51" in marker_names:
        print("🛠 Remplacement de 'RM5' par 'RM51'")
        positions["RM5"] = positions["RM51"]

    left_heel_z = positions["LHEE"][heel_strike_frame, 2]
    right_heel_z = positions["RHEE"][heel_strike_frame, 2]
    side = "L" if left_heel_z < right_heel_z else "R"
    oppside = "R" if side == "L" else "L"

    COM_x = np.mean([positions[m][:, 0] for m in ['LPSI', 'RPSI', 'LASI', 'RASI']], axis=0)
    COM_y = np.mean([positions[m][:, 1] for m in ['LPSI', 'RPSI', 'LASI', 'RASI']], axis=0)
    COM_z = np.mean([positions[m][:, 2] for m in ['LPSI', 'RPSI', 'LASI', 'RASI']], axis=0)

    velocities_x = np.diff(COM_x) * 100
    velocities_y = np.diff(COM_y) * 100

    g = 9810
    l_z = np.abs(positions[f"{side}ANK"][:-1, 2] - COM_z[:-1])
    k = np.sqrt(g / (l_z + 1e-6))

    xCOM_x = COM_x[:-1] + velocities_x / k
    xCOM_y = COM_y[:-1] + velocities_y / k

    directionML = -1 if (positions[f"{side}M5"][0, 0] - positions[f"{oppside}M5"][0, 0]) < 0 else 1

    HEEL_y = positions[f"{side}HEE"][:, 1]
    TOE_y = positions[f"{side}TOE"][:, 1]
    ANKLE_x = positions[f"{side}ANK"][:, 0]
    M5_x = positions[f"{side}M5"][:, 0]

    MoS_AP_heel = HEEL_y[heel_strike_frame:heel_off_frame] - xCOM_y[heel_strike_frame:heel_off_frame]
    MoS_AP_toe = TOE_y[heel_off_frame:toe_off_frame] - xCOM_y[heel_off_frame:toe_off_frame]
    MoS_AP_Mean = np.mean([np.mean(MoS_AP_heel), np.mean(MoS_AP_toe)])

    MoS_ML_ankle = (ANKLE_x[heel_strike_frame:heel_off_frame] - xCOM_x[heel_strike_frame:heel_off_frame]) * directionML
    MoS_ML_m5 = (M5_x[heel_off_frame:toe_off_frame] - xCOM_x[heel_off_frame:toe_off_frame]) * directionML
    mos_ml_total = np.concatenate([MoS_ML_ankle, MoS_ML_m5])

    return {
        "MoS_AP_Heel_Min": np.min(MoS_AP_heel),
        "MoS_AP_Heel_Max": np.max(MoS_AP_heel),
        "MoS_AP_Heel_Mean": np.mean(MoS_AP_heel),
        "MoS_AP_Heel_SD": np.std(MoS_AP_heel),
        "MoS_AP_Toe_Min": np.min(MoS_AP_toe),
        "MoS_AP_Toe_Max": np.max(MoS_AP_toe),
        "MoS_AP_Toe_Mean": np.mean(MoS_AP_toe),
        "MoS_AP_Toe_SD": np.std(MoS_AP_toe),
        "MoS_AP_Mean": MoS_AP_Mean,
        "MoS_ML_Ankle_Min": np.min(MoS_ML_ankle),
        "MoS_ML_Ankle_Max": np.max(MoS_ML_ankle),
        "MoS_ML_Ankle_Mean": np.mean(MoS_ML_ankle),
        "MoS_ML_Ankle_SD": np.std(MoS_ML_ankle),
        "MoS_ML_M5_Min": np.min(MoS_ML_m5),
        "MoS_ML_M5_Max": np.max(MoS_ML_m5),
        "MoS_ML_M5_Mean": np.mean(MoS_ML_m5),
        "MoS_ML_M5_SD": np.std(MoS_ML_m5),
        "MoS_ML_Min": np.min(mos_ml_total),
        "MoS_ML_Max": np.max(mos_ml_total),
        "MoS_ML_Mean": np.mean(mos_ml_total),
        "MoS_ML_SD": np.std(mos_ml_total),
        "MoS_Heel_Strike_AP": MoS_AP_heel[0],
        "MoS_Heel_Strike_ML": mos_ml_total[0]
    }

# === Chargement events_CTL_27.csv ===
df = pd.read_csv("/Users/bastin/events_CTL_38.csv")
df.columns = df.columns.str.strip()
df = df[df['Sujet'] == 'CTL_38']

# Colonnes de résultats
mos_columns = [
    "MoS_AP_Heel_Min", "MoS_AP_Heel_Max", "MoS_AP_Heel_Mean", "MoS_AP_Heel_SD",
    "MoS_AP_Toe_Min", "MoS_AP_Toe_Max", "MoS_AP_Toe_Mean", "MoS_AP_Toe_SD", "MoS_AP_Mean",
    "MoS_ML_Ankle_Min", "MoS_ML_Ankle_Max", "MoS_ML_Ankle_Mean", "MoS_ML_Ankle_SD",
    "MoS_ML_M5_Min", "MoS_ML_M5_Max", "MoS_ML_M5_Mean", "MoS_ML_M5_SD",
    "MoS_ML_Min", "MoS_ML_Max", "MoS_ML_Mean", "MoS_ML_SD",
    "MoS_Heel_Strike_AP", "MoS_Heel_Strike_ML"
]
for col in mos_columns:
    df[col] = np.nan

# Répertoire des fichiers
base_dir = "/Users/bastin/Desktop/CTL_38/"
c3d_files = {f"Plat_{i}": f"CTL_38_Plat_Lunette_0{i}.c3d" for i in range(1, 5)}
c3d_files.update({f"Mixte_{i}": f"CTL_38_Mixte_0{i}.c3d" for i in range(1, 5)})

# Boucle de traitement
for idx, row in df.iterrows():
    surface = row["Surface"]
    essai = int(row["Essai"])
    key = f"{surface}_{essai}"
    if key not in c3d_files:
        print(f"⚠️ Fichier C3D introuvable pour {key}")
        continue
    c3d_path = os.path.join(base_dir, c3d_files[key])
    try:
        hs = int(row["Pose Talon"])
        ho = int(row["Décollement Talon"])
        to = int(row["Décollement Pied"])
        mos = calculate_MoS(c3d_path, hs, ho, to)
        for col in mos_columns:
            df.at[idx, col] = mos[col]
    except Exception as e:
        print(f"❌ Erreur ligne {idx} ({key}) : {e}")

# Sauvegarde
output_filename = "/Users/bastin/Desktop/MoS_results_CTL_38.csv"
df.to_csv(output_filename, index=False)
print(f"✅ Fichier généré : {output_filename}")
