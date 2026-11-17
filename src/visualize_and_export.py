#!/usr/bin/env python3
"""
visualize_and_export.py
Generates:
1. A 3-panel high-DPI mission operations chart (Global J2 Ground Tracks, 24-Hr Station
   Contact & Contention Gantt Timeline, and Peak X-Band Link Margin / Doppler Profile).
2. An OGC KML 2.2 3D orbital trajectory and ground-station overlay for Google Earth Pro.
"""

import os
import json
import sqlite3
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = os.path.join(BASE_DIR, "data", "raw")
PROC_DIR = os.path.join(BASE_DIR, "data", "processed")
DB_PATH = os.path.join(PROC_DIR, "space_ops_schedule.db")
PNG_OUT = os.path.join(PROC_DIR, "leo_mission_schedule_analysis.png")
KML_OUT = os.path.join(PROC_DIR, "leo_orbital_tracks_and_stations.kml")


def split_antimeridian(lons: list, lats: list):
    """Inserts NaNs across +/-180 deg longitude jumps to prevent horizontal map wrap lines."""
    out_lon, out_lat = [lons[0]], [lats[0]]
    for i in range(1, len(lons)):
        if abs(lons[i] - lons[i - 1]) > 180.0:
            out_lon.append(np.nan)
            out_lat.append(np.nan)
        out_lon.append(lons[i])
        out_lat.append(lats[i])
    return out_lon, out_lat


def export_visuals_and_kml():
    with open(os.path.join(PROC_DIR, "orbital_tracks.json"), "r", encoding="utf-8") as f:
        tracks = json.load(f)
    with open(os.path.join(PROC_DIR, "best_pass_profile.json"), "r", encoding="utf-8") as f:
        best_pass = json.load(f)
    with open(os.path.join(RAW_DIR, "ground_stations.json"), "r", encoding="utf-8") as f:
        stations = json.load(f)

    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    cur.execute("""
        SELECT station_id, sat_id, aos_sec, los_sec, schedule_status
        FROM contact_windows
        ORDER BY aos_sec ASC;
    """)
    contacts = cur.fetchall()
    conn.close()

    sat_colors = {
        "VANGUARD-SAR-1": "#2563eb",
        "VANGUARD-EO-2": "#16a34a",
        "SENTINEL-SIG-3": "#d97706",
        "PATHFINDER-TD-4": "#7c3aed"
    }

    # 1. Generate 3-Panel Space Operations Dashboard
    fig = plt.figure(figsize=(15, 10.5), dpi=220)
    fig.suptitle(
        "LEO Constellation J2 Pass Scheduler, X-Band Link Budget & Dish Deconfliction (24-Hr Horizon)",
        fontsize=13, fontweight="bold"
    )

    # Panel A: Global Sub-Satellite Ground Tracks (100-Min Window) & Ground Network
    ax_map = fig.add_subplot(2, 2, 1)
    for sat_id, trk in tracks.items():
        clean_lon, clean_lat = split_antimeridian(trk["lon_deg"], trk["lat_deg"])
        ax_map.plot(
            clean_lon, clean_lat,
            color=sat_colors.get(sat_id, "#334155"),
            linewidth=1.5, label=f"{sat_id} (P{trk['priority_tier']})"
        )

    for gs in stations:
        ax_map.scatter(gs["lon_deg"], gs["lat_deg"], c="#dc2626", marker="^", s=85, zorder=5)
        ax_map.text(
            gs["lon_deg"] + 3.0, gs["lat_deg"] - 4.5, gs["station_id"].replace("GS-", ""),
            fontsize=8, fontweight="bold", color="#991b1b"
        )

    ax_map.axhline(0.0, color="#94a3b8", linestyle="--", linewidth=0.7)
    ax_map.set_xlim(-180, 180)
    ax_map.set_ylim(-90, 90)
    ax_map.set_xlabel("Longitude (deg WGS84)")
    ax_map.set_ylabel("Latitude (deg WGS84)")
    ax_map.set_title("Constellation J2 Sub-Satellite Ground Tracks (100-Min Slice) & Ground Sites")
    ax_map.grid(True, linestyle=":", alpha=0.6)
    ax_map.legend(loc="lower left", fontsize=7.5)

    # Panel B: Peak Pass X-Band RF Link Margin & Doppler Curve
    ax_rf = fig.add_subplot(2, 2, 2)
    t_rel = best_pass["rel_time_min"]
    ax_rf.plot(t_rel, best_pass["link_margin_db"], color="#16a34a", linewidth=2.0, label="X-Band Link Margin (dB)")
    ax_rf.plot(t_rel, [e / 5.0 for e in best_pass["el_deg"]], color="#0f172a", linestyle="--",
               linewidth=1.2, label="Elevation / 5 (deg)")
    ax_rf.axhline(3.0, color="#dc2626", linestyle=":", linewidth=1.1, label="3.0 dB Min Margin Cutoff")

    ax_dopp = ax_rf.twinx()
    ax_dopp.plot(t_rel, best_pass["doppler_khz"], color="#2563eb", linewidth=1.4, alpha=0.8, label="8.2 GHz Doppler (kHz)")
    ax_dopp.set_ylabel("X-Band Doppler Shift (kHz)", color="#2563eb")
    ax_dopp.tick_params(axis="y", labelcolor="#2563eb")

    ax_rf.set_title(
        f"Dynamic RF Link Profile: {best_pass['sat_id']} over {best_pass['station_id']} "
        f"(Max El {best_pass['max_el_deg']} deg)"
    )
    ax_rf.set_xlabel("Pass Elapsed Time Since AOS (minutes)")
    ax_rf.set_ylabel("Eb/N0 Link Margin (dB)  |  Scaled El")
    ax_rf.grid(True, linestyle=":", alpha=0.6)
    l1, lb1 = ax_rf.get_legend_handles_labels()
    l2, lb2 = ax_dopp.get_legend_handles_labels()
    ax_rf.legend(l1 + l2, lb1 + lb2, loc="lower center", fontsize=7.5, ncol=2)

    # Panel C: 24-Hour Multi-Station Contact Allocation & Dish Contention Gantt
    ax_gantt = fig.add_subplot(2, 1, 2)
    st_ids = [gs["station_id"] for gs in stations]
    y_map = {st_id: idx for idx, st_id in enumerate(st_ids)}

    for st_id, sat_id, aos_s, los_s, status in contacts:
        y_val = y_map[st_id]
        start_hr = aos_s / 3600.0
        width_hr = max((los_s - aos_s) / 3600.0, 0.08)  # Min visual width on 24-hr axis
        if status == "SCHEDULED":
            bar_color = sat_colors.get(sat_id, "#16a34a")
            y_offset = 0.16
            alpha_val = 0.92
        elif status == "REJECTED_DISH_CONFLICT":
            bar_color = "#dc2626"
            y_offset = -0.18
            alpha_val = 0.85
        else:
            bar_color = "#94a3b8"
            y_offset = -0.18
            alpha_val = 0.55

        ax_gantt.broken_barh([(start_hr, width_hr)], (y_val + y_offset - 0.13, 0.26),
                             facecolors=bar_color, alpha=alpha_val)

    ax_gantt.set_yticks(range(len(st_ids)))
    ax_gantt.set_yticklabels(st_ids)
    ax_gantt.set_xlim(0.0, 24.0)
    ax_gantt.set_xticks(range(0, 25, 2))
    ax_gantt.set_xlabel("Mission Elapsed Time Since Epoch 2026-11-17 12:00 UTC (Hours)")
    ax_gantt.set_title("24-Hour Ground Station Contact Timeline (Upper Track: Scheduled | Lower Track: Dropped Dish Conflict / Low Margin)")
    ax_gantt.grid(True, axis="x", linestyle=":", alpha=0.6)

    legend_patches = [
        Patch(facecolor="#2563eb", label="Scheduled: VANGUARD-SAR-1 (P1)"),
        Patch(facecolor="#16a34a", label="Scheduled: VANGUARD-EO-2 (P2)"),
        Patch(facecolor="#d97706", label="Scheduled: SENTINEL-SIG-3 (P1)"),
        Patch(facecolor="#7c3aed", label="Scheduled: PATHFINDER-TD-4 (P3)"),
        Patch(facecolor="#dc2626", label="Dropped: Single-Dish Contention"),
        Patch(facecolor="#94a3b8", label="Dropped: X-Band Margin < 3.0 dB")
    ]
    ax_gantt.legend(handles=legend_patches, loc="upper right", fontsize=8, ncol=3)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(PNG_OUT)
    plt.close(fig)

    # 2. Export 3D Orbital Trajectory & Ground Network KML for Google Earth Pro
    kml_Colors = {
        "VANGUARD-SAR-1": "ffff6325",
        "VANGUARD-EO-2": "ff4aa316",
        "SENTINEL-SIG-3": "ff0677d9",
        "PATHFINDER-TD-4": "ffed3a7c"
    }

    kml_parts = ["""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
<Document>
    <name>LEO Constellation 3D J2 Orbits &amp; Ground Station Network</name>"""]

    for gs in stations:
        kml_parts.append(f"""
    <Placemark>
        <name>{gs['station_id']} ({gs['station_name']})</name>
        <description>Min El: {gs['min_el_deg']} deg | G/T: {gs['g_over_t_dbk']} dB/K</description>
        <Point>
            <altitudeMode>clampToGround</altitudeMode>
            <coordinates>{gs['lon_deg']:.5f},{gs['lat_deg']:.5f},0</coordinates>
        </Point>
    </Placemark>""")

    for sat_id, trk in tracks.items():
        coord_tokens = [
            f"{lon:.5f},{lat:.5f},{alt:.1f}"
            for lon, lat, alt in zip(trk["lon_deg"], trk["lat_deg"], trk["alt_m"])
        ]
        kml_parts.append(f"""
    <Placemark>
        <name>{sat_id} ({trk['mission_role']})</name>
        <Style>
            <LineStyle><color>{kml_Colors.get(sat_id, 'ff00ffff')}</color><width>3</width></LineStyle>
        </Style>
        <LineString>
            <extrude>0</extrude>
            <tessellate>0</tessellate>
            <altitudeMode>absolute</altitudeMode>
            <coordinates>{" ".join(coord_tokens)}</coordinates>
        </LineString>
    </Placemark>""")

    kml_parts.append("\n</Document>\n</kml>\n")
    with open(KML_OUT, mode="w", encoding="utf-8") as f:
        f.write("".join(kml_parts))

    print(f"[+] Exported 3-panel mission analysis PNG : {PNG_OUT}")
    print(f"[+] Exported 3D Google Earth Pro KML      : {KML_OUT}")


if __name__ == "__main__":
    export_visuals_and_kml()