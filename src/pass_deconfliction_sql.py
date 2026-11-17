#!/usr/bin/env python3
"""
pass_deconfliction_sql.py
Loads raw orbital visibility passes into an embedded SQLite mission database, detects
single-dish antenna contention (including 60-second slew/re-lock turnaround buffers),
executes priority-weighted greedy deconfliction, and outputs SQL mission analytics.
"""

import os
import csv
import sqlite3
from datetime import datetime, timedelta

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC_DIR = os.path.join(BASE_DIR, "data", "processed")
DB_PATH = os.path.join(PROC_DIR, "space_ops_schedule.db")
CSV_OUT = os.path.join(PROC_DIR, "deconflicted_contact_schedule.csv")
REPORT_OUT = os.path.join(PROC_DIR, "mission_schedule_report.txt")

EPOCH_START = datetime(2026, 11, 17, 12, 0, 0)
MIN_PEAK_MARGIN_DB = 3.0  # Minimum required peak X-band safety margin


def sec_to_utc_str(sec: float) -> str:
    return (EPOCH_START + timedelta(seconds=float(sec))).strftime("%Y-%m-%d %H:%M:%SZ")


def format_sql_table(headers, rows) -> str:
    if not rows:
        return "  (No records)\n"
    str_rows = [[str(cell) for cell in r] for r in rows]
    widths = [max(len(str(h)), max(len(r[i]) for r in str_rows)) for i, h in enumerate(headers)]
    hdr_line = "  ".join(str(h).ljust(widths[i]) for i, h in enumerate(headers))
    sep_line = "  ".join("-" * widths[i] for i in range(len(headers)))
    body = "\n".join("  ".join(r[i].ljust(widths[i]) for i, _ in enumerate(headers)) for r in str_rows)
    return f"{hdr_line}\n{sep_line}\n{body}"


def build_and_deconflict_schedule(raw_passes: list, stations: list, specs: list):
    os.makedirs(PROC_DIR, exist_ok=True)
    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)

    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()

    cur.execute("""
        CREATE TABLE ground_stations (
            station_id TEXT PRIMARY KEY,
            station_name TEXT NOT NULL,
            lat_deg REAL NOT NULL,
            lon_deg REAL NOT NULL,
            min_el_deg REAL NOT NULL,
            g_over_t_dbk REAL NOT NULL,
            turnaround_slew_sec REAL NOT NULL
        );
    """)

    cur.execute("""
        CREATE TABLE spacecraft (
            sat_id TEXT PRIMARY KEY,
            norad_id INTEGER UNIQUE NOT NULL,
            mission_role TEXT NOT NULL,
            priority_tier INTEGER NOT NULL,
            eirp_dbw REAL NOT NULL,
            downlink_rate_mbps REAL NOT NULL
        );
    """)

    cur.execute("""
        CREATE TABLE contact_windows (
            pass_id TEXT PRIMARY KEY,
            sat_id TEXT NOT NULL,
            station_id TEXT NOT NULL,
            priority_tier INTEGER NOT NULL,
            aos_sec REAL NOT NULL,
            los_sec REAL NOT NULL,
            aos_utc TEXT NOT NULL,
            los_utc TEXT NOT NULL,
            duration_sec REAL NOT NULL,
            max_el_deg REAL NOT NULL,
            min_slant_km REAL NOT NULL,
            peak_link_margin_db REAL NOT NULL,
            mean_link_margin_db REAL NOT NULL,
            max_doppler_khz REAL NOT NULL,
            est_data_vol_gb REAL NOT NULL,
            schedule_status TEXT NOT NULL,
            conflict_with_pass_id TEXT,
            FOREIGN KEY(sat_id) REFERENCES spacecraft(sat_id),
            FOREIGN KEY(station_id) REFERENCES ground_stations(station_id)
        );
    """)

    for gs in stations:
        cur.execute(
            "INSERT INTO ground_stations VALUES (?, ?, ?, ?, ?, ?, ?);",
            (gs["station_id"], gs["station_name"], gs["lat_deg"], gs["lon_deg"],
             gs["min_el_deg"], gs["g_over_t_dbk"], gs["turnaround_slew_sec"])
        )

    for sp in specs:
        cur.execute(
            "INSERT INTO spacecraft VALUES (?, ?, ?, ?, ?, ?);",
            (sp["sat_id"], sp["norad_id"], sp["mission_role"],
             sp["priority_tier"], sp["eirp_dbw"], sp["downlink_rate_mbps"])
        )

    for p in raw_passes:
        init_status = "CANDIDATE" if p["peak_link_margin_db"] >= MIN_PEAK_MARGIN_DB else "REJECTED_LOW_MARGIN"
        cur.execute(
            "INSERT INTO contact_windows VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);",
            (
                p["pass_id"], p["sat_id"], p["station_id"], p["priority_tier"],
                p["aos_sec"], p["los_sec"], sec_to_utc_str(p["aos_sec"]), sec_to_utc_str(p["los_sec"]),
                p["duration_sec"], p["max_el_deg"], p["min_slant_km"],
                p["peak_link_margin_db"], p["mean_link_margin_db"],
                p["max_doppler_khz"], p["est_data_vol_gb"], init_status, None
            )
        )

    # Priority + Data-Volume Greedy Deconfliction with Dish Slew Turnaround Buffer
    slew_lookup = {gs["station_id"]: gs["turnaround_slew_sec"] for gs in stations}
    cur.execute("""
        SELECT pass_id, station_id, aos_sec, los_sec
        FROM contact_windows
        WHERE schedule_status = 'CANDIDATE'
        ORDER BY priority_tier ASC, est_data_vol_gb DESC, max_el_deg DESC;
    """)
    candidates = cur.fetchall()

    scheduled_by_station = {gs["station_id"]: [] for gs in stations}

    for pass_id, st_id, aos_s, los_s in candidates:
        slew_buf = slew_lookup[st_id]
        conflict_id = None
        for sched_id, sched_aos, sched_los in scheduled_by_station[st_id]:
            if (aos_s < sched_los + slew_buf) and (los_s + slew_buf > sched_aos):
                conflict_id = sched_id
                break

        if conflict_id is None:
            scheduled_by_station[st_id].append((pass_id, aos_s, los_s))
            cur.execute(
                "UPDATE contact_windows SET schedule_status = 'SCHEDULED' WHERE pass_id = ?;",
                (pass_id,)
            )
        else:
            cur.execute(
                "UPDATE contact_windows SET schedule_status = 'REJECTED_DISH_CONFLICT', "
                "conflict_with_pass_id = ? WHERE pass_id = ?;",
                (conflict_id, pass_id)
            )

    conn.commit()

    # Export Full Deconflicted Contact Deck to CSV
    cur.execute("SELECT * FROM contact_windows ORDER BY aos_sec ASC;")
    col_names = [desc[0] for desc in cur.description]
    all_rows = cur.fetchall()
    with open(CSV_OUT, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(col_names)
        writer.writerows(all_rows)

    # Analytical SQL Summaries
    cur.execute("""
        SELECT
            s.sat_id AS Spacecraft,
            'PRI-' || s.priority_tier AS Pri,
            COUNT(c.pass_id) AS Raw_Passes,
            SUM(CASE WHEN c.schedule_status = 'SCHEDULED' THEN 1 ELSE 0 END) AS Scheduled,
            SUM(CASE WHEN c.schedule_status = 'REJECTED_DISH_CONFLICT' THEN 1 ELSE 0 END) AS Dish_Conflicts,
            SUM(CASE WHEN c.schedule_status = 'REJECTED_LOW_MARGIN' THEN 1 ELSE 0 END) AS Low_RF_Margin,
            ROUND(SUM(CASE WHEN c.schedule_status = 'SCHEDULED' THEN c.est_data_vol_gb ELSE 0 END), 1) AS Downlink_GB
        FROM spacecraft s
        LEFT JOIN contact_windows c ON s.sat_id = c.sat_id
        GROUP BY s.sat_id
        ORDER BY s.priority_tier ASC, Downlink_GB DESC;
    """)
    sat_headers = [d[0] for d in cur.description]
    sat_rows = cur.fetchall()

    cur.execute("""
        SELECT
            g.station_id AS Station,
            COUNT(c.pass_id) AS Visibility_Hits,
            SUM(CASE WHEN c.schedule_status = 'SCHEDULED' THEN 1 ELSE 0 END) AS Contacts_Allocated,
            SUM(CASE WHEN c.schedule_status = 'REJECTED_DISH_CONFLICT' THEN 1 ELSE 0 END) AS Contention_Drops,
            ROUND(SUM(CASE WHEN c.schedule_status = 'SCHEDULED' THEN c.duration_sec ELSE 0 END) / 864.0, 1) || '%' AS Duty_Cycle,
            ROUND(AVG(CASE WHEN c.schedule_status = 'SCHEDULED' THEN c.peak_link_margin_db END), 2) AS Avg_Peak_Margin_dB
        FROM ground_stations g
        LEFT JOIN contact_windows c ON g.station_id = c.station_id
        GROUP BY g.station_id
        ORDER BY Contacts_Allocated DESC;
    """)
    gs_headers = [d[0] for d in cur.description]
    gs_rows = cur.fetchall()

    cur.execute("""
        SELECT
            c.pass_id AS Dropped_Pass,
            c.sat_id AS Dropped_Sat,
            c.station_id AS Station,
            SUBSTR(c.aos_utc, 12, 8) AS AOS_UTC,
            c.conflict_with_pass_id AS Winning_Pass,
            w.sat_id AS Winning_Sat,
            c.max_el_deg AS Max_El
        FROM contact_windows c
        JOIN contact_windows w ON c.conflict_with_pass_id = w.pass_id
        WHERE c.schedule_status = 'REJECTED_DISH_CONFLICT'
        ORDER BY c.aos_sec ASC
        LIMIT 8;
    """)
    conf_headers = [d[0] for d in cur.description]
    conf_rows = cur.fetchall()

    report_str = "\n".join([
        "=" * 86,
        "24-HOUR LEO CONSTELLATION PASS SCHEDULER & X-BAND LINK DECONFLICTION REPORT",
        "=" * 86,
        "\n[1] SPACECRAFT DOWNLINK FULFILLMENT & CONTENTION SUMMARY:",
        format_sql_table(sat_headers, sat_rows),
        "\n[2] GROUND STATION ANTENNA UTILIZATION & RF LINK PERFORMANCE:",
        format_sql_table(gs_headers, gs_rows),
        "\n[3] SAMPLE SINGLE-DISH CONTENTION RESOLUTIONS (60-SEC SLEW BUFFER ENFORCED):",
        format_sql_table(conf_headers, conf_rows),
        "=" * 86
    ])

    print(report_str)
    with open(REPORT_OUT, mode="w", encoding="utf-8") as f:
        f.write(report_str + "\n")

    conn.close()
    print(f"[+] Saved SQLite schedule database : {DB_PATH}")
    print(f"[+] Exported contact deck CSV      : {CSV_OUT}")


if __name__ == "__main__":
    from src.orbit_link_propagator import propagate_and_compute_passes
    passes, g_stations, s_specs = propagate_and_compute_passes()
    build_and_deconflict_schedule(passes, g_stations, s_specs)