#!/usr/bin/env python3
"""
tle_and_stations.py
Generates valid 69-character NORAD Two-Line Element (TLE) sets with verified Modulo-10
checksums for a 4-spacecraft LEO constellation (including a close tandem SAR/EO
tip-and-cue pair that induces realistic ground-station dish contention) and a 4-site
global tracking network.
"""

import os
import json

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = os.path.join(BASE_DIR, "data", "raw")


def compute_tle_checksum(line_68: str) -> int:
    """Computes standard NORAD Modulo-10 checksum over the first 68 TLE columns."""
    if len(line_68) != 68:
        raise ValueError(f"Expected 68-char TLE body, got {len(line_68)} chars.")
    total = 0
    for ch in line_68:
        if ch.isdigit():
            total += int(ch)
        elif ch == "-":
            total += 1
    return total % 10


def format_tle_pair(
    norad_id: int,
    intl_desig: str,
    epoch_yr: int,
    epoch_day: float,
    inc_deg: float,
    raan_deg: float,
    ecc: float,
    argp_deg: float,
    ma_deg: float,
    mean_motion_rev_day: float,
    rev_num: int = 1840
):
    """Formats exact 69-character NORAD TLE Line 1 and Line 2 strings."""
    l1_body = (
        f"1 {norad_id:05d}U {intl_desig:<8s} "
        f"{epoch_yr:02d}{epoch_day:012.8f}  .00001200  00000-0  52000-4 0  999"
    )
    ecc_str = f"{int(round(ecc * 1.0e7)):07d}"
    l2_body = (
        f"2 {norad_id:05d} {inc_deg:8.4f} {raan_deg:8.4f} {ecc_str} "
        f"{argp_deg:8.4f} {ma_deg:8.4f} {mean_motion_rev_day:11.8f}{rev_num:5d}"
    )
    line1 = f"{l1_body}{compute_tle_checksum(l1_body)}"
    line2 = f"{l2_body}{compute_tle_checksum(l2_body)}"
    return line1, line2


def generate_mission_inputs():
    os.makedirs(RAW_DIR, exist_ok=True)

    # 4-Satellite Unclassified LEO Constellation
    # Note: VANGUARD-SAR-1 and VANGUARD-EO-2 fly a co-planar/tandem tip-and-cue formation
    # (~2.2 deg separation), intentionally triggering single-dish contention at high latitudes.
    constellation = [
        {
            "sat_id": "VANGUARD-SAR-1",
            "norad_id": 61001,
            "intl_desig": "26041A",
            "mission_role": "All-Weather X-Band SAR",
            "priority_tier": 1,
            "eirp_dbw": 16.5,
            "downlink_rate_mbps": 150.0,
            "inc_deg": 97.5200,
            "raan_deg": 48.0000,
            "ecc": 0.0011500,
            "argp_deg": 92.0000,
            "ma_deg": 105.0000,
            "mean_motion": 15.12540000  # ~530 km Sun-Synchronous Orbit
        },
        {
            "sat_id": "VANGUARD-EO-2",
            "norad_id": 61002,
            "intl_desig": "26041B",
            "mission_role": "High-Res Optical Tip-and-Cue",
            "priority_tier": 2,
            "eirp_dbw": 15.5,
            "downlink_rate_mbps": 150.0,
            "inc_deg": 97.5200,
            "raan_deg": 49.2000,
            "ecc": 0.0011800,
            "argp_deg": 92.0000,
            "ma_deg": 107.4000,
            "mean_motion": 15.12540000  # Trailing ~2.4 deg behind SAR-1
        },
        {
            "sat_id": "SENTINEL-SIG-3",
            "norad_id": 61003,
            "intl_desig": "26055A",
            "mission_role": "Mid-Inclination Maritime RF",
            "priority_tier": 1,
            "eirp_dbw": 14.8,
            "downlink_rate_mbps": 150.0,
            "inc_deg": 51.6400,
            "raan_deg": 135.5000,
            "ecc": 0.0008500,
            "argp_deg": 45.0000,
            "ma_deg": 210.0000,
            "mean_motion": 15.21980000  # ~500 km Mid-Inclination Orbit
        },
        {
            "sat_id": "PATHFINDER-TD-4",
            "norad_id": 61004,
            "intl_desig": "26062C",
            "mission_role": "Polar Tech-Demo CubeSat",
            "priority_tier": 3,
            "eirp_dbw": 11.2,
            "downlink_rate_mbps": 150.0,
            "inc_deg": 88.2000,
            "raan_deg": 48.6000,
            "ecc": 0.0019000,
            "argp_deg": 180.0000,
            "ma_deg": 106.0000,
            "mean_motion": 14.91850000  # ~595 km Polar Orbit (marginal EIRP at low el)
        }
    ]

    tle_lines = []
    sat_specs = []
    for sat in constellation:
        l1, l2 = format_tle_pair(
            norad_id=sat["norad_id"],
            intl_desig=sat["intl_desig"],
            epoch_yr=26,
            epoch_day=321.50000000,  # 2026-11-17 12:00:00 UTC
            inc_deg=sat["inc_deg"],
            raan_deg=sat["raan_deg"],
            ecc=sat["ecc"],
            argp_deg=sat["argp_deg"],
            ma_deg=sat["ma_deg"],
            mean_motion_rev_day=sat["mean_motion"]
        )
        tle_lines.extend([sat["sat_id"], l1, l2])
        sat_specs.append({
            "sat_id": sat["sat_id"],
            "norad_id": sat["norad_id"],
            "mission_role": sat["mission_role"],
            "priority_tier": sat["priority_tier"],
            "eirp_dbw": sat["eirp_dbw"],
            "downlink_rate_mbps": sat["downlink_rate_mbps"]
        })

    tle_path = os.path.join(RAW_DIR, "constellation.tle")
    with open(tle_path, mode="w", encoding="utf-8") as f:
        f.write("\n".join(tle_lines) + "\n")

    specs_path = os.path.join(RAW_DIR, "spacecraft_specs.json")
    with open(specs_path, mode="w", encoding="utf-8") as f:
        json.dump(sat_specs, f, indent=4)

    ground_stations = [
        {
            "station_id": "GS-SVALBARD",
            "station_name": "Svalbard Polar Station (Norway)",
            "lat_deg": 78.2298,
            "lon_deg": 15.4078,
            "alt_m": 460.0,
            "min_el_deg": 5.0,
            "g_over_t_dbk": 29.5,
            "turnaround_slew_sec": 60.0
        },
        {
            "station_id": "GS-WALLOPS",
            "station_name": "Wallops Flight Facility (Virginia, US)",
            "lat_deg": 37.9402,
            "lon_deg": -75.4664,
            "alt_m": 12.0,
            "min_el_deg": 7.0,
            "g_over_t_dbk": 28.5,
            "turnaround_slew_sec": 60.0
        },
        {
            "station_id": "GS-MAHIA",
            "station_name": "Mahia / Auckland Tracking Site (NZ)",
            "lat_deg": -39.2608,
            "lon_deg": 177.8659,
            "alt_m": 115.0,
            "min_el_deg": 6.0,
            "g_over_t_dbk": 27.8,
            "turnaround_slew_sec": 60.0
        },
        {
            "station_id": "GS-GUAM",
            "station_name": "Andersen / Indo-Pacific Node (Guam)",
            "lat_deg": 13.5840,
            "lon_deg": 144.9245,
            "alt_m": 162.0,
            "min_el_deg": 8.0,
            "g_over_t_dbk": 26.2,
            "turnaround_slew_sec": 60.0
        }
    ]

    gs_path = os.path.join(RAW_DIR, "ground_stations.json")
    with open(gs_path, mode="w", encoding="utf-8") as f:
        json.dump(ground_stations, f, indent=4)

    print(f"[+] Generated verified NORAD 3LE file : {tle_path} ({len(constellation)} spacecraft)")
    print(f"[+] Generated ground station network  : {gs_path} ({len(ground_stations)} sites)")


if __name__ == "__main__":
    generate_mission_inputs()