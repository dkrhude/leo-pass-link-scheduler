#!/usr/bin/env python3
"""
orbit_link_propagator.py
Validates NORAD TLE checksums, propagates 24-hour LEO trajectories under secular WGS84
J2 oblateness perturbations, transforms ECI -> ECEF -> Topocentric SEZ horizon frames,
and computes 8.2 GHz X-Band downlink margins and Doppler shift profiles.
"""

import os
import json
import math
import numpy as np
from src.tle_and_stations import compute_tle_checksum

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = os.path.join(BASE_DIR, "data", "raw")
PROC_DIR = os.path.join(BASE_DIR, "data", "processed")

# WGS84 & Physical Astrodynamics Constants
MU_EARTH_KM3_S2 = 398600.4418
R_EARTH_KM = 6378.137
ECC_SQ_WGS84 = 0.00669437999014
J2_EARTH = 1.08262668e-3
OMEGA_EARTH_RAD_S = 7.2921150e-5
C_LIGHT_MPS = 299792458.0
BOLTZMANN_DBW_K_HZ = -228.6
XBAND_FREQ_HZ = 8.2e9
REQ_EB_N0_DB = 10.5  # QPSK/8PSK LDPC threshold + implementation loss


def parse_tle_file(tle_path: str):
    """Parses 3-line NORAD TLE sets and verifies Modulo-10 checksum integrity."""
    with open(tle_path, mode="r", encoding="utf-8") as f:
        raw_lines = [ln.rstrip("\r\n") for ln in f if ln.strip()]

    satellites = []
    for idx in range(0, len(raw_lines), 3):
        name = raw_lines[idx].strip()
        l1 = raw_lines[idx + 1]
        l2 = raw_lines[idx + 2]

        for line_label, ln in (("Line 1", l1), ("Line 2", l2)):
            expected = compute_tle_checksum(ln[:68])
            actual = int(ln[68])
            if expected != actual:
                raise ValueError(f"TLE checksum mismatch on {name} {line_label}: {actual} != {expected}")

        satellites.append({
            "sat_id": name,
            "norad_id": int(l1[2:7]),
            "epoch_yr": 2000 + int(l1[18:20]),
            "epoch_day": float(l1[20:32]),
            "inc_deg": float(l2[8:16]),
            "raan_deg": float(l2[17:25]),
            "ecc": float("0." + l2[26:33].strip()),
            "argp_deg": float(l2[34:42]),
            "ma_deg": float(l2[43:51]),
            "mean_motion_rev_day": float(l2[52:63])
        })
    return satellites


def gmst_at_epoch_rad(year: int, day_of_year: float) -> float:
    """Computes Greenwich Mean Sidereal Time (radians) at TLE epoch."""
    # Julian Date at Jan 0.0 of given year + day_of_year
    a = (14 - 1) // 12
    y = year + 4800 - a
    m = 1 + 12 * a - 3
    jdn_jan1 = 1 + ((153 * m + 2) // 5) + 365 * y + (y // 4) - (y // 100) + (y // 400) - 32045
    jd = (jdn_jan1 - 0.5) + (day_of_year - 1.0)
    t_ut1 = (jd - 2451545.0) / 36525.0
    gmst_deg = (
        280.46061837
        + 360.98564736629 * (jd - 2451545.0)
        + 0.000387933 * (t_ut1 ** 2)
    ) % 360.0
    return math.radians(gmst_deg)


def propagate_j2_ecef(tle_rec: dict, t_sec: np.ndarray):
    """Vectorized Keplerian + secular J2 propagator returning ECEF (km) and sub-sat coords."""
    n_rad_s = tle_rec["mean_motion_rev_day"] * (2.0 * math.pi / 86400.0)
    a_km = (MU_EARTH_KM3_S2 / (n_rad_s ** 2)) ** (1.0 / 3.0)
    ecc = tle_rec["ecc"]
    inc_rad = math.radians(tle_rec["inc_deg"])
    p_km = a_km * (1.0 - ecc ** 2)

    j2_term = 1.5 * n_rad_s * J2_EARTH * ((R_EARTH_KM / p_km) ** 2)
    raan_dot = -j2_term * math.cos(inc_rad)
    argp_dot = 0.5 * j2_term * (5.0 * (math.cos(inc_rad) ** 2) - 1.0)
    ma_dot = n_rad_s + 0.5 * j2_term * math.sqrt(1.0 - ecc ** 2) * (3.0 * (math.cos(inc_rad) ** 2) - 1.0)

    raan = math.radians(tle_rec["raan_deg"]) + raan_dot * t_sec
    argp = math.radians(tle_rec["argp_deg"]) + argp_dot * t_sec
    mean_anom = np.mod(math.radians(tle_rec["ma_deg"]) + ma_dot * t_sec, 2.0 * np.pi)

    # Vectorized Newton-Raphson Kepler Equation Solver (M = E - e*sin(E))
    ecc_anom = mean_anom.copy()
    for _ in range(8):
        ecc_anom -= (ecc_anom - ecc * np.sin(ecc_anom) - mean_anom) / (1.0 - ecc * np.cos(ecc_anom))

    true_anom = 2.0 * np.arctan2(
         math.sqrt(1.0 + ecc) * np.sin(0.5 * ecc_anom),
         math.sqrt(1.0 - ecc) * np.cos(0.5 * ecc_anom)
    )
    r_km = a_km * (1.0 - ecc * np.cos(ecc_anom))
    u = argp + true_anom

    x_eci = r_km * (np.cos(raan) * np.cos(u) - np.sin(raan) * np.sin(u) * math.cos(inc_rad))
    y_eci = r_km * (np.sin(raan) * np.cos(u) + np.cos(raan) * np.sin(u) * math.cos(inc_rad))
    z_eci = r_km * (np.sin(u) * math.sin(inc_rad))

    gmst0 = gmst_at_epoch_rad(tle_rec["epoch_yr"], tle_rec["epoch_day"])
    gmst = gmst0 + OMEGA_EARTH_RAD_S * t_sec

    x_ecef = np.cos(gmst) * x_eci + np.sin(gmst) * y_eci
    y_ecef = -np.sin(gmst) * x_eci + np.cos(gmst) * y_eci
    z_ecef = z_eci

    sub_lat_deg = np.degrees(np.arcsin(np.clip(z_ecef / r_km, -1.0, 1.0)))
    sub_lon_deg = np.degrees(np.arctan2(y_ecef, x_ecef))
    alt_km = r_km - R_EARTH_KM

    return x_ecef, y_ecef, z_ecef, sub_lat_deg, sub_lon_deg, alt_km


def station_geodetic_to_ecef_km(lat_deg: float, lon_deg: float, alt_m: float):
    """Converts WGS84 geodetic ground-station coordinates to ECEF (km)."""
    lat_r = math.radians(lat_deg)
    lon_r = math.radians(lon_deg)
    alt_km = alt_m / 1000.0
    n_prime = R_EARTH_KM / math.sqrt(1.0 - ECC_SQ_WGS84 * (math.sin(lat_r) ** 2))
    xs = (n_prime + alt_km) * math.cos(lat_r) * math.cos(lon_r)
    ys = (n_prime + alt_km) * math.cos(lat_r) * math.sin(lon_r)
    zs = (n_prime * (1.0 - ECC_SQ_WGS84) + alt_km) * math.sin(lat_r)
    return xs, ys, zs


def propagate_and_compute_passes():
    os.makedirs(PROC_DIR, exist_ok=True)
    tles = parse_tle_file(os.path.join(RAW_DIR, "constellation.tle"))
    with open(os.path.join(RAW_DIR, "spacecraft_specs.json"), "r", encoding="utf-8") as f:
        specs_map = {s["sat_id"]: s for s in json.load(f)}
    with open(os.path.join(RAW_DIR, "ground_stations.json"), "r", encoding="utf-8") as f:
        stations = json.load(f)

    dt = 10.0
    t_sec = np.arange(0.0, 86400.0 + dt, dt)  # 24-hour scheduling horizon

    raw_passes = []
    orbit_tracks_export = {}
    best_pass_profile = {"max_el_deg": -90.0}
    pass_counter = 1

    for tle in tles:
        sat_id = tle["sat_id"]
        spec = specs_map[sat_id]
        x_e, y_e, z_e, sub_lat, sub_lon, alt_km = propagate_j2_ecef(tle, t_sec)

        # Save 100-minute (~1 full orbit, 600 steps) ground/3D track for plotting & KML
        orbit_tracks_export[sat_id] = {
            "priority_tier": spec["priority_tier"],
            "mission_role": spec["mission_role"],
            "lat_deg": np.round(sub_lat[:600], 4).tolist(),
            "lon_deg": np.round(sub_lon[:600], 4).tolist(),
            "alt_m": np.round(alt_km[:600] * 1000.0, 1).tolist()
        }

        for gs in stations:
            xs, ys, zs = station_geodetic_to_ecef_km(gs["lat_deg"], gs["lon_deg"], gs["alt_m"])
            lat_r = math.radians(gs["lat_deg"])
            lon_r = math.radians(gs["lon_deg"])

            dx = x_e - xs
            dy = y_e - ys
            dz = z_e - zs

            south = math.sin(lat_r) * math.cos(lon_r) * dx + math.sin(lat_r) * math.sin(lon_r) * dy - math.cos(lat_r) * dz
            east = -math.sin(lon_r) * dx + math.cos(lon_r) * dy
            zenith = math.cos(lat_r) * math.cos(lon_r) * dx + math.cos(lat_r) * math.sin(lon_r) * dy + math.sin(lat_r) * dz

            slant_km = np.sqrt(south ** 2 + east ** 2 + zenith ** 2)
            el_deg = np.degrees(np.arcsin(np.clip(zenith / slant_km, -1.0, 1.0)))
            range_rate_mps = np.gradient(slant_km * 1000.0, dt)
            doppler_khz = -(range_rate_mps / C_LIGHT_MPS) * (XBAND_FREQ_HZ / 1000.0)

            # X-Band (8.2 GHz) Downlink RF Link Budget
            fspl_db = 20.0 * np.log10((4.0 * math.pi * slant_km * 1000.0 * XBAND_FREQ_HZ) / C_LIGHT_MPS)
            sin_el = np.maximum(np.sin(np.radians(np.maximum(el_deg, 2.0))), 0.05)
            atm_loss_db = (0.38 / sin_el) + 0.50
            c_n0_dbhz = spec["eirp_dbw"] - fspl_db - atm_loss_db + gs["g_over_t_dbk"] - BOLTZMANN_DBW_K_HZ
            eb_n0_db = c_n0_dbhz - 10.0 * math.log10(spec["downlink_rate_mbps"] * 1.0e6)
            link_margin_db = eb_n0_db - REQ_EB_N0_DB

            in_view = el_deg >= gs["min_el_deg"]
            padded = np.concatenate(([False], in_view, [False])).astype(int)
            diffs = np.diff(padded)
            starts = np.where(diffs == 1)[0]
            ends = np.where(diffs == -1)[0] - 1

            for s_idx, e_idx in zip(starts, ends):
                dur_sec = float((e_idx - s_idx) * dt)
                if dur_sec < 60.0:
                    continue

                sl = slice(s_idx, e_idx + 1)
                peak_el = float(np.max(el_deg[sl]))
                peak_margin = float(np.max(link_margin_db[sl]))
                mean_margin = float(np.mean(link_margin_db[sl]))
                min_slant = float(np.min(slant_km[sl]))
                max_dopp = float(np.max(np.abs(doppler_khz[sl])))

                # Usable downlink seconds where link margin >= 0 dB
                closed_steps = int(np.sum(link_margin_db[sl] >= 0.0))
                data_vol_gb = (closed_steps * dt * spec["downlink_rate_mbps"]) / 8000.0

                pass_id = f"PASS-{pass_counter:04d}"
                pass_counter += 1

                raw_passes.append({
                    "pass_id": pass_id,
                    "sat_id": sat_id,
                    "station_id": gs["station_id"],
                    "priority_tier": spec["priority_tier"],
                    "aos_sec": float(t_sec[s_idx]),
                    "los_sec": float(t_sec[e_idx]),
                    "duration_sec": dur_sec,
                    "max_el_deg": round(peak_el, 2),
                    "min_slant_km": round(min_slant, 1),
                    "peak_link_margin_db": round(peak_margin, 2),
                    "mean_link_margin_db": round(mean_margin, 2),
                    "max_doppler_khz": round(max_dopp, 1),
                    "est_data_vol_gb": round(data_vol_gb, 3)
                })

                if peak_el > best_pass_profile["max_el_deg"] and spec["priority_tier"] == 1:
                    best_pass_profile = {
                        "pass_id": pass_id,
                        "sat_id": sat_id,
                        "station_id": gs["station_id"],
                        "max_el_deg": round(peak_el, 2),
                        "rel_time_min": np.round((t_sec[sl] - t_sec[s_idx]) / 60.0, 2).tolist(),
                        "el_deg": np.round(el_deg[sl], 2).tolist(),
                        "link_margin_db": np.round(link_margin_db[sl], 2).tolist(),
                        "doppler_khz": np.round(doppler_khz[sl], 2).tolist()
                    }

    with open(os.path.join(PROC_DIR, "orbital_tracks.json"), "w", encoding="utf-8") as f:
        json.dump(orbit_tracks_export, f)
    with open(os.path.join(PROC_DIR, "best_pass_profile.json"), "w", encoding="utf-8") as f:
        json.dump(best_pass_profile, f)

    print(f"[+] Propagated 24-hr J2 orbits & detected {len(raw_passes)} geometric station passes.")
    return raw_passes, stations, list(specs_map.values())


if __name__ == "__main__":
    propagate_and_compute_passes()