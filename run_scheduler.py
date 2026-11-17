#!/usr/bin/env python3
"""
run_scheduler.py
Single-command orchestrator for the LEO Pass & X-Band Link Scheduler.
"""

from src.tle_and_stations import generate_mission_inputs
from src.orbit_link_propagator import propagate_and_compute_passes
from src.pass_deconfliction_sql import build_and_deconflict_schedule
from src.visualize_and_export import export_visuals_and_kml

if __name__ == "__main__":
    print("======================================================================================")
    print("INITIALIZING LEO CONSTELLATION J2 PASS SCHEDULER & X-BAND LINK DECONFLICTION PIPELINE")
    print("======================================================================================\n")
    generate_mission_inputs()
    raw_passes, stations, specs = propagate_and_compute_passes()
    build_and_deconflict_schedule(raw_passes, stations, specs)
    export_visuals_and_kml()
    print("\n[SUCCESS] LEO Pass & Link Scheduler complete. Artifacts saved to data/processed/.")