#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""生成 CT 扫描数据（解析体模 → 投影 → 计数）。同时用于测试变体。"""
import json
import math

# 多能谱前向：透射率 T = Σ w_k·exp(-f_k·P)，P 是 70 keV 下的线积分（µ70·长度）
# f_k 是同一材料在 E_k 下相对 70 keV 的衰减比；材料卡里带能谱与 f，各材料共用同一条能量依赖
SPECTRUM = [(70.0, 0.60, 1.00), (90.0, 0.30, 0.71), (110.0, 0.10, 0.55)]
REF_ENERGY = 70.0


def transmission(p):
    return sum(w * math.exp(-f * p) for _kev, w, f in SPECTRUM)


DEFAULT_PHANTOM = [
    {"x": 0.0, "y": 0.0, "r": 60.0, "mu": 0.50},
    {"x": -20.0, "y": 15.0, "r": 22.0, "mu": 0.62},
    {"x": 25.0, "y": -10.0, "r": 9.0, "mu": 0.05},
]


def chord(radius, dist):
    if abs(dist) >= radius:
        return 0.0
    return 2.0 * math.sqrt(radius * radius - dist * dist)


def line_integral(phantom, s, theta_deg):
    """µ 单位 1/cm，长度 mm：Σ µ[1/cm]·弦长[mm]/10 = 无量纲线积分。"""
    th = math.radians(theta_deg)
    c, sn = math.cos(th), math.sin(th)
    bg = phantom[0]["mu"]
    total = bg * chord(phantom[0]["r"], s - (phantom[0]["x"] * c + phantom[0]["y"] * sn)) / 10.0
    for disc in phantom[1:]:
        dc = disc["x"] * c + disc["y"] * sn
        total += (disc["mu"] - bg) * chord(disc["r"], s - dc) / 10.0
    return total


def source_drift(index, amplitude=0.03, secondary=0.015):
    return 1.0 + amplitude * math.sin(2 * math.pi * index / 37.0) + secondary * math.cos(2 * math.pi * index / 11.0)


def make_scan(phantom=None, size=80, px=2.0, channels=121, pitch=2.0, n_angles=180,
              i0=20000.0, dark=50.0, bad_channels=(40, 41, 118), silent_dead=(77,),
              drift_amp=0.03, drift_2=0.015, rois=None, inspect_r=None, thr=0.20,
              min_area=30.0, matrix_mu=None, center_offset_mm=0.0):
    phantom = phantom or DEFAULT_PHANTOM
    if inspect_r is None:
        inspect_r = phantom[0]["r"] - 2.0
    angles = [a * 180.0 / n_angles for a in range(n_angles)]
    s0 = -(channels - 1) / 2.0 * pitch
    per_angle, counts = [], []
    for a in range(n_angles):
        i0a = i0 * source_drift(a, drift_amp, drift_2)
        per_angle.append(round(i0a, 6))
        row = []
        for k in range(channels):
            s = s0 + k * pitch
            p = line_integral(phantom, s - center_offset_mm, angles[a])
            cnt = int(round(i0a * transmission(p) + dark))
            if k in bad_channels or k in silent_dead:
                cnt = 0
            row.append(cnt)
        counts.append(row)
    rois = rois or [
        {"name": "matrix", "x_mm": 0.0, "y_mm": -35.0, "r_mm": 8.0},
        {"name": "insert", "x_mm": -20.0, "y_mm": 15.0, "r_mm": 8.0},
        {"name": "void", "x_mm": 25.0, "y_mm": -10.0, "r_mm": 3.0},
    ]
    return {
        "instrument": {"source_photons": i0, "dark_counts": dark},
        "geometry": {
            "image_size": size, "pixel_size_mm": px, "channels": channels,
            "channel_pitch_mm": pitch, "angles_deg": angles,
        },
        "projections": {"per_angle_source": per_angle, "counts": counts},
        "inspection": {
            "radius_mm": inspect_r,
            "defect_mu_threshold_per_cm": thr,
            "min_defect_area_mm2": min_area,
        },
        "rois": rois,
        "material_card": {
            "reference_energy_kev": REF_ENERGY,
            "matrix_mu_per_cm": matrix_mu if matrix_mu is not None else phantom[0]["mu"],
            "spectrum": [{"kev": kev, "weight": w, "mu_relative": f} for kev, w, f in SPECTRUM],
        },
        "detector": {"bad_channels": sorted(set(bad_channels))},
    }


if __name__ == "__main__":
    import os
    import sys

    out = sys.argv[1] if len(sys.argv) > 1 else "scan.json"
    data = make_scan()
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        json.dump(data, f, ensure_ascii=False, separators=(",", ":"))
        f.write("\n")
    print("wrote %s (%d bytes)" % (out, os.path.getsize(out)))
