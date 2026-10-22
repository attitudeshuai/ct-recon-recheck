#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""判分输出：11 条行为测试，按重要程度加权，权重和 = 1。

每条测试自己生成一份扫描数据（解析体模 → 投影 → 计数），驱动被测程序的命令行入口，
把它写出的报告与体模真值比对。真值由测试按几何解析式独立复算，不读被测代码，
也不与参考实现的输出做比对。
"""
import json
import math
import os
import subprocess
import sys
import tempfile

APP = os.environ.get("CT_APP_DIR", "/app")
ENTRY = os.path.join(APP, "src", "recon.py")

WEIGHTS = {
    "test_roi_mu_matches_phantom_on_base_scan": 0.15,          # 核心：三个采区的 µ
    "test_roi_mu_matches_second_phantom_geometry": 0.12,       # 换几何换材料，防写死
    "test_dark_current_large_offset_recovered": 0.10,          # 暗电流口径
    "test_dead_channel_recovery_with_unlisted_zero": 0.10,     # 失效道（含自检漏报）
    "test_defect_area_and_center_match_phantom": 0.11,         # 缺陷面积与质心
    "test_no_spurious_defect_outside_inspection_radius": 0.07,  # 检视范围边界
    "test_repaired_channels_listed_exactly": 0.07,             # 补道清单如实
    "test_invalid_scan_rejected_without_output": 0.10,         # 非法输入拒收
    "test_report_schema_and_rounding": 0.07,                   # 输出契约
    "test_rerun_produces_identical_bytes": 0.06,               # 重复运行一致
    "test_failed_run_keeps_previous_report": 0.05,             # 失败不污染旧报告
}

# 多能谱：材料卡给能谱与相对衰减比，透射率 T = Σ w·exp(-f·P)，P 是 70 keV 下的线积分
SPECTRUM = [(70.0, 0.60, 1.00), (90.0, 0.30, 0.71), (110.0, 0.10, 0.55)]
REF_ENERGY = 70.0


def _transmission(p):
    return sum(w * math.exp(-f * p) for _k, w, f in SPECTRUM)


BASE_PHANTOM = [
    {"x": 0.0, "y": 0.0, "r": 60.0, "mu": 0.50},
    {"x": -20.0, "y": 15.0, "r": 22.0, "mu": 0.62},
    {"x": 25.0, "y": -10.0, "r": 9.0, "mu": 0.05},
]
SECOND_PHANTOM = [
    {"x": 0.0, "y": 0.0, "r": 55.0, "mu": 0.44},
    {"x": 15.0, "y": 25.0, "r": 15.0, "mu": 0.60},
    {"x": -30.0, "y": -20.0, "r": 7.0, "mu": 0.15},
]
BASE_ROIS = [
    {"name": "matrix", "x_mm": 0.0, "y_mm": -35.0, "r_mm": 8.0, "mu": 0.50, "rel": 0.02},
    {"name": "insert", "x_mm": -20.0, "y_mm": 15.0, "r_mm": 8.0, "mu": 0.62, "rel": 0.02},
    {"name": "void", "x_mm": 25.0, "y_mm": -10.0, "r_mm": 3.0, "mu": 0.05, "abs": 0.005},
]
SECOND_ROIS = [
    {"name": "matrix", "x_mm": 0.0, "y_mm": -30.0, "r_mm": 7.0, "mu": 0.44, "rel": 0.02},
    {"name": "insert", "x_mm": 15.0, "y_mm": 25.0, "r_mm": 6.0, "mu": 0.60, "rel": 0.02},
    {"name": "void", "x_mm": -30.0, "y_mm": -20.0, "r_mm": 2.5, "mu": 0.15, "abs": 0.005},
]


# ---------------------------------------------------------------- 数据生成
def _chord(radius, dist):
    if abs(dist) >= radius:
        return 0.0
    return 2.0 * math.sqrt(radius * radius - dist * dist)


def _line_integral(phantom, s, theta_deg):
    th = math.radians(theta_deg)
    c, sn = math.cos(th), math.sin(th)
    bg = phantom[0]["mu"]
    total = bg * _chord(phantom[0]["r"], s - (phantom[0]["x"] * c + phantom[0]["y"] * sn)) / 10.0
    for disc in phantom[1:]:
        total += (disc["mu"] - bg) * _chord(disc["r"], s - (disc["x"] * c + disc["y"] * sn)) / 10.0
    return total


def make_scan(phantom, rois, size=80, px=2.0, channels=121, pitch=2.0, n_angles=180,
              i0=20000.0, dark=50.0, bad=(), silent=(), drift=0.03, inspect_r=None,
              center_offset=0.0):
    angles = [a * 180.0 / n_angles for a in range(n_angles)]
    s0 = -(channels - 1) / 2.0 * pitch
    per_angle, counts = [], []
    for a in range(n_angles):
        i0a = i0 * (1.0 + drift * math.sin(2 * math.pi * a / 37.0)
                    + 0.5 * drift * math.cos(2 * math.pi * a / 11.0))
        per_angle.append(round(i0a, 6))
        row = []
        for k in range(channels):
            depth = _line_integral(phantom, s0 + k * pitch - center_offset, angles[a])
            cnt = int(round(i0a * _transmission(depth) + dark))
            if k in bad or k in silent:
                cnt = 0
            row.append(cnt)
        counts.append(row)
    return {
        "instrument": {"source_photons": i0, "dark_counts": dark},
        "geometry": {"image_size": size, "pixel_size_mm": px, "channels": channels,
                     "channel_pitch_mm": pitch, "angles_deg": angles},
        "projections": {"per_angle_source": per_angle, "counts": counts},
        "inspection": {
            "radius_mm": inspect_r if inspect_r is not None else phantom[0]["r"] - 2.0,
            "defect_mu_threshold_per_cm": 0.20,
            "min_defect_area_mm2": 30.0,
        },
        "rois": [{"name": r["name"], "x_mm": r["x_mm"], "y_mm": r["y_mm"], "r_mm": r["r_mm"]} for r in rois],
        "material_card": {
            "reference_energy_kev": REF_ENERGY,
            "matrix_mu_per_cm": phantom[0]["mu"],
            "spectrum": [{"kev": kev, "weight": w, "mu_relative": f} for kev, w, f in SPECTRUM],
        },
        "detector": {"bad_channels": sorted(set(bad))},
    }


# ---------------------------------------------------------------- 驱动被测程序
def _run(scan, out_path, timeout=240):
    cmd = [sys.executable, ENTRY, "--scan", scan, "--out", out_path]
    return subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                          errors="replace", timeout=timeout)


def run_scan(scan_obj):
    """写扫描 → 跑程序 → 读回报告；返回 (report 或 None, 上下文)。"""
    tmp = tempfile.mkdtemp(prefix="ct-case-")
    scan_path = os.path.join(tmp, "scan.json")
    out_path = os.path.join(tmp, "report.json")
    with open(scan_path, "w", encoding="utf-8") as f:
        json.dump(scan_obj, f, ensure_ascii=False)
    proc = _run(scan_path, out_path)
    report = None
    if os.path.exists(out_path):
        try:
            with open(out_path, encoding="utf-8") as f:
                report = json.load(f)
        except ValueError:
            report = None
    return report, {"scan": scan_path, "out": out_path, "proc": proc, "dir": tmp}


def _roi_ok(report, rois):
    if not isinstance(report, dict):
        return False
    got = {r.get("name"): r.get("mu_per_cm") for r in report.get("roi") or []}
    for roi in rois:
        if roi["name"] not in got:
            return False
        v = float(got[roi["name"]])
        if not math.isfinite(v):
            return False
        if "rel" in roi:
            if abs(v - roi["mu"]) > roi["rel"] * roi["mu"]:
                return False
        elif abs(v - roi["mu"]) > roi["abs"]:
            return False
    return True


# ---------------------------------------------------------------- 测试
def test_roi_mu_matches_phantom_on_base_scan():
    """题面原文：「roi 按扫描卡里 rois 的顺序输出」「按 rois 里给的采区计算平均 µ」。

    基体与夹杂按相对误差 2% 卡，空洞按绝对误差 0.005 /cm 卡。数据是多能谱的，
    不做线性化时这条路会整体偏低约两成。
    """
    report, _ = run_scan(make_scan(BASE_PHANTOM, BASE_ROIS, bad=(40, 41, 118), silent=(77,), center_offset=-3.4))
    return _roi_ok(report, BASE_ROIS)


def test_roi_mu_matches_second_phantom_geometry():
    """题面原文：「复检时程序会被喂别的批次，几何、材料、道数与角度都和前一批不同」。

    换一套几何（道数、角度覆盖、源强、暗电流都变）与材料后，三个采区仍要落在容差内。
    """
    report, _ = run_scan(make_scan(SECOND_PHANTOM, SECOND_ROIS, channels=121, pitch=2.0,
                                   n_angles=180, i0=15000.0, dark=120.0, bad=(55,), silent=(30,), center_offset=2.6))
    return _roi_ok(report, SECOND_ROIS)


def test_dark_current_large_offset_recovered():
    """题面原文：「计数按 I = I0(θ)·exp(-∫µ dl) + 暗电流 采集」。

    暗电流抬到源强的百分之四后，不扣暗电流的实现采区会明显偏低；采区仍须落在 2% 容差内。
    """
    report, _ = run_scan(make_scan(BASE_PHANTOM, BASE_ROIS, dark=800.0, bad=(40,), silent=(77,), center_offset=-3.4))
    return _roi_ok(report, BASE_ROIS)


def test_dead_channel_recovery_with_unlisted_zero():
    """题面原文：「自检也可能漏报，计数为 0 的道一律按失效道处理」。

    自检清单为空而存在 4 个零计数道时，采区仍要合格。
    """
    report, _ = run_scan(make_scan(BASE_PHANTOM, BASE_ROIS, bad=(), silent=(20, 21, 22, 90), center_offset=-3.4))
    return _roi_ok(report, BASE_ROIS)


def test_defect_area_and_center_match_phantom():
    """题面原文：「defects 是检视范围内 µ 低于阈值的连通块」「面积按 mm² 计」「按面积从大到小排列」。

    体模里只有一个 9 mm 半径的空洞：恰好检出 1 个，质心偏差不超过 2 mm，面积在真值 ±12% 内。
    """
    report, _ = run_scan(make_scan(BASE_PHANTOM, BASE_ROIS, bad=(40, 41, 118), silent=(77,), center_offset=-3.4))
    if not isinstance(report, dict):
        return False
    defects = report.get("defects")
    if not isinstance(defects, list) or len(defects) != 1:
        return False
    truth_area = math.pi * 9.0 * 9.0
    try:
        cx, cy, area = float(defects[0]["x_mm"]), float(defects[0]["y_mm"]), float(defects[0]["area_mm2"])
    except (KeyError, TypeError, ValueError):
        return False
    return (math.hypot(cx - 25.0, cy + 10.0) <= 2.0
            and abs(area - truth_area) <= 0.12 * truth_area)


def test_no_spurious_defect_outside_inspection_radius():
    """题面原文：「在 inspection.radius_mm 检视范围内找内部缺陷」。

    检视半径收到 52 mm（样品半径 60 mm）时，样品外的空气不得被报成缺陷，报告里只应有那个空洞。
    """
    report, _ = run_scan(make_scan(BASE_PHANTOM, BASE_ROIS, bad=(40, 41, 118), silent=(77,),
                                   inspect_r=52.0, center_offset=-3.4))
    if not isinstance(report, dict):
        return False
    defects = report.get("defects")
    if not isinstance(defects, list) or len(defects) != 1:
        return False
    try:
        area = float(defects[0]["area_mm2"])
    except (KeyError, TypeError, ValueError):
        return False
    return 150.0 <= area <= 350.0


def test_repaired_channels_listed_exactly():
    """题面原文：「repaired_channels 是本次判为失效并补过的道号，升序去重」。

    自检列了 12 与 60，另有一道 90 只有零计数：三者都要列出，且升序去重。
    """
    report, _ = run_scan(make_scan(BASE_PHANTOM, BASE_ROIS, bad=(12, 60), silent=(90,), center_offset=-3.4))
    if not isinstance(report, dict):
        return False
    got = (report.get("meta") or {}).get("repaired_channels")
    return got == [12, 60, 90]


def test_invalid_scan_rejected_without_output():
    """题面原文：「扫描文件不合法时以退出码 2 结束，往 stderr 写一句能看懂的原因，并且不产出报告」。

    四种坏输入逐个验证：退出码 2、无报告文件、无遗留临时文件、stderr 非空。
    """
    good = make_scan(BASE_PHANTOM, BASE_ROIS, bad=(40,), silent=(), center_offset=-3.4)
    missing = json.loads(json.dumps(good))
    missing.pop("projections")
    wrong_ch = json.loads(json.dumps(good))
    wrong_ch["geometry"]["channels"] = 100
    bad_instrument = json.loads(json.dumps(good))
    bad_instrument["instrument"]["dark_counts"] = bad_instrument["instrument"]["source_photons"]
    bad_geom = json.loads(json.dumps(good))
    bad_geom["geometry"]["pixel_size_mm"] = -1.0
    for scan in (missing, wrong_ch, bad_instrument, bad_geom):
        tmp = tempfile.mkdtemp(prefix="ct-bad-")
        scan_path = os.path.join(tmp, "scan.json")
        out_path = os.path.join(tmp, "report.json")
        with open(scan_path, "w", encoding="utf-8") as f:
            json.dump(scan, f, ensure_ascii=False)
        proc = _run(scan_path, out_path)
        if proc.returncode != 2 or os.path.exists(out_path):
            return False
        if not (proc.stderr or "").strip():
            return False
        if [n for n in os.listdir(tmp) if n.endswith(".tmp")]:
            return False
    return True


def test_report_schema_and_rounding():
    """题面原文：「报告是一个 UTF-8 编码的 JSON 对象，结构固定」「mu_per_cm 保留三位小数」。

    校验顶层键、采区顺序与名称、meta 计数、数值位数，以及文件里不出现 NaN 或 Infinity。
    """
    scan_obj = make_scan(BASE_PHANTOM, BASE_ROIS, bad=(40, 41, 118), silent=(77,))
    report, ctx = run_scan(scan_obj)
    if not isinstance(report, dict):
        return False
    if set(report) != {"geometry", "roi", "defects", "meta"}:
        return False
    geo = report.get("geometry") or {}
    if geo.get("image_size") != 80 or float(geo.get("pixel_size_mm", 0)) != 2.0:
        return False
    if [r.get("name") for r in report.get("roi") or []] != [r["name"] for r in BASE_ROIS]:
        return False
    for r in report["roi"]:
        mu = float(r["mu_per_cm"])
        if not math.isfinite(mu) or abs(mu - round(mu, 3)) > 1e-12:
            return False
        if not isinstance(r.get("pixels"), int) or r["pixels"] <= 0:
            return False
    meta = report.get("meta") or {}
    if meta.get("n_angles") != 180 or meta.get("n_channels") != 121:
        return False
    raw = open(ctx["out"], encoding="utf-8").read()
    return "NaN" not in raw and "Infinity" not in raw


def test_rerun_produces_identical_bytes():
    """题面原文：「同一份扫描重复运行，报告逐字节一致」。

    连跑两次，两个报告文件的字节完全相同。
    """
    scan_obj = make_scan(BASE_PHANTOM, BASE_ROIS, bad=(40, 41, 118), silent=(77,))
    first, ctx1 = run_scan(scan_obj)
    if first is None:
        return False
    body1 = open(ctx1["out"], "rb").read()
    second, ctx2 = run_scan(scan_obj)
    if second is None:
        return False
    return body1 == open(ctx2["out"], "rb").read()


def test_failed_run_keeps_previous_report():
    """题面原文：「运行失败时不能覆盖上一次已经写好的报告」。

    先跑一次合格扫描留下报告，再用道数不匹配的坏输入跑一次：退出码 2，旧报告字节不变。
    """
    good = make_scan(BASE_PHANTOM, BASE_ROIS, bad=(40,), silent=(), center_offset=-3.4)
    tmp = tempfile.mkdtemp(prefix="ct-keep-")
    scan_path = os.path.join(tmp, "scan.json")
    out_path = os.path.join(tmp, "report.json")
    with open(scan_path, "w", encoding="utf-8") as f:
        json.dump(good, f, ensure_ascii=False)
    first = _run(scan_path, out_path)
    if first.returncode != 0 or not os.path.exists(out_path):
        return False
    before = open(out_path, "rb").read()
    broken = json.loads(json.dumps(good))
    broken["geometry"]["channels"] = 100
    with open(scan_path, "w", encoding="utf-8") as f:
        json.dump(broken, f, ensure_ascii=False)
    second = _run(scan_path, out_path)
    if second.returncode != 2:
        return False
    return open(out_path, "rb").read() == before


def main():
    total = 0.0
    for name, weight in WEIGHTS.items():
        fn = globals().get(name)
        try:
            passed = bool(fn and fn())
        except Exception as e:  # noqa: BLE001 —— 任何异常都按不通过计
            print("%s: 异常 %r" % (name, e), file=sys.stderr)
            passed = False
        total += weight * passed
        print("%s: %s (权重 %.2f)" % (name, "PASS" if passed else "FAIL", weight), file=sys.stderr)
    print("%.6f" % total)


if __name__ == "__main__":
    main()
