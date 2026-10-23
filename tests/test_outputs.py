# -*- coding: utf-8 -*-
"""sxw-002 判分测试（pytest 版）：每条测试自己造扫描数据，断言写在测试函数体内。

期望值来自解析体模的真值（每条扫描由测试自己按衰减定律正向生成），
不读被测代码，也不与参考解的输出比对。权重见 conftest.py 的 WEIGHTS。
"""
import json
import math
import os
import subprocess
import sys
import tempfile

APP = os.environ.get("CT_APP_DIR", "/app")
ENTRY = os.path.join(APP, "src", "recon.py")

SPECTRUM = [(70.0, 0.60, 1.00), (90.0, 0.30, 0.71), (110.0, 0.10, 0.55)]
REF_ENERGY = 70.0

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
    {"name": "matrix", "x_mm": 0.0, "y_mm": -35.0, "r_mm": 8.0, "mu": 0.50, "tol": 0.02},
    {"name": "insert", "x_mm": -20.0, "y_mm": 15.0, "r_mm": 8.0, "mu": 0.62, "tol": 0.02},
    {"name": "void", "x_mm": 25.0, "y_mm": -10.0, "r_mm": 3.0, "mu": 0.05, "abs": 0.012},
]
SECOND_ROIS = [
    {"name": "matrix", "x_mm": 0.0, "y_mm": -30.0, "r_mm": 7.0, "mu": 0.44, "tol": 0.02},
    {"name": "insert", "x_mm": 15.0, "y_mm": 25.0, "r_mm": 6.0, "mu": 0.60, "tol": 0.02},
    {"name": "void", "x_mm": -30.0, "y_mm": -20.0, "r_mm": 2.5, "mu": 0.15, "abs": 0.012},
]


# ---------------------------------------------------------------- 造扫描（独立实现）
def _transmission(p):
    return sum(w * math.exp(-f * p) for _k, w, f in SPECTRUM)


def _chord(radius, dist):
    return 0.0 if abs(dist) >= radius else 2.0 * math.sqrt(radius * radius - dist * dist)


def _line_integral(phantom, s, theta_deg):
    th = math.radians(theta_deg)
    c, sn = math.cos(th), math.sin(th)
    bg = phantom[0]["mu"]
    total = bg * _chord(phantom[0]["r"], s - (phantom[0]["x"] * c + phantom[0]["y"] * sn)) / 10.0
    for disc in phantom[1:]:
        total += (disc["mu"] - bg) * _chord(disc["r"], s - (disc["x"] * c + disc["y"] * sn)) / 10.0
    return total


# 低对比度工况：夹杂只比基体高 2.4%，靠单能谱近似或漏做居中校正会把它抹平
LOW_PHANTOM = [
    {"x": 0.0, "y": 0.0, "r": 60.0, "mu": 0.50},
    {"x": -18.0, "y": 12.0, "r": 26.0, "mu": 0.512},
]
LOW_ROIS = [
    {"name": "matrix", "x_mm": 10.0, "y_mm": -30.0, "r_mm": 8.0},
    {"name": "insert", "x_mm": -18.0, "y_mm": 12.0, "r_mm": 10.0},
]

def _scatter(a, x, i0a):
    """每帧一层缓变的散射本底：按道号是二次多项式，幅度是源强的几个百分点。"""
    amp = 0.02 * (1.0 + 0.6 * math.sin(2 * math.pi * a / 13.0)
                  + 0.3 * math.cos(2 * math.pi * a / 29.0))
    return i0a * amp * (0.6 + 0.35 * x + 0.45 * x * x)


def make_scan(phantom, rois, size=64, px=2.0, channels=81, pitch=2.0, n_angles=90,
              i0=20000.0, dark=50.0, bad=(), silent=(), drift=0.03, inspect_r=None,
              center_offset=0.0, air_skip=8):
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
            x = (s0 + k * pitch) / (max(abs(s0), abs(s0 + (channels - 1) * pitch)))
            cnt = int(round(i0a * _transmission(depth) + dark + _scatter(a, x, i0a)))
            if k in bad or k in silent:
                cnt = 0
            row.append(cnt)
        counts.append(row)
    return {
        "instrument": {"source_photons": i0, "dark_counts": dark,
                       "air_skip_channels": air_skip},
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


def _run(scan_path, out_path, timeout=240):
    cmd = [sys.executable, ENTRY, "--scan", scan_path, "--out", out_path]
    return subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                          errors="replace", timeout=timeout)


def run_scan(scan_obj, out_path=None):
    tmp = tempfile.mkdtemp(prefix="ct-case-")
    scan_path = os.path.join(tmp, "scan.json")
    target = out_path or os.path.join(tmp, "report.json")
    with open(scan_path, "w", encoding="utf-8") as f:
        json.dump(scan_obj, f, ensure_ascii=False)
    proc = _run(scan_path, target)
    report = None
    if os.path.exists(target):
        try:
            with open(target, encoding="utf-8") as f:
                report = json.load(f)
        except ValueError:
            report = None
    return report, {"scan": scan_path, "out": target, "proc": proc, "dir": tmp}


def _roi_values(report):
    return {r["name"]: r["mu_per_cm"] for r in report["roi"]}


# ---------------------------------------------------------------- 测试
def test_roi_mu_matches_phantom_on_base_scan():
    """题面原文：「采区均值与缺陷面积分别按三位和一位小数给」

    基体与夹杂按相对 2% 卡，空洞按绝对 0.005 /cm 卡。
    """
    report, _ = run_scan(make_scan(BASE_PHANTOM, BASE_ROIS, bad=(27, 28, 79), silent=(52,)))
    assert isinstance(report, dict), "没有产出报告"
    got = _roi_values(report)
    assert "matrix" in got, "报告缺采区 matrix"
    assert "insert" in got, "报告缺采区 insert"
    assert "void" in got, "报告缺采区 void"
    matrix = float(got["matrix"])
    insert = float(got["insert"])
    void = float(got["void"])
    assert abs(matrix - 0.50) <= 0.012 * 0.50, "基体偏差 %.2f%% 超容差" % ((matrix - 0.50) / 0.50 * 100)
    assert abs(insert - 0.62) <= 0.012 * 0.62, "夹杂偏差 %.2f%% 超容差" % ((insert - 0.62) / 0.62 * 100)
    assert abs(void - 0.05) <= 0.012, "空洞偏差 %.4f 超过绝对容差" % (void - 0.05)


def test_roi_mu_matches_second_phantom_geometry():
    """题面原文：「复检时换的批次连几何、材料、道数与角度都不一样」。

    换几何与材料后，三个采区仍要落在容差内。
    """
    report, _ = run_scan(make_scan(SECOND_PHANTOM, SECOND_ROIS, channels=81, pitch=2.0,
                                   n_angles=90, i0=15000.0, dark=120.0, bad=(37,),
                                   silent=(30,), center_offset=2.6))
    assert isinstance(report, dict), "没有产出报告"
    got = _roi_values(report)
    assert "matrix" in got, "报告缺采区 matrix"
    assert "insert" in got, "报告缺采区 insert"
    assert "void" in got, "报告缺采区 void"
    matrix = float(got["matrix"])
    insert = float(got["insert"])
    void = float(got["void"])
    assert abs(matrix - 0.44) <= 0.012 * 0.44, "第二批基体偏差 %.2f%%" % ((matrix - 0.44) / 0.44 * 100)
    assert abs(insert - 0.60) <= 0.012 * 0.60, "第二批夹杂偏差 %.2f%%" % ((insert - 0.60) / 0.60 * 100)
    assert abs(void - 0.15) <= 0.012, "第二批空洞偏差 %.4f" % (void - 0.15)


def test_dark_current_large_offset_recovered():
    """题面原文：「从计数恢复射线穿过样品的线积分」。

    暗电流抬到源强的百分之四后，采区仍要落在容差内。
    """
    report, _ = run_scan(make_scan(BASE_PHANTOM, BASE_ROIS, dark=800.0, bad=(27,),
                                   silent=(77,), center_offset=-3.4))
    assert isinstance(report, dict), "没有产出报告"
    got = _roi_values(report)
    assert abs(float(got["matrix"]) - 0.50) <= 0.012 * 0.50, "暗电流工况下基体偏差过大"
    assert abs(float(got["insert"]) - 0.62) <= 0.012 * 0.62, "暗电流工况下夹杂偏差过大"
    assert abs(float(got["void"]) - 0.05) <= 0.005, "暗电流工况下空洞偏差过大"


def test_low_contrast_insert_resolved():
    """题面原文：「均匀材料上的重建值应当等于该材料在参考能量下的线性衰减系数」。

    夹杂只比基体高 2.4%（0.512 对 0.500）：夹杂的衰减系数要落在 0.8% 以内，
    两区之差要保持在 0.009 以上；同时它远高于缺陷阈值，不得被报成缺陷。
    """
    report, _ = run_scan(make_scan(LOW_PHANTOM, LOW_ROIS, bad=(27,), silent=(52,)))
    assert isinstance(report, dict), "没有产出报告"
    got = _roi_values(report)
    matrix = float(got["matrix"])
    insert = float(got["insert"])
    assert abs(insert - 0.512) <= 0.008 * 0.512, \
        "低对比度夹杂偏差 %.4f 超容差" % (insert - 0.512)
    assert abs(matrix - 0.500) <= 0.008 * 0.500, \
        "低对比度基体偏差 %.4f 超容差" % (matrix - 0.500)
    assert insert - matrix >= 0.009, \
        "两区对比度被抹平：insert-matrix = %.4f（真值 0.012）" % (insert - matrix)
    assert len(report["defects"]) == 0, \
        "低对比度夹杂被误报成缺陷：%r" % (report["defects"],)


def test_dead_channel_recovery_with_unlisted_zero():
    """题面原文：「计数为 0 的道一律按失效道处理」。

    自检清单为空而存在 4 个零计数道时，采区仍要合格。
    """
    report, _ = run_scan(make_scan(BASE_PHANTOM, BASE_ROIS, bad=(), silent=(20, 21, 22, 60),
                                   center_offset=-3.4))
    assert isinstance(report, dict), "没有产出报告"
    got = _roi_values(report)
    assert abs(float(got["matrix"]) - 0.50) <= 0.012 * 0.50
    assert abs(float(got["insert"]) - 0.62) <= 0.012 * 0.62


def test_defect_area_and_center_match_phantom():
    """题面原文：「是检视范围内 µ 低于」。

    只有一个 9 mm 半径的空洞：恰好检出 1 个，质心不超过 2 mm，面积在真值 ±12% 内。
    """
    report, _ = run_scan(make_scan(BASE_PHANTOM, BASE_ROIS, bad=(27, 28, 79), silent=(52,)))
    assert isinstance(report, dict), "没有产出报告"
    defects = report["defects"]
    assert len(defects) == 1, "检出缺陷数 %d，应为 1" % len(defects)
    d = defects[0]
    truth_area = math.pi * 81.0
    assert math.hypot(float(d["x_mm"]) - 25.0, float(d["y_mm"]) + 10.0) <= 1.5, "缺陷质心偏差过大"
    assert abs(float(d["area_mm2"]) - truth_area) <= 0.12 * truth_area, "缺陷面积偏差过大"


def test_no_spurious_defect_outside_inspection_radius():
    """题面原文：「检视范围内找内部缺陷」。

    检视半径收到 52 mm（样品半径 60 mm）时，样品外的空气不得被报成缺陷。
    """
    report, _ = run_scan(make_scan(BASE_PHANTOM, BASE_ROIS, bad=(27, 28, 79), silent=(52,),
                                   inspect_r=52.0))
    assert isinstance(report, dict), "没有产出报告"
    assert len(report["defects"]) == 1, "检视范围外被报成了缺陷"
    assert 150.0 <= float(report["defects"][0]["area_mm2"]) <= 350.0, "缺陷面积不合理"


def test_repaired_channels_listed_exactly():
    """题面原文：「是本次判为失效并补过的道号，升序去重」。

    自检列了 12 与 40、另有一道 60 只有零计数：三者都要列出。
    """
    report, _ = run_scan(make_scan(BASE_PHANTOM, BASE_ROIS, bad=(12, 40), silent=(60,),
                                   center_offset=-3.4))
    assert isinstance(report, dict), "没有产出报告"
    assert report["meta"]["repaired_channels"] == [12, 40, 60], "补道清单不对：%r" % (
        report["meta"]["repaired_channels"],)


def test_invalid_scan_rejected_without_output():
    """题面原文：「并在 stderr 留一句人话，报告文件一个字都不写」。

    六种坏输入（含边界值）逐个验证：退出码 2、无报告、无临时文件、无日志/小结副作用、
    已有报告字节不变、stderr 非空。
    """
    good = make_scan(BASE_PHANTOM, BASE_ROIS, bad=(27,), silent=(), center_offset=-3.4)
    missing = json.loads(json.dumps(good))
    missing.pop("projections")
    wrong_ch = json.loads(json.dumps(good))
    wrong_ch["geometry"]["channels"] = 100
    bad_instrument = json.loads(json.dumps(good))
    bad_instrument["instrument"]["dark_counts"] = bad_instrument["instrument"]["source_photons"]
    zero_px = json.loads(json.dumps(good))
    zero_px["geometry"]["pixel_size_mm"] = 0.0
    tiny_image = json.loads(json.dumps(good))
    tiny_image["geometry"]["image_size"] = 3
    bad_channel = json.loads(json.dumps(good))
    bad_channel["detector"]["bad_channels"] = [good["geometry"]["channels"]]

    tmp = tempfile.mkdtemp(prefix="ct-neg-")
    keep_path = os.path.join(tmp, "report.json")
    report, _ = run_scan(good, out_path=keep_path)
    assert isinstance(report, dict), "前置的正常扫描没跑通"
    before = open(keep_path, "rb").read()

    for scan in (missing, wrong_ch, bad_instrument, zero_px, tiny_image, bad_channel):
        case = tempfile.mkdtemp(prefix="ct-bad-")
        sp = os.path.join(case, "scan.json")
        op = os.path.join(case, "report.json")
        with open(sp, "w", encoding="utf-8") as f:
            json.dump(scan, f, ensure_ascii=False)
        proc = _run(sp, op)
        assert proc.returncode == 2, "坏输入退出码是 %r，应为 2" % proc.returncode
        assert not os.path.exists(op), "坏输入却写出了报告"
        assert (proc.stderr or "").strip(), "坏输入没有给出可读原因"
        leftovers = [n for n in os.listdir(case) if n != "scan.json"]
        assert not leftovers, "坏输入留下了副作用文件：%r" % leftovers

    bad_input = os.path.join(tmp, "bad.json")
    with open(bad_input, "w", encoding="utf-8") as f:
        json.dump(wrong_ch, f, ensure_ascii=False)
    proc = _run(bad_input, keep_path)
    assert proc.returncode == 2
    assert open(keep_path, "rb").read() == before, "坏输入覆盖了上一次的报告"


def test_report_schema_and_rounding():
    """题面原文：「结论写成一个 UTF-8 的 JSON 对象，字段与层次如下」

    校验顶层键、采区顺序、位数、meta 计数与 JSON 合法性。
    """
    scan_obj = make_scan(BASE_PHANTOM, BASE_ROIS, bad=(27, 28, 79), silent=(52,))
    report, ctx = run_scan(scan_obj)
    assert isinstance(report, dict)
    assert set(report) == {"geometry", "roi", "defects", "meta"}, "顶层键不对：%r" % (set(report),)
    assert report["geometry"]["image_size"] == 64
    assert abs(float(report["geometry"]["pixel_size_mm"]) - 2.0) <= 1e-9, "像素尺寸不是 2.0"
    assert [r["name"] for r in report["roi"]] == [r["name"] for r in BASE_ROIS], "采区顺序不对"
    for r in report["roi"]:
        mu = float(r["mu_per_cm"])
        assert math.isfinite(mu)
        assert abs(mu - round(mu, 3)) <= 1e-12, "采区 %s 不是三位小数" % r["name"]
        assert isinstance(r["pixels"], int) and r["pixels"] > 0
    assert report["meta"]["n_angles"] == 90
    assert report["meta"]["n_channels"] == 81
    raw = open(ctx["out"], encoding="utf-8").read()
    assert "NaN" not in raw and "Infinity" not in raw


def test_rerun_produces_identical_bytes():
    """题面原文：「同一份扫描再跑一遍，结论文件要一个字节都不差」。

    连跑两次，两个报告文件字节完全相同。
    """
    scan_obj = make_scan(BASE_PHANTOM, BASE_ROIS, bad=(27, 28, 79), silent=(52,))
    first, ctx1 = run_scan(scan_obj)
    assert isinstance(first, dict), "第一次就没跑通"
    body1 = open(ctx1["out"], "rb").read()
    second, ctx2 = run_scan(scan_obj)
    assert isinstance(second, dict), "第二次没跑通"
    assert body1 == open(ctx2["out"], "rb").read(), "两次运行的报告字节不一致"


def test_failed_run_keeps_previous_report():
    """题面原文：「这次跑失败，上一版报告要原封不动地留着」。

    先跑合格扫描留下报告，再用道数不匹配的坏输入跑：退出码 2，旧报告字节不变。
    """
    good = make_scan(BASE_PHANTOM, BASE_ROIS, bad=(27,), silent=(), center_offset=-3.4)
    tmp = tempfile.mkdtemp(prefix="ct-keep-")
    sp = os.path.join(tmp, "scan.json")
    op = os.path.join(tmp, "report.json")
    with open(sp, "w", encoding="utf-8") as f:
        json.dump(good, f, ensure_ascii=False)
    first = _run(sp, op)
    assert first.returncode == 0 and os.path.exists(op), "前置的正常扫描没跑通"
    before = open(op, "rb").read()
    broken = json.loads(json.dumps(good))
    broken["geometry"]["channels"] = 100
    with open(sp, "w", encoding="utf-8") as f:
        json.dump(broken, f, ensure_ascii=False)
    second = _run(sp, op)
    assert second.returncode == 2
    assert open(op, "rb").read() == before, "失败的那次覆盖了旧报告"


def test_zero_reference_energy_rejected():
    """题面原文：「取 0 或负数同样算参数越界」。

    把参考能量设成 0（正好卡在边界）：必须退出码 2 且不产出报告。
    """
    data = make_scan(BASE_PHANTOM, BASE_ROIS, bad=(27,), silent=(52,))
    data["material_card"]["reference_energy_kev"] = 0.0
    tmp = tempfile.mkdtemp(prefix="ct-ref0-")
    spath = os.path.join(tmp, "scan.json")
    opath = os.path.join(tmp, "report.json")
    with open(spath, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False)
    proc = _run(spath, opath)
    assert proc.returncode == 2, "参考能量为 0 时退出码应为 2，实际 %r" % proc.returncode
    assert not os.path.exists(opath), "参考能量为 0 时不该产出报告"