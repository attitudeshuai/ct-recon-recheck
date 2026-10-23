#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""平行束 CT 截面重建与缺陷复核（参考实现）。

命令行：
    python3 src/recon.py --scan data/scan.json --out out/report.json

口径：
  · 计数 → 透射率要先扣暗电流、按每角度源强归一，再用材料卡的能谱做多能谱线性化反演：
    T = (I - dark)/(I0(θ) - dark)，T = Σ w_k·exp(-f_k·P)，反解出 70 keV 下的线积分 P；
  · 滤波反投影：Ram-Lak 核，反投影权重 2·Δθ（均匀圆盘定标，见 solution/README.md）；
  · 重建值是线性衰减系数 µ（1/cm），基体标称值见材料卡；
  · 缺陷 = 检视半径内 µ 低于阈值的连通块（8 邻接），面积按 mm² 报；
  · 输出先写临时文件再原子替换；任何非法输入以退出码 2 结束且不产出报告。
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import tempfile

EXIT_INVALID = 2


# ---------------------------------------------------------------- 输入校验
def _fail(msg):
    sys.stderr.write("非法输入：%s\n" % msg)
    raise SystemExit(EXIT_INVALID)


def _get(obj, key, kind, where):
    if not isinstance(obj, dict) or key not in obj:
        _fail("%s 缺字段 %s" % (where, key))
    val = obj[key]
    if kind is float:
        if not isinstance(val, (int, float)) or isinstance(val, bool):
            _fail("%s.%s 不是数值" % (where, key))
        return float(val)
    if kind is int:
        if not isinstance(val, int) or isinstance(val, bool):
            _fail("%s.%s 不是整数" % (where, key))
        return val
    if not isinstance(val, kind):
        _fail("%s.%s 类型不对" % (where, key))
    return val


def load_scan(path):
    try:
        with open(path, encoding="utf-8") as f:
            raw = json.load(f)
    except FileNotFoundError:
        _fail("读不到扫描文件 %s" % path)
    except ValueError as e:
        _fail("扫描文件不是合法 JSON：%s" % e)

    geo = _get(raw, "geometry", dict, "scan")
    inst = _get(raw, "instrument", dict, "scan")
    proj = _get(raw, "projections", dict, "scan")
    insp = _get(raw, "inspection", dict, "scan")

    size = _get(geo, "image_size", int, "geometry")
    px = _get(geo, "pixel_size_mm", float, "geometry")
    n_ch = _get(geo, "channels", int, "geometry")
    pitch = _get(geo, "channel_pitch_mm", float, "geometry")
    angles = _get(geo, "angles_deg", list, "geometry")
    if size < 4:
        _fail("geometry.image_size 太小")
    if not (px > 0) or not (pitch > 0):
        _fail("geometry 的像素尺寸与道间距必须为正")
    if len(angles) < 3:
        _fail("geometry.angles_deg 少于 3 个角度")
    for a in angles:
        if not isinstance(a, (int, float)) or isinstance(a, bool) or not math.isfinite(float(a)):
            _fail("geometry.angles_deg 里有非数值")
    step = 180.0 / len(angles)
    for i, a in enumerate(angles):
        if abs(float(a) - i * step) > 0.01:
            _fail("geometry.angles_deg 不是从 0 起、步长 180/N 的均匀半周采样")

    i0 = _get(inst, "source_photons", float, "instrument")
    dark = _get(inst, "dark_counts", float, "instrument")
    if not (i0 > 0):
        _fail("instrument.source_photons 必须为正")
    if dark < 0:
        _fail("instrument.dark_counts 不能为负")
    if not (i0 > dark):
        _fail("instrument.source_photons 不大于 dark_counts")

    per_angle = _get(proj, "per_angle_source", list, "projections")
    counts = _get(proj, "counts", list, "projections")
    n_ang = len(angles)
    if len(per_angle) != n_ang:
        _fail("per_angle_source 有 %d 个值，和 angles_deg 的 %d 个角度对不上" % (len(per_angle), n_ang))
    if len(counts) != n_ang:
        _fail("counts 有 %d 行，和 angles_deg 的 %d 个角度对不上" % (len(counts), n_ang))
    for i, v in enumerate(per_angle):
        if not isinstance(v, (int, float)) or isinstance(v, bool):
            _fail("per_angle_source[%d] 不是数值" % i)
        if not (float(v) > dark):
            _fail("per_angle_source[%d] 不大于暗电流" % i)
    for r, row in enumerate(counts):
        if not isinstance(row, list) or len(row) != n_ch:
            _fail("counts 第 %d 行不是 %d 道" % (r + 1, n_ch))
        for c, v in enumerate(row):
            if not isinstance(v, (int, float)) or isinstance(v, bool):
                _fail("counts[%d][%d] 不是数值" % (r, c))
            if float(v) < 0:
                _fail("counts[%d][%d] 为负" % (r, c))
            if float(v) - dark > float(per_angle[r]) * 1.10 + 100.0:
                _fail("counts[%d][%d] 扣掉暗电流后仍大于该角度源强（曝光异常）" % (r, c))

    r_insp = _get(insp, "radius_mm", float, "inspection")
    thr = _get(insp, "defect_mu_threshold_per_cm", float, "inspection")
    min_area = _get(insp, "min_defect_area_mm2", float, "inspection")
    if not (r_insp > 0) or min_area < 0:
        _fail("inspection 的检视半径必须为正、最小缺陷面积不能为负")

    rois = _get(raw, "rois", list, "scan")
    if not rois:
        _fail("rois 为空")
    clean_rois = []
    for i, r in enumerate(rois):
        if not isinstance(r, dict):
            _fail("rois[%d] 不是对象" % i)
        name = _get(r, "name", str, "rois[%d]" % i)
        cx = _get(r, "x_mm", float, "rois[%d]" % i)
        cy = _get(r, "y_mm", float, "rois[%d]" % i)
        rr = _get(r, "r_mm", float, "rois[%d]" % i)
        if not (rr > 0):
            _fail("rois[%d].r_mm 必须为正" % i)
        clean_rois.append({"name": name, "x_mm": cx, "y_mm": cy, "r_mm": rr})

    card = _get(raw, "material_card", dict, "scan")
    if len(angles) == 0:
        _fail("geometry.angles_deg 是空的")
    spectrum_raw = _get(card, "spectrum", list, "material_card")
    if len(spectrum_raw) < 2:
        _fail("material_card.spectrum 至少要两个能量点")
    spectrum = []
    for i, item in enumerate(spectrum_raw):
        if not isinstance(item, dict):
            _fail("material_card.spectrum[%d] 不是对象" % i)
        kev = _get(item, "kev", float, "spectrum[%d]" % i)
        weight = _get(item, "weight", float, "spectrum[%d]" % i)
        rel = _get(item, "mu_relative", float, "spectrum[%d]" % i)
        if not (kev > 0) or not (weight > 0) or not (rel > 0):
            _fail("material_card.spectrum[%d] 的 kev/weight/mu_relative 都要为正" % i)
        spectrum.append((kev, weight, rel))
    if abs(sum(w for _k, w, _f in spectrum) - 1.0) > 0.05:
        _fail("material_card.spectrum 的权重和应接近 1")
    ref_kev = _get(card, "reference_energy_kev", float, "material_card")
    if ref_kev <= 0:
        _fail("material_card.reference_energy_kev 必须为正")

    detector = raw.get("detector") or {}
    bad = detector.get("bad_channels") or []
    clean_bad = []
    for c in bad:
        if not isinstance(c, int) or isinstance(c, bool) or c < 0 or not (c < n_ch):
            _fail("detector.bad_channels 里有越界道号")
        clean_bad.append(c)

    return {
        "size": size, "px": px, "n_ch": n_ch, "pitch": pitch,
        "angles_deg": [float(a) for a in angles], "i0": i0, "dark": dark,
        "per_angle": [float(v) for v in per_angle],
        "counts": [[float(v) for v in row] for row in counts],
        "inspect_r": r_insp, "thr": thr, "min_area": min_area,
        "rois": clean_rois, "bad": sorted(set(clean_bad)),
        "air_skip": int(inst.get("air_skip_channels", 0)),
        "spectrum": spectrum, "ref_kev": ref_kev,
    }


# ---------------------------------------------------------------- 多能谱线性化
def _transmission_of_p(spectrum, p):
    return sum(w * math.exp(-f * p) for _kev, w, f in spectrum)


def build_inverse(spectrum, p_max=40.0, n=20001):
    """建 P ↔ T 的单调表；T 随 P 递减，反解用二分 + 线性插值。"""
    ps = [p_max * i / (n - 1) for i in range(n)]
    ts = [_transmission_of_p(spectrum, v) for v in ps]
    return ps, ts


def optical_depth(ts, ps, t):
    """由透射率反解 70 keV 线积分：单能假设下这是一条过原点的直线，多能谱下是条曲线。"""
    if not (t < ts[0]):
        return 0.0
    if t <= ts[-1]:
        return ps[-1]
    lo, hi = 0, len(ts) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ts[mid] > t:
            lo = mid
        else:
            hi = mid
    t0, t1 = ts[lo], ts[hi]
    if not (t1 > t0):
        return ps[lo]
    return ps[lo] + (t0 - t) / (t0 - t1) * (ps[hi] - ps[lo])


# ---------------------------------------------------------------- 投影 → 线积分
def _solve_linear(a, b):
    n = len(a)
    m = [a[i][:] + [b[i]] for i in range(n)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(m[r][col]))
        if abs(m[piv][col]) < 1e-14:
            return None
        m[col], m[piv] = m[piv], m[col]
        pv = m[col][col]
        for k in range(col, n + 1):
            m[col][k] /= pv
        for r in range(n):
            if r == col:
                continue
            fac = m[r][col]
            if fac:
                for k in range(col, n + 1):
                    m[r][k] -= fac * m[col][k]
    return [m[i][n] for i in range(n)]


def _air_background(row, dark, air_skip, bad, i0):
    """两端空气道只看得到「源强 + 散射本底」，据此拟二次多项式并返回逐道本底。"""
    n_ch = len(row)
    xs, ys = [], []
    for k in list(range(air_skip)) + list(range(n_ch - air_skip, n_ch)):
        if k in bad:
            continue
        v = row[k] - dark - i0
        if v <= -0.5 * i0:
            continue
        xs.append(k - (n_ch - 1) / 2.0)
        ys.append(v)
    if len(xs) < 4:
        return [0.0] * n_ch
    xtx = [[0.0] * 3 for _ in range(3)]
    xty = [0.0] * 3
    for x, y in zip(xs, ys):
        basis = [1.0, x, x * x]
        for i in range(3):
            xty[i] += basis[i] * y
            for j in range(3):
                xtx[i][j] += basis[i] * basis[j]
    coef = _solve_linear(xtx, xty)
    if coef is None:
        return [0.0] * n_ch
    out = []
    for k in range(n_ch):
        x = k - (n_ch - 1) / 2.0
        out.append(max(0.0, coef[0] + coef[1] * x + coef[2] * x * x))
    return out


def to_line_integrals(scan):
    """先按空气道扣散射本底与暗电流，再按每角度源强归一；坏道用两侧好道插值。"""
    dark = scan["dark"]
    n_ch = scan["n_ch"]
    air_skip = scan.get("air_skip", 0)
    bad_set = set(scan["bad"])
    ps, ts = build_inverse(scan["spectrum"])
    p, repaired = [], set()
    for a, row in enumerate(scan["counts"]):
        i0 = scan["per_angle"][a]
        bg = _air_background(row, dark, air_skip, bad_set, i0) if air_skip > 1 else [0.0] * n_ch
        row = [row[k] - bg[k] for k in range(n_ch)]
        line = []
        for k, cnt in enumerate(row):
            denom = cnt - dark
            if denom <= 0.0:          # 该道没有有效信号：先占位，后面统一插值
                line.append(None)
                repaired.add(k)
                continue
            line.append(optical_depth(ts, ps, min(1.0, denom / (i0 - dark))))
        p.append(line)
    repaired |= set(scan["bad"])
    for k in scan["bad"]:
        for a in range(len(p)):
            p[a][k] = None
    for a in range(len(p)):
        for k in range(n_ch):
            if p[a][k] is not None:
                continue
            lo, hi = k - 1, k + 1
            while not (lo < 0) and p[a][lo] is None:
                lo -= 1
            while hi < n_ch and p[a][hi] is None:
                hi += 1
            if lo < 0 and not (hi < n_ch):
                p[a][k] = 0.0
            elif lo < 0:
                p[a][k] = p[a][hi]
            elif not (hi < n_ch):
                p[a][k] = p[a][lo]
            else:
                w = (k - lo) / float(hi - lo)
                p[a][k] = p[a][lo] * (1.0 - w) + p[a][hi] * w
    return p, sorted(repaired)


# ---------------------------------------------------------------- 转台轴偏心
def estimate_axis_offset(p, scan):
    """用投影质心的角度均值估转台轴偏心（样品质心投影在整周上正负抵消，均值即偏心量）。"""
    ds = scan["pitch"]
    s0 = -(scan["n_ch"] - 1) / 2.0 * ds
    total, count = 0.0, 0
    for row in p:
        weight = sum(row)
        if weight <= 0.0:
            continue
        centroid = sum((s0 + k * ds) * v for k, v in enumerate(row)) / weight
        total += centroid
        count += 1
    return total / count if count else 0.0


# ---------------------------------------------------------------- 滤波反投影
def ramp_kernel(n_ch, ds):
    """Ram-Lak：h(0)=1/(4Δs²)，奇数位 −1/(π²n²Δs²)，偶数位 0。"""
    half = n_ch - 1
    h = [0.0] * (2 * half + 1)
    h[half] = 1.0 / (4.0 * ds * ds)
    for n in range(1, half + 1):
        if n % 2 == 1:
            v = -1.0 / (math.pi * math.pi * n * n * ds * ds)
            h[half + n] = v
            h[half - n] = v
    return h, half


def reconstruct(p, scan, axis_offset=0.0):
    n_ch, ds = scan["n_ch"], scan["pitch"]
    h, half = ramp_kernel(n_ch, ds)
    s0 = -(n_ch - 1) / 2.0 * ds
    n_ang = len(p)
    dtheta = math.pi / n_ang
    # 2·Δθ 是离散反投影的定标：单色均匀圆盘重建值必须等于其衰减系数
    weight = 2.0 * dtheta
    size, px = scan["size"], scan["px"]
    img = [[0.0] * size for _ in range(size)]
    for a in range(n_ang):
        row = p[a]
        q = [0.0] * n_ch
        for k in range(n_ch):
            acc = 0.0
            lo = max(0, k - half)
            hi = min(n_ch - 1, k + half)
            for kk in range(lo, hi + 1):
                acc += row[kk] * h[half + (k - kk)]
            q[k] = acc
        th = math.radians(scan["angles_deg"][a])
        c, sn = math.cos(th), math.sin(th)
        for iy in range(size):
            y = (iy - (size - 1) / 2.0) * px
            base = y * sn - s0
            for ix in range(size):
                x = (ix - (size - 1) / 2.0) * px
                t = (x * c + base + axis_offset) / ds
                k0 = int(math.floor(t))
                if k0 < 0 or not (k0 < n_ch - 1):
                    continue
                f = t - k0
                img[iy][ix] += weight * ((1.0 - f) * q[k0] + f * q[k0 + 1])
    return img


# ---------------------------------------------------------------- 汇总
def to_per_cm(img):
    """内部按 1/mm 计算，输出统一换成 1/cm。"""
    return [[v * 10.0 for v in row] for row in img]


def roi_report(img_cm, scan):
    size, px = scan["size"], scan["px"]
    out = []
    for roi in scan["rois"]:
        vals = []
        for iy in range(size):
            y = (iy - (size - 1) / 2.0) * px
            dy2 = (y - roi["y_mm"]) ** 2
            for ix in range(size):
                x = (ix - (size - 1) / 2.0) * px
                if (x - roi["x_mm"]) ** 2 + dy2 <= roi["r_mm"] ** 2:
                    vals.append(img_cm[iy][ix])
        mu = sum(vals) / len(vals) if vals else 0.0
        out.append({"name": roi["name"], "mu_per_cm": round(mu, 3), "pixels": len(vals)})
    return out


def defect_report(img_cm, scan):
    size, px = scan["size"], scan["px"]
    r_insp, thr, min_area = scan["inspect_r"], scan["thr"], scan["min_area"]
    seen = [[False] * size for _ in range(size)]

    def inside(ix, iy):
        x = (ix - (size - 1) / 2.0) * px
        y = (iy - (size - 1) / 2.0) * px
        return x * x + y * y <= r_insp * r_insp

    out = []
    for iy in range(size):
        for ix in range(size):
            if seen[iy][ix] or not inside(ix, iy) or not (img_cm[iy][ix] < thr):
                continue
            stack, comp = [(ix, iy)], []
            seen[iy][ix] = True
            while stack:
                cx, cy = stack.pop()
                comp.append((cx, cy))
                for dx in (-1, 0, 1):
                    for dy in (-1, 0, 1):
                        nx, ny = cx + dx, cy + dy
                        if 0 <= nx < size and 0 <= ny < size and not seen[ny][nx] \
                                and inside(nx, ny) and img_cm[ny][nx] < thr:
                            seen[ny][nx] = True
                            stack.append((nx, ny))
            area = len(comp) * px * px
            if area < min_area:
                continue
            gx = sum(c[0] for c in comp) / float(len(comp))
            gy = sum(c[1] for c in comp) / float(len(comp))
            out.append({
                "x_mm": round((gx - (size - 1) / 2.0) * px, 2),
                "y_mm": round((gy - (size - 1) / 2.0) * px, 2),
                "area_mm2": round(area, 1),
                "min_mu_per_cm": round(min(img_cm[c[1]][c[0]] for c in comp), 3),
            })
    return sorted(out, key=lambda d: (-d["area_mm2"], d["x_mm"], d["y_mm"]))


def build_report(img_cm, scan, repaired):
    return {
        "geometry": {"image_size": scan["size"], "pixel_size_mm": scan["px"]},
        "roi": roi_report(img_cm, scan),
        "defects": defect_report(img_cm, scan),
        "meta": {
            "n_angles": len(scan["angles_deg"]),
            "n_channels": scan["n_ch"],
            "repaired_channels": repaired,
        },
    }


def atomic_write(path, obj):
    out_dir = os.path.dirname(os.path.abspath(path))
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    body = json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=False) + "\n"
    fd, tmp = tempfile.mkstemp(prefix=".recon-", suffix=".tmp", dir=out_dir or ".")
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as f:
            f.write(body)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def main(argv=None):
    ap = argparse.ArgumentParser(description="平行束 CT 重建与缺陷复核")
    ap.add_argument("--scan", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    scan = load_scan(args.scan)
    p, repaired = to_line_integrals(scan)
    axis_offset = estimate_axis_offset(p, scan)
    img = to_per_cm(reconstruct(p, scan, axis_offset))
    report = build_report(img, scan, repaired)
    atomic_write(args.out, report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
