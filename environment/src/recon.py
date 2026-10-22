#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""工业 CT 截面重建与缺陷复核。

读扫描站导出的 JSON，重建截面衰减系数，输出采区平均值与缺陷清单。

用法：
    python3 src/recon.py --scan data/scan.json --out out/report.json
"""
import argparse
import json
import math
import os


def load_scan(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def to_line_integrals(scan):
    """计数转线积分：每个角度做一次对数变换。"""
    source = scan["instrument"]["source_photons"]
    lines = []
    for row in scan["projections"]["counts"]:
        lines.append([-math.log(count / source) for count in row])
    return lines


def ramp_kernel(n_ch, pitch):
    """Ram-Lak 滤波核。"""
    half = n_ch - 1
    kernel = [0.0] * (2 * half + 1)
    kernel[half] = 1.0 / (4.0 * pitch * pitch)
    for n in range(1, half + 1):
        if n % 2 == 1:
            value = -1.0 / (math.pi * math.pi * n * n * pitch * pitch)
            kernel[half + n] = value
            kernel[half - n] = value
    return kernel, half


def reconstruct(lines, scan):
    geo = scan["geometry"]
    n_ch = geo["channels"]
    pitch = geo["channel_pitch_mm"]
    size = geo["image_size"]
    px = geo["pixel_size_mm"]
    n_ang = len(lines)

    kernel, half = ramp_kernel(n_ch, pitch)
    s0 = -(n_ch - 1) / 2.0 * pitch
    dtheta = math.pi / n_ang

    image = [[0.0] * size for _ in range(size)]
    for a in range(n_ang):
        row = lines[a]
        filtered = [0.0] * n_ch
        for k in range(n_ch):
            acc = 0.0
            for j in range(-half, half + 1):
                kk = k - j
                if 0 <= kk < n_ch:
                    acc += row[kk] * kernel[half + j]
            filtered[k] = acc

        theta = math.radians(geo["angles_deg"][a])
        cos_t, sin_t = math.cos(theta), math.sin(theta)
        for iy in range(size):
            y = (iy - (size - 1) / 2.0) * px
            for ix in range(size):
                x = (ix - (size - 1) / 2.0) * px
                t = (x * cos_t + y * sin_t - s0) / pitch
                k0 = int(math.floor(t))
                if 0 <= k0 < n_ch - 1:
                    frac = t - k0
                    image[iy][ix] += dtheta * ((1.0 - frac) * filtered[k0] + frac * filtered[k0 + 1])

    for iy in range(size):
        for ix in range(size):
            image[iy][ix] *= 10.0
    return image


def summarize(image, scan):
    geo = scan["geometry"]
    size = geo["image_size"]
    px = geo["pixel_size_mm"]

    roi_rows = []
    for roi in scan["rois"]:
        values = []
        for iy in range(size):
            y = (iy - (size - 1) / 2.0) * px
            for ix in range(size):
                x = (ix - (size - 1) / 2.0) * px
                if (x - roi["x_mm"]) ** 2 + (y - roi["y_mm"]) ** 2 <= roi["r_mm"] ** 2:
                    values.append(image[iy][ix])
        mean = sum(values) / len(values)
        roi_rows.append({"name": roi["name"], "mu_per_cm": round(mean, 3), "pixels": len(values)})

    threshold = scan["inspection"]["defect_mu_threshold_per_cm"]
    seen = [[False] * size for _ in range(size)]
    defects = []
    for iy in range(size):
        for ix in range(size):
            if seen[iy][ix] or image[iy][ix] >= threshold:
                continue
            stack = [(ix, iy)]
            seen[iy][ix] = True
            cells = []
            while stack:
                cx, cy = stack.pop()
                cells.append((cx, cy))
                for dx in (-1, 0, 1):
                    for dy in (-1, 0, 1):
                        nx, ny = cx + dx, cy + dy
                        if 0 <= nx < size and 0 <= ny < size and not seen[ny][nx] \
                                and image[ny][nx] < threshold:
                            seen[ny][nx] = True
                            stack.append((nx, ny))
            if len(cells) < scan["inspection"]["min_defect_area_mm2"]:
                continue
            mx = sum(c[0] for c in cells) / float(len(cells))
            my = sum(c[1] for c in cells) / float(len(cells))
            defects.append({
                "x_mm": round((mx - (size - 1) / 2.0) * px, 2),
                "y_mm": round((my - (size - 1) / 2.0) * px, 2),
                "area_mm2": len(cells),
                "min_mu_per_cm": round(min(image[c[1]][c[0]] for c in cells), 3),
            })
    defects.sort(key=lambda item: -item["area_mm2"])

    return {
        "geometry": {"image_size": size, "pixel_size_mm": px},
        "roi": roi_rows,
        "defects": defects,
        "meta": {
            "n_angles": len(geo["angles_deg"]),
            "n_channels": geo["channels"],
            "repaired_channels": list(scan.get("detector", {}).get("bad_channels") or []),
        },
    }


def main():
    parser = argparse.ArgumentParser(description="工业 CT 截面重建与缺陷复核")
    parser.add_argument("--scan", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    scan = load_scan(args.scan)
    lines = to_line_integrals(scan)
    image = reconstruct(lines, scan)
    report = summarize(image, scan)

    out_dir = os.path.dirname(os.path.abspath(args.out))
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=2)
        f.write("\n")


if __name__ == "__main__":
    main()
