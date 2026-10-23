# -*- coding: utf-8 -*-
"""pytest 汇总钩子：按权重把测试结果写成 reward。

口径：每条测试占一个权重，权重和 = 1（见 WEIGHTS）；只有 PASS 才拿到该条权重；
不在 WEIGHTS 里的测试名一律不给分（防止误加测试把总分抬高）。
reward 写在 /logs/verifier/reward.txt（本地调试可用 CT_REWARD_PATH 覆盖）。
"""
import os

WEIGHTS = {
    "test_roi_mu_matches_phantom_on_base_scan": 0.13,
    "test_roi_mu_matches_second_phantom_geometry": 0.10,
    "test_low_contrast_insert_resolved": 0.12,
    "test_dark_current_large_offset_recovered": 0.09,
    "test_dead_channel_recovery_with_unlisted_zero": 0.09,
    "test_defect_area_and_center_match_phantom": 0.11,
    "test_no_spurious_defect_outside_inspection_radius": 0.06,
    "test_repaired_channels_listed_exactly": 0.06,
    "test_invalid_scan_rejected_without_output": 0.09,
    "test_report_schema_and_rounding": 0.06,
    "test_rerun_produces_identical_bytes": 0.05,
    "test_failed_run_keeps_previous_report": 0.04,
}

_RESULTS = {}


def pytest_runtest_logreport(report):
    if report.when in ("setup", "call") and report.failed:
        _RESULTS[report.nodeid.split("::")[-1]] = False
    elif report.when == "call" and report.passed:
        _RESULTS.setdefault(report.nodeid.split("::")[-1], True)


def pytest_sessionfinish(session, exitstatus):
    # 完成报告校验：必须有每一条测试的结果，缺任何一条都按 INCOMPLETE 记 0
    missing = [n for n in WEIGHTS if n not in _RESULTS]
    total = 0.0
    if missing:
        print("[INCOMPLETE] 这些测试没有跑完：%s；reward 记 0" % ", ".join(missing))
    else:
        for name, weight in WEIGHTS.items():
            if _RESULTS.get(name) is True:
                total += weight
    path = os.environ.get("CT_REWARD_PATH", "/logs/verifier/reward.txt")
    out_dir = os.path.dirname(os.path.abspath(path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("%.6f\n" % total)
