#!/bin/bash
# 评分脚本：pytest 跑 /tests 下的测试，权重汇总由 tests/conftest.py 写 reward。
# reward 只写 /logs/verifier/reward.txt；开头先覆盖为 0，崩溃/超时/信号都归零。
set -u

mkdir -p /logs/verifier
echo 0 > /logs/verifier/reward.txt
DONE=0
trap '[ "$DONE" = 1 ] || echo 0 > /logs/verifier/reward.txt' EXIT
trap 'exit 1' INT TERM HUP

cd /app

# 判分工具不预装进镜像（B5）：pytest 随包放在 /tests/wheels，验证阶段离线安装
pip install --no-index --find-links /tests/wheels --quiet pytest > /logs/verifier/pip.log 2>&1 || true

timeout 240 python3 -m pytest /tests -q -p no:cacheprovider -o addopts="" \
    --junitxml=/logs/verifier/pytest.xml > /logs/verifier/pytest.log 2>&1 || true

# 完成报告校验：没有结构化结果（pytest 没跑起来 / 被跳过）就按 INCOMPLETE 记 0
if [ ! -s /logs/verifier/pytest.xml ]; then
  echo 0 > /logs/verifier/reward.txt
fi

DONE=1
chmod 0444 /logs/verifier/reward.txt 2>/dev/null || true
