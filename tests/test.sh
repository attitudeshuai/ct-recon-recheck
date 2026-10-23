#!/bin/bash
# 评分脚本：pytest 跑 /tests 下的测试，权重汇总由 tests/conftest.py 写 reward。
# reward 只写 /logs/verifier/reward.txt；开头先覆盖为 0，崩溃/超时/信号都归零。
set -u

# 评分脚本：先把评分环境隔离（不让被测目录里同名模块劫持 import），再跑测试。
#   · cd / 且 python3 -P：不让工作目录进 sys.path，挡住 /app/pytest.py 这类同名劫持
#   · env -u PYTHON*：挡住 PYTHONPATH / sitecustomize 注入
#   · 解释器走绝对路径：挡住 PATH 上的同名 python3

mkdir -p /logs/verifier
rm -f /logs/verifier/pytest.xml /logs/verifier/reward.txt  # 先清旧报告，挡预置伪造的 junit
echo 0 > /logs/verifier/reward.txt
DONE=0
trap '[ "$DONE" = 1 ] || echo 0 > /logs/verifier/reward.txt' EXIT
trap 'exit 1' INT TERM HUP

cd /

# 判分工具不预装进镜像（B5）：pytest 随包放在 /tests/wheels，验证阶段离线安装
/usr/local/bin/python3 -m pip install --no-index --find-links /tests/wheels --quiet pytest > /logs/verifier/pip.log 2>&1 || true

timeout 240 env -u PYTHONPATH -u PYTHONHOME -u PYTHONSTARTUP /usr/local/bin/python3 -P -m pytest /tests -q -p no:cacheprovider -o addopts="" \
    --junitxml=/logs/verifier/pytest.xml > /logs/verifier/pytest.log 2>&1 || true

# 完成报告校验：没有结构化结果（pytest 没跑起来 / 被跳过）就按 INCOMPLETE 记 0
if [ ! -s /logs/verifier/pytest.xml ]; then
  echo 0 > /logs/verifier/reward.txt
fi

DONE=1
chmod 0444 /logs/verifier/reward.txt 2>/dev/null || true
