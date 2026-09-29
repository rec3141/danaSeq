#!/usr/bin/env bash
# Print the training device for a PyTorch binner: "cuda" or "cpu".
#   torch_device.sh <env suffix> <auto|true|false>
# The env suffix names the tool's env (vamb -> dana-mag-vamb). The check runs
# that env's own Python, since the torch build (CPU or CUDA) differs per env
# and the python3 on PATH may have no torch at all. With "true", a missing GPU
# is reported loudly and the tool falls back to CPU instead of failing.
set -uo pipefail
env_name="dana-mag-$1"; want="${2:-auto}"

case "$want" in
    false|off|cpu) echo cpu; exit 0 ;;
    true|on|auto)  ;;
    *) echo "[WARNING] torch_device.sh: --gpu must be auto, true or false (got '$want'); using auto" >&2 ;;
esac

py=""
for cand in "${CONDA_PREFIX:-}" "/opt/conda/envs/$env_name"; do
    [[ -n "$cand" && -x "$cand/bin/python" && "$(basename "$cand")" == "$env_name" ]] && { py="$cand/bin/python"; break; }
done
[[ -z "$py" ]] && py=$(command -v python3 || true)

status=$("$py" - <<'PY' 2>&1
import torch
if torch.version.cuda is None:
    print("cpu-build")
elif not torch.cuda.is_available():
    print("no-device cuda-" + torch.version.cuda)
else:
    print("cuda " + torch.cuda.get_device_name(0))
PY
) || status="no-torch"

case "$status" in
    cuda\ *)
        echo "[INFO] $env_name: training on GPU (${status#cuda })" >&2
        echo cuda ;;
    *)
        if [[ "$want" == true || "$want" == on ]]; then
            echo "[WARNING] $env_name: GPU requested but not usable ($status via $py); training on CPU" >&2
        else
            echo "[INFO] $env_name: training on CPU ($status)" >&2
        fi
        echo cpu ;;
esac
