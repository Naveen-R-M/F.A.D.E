# app/hpc/runner.py
from __future__ import annotations
import os, shlex, time
from typing import Dict, Any
from dotenv import load_dotenv
from .ssh_client import SSHCredentials, HPCSSH

load_dotenv()

# ------- SSH / HPC config -------
HPC_HOST       = os.getenv("HPC_HOST", "")
HPC_PORT       = int(os.getenv("HPC_PORT", "22"))
HPC_USER       = os.getenv("HPC_USER", "")
SSH_KEY        = os.getenv("SSH_PRIVATE_KEY_PATH") or None
SSH_PASSPHRASE = os.getenv("SSH_KEY_PASSPHRASE") or None
HPC_PASSWORD   = os.getenv("HPC_PASSWORD") or None

# Where your unchanged script lives
HPC_FADE_BASE  = os.getenv("HPC_FADE_BASE", "/scratch/rajagopalmohanraj.n/F.A.D.E")

# *** NEW: Fixed log root (can be overridden by .env: HPC_LOG_ROOT=...) ***
HPC_LOG_ROOT   = os.getenv(
    "HPC_LOG_ROOT",
    "/scratch/rajagopalmohanraj.n/F.A.D.E/backend/nextflow/logs"
)

# Optional: enable bash -x (trace) when invoking the script
HPC_TRACE      = os.getenv("HPC_TRACE", "0").lower() in ("1", "true", "yes")

def _bash_login(cmd: str) -> str:
    """Run remote commands in a login shell so env/modules are available."""
    return f"/bin/bash -lc {shlex.quote(cmd)}"

def _ts() -> str:
    return time.strftime("%Y%m%d_%H%M%S")

def submit_fade_job(job_id: str, prompt: str) -> Dict[str, Any]:
    """
    Fire-and-forget background submission. Writes logs under HPC_LOG_ROOT.
    Does NOT modify your run_fade_nextflow.sh; passes the prompt as a positional arg.
    """
    creds = SSHCredentials(
        host=HPC_HOST, user=HPC_USER, port=HPC_PORT,
        key_path=SSH_KEY, passphrase=SSH_PASSPHRASE, password=HPC_PASSWORD,
        strict_host_key_checking=True,
    )

    # Per-job log dir under the requested scratch path
    run_tag   = f"{_ts()}_{job_id[:8]}"
    run_dir   = f"{HPC_LOG_ROOT.rstrip('/')}/{run_tag}"
    script    = f"{HPC_FADE_BASE.rstrip('/')}/backend/run_fade_nextflow.sh"

    submitlog = f"{run_dir}/submit.log"
    runlog    = f"{run_dir}/run.log"
    pidfile   = f"{run_dir}/run.pid"

    # Persist the request for provenance
    req_json = (
        "{\\n"
        f"  \\\"job_id\\\": \\\"{job_id}\\\",\\n"
        f"  \\\"prompt\\\": \\\"{prompt.replace('\"','\\\\\\\"')}\\\",\\n"
        f"  \\\"created_at\\\": \\\"{_ts()}\\\"\\n"
        "}"
    )

    quoted_script = shlex.quote(script)
    quoted_prompt = shlex.quote(prompt)
    # Background start (nohup) with optional -x tracing
    start_line = (
        f"nohup {quoted_script} {quoted_prompt} >> {shlex.quote(runlog)} 2>&1 & echo $! > {shlex.quote(pidfile)}"
        if not HPC_TRACE else
        f"nohup bash -x {quoted_script} {quoted_prompt} >> {shlex.quote(runlog)} 2>&1 & echo $! > {shlex.quote(pidfile)}"
    )

    remote = f"""
set -euo pipefail

# Ensure log root and per-job dir
mkdir -p {shlex.quote(run_dir)}
echo "[CHECK] Host: $(hostname)"                         |& tee -a {shlex.quote(submitlog)}
echo "[CHECK] User: $(whoami)"                           |& tee -a {shlex.quote(submitlog)}
echo "[CHECK] CWD:  $(pwd)"                              |& tee -a {shlex.quote(submitlog)}
echo "[CHECK] PATH: $PATH"                               |& tee -a {shlex.quote(submitlog)}
echo "[INFO] Log dir: {run_dir}"                         |& tee -a {shlex.quote(submitlog)}

# Verify script location
if [ ! -f {shlex.quote(script)} ]; then
  echo "[ERR] Script missing: {script}"                  |& tee -a {shlex.quote(submitlog)}
  echo "[HINT] Check HPC_FADE_BASE or script path"       |& tee -a {shlex.quote(submitlog)}
  exit 44
fi
chmod +x {shlex.quote(script)} || true
echo "[OK] Script: {script}"                             |& tee -a {shlex.quote(submitlog)}
ls -l {shlex.quote(script)}                              |& tee -a {shlex.quote(submitlog)} || true

# Save request.json
cat > {shlex.quote(run_dir)}/request.json <<'REQ'
{req_json}
REQ
echo "[OK] wrote request.json"                           |& tee -a {shlex.quote(submitlog)}

# Announce command (prompt redacted in echo)
echo "[RUN] {script} '<omitted>' (background via nohup)" |& tee -a {shlex.quote(submitlog)}

# Start in background and capture PID
cd {shlex.quote(run_dir)}
{start_line}
echo "[OK] started; PID=$(cat {shlex.quote(pidfile)})"   |& tee -a {shlex.quote(submitlog)}
"""

    with HPCSSH(creds) as ssh:
        code, out, err = ssh.run(_bash_login(remote), get_pty=True, timeout=60)
        ok = (code == 0)
        return {
            "ok": ok,
            "job_id": job_id,
            "cwd": run_dir,  # <- the log directory you asked for
            "command": f"{script} '<omitted>' (background via nohup)",
            "exit_code": code,
            "stdout": out[-8000:],
            "stderr": err[-8000:],
        }
