# app/hpc/test_connection.py
from __future__ import annotations
import os
from dotenv import load_dotenv

# Import the SSH client (works whether run as module or from the folder)
try:
    from .ssh_client import SSHCredentials, HPCSSH   # python -m app.hpc.test_connection
except Exception:
    from ssh_client import SSHCredentials, HPCSSH    # cd app/hpc && python test_connection.py

load_dotenv()  # read .env at project root if present

# --- HPC connection config ---
HPC_HOST = os.getenv("HPC_HOST", "")
HPC_PORT = int(os.getenv("HPC_PORT", "22"))
HPC_USER = os.getenv("HPC_USER", "")
HPC_SCRATCH_DIR = os.getenv("HPC_SCRATCH_DIR", "/scratch/rajagopalmohanraj.n")

SSH_PRIVATE_KEY_PATH = os.getenv("SSH_PRIVATE_KEY_PATH") or None
SSH_KEY_PASSPHRASE = os.getenv("SSH_KEY_PASSPHRASE") or None
HPC_PASSWORD = os.getenv("HPC_PASSWORD") or None

BASTION_HOST = os.getenv("BASTION_HOST") or None
BASTION_USER = os.getenv("BASTION_USER") or None
BASTION_PORT = int(os.getenv("BASTION_PORT", "22"))

# Script path (symlink or real file)
HPC_SCRIPT_PATH = os.getenv(
    "HPC_SCRIPT_PATH",
    "/scratch/rajagopalmohanraj.n/F.A.D.E/run_fade_nextflow.sh"
)

# Increase if your script runs the pipeline in foreground
HPC_SUBMIT_TIMEOUT = float(os.getenv("HPC_SUBMIT_TIMEOUT", "300"))

QUERY = "Find molecules targeting EGFR for lung cancer"

def build_creds():
    return SSHCredentials(
        host=HPC_HOST,
        user=HPC_USER,
        port=HPC_PORT,
        key_path=SSH_PRIVATE_KEY_PATH,
        passphrase=SSH_KEY_PASSPHRASE,
        password=HPC_PASSWORD,
        bastion_host=BASTION_HOST,
        bastion_user=BASTION_USER,
        bastion_port=BASTION_PORT,
        strict_host_key_checking=os.getenv("STRICT_HOST_KEY_CHECKING", "true").lower() == "true",
        known_hosts_path=os.getenv("KNOWN_HOSTS_PATH") or None,
        connect_timeout=float(os.getenv("HPC_CONNECT_TIMEOUT", "15")),
        keepalive_secs=int(os.getenv("HPC_KEEPALIVE_SECS", "30")),
    )

def main():
    creds = build_creds()

    # --- Preflight (resolve symlink, ensure executable, require sibling .env) ---
    preflight_cmd = r"""bash -lc '
        set -Eeuo pipefail
        tgt="$(readlink -f "$SCRIPT" 2>/dev/null || echo "$SCRIPT")"
        if [ ! -e "$tgt" ]; then
            echo "MISSING: $tgt"
            exit 2
        fi
        [ -x "$tgt" ] || chmod +x "$tgt" || true
        dir="$(dirname "$tgt")"
        echo "Preflight OK: script=$tgt env=$HPC_SCRATCH_DIR/F.A.D.E/nextflow/.env"
    '"""

    with HPCSSH(creds) as ssh:
        pcode, pout, perr = ssh.run(
            preflight_cmd,
            env={"SCRIPT": HPC_SCRIPT_PATH, "HPC_SCRATCH_DIR": HPC_SCRATCH_DIR},
            get_pty=True,
            timeout=60,
        )
        print("preflight exit:", pcode)
        print("preflight stdout:\n", pout)
        print("preflight stderr:\n", perr)
        if pcode != 0:
            return

    # --- Execute: cd to script dir and invoke with the query ---
    run_cmd = r"""bash -lc '
        set -Eeuo pipefail
        tgt="$(readlink -f "$SCRIPT" 2>/dev/null || echo "$SCRIPT")"
        dir="$(dirname "$tgt")"
        cd "$dir"
        "$tgt" "$ARG"
    '"""

    with HPCSSH(creds) as ssh:
        code, out, err = ssh.run(
            run_cmd,
            env={"SCRIPT": HPC_SCRIPT_PATH, "ARG": QUERY},
            get_pty=True,
            timeout=HPC_SUBMIT_TIMEOUT,
        )

    print("exit:", code)
    print("stdout:\n", out)
    print("stderr:\n", err)

if __name__ == "__main__":
    main()
