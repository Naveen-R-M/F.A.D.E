# hpc/ssh_client.py
from __future__ import annotations
import os
import shlex
from dataclasses import dataclass
from typing import Optional, Dict, Tuple
import paramiko
from paramiko.proxy import ProxyCommand

def _load_pkey(path: str, passphrase: Optional[str]) -> paramiko.PKey:
    """
    Load an SSH private key file, trying common key types in order.
    """
    for key_cls in (paramiko.Ed25519Key, paramiko.RSAKey, paramiko.ECDSAKey):
        try:
            return key_cls.from_private_key_file(path, password=passphrase)
        except Exception:
            continue
    raise ValueError(f"Unsupported or unreadable SSH key at: {path}")

@dataclass
class SSHCredentials:
    host: str
    user: str
    port: int = 22
    key_path: Optional[str] = None       # ~/.ssh/id_ed25519 (recommended)
    password: Optional[str] = None       # only if key not used
    passphrase: Optional[str] = None     # if the private key is encrypted
    bastion_host: Optional[str] = None
    bastion_user: Optional[str] = None
    bastion_port: int = 22
    strict_host_key_checking: bool = True
    known_hosts_path: Optional[str] = None
    connect_timeout: float = 15.0
    keepalive_secs: int = 30

class HPCSSH:
    """
    Minimal SSH client for HPC interactions (exec + SFTP).
    Use as a context manager to ensure the connection closes cleanly.
    """
    def __init__(self, creds: SSHCredentials):
        self.creds = creds
        self.client: Optional[paramiko.SSHClient] = None

    def __enter__(self) -> "HPCSSH":
        self.client = self._connect()
        return self

    def __exit__(self, exc_type, exc, tb):
        try:
            if self.client:
                self.client.close()
        finally:
            self.client = None

    def _connect(self) -> paramiko.SSHClient:
        c = self.creds
        client = paramiko.SSHClient()

        # Host key policy
        if c.strict_host_key_checking:
            if c.known_hosts_path and os.path.exists(c.known_hosts_path):
                client.load_host_keys(c.known_hosts_path)
            client.load_system_host_keys()
            client.set_missing_host_key_policy(paramiko.RejectPolicy())
        else:
            # For first-run/dev only. Prefer verifying host keys!
            client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

        # Optional bastion/jump host using local OpenSSH ProxyCommand
        proxy = None
        if c.bastion_host:
            # Requires the local `ssh` binary on PATH
            bastion_user = c.bastion_user or c.user
            proxy_cmd = f"ssh -W {c.host}:{c.port} {bastion_user}@{c.bastion_host} -p {c.bastion_port}"
            proxy = ProxyCommand(proxy_cmd)

        # Auth preferences: key first, else password
        pkey = _load_pkey(c.key_path, c.passphrase) if c.key_path else None

        client.connect(
            hostname=c.host,
            username=c.user,
            port=c.port,
            pkey=pkey,
            password=(None if pkey else c.password),
            timeout=c.connect_timeout,
            sock=proxy,
            allow_agent=True,           # use ssh-agent/Pageant if available
            look_for_keys=bool(not pkey and not c.password)  # allow ~/.ssh keys
        )

        # Keepalive to prevent idle drops
        transport = client.get_transport()
        if transport and c.keepalive_secs:
            transport.set_keepalive(c.keepalive_secs)

        return client

    def run(
        self,
        command: str,
        *,
        cwd: Optional[str] = None,
        env: Optional[Dict[str, str]] = None,
        get_pty: bool = False,
        timeout: Optional[float] = None,
    ) -> Tuple[int, str, str]:
        """
        Execute a remote command. Returns (exit_code, stdout, stderr).
        """
        if not self.client:
            raise RuntimeError("SSH client is not connected")

        full_cmd = command
        if cwd:
            full_cmd = f"cd {shlex.quote(cwd)} && {command}"
        if env:
            exports = " ".join(f"{k}={shlex.quote(v)}" for k, v in env.items())
            full_cmd = f"{exports} {full_cmd}"

        stdin, stdout, stderr = self.client.exec_command(full_cmd, get_pty=get_pty, timeout=timeout)
        out = stdout.read().decode("utf-8", errors="replace")
        err = stderr.read().decode("utf-8", errors="replace")
        code = stdout.channel.recv_exit_status()
        return code, out, err

    def put(self, local_path: str, remote_path: str):
        """Upload a file via SFTP."""
        if not self.client:
            raise RuntimeError("SSH client is not connected")
        sftp = self.client.open_sftp()
        try:
            sftp.put(local_path, remote_path)
        finally:
            sftp.close()

    def get(self, remote_path: str, local_path: str):
        """Download a file via SFTP."""
        if not self.client:
            raise RuntimeError("SSH client is not connected")
        sftp = self.client.open_sftp()
        try:
            sftp.get(remote_path, local_path)
        finally:
            sftp.close()
