"""
SSH client for connecting to HPC cluster.

This module provides SSH connectivity for running remote commands
and transferring files to/from the HPC cluster.
"""

import logging
import paramiko
import os
from pathlib import Path
from typing import Optional, Tuple, List
import time

from fade.config import config
from fade.utils import get_logger

logger = get_logger("tools.ssh_client")


class HPCSSHClient:
    """SSH client for HPC cluster operations."""
    
    def __init__(self):
        """Initialize SSH client with credentials from config."""
        self.host = os.getenv("HPC_HOST", "explorer.northeastern.edu")
        self.port = int(os.getenv("HPC_PORT", "22"))
        self.username = os.getenv("HPC_USER")
        self.password = os.getenv("HPC_PASSWORD")
        self.ssh_key = os.getenv("HPC_SSH_KEY")
        self.passphrase = os.getenv("HPC_SSH_PASSPHRASE")
        
        self.client = None
        self.sftp = None
        
    def connect(self) -> bool:
        """
        Establish SSH connection to HPC cluster.
        
        Returns:
            True if connection successful
        """
        try:
            self.client = paramiko.SSHClient()
            self.client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            
            connect_kwargs = {
                "hostname": self.host,
                "port": self.port,
                "username": self.username,
                "timeout": 30
            }
            
            # Use SSH key if available, otherwise use password
            if self.ssh_key and Path(self.ssh_key).exists():
                key = paramiko.RSAKey.from_private_key_file(
                    self.ssh_key,
                    password=self.passphrase
                )
                connect_kwargs["pkey"] = key
                logger.info(f"Connecting to {self.host} with SSH key")
            else:
                connect_kwargs["password"] = self.password
                logger.info(f"Connecting to {self.host} with password")
            
            self.client.connect(**connect_kwargs)
            self.sftp = self.client.open_sftp()
            
            logger.info(f"Successfully connected to {self.host}")
            return True
            
        except Exception as e:
            logger.error(f"Failed to connect to HPC: {e}")
            return False
    
    def execute_command(self, command: str, timeout: int = 60) -> Tuple[str, str, int]:
        """
        Execute command on HPC cluster.
        
        Args:
            command: Command to execute
            timeout: Command timeout in seconds
            
        Returns:
            Tuple of (stdout, stderr, exit_code)
        """
        if not self.client:
            if not self.connect():
                return "", "Not connected to HPC", 1
        
        try:
            logger.debug(f"Executing: {command}")
            stdin, stdout, stderr = self.client.exec_command(
                command,
                timeout=timeout
            )
            
            # Wait for command to complete
            exit_status = stdout.channel.recv_exit_status()
            
            stdout_text = stdout.read().decode('utf-8')
            stderr_text = stderr.read().decode('utf-8')
            
            if exit_status != 0:
                logger.warning(f"Command exited with status {exit_status}")
                logger.debug(f"stderr: {stderr_text}")
            
            return stdout_text, stderr_text, exit_status
            
        except Exception as e:
            logger.error(f"Command execution failed: {e}")
            return "", str(e), 1
    
    def upload_file(self, local_path: str, remote_path: str) -> bool:
        """
        Upload file to HPC cluster.
        
        Args:
            local_path: Local file path
            remote_path: Remote destination path
            
        Returns:
            True if upload successful
        """
        if not self.sftp:
            if not self.connect():
                return False
        
        try:
            logger.info(f"Uploading {local_path} to {remote_path}")
            self.sftp.put(local_path, remote_path)
            return True
            
        except Exception as e:
            logger.error(f"File upload failed: {e}")
            return False
    
    def download_file(self, remote_path: str, local_path: str) -> bool:
        """
        Download file from HPC cluster.
        
        Args:
            remote_path: Remote file path
            local_path: Local destination path
            
        Returns:
            True if download successful
        """
        if not self.sftp:
            if not self.connect():
                return False
        
        try:
            logger.info(f"Downloading {remote_path} to {local_path}")
            self.sftp.get(remote_path, local_path)
            return True
            
        except Exception as e:
            logger.error(f"File download failed: {e}")
            return False
    
    def list_directory(self, remote_path: str) -> List[str]:
        """
        List contents of remote directory.
        
        Args:
            remote_path: Remote directory path
            
        Returns:
            List of filenames
        """
        if not self.sftp:
            if not self.connect():
                return []
        
        try:
            return self.sftp.listdir(remote_path)
        except Exception as e:
            logger.error(f"Failed to list directory: {e}")
            return []
    
    def mkdir_remote(self, remote_path: str) -> bool:
        """
        Create directory on remote server.
        
        Args:
            remote_path: Remote directory path
            
        Returns:
            True if successful
        """
        # Use mkdir -p to create parent directories if needed
        stdout, stderr, exit_code = self.execute_command(f"mkdir -p {remote_path}")
        return exit_code == 0
    
    def file_exists(self, remote_path: str) -> bool:
        """
        Check if remote file exists.
        
        Args:
            remote_path: Remote file path
            
        Returns:
            True if file exists
        """
        stdout, stderr, exit_code = self.execute_command(f"test -e {remote_path}")
        return exit_code == 0
    
    def close(self):
        """Close SSH connection."""
        try:
            if self.sftp:
                self.sftp.close()
            if self.client:
                self.client.close()
            logger.info("SSH connection closed")
        except:
            # Ignore errors during shutdown
            pass
    
    def __del__(self):
        """Cleanup on deletion."""
        try:
            self.close()
        except:
            # Ignore errors during shutdown
            pass


# Singleton instance
_hpc_client = None

def get_hpc_client() -> HPCSSHClient:
    """Get or create HPC SSH client singleton."""
    global _hpc_client
    if _hpc_client is None:
        _hpc_client = HPCSSHClient()
    return _hpc_client
