import subprocess


def run_cmd(cmd):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout_str, stderr_str = p.communicate()
    return_code = p.returncode
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd, stdout_str, stderr_str)
    return stdout_str, stderr_str
