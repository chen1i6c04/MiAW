import subprocess
from subprocess import CalledProcessError


class Option:
    def __init__(self, names, checker_function=None):
        self.names = names
        self.checker_function = checker_function
        self.value = None

    def __str__(self):
        return f"{self.names[0]} {self.value}"


class Switch:
    def __init__(self, names, checker_function=lambda x: isinstance(x, bool)):
        self.names = names
        self.checker_function = checker_function
        self.value = None

    def __str__(self):
        if self.value:
            return f"{self.names[0]}"
        else:
            return ""


class AbstractCommandline:
    def __init__(self, cmd, *args, **kwargs):
        self.program_name = cmd
        self.dispatch = args
        for key, value in kwargs.items():
            self.set_parameter(key, value)

    def _check_value(self, value, name, check_function):
        if check_function:
            is_good = check_function(value)
            if not is_good:
                raise ValueError(f"Invalid parameter value {value} for parameter {name}")

    def set_parameter(self, name, value=None):
        for parameter in self.parameters:
            if name in parameter.names:
                if not isinstance(parameter, Switch):
                    if value:
                        self._check_value(value, name, parameter.checker_function)
                        parameter.value = value
                else:
                    if value:
                        self._check_value(value, name, parameter.checker_function)
                        parameter.value = True

    def __str__(self):
        commandline = [self.program_name]
        for parameter in self.parameters:
            if parameter.value is not None:
                commandline.append(str(parameter))
        for parameter in self.dispatch:
            commandline.append(parameter)
        return ' '.join(commandline)

    def __call__(self):
        cmd = str(self)
        child_process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        stdout_str, stderr_str = child_process.communicate()
        return_code = child_process.returncode
        if return_code:
            raise CalledProcessError(return_code, str(self), stdout_str, stderr_str)
        return stdout_str, stderr_str
