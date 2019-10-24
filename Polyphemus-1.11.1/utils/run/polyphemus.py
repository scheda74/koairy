# Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
#     Author(s): Vivien Mallet
#
# This file is part of the air quality modeling system Polyphemus.
#
# Polyphemus is developed in the INRIA - ENPC joint project-team CLIME and in
# the ENPC - EDF R&D joint laboratory CEREA.
#
# Polyphemus is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Polyphemus is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
#
# For more information, visit the Polyphemus web site:
#      http://cerea.enpc.fr/polyphemus/


"""\package polyphemus

This module provides facilities to launch Polyphemus programs.
"""


try:
    from network import Network
except:
    ## The module 'str' if the module 'network' is not available.
    Network = str


import sys


##############
# POLYPHEMUS #
##############


class Polyphemus:
    """This class manages the execution of several Polyphemus programs.
    """


    def __init__(self, net = None):
        """Initializes the network.

        \param net the network over which the simulations should be launched.
        """
        ## The list of programs.
        self.program_list = []

        ## A string.
        self.log_separator = "-" * 78 + "\n\n"

        ## The list of processes.
        self.process = []

        ## A network.Network instance.
        self.net = net


    def SetNetwork(self, net):
        """Sets the network.

        \param net the network over which the simulations should be launched.
        """
        self.net = net


    def AddProgram(self, program):
        """Adds a program.

        \param program the program to be added.
        """
        if isinstance(program, str):
            self.program_list.append(Program(program))
        else:
            self.program_list.append(program)
        def compare(x, y):
            if x.group < y.group:
                return -1
            elif x.group > y.group:
                return 1
            else:
                return 0
        self.program_list.sort(compare)


    def Clear(self):
        """Clears the program list.
        """
        self.program_list = []


    def Run(self, log = sys.stdout):
        """Executes the set of programs on localhost.
        """
        for program in self.program_list:
            print "Program name: ", program.name.split("/")[-1]
            program.Run(log)
            log.write(self.log_separator)
            if program.status != 0:
                raise Exception, "Program \"" + program.basename \
                      + "\" failed  (status " + str(program.status) + ")."


    def RunNetwork(self, log = sys.stdout, delay = 30):
        """Executes the set of programs on the network.

        \param delay the minimum period of time between the launch of two
        programs. Unit: seconds.
        """
        import time, commands
        self.process = []
        host = []
        beg_time = []
        ens_time = ["" for i in range(len(self.program_list))]
        # Index of the first program from current group.
        i_group = 0
        count_program = 1
        count_host = 0
        # List of available hosts (tuple: (hostname, Ncpu))
        host_available = self.net.GetAvailableHosts()
        Ncpu_list = [x[1] for x in host_available]
        # Sum of available cpus.
        Ncpu = sum(Ncpu_list)
        cpu_cumsum = [sum(Ncpu_list[0:x + 1]) for x in range(len(Ncpu_list))]
        # Copies and replaces for configuration files.
        for i in range(len(self.program_list)):
            program = self.program_list[i]
            program.config.Proceed()
        # Program runs on Network.
        for i in range(len(self.program_list)):
            program = self.program_list[i]
            if i > i_group and program.group != self.program_list[i-1].group:
                # If any process from the previous group is still up.
                while min([x.poll() for x in self.process[i_group:]]) == -1:
                    time.sleep(delay)
                i_group = i
                count_program = 1
                count_host = 0
                host_available = self.net.GetAvailableHosts()
                Ncpu_list = [x[1] for x in host_available]
                Ncpu = sum(Ncpu_list)
                cpu_cumsum = [sum(Ncpu_list[0:x + 1]) for x \
                              in range(len(Ncpu_list))]
            # If all hosts are busy.
            if count_program > Ncpu:
                time.sleep(70.)
                host_available = self.net.GetAvailableHosts()
                Ncpu_list = [x[1] for x in host_available]
                Ncpu = sum(Ncpu_list)
                cpu_cumsum = [sum(Ncpu_list[0:x + 1]) for x \
                              in range(len(Ncpu_list))]
                count_host = 0
                count_program = 1
                while Ncpu == 0:
                    time.sleep(60.)
                    host_available = self.net.GetAvailableHosts()
                    Ncpu_list = [x[1] for x in host_available]
                    Ncpu = sum(Ncpu_list)
                cpu_cumsum = [sum(Ncpu_list[0:x + 1]) for x \
                              in range(len(Ncpu_list))]

            # Changes host.
            if count_program > cpu_cumsum[count_host]:
                count_host += 1

            current_host = host_available[count_host][0]
            print "Program: ", program.name.split("/")[-1], \
                  " - Available host: ", current_host
            p = self.net.LaunchBG(program.Command(), host = current_host)
            self.process.append(p)
            count_program += 1
            host.append(current_host)
            beg_time.append(time.asctime())

            for j in range(i):
                if ens_time[j] == "" and self.process[i].poll() != 1:
                    ens_time[j] = time.asctime()

            # Checks process status.
            for j in range(i):
                if self.process[j].poll() != -1 and self.process[j].wait() != 0:
                    raise Exception, "The command: \"" + self.process[j].cmd \
                          + "\" does not work.\n" \
                          + "status: " + str(self.process[j].wait()) + ".\n" \
                          + "Error message: " \
                          + commands.getoutput("cat " +
                                               self.process[j].cmd.split()[-1])

        # Waits for the latest programs.
        while min([x.poll() for x in self.process[i_group:]]) == -1:
            time.sleep(delay)
            for j in range(len(self.program_list)):
                if ens_time[j] == "" and self.process[i].poll() != 1:
                    ens_time[j] = time.asctime()

        i_group = 0
        for i in range(len(self.program_list)):
            program = self.program_list[i]
            log_str = ""
            # New group ?
            if i > i_group and program.group != self.program_list[i-1].group:
                log_str += ("### GROUP " + str(program.group)
                            + " ###").center(78) + "\n" + self.log_separator
                i_group = i

            log_str += program.Command() + "\n" \
                       + "\nStatus: " + str(self.process[i].poll()) + "\n" \
                       + "Hostname: " + str(host[i]) + "\n" \
                       + "Started at " + str(beg_time[i]) + "\n" \
                       + "Ended approximatively at " + str(ens_time[i]) \
                       + "\n" + self.log_separator
            log.write(log_str)


    def Try(self, log = sys.stdout):
        """Performs a dry run.
        """
        for program in self.program_list:
            program.Try(log)
            log.Write(self.log_separator)
            if program.status != 0:
                raise Exception, "Program \"" + program.basename \
                      + "\" failed (status " + str(program.status) + ")."


###########
# PROGRAM #
###########


class Program:
    """This class manages a program associated with configuration files.
    """


    def __init__(self, name = None, config = None, arguments_format = " %a", group = 0):
        """Full initialization.

        \param name the program name.
        \param config the program configuration.
        \param arguments_format the arguments_format of the arguments, where "%a" is replaced with
        the configuration files.
        \param group the group index.
        """
        if config is not None:
            ## A configuration file or a polyphemus.Configuration instance.
            self.config = config
        else:
            self.config = Configuration()

        ## The path of the program.
        self.name = name

        import os

        ## The basename of the program.
        self.basename = os.path.basename(name)

        ## A string.
        self.exec_path = "./"

        ## The format of the arguments.
        self.arguments_format = arguments_format

        ## A string.
        self.priority = "0"

        ## None.
        self.output_directory = None

        ## The group index.
        self.group = group

        ## The status of the program.
        self.status = None

    def Run(self, log = sys.stdout):
        """Executes the program.
        """
        self.config.Proceed()
        from subprocess import Popen, PIPE
        p = Popen([self.name, self.arguments_format], stdout=PIPE, bufsize=1)
        for line in iter(p.stdout.readline, b''):
            if line is not None:
                log.write(line)
        p.communicate()
        self.status = p.returncode

    def Command(self):
        """Returns the command to launch the program.

        The program must be ready to be launched.
        """
        if not self.IsReady():
            raise Exception, "Program \"" + self.name + "\" is not ready."
        arguments_format = self.arguments_format[:]
        command = "nice time " + self.name \
                  + arguments_format.replace("%a", self.config.GetArgument())
        return command


    def SetConfiguration(self, config, mode = "random", path = None,
                         replacement = None, additional_file_list = []):
        """Sets the program configuration files.

        \param config the configuration files associated with the program.
        \param mode the copy mode. Mode "raw" just copies to the target path,
        while mode "random" appends a random string at the end of the file
        name. Mode "random_path" appends a random directory in the path. This
        entry is useless if "config" is a Configuration instance.
        \param path the path where configuration files should be copied. If
        set to None, then the temporary directory "/tmp" is used. This entry
        is useless if "config" is a Configuration instance.
        \param replacement the map of replaced strings and the replacement
        values. This entry is useless if "config" is a Configuration
        instance.
        \param additional_file_list an additional configuration file or a
        list of additional configuration files to be managed. Just like
        primary configuration files, they are subject to the replacements and
        copies, but are not considered as program arguments.
        """
        if isinstance(config, list) or isinstance(config, str):
            if mode is None:
                mode = "tmp"
            self.config = Configuration(config, mode, path,
                                        additional_file_list)
            if replacement is not None:
                self.config.SetReplacementMap(replacement)
        else:
            self.config = config


    def Try(self, log = sys.stdout):
        """Performs a dry run.
        """
        self.config.Proceed()
        import os
        if os.path.isfile(self.name):
            self.status = 0
        else:
            self.status = 1
        arguments_format = self.arguments_format[:]
        command = self.name + arguments_format.replace("%a", self.config.GetArgument())
        log.write( "Running program \"" + self.basename + "\":\n" \
                   + "   " + command )


    def IsReady(self):
        """Checks whether the program can be launched.

        @return True if the program can be executed, False otherwise.
        """
        return self.config.IsReady()


#################
# CONFIGURATION #
#################


class Configuration:
    """This class manages configuration files.

    It proceeds replacements in the files and makes copies of the files.
    """


    def __init__(self, file_list = [], mode = "random", path = None,
                 additional_file_list = []):
        """Initialization of configuration information.

        \param file_list the configuration file or the list of configuration
        files to be managed.
        \param mode the copy mode. Mode "raw" just copies to the target path,
        while mode "random" appends a random string at the end of the file
        name. Mode "random_path" appends a random directory in the path.
        \param path the path where configuration files should be copied. If
        set to None, then the temporary directory "/tmp" is used.
        \param additional_file_list an additional configuration file or a
        list of additional configuration files to be managed. Just like
        primary configuration files, they are subject to the replacements and
        copies, but are not considered as program arguments.
        """
        import os
        if isinstance(file_list, str):
            ## The configuration file or list of configuration to be managed.
            self.raw_file_list = [file_list]
        else:
            self.raw_file_list = file_list

        ## The number of configuration files.
        self.Narg = len(self.raw_file_list)

        if isinstance(additional_file_list, str):
            self.raw_file_list.append(additional_file_list)
        else:
            self.raw_file_list += additional_file_list
        for f in self.raw_file_list:
            if not os.path.isfile(f):
                raise Exception, "Unable to find \"" + f + "\"."

        ## The list of replaced configuration files.
        self.file_list = []

        ## Are the configuration files ready to use?
        self.ready = False

        self.SetMode(mode)
        self.SetPath(path)

        ## The map of replaced strings and the replacement values.
        self.config = {}


    def SetMode(self, mode = "random"):
        """Sets the copy mode.
        \param mode the copy mode. Mode "raw" just copies to the target path,
        while mode "random" appends a random string at the end of the file
        name. Mode "random_path" appends a random directory in the path.
        """
        if mode not in ["random", "random_path", "raw"]:
            raise Exception, "Mode \"" + str(mode) + "\" is not supported."

        ## The copy mode.
        self.mode = mode


    def SetPath(self, path):
        """Sets the path.
        \param path the path where configuration fields should be copied. If
        set to None, then the temporary directory "/tmp" is used.
        """
        if path is not None:
            ## The path where configuration files should be copied.
            self.path = path
        else:
            self.path = "/tmp/"


    def IsReady(self):
        """Tests whether the configuration files are ready for use.

        @return True if the configuration files are ready for use, False
        otherwise.
        """
        return self.ready or self.raw_file_list == []


    def GetReplacementMap(self):
        """Returns the map of replaced strings and the replacement values.

        @return the map of replaced strings and the replacement values.
        """
        return self.config


    def SetReplacementMap(self, config):
        """Sets the map of replaced strings and the replacement values.

        \param config the map of replaced strings and the replacement
        values.
        """
        self.config = config


    def SetConfiguration(self, config, mode = "random", path = None):
        """Initialization of configuration information, except file names.

        \param config the map of replaced strings and the replacement
        values.
        \param mode the copy mode. Mode "raw" just copies to the target path,
        while mode "random" appends a random string at the end of the file
        name. Mode "random_path" appends a random directory in the path.
        \param path the path where configuration fields should be copied. If
        set to None, then the temporary directory "/tmp" is used.
        """
        self.SetMode(mode)
        self.SetPath(path)
        self.SetReplacementMap(config)
        self.Proceed()


    def Proceed(self):
        """Proceeds replacement in configuration files and copy them.
        """
        import os, shutil, fileinput
        self.file_list = []
        if self.mode == "random_path" and self.raw_file_list is not []:
            import tempfile
            random_path = tempfile.mkdtemp(prefix = self.path)
        for f in self.raw_file_list:
            if self.mode == "raw":
                if os.path.dirname(f) == self.path:
                    raise Exception, "Error: attempt to overwrite" \
                          + " the raw configuration file \"" + f + "\"."
                name = os.path.join(self.path, os.path.basename(f))
                shutil.copy(f, name)
            elif self.mode == "random":
                name = os.path.join(self.path, os.path.basename(f))
                _, name = tempfile.mkstemp(prefix = name + "-")
                shutil.copy(f, name)
            elif self.mode == "random_path":
                name = os.path.join(random_path, os.path.basename(f))
                shutil.copy(f, name)
            self.file_list.append(name)

        if self.file_list != []:
            for line in fileinput.input(self.file_list, 1):
                new_line = line
                for i in self.config.keys():
                    new_line = new_line.replace(str(i), str(self.config[i]))
                if self.mode == "random_path":
                    new_line = new_line.replace("%random_path%", random_path)
                print new_line,
            fileinput.close()
        self.ready = True


    def GetRawFileList(self):
        """Returns the list of reference (or raw) configuration files.

        @return the list of reference (or raw) configuration files.
        """
        return self.raw_file_list


    def SetRawFileList(self, file_list):
        """Sets the list of reference (or raw) configuration files.

        \param file_list the list of reference (or raw) configuration files.
        """
        self.raw_file_list = file_list
        self.ready = False


    def Clear(self):
        """Clears all, including configuration file names.
        """
        self.raw_file_list = []
        self.file_list = []
        self.ready = False
        self.mode = "random"
        self.path = "/tmp/"
        self.config = {}


    def GetArgument(self):
        """Returns the list of program arguments.

        @return the list of program arguments aggregated in a string (and
        split by an empty space).
        """
        if self.IsReady():
            return " ".join(self.file_list[:self.Narg])
        else:
            raise Exception, "Not ready."
