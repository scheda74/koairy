# Copyright (C) 2005-2009 INRIA
# Authors: Vivien Mallet, Damien Garaud
#
# This file is part of the air quality modeling system Polyphemus. It provides
# facilities to deal with computations over a Linux network.
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

"""\package network

Provides class 'Network' designed to managed processes launched over the
network.

The tools 'dsh' must be installed on your computer. It can be downloaded in
this website:
http://www.netfort.gr.jp/~dancer/software/downloads/list.cgi

You have to create a file with the list of hosts in this directory: \code
~/.dsh/group/your_file \endcode
"""

import os, types, commands, socket, popen2

## The ssh command.
ssh = 'ssh '


###########################
# MISCELLANEOUS FUNCTIONS #
###########################


def IsNum(str_value):
    """Is the string a float?
    \param str_value a string
    @return True if 'str_value' is a float, False otherwise.
    """
    is_num = True
    try:
        float(str_value)
    except ValueError:
        is_num = False
    return is_num


def IsInt(str_value):
    """Tests whether a string is an integer.
    \param str_value a string to be tested.
    @return True if 'str_value' is a integer, False otherwise.
    """
    is_num = True
    try:
        int(str_value)
    except ValueError:
        is_num = False
    return is_num


def ToNum(str_value):
    """Converts a string to a number.
    \param str_value a string to be converted.
    @return The number represented by 'str_value'.
    """
    if IsInt(str_value):
        return int(str_value)
    elif IsNum(str_value):
        return float(str_value)
    else:
        raise Exception, "\"" + str_value + "\" is not a number."


def remove(files):
    """Removes one or more files and/or directories."""
    if isinstance(files, str):
        files = [files]
    if not isinstance(files, list):
        raise ValueError
    for path in files:
        if os.path.isdir(path):
            import shutil
            shutil.rmtree(path)
        elif os.path.isfile(path):
            os.remove(path)


########
# HOST #
########


class Host:
    """Dedicated to host management."""


    def __init__(self, host = socket.gethostname()):
        """Initializes the name of the host and the number of CPUs.

        \param host the name of the host.
        """
        command_name = "cat /proc/cpuinfo | grep ^processor | wc -l"
        if type(host) == types.TupleType or type(host) == types.ListType:
            ## The name of the host.
            self.name = str(host[0])
            if len(host) == 1:
                (_, o) = commands.getstatusoutput(ssh + self.name
                                                  + ' ' + command_name)
                ## Number of CPU.
                self.cpu = ToNum(o)
            else:
                self.cpu = int(host[1])
        else:
            self.name = str(host)
            (_, o) = commands.getstatusoutput(ssh + self.name
                                              + ' ' + command_name)
            self.cpu = ToNum(o)


    def LoadAverage(self):
        """Returns the load average.

        @return a tuple of three load averages.
        """
        if self.name == socket.gethostname():
            return os.getloadavg()
        else:
            command = ssh + self.name + ' uptime'
            (s, o) = commands.getstatusoutput(command)
            if s == 0:
                return tuple([ToNum(elt.strip(',')) for elt in \
                                  o.split()[-3:]])
            else:
                return (9999, 9999, 9999)


    def LoadAverageFile(self):
        """Returns the load average in a file.

        @return the file where are the load averages.
        """
        import warnings
        warnings.simplefilter("ignore")
        tmp_file = os.tempnam("/tmp/")
        load_average = self.LoadAverage()
        sfile = open(tmp_file, 'w')
        sfile.write(' '.join([str(x) for x in load_average]))
        sfile.close()
        return tmp_file


    def LaunchInt(self, command):
        """Launches a command in interactive mode (using os.system).

        \param command the name of the command.
        @return the status of the command.
        """
        return os.system(command)


    def LaunchFG(self, command):
        """Launches a command in the foreground.

        \param command the name of the command.
        @return the output and the status of the command.
        """
        return commands.getstatusoutput(command)


    def LaunchBG(self, command):
        """Launches a command in the background and returns a Popen4 object.

        \param command the name of the command.
        @return a Popen4 object.
        """
        # Output is redirected.
        command = "( " + command + "; ) &> /dev/null"
        return popen2.Popen4(command)


    def LaunchWait(self, command, ltime, wait = 0.1):
        """Launches a command in the foreground and waits for its output for a
        given time after which the process is killed.

        \param command the name of the command.
        \param ltime the limit time.
        \param wait the waiting time.
        @return the output and the status of the command.
        """
        import time
        # Output is redirected.
        file_index = 0
        while os.path.exists("/tmp/output-" + str(file_index)):
            file_index += 1
        file_tmp = open("/tmp/output-" + str(file_index), 'w')
        file_tmp.close()
        file_name = "/tmp/output-" + str(file_index)
        os.chmod(file_name, 0600)
        command = "( " + command + "; ) &> " + file_name
        p = popen2.Popen4(command)
        t = time.time()
        while (p.poll() == -1 and time.time() - t < ltime):
            time.sleep(wait)
        s = p.poll()
        if s == -1:
            try:
                os.kill(p.pid, 9)
            except:
                pass
        o = open("/tmp/output-" + str(file_index), 'r').read()
        os.remove(file_name)
        return (s, o)


###########
# NETWORK #
###########


class Network:
    """Dedicated to network management."""


    def __init__(self, group = 'all'):
        """Initiliazes the list of hosts.

        \param group the file name where there is a list of hosts.
        """
        dsh_dir = os.environ['HOME'] + '/.dsh/'
        ## A list of Host instances.
        self.hosts = []
        # Checking files and directories.
        if not os.path.isdir(dsh_dir):
            raise Exception, "You must install 'dsh' and create " \
                + "a directory '.dsh/group' in your home."
        if not os.path.isfile(dsh_dir + 'group/' + group):
            raise Exception, "You must create a file with the list "\
                + "of your hosts. Unable to find \"" \
                + dsh_dir + "group/" + group + "\"."
        command = "dsh -g " + group  +" -r " + ssh + " -Mc " \
            +"'cat /proc/cpuinfo | grep ^processor | wc -l'"
        s, o = commands.getstatusoutput(command)
        # If the command does not work.
        if s != 0:
            raise Exception, "The command \"" + command + "\" did not work. " \
                + "The command returns:\n\n" \
                + o
        o = o.split()
        # A loop to store the name and cpu hosts.
        for i in range(0, len(o), 2):
            self.hosts.append(Host((o[i].strip(':'), o[i+1].strip(':'))))
        ## The name of the file where there is the hosts list.
        self.group = group


    def PrintHostNames(self):
        """Prints the name of hosts.
        """
        for host in self.hosts: print host.name


    def GetHostNames(self):
        """Returns the hosts names.
        @return the list of hosts.
        """
        l = []
        for host in self.hosts: l.append(host.name)
        return l


    def LaunchInt(self, command, host = socket.gethostname()):
        """Launches a command in interactive mode (using os.system).

        \param command the name of the command.
        \param host the name of the host.
        @return the status of the command.
        """
        if not isinstance(host, Host):
            host = Host(host)
        if host.name == socket.gethostname():
            return os.system(command)
        else:
            return os.system(ssh + host.name + " \"" + command + "\"")


    def LaunchFG(self, command, host = socket.gethostname()):
        """Launches a command in the foreground.

        \param command the name of the command.
        \param host the name of the host.
        @return the output and the status of the command.
        """
        if not isinstance(host, Host):
            host = Host(host)
        if host.name == socket.gethostname():
            (s, o) = commands.getstatusoutput(command)
            # # "commands.getstatusoutput" removes the last '\n' if it exists!
            # # Assuming that there was a line break:
            # o += "\n"
        else:
            (s, o) = commands.getstatusoutput(ssh + host.name \
                                              + " \"" + command + "\"")
        # Note: "commands.getstatusoutput" removes the last '\n' if it exists!
        return (s, o)


    def LaunchWait(self, command, ltime, wait = 0.1,
                   host = socket.gethostname()):
        """Launches a command in the foreground and waits for its output for a given
        time after which the process is killed.

        \param command the name of the command.
        \param ltime the limit time.
        \param wait the waiting time.
        \param host the name of the host.
        @return the output and the status of the command.
        """
        import time
        if not isinstance(host, Host):
            host = Host(host)
        # Output is redirected.
        file_index = 0
        while os.path.exists("/tmp/output-" + str(file_index)):
            file_index += 1
        file_tmp = open("/tmp/output-" + str(file_index), 'w')
        file_tmp.close()
        file_name = "/tmp/output-" + str(file_index)
        os.chmod(file_name, 0600)
        if host.name == socket.gethostname():
            command = "( " + command + "; ) &> " + file_name
        else:
            command = "( " + ssh + " " + host.name + " \"" \
                      + command + "\"; ) &> " + file_name
        p = popen2.Popen4(command)
        t = time.time()
        while (p.poll() == -1 and time.time() - t < ltime):
            time.sleep(wait)
        s = p.poll()
        if s == -1:
            try:
                os.kill(p.pid, 9)
            except:
                pass
        o = open("/tmp/output-" + str(file_index), 'r').read()
        os.remove(file_name)
        return (s, o)


    def LaunchBGTest(self):
        """Checks the LaunchBG function on all hosts."""
        command = "/bin/ls"
        proc = []
        for h in self.hosts:
            proc.append(self.LaunchBG(command, h))
        for p in proc:
            if p.poll() != 0 and p.poll() != -1:
                raise Exception, "An error appears during the launching of " \
                      + "the command \"" + command + "\" on host \"" \
                      + h.name + "\"."


    def LaunchBG(self, command, host = socket.gethostname()):
        """Launches a command in the background and returns a Popen4 object.

        \param command the name of the command.
        \param host the name of the host.
        @return a Popen4 object.
        """
        if not isinstance(host, Host):
            host = Host(host)
        # Output is redirected.
        command = "( " + command + "; ) &> /dev/null"
        if host.name == socket.gethostname():
            return popen2.Popen4(command)
        else:
            return popen2.Popen4(ssh + host.name + " \"" + command + "\"")


    def LaunchScreen(self, command, host = socket.gethostname()):
        """Launches a command in a screen.

        \param command the name of the command.
        \param host the name of the host.
        """
        if not isinstance(host, Host):
            host = Host(host)
        if host.name != socket.gethostname():
            command = ssh + host.name + " screen -d -m " + command
        else:
            command = "screen -d -m " + command
        os.system(command)


    def SendMail(self, subject, fromaddr, toaddr, msg = ""):
        """Sends a mail.
        """
        from email.MIMEText import MIMEText
        msg = MIMEText(msg)
        msg['Subject'] = subject
        msg['From'] = fromaddr
        msg['To'] = toaddr
        import smtplib
        server = smtplib.SMTP('localhost')
        server.sendmail(fromaddr, toaddr, msg.as_string())
        server.quit()


    def SendMailAttach(self, subject, attachments, fromaddr, toaddr,
                       msg = ""):
        """Sends a mail with attachments.
        """
        from email.MIMEMultipart import MIMEMultipart
        from email.MIMEText import MIMEText
        tmp = MIMEMultipart()
        tmp.attach(MIMEText(msg))
        msg = tmp
        msg['Subject'] = subject
        msg['From'] = fromaddr
        msg['To'] = toaddr
        if type(attachments) == types.StringType:
            attached_file = open(attachments, 'rb')
            msg.attach(MIMEText(attached_file.read()))
            attached_file.close()
        else:
            for attachment in attachments:
                attached_file = open(attachment, 'rb')
                msg.attach(MIMEText(attached_file.read()))
                attached_file.close()
        import smtplib
        server = smtplib.SMTP('localhost')
        server.sendmail(fromaddr, toaddr, msg.as_string())
        server.quit()


    def GetLoadAverages(self):
        """Returns the load averages.

        @return a list of hosts with the load averages for the past 1 minute.
        """
        result = []
        command = "dsh -r " + ssh + " -Mc -g "+ self.group \
            + " 'echo `hostname` `uptime`'"
        status, out = commands.getstatusoutput(command)
        # If the command does not work.
        if status != 0:
            raise Exception, "Unable to launch the command: \"" + command \
                  + "\"."
        # A loop for every hosts to pick the load averages.
        for line in out.split("\n"):
            name = line.split()[0].strip(':')
            load_average = ToNum(line.split()[-3].strip(','))
            result.append((name, load_average))
        return result


    def GetAvailableHosts(self):
        """Returns the available hosts.

        @return a list of hosts with the number of available cpu.
        """
        # The network load averages.
        host_avr = self.GetLoadAverages()
        result = []
        # A loop for every hosts to pick up the number of available cpu.
        for host in self.hosts:
            index = [x[0] for x in host_avr].index(host.name)
            avr = host_avr[index][1]
            if avr < host.cpu - 0.5:
                result.append((host.name, int(host.cpu - avr + 0.5)))
        return result


    def GetAvailableHost(self, load_limit = 0.3):
        """Return an available host.

        \param load_limit the limit time to load.
        @return an available host.
        """
        result = ""
        host_list = self.GetLoadAverages()
        for host in host_list:
            if host[1] <= load_limit:
                result = host[0]
        if len(result) == 0:
            load_list = [x[1] for x in host_list]
            load_min = min(load_list)
            result = host[load_list.index(load_min)]
        return result


    def GetRandomAvailableHost(self, load_limit = 0.3):
        """Returns a random available host.

        \param load_limit the limit time to load.
        @return the name of a random available host.
        """
        host_available = []
        for host in self.hosts:
            l = host.LoadAverage()[0]
            if l < load_limit:
                host_available.append(host)
        if len(host_available) != 0:
            import random
            return host_available[random.randint(0, len(host_available) - 1)].name
        else:
            return self.hosts[0].name


if __name__ == "__main__":
    ## A network.Network instance.
    net = Network()
    print "# Name \t # CPU"
    for host in net.hosts:
        print(host.name + '\t' + str(host.cpu))
