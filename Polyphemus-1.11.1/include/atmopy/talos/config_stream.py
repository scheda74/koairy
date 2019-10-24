# Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
#     Author(s): Vivien Mallet, Vincent Picavet
#
# This file is part of AtmoPy library, a tool for data processing and
# visualization in atmospheric sciences.
#
# AtmoPy is developed in the INRIA - ENPC joint project-team CLIME and in
# the ENPC - EDF R&D joint laboratory CEREA.
#
# AtmoPy is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# AtmoPy is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
#
# For more information, visit the AtmoPy home page:
#     http://cerea.enpc.fr/polyphemus/atmopy.html


import commands
import datetime
import miscellaneous
import os


class ConfigStream:
    """Manages a configuration file. The ConfigStream class provides
    an interface to the Talos config file system, through calls to
    the extract_configuration program."""

    def __init__(self, file):
        import os.path
        self.filename = file
        self.extract = os.path.dirname(os.path.abspath(__file__)) \
                       + "/extract_configuration"
        if os.name == "nt":
            self.extract += ".exe"

        # Tests that "extract_configuration" has been built.
        if not os.path.isfile(self.extract):
            raise Exception, \
                '"' + self.extract + '" does not exist (Did you build it?).'

        # Tests that "extract_configuration" is okay and that the configuration
        # file is accessible.
        self.GetOutput("-t " + file)

    def GetOutput(self, command):
        """ Calls external program extract_configuration (which must
        be in the PATH), and returns output of execution on success.
        Raises exception on failure.

        @type command: string
        @param command: Options string to pass to extract_configuration
        program.

        @rtype: string
        @return: extract_configuration Output and error messages, or
        launch an exception if an error occured.
        """
        cmd = " ".join([self.extract, self.filename, command])
        error = None
        if os.name == "nt":
            import popen2
            o, w, e = popen2.popen3(cmd)
            errors = e.readlines()
            output = o.readlines()
            e.close()
            o.close()
            w.close()
            if len(errors) == 0:
                return "".join(output)[:-1]
            else:
                error = "".join(errors)[:-1]
        else:
            (s, o) = commands.getstatusoutput(cmd)
            if s == 0:
                return o
            else:
                error = o
        if error:
            raise Exception, "While running \"" + cmd + "\",\n" \
                + "the following errors occurred:\n" + error

    def ListSections(self):
        """
        Lists all sections of config file.

        @rtype: string
        @return: String containing all sections of config file.
        """
        return self.GetOutput("")

    def ListAll(self):
        """
        Lists everything in config file.

        @rtype: string
        @return: String containing everything in config file.
        """
        return self.GetOutput("-ll")

    def ListSectionLines(self, section):
        """
        Lists everything in specified section of config file.

        @type section: string
        @param section: Name of the section to list elements from.

        @rtype: string
        @return: String containing everything in given section of config file.
        """
        return self.GetOutput("-s " + section)

    def GetElement(self, element, section = "", type = "String"):
        """
        Returns the value of a given field.

        @type element: string
        @param element: The element (field) to be found in the configuration
        file.
        @type section: string
        @param section: The section in which the element should be found, if
        any.
        @type type: string
        @param type: The type of the element to be returned. It could be:
        'Int', 'IntList', 'IntSection', 'Bool', 'Float', 'FloatList',
        'FloatSection', 'String',  'StringList', 'StringSection', 'Num',
        'NumList' or 'NumSection'. 'Num' means 'Int' or 'Float'.

        @rtype: given by 'type'
        @return: The value of the field in the configuration file.

        @warning: In case the output is a list, the section is not taken into
        account, which means that two lists against the same field name should
        not be introduced in the configuration file, even in two different
        sections.
        """
        if type == "String":
            return self.GetString(element, section)
        elif type == "StringList":
            return self.ListSectionLines(element).split('\n')[0].split()
        elif type == "StringSection":
            return self.ListSectionLines(element).split('\n')
        elif type == "Num":
            return self.GetNum(element, section)
        elif type == "NumList":
            return [miscellaneous.to_num(x) for x in \
                    self.ListSectionLines(element).split('\n')[0].split()]
        elif type == "NumSection":
            return [miscellaneous.to_num(x) for x in \
                    self.ListSectionLines(element).split('\n')]
        elif type == "Int":
            return self.GetInt(element, section)
        elif type == "IntList":
            return [int(x) for x in \
                    self.ListSectionLines(element).split('\n')[0].split()]
        elif type == "IntSection":
            return [int(x) for x in \
                    self.ListSectionLines(element).split('\n')]
        elif type == "Bool":
            return self.GetBool(element, section)
        elif type == "Float":
            return self.GetFloat(element, section)
        elif type == "FloatList":
            return [float(x) for x in \
                    self.ListSectionLines(element).split('\n')[0].split()]
        elif type == "FloatSection":
            return [float(x) for x in \
                    self.ListSectionLines(element).split('\n')]
        elif type == "DateTime":
            return self.GetDateTime(element, section)
        elif type == "DateTimeList":
            return [self.StringToDateTime(x) for x in \
                    self.ListSectionLines(element).split('\n')[0].split()]
        elif type == "DateTimeSection":
            return [self.StringToDateTime(x) for x in \
                    self.ListSectionLines(element).split('\n')]
        else:
            raise Exception, "Type \"" + type + "\" is unknown."

    def GetString(self, element, section = ""):
        """
        Returns the value (string) of a given field.

        @type element: string
        @param element: The element (field) to be found in the configuration
        file.
        @type section: string
        @param section: The section in which the element should be found, if
        any.

        @rtype: string
        @return: The value of the field in the configuration file.
        """
        if section == "":
            return self.GetOutput(element)
        else:
            return self.GetOutput("-s " + section + " " + element)

    def GetNum(self, element, section = ""):
        """
        Returns the value (number: int or float) of a given field.

        @type element: string
        @param element: The element (field) to be found in the configuration
        file.
        @type section: string
        @param section: The section in which the element should be found, if
        any.

        @rtype: int (preferred) or float
        @return: The value of the field in the configuration file.
        """
        if section == "":
            str = self.GetOutput(element)
        else:
            str = self.GetOutput("-s " + section + " " \
                                 + element)
        return miscellaneous.to_num(str)

    def GetFloat(self, element, section = ""):
        """
        Returns specified element value in given section, as float.

        @type element: string
        @param element: The element (field) to be found in the configuration
        file.
        @type section: string
        @param section: The section in which the element should be found, if
        any.

        @rtype: float
        @return: The specified element value in given section, as float.
        """
        if section == "":
            return float(self.GetOutput(element))
        else:
            return float(self.GetOutput("-s " + section + " " \
                                        + element))

    def GetInt(self, element, section = ""):
        """
        Returns specified element value in given section, as integer

        @type element: string
        @param element: The element (field) to be found in the configuration
        file.
        @type section: string
        @param section: The section in which the element should be found, if
        any.

        @rtype: int
        @return: The specified element value in given section, as integer.
        """
        if section == "":
            return int(self.GetOutput(element))
        else:
            return int(self.GetOutput("-s " + section + " " \
                                      + element))

    def GetBool(self, element, section = ""):
        """
        Returns the value (Boolean) of a given field.

        @type element: string
        @param element: The element (field) to be found in the configuration
        file.
        @type section: string
        @param section: The section in which the element should be found, if
        any.

        @rtype: Boolean
        @return: The value of the field in the configuration file.
        """
        if section == "":
            elt = self.GetOutput(element)
        else:
            elt = self.GetOutput("-s " + section + " " + element)
        elt = elt.lower()
        if elt == "true" or elt == "t" or elt == "y" or elt == "yes":
            return True
        elif elt == "false" or elt == "f" or elt == "n" or elt == "no":
            return False
        else:
            raise Exception, "Field \"" + element + "\" is not a Boolean " \
                  "in " + self.filename + "."

    def GetDateTime(self, element, section = ""):
        """
        Returns the value (datetime) of a given field.

        @type element: string
        @param element: The element (field) to be found in the configuration
        file.
        @type section: string
        @param section: The section in which the element should be found, if
        any.

        @rtype: datetime
        @return: The value of the field in the configuration file.

        @note: The field value should be in format YYYY, YYYYMM, YYYYMMDD,
        YYYYMMDDHH, YYYYMMDDHHMM or YYYYMMDDHHMMSS. Delimiters (any character
        except numbers) can be added around the month, the day, etc.
        """
        if section == "":
            ret_str = self.GetOutput(element)
        else:
            ret_str = self.GetOutput("-s " + section + " " + element)
        return self.StringToDateTime(ret_str)

    def StringToDateTime(self, str):
        """
        Converts a string into a datetime object.

        @type str: string
        @param str: String to be converted. It must be in format YYYY, YYYYMM,
        YYYYMMDD, YYYYMMDDHH, YYYYMMDDHHMM or YYYYMMDDHHMMSS. Delimiters (any
        character except numbers) can be added around the month, the day, etc.

        @rtype: datetime
        @return: The datetime object corresponding to 'str'.
        """
        # First filters useless characters.
        str = [x for x in str if miscellaneous.is_num(x)]
        str = reduce(lambda x, y: x + y, str)

        year = int(str[0:4])
        if len(str) > 5:
            month = int(str[4:6])
        else:
            month = 1
        if len(str) > 7:
            day = int(str[6:8])
        else:
            day = 1
        if len(str) > 9:
            hour = int(str[8:10])
        else:
            hour = 0
        if len(str) > 11:
            minute = int(str[10:12])
        else:
            minute = 0
        if len(str) > 13:
            sec = int(str[12:14])
        else:
            sec = 0
        return datetime.datetime(year, month, day, hour, minute, sec)
