# Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
#     Author(s): Vivien Mallet
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


import config_stream


class Config:
    """
    Instances of Config store a complete configuration extracted from a
    configuration file.
    """

    def __init__(self, filename, additional_content = [], new_content = [],
                 show_error = False):
        """
        Config constructor. It reads a set of attributes in a configuration
        file.

        @type filename: string
        @param filename: The name of the configuration file.
        @type additional_content: list of tuples of strings
        @param additional_content: Description of attributes (in addition to
        the default attributes) read in the configuration file. There are
        three or four elements in each tuple:
           0. the field name or section name (in case a whole section is put
           in the attribute), to be read in the configuration file;
           1. the section name (in which the field lies), discarded in case a
           whole section is read;
           2. the name of the attribute, optional (default: the field name);
           3. the type of the attribute: 'Int', 'IntList', 'IntSection',
           'Bool', 'Float', 'FloatList', 'FloatSection', 'String',
           'StringList, 'StringSection', 'Num, 'NumList' or 'NumSection',
           where 'Num' means 'Int' or 'Float', and 'Section' means that a
           whole section is read and returned in a list made of the section
           lines.
        @type new_content: list of tuples of strings
        @param new_content: Description of all attributes. It overwrites the
        default attributes.
        @type show_error: Boolean
        @param show_error: True if an exception may be launched when an error
        occurs, False otherwise.
        """
        self.filename = filename
        self.stream = config_stream.ConfigStream(self.filename)
        self.content = [("t_min", "", "DateTime"), \
                            ("x_min", "", "Float"), \
                            ("y_min", "", "Float"), \
                            ("time_angle_min", "", "Float"), \
                            ("day_min", "", "Float"), \
                            ("Delta_t", "", "Float"), \
                            ("Delta_x", "", "Float"), \
                            ("Delta_y", "", "Float"), \
                            ("Delta_time_angle", "", "Float"), \
                            ("Delta_t", "", "Float"), \
                            ("Nt", "", "Int"), ("Nx", "", "Int"), \
                            ("Ny", "", "Int"), ("Nz", "", "Int"), \
                            ("Ntheta", "", "Int"), ("Ndays", "", "Int"), \
                            ("[levels_heights]", "",
                             "levels", "IntList"), \
                            ("[Jlevels_height]", "", "Jlevels", "IntList"), \
                            ("file", "", "input_file", "String"), \
                            ("species", "", "String"), \
                            ("obs_dir", "", "String"), \
                            ("station_file", "", "String"), \
                            ("station_file_type", "", "String"), \
                            ("type", "[output]", "String"), \
                            ("terminal", "[output]", "String"), \
                            ("t_range", "[output]", "DateTimeList"), \
                            ("x_range", "[output]", "NumList"), \
                            ("y_range", "[output]", "NumList"), \
                            ("station", "[output]", "String"), \
                            ("file", "[output]", "output_file", "String"), \
                            ("[species_list]", "", "species_list",
                             "StringSection"), \
                            ("[dir_list]", "", "dir_list", "StringSection"), \
                            ("[file_list]", "", "file_list", "StringSection")]
        self.content.extend(additional_content)
        if len(new_content) != 0:
            self.content = new_content[:]
        self.show_error = show_error
        for x in self.content:
            self.SetAttribute(x)
        self.SetMetaAttributes()

    def SetAttribute(self, x):
        """
        Sets an attribute based on its value in the configuration file. If the
        corresponding field does not appear in the configuration files, the
        attribute is not created.

        @type x: list of lists of strings
        @param x: Each element of x contains three or four elements:
           0. the field name or section name (in case a whole section is put
           in the attribute), to be read in the configuration file;
           1. the section name (in which the field lies), discarded in case a
           whole section is read;
           2. the name of the attribute, optional (default: the field name);
           3. the type of the attribute: 'Int', 'IntList', 'IntSection',
           'Bool', 'Float', 'FloatList', 'FloatSection', 'String',
           'StringList, 'StringSection', 'Num, 'NumList' or 'NumSection',
           where 'Num' means 'Int' or 'Float', and 'Section' means that a
           whole section is read and returned in a list made of the section
           lines.
        """
        try:
            val = self.stream.GetElement(x[0], section = x[1], type = x[-1])
            if len(x) == 4:
                setattr(self, x[2], val)
            else:
                setattr(self, x[0], val)
        except Exception, e:
            if self.show_error:
                raise

    def SetMetaAttributes(self):
        """
        Adds meta-attributes based on the primary attributes. The
        meta-attributes are: origin = (t_min, y_min, x_min); Delta = (Delta_t,
        Delta_y, Delta_x); shape = (Nt, Ny, Nx).

        @note: A meta-attribute is not created in case one of its primary
        attribute is missing. But no exception is raised.
        """
        try:
            self.origin = (self.t_min, self.y_min, self.x_min)
        except:
            pass
        try:
            self.Delta = (self.Delta_t, self.Delta_y, self.Delta_x)
        except:
            pass
        try:
            self.shape = (self.Nt, self.Ny, self.Nx)
        except:
            pass


def write_configuration_file(config, filename):
    """
    Writes a configuration file.

    @type config: list of list.
    @param config: The configuration is a list of fields which are described
    by a list with two or three element:
       0. The field name (string);
       1. The section (string, optional);
       2. The value of the field (any type that can be converted to a string).
    @type filename: string
    @param filename: configuration-file name.
    """
    # First sorts all fields per section.
    sorted_config = {"": []}
    for field in config:
        # Determines the delimiter: semi-colon for strings, equal sign
        # for other types.
        output = field[-1]
        if isinstance(output, str):
            output = ": " + output
        else:
            output = " = " + str(output)
        output = field[0] + output
        # If no section is specified.
        if len(field) == 2:
            sorted_config[""] = output
        else:
            if sorted_config.has_key(field[1]):
                sorted_config[field[1]].append(output)
            else:
                sorted_config[field[1]] = [output]

    f = open(filename, 'w')

    # Writes fields outside of sections.
    for field in sorted_config.pop(""):
        f.write(field + "\n")

    # Writes sections.
    for section in sorted_config:
        f.write("\n\n" + section + "\n\n")
        for field in sorted_config[section]:
            f.write(field + "\n")

    f.close()
