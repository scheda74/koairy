polyphemus_path = "/path/to/polyphemus/"

-- This path is necessary for the script "run.py" to run.
model_configuration_file = polyphemus_path .. "processing/chimere/chimere.cfg"

output_directory = "/path/to/the/output/results/"


------------------------------------ RUN -------------------------------------


run = {

   -- Please, provide all paths as absolute paths.

   main_program_path = polyphemus_path .. "processing/chimere/chimere",
   main_configuration_file = polyphemus_path .. "processing/chimere/chimere.cfg",
   compilation_option = "",
   chimere_script = polyphemus_path .. "processing/chimere/chimere.sh",
   lammpi_fortran_compiler = "g95",
   polyphemus_path = polyphemus_path,
   verdandi_path = polyphemus_path .. "include/verdandi",

}


---------------------------------- METHOD ------------------------------------


-- Forward simulation.
forward = {

   model = {

      configuration_file = model_configuration_file

   },

   display = {

      show_iteration = true,
      show_time = true

   },

   output_saver = {

      variable_list = {},
      file = ""

   },

   output = {

      configuration = output_directory .. "/forward.lua",
      log = output_directory .. "/forward.log"

   }

}


-- Simulation with assimilation using optimal interpolation.
optimal_interpolation = {

   model = {

      configuration_file = model_configuration_file

   },

   -- Computation mode for BLUE: "vector" or "matrix".
   BLUE_computation = "vector",

   data_assimilation = {

      analyze_first_step = false,

   },

   display = {

      show_iteration = true,
      show_time = true

   },

   output_saver = {

      variable_list = {"state_forecast", "state_analysis"},
      file = output_directory .. "oi-%{name}.%{extension}",
      -- time = "step 3600"

   },

   output = {

      configuration = output_directory .. "oi.lua",
      log = output_directory .. "oi.log"

   }

}


ground_network_observation = {

  option = {

     -- Should all observations should be loaded in memory?
     load_all_observation = true

  },

  network = {

     network_type = "default",
     stations = "/path/to/observations/O3/bdqa/stations",
     -- Path to the station data file(s). The special markup "&l" is
     -- replaced with the short name of the station.
     station_data_file = "/path/to/observations/O3/bdqa-v2/&l",

     frequency = "hourly",

     unknown_values = {-999.}

  },

  domain = {

     date_begin = "2003-07-30_00",
     date_end = "2003-08-03_00",
     Delta_t = 3600,
     x_min = -14.,
     Delta_x = 0.5,
     Nx = 79,
     y_min = 35.,
     Delta_y = 0.5,
     Ny = 47

  },

  error = {

     -- Variance of observational errors.
     variance = 100.

  }

}