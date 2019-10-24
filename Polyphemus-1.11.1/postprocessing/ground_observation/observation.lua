observation_directory = "/stockage/polyphemus/common/obs/2001-v2/O3/emep/"

ground_network_observation = {

   option = {

      -- Should all observations should be loaded in memory?
      load_all_observation = true

   },

   network = {

      network_type = "default",
      stations = observation_directory .. "stations-default_format",
      -- Path to the station data file(s). The special markup "&l" is
      -- replaced with the short name of the station.
      station_data_file = observation_directory .. "&l",

      frequency = "hourly",

      unknown_values = {-999.}

   },

   domain = {

      date_begin = "2001-01-01_01",
      date_end = "2001-12-01_00",
      Delta_t = 3600,

      x_min = -10.5,
      Delta_x = 0.5,
      Nx = 67,
      y_min = 35.,
      Delta_y = 0.5,
      Ny = 46

   },

   error = {

      -- Variance of observational errors.
      variance = 49.

   }

}
