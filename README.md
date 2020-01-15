# [KoAirY] Documentation
###### Documentation for KoAiry - an emission visualization tool based on SUMO, FastAPI and React
________
The developed software is compatible with **macOS**, **Windows 10** and **Linux Debian 18.04**. Used software licenses are **GPLv3** and **EPLv2** (SUMO).

## Prerequisites
In order to run a SUMO simulation, either XLaunch for Windows or XQuartz for macOS has to be installed.

## Installation


1. Start the Docker daemon on your operating system.
2. Go to the root of the repository *koairy/* and enter the command  ```docker create network reverse_proxy```
3. Execute command ```docker-compose build``` *Note: Due to the building of SUMO, this might take several minutes.*
4. Optional: In the repository go to path *koairy-frontend/* and enter the command ```docker build -t koairy-frontend .``` 
5. Optional: In the repository go to path *koairy-backend/* and enter the command ```docker build -t koairy-backend .``` 
6. Go to the root of the repository and enter the command ```docker-compose up```.
7. In your web browser open ```0.0.0.0:3000/koairy``` in order to see the KoAiry user interface.

If problems during building or running SUMO arise, check the SUMO documentation (*https://sumo.dlr.de/docs/Downloads.html*).

## Additional Information
### Ports
All subsystems rely on the following ports being available:

* MongoDB: port 27017
* KoAiry backend: port 8000
* SUMO TraCi: port 4040
* KoAiry frontend: port 3000

Make sure these ports are free as otherwise problems before or during usage may arise.

### Data

#### Simulation Inputs
General simulation files like the road network or additional configuration files can be found in *koairy-backend/tools/simulation/data/traffic-input/*

#### Simulation Emission Output
XML files containing the emission outputs of simulations can be found in *koairy-backend/tools/simulation/data/emission-output/*

#### MongoDB data volume
If the volume of the database needs to be exported or imported, the docker volume file path is *koairy-backend:/data/db*.

### Known Issue
#### Traffic Sensor API
Unfortunately, the backend service for fetching Bremicker's traffic sensor data is sometimes unreliable. For some reason response objects differ from request to request which leads to errors in formatting and saving. The easiest solution is to just repeatedly try to fetch the data until it is being saved to the database. After that, no communication with the sensor backend takes place anymore and data fetching should run smoothly.

#### Wind Data
As no reliable, free and performant API for wind speed and direction data was found, said data is being provided by static files located at path \\ *koairy-backend/tools/predictor/data/weather/*. \\ If up-to-date wind conditions are needed, the platform of *Deutscher Wetterdienst* (DWD) (*https://opendata.dwd.de/test/CDC/observations_germany/climate/hourly/wind/recent/*) provides various types of data sets. For sensor measurments near Kirchheim, the files of sensor id *03379} have to be downloaded into the path specified above.


# Development with macOS

## Environment Variables
You need to set the following variables in order to start the server:
On MacOS:
```
export SUMO_HOME="/usr/local/opt/sumo/share/sumo"
```

## Installation
building a docker image:

docker build -t image_name .

NOTE: if SUMO fails to compile, e.g.: “internal compiler error, killed program cc1plus” try as follows:
Your virtual machine needs more RAM! At least 4096MB. (Happens on MacOS as default is 2GB)
1) docker-machine stop
2) VBOXMan --memory 4096
3) docker-machine start

## Useful Links

## References

## Important Commands / Common Errors:
```.shell script
vboxmanage modifyvm api --memory 4096
vboxmanage modifyvm api --natpf1 "backend,tcp,,8000,,8000"
vboxmanage modifyvm api --natpf1 "frontend,tcp,,3000,,3000"
vboxmanage modifyvm api --natpf1 "db,tcp,,27017,,27017"
```

### No space left on device error
Remove images that are not linked to a container
```bash
docker rmi -f $(docker images | grep '^<none>' | awk '{print $3}')
```


Remove old and exited container (CAUTION)
```bash
docker rm $(docker ps -a -f status=exited -q)
```

### MonoDB authentication fails
If authenticaton for MonoDB fails you probably need to create a database user and set user and passord.
Hence, you first need to create an admin user.
```
db.createUser({user: "admin", pwd: "example", roles: [{role: "root", db: "admin"}]})
```
After that you need to create a user 'root' for the 'mongo-db' database.
```
db.auth('admin', 'example')
use admin
db.createUser({user: "root", pwd: "example", roles: [{role: "readWrite", db: "mongo-db"}]})
```

Accessing docker shell in MongoDB container:
```shell script
docker exec -it mongo-db mongo --username root --password example --authenticationDatabase admin
```