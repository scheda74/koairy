# [EM-ViZ] Documentation
## Documentation for EM-ViZ - an emission visualization tool based on SUMO, FastAPI and React

## Starting the Backend

# Environment Variables
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
use mongo-db
db.createUser({user: "root", pwd: "example", roles: [{role: "readWrite", db: "mongo-db"}]})
```
