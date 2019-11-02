# [EM-ViZ] Documentation
## Documentation for EM-ViZ - an emission visualization tool based on SUMO, FastAPI and React

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
