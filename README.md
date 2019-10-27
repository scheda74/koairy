# [EM-ViZ] Documentation
## Documentation for EM-ViZ - an emission visualization tool based on SUMO, FastAPI and React

## Installation

## Useful Links

## References

## Important Commands:

Remove images that are not linked to a container
```bash
docker rmi -f $(docker images | grep '^<none>' | awk '{print $3}')
```

Remove old and exited container (CAUTION)
```bash
docker rm $(docker ps -a -f status=exited -q)
```