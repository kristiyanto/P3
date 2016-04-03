
# List exisiting docker engines
docker-machine ls
# Dowload and create swarm 
docker pull swarm
sid=$(docker run swarm create)
# Create the swarm master
docker-machine create -d virtualbox --swarm --swarm-master --swarm-discovery token://$sid swarm-master
# Create the nodes. Repeat accordingly.
docker-machine create -d virtualbox --swarm --swarm-discovery token://$sid swarm-node-01
docker-machine create -d virtualbox --swarm --swarm-discovery token://$sid swarm-node-02
# List of exsisting swarm
docker run swarm list token://$sid
# Get into the swarm environtment
eval `docker-machine env --swarm swarm-master`
docker ps

