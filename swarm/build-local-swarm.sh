docker-machine ls
docker pull swarm
sid=$(docker run swarm create)
docker-machine create -d virtualbox --swarm --swarm-master --swarm-discovery token://$sid swarm-master
docker-machine create -d virtualbox --swarm --swarm-discovery token://$sid swarm-node-01
docker-machine create -d virtualbox --swarm --swarm-discovery token://$sid swarm-node-02

`docker-machine env --swarm swarm-master`
docker ps
