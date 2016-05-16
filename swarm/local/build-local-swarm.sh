# List exisiting docker engines
docker-machine ls
# Dowload and create swarm 
docker pull swarm
sid=$(docker run swarm create)

echo "Number of nodes"

read key
if key < 1 then
	echo "Invalid input. Quitting."
fi
# Create the swarm master


docker-machine create -d virtualbox --virtualbox-memory "4096" --swarm --swarm-master --swarm-discovery token://$sid swarm-master
for i in `seq 1 $key`;
do
# Create the nodes. Repeat accordingly.
	docker-machine create -d virtualbox --virtualbox-memory "4096" --swarm --swarm-discovery token://$sid swarm-node-01
	docker-machine create -d virtualbox --virtualbox-memory "4096" --swarm --swarm-discovery token://$sid swarm-node-02
done	

# List of exsisting swarm
docker run --rm swarm list token://$sid
# Get into the swarm environtment
eval `docker-machine env --swarm swarm-master`
docker ps

