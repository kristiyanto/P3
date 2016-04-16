if [ -z "$AWS_ACCESS_KEY_ID"  ] || [ -z "$AWS_SECRET_ACCESS_KEY" ] || [ -z "$AWS_VPC_ID" || node]; then
  echo "set vars AWS_ACCESS_KEY_ID , AWS_SECRET_ACCESS_KEY and AWS_VPC_ID"
  exit 1;
fi
echo "Creating Multi Host Keystore"
docker-machine create --driver amazonec2 --amazonec2-access-key $AWS_ACCESS_KEY_ID --amazonec2-secret-key $AWS_SECRET_ACCESS_KEY --amazonec2-vpc-id $AWS_VPC_ID --engine-opt dns=8.8.8.8 aws-mh-keystore
eval "$(docker-machine env aws-mh-keystore)"
echo "Starting Consul at Keystore Machine"
docker run -d -p "8500:8500" -h "consul"  progrium/consul -server -bootstrap
echo "Now its time to log in https://console.aws.amazon.com/ec2/ and setup the 'docker-machine' group inbound rules."
echo "DON'T PROCEED UNTIL inbound rules is configured"
echo "If the 'docker-machine' group is configured, press any key to continue..."
read key
echo "Creating Swarm master ..."
docker-machine create --driver amazonec2 --amazonec2-access-key $AWS_ACCESS_KEY_ID --amazonec2-secret-key $AWS_SECRET_ACCESS_KEY --amazonec2-vpc-id $AWS_VPC_ID --engine-opt dns=8.8.8.8 --swarm --swarm-master --swarm-strategy "spread" --swarm-discovery="consul://$(docker-machine ip aws-mh-keystore):8500" --engine-opt="cluster-store=consul://$(docker-machine ip aws-mh-keystore):8500" --engine-opt="cluster-advertise=eth0:2376" p3-swarm-master
echo "Enter the number of nodes"
read key
if [key > 1]; then
	for i in `seq 1 $key`;
	do
		echo "Creating Swarm node $key"
		docker-machine create --driver amazonec2 --amazonec2-access-key $AWS_ACCESS_KEY_ID --amazonec2-secret-key $AWS_SECRET_ACCESS_KEY --amazonec2-vpc-id $AWS_VPC_ID --engine-opt dns=8.8.8.8 --swarm --swarm-discovery="consul://$(docker-machine ip aws-mh-keystore):8500" --engine-opt="cluster-store=consul://$(docker-machine ip aws-mh-keystore):8500" --engine-opt="cluster-advertise=eth0:2376"  p3-node$key
	done
else
	echo "Invalid key. Creating 1 node."
		docker-machine create --driver amazonec2 --amazonec2-access-key $AWS_ACCESS_KEY_ID --amazonec2-secret-key $AWS_SECRET_ACCESS_KEY --amazonec2-vpc-id $AWS_VPC_ID --engine-opt dns=8.8.8.8 --swarm --swarm-discovery="consul://$(docker-machine ip aws-mh-keystore):8500" --engine-opt="cluster-store=consul://$(docker-machine ip aws-mh-keystore):8500" --engine-opt="cluster-advertise=eth0:2376"  p3-node
fi
