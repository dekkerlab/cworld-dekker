set -ex

USERNAME=dekkerlab
IMAGE=cworld_env

cp -r ../* ./
#cp ../VERSION ./VERSION
#cp ../cworld_environment.yml ./cworld_environment.yml

function cleanup {
    rm  ./VERSION
    rm  ./cworld_environment.yml
}

trap cleanup EXIT

docker build -t $USERNAME/$IMAGE:latest . 
