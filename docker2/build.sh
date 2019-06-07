set -ex

USERNAME=dekkerlab
IMAGE=cworld_bioperl

cp -r ../lib ./
cp -r ../scripts ./
cp ../Build.PL ./Build.PL
cp ../VERSION ./VERSION
cp ../cworld_environment.yml ./cworld_environment.yml
cp ../MANIFEST ./MANIFEST

function cleanup {
    rm  ./VERSION
    rm  ./MANIFEST
    rm  -r ./lib
    rm  -r ./scripts
    rm  ./Build.PL
    rm  ./cworld_environment.yml
}

trap cleanup EXIT

docker build -t $USERNAME/$IMAGE:latest . 
