FROM ubuntu:20.04

RUN apt-get update -y && \
    apt-get install -y autoconf automake m4 python3 python3-pip \
        libhts-dev libboost-math-dev libboost-random-dev zlib1g-dev && \ 
    DEBIAN_FRONTEND=noninteractive apt-get install -y pkg-config 

SHELL ["/bin/bash", "-c"]

COPY requirements.txt .
RUN pip3 install -r requirements.txt

COPY . ./usr/app

# GitHub Actions chimes in here and sets docker's WORKDIR=${GITHUB_WORKSPACE}
# https://docs.github.com/en/actions/creating-actions/dockerfile-support-for-github-actions#workdir

ENTRYPOINT ["./entrypoint.sh"]

CMD pip install pytest && cd ./tests/b2w && pytest