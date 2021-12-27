FROM python:3.10

ENV POETRY_VERSION=1.1.11

RUN apt-get update -y && \
    apt-get install -y libhts-dev libboost-random-dev

RUN pip install "poetry==$POETRY_VERSION"

COPY . /usr/app/

# GitHub Actions chimes in here and sets docker's WORKDIR=${GITHUB_WORKSPACE}
# https://docs.github.com/en/actions/creating-actions/dockerfile-support-for-github-actions#workdir

ENTRYPOINT ["./entrypoint.sh"]

#CMD pip install pytest && cd ./tests && pytest