FROM python:3.10

ENV POETRY_VERSION=1.2.2
# Keeps Python from generating .pyc files in the container
ENV PYTHONDONTWRITEBYTECODE 1
# Turns off buffering for easier container logging
ENV PYTHONUNBUFFERED 1

RUN apt-get update -y && \
    apt-get install -y libhts-dev libboost-random-dev

RUN pip install "poetry==$POETRY_VERSION"

COPY pyproject.toml poetry.lock /usr/app/

# https://stackoverflow.com/questions/53835198/integrating-python-poetry-with-docker
RUN cd /usr/app && poetry install --no-interaction --no-ansi --no-root

# GitHub Actions chimes in here and sets docker's WORKDIR=${GITHUB_WORKSPACE}
# https://docs.github.com/en/actions/creating-actions/dockerfile-support-for-github-actions#workdir
# poetry install --only-root would be more elegant but does not work in Github Actions
CMD poetry install --no-interaction --no-ansi && cd ./tests && poetry run pytest