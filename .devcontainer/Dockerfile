FROM python:3.10.14-bookworm
COPY requirements.txt .
ARG USERNAME=sahadeva
ARG USER_UID=1001
ARG USER_GID=$USER_UID
RUN apt-get update \
    && groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    && pip install --no-cache-dir -r requirements.txt
SHELL ["/bin/bash", "-c"]
USER $USERNAME