FROM nvidia/cuda:12.0.1-cudnn8-devel-ubuntu22.04

ARG USERNAME=user
ARG PROJECT_NAME=ic

RUN useradd -m $USERNAME \
    && apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y coreutils sudo git gcc g++ make cmake python3 \
    && usermod -aG sudo user \
    && echo "$USERNAME ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers

USER $USERNAME

RUN mkdir /home/$USERNAME/$PROJECT_NAME
RUN mkdir -p /home/$USERNAME/.local/bin

WORKDIR /home/$USERNAME/$PROJECT_NAME
