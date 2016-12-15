FROM ubuntu
WORKDIR /home/
RUN apt-get update
RUN apt-get install -y \
build-essential \
wget \
libgts-dev
ADD \
fast-poisson-disk.makefile \
fast-poisson-disk_make_linux.patch \
/home/
RUN make -f fast-poisson-disk.makefile
RUN apt-get install -y \
python-pil \
python-scipy
ADD \
cellgen.py \
/home/
VOLUME /data
ENTRYPOINT python cellgen.py /data
