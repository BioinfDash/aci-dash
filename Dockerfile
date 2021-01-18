FROM quay.io/mariusdieckmann/aci-dash-baseimage:v0.1.7

#USER aci-dash
COPY aci_dash aci_dash
COPY utility utility
COPY .pyup.yml .pyup.yml
COPY .git .git

RUN export VERSION="git describe --tags --long"
ENV VERSION=$VERSION
RUN export PATH=~/home/aci-dash/.local/bin:$PATH

CMD gunicorn aci_dash.app:server --bind 0.0.0.0:8080 --timeout 300
