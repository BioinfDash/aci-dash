FROM alpine

RUN apk update
RUN apk upgrade
RUN apk add python3 py3-pip py3-pandas git nodejs npm

WORKDIR /aci-dash
RUN adduser -D -s /bin/bash aci-dash
RUN pip install --upgrade pip
RUN npm install npm -g
RUN npm install plotly.js-dist

RUN pip install --no-cache-dir dash.ly --upgrade
RUN pip install gunicorn wheel
#USER aci-dash
COPY aci_dash aci_dash
COPY utility utility
COPY .pyup.yml .pyup.yml
COPY .git .git

RUN VERSION="$(git describe --tags --long)"
ENV VERSION=$VERSION

RUN export PATH=~/home/aci-dash/.local/bin:$PATH

COPY requirements.txt requirements_aci_dash.txt
RUN pip install --no-cache-dir -r requirements_aci_dash.txt

CMD gunicorn aci_dash.app:server --bind 0.0.0.0:8080 --timeout 300
