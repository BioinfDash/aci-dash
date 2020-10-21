FROM python:slim-buster

RUN apt-get -y update
RUN apt-get -y upgrade
RUN apt-get -y install build-essential

RUN apt -y install nodejs
RUN apt -y install npm
RUN apt -y install git
RUN apt-get -y install gunicorn
RUN apt-get -y install python3-pandas

RUN useradd -ms /bin/bash aci-dash
RUN npm install npm@latest -g
RUN npm install plotly.js-dist


#USER aci-dash
COPY aci_dash aci_dash
COPY utility utility
COPY requirements.txt requirements.txt
COPY .pyup.yml .pyup.yml
COPY .git .git

RUN export VERSION="git describe --tags --long"
RUN pip install dash.ly --upgrade
RUN export PATH=~/home/aci-dash/.local/bin:$PATH

RUN pip install -r requirements.txt

CMD gunicorn aci_dash.app:server --bind 0.0.0.0:8080 --timeout 300