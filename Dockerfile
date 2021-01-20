FROM quay.io/mariusdieckmann/aci-dash-baseimage:v0.1.8

#USER aci-dash
COPY aci_dash aci_dash
COPY utility utility
COPY .pyup.yml .pyup.yml
COPY .git .git

ENV VERSION="$(git describe --tags --long)"
RUN export PATH=~/home/aci-dash/.local/bin:$PATH

COPY requirements.txt requirements_aci_dash.txt
RUN pip install --no-cache-dir -r requirements_aci_dash.txt

CMD gunicorn aci_dash.app:server --bind 0.0.0.0:8080 --timeout 300
