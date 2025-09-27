FROM alpine:latest
WORKDIR /aci-dash
#USER aci-dash
COPY aci_dash aci_dash
COPY utility utility
COPY .pyup.yml .git ./
COPY requirements.txt requirements_aci_dash.txt
RUN apk update \
    && apk upgrade \
    && apk add python3 py3-pip py3-pandas git nodejs npm \
    && adduser -D -s /bin/bash aci-dash \
    && pip install --upgrade pip \
    && npm install npm -g \
    && npm install plotly.js-dist \
    && pip install --no-cache-dir --break-system-packages dash.ly --upgrade \
    && pip install --break-system-packages gunicorn wheel \
    && pip install --no-cache-dir --break-system-packages -r requirements_aci_dash.txt
CMD ["gunicorn", "aci_dash.app:server", "--bind", "0.0.0.0:8080", "--timeout", "300"]
ENV VERSION="$(git describe --tags --long)"
ENV PATH="~/home/aci-dash/.local/bin:${PATH}"
