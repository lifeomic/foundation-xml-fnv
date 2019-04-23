FROM python:2.7-jessie

LABEL name "foundation-xml-fnv"
LABEL version "1.0.0"
LABEL maintainer "LifeOmic <development@lifeomic.com>"

RUN mkdir -p /opt/app
WORKDIR /opt/app
COPY . /opt/app
RUN pip install -r requirements.txt

ENTRYPOINT ["python", "/opt/app/src/convert.py"]
