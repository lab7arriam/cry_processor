FROM python:3-alpine

RUN apk add --no-cache make cmake gcc build-base && \
	pip install biopython && \
	mkdir /cry_processor

ENV PATH=$PATH:/cry_processor

WORKDIR /data

ADD . / /cry_processor/

RUN chmod 775 /cry_processor/cry_processor.py

CMD cry_processor.py
