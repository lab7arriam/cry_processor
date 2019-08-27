FROM python:3-alpine

RUN apk add --no-cache make cmake gcc build-base && \
	pip install biopython && \
	mkdir /cry_processor

ADD . / /cry_processor/

ENV PATH=$PATH:/cry_processor

CMD cry_processor.py
